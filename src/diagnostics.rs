use crate::{Args, Result, STATES, State, reference::Reference, vcf::Record};
use serde::Serialize;
use serde_json::{Value, json};

#[derive(Serialize)]
pub struct IntervalExample {
    pub contig: String,
    pub start: u64,
    pub end: u64,
    pub reason: &'static str,
}

#[derive(Serialize)]
pub struct Diagnostics {
    pub max_examples_per_reason: usize,
    pub reference_mismatch_records_overlapping_scope: u64,
    pub record_examples: Vec<Value>,
    pub interval_examples: Vec<IntervalExample>,
    #[serde(skip)]
    record_counts: [usize; STATES.len()],
    #[serde(skip)]
    missing_min_dp_examples: usize,
    #[serde(skip)]
    interval_indices: [Vec<usize>; STATES.len()],
}

impl Diagnostics {
    pub fn new(limit: usize) -> Self {
        Self {
            max_examples_per_reason: limit,
            reference_mismatch_records_overlapping_scope: 0,
            record_examples: Vec::new(),
            interval_examples: Vec::new(),
            record_counts: [0; STATES.len()],
            missing_min_dp_examples: 0,
            interval_indices: std::array::from_fn(|_| Vec::new()),
        }
    }

    pub fn needs_record(&self, r: &Record<'_>) -> bool {
        r.state == State::ReferenceMismatch
            || (!r.state.callable()
                && self.record_counts[r.state as usize] < self.max_examples_per_reason)
            || (r.block
                && r.min_dp.is_none()
                && self.missing_min_dp_examples < self.max_examples_per_reason)
    }

    pub fn record(
        &mut self,
        r: &Record<'_>,
        line: u64,
        args: &Args,
        reference: &Reference,
    ) -> Result<()> {
        if r.state == State::ReferenceMismatch {
            self.reference_mismatch_records_overlapping_scope += 1;
        }
        let primary = !r.state.callable()
            && self.record_counts[r.state as usize] < self.max_examples_per_reason;
        let missing = r.block
            && r.min_dp.is_none()
            && self.missing_min_dp_examples < self.max_examples_per_reason;
        if !primary && !missing {
            return Ok(());
        }
        let depth_field = if r.block && !(args.allow_block_dp && r.min_dp.is_none()) {
            "MIN_DP"
        } else {
            "DP"
        };
        let mismatch = if r.state == State::ReferenceMismatch {
            reference.first_mismatch(r.contig, r.start, r.reference)?
                .map(|(position, vcf, fasta)| json!({"position": position, "vcf_base": vcf, "reference_base": fasta}))
        } else {
            None
        };
        let observed = json!({
            "GT": preview(r.genotype), "REF": preview(r.reference), "ALT": preview(r.alternate),
            "FILTER": preview(r.filter), "FT": r.genotype_filter.map(preview),
            "DP": r.dp, "MIN_DP": r.min_dp, "GQ": r.gq,
            "reference_block": r.block, "depth_field_used": depth_field,
            "first_reference_mismatch": mismatch
        });
        if primary {
            self.record_examples.push(json!({
                "contig": reference.contigs[r.contig].name, "input_contig": r.chrom,
                "start": r.start, "end": r.end, "line": line,
                "reason": r.state.name(), "rule": rule(r.state, depth_field, args), "observed": observed
            }));
            self.record_counts[r.state as usize] += 1;
        }
        if missing {
            self.record_examples.push(json!({
                "contig": reference.contigs[r.contig].name, "input_contig": r.chrom,
                "start": r.start, "end": r.end, "line": line,
                "reason": "reference_blocks_missing_MIN_DP",
                "rule": "Minimum block depth requires MIN_DP. DP fallback, when enabled, cannot recover that minimum.",
                "observed": observed
            }));
            self.missing_min_dp_examples += 1;
        }
        Ok(())
    }

    pub fn interval(&mut self, contig: &str, start: u64, end: u64, state: State) {
        if state.callable() || self.max_examples_per_reason == 0 {
            return;
        }
        let indices = &mut self.interval_indices[state as usize];
        if let Some(&i) = indices.last() {
            let previous = &mut self.interval_examples[i];
            if previous.contig == contig && previous.end == start {
                previous.end = end;
                return;
            }
        }
        if indices.len() < self.max_examples_per_reason {
            indices.push(self.interval_examples.len());
            self.interval_examples.push(IntervalExample {
                contig: contig.into(),
                start,
                end,
                reason: state.name(),
            });
        }
    }
}

fn preview(value: &str) -> String {
    let mut text: String = value.chars().take(80).collect();
    if text.len() < value.len() {
        text.push('…');
    }
    text
}

fn rule(state: State, depth: &str, args: &Args) -> String {
    match state {
        State::LowDepth => format!("{depth} must be at least {}.", args.min_dp),
        State::LowGq => format!("GQ must be at least {}.", args.min_gq),
        State::QualityMissing => format!("Both {depth} and GQ must be present."),
        State::ReferenceMismatch => "REF must agree with the supplied FASTA (ignoring letter case).".into(),
        State::ReferenceAmbiguous => "The REF allele must contain only A, C, G or T.".into(),
        State::NoCall | State::PartialNoCall => "Every GT allele must be called.".into(),
        State::UnsupportedPloidy => "The genotype must be haploid or diploid.".into(),
        State::GenotypeFiltered => "FORMAT/FT must be PASS, missing or absent.".into(),
        State::FilterNotAssessed => "FILTER=. requires a recognized reference-block convention or explicit --allow-unfiltered.".into(),
        State::SiteFiltered => "Site filters must pass the selected caller policy.".into(),
        State::UnsupportedAllele => "Called alleles must be literal A/C/G/T alleles; a reference block must be hom-ref.".into(),
        State::ComplexVariant => "Indels, symbolic alleles and extended variant spans remain unresolved under this policy.".into(),
        _ => "See the final interval status and documented callability policy.".into(),
    }
}
