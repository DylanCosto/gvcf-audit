mod compare;
mod compare_bed;
mod compare_report;
mod diagnostics;
mod genes;
mod intervals;
mod reference;
mod report;
mod vcf;

use clap::Parser;
use flate2::read::MultiGzDecoder;
use intervals::{Output, Scope, Sweep};
use reference::Reference;
use serde::Serialize;
use std::{
    collections::{BTreeMap, HashSet},
    fs::File,
    io::{BufRead, BufReader, Read},
    path::{Path, PathBuf},
    time::{Instant, SystemTime, UNIX_EPOCH},
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;
fn fail<T>(message: impl Into<String>) -> Result<T> {
    Err(message.into().into())
}

#[derive(Parser, Serialize)]
#[command(
    version,
    about = "Explain callability and missing evidence in a gVCF",
    after_help = "Compare existing audits with: gvcf-audit compare --help\n\nOutputs: callable.bed, unresolved.bed, regions.tsv, report.json and report.html; --gene-targets adds genes.tsv and exons.tsv.\nCallability is a quality policy for small-variant genotype evidence, not a guarantee of variant detection.\nThe audit scans the entire input, even with selected targets; compressed FASTA and BCF are not supported."
)]
pub struct Args {
    /// gVCF or VCF, plain text, gzip or BGZF
    #[arg(long)]
    gvcf: PathBuf,
    /// Uncompressed FASTA with matching .fai
    #[arg(long)]
    reference: PathBuf,
    /// Optional 0-based, half-open BED; without targets, covers every FASTA contig
    #[arg(long)]
    bed: Option<PathBuf>,
    /// Annotated targets TSV: contig, start, end, gene, exon (header required; replaces --bed)
    #[arg(long, conflicts_with = "bed")]
    gene_targets: Option<PathBuf>,
    /// A new output directory (existing directories are never overwritten)
    #[arg(long)]
    out: PathBuf,
    /// Sample name; required for multi-sample input
    #[arg(long)]
    sample: Option<String>,
    /// Minimum DP, or MIN_DP for reference blocks
    #[arg(long, default_value_t = 10.0)]
    min_dp: f64,
    /// Minimum reported genotype quality
    #[arg(long, default_value_t = 20.0)]
    min_gq: f64,
    /// Header detection by default; explicit deepvariant also enables hom-ref RefCall adaptation
    #[arg(long, default_value = "auto", value_parser = ["auto", "deepvariant", "gatk", "dragen", "generic"])]
    caller: String,
    /// Accept FILTER=. (not assessed) in addition to PASS
    #[arg(long)]
    allow_unfiltered: bool,
    /// Use DP when a reference block has no MIN_DP; DP may not be the block minimum
    #[arg(long)]
    allow_block_dp: bool,
    /// Exit with status 2 if fewer than this percent of requested bases are callable
    #[arg(long)]
    min_callable_percent: Option<f64>,
    /// Exit with status 2 if requested reference-mismatch bases exceed this limit
    #[arg(long)]
    max_reference_mismatch_bases: Option<u64>,
    /// Maximum mismatching records overlapping the scope, including masked or overlapping records
    #[arg(long)]
    max_reference_mismatch_records: Option<u64>,
    /// Maximum examples per reason and example type; 0 disables examples
    #[arg(long, default_value_t = 5, value_parser = clap::value_parser!(u32).range(0..=100))]
    max_examples: u32,
    /// Suppress progress and completion messages (errors are still printed)
    #[arg(long)]
    quiet: bool,
}

#[derive(Clone, Copy, PartialEq, Eq)]
#[repr(usize)]
pub enum State {
    CallableReference,
    CallableVariant,
    NoRecord,
    NoCall,
    PartialNoCall,
    SiteFiltered,
    GenotypeFiltered,
    QualityMissing,
    LowDepth,
    LowGq,
    UnsupportedPloidy,
    UnsupportedAllele,
    ComplexVariant,
    OverlappingRecords,
    ReferenceMismatch,
    ReferenceAmbiguous,
    FilterNotAssessed,
}
pub const STATES: [State; 17] = [
    State::CallableReference,
    State::CallableVariant,
    State::NoRecord,
    State::NoCall,
    State::PartialNoCall,
    State::SiteFiltered,
    State::GenotypeFiltered,
    State::QualityMissing,
    State::LowDepth,
    State::LowGq,
    State::UnsupportedPloidy,
    State::UnsupportedAllele,
    State::ComplexVariant,
    State::OverlappingRecords,
    State::ReferenceMismatch,
    State::ReferenceAmbiguous,
    State::FilterNotAssessed,
];
impl State {
    pub fn name(self) -> &'static str {
        [
            "callable_reference",
            "callable_variant",
            "no_record",
            "no_call",
            "partial_no_call",
            "site_filtered",
            "genotype_filtered",
            "quality_missing",
            "low_depth",
            "low_gq",
            "unsupported_ploidy",
            "unsupported_allele",
            "complex_variant",
            "overlapping_records",
            "reference_mismatch",
            "reference_ambiguous",
            "filter_not_assessed",
        ][self as usize]
    }
    pub fn callable(self) -> bool {
        matches!(self, Self::CallableReference | Self::CallableVariant)
    }
}

pub fn reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut probe = File::open(path)?;
    let mut magic = [0u8; 3];
    let n = probe.read(&mut magic)?;
    if n == 3 && magic == *b"BCF" {
        return fail("BCF is not supported; convert to VCF first");
    }
    let f = File::open(path)?;
    if magic[..2] == [31, 139] {
        Ok(Box::new(BufReader::with_capacity(
            1024 * 1024,
            MultiGzDecoder::new(BufReader::with_capacity(1024 * 1024, f)),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(1024 * 1024, f)))
    }
}

fn audit(args: &Args, dir: &Path) -> Result<bool> {
    let begin = Instant::now();
    let reference = Reference::open(&args.reference)?;
    let scope = Scope::load(
        args.bed.as_deref(),
        args.gene_targets.as_deref(),
        &reference,
    )?;
    let mut output = Output::new(scope, dir, args.max_examples as usize)?;
    let mut input = reader(&args.gvcf)?;
    let mut h = vcf::Header::default();
    let mut line = String::new();
    let mut line_no = 0u64;
    let mut header_done = false;
    let mut records = 0u64;
    let mut record_states = BTreeMap::new();
    let mut adaptations = BTreeMap::new();
    let mut seen = HashSet::new();
    let mut sweep: Option<Sweep> = None;
    loop {
        line.clear();
        if input.read_line(&mut line)? == 0 {
            break;
        }
        line_no += 1;
        let s = line.trim_end_matches(['\r', '\n']);
        if line_no == 1 && !s.starts_with("##fileformat=VCFv4.") {
            return fail("Input must start with a VCFv4 fileformat header");
        }
        if !header_done {
            header_done = h
                .line(s, args, &reference)
                .map_err(|e| format!("Line {line_no}: {e}"))?;
            continue;
        }
        if s.is_empty() {
            continue;
        }
        if s.starts_with('#') {
            return fail(format!("Line {line_no}: header inside VCF records"));
        }
        let r = vcf::parse(s, &h, args, &reference, &mut adaptations)
            .map_err(|e| format!("Line {line_no}: {e}"))?;
        if output.diagnostics.needs_record(&r) && output.scope.overlaps(r.contig, r.start, r.end) {
            output.diagnostics.record(&r, line_no, args, &reference)?;
        }
        if sweep.as_ref().is_none_or(|sw| sw.contig != r.contig) {
            if let Some(mut previous) = sweep.take() {
                previous.advance(reference.contigs[previous.contig].length, &mut output)?;
            }
            if !seen.insert(r.contig) {
                return fail(format!(
                    "Line {line_no}: contig {} appears in multiple groups",
                    r.chrom
                ));
            }
            sweep = Some(Sweep::new(r.contig, r.chrom.into(), &output, &reference)?);
        }
        let sw = sweep.as_mut().unwrap();
        if sw.raw_name != r.chrom {
            return fail(format!(
                "Line {line_no}: two contig names resolve to the same reference"
            ));
        }
        sw.add(r.start, r.end, r.state, &mut output)
            .map_err(|e| format!("Line {line_no}: {e}"))?;
        records += 1;
        *record_states
            .entry(r.state.name().to_string())
            .or_insert(0u64) += 1;
        if !args.quiet && records.is_multiple_of(5_000_000) {
            eprintln!(
                "Read {records} records ({:.1}s)",
                begin.elapsed().as_secs_f64()
            );
        }
    }
    if !header_done {
        return fail("No #CHROM sample header found");
    }
    if let Some(mut sw) = sweep {
        sw.advance(reference.contigs[sw.contig].length, &mut output)?;
    }
    for (id, c) in reference.contigs.iter().enumerate() {
        if !seen.contains(&id) {
            Sweep::new(id, c.name.clone(), &output, &reference)?.advance(c.length, &mut output)?;
        }
    }
    output.finish(dir)?;
    let mut warnings = Vec::<String>::new();
    if h.caller == "generic" {
        warnings.push("Caller was not uniquely identified; generic rules were used.".into());
    }
    if h.declared_contigs == 0 {
        warnings.push("The VCF has no contig declarations; reference lengths could not be checked against its header.".into());
    }
    if !h.unmapped_header_contigs.is_empty() {
        warnings.push(format!("{} header contigs are absent from the reference and have no records. They were not audited; their names are listed in JSON. Every record and requested region must have a matching reference contig.", h.unmapped_header_contigs.len()));
    }
    if args.bed.is_none() && args.gene_targets.is_none() {
        warnings.push("Scope includes every FASTA contig, including contigs absent from the VCF. Use --bed for a selected assay or region set.".into());
    }
    if args.allow_block_dp {
        warnings.push("DP fallback is enabled for blocks without MIN_DP; DP is not guaranteed to be the minimum depth across the block.".into());
    }
    if h.reblocked {
        warnings.push("ReblockGVCF processing was detected. Reblocking can remove MIN_DP and reduce GQ resolution; the audit uses only the quality fields retained in this file.".into());
    }
    if !args.allow_block_dp
        && adaptations
            .get("reference_blocks_missing_MIN_DP")
            .copied()
            .unwrap_or(0)
            > 0
    {
        warnings.push("Some reference blocks lack MIN_DP. Their minimum depth cannot be checked and otherwise passing blocks remain quality_missing. --allow-block-dp explicitly uses DP instead, which may not be a minimum.".into());
    }
    if args.allow_unfiltered {
        warnings.push(
            "FILTER=. is accepted even though site filters were not reported as passed.".into(),
        );
    }
    if record_states
        .get("reference_mismatch")
        .copied()
        .unwrap_or(0)
        > 0
    {
        warnings.push("REF mismatches were found. Check the reference assembly and input provenance before interpreting this report.".into());
    }
    if records == 0 {
        warnings.push(
            "The input contains no records. All requested bases lack genotype evidence.".into(),
        );
    }
    let quality = report::QualityChecks::evaluate(args, &output);
    let passed = quality.passed;
    report::write(
        dir,
        args,
        &h,
        &output,
        records,
        record_states,
        adaptations,
        warnings,
        begin.elapsed().as_secs_f64(),
        quality,
    )?;
    let called = output.counts[0] + output.counts[1];
    if !args.quiet {
        eprintln!(
            "{records} records; {called}/{} requested bases callable ({:.2}%). Report: {}",
            output.scope.bases(),
            100.0 * called as f64 / output.scope.bases() as f64,
            args.out.join("report.html").display()
        );
    }
    Ok(passed)
}

fn run(args: &Args) -> Result<bool> {
    if !args.min_dp.is_finite()
        || !args.min_gq.is_finite()
        || args.min_dp < 0.0
        || args.min_gq < 0.0
    {
        return fail("Depth and genotype-quality thresholds must be finite and nonnegative");
    }
    if args
        .min_callable_percent
        .is_some_and(|p| !p.is_finite() || !(0.0..=100.0).contains(&p))
    {
        return fail("--min-callable-percent must be between 0 and 100");
    }
    publish(&args.out, |dir| audit(args, dir))
}

fn publish(out: &Path, operation: impl FnOnce(&Path) -> Result<bool>) -> Result<bool> {
    if out.exists() {
        return fail("Output directory already exists; choose a new --out directory");
    }
    let parent = out
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    std::fs::create_dir_all(parent)?;
    let staging = parent.join(format!(
        ".gvcf-audit-{}-{}",
        std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
    ));
    std::fs::create_dir(&staging)?;
    let result = operation(&staging).and_then(|passed| {
        if out.exists() {
            return fail("Output directory appeared during the audit; refusing to replace it");
        }
        std::fs::rename(&staging, out)?;
        Ok(passed)
    });
    if result.is_err() {
        let _ = std::fs::remove_dir_all(&staging);
    }
    result
}

fn parsed<T>(result: std::result::Result<T, clap::Error>) -> T {
    result.unwrap_or_else(|error| {
        let code = if error.use_stderr() { 1 } else { 0 };
        let _ = error.print();
        std::process::exit(code);
    })
}

fn main() {
    let result = if std::env::args_os().nth(1).is_some_and(|s| s == "compare") {
        let args = parsed(compare::Args::try_parse_from(
            std::iter::once(std::ffi::OsString::from("gvcf-audit compare"))
                .chain(std::env::args_os().skip(2)),
        ));
        publish(&args.out, |dir| compare::run(&args, dir))
    } else {
        run(&parsed(Args::try_parse()))
    };
    match result {
        Ok(true) => {}
        Ok(false) => {
            eprintln!("Quality limits were not met. See the completed report for details.");
            std::process::exit(2);
        }
        Err(e) => {
            eprintln!("Error: {e}");
            std::process::exit(1);
        }
    }
}
