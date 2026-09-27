use crate::{
    Result, STATES,
    compare_bed::{BedIndex, Partition},
    fail,
};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::{collections::BTreeMap, fs::File, io::BufReader, path::Path};

pub(crate) type Counts = [u64; STATES.len()];
pub(crate) type Scope = BTreeMap<String, Vec<(u64, u64)>>;
pub(crate) type Key = (String, String, Option<String>);
pub(crate) type GroupMap<'a> = BTreeMap<Key, (&'a Group, Vec<(u64, u64)>)>;

#[derive(Deserialize, Serialize, PartialEq)]
pub(crate) struct Provenance {
    pub(crate) path: String,
    pub(crate) bytes: u64,
    pub(crate) modified_unix_seconds: u64,
}
#[derive(Deserialize, Serialize, PartialEq)]
pub(crate) struct Inputs {
    pub(crate) reference: Provenance,
    pub(crate) reference_index: Provenance,
}
#[derive(Deserialize, Serialize)]
pub(crate) struct Header {
    pub(crate) sample: String,
    pub(crate) caller: String,
    #[serde(default)]
    pub(crate) reblocked: bool,
}
#[derive(Deserialize, Serialize, PartialEq)]
pub(crate) struct Options {
    pub(crate) min_dp: f64,
    pub(crate) min_gq: f64,
    pub(crate) allow_block_dp: bool,
    pub(crate) allow_unfiltered: bool,
}
#[derive(Deserialize)]
pub(crate) struct Region {
    pub(crate) contig: String,
    pub(crate) start: u64,
    pub(crate) end: u64,
}
#[derive(Deserialize)]
pub(crate) struct Contig {
    pub(crate) contig: String,
    pub(crate) requested_bases: u64,
    pub(crate) callable_bases: u64,
    pub(crate) bases_by_state: BTreeMap<String, u64>,
}
#[derive(Deserialize)]
pub(crate) struct Group {
    pub(crate) contig: String,
    pub(crate) gene: String,
    pub(crate) exon: Option<String>,
    pub(crate) intervals: Vec<(u64, u64)>,
    pub(crate) bases: u64,
    pub(crate) callable_bases: u64,
    pub(crate) bases_by_state: BTreeMap<String, u64>,
}
#[derive(Deserialize)]
pub(crate) struct Audit {
    pub(crate) schema: String,
    pub(crate) version: String,
    pub(crate) status: String,
    pub(crate) inputs: Inputs,
    pub(crate) header: Header,
    pub(crate) options: Options,
    pub(crate) policy: Value,
    pub(crate) quality_checks: Value,
    pub(crate) warnings: Vec<String>,
    pub(crate) requested_bases: u64,
    pub(crate) callable_bases: u64,
    pub(crate) bases_by_state: BTreeMap<String, u64>,
    pub(crate) regions: Vec<Region>,
    pub(crate) contigs: Vec<Contig>,
    #[serde(default)]
    pub(crate) genes: Vec<Group>,
    #[serde(default)]
    pub(crate) exons: Vec<Group>,
}

pub(crate) fn counts(map: &BTreeMap<String, u64>, total: u64, callable: u64) -> Result<Counts> {
    let mut out = [0; STATES.len()];
    if map.len() != STATES.len() {
        return fail("Report contains missing or unsupported states");
    }
    for (i, state) in STATES.iter().enumerate() {
        out[i] = *map.get(state.name()).ok_or("Report state is missing")?;
    }
    if total > i64::MAX as u64
        || out.iter().try_fold(0u64, |s, &n| s.checked_add(n)) != Some(total)
        || out[0].checked_add(out[1]) != Some(callable)
    {
        return fail("Report counts are inconsistent or too large");
    }
    Ok(out)
}
pub(crate) fn merge(intervals: &mut Vec<(u64, u64)>) -> Result<()> {
    intervals.sort_unstable();
    let mut merged: Vec<(u64, u64)> = Vec::new();
    for &(s, e) in intervals.iter() {
        if s >= e {
            return fail("Report contains an empty or reversed interval");
        }
        if let Some(last) = merged.last_mut().filter(|last| s <= last.1) {
            last.1 = last.1.max(e);
        } else {
            merged.push((s, e));
        }
    }
    *intervals = merged;
    Ok(())
}
pub(crate) fn bases(intervals: &[(u64, u64)]) -> Result<u64> {
    intervals.iter().try_fold(0u64, |n, &(s, e)| {
        n.checked_add(e - s)
            .ok_or_else(|| "Target length overflow".into())
    })
}
impl Audit {
    pub(crate) fn validate_beds(&self, dir: &Path, scope: &Scope) -> Result<()> {
        let contigs = self.contig_counts(scope)?;
        let groups = self.groups()?;
        let mut targets: BTreeMap<String, Vec<(u64, u64, usize)>> = BTreeMap::new();
        let mut expected = Vec::new();
        for ((contig, _, _), (group, intervals)) in groups {
            if !scope.contains_key(&contig) {
                return fail("Gene/exon contig is outside the scope");
            }
            for (start, end) in intervals {
                targets
                    .entry(contig.clone())
                    .or_default()
                    .push((start, end, expected.len()));
            }
            expected.push(counts(
                &group.bases_by_state,
                group.bases,
                group.callable_bases,
            )?);
        }
        for intervals in targets.values_mut() {
            intervals.sort_unstable();
        }
        let beds = [
            BedIndex::open(&dir.join("callable.bed"), true, scope)?,
            BedIndex::open(&dir.join("unresolved.bed"), false, scope)?,
        ];
        let mut observed = vec![[0u64; STATES.len()]; expected.len()];
        for (contig, scope_intervals) in scope {
            let mut partition = Partition::new(&beds, contig, scope_intervals)?;
            let intervals = targets.get(contig).map(Vec::as_slice).unwrap_or(&[]);
            let mut next = 0;
            let mut active = Vec::new();
            let mut total = [0u64; STATES.len()];
            while let Some(segment) = partition.next()? {
                total[segment.state] += segment.end - segment.start;
                while next < intervals.len() && intervals[next].0 < segment.end {
                    active.push(intervals[next]);
                    next += 1;
                }
                active.retain(|&(_, end, _)| end > segment.start);
                for &(start, end, i) in &active {
                    observed[i][segment.state] += end.min(segment.end) - start.max(segment.start);
                }
            }
            if total != contigs[contig] {
                return fail("BED states disagree with reported contig counts");
            }
        }
        if observed != expected {
            return fail("BED states disagree with reported gene/exon counts");
        }
        Ok(())
    }

    pub(crate) fn load(dir: &Path) -> Result<Self> {
        let report: Self =
            serde_json::from_reader(BufReader::new(File::open(dir.join("report.json"))?))?;
        if report.schema != "gvcf-audit-report-v1"
            || report.status != "complete"
            || report.policy["id"] != "small-variant-evidence-v1"
        {
            return fail("Expected complete reports with the supported v1 callability policy");
        }
        counts(
            &report.bases_by_state,
            report.requested_bases,
            report.callable_bases,
        )?;
        if report.requested_bases == 0
            || report.header.sample.is_empty()
            || report.header.caller.is_empty()
            || report.quality_checks["passed"].as_bool().is_none()
            || !report.options.min_dp.is_finite()
            || !report.options.min_gq.is_finite()
            || report.options.min_dp < 0.0
            || report.options.min_gq < 0.0
        {
            return fail("Invalid report metadata");
        }
        Ok(report)
    }
    pub(crate) fn scope(&self) -> Result<Scope> {
        let mut scope: Scope = BTreeMap::new();
        for r in &self.regions {
            if r.contig.is_empty() {
                return fail("Report contig is empty");
            }
            scope
                .entry(r.contig.clone())
                .or_default()
                .push((r.start, r.end));
        }
        for intervals in scope.values_mut() {
            merge(intervals)?;
        }
        let total = scope.values().try_fold(0u64, |n, v| -> Result<u64> {
            n.checked_add(bases(v)?)
                .ok_or_else(|| "Target length overflow".into())
        })?;
        if total != self.requested_bases {
            return fail("Requested bases do not match the report's target union");
        }
        Ok(scope)
    }
    pub(crate) fn contig_counts(&self, scope: &Scope) -> Result<BTreeMap<String, Counts>> {
        let mut result = BTreeMap::new();
        for c in &self.contigs {
            let intervals = scope
                .get(&c.contig)
                .ok_or("Report contig is outside scope")?;
            if bases(intervals)? != c.requested_bases || result.contains_key(&c.contig) {
                return fail("Inconsistent report contig totals");
            }
            result.insert(
                c.contig.clone(),
                counts(&c.bases_by_state, c.requested_bases, c.callable_bases)?,
            );
        }
        if result.len() != scope.len() {
            return fail("Report contig totals are missing");
        }
        let mut total = [0; STATES.len()];
        for row in result.values() {
            for i in 0..STATES.len() {
                total[i] += row[i];
            }
        }
        if total
            != counts(
                &self.bases_by_state,
                self.requested_bases,
                self.callable_bases,
            )?
        {
            return fail("Global and contig counts disagree");
        }
        Ok(result)
    }
    pub(crate) fn groups(&self) -> Result<GroupMap<'_>> {
        let mut result = BTreeMap::new();
        for (rows, is_exon) in [(&self.genes, false), (&self.exons, true)] {
            for g in rows {
                if g.gene.is_empty() || g.exon.is_some() != is_exon || g.exon.as_deref() == Some("")
                {
                    return fail("Invalid gene/exon identifiers");
                }
                let mut intervals = g.intervals.clone();
                merge(&mut intervals)?;
                if intervals.is_empty() || bases(&intervals)? != g.bases {
                    return fail("Group target lengths disagree");
                }
                counts(&g.bases_by_state, g.bases, g.callable_bases)?;
                let key = (g.contig.clone(), g.gene.clone(), g.exon.clone());
                if result.insert(key, (g, intervals)).is_some() {
                    return fail("Duplicate gene/exon group");
                }
            }
        }
        Ok(result)
    }
}
