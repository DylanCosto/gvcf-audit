use crate::{
    Result, STATES,
    audit_report::{Audit, Counts, Key, counts},
    compare_bed::{BedIndex, Partition},
    compare_report, fail,
};
use clap::Parser;
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

#[derive(Parser)]
#[command(
    name = "gvcf-audit compare",
    version,
    about = "Compare two complete audit directories without rescanning their gVCFs"
)]
pub struct Args {
    #[arg(long)]
    before: PathBuf,
    #[arg(long)]
    after: PathBuf,
    /// A new directory for comparison.json, changes.bed, gene/exon TSVs and report.html
    #[arg(long)]
    pub out: PathBuf,
    /// Compare different thresholds or caller policies; record the differences prominently
    #[arg(long)]
    allow_policy_change: bool,
    /// Compare different sample names; this is not a within-sample pipeline comparison
    #[arg(long)]
    allow_different_samples: bool,
    /// Assert the same reference sequence was used when recorded reference metadata differ
    #[arg(long)]
    assume_same_reference: bool,
}

struct GroupChange {
    key: Key,
    intervals: Vec<(u64, u64)>,
    before: Counts,
    after: Counts,
    observed_before: Counts,
    observed_after: Counts,
    gained: u64,
    lost: u64,
}
fn named(values: &Counts) -> BTreeMap<&'static str, u64> {
    STATES
        .iter()
        .zip(values)
        .map(|(s, &v)| (s.name(), v))
        .collect()
}
fn summary(before: &Counts, after: &Counts, gained: u64, lost: u64) -> Value {
    let total: u64 = before.iter().sum();
    let old = before[0] + before[1];
    let new = after[0] + after[1];
    json!({"bases":total,"before_callable_bases":old,"after_callable_bases":new,
        "gained_callable_bases":gained,"lost_callable_bases":lost,"net_callable_bases":new as i64-old as i64,
        "before_callable_fraction":old as f64/total as f64,"after_callable_fraction":new as f64/total as f64,
        "callable_percentage_point_change":100.0*(new as f64-old as f64)/total as f64,
        "before_bases_by_state":named(before),"after_bases_by_state":named(after),
        "bases_by_state_delta":STATES.iter().enumerate().map(|(i,s)|(s.name(),after[i] as i64-before[i] as i64)).collect::<BTreeMap<_,_>>()})
}
fn compatibility(args: &Args, before: &Audit, after: &Audit) -> Result<Vec<String>> {
    let mut warnings = vec!["Reference sequence identity was not verified: reports contain paths, sizes and modification times, not sequence checksums.".into()];
    if before.header.sample != after.header.sample {
        if !args.allow_different_samples {
            return fail(
                "Sample names differ; use --allow-different-samples only for an intentional cross-sample comparison",
            );
        }
        warnings.push(format!(
            "Different samples: {} → {}. Changes cannot be attributed to the pipeline alone.",
            before.header.sample, after.header.sample
        ));
    }
    if before.inputs.reference != after.inputs.reference
        || before.inputs.reference_index != after.inputs.reference_index
    {
        if !args.assume_same_reference {
            return fail(
                "Recorded reference metadata differ; confirm both audits used identical reference sequences before using --assume-same-reference",
            );
        }
        warnings.push("Reference metadata differ. The user asserted identical reference sequences with --assume-same-reference.".into());
    }
    let mut differences = Vec::new();
    let a = serde_json::to_value(&before.options)?;
    let b = serde_json::to_value(&after.options)?;
    for field in ["min_dp", "min_gq", "allow_block_dp", "allow_unfiltered"] {
        if a[field] != b[field] {
            differences.push(format!("{field}: {} → {}", a[field], b[field]));
        }
    }
    if before.header.caller != after.header.caller {
        differences.push(format!(
            "caller: {} → {}",
            before.header.caller, after.header.caller
        ));
    }
    if before.policy != after.policy {
        differences.push("policy definitions differ".into());
    }
    if !differences.is_empty() {
        if !args.allow_policy_change {
            return fail(format!(
                "Audit rules differ ({}); use --allow-policy-change for an intentional policy comparison",
                differences.join("; ")
            ));
        }
        warnings.push(format!(
            "Rules changed: {}. These results include the effect of those changes.",
            differences.join("; ")
        ));
    }
    if before.header.reblocked != after.header.reblocked {
        warnings.push(
            "Reblocking status changed; the inputs may retain different quality evidence.".into(),
        );
    }
    for (label, audit) in [("Before", before), ("After", after)] {
        if audit.quality_checks["passed"] == false {
            warnings.push(format!("{label} audit failed its quality limits."));
        }
        warnings.extend(audit.warnings.iter().map(|w| format!("{label} audit: {w}")));
    }
    Ok(warnings)
}

pub fn run(args: &Args, dir: &Path) -> Result<bool> {
    let mut before = Audit::load(&args.before)?;
    let mut after = Audit::load(&args.after)?;
    let warnings = compatibility(args, &before, &after)?;
    let scope = before.scope()?;
    if scope != after.scope()? {
        return fail(
            "Target unions differ; audit both inputs with the same targets before comparing",
        );
    }
    let old_contigs = before.contig_counts(&scope)?;
    let new_contigs = after.contig_counts(&scope)?;
    let old_groups = before.groups()?;
    let new_groups = after.groups()?;
    if old_groups.keys().ne(new_groups.keys()) {
        return fail("Gene/exon identifiers differ; use the same annotation for both audits");
    }
    let mut groups = Vec::new();
    let mut group_intervals: BTreeMap<String, Vec<(u64, u64, usize)>> = BTreeMap::new();
    for (key, (old, intervals)) in old_groups {
        let (new, new_intervals) = &new_groups[&key];
        if &intervals != new_intervals {
            return fail(
                "Gene/exon target intervals differ; use the same annotation for both audits",
            );
        }
        if !scope.contains_key(&key.0) {
            return fail("Gene/exon contig is outside the scope");
        }
        let index = groups.len();
        for &(s, e) in &intervals {
            group_intervals
                .entry(key.0.clone())
                .or_default()
                .push((s, e, index));
        }
        groups.push(GroupChange {
            key,
            intervals,
            before: counts(&old.bases_by_state, old.bases, old.callable_bases)?,
            after: counts(&new.bases_by_state, new.bases, new.callable_bases)?,
            observed_before: [0; STATES.len()],
            observed_after: [0; STATES.len()],
            gained: 0,
            lost: 0,
        });
    }
    drop(new_groups);
    before.genes.clear();
    before.exons.clear();
    before.regions.clear();
    after.genes.clear();
    after.exons.clear();
    after.regions.clear();
    for entries in group_intervals.values_mut() {
        entries.sort_unstable();
    }
    let indices = |p: &Path| -> Result<[BedIndex; 2]> {
        Ok([
            BedIndex::open(&p.join("callable.bed"), true, &scope)?,
            BedIndex::open(&p.join("unresolved.bed"), false, &scope)?,
        ])
    };
    let old_beds = indices(&args.before)?;
    let new_beds = indices(&args.after)?;
    let mut changes = BufWriter::new(File::create(dir.join("changes.bed"))?);
    let mut transitions = [[0u64; STATES.len()]; STATES.len()];
    let mut contigs = Vec::new();
    let mut gained_total = 0;
    let mut lost_total = 0;
    for (contig, intervals) in &scope {
        let mut a = Partition::new(&old_beds, contig, intervals)?;
        let mut b = Partition::new(&new_beds, contig, intervals)?;
        let mut left = a.next()?;
        let mut right = b.next()?;
        let mut old = [0; STATES.len()];
        let mut new = [0; STATES.len()];
        let mut gained = 0;
        let mut lost = 0;
        let targets = group_intervals
            .get(contig)
            .map(Vec::as_slice)
            .unwrap_or(&[]);
        let mut next = 0;
        let mut active = Vec::new();
        while let (Some(x), Some(y)) = (&mut left, &mut right) {
            if x.start != y.start {
                return fail("Audit BED partitions disagree");
            }
            let start = x.start;
            let end = x.end.min(y.end);
            let n = end - start;
            old[x.state] += n;
            new[y.state] += n;
            transitions[x.state][y.state] += n;
            let gain =
                usize::from(!STATES[x.state].callable() && STATES[y.state].callable()) as u64;
            let loss =
                usize::from(STATES[x.state].callable() && !STATES[y.state].callable()) as u64;
            gained += gain * n;
            lost += loss * n;
            if x.state != y.state {
                writeln!(
                    changes,
                    "{contig}\t{start}\t{end}\t{}->{}",
                    STATES[x.state].name(),
                    STATES[y.state].name()
                )?;
            }
            while next < targets.len() && targets[next].0 < end {
                active.push(targets[next]);
                next += 1;
            }
            active.retain(|&(_, e, _)| e > start);
            for &(s, e, i) in &active {
                let size = e.min(end) - s.max(start);
                let g = &mut groups[i];
                g.observed_before[x.state] += size;
                g.observed_after[y.state] += size;
                g.gained += gain * size;
                g.lost += loss * size;
            }
            if x.end == end {
                left = a.next()?;
            } else {
                x.start = end;
            }
            if y.end == end {
                right = b.next()?;
            } else {
                y.start = end;
            }
        }
        if left.is_some()
            || right.is_some()
            || old != old_contigs[contig]
            || new != new_contigs[contig]
        {
            return fail("BED states disagree with reported contig counts");
        }
        gained_total += gained;
        lost_total += lost;
        let mut row = summary(&old, &new, gained, lost);
        row["contig"] = json!(contig);
        contigs.push(row);
    }
    changes.flush()?;
    let mut genes = Vec::new();
    let mut exons = Vec::new();
    for g in groups {
        if g.before != g.observed_before || g.after != g.observed_after {
            return fail("BED states disagree with reported gene/exon counts");
        }
        let mut row = summary(&g.before, &g.after, g.gained, g.lost);
        row["contig"] = json!(g.key.0);
        row["gene"] = json!(g.key.1);
        row["intervals"] = json!(g.intervals);
        if let Some(exon) = g.key.2 {
            row["exon"] = json!(exon);
            exons.push(row);
        } else {
            genes.push(row);
        }
    }
    let source = |p: &Path, r: &Audit| -> Result<Value> {
        Ok(
            json!({"directory":std::fs::canonicalize(p)?.to_string_lossy(),"version":r.version,"header":r.header,"reference":r.inputs,"options":r.options,"quality_checks":r.quality_checks}),
        )
    };
    let transition_rows: Vec<_> = transitions
        .iter()
        .enumerate()
        .flat_map(|(i, row)| {
            row.iter().enumerate().filter(|(_, n)| **n > 0).map(move |(j, &n)| json!({"before":STATES[i].name(),"after":STATES[j].name(),"bases":n}))
        })
        .collect();
    let mut result = json!({"schema":"gvcf-audit-comparison-v1","version":env!("CARGO_PKG_VERSION"),"status":"complete",
        "before":source(&args.before,&before)?,"after":source(&args.after,&after)?,"warnings":warnings,
        "compatibility":{"same_target_union":true,"same_group_targets":true,"reference_identity_verified":false,
            "allow_policy_change":args.allow_policy_change,"allow_different_samples":args.allow_different_samples,"assume_same_reference":args.assume_same_reference},
        "summary":summary(&counts(&before.bases_by_state,before.requested_bases,before.callable_bases)?,&counts(&after.bases_by_state,after.requested_bases,after.callable_bases)?,gained_total,lost_total),
        "contigs":contigs,"transitions":transition_rows,
        "coordinates":"0-based half-open; changes.bed columns: contig, start, end, before_state->after_state"});
    result["genes"] = Value::Array(genes);
    result["exons"] = Value::Array(exons);
    compare_report::write(dir, &result)?;
    eprintln!(
        "{gained_total} bases gained callability; {lost_total} lost. Report: {}",
        args.out.join("report.html").display()
    );
    Ok(true)
}
