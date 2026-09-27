use crate::{
    Result,
    audit_report::{Audit, Key, Scope},
    fail,
};
use clap::{Parser, ValueEnum};
use serde::Serialize;
use serde_json::{Value, json};
use std::{
    collections::{BTreeMap, HashSet},
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

#[derive(Clone, Copy, PartialEq, Serialize, ValueEnum)]
#[serde(rename_all = "lowercase")]
enum Level {
    Gene,
    Exon,
}

#[derive(Parser, Serialize)]
#[command(
    name = "gvcf-audit summarize",
    version,
    about = "Summarize gene or exon callability across saved audits"
)]
pub struct Args {
    /// Audit directories containing report.json, callable.bed and unresolved.bed
    #[arg(value_name = "AUDIT_DIR")]
    audits: Vec<PathBuf>,
    /// Optional file of audit directory paths, one per line, relative to this file
    #[arg(long)]
    audit_list: Option<PathBuf>,
    /// A new directory for cohort.json, sample/group summaries and callability matrices
    #[arg(long)]
    pub out: PathBuf,
    /// Summarize genes (default) or individual exon identifiers
    #[arg(long, value_enum, default_value = "gene")]
    level: Level,
    /// Count samples below this callability percentage for each gene or exon
    #[arg(long, default_value_t = 95.0)]
    min_callable_percent: f64,
    /// Assert identical reference sequences when the recorded reference metadata differ
    #[arg(long)]
    assume_same_reference: bool,
}

struct Group {
    intervals: Vec<(u64, u64)>,
    bases: u64,
    called: Vec<u64>,
}

fn identifier(value: &str) -> Result<()> {
    if value.chars().any(char::is_control) {
        return fail("Sample and annotation identifiers must not contain control characters");
    }
    Ok(())
}

fn cell(value: &str) -> String {
    if value.contains('"') {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

fn inputs(args: &Args) -> Result<Vec<PathBuf>> {
    let mut paths = args.audits.clone();
    if let Some(list) = &args.audit_list {
        let parent = list.parent().unwrap_or(Path::new("."));
        for line in BufReader::new(File::open(list)?).lines() {
            let line = line?;
            if !line.trim().is_empty() && !line.starts_with('#') {
                paths.push(parent.join(line));
            }
        }
    }
    if paths.is_empty() {
        return fail("Supply audit directories or --audit-list");
    }
    Ok(paths)
}

pub fn run(args: &Args, dir: &Path) -> Result<bool> {
    if !args.min_callable_percent.is_finite() || !(0.0..=100.0).contains(&args.min_callable_percent)
    {
        return fail("--min-callable-percent must be between 0 and 100");
    }
    let mut groups: BTreeMap<Key, Group> = BTreeMap::new();
    let mut scope: Option<Scope> = None;
    let mut reference = Value::Null;
    let mut policy = Value::Null;
    let mut names = HashSet::new();
    let mut paths_seen = HashSet::new();
    let mut samples = Vec::new();
    let mut warnings = vec!["Reference sequence identity was not verified; saved reports contain file metadata, not sequence checksums.".to_string()];
    let mut reblocked = None;
    for path in inputs(args)? {
        let path = fs::canonicalize(path)?;
        if !paths_seen.insert(path.clone()) {
            return fail(format!("Duplicate audit directory: {}", path.display()));
        }
        let audit = Audit::load(&path).map_err(|e| format!("{}: {e}", path.display()))?;
        if audit.quality_checks["requested"].as_bool().is_none() {
            return fail("Report is missing its quality-gate request status");
        }
        identifier(&audit.header.sample)?;
        if !names.insert(audit.header.sample.clone()) {
            return fail(format!("Duplicate sample name: {}", audit.header.sample));
        }
        let current_scope = audit.scope()?;
        let current_reference = serde_json::to_value(&audit.inputs)?;
        let current_policy = json!({"options": audit.options, "caller": audit.header.caller, "policy": audit.policy});
        if let Some(scope) = &scope {
            if scope != &current_scope {
                return fail(format!("{}: target unions differ", audit.header.sample));
            }
            if policy != current_policy {
                return fail(format!(
                    "{}: quality rules or caller policies differ; use matching policies for a cohort",
                    audit.header.sample
                ));
            }
            if reference != current_reference {
                if !args.assume_same_reference {
                    return fail(
                        "Reference metadata differ; use --assume-same-reference only after confirming identical reference sequences",
                    );
                }
                warnings.push(format!("{}: reference metadata differ; the user asserted identical reference sequences.", audit.header.sample));
            }
            if reblocked != Some(audit.header.reblocked) {
                warnings.push(format!(
                    "{}: reblocking status differs from the first sample.",
                    audit.header.sample
                ));
            }
        } else {
            scope = Some(current_scope.clone());
            reference = current_reference;
            policy = current_policy;
            reblocked = Some(audit.header.reblocked);
        }
        audit
            .validate_beds(&path, &current_scope)
            .map_err(|e| format!("{}: {e}", path.display()))?;
        let selected: BTreeMap<_, _> = audit
            .groups()?
            .into_iter()
            .filter(|(key, _)| key.2.is_some() == (args.level == Level::Exon))
            .collect();
        if selected.is_empty() {
            return fail("No gene/exon summaries found; generate audits with --gene-targets first");
        }
        if !samples.is_empty() && groups.keys().ne(selected.keys()) {
            return fail(format!(
                "{}: gene/exon identifiers differ",
                audit.header.sample
            ));
        }
        for (key, (group, intervals)) in selected {
            identifier(&key.0)?;
            identifier(&key.1)?;
            if let Some(exon) = &key.2 {
                identifier(exon)?;
            }
            let row = groups.entry(key).or_insert_with(|| Group {
                intervals: intervals.clone(),
                bases: group.bases,
                called: Vec::new(),
            });
            if row.intervals != intervals || row.bases != group.bases {
                return fail(format!(
                    "{}: gene/exon target intervals differ",
                    audit.header.sample
                ));
            }
            row.called.push(group.callable_bases);
        }
        samples.push(json!({"sample":audit.header.sample, "directory":path, "version":audit.version,
            "requested_bases":audit.requested_bases,"callable_bases":audit.callable_bases,
            "callable_percent":100.0 * audit.callable_bases as f64 / audit.requested_bases as f64,
            "quality_checks":audit.quality_checks,"warnings":audit.warnings,"reference":audit.inputs,
            "reblocked":audit.header.reblocked,"bases_by_state":audit.bases_by_state}));
    }
    write_tables(dir, args, &samples, &groups)?;
    let report = json!({"schema":"gvcf-audit-cohort-v1", "version":env!("CARGO_PKG_VERSION"), "status":"complete",
        "level":args.level,"sample_count":samples.len(),"group_count":groups.len(),"samples":samples,
        "policy":policy,"min_callable_percent":args.min_callable_percent,"warnings":warnings,
        "compatibility":{"same_target_union":true,"same_group_targets":true,"same_policy":true,
            "reference_identity_verified":false,"assume_same_reference":args.assume_same_reference},
        "definition":"Callability describes usable genotype evidence in the gVCF, not biological accuracy. The percentage threshold counts samples per group; it does not set an exit-code gate.",
        "files":{"percent":"callable_percent.tsv","bases":"callable_bases.tsv","groups":"groups.tsv","samples":"samples.tsv"}});
    fs::write(dir.join("cohort.json"), serde_json::to_vec_pretty(&report)?)?;
    eprintln!(
        "{} samples; {} groups. Summary: {}",
        samples.len(),
        groups.len(),
        args.out.display()
    );
    Ok(true)
}

fn write_tables(
    dir: &Path,
    args: &Args,
    samples: &[Value],
    groups: &BTreeMap<Key, Group>,
) -> Result<()> {
    let mut percents = BufWriter::new(File::create(dir.join("callable_percent.tsv"))?);
    let mut bases = BufWriter::new(File::create(dir.join("callable_bases.tsv"))?);
    let mut summary = BufWriter::new(File::create(dir.join("groups.tsv"))?);
    let mut sample_table = BufWriter::new(File::create(dir.join("samples.tsv"))?);
    let key_header = if args.level == Level::Exon {
        "contig\tgene\texon\trequested_bases"
    } else {
        "contig\tgene\trequested_bases"
    };
    for file in [&mut percents, &mut bases] {
        write!(file, "{key_header}")?;
        for sample in samples {
            write!(
                file,
                "\t{}",
                cell(&format!("sample:{}", sample["sample"].as_str().unwrap()))
            )?;
        }
        writeln!(file)?;
    }
    writeln!(
        summary,
        "{key_header}\tmin_callable_percent\tmean_callable_percent\tmax_callable_percent\tsamples_below_threshold\tsamples_fully_callable"
    )?;
    for (key, group) in groups {
        let mut label = format!("{}\t{}", cell(&key.0), cell(&key.1));
        if let Some(exon) = &key.2 {
            label.push_str(&format!("\t{}", cell(exon)));
        }
        label.push_str(&format!("\t{}", group.bases));
        write!(percents, "{label}")?;
        write!(bases, "{label}")?;
        let percent = |n: u64| 100.0 * n as f64 / group.bases as f64;
        for &count in &group.called {
            write!(percents, "\t{:.6}", percent(count))?;
            write!(bases, "\t{count}")?;
        }
        writeln!(percents)?;
        writeln!(bases)?;
        let sum: u128 = group.called.iter().map(|&n| n as u128).sum();
        let mean = 100.0 * sum as f64 / (group.bases as f64 * samples.len() as f64);
        let below = group
            .called
            .iter()
            .filter(|&&n| percent(n) < args.min_callable_percent)
            .count();
        let full = group.called.iter().filter(|&&n| n == group.bases).count();
        writeln!(
            summary,
            "{label}\t{:.6}\t{mean:.6}\t{:.6}\t{below}\t{full}",
            percent(*group.called.iter().min().unwrap()),
            percent(*group.called.iter().max().unwrap())
        )?;
    }
    writeln!(
        sample_table,
        "sample\trequested_bases\tcallable_bases\tcallable_percent\tquality_gates"
    )?;
    for sample in samples {
        let quality = &sample["quality_checks"];
        let gate = if quality["requested"] == false {
            "not requested"
        } else if quality["passed"] == true {
            "passed"
        } else {
            "failed"
        };
        writeln!(
            sample_table,
            "{}\t{}\t{}\t{:.6}\t{gate}",
            cell(sample["sample"].as_str().unwrap()),
            sample["requested_bases"],
            sample["callable_bases"],
            sample["callable_percent"].as_f64().unwrap()
        )?;
    }
    for file in [&mut percents, &mut bases, &mut summary, &mut sample_table] {
        file.flush()?;
    }
    Ok(())
}
