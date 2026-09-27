use crate::{Result, STATES, fail, intervals::Scope, reader, reference::Reference};
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufRead, BufWriter, Write},
    path::Path,
};

pub struct Group {
    pub contig: String,
    pub gene: String,
    pub exon: Option<String>,
    pub intervals: Vec<(u64, u64)>,
    pub bases_by_state: [u64; STATES.len()],
}

type GroupKey = (usize, String, Option<String>);

pub fn load(path: &Path, scope: &mut Scope, reference: &Reference) -> Result<()> {
    let mut header = false;
    let mut groups: BTreeMap<GroupKey, Vec<(u64, u64)>> = BTreeMap::new();
    for (line_no, line) in reader(path)?.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        if !header {
            if fields != ["contig", "start", "end", "gene", "exon"] {
                return fail(
                    "Gene targets require the tab-separated header: contig, start, end, gene, exon",
                );
            }
            header = true;
            continue;
        }
        if fields.len() != 5 || fields[3..].iter().any(|s| s.trim().is_empty() || *s == ".") {
            return fail(format!(
                "Gene targets line {}: expected five columns and nonempty gene/exon identifiers",
                line_no + 1
            ));
        }
        let id = reference.resolve(fields[0])?;
        let start: u64 = fields[1].parse()?;
        let end: u64 = fields[2].parse()?;
        if start >= end || end > reference.contigs[id].length {
            return fail(format!(
                "Gene targets line {}: empty or out-of-bounds interval",
                line_no + 1
            ));
        }
        let gene = fields[3].to_string();
        let exon = fields[4].to_string();
        let i = scope.add(
            id,
            reference.contigs[id].name.clone(),
            start,
            end,
            format!("{gene}/{exon}"),
        );
        scope.regions[i].gene = Some(gene.clone());
        scope.regions[i].exon = Some(exon.clone());
        for exon_key in [None, Some(exon)] {
            groups
                .entry((id, gene.clone(), exon_key))
                .or_default()
                .push((start, end));
        }
    }
    if scope.regions.is_empty() {
        return fail("The gene targets file contains no intervals");
    }
    for ((id, gene, exon), mut intervals) in groups {
        intervals.sort_unstable();
        let mut merged: Vec<(u64, u64)> = Vec::new();
        for (start, end) in intervals {
            if let Some(last) = merged.last_mut().filter(|last| start <= last.1) {
                last.1 = last.1.max(end);
            } else {
                merged.push((start, end));
            }
        }
        let index = scope.groups.len();
        for &(start, end) in &merged {
            scope.group_intervals[id].push((start, end, index));
        }
        scope.groups.push(Group {
            contig: reference.contigs[id].name.clone(),
            gene,
            exon,
            intervals: merged,
            bases_by_state: [0; STATES.len()],
        });
    }
    for intervals in &mut scope.group_intervals {
        intervals.sort_unstable();
    }
    Ok(())
}

pub fn summaries(groups: &[Group]) -> (Vec<Value>, Vec<Value>) {
    let mut exon_counts: BTreeMap<(&str, &str), (usize, usize)> = BTreeMap::new();
    for g in groups.iter().filter(|g| g.exon.is_some()) {
        let counts = exon_counts.entry((&g.contig, &g.gene)).or_default();
        counts.0 += 1;
        counts.1 += usize::from(g.bases_by_state[2..].iter().sum::<u64>() > 0);
    }
    let mut genes = Vec::new();
    let mut exons = Vec::new();
    for g in groups {
        let bases: u64 = g.intervals.iter().map(|(s, e)| e - s).sum();
        let callable = g.bases_by_state[0] + g.bases_by_state[1];
        let states: BTreeMap<_, _> = STATES
            .iter()
            .zip(&g.bases_by_state)
            .map(|(s, &n)| (s.name(), n))
            .collect();
        let mut main_reason = None;
        let mut most = 0;
        for (state, &n) in STATES.iter().zip(&g.bases_by_state) {
            if !state.callable() && n > most {
                main_reason = Some(state.name());
                most = n;
            }
        }
        let mut row = json!({"contig": g.contig, "gene": g.gene, "intervals": g.intervals,
            "bases": bases, "callable_bases": callable, "unresolved_bases": bases-callable,
            "callable_fraction": callable as f64 / bases as f64, "main_reason": main_reason, "bases_by_state": states});
        if let Some(exon) = &g.exon {
            row["exon"] = json!(exon);
            exons.push(row);
        } else {
            let (total, unresolved) = exon_counts[&(g.contig.as_str(), g.gene.as_str())];
            row["targeted_exons"] = json!(total);
            row["unresolved_exons"] = json!(unresolved);
            genes.push(row);
        }
    }
    (genes, exons)
}

pub fn write(dir: &Path, groups: &[Group]) -> Result<()> {
    for g in groups {
        if g.bases_by_state.iter().sum::<u64>()
            != g.intervals.iter().map(|(s, e)| e - s).sum::<u64>()
        {
            return fail("Internal gene/exon accounting error");
        }
    }
    let (genes, exons) = summaries(groups);
    for (filename, rows, is_exon) in [("genes.tsv", genes, false), ("exons.tsv", exons, true)] {
        let mut f = BufWriter::new(File::create(dir.join(filename))?);
        write!(
            f,
            "contig\tgene\t{}\tbases\tcallable_bases\tunresolved_bases\tcallable_percent\tmain_reason",
            if is_exon {
                "exon\tintervals"
            } else {
                "targeted_exons\tunresolved_exons"
            }
        )?;
        for s in STATES {
            write!(f, "\t{}", s.name())?;
        }
        writeln!(f)?;
        for r in rows {
            write!(
                f,
                "{}\t{}\t",
                r["contig"].as_str().unwrap(),
                r["gene"].as_str().unwrap()
            )?;
            if is_exon {
                let intervals = r["intervals"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| format!("{}-{}", v[0], v[1]))
                    .collect::<Vec<_>>()
                    .join(",");
                write!(f, "{}\t{intervals}", r["exon"].as_str().unwrap())?;
            } else {
                write!(f, "{}\t{}", r["targeted_exons"], r["unresolved_exons"])?;
            }
            write!(
                f,
                "\t{}\t{}\t{}\t{:.4}\t{}",
                r["bases"],
                r["callable_bases"],
                r["unresolved_bases"],
                100.0 * r["callable_fraction"].as_f64().unwrap(),
                r["main_reason"].as_str().unwrap_or(".")
            )?;
            for s in STATES {
                write!(f, "\t{}", r["bases_by_state"][s.name()])?;
            }
            writeln!(f)?;
        }
        f.flush()?;
    }
    Ok(())
}
