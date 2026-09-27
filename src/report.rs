use crate::{Args, Result, STATES, State, genes, intervals::Output, vcf::Header};
use serde::Serialize;
use serde_json::{Value, json};
use std::{collections::BTreeMap, path::Path, time::UNIX_EPOCH};

pub(crate) fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}

fn file_info(path: &Path) -> Result<Value> {
    let meta = std::fs::metadata(path)?;
    Ok(
        json!({"path": std::fs::canonicalize(path)?.to_string_lossy(), "bytes": meta.len(), "modified_unix_seconds": meta.modified()?.duration_since(UNIX_EPOCH)?.as_secs()}),
    )
}

fn counts(values: &[u64]) -> BTreeMap<&'static str, u64> {
    STATES
        .iter()
        .zip(values)
        .map(|(s, &v)| (s.name(), v))
        .collect()
}

#[derive(Serialize)]
pub struct QualityChecks {
    pub passed: bool,
    requested: bool,
    failures: Vec<String>,
}

impl QualityChecks {
    pub fn evaluate(args: &Args, out: &Output) -> Self {
        let mut failures = Vec::new();
        let percent = 100.0 * (out.counts[0] + out.counts[1]) as f64 / out.scope.bases() as f64;
        if let Some(minimum) = args.min_callable_percent
            && percent < minimum
        {
            failures.push(format!(
                "Callable bases: {percent:.6}%, required at least {minimum}%"
            ));
        }
        let mismatches = out.counts[State::ReferenceMismatch as usize];
        if let Some(maximum) = args.max_reference_mismatch_bases
            && mismatches > maximum
        {
            failures.push(format!(
                "Reference-mismatch bases: {mismatches}, allowed at most {maximum}"
            ));
        }
        if let Some(maximum) = args.max_reference_mismatch_records {
            let records = out.diagnostics.reference_mismatch_records_overlapping_scope;
            if records > maximum {
                failures.push(format!("Reference-mismatch records overlapping scope: {records}, allowed at most {maximum}"));
            }
        }
        Self {
            passed: failures.is_empty(),
            requested: args.min_callable_percent.is_some()
                || args.max_reference_mismatch_bases.is_some()
                || args.max_reference_mismatch_records.is_some(),
            failures,
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn write(
    dir: &Path,
    args: &Args,
    h: &Header,
    out: &Output,
    records: u64,
    record_states: BTreeMap<String, u64>,
    adaptations: BTreeMap<String, u64>,
    warnings: Vec<String>,
    seconds: f64,
    quality: QualityChecks,
) -> Result<()> {
    let total = out.scope.bases();
    let called = out.counts[0] + out.counts[1];
    let percent = 100.0 * called as f64 / total as f64;
    let scope = if args.gene_targets.is_some() {
        "Union of supplied gene/exon targets"
    } else if args.bed.is_some() {
        "Union of supplied BED intervals"
    } else {
        "All reference contigs"
    };
    let mut regions = Vec::new();
    for r in &out.scope.regions {
        let size = r.end - r.start;
        let c = r.bases_by_state[0] + r.bases_by_state[1];
        let mut row = json!({"contig": r.contig, "start": r.start, "end": r.end, "name": r.name, "bases": size, "callable_bases": c, "callable_fraction": c as f64 / size as f64, "bases_by_state": counts(&r.bases_by_state)});
        if let Some(gene) = &r.gene {
            row["gene"] = json!(gene);
            row["exon"] = json!(r.exon);
        }
        regions.push(row);
    }
    let mut contigs = Vec::new();
    for (id, values) in out.contig_counts.iter().enumerate() {
        let bases: u64 = values.iter().sum();
        if bases == 0 {
            continue;
        }
        let name = &out.scope.regions[out.scope.by_contig[id][0]].contig;
        let callable = values[0] + values[1];
        contigs.push(
            json!({"contig": name, "requested_bases": bases, "callable_bases": callable,
            "callable_fraction": callable as f64 / bases as f64, "bases_by_state": counts(values)}),
        );
    }
    let definition = "Callable means the input supplies haploid or diploid small-variant genotype evidence passing the stated filters, depth and GQ thresholds. It does not establish absence of every variant type or correctness of the caller. Indel/structural spans and overlapping records are conservatively unresolved.";
    let (gene_summaries, exon_summaries) = genes::summaries(&out.scope.groups);
    let mismatching_records = record_states
        .get("reference_mismatch")
        .copied()
        .unwrap_or(0);
    let mismatch_summary = json!({"all_input_records": mismatching_records,
        "records_overlapping_scope": out.diagnostics.reference_mismatch_records_overlapping_scope,
        "primary_bases_in_scope": out.counts[State::ReferenceMismatch as usize]});
    let data = json!({
        "schema": "gvcf-audit-report-v1", "version": env!("CARGO_PKG_VERSION"), "status": "complete",
        "inputs": {"gvcf": file_info(&args.gvcf)?, "reference": file_info(&args.reference)?, "reference_index": file_info(Path::new(&format!("{}.fai", args.reference.display())))?, "bed": args.bed.as_deref().map(file_info).transpose()?, "gene_targets": args.gene_targets.as_deref().map(file_info).transpose()?},
        "checksums_computed": false, "header": h, "options": args, "definition": definition,
        "policy": {"id": "small-variant-evidence-v1", "block_depth": "MIN_DP, with DP only when explicitly allowed", "reference_validation": "header lengths where supplied, record REF alleles and FASTA sequence layout; assembly identity is not inferred", "haploid_calls": "accepted without dosage conversion", "overlap": "two or more records unresolved, including split multiallelic records", "primary_reason_order": ["reference_ambiguous (base mask)", "overlapping_records", "reference_mismatch", "reference_ambiguous (REF allele)", "no_call", "partial_no_call", "unsupported_ploidy", "genotype_filtered", "filter_not_assessed or site_filtered", "quality_missing", "low_depth", "low_gq", "unsupported_allele", "complex_variant", "callable_reference or callable_variant", "no_record when absent"]},
        "scope": scope, "requested_bases": total, "callable_bases": called, "callable_fraction": called as f64 / total as f64,
        "bases_by_state": counts(&out.counts), "records_scanned": records, "record_states_before_overlap_and_reference_mask": record_states,
        "adaptations": adaptations, "warnings": warnings, "elapsed_seconds": seconds, "regions": regions,
        "contigs": contigs, "quality_checks": quality,
        "diagnostics": out.diagnostics, "reference_mismatches": mismatch_summary,
        "genes": gene_summaries, "exons": exon_summaries,
        "coordinates": "0-based half-open; overlapping input rows are independent, overall totals use their union"
    });
    std::fs::write(dir.join("report.json"), serde_json::to_vec_pretty(&data)?)?;
    let mut states_html = String::new();
    for state in STATES {
        let n = out.counts[state as usize];
        if n == 0 {
            continue;
        }
        let pct = 100.0 * n as f64 / total as f64;
        let class = if state.callable() { "good" } else { "other" };
        states_html.push_str(&format!("<tr><td>{}</td><td class=number>{n}</td><td class=number>{pct:.2}%</td><td><div class=track><div class='{class}' style='width:{pct:.4}%'></div></div></td></tr>", state.name().replace('_', " ")));
    }
    let warnings_html = warnings
        .iter()
        .map(|s| format!("<li>{}</li>", escape(s)))
        .collect::<String>();
    let region_rows: Vec<_> = out
        .scope
        .regions
        .iter()
        .map(|r| {
            json!([
                r.contig,
                r.start,
                r.end,
                r.name,
                r.bases_by_state[0] + r.bases_by_state[1]
            ])
        })
        .collect();
    let region_data = serde_json::to_string(&region_rows)?.replace('<', "\\u003c");
    let gene_data = serde_json::to_string(&gene_summaries)?.replace('<', "\\u003c");
    let exon_data = serde_json::to_string(&exon_summaries)?.replace('<', "\\u003c");
    let gene_html = if args.gene_targets.is_some() {
        format!(
            "{}{}",
            group_table(
                "gene",
                "Genes",
                "Gene, contig",
                "<th>Gene</th><th>Contig</th><th>Targeted exons</th><th>Unresolved exons</th>"
            ),
            group_table(
                "exon",
                "Exons",
                "Gene, exon, contig",
                "<th>Gene</th><th>Exon</th><th>Intervals</th>"
            )
        )
    } else {
        String::new()
    };
    let diagnostics_html = diagnostics_html(out, mismatching_records)?;
    let region_script = include_str!("report.js");
    let sample = escape(&h.sample);
    let caller = escape(&h.caller);
    let gvcf = escape(&args.gvcf.display().to_string());
    let version = env!("CARGO_PKG_VERSION");
    let quality_html = if !quality.requested {
        "No pass/fail limits were requested.".to_string()
    } else if quality.passed {
        "The requested quality limits passed.".to_string()
    } else {
        format!(
            "Quality limits failed: {}. The audit completed successfully; review the results below.",
            escape(&quality.failures.join("; "))
        )
    };
    let html = format!(
        r#"<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>gVCF audit · {sample}</title>
<style>
:root{{color-scheme:light;font-family:system-ui,-apple-system,sans-serif;color:#17313c;background:#f4f7f8}}*{{box-sizing:border-box}}body{{margin:0}}main{{max-width:1100px;margin:auto;padding:40px 24px 70px}}header{{border-top:5px solid #176c61;padding-top:24px}}.eyebrow{{text-transform:uppercase;letter-spacing:.12em;font-size:12px;color:#4d696c}}h1{{font-size:36px;margin:8px 0 12px;overflow-wrap:anywhere}}h2{{font-size:21px;margin:0 0 16px}}p{{line-height:1.6}}.muted{{color:#516b74;font-size:14px}}.cards{{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;margin:28px 0}}.card,section{{background:white;border:1px solid #dce6e8;border-radius:12px;padding:22px}}.big{{font-size:34px;font-weight:700;margin:8px 0}}.green{{color:#176c61}}td{{overflow-wrap:anywhere;max-width:380px}}details{{margin:12px 0}}summary{{cursor:pointer}}input,select,button{{font:inherit;padding:6px}}section{{margin-top:20px;overflow-x:auto}}table{{width:100%;border-collapse:collapse;font-size:14px}}th{{text-align:left;color:#526971;font-weight:600}}td,th{{padding:11px 8px;border-bottom:1px solid #e6edef}}td:first-child{{padding-left:0}}.number{{text-align:right;font-variant-numeric:tabular-nums}}.track{{width:130px;height:8px;border-radius:8px;background:#edf1f2;overflow:hidden}}.good{{height:100%;background:#208373}}.other{{height:100%;background:#d59b4c}}.notice{{border-left:4px solid #d59b4c}}li{{line-height:1.6;margin:8px 0}}a{{color:#156c62}}.path{{overflow-wrap:anywhere}}footer{{margin-top:25px;color:#516b74;font-size:13px}}@media(max-width:640px){{.cards{{grid-template-columns:1fr}}main{{padding:22px 14px}}h1{{font-size:28px}}.track{{width:65px}}}}
</style><main><header><div class=eyebrow>gvcf-audit · {version} · completed</div><h1>{sample}</h1><p>Where this sample has usable genotype evidence, and why other regions remain unresolved.</p><p class=muted>{scope} · {caller} rules · minimum depth {dp} · minimum GQ {gq}</p><p>{quality_html}</p></header>
<div class=cards><div class=card><div>Callable bases</div><div class="big green">{percent:.2}%</div><div class=muted>{called} of {total} requested bases</div></div><div class=card><div>Unresolved bases</div><div class=big>{unresolved}</div><div class=muted>Each assigned a primary reason</div></div><div class=card><div>Records read</div><div class=big>{records}</div><div class=muted>Full-file scan · {seconds:.2} seconds</div></div></div>
<section><h2>Evidence by status</h2><table><thead><tr><th>Status</th><th class=number>Bases</th><th class=number>Share</th><th>Distribution</th></tr></thead><tbody>{states_html}</tbody></table><p class=muted>Statuses are mutually exclusive. Reference ambiguity overrides record overlap, which overrides a record’s primary reason.</p></section>
<section class=notice><h2>How to interpret this report</h2><p>{definition}</p><ul>{warnings_html}</ul></section>
{diagnostics_html}{gene_html}
<section><h2>Requested regions</h2><p><label>Find a region <input id="region-search" type="search" placeholder="Name or contig"></label> <label>Sort <select id="region-sort"><option value="input">Input order</option><option value="lowest">Lowest callability</option><option value="highest">Highest callability</option></select></label></p><table><thead><tr><th>Name</th><th>Interval</th><th class=number>Bases</th><th class=number>Callable</th></tr></thead><tbody id="region-rows"></tbody></table><p class=muted><span id="region-count"></span> <button id="region-prev" type="button">Previous</button> <button id="region-next" type="button">Next</button></p><noscript><p>Enable JavaScript to explore regions, or open the All regions TSV download.</p></noscript><p class=muted>Coordinates are 0-based, half-open. Overlapping input rows are reported independently; overall totals count their union once.</p></section>
<section><h2>Downloads and provenance</h2><p><a href="callable.bed">Callable BED</a> · <a href="unresolved.bed">Unresolved BED</a> · <a href="regions.tsv">All regions TSV</a> · <a href="report.json">Full JSON report</a></p><p class="muted path">Input: {gvcf}</p><p class=muted>Input paths, sizes, modification times, caller adaptations and policy options are recorded in JSON. No checksums were computed. This report is self-contained and makes no network requests.</p></section><footer>Generated by gvcf-audit. Review the input caller’s limitations alongside these quality thresholds.</footer></main><script id="region-data" type="application/json">{region_data}</script><script id="gene-data" type="application/json">{gene_data}</script><script id="exon-data" type="application/json">{exon_data}</script><script>{region_script}</script></html>"#,
        dp = args.min_dp,
        gq = args.min_gq,
        unresolved = total - called
    );
    std::fs::write(dir.join("report.html"), html)?;
    Ok(())
}

fn group_table(id: &str, title: &str, placeholder: &str, headers: &str) -> String {
    format!(
        r#"<section><h2>{title}</h2><p><label>Find {id} <input id="{id}-search" type="search" placeholder="{placeholder}"></label> <label>Sort <select id="{id}-sort"><option value="lowest">Lowest callability</option><option value="highest">Highest callability</option><option value="input">Identifier order</option></select></label></p><table><thead><tr>{headers}<th>Bases</th><th>Unresolved</th><th>Callable</th><th>Main reason</th></tr></thead><tbody id="{id}-rows"></tbody></table><p class=muted><span id="{id}-count"></span> <button id="{id}-prev" type="button">Previous</button> <button id="{id}-next" type="button">Next</button></p><p class=muted>Only supplied exon targets are summarized; introns and other exons are outside this scope. Overlaps count once within each group. Different genes or exon identifiers may overlap; do not sum their totals. Exons with any unresolved base are counted as unresolved.</p><p><a href="{id}s.tsv">Download all {title} TSV</a></p><noscript>Enable JavaScript to explore this table, or use the TSV download.</noscript></section>"#
    )
}

fn diagnostics_html(out: &Output, all_records: u64) -> Result<String> {
    let d = &out.diagnostics;
    let mut text = format!(
        "<section><h2>Reference agreement</h2><p><strong>{all_records}</strong> mismatching records in the full input; <strong>{}</strong> overlap the requested scope; <strong>{}</strong> requested bases have reference mismatch as their primary status.</p><p class=muted>Ambiguous reference bases and overlapping records take precedence in base counts. The record limit includes those mismatches; the base limit counts only the primary status. A record is in scope when any part of its span overlaps a requested interval.</p></section>",
        d.reference_mismatch_records_overlapping_scope,
        out.counts[State::ReferenceMismatch as usize]
    );
    text.push_str(&format!("<section><h2>Examples to investigate</h2><p class=muted>Up to {} examples per reason and type, in input order. Coordinates are 0-based, half-open. Record examples overlap the scope but may extend beyond it. They show a primary record reason before overlap and reference masking, or missing MIN_DP; interval examples show the final base status. Examples are not an exhaustive list of problems.</p>", d.max_examples_per_reason));
    if d.record_examples.is_empty() && d.interval_examples.is_empty() {
        text.push_str("<p>No examples were retained. Check the counts above; examples can be disabled with --max-examples 0.</p>");
    }
    if !d.record_examples.is_empty() {
        text.push_str("<details><summary>Record details</summary><table><thead><tr><th>Reason</th><th>Record / input line</th><th>Observed values</th><th>Rule</th></tr></thead><tbody>");
        for r in &d.record_examples {
            let o = &r["observed"];
            let mut observed = [
                "GT",
                "REF",
                "ALT",
                "FILTER",
                "FT",
                "DP",
                "MIN_DP",
                "GQ",
                "depth_field_used",
            ]
            .iter()
            .map(|key| {
                let value = &o[*key];
                let label = if value.is_null() {
                    "missing".into()
                } else {
                    value
                        .as_str()
                        .map(str::to_string)
                        .unwrap_or_else(|| value.to_string())
                };
                format!("{key}={label}")
            })
            .collect::<Vec<_>>()
            .join("; ");
            let m = &o["first_reference_mismatch"];
            if !m.is_null() {
                observed.push_str(&format!(
                    "; first mismatch at {}: VCF={}, reference={}",
                    m["position"],
                    m["vcf_base"].as_str().unwrap(),
                    m["reference_base"].as_str().unwrap()
                ));
            }
            text.push_str(&format!(
                "<tr><td>{}</td><td>{}:{}–{}<br>line {}</td><td>{}</td><td>{}</td></tr>",
                escape(r["reason"].as_str().unwrap()),
                escape(r["contig"].as_str().unwrap()),
                r["start"],
                r["end"],
                r["line"],
                escape(&observed),
                escape(r["rule"].as_str().unwrap())
            ));
        }
        text.push_str("</tbody></table></details>");
    }
    if !d.interval_examples.is_empty() {
        text.push_str("<details><summary>Unresolved intervals</summary><table><thead><tr><th>Reason</th><th>Interval</th></tr></thead><tbody>");
        for e in &d.interval_examples {
            text.push_str(&format!(
                "<tr><td>{}</td><td>{}:{}–{}</td></tr>",
                e.reason,
                escape(&e.contig),
                e.start,
                e.end
            ));
        }
        text.push_str("</tbody></table></details>");
    }
    text.push_str("</section>");
    Ok(text)
}
