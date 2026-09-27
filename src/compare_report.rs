use crate::{Result, STATES, report::escape};
use serde_json::Value;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

pub fn write(dir: &Path, data: &Value) -> Result<()> {
    let mut json = BufWriter::new(File::create(dir.join("comparison.json"))?);
    serde_json::to_writer_pretty(&mut json, data)?;
    json.flush()?;
    for kind in ["genes", "exons"] {
        let mut f = BufWriter::new(File::create(dir.join(format!("{kind}.tsv")))?);
        write!(f, "contig\tgene")?;
        if kind == "exons" {
            write!(f, "\texon")?;
        }
        writeln!(
            f,
            "\tbases\tbefore_callable_bases\tafter_callable_bases\tgained_callable_bases\tlost_callable_bases\tnet_callable_bases\tbefore_callable_percent\tafter_callable_percent\tcallable_percentage_point_change"
        )?;
        for r in data[kind].as_array().unwrap() {
            write!(
                f,
                "{}\t{}",
                r["contig"].as_str().unwrap(),
                r["gene"].as_str().unwrap()
            )?;
            if kind == "exons" {
                write!(f, "\t{}", r["exon"].as_str().unwrap())?;
            }
            for key in [
                "bases",
                "before_callable_bases",
                "after_callable_bases",
                "gained_callable_bases",
                "lost_callable_bases",
                "net_callable_bases",
            ] {
                write!(f, "\t{}", r[key])?;
            }
            writeln!(
                f,
                "\t{:.4}\t{:.4}\t{:.4}",
                100.0 * r["before_callable_fraction"].as_f64().unwrap(),
                100.0 * r["after_callable_fraction"].as_f64().unwrap(),
                r["callable_percentage_point_change"].as_f64().unwrap()
            )?;
        }
        f.flush()?;
    }
    let warnings = data["warnings"]
        .as_array()
        .unwrap()
        .iter()
        .map(|w| format!("<li>{}</li>", escape(w.as_str().unwrap())))
        .collect::<String>();
    let s = &data["summary"];
    let mut states = String::new();
    for state in STATES {
        let key = state.name();
        states.push_str(&format!(
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
            key.replace('_', " "),
            s["before_bases_by_state"][key],
            s["after_bases_by_state"][key],
            s["bases_by_state_delta"][key]
        ));
    }
    let mut tables = String::new();
    for (kind, title) in [("genes", "Genes"), ("exons", "Exons")] {
        tables.push_str(&format!(r#"<section><h2>{title}</h2><label>Search <input id="{kind}-search" type="search" placeholder="Gene, exon or contig"></label> <label>Sort <select id="{kind}-sort"><option value="lost_callable_bases">Most bases lost</option><option value="gained_callable_bases">Most bases gained</option><option value="net_callable_bases">Largest net gain</option></select></label><table><thead><tr><th>Gene / exon</th><th>Contig</th><th>Before</th><th>After</th><th>Gained</th><th>Lost</th><th>Net change</th></tr></thead><tbody id="{kind}-rows"></tbody></table><p><span id="{kind}-count"></span> <button id="{kind}-prev">Previous</button> <button id="{kind}-next">Next</button></p><a href="{kind}.tsv">Download {kind} TSV</a></section>"#));
    }
    let mut table_data = serde_json::Map::new();
    for kind in ["genes", "exons"] {
        let rows: Vec<_> = data[kind]
            .as_array()
            .unwrap()
            .iter()
            .map(|r| {
                let mut row = serde_json::Map::new();
                for key in [
                    "gene",
                    "exon",
                    "contig",
                    "before_callable_fraction",
                    "after_callable_fraction",
                    "gained_callable_bases",
                    "lost_callable_bases",
                    "net_callable_bases",
                ] {
                    row.insert(key.into(), r[key].clone());
                }
                Value::Object(row)
            })
            .collect();
        table_data.insert(kind.into(), Value::Array(rows));
    }
    let rows = serde_json::to_string(&table_data)?.replace('<', "\\u003c");
    let html = format!(
        r#"<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>gvcf-audit comparison</title><style>
body{{font:16px/1.55 system-ui,sans-serif;color:#172b38;background:#f4f7f8;margin:0}}main{{max-width:1120px;margin:auto;padding:36px 24px}}h1{{font-size:34px;margin-bottom:4px}}h2{{margin-top:0}}section{{background:white;border:1px solid #d7e1e5;border-radius:10px;padding:22px;margin:22px 0;overflow:auto}}table{{width:100%;border-collapse:collapse;font-variant-numeric:tabular-nums}}th,td{{padding:9px;text-align:right;border-bottom:1px solid #e4e9ec}}th:first-child,td:first-child{{text-align:left;overflow-wrap:anywhere}}input,select,button{{font:inherit;padding:7px;border:1px solid #bccbd2;border-radius:5px}}input,select{{margin:0 12px 14px 0}}button{{cursor:pointer;background:#fff}}button:disabled{{opacity:.4;cursor:default}}a{{color:#146b83}}.muted{{color:#536975}}.cards{{display:flex;gap:24px;flex-wrap:wrap}}.cards strong{{font-size:28px;display:block}}.notice{{background:#fff9e9;border-color:#e4d59e}}.path{{overflow-wrap:anywhere}}</style>
<main><p class="muted">gvcf-audit · comparison</p><h1>{before} → {after}</h1><p>Changes in reported genotype evidence on identical target intervals.</p><section><div class="cards"><div><strong>{gained}</strong>bases gained callability</div><div><strong>{lost}</strong>bases lost callability</div><div><strong>{net}</strong>net change</div></div><p>Before: {old:.2}% callable · After: {new:.2}% callable · {bases} requested bases</p><p class="muted">Gained and lost count actual positions, even when the net change is zero. Changes between two callable states do not count as gains or losses. Gene/exon targets can overlap, so their totals must not be added together.</p></section>
<section class="notice"><h2>Comparison conditions</h2><ul>{warnings}</ul><p>These changes describe the caller's reported evidence, not a change in biological accuracy.</p></section>{tables}
<section><h2>Reasons before and after</h2><table><thead><tr><th>State</th><th>Before bases</th><th>After bases</th><th>Net change</th></tr></thead><tbody>{states}</tbody></table><p>Reason totals can cancel across positions. The changed intervals and transition counts preserve those differences.</p></section>
<section><h2>Files and provenance</h2><p><a href="changes.bed">Changed intervals</a> · <a href="comparison.json">Full comparison JSON</a></p><p class="path">Before: {before_path}</p><p class="path">After: {after_path}</p><p class="muted">BED coordinates are 0-based, half-open. No gVCF or FASTA was rescanned. Keep these files together. This report works offline.</p><noscript>Enable JavaScript to explore genes and exons, or use the TSV downloads.</noscript></section></main><script id="comparison-data" type="application/json">{rows}</script><script>{script}</script></html>"#,
        before = escape(data["before"]["header"]["sample"].as_str().unwrap()),
        after = escape(data["after"]["header"]["sample"].as_str().unwrap()),
        gained = s["gained_callable_bases"],
        lost = s["lost_callable_bases"],
        net = s["net_callable_bases"],
        bases = s["bases"],
        old = 100.0 * s["before_callable_fraction"].as_f64().unwrap(),
        new = 100.0 * s["after_callable_fraction"].as_f64().unwrap(),
        before_path = escape(data["before"]["directory"].as_str().unwrap()),
        after_path = escape(data["after"]["directory"].as_str().unwrap()),
        script = include_str!("compare.js")
    );
    std::fs::write(dir.join("report.html"), html)?;
    Ok(())
}
