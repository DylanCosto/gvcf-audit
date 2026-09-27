mod common;

use common::{Workspace, fixture};
use flate2::{Compression, write::GzEncoder};
use serde_json::{Value, json};
use std::{fs, io::Write};

#[test]
fn audit_matches_hand_calculated_intervals_and_counts() {
    let t = Workspace::core();
    t.audit("audit", &[], 0);
    for name in ["callable.bed", "unresolved.bed"] {
        assert_eq!(
            t.read(&format!("audit/{name}")),
            fs::read_to_string(fixture(name)).unwrap()
        );
    }
    let actual = t.json("audit/report.json");
    let expected: Value =
        serde_json::from_str(&fs::read_to_string(fixture("summary.json")).unwrap()).unwrap();
    for (key, value) in expected.as_object().unwrap() {
        assert_eq!(&actual[key], value, "{key}");
    }
    assert_eq!(actual["header"]["caller"], "gatk");
    assert_eq!(
        actual["reference_mismatches"],
        json!({
            "all_input_records": 1, "records_overlapping_scope": 1, "primary_bases_in_scope": 1
        })
    );
}

#[test]
fn raw_gatk_filter_is_visible_and_override_is_explicit() {
    let t = Workspace::core();
    let result = t.audit("default", &[], 0);
    assert!(String::from_utf8_lossy(&result.stderr).contains("filter_not_assessed"));
    assert!(
        t.read("default/report.html")
            .contains("Raw HaplotypeCaller")
    );
    assert!(
        t.json("default/report.json")["warnings"]
            .as_array()
            .unwrap()
            .iter()
            .any(|v| v.as_str().unwrap().contains("--allow-unfiltered"))
    );
    let quiet = t.audit("quiet", &["--quiet"], 0);
    assert!(quiet.stderr.is_empty());
    assert!(t.read("quiet/report.html").contains("Raw HaplotypeCaller"));
    t.audit("allowed", &["--allow-unfiltered"], 0);
    let allowed = t.json("allowed/report.json");
    assert_eq!(allowed["callable_bases"], 11);
    assert_eq!(allowed["bases_by_state"]["filter_not_assessed"], 0);
    assert_eq!(allowed["bases_by_state"]["low_depth"], 5);
    assert_eq!(allowed["bases_by_state"]["quality_missing"], 6);
}

#[test]
fn gene_and_exon_totals_count_overlaps_once() {
    let t = Workspace::core();
    t.audit("audit", &["--gene-targets", "targets.tsv"], 0);
    let d = t.json("audit/report.json");
    assert_eq!(d["requested_bases"], 30);
    assert_eq!(d["callable_bases"], 10);
    let genes: Vec<_> = d["genes"]
        .as_array()
        .unwrap()
        .iter()
        .map(|g| json!([g["gene"], g["bases"], g["callable_bases"]]))
        .collect();
    assert_eq!(genes, vec![json!(["G", 22, 10]), json!(["H", 15, 1])]);
    let exons: Vec<_> = d["exons"]
        .as_array()
        .unwrap()
        .iter()
        .map(|e| json!([e["exon"], e["bases"], e["callable_bases"]]))
        .collect();
    assert_eq!(
        exons,
        vec![
            json!(["E1", 12, 6]),
            json!(["E2", 12, 5]),
            json!(["E3", 15, 1])
        ]
    );
}

#[test]
fn quality_gate_keeps_the_report_but_input_errors_do_not() {
    let t = Workspace::core();
    t.audit("pass", &["--min-callable-percent", "20"], 0);
    t.audit("fail", &["--min-callable-percent", "20.001"], 2);
    assert_eq!(
        t.json("fail/report.json")["quality_checks"]["passed"],
        false
    );
    t.write(
        "sample.g.vcf",
        &t.read("sample.g.vcf").replace("END=5", "END=0"),
    );
    t.audit("bad", &[], 1);
    assert!(!t.0.join("bad").exists());
    let before = t.read("pass/report.json");
    t.audit("pass", &[], 1);
    assert_eq!(t.read("pass/report.json"), before);
}

#[test]
fn gzip_matches_plain_and_truncation_leaves_no_report() {
    let t = Workspace::core();
    t.audit("plain", &[], 0);
    let mut gzip = GzEncoder::new(Vec::new(), Compression::default());
    gzip.write_all(t.read("sample.g.vcf").as_bytes()).unwrap();
    let bytes = gzip.finish().unwrap();
    fs::write(t.0.join("sample.g.vcf"), &bytes).unwrap();
    t.audit("compressed", &[], 0);
    for name in ["callable.bed", "unresolved.bed", "regions.tsv"] {
        assert_eq!(
            t.read(&format!("plain/{name}")),
            t.read(&format!("compressed/{name}"))
        );
    }
    fs::write(t.0.join("sample.g.vcf"), &bytes[..bytes.len() - 5]).unwrap();
    t.audit("truncated", &[], 1);
    assert!(!t.0.join("truncated").exists());
}

#[test]
fn stale_fasta_index_is_rejected() {
    let t = Workspace::core();
    t.write(
        "reference.fa.fai",
        &t.read("reference.fa.fai").replace("\t6\t", "\t7\t"),
    );
    t.audit("bad", &[], 1);
    assert!(!t.0.join("bad").exists());
}
