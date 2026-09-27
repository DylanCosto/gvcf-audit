mod common;

use common::Workspace;
use serde_json::json;
use std::fs;

fn pair() -> Workspace {
    let t = Workspace::core();
    t.audit("a", &["--gene-targets", "targets.tsv"], 0);
    t.write(
        "sample.g.vcf",
        &t.read("sample.g.vcf")
            .replace("\tS\n", "\tT\n")
            .replace("0/0:2:30", "0/0:20:30"),
    );
    t.audit(
        "b",
        &[
            "--gene-targets",
            "targets.tsv",
            "--min-callable-percent",
            "90",
        ],
        2,
    );
    t
}

#[test]
fn cohort_matrices_match_hand_calculated_genes_and_exons() {
    let t = pair();
    // Saved audits remain usable after their original gVCF and FASTA are removed.
    fs::remove_file(t.0.join("sample.g.vcf")).unwrap();
    fs::remove_file(t.0.join("reference.fa")).unwrap();
    t.run(
        &[
            "summarize",
            "a",
            "b",
            "--min-callable-percent",
            "50",
            "--out",
            "cohort",
        ],
        0,
    );
    assert_eq!(
        t.read("cohort/callable_bases.tsv"),
        "contig\tgene\trequested_bases\tsample:S\tsample:T\nchr1\tG\t22\t10\t15\nchr1\tH\t15\t1\t1\n"
    );
    assert_eq!(
        t.read("cohort/callable_percent.tsv"),
        "contig\tgene\trequested_bases\tsample:S\tsample:T\nchr1\tG\t22\t45.454545\t68.181818\nchr1\tH\t15\t6.666667\t6.666667\n"
    );
    assert!(
        t.read("cohort/groups.tsv")
            .contains("chr1\tG\t22\t45.454545\t56.818182\t68.181818\t1\t0\n")
    );
    assert!(
        t.read("cohort/samples.tsv")
            .contains("T\t30\t15\t50.000000\tfailed\n")
    );
    assert_eq!(t.json("cohort/cohort.json")["sample_count"], 2);
    fs::create_dir(t.0.join("lists")).unwrap();
    t.write("lists/audits.txt", "# Saved audits\n../b\n../a\n");
    t.run(
        &[
            "summarize",
            "--audit-list",
            "lists/audits.txt",
            "--level",
            "exon",
            "--out",
            "exons",
        ],
        0,
    );
    assert_eq!(
        t.read("exons/callable_bases.tsv"),
        "contig\tgene\texon\trequested_bases\tsample:T\tsample:S\nchr1\tG\tE1\t12\t11\t6\nchr1\tG\tE2\t12\t5\t5\nchr1\tH\tE3\t15\t1\t1\n"
    );
}

#[test]
fn incompatible_or_damaged_audits_leave_no_cohort() {
    let t = pair();
    let original = t.json("b/report.json");
    let mut variants = Vec::new();
    let mut duplicate = original.clone();
    duplicate["header"]["sample"] = json!("S");
    variants.push(duplicate);
    let mut policy = original.clone();
    policy["options"]["min_dp"] = json!(11);
    variants.push(policy);
    let mut count = original.clone();
    count["genes"][0]["callable_bases"] = json!(14);
    variants.push(count);
    let mut interval = original.clone();
    interval["genes"][0]["intervals"][0] = json!([1, 23]);
    variants.push(interval);
    let mut missing = original.clone();
    missing["genes"] = json!([]);
    variants.push(missing);
    let mut missing_gate = original.clone();
    missing_gate["quality_checks"]["requested"] = json!(null);
    variants.push(missing_gate);
    for (i, report) in variants.iter().enumerate() {
        t.write("b/report.json", &report.to_string());
        let out = format!("bad-{i}");
        t.run(&["summarize", "a", "b", "--out", &out], 1);
        assert!(!t.0.join(out).exists());
    }
    t.write("b/report.json", &original.to_string());
    t.write("b/callable.bed", "");
    t.run(&["summarize", "a", "b", "--out", "damaged"], 1);
    assert!(!t.0.join("damaged").exists());
    t.run(&["summarize", "a", "a", "--out", "duplicate"], 1);
    assert!(!t.0.join("duplicate").exists());
}

#[test]
fn reference_override_is_explicit_and_recorded() {
    let t = pair();
    let mut report = t.json("b/report.json");
    report["inputs"]["reference"]["path"] = json!("/relocated/reference.fa");
    t.write("b/report.json", &report.to_string());
    t.run(&["summarize", "a", "b", "--out", "refused"], 1);
    assert!(!t.0.join("refused").exists());
    t.run(
        &[
            "summarize",
            "a",
            "b",
            "--assume-same-reference",
            "--out",
            "accepted",
        ],
        0,
    );
    let result = t.json("accepted/cohort.json");
    assert_eq!(result["compatibility"]["assume_same_reference"], true);
    assert_eq!(
        result["compatibility"]["reference_identity_verified"],
        false
    );
    assert!(
        result["warnings"]
            .as_array()
            .unwrap()
            .iter()
            .any(|w| w.as_str().unwrap().contains("T: reference metadata differ"))
    );
}

#[test]
fn five_hundred_saved_samples_keep_distinct_matrix_columns() {
    let t = Workspace::core();
    t.audit("original", &["--gene-targets", "targets.tsv"], 0);
    let mut report = t.json("original/report.json");
    let mut list = String::new();
    for i in 0..500 {
        let name = format!("S{i:03}");
        fs::create_dir(t.0.join(&name)).unwrap();
        for file in ["callable.bed", "unresolved.bed"] {
            fs::copy(t.0.join("original").join(file), t.0.join(&name).join(file)).unwrap();
        }
        report["header"]["sample"] = json!(name);
        t.write(&format!("{name}/report.json"), &report.to_string());
        list.push_str(&format!("{name}\n"));
    }
    t.write("audits.txt", &list);
    t.run(
        &["summarize", "--audit-list", "audits.txt", "--out", "cohort"],
        0,
    );
    let table = t.read("cohort/callable_bases.tsv");
    let rows: Vec<Vec<_>> = table
        .lines()
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 3);
    assert_eq!(rows[0].len(), 503);
    assert_eq!(rows[0][3], "sample:S000");
    assert_eq!(rows[0][502], "sample:S499");
    assert!(rows[1][3..].iter().all(|&n| n == "10"));
    assert!(rows[2][3..].iter().all(|&n| n == "1"));
    assert_eq!(t.json("cohort/cohort.json")["sample_count"], 500);
}
