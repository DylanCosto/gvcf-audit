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
fn multiqc_summary_preserves_counts_names_and_gate_status() {
    let t = Workspace::core();
    let original = t.read("sample.g.vcf");
    for (i, sample) in ["001", "#sample", "sample \"quoted\""].iter().enumerate() {
        t.write(
            "sample.g.vcf",
            &original.replace("\tS\n", &format!("\t{sample}\n")),
        );
        let out = format!("audit-{i}");
        t.audit(&out, &[], 0);
        let text = t.read(&format!("{out}/gvcf_audit_mqc.tsv"));
        assert!(text.contains("# plot_type: generalstats\n"));
        assert_eq!(
            text.lines().last().unwrap(),
            format!(
                "\"{sample}\"\t20.000000\t50\t46.000000\t2.000000\t\"not requested\"\t10\t20\t\"gatk\"\tfalse\tfalse"
            )
        );
    }
    t.audit("pass", &["--min-callable-percent", "20"], 0);
    t.audit("fail", &["--min-callable-percent", "21", "--quiet"], 2);
    assert!(t.read("pass/gvcf_audit_mqc.tsv").contains("\t\"passed\"\t"));
    assert!(t.read("fail/gvcf_audit_mqc.tsv").contains("\t\"failed\"\t"));
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

#[test]
fn soft_masked_reference_bases_are_callable() {
    let t = Workspace::new();
    t.write("reference.fa", ">chr1\nacgTacgtn\n");
    t.write("reference.fa.fai", "chr1\t9\t6\t9\t10\n");
    t.write(
        "sample.g.vcf",
        concat!(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=9>\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
            "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=1\tGT:MIN_DP:GQ\t0/0:20:30\n",
            "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DP:GQ\t0/1:20:30\n",
            "chr1\t3\t.\tG\t<NON_REF>\t.\tPASS\tEND=9\tGT:MIN_DP:GQ\t0/0:20:30\n",
        ),
    );
    t.audit("audit", &[], 0);
    assert_eq!(
        t.read("audit/callable.bed"),
        "chr1\t0\t1\tcallable_reference\nchr1\t1\t2\tcallable_variant\nchr1\t2\t8\tcallable_reference\n"
    );
    assert_eq!(
        t.read("audit/unresolved.bed"),
        "chr1\t8\t9\treference_ambiguous\n"
    );
    let report = t.json("audit/report.json");
    assert_eq!(report["callable_bases"], 8);
    assert_eq!(report["reference_mismatches"]["all_input_records"], 0);
}

#[test]
fn deletion_and_spanning_deletion_stay_unresolved() {
    let t = Workspace::new();
    t.write("reference.fa", ">chr1\nACGTA\n");
    t.write("reference.fa.fai", "chr1\t5\t6\t5\t6\n");
    t.write(
        "sample.g.vcf",
        concat!(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=5>\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
            "chr1\t1\t.\tACG\tA\t.\tPASS\t.\tGT:DP:GQ\t0/1:20:30\n",
            "chr1\t2\t.\tC\t*\t.\tPASS\t.\tGT:DP:GQ\t0/1:20:30\n",
            "chr1\t4\t.\tT\t<NON_REF>\t.\tPASS\tEND=5\tGT:MIN_DP:GQ\t0/0:20:30\n",
        ),
    );
    t.audit("audit", &[], 0);
    assert_eq!(
        t.read("audit/callable.bed"),
        "chr1\t3\t5\tcallable_reference\n"
    );
    assert_eq!(
        t.read("audit/unresolved.bed"),
        "chr1\t0\t1\tcomplex_variant\nchr1\t1\t2\toverlapping_records\nchr1\t2\t3\tcomplex_variant\n"
    );
    let report = t.json("audit/report.json");
    assert_eq!(report["callable_bases"], 2);
    assert_eq!(
        report["record_states_before_overlap_and_reference_mask"]["unsupported_allele"],
        1
    );
}
