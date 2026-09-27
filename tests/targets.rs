mod common;

use common::Workspace;
use std::{fs, path::Path, process::Command};

fn examples() -> Workspace {
    let t = Workspace::new();
    for name in ["genes.gtf", "genes.gff3"] {
        t.copy(
            &Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("examples")
                .join(name),
            name,
        );
    }
    t
}

#[test]
fn both_annotation_formats_match_saved_targets() {
    let t = examples();
    let expected = include_str!("fixtures/annotation-targets.tsv");
    for format in ["gtf", "gff3"] {
        let out = format!("{format}.tsv");
        t.run(
            &[
                "targets",
                "--annotation",
                &format!("genes.{format}"),
                "--out",
                &out,
            ],
            0,
        );
        assert_eq!(t.read(&out), expected);
    }
}

#[test]
fn forward_parents_shared_exons_and_negative_single_base_coordinates() {
    let t = Workspace::new();
    t.write("input.gff3", "##gff-version 3\nchr1\t.\texon\t7\t7\t.\t-\t.\tParent=T,H\nchr1\t.\tmRNA\t1\t10\t.\t-\t.\tID=T;Parent=G\nchr1\t.\tgene\t1\t10\t.\t-\t.\tID=G\nchr1\t.\tgene\t1\t10\t.\t-\t.\tID=H\n");
    t.run(
        &[
            "targets",
            "--annotation",
            "input.gff3",
            "--out",
            "targets.tsv",
        ],
        0,
    );
    assert_eq!(
        t.read("targets.tsv"),
        "contig\tstart\tend\tgene\texon\nchr1\t6\t7\tG\tgvcf-audit:exon:6-7:-\nchr1\t6\t7\tH\tgvcf-audit:exon:6-7:-\n"
    );
}

#[test]
fn missing_parent_and_missing_gene_filter_leave_no_output() {
    let t = examples();
    t.write(
        "broken.gff3",
        "chr1\t.\texon\t1\t10\t.\t+\t.\tParent=missing\n",
    );
    t.run(
        &[
            "targets",
            "--annotation",
            "broken.gff3",
            "--out",
            "broken.tsv",
        ],
        1,
    );
    t.run(
        &[
            "targets",
            "--annotation",
            "genes.gtf",
            "--gene",
            "missing",
            "--out",
            "missing.tsv",
        ],
        1,
    );
    assert!(!t.0.join("broken.tsv").exists());
    assert!(!t.0.join("missing.tsv").exists());
}

#[test]
fn concurrent_converters_never_replace_completed_output() {
    let t = examples();
    let mut command = Command::new(env!("CARGO_BIN_EXE_gvcf-audit"));
    command
        .current_dir(&t.0)
        .args([
            "targets",
            "--annotation",
            "genes.gtf",
            "--out",
            "targets.tsv",
            "--quiet",
        ])
        .stderr(std::process::Stdio::null());
    let mut a = command.spawn().unwrap();
    let mut b = command.spawn().unwrap();
    let mut codes = [a.wait().unwrap().code(), b.wait().unwrap().code()];
    codes.sort();
    assert_eq!(codes, [Some(0), Some(1)]);
    assert_eq!(
        t.read("targets.tsv"),
        include_str!("fixtures/annotation-targets.tsv")
    );
    assert_eq!(fs::read_dir(&t.0).unwrap().count(), 3);
}
