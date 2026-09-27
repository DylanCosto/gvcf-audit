mod common;

use common::Workspace;

#[test]
fn comparison_tracks_exact_gains_and_requires_policy_override() {
    let t = Workspace::core();
    t.audit("before", &[], 0);
    t.audit("after", &["--allow-unfiltered", "--allow-block-dp"], 0);
    t.run(
        &[
            "compare", "--before", "before", "--after", "after", "--out", "refused",
        ],
        1,
    );
    assert!(!t.0.join("refused").exists());
    t.run(
        &[
            "compare",
            "--before",
            "before",
            "--after",
            "after",
            "--out",
            "comparison",
            "--allow-policy-change",
        ],
        0,
    );
    let d = t.json("comparison/comparison.json");
    assert_eq!(d["summary"]["gained_callable_bases"], 7);
    assert_eq!(d["summary"]["lost_callable_bases"], 0);
    assert_eq!(
        t.read("comparison/changes.bed"),
        "chr1\t10\t11\tfilter_not_assessed->callable_variant\nchr1\t24\t30\tquality_missing->callable_reference\n"
    );
    t.run(
        &[
            "compare", "--before", "before", "--after", "before", "--out", "same",
        ],
        0,
    );
    assert!(t.read("same/changes.bed").is_empty());
}

#[test]
fn incomplete_bed_partition_is_rejected() {
    let t = Workspace::core();
    t.audit("before", &[], 0);
    t.audit("after", &[], 0);
    t.write("after/callable.bed", "");
    t.run(
        &[
            "compare", "--before", "before", "--after", "after", "--out", "bad",
        ],
        1,
    );
    assert!(!t.0.join("bad").exists());
}
