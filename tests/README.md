# Running the checks

From a fresh checkout, with Rust 1.88 or later:

```sh
cargo test --locked
cargo clippy --locked --all-targets -- -D warnings
cargo fmt --check
```

The Rust integration tests run the actual command-line binary. They use small,
synthetic files in `fixtures/` and the public `examples/`, with no downloaded data,
Python, caller installations or old release binaries required.

The independent Python checks use only Python 3.10+ and its standard library:

```sh
cargo build --locked
python3 tests/oracles/run.py target/debug/gvcf-audit
```

Each suite creates and cleans up its own temporary directory. The command works
from any working directory when given absolute paths. Do not use Python `-O`,
which disables assertions; the runner refuses it. Failures return a nonzero exit
code. Both layers run in CI on stable Rust and Rust 1.88, including a release build.
Push CI runs only on `main`; pull requests run their own checks.

A separate weekly job runs the synthetic suites with a fresh seed each Monday at
04:17 UTC. It prints the seed and replay command before starting, and repeats the
command on failure. To try new cases locally or reproduce a weekly failure:

```sh
python3 tests/oracles/run.py target/debug/gvcf-audit --random-seed
python3 tests/oracles/run.py target/debug/gvcf-audit --seed 12345
```

The override applies to each randomized suite. Omitting both flags retains the
fixed suite defaults. No additional Python dependencies are needed.

## What is covered

- Audit state precedence, gaps, overlaps, reference ambiguity/mismatches, contig
  aliases, samples, caller adaptations, GQ/depth and filter rules.
- Reference index validation, compressed/truncated input, HTML escaping, quality
  gates, errors that must leave no report, and existing output protection.
- Lowercase reference bases (including ambiguous `n`) and a deletion followed by
  a `*` spanning-deletion record, with the overlap remaining unresolved.
- Gene/exon interval unions and bounded diagnostic examples.
- Exact comparison gains/losses, equal totals at different positions, metadata
  guards, malformed partitions, group totals and changed-interval exports.
- GTF/GFF3 coordinates, attributes, shared exons, parent graphs, selectors,
  missing IDs, compressed input and concurrent output protection.

`oracles/` preserves the latest synthetic checks used during development, with
workspace-specific paths and dependencies on old binaries removed. Fixed seeds
make the randomized cases repeatable: audit 892, genes 71942, comparison 40521,
annotation conversion 5021. These use per-base sets/counts and independently
constructed interval models, rather than the Rust implementation's interval sweep.
Comparison oracles check the comparison against the two input audit partitions;
the audit and gene oracles separately check classification and grouping.

The large public-genome runs documented in `docs/validation.md` remain historical
validation evidence. They are not rerun by this suite and are not prerequisites.
Passing these checks is not proof of biological accuracy or every caller version.

## Saved expected outputs

`fixtures/core/` uses a 40-base chr1 and an unobserved 10-base chr2. Its BEDs and
summary were calculated by hand. The default has 10 callable bases: nine reference
bases and one passing variant. One raw GATK variant has an unset filter, six block
bases lack MIN_DP, and overlapping records remain unresolved. Explicitly allowing
unset filters and block DP adds exactly seven callable bases.

The overlapping gene targets cover 30 bases in total. Gene G covers 22 bases with
10 callable; H covers 15 with one callable. Exons E1/E2/E3 cover 12/12/15 bases
with 6/5/1 callable. Repeating E1 must not inflate these totals.

`fixtures/annotation-targets.tsv` contains the expected conversion of the example
GTF/GFF3: subtract one from each start, retain each end and preserve gene/exon IDs.
It was independently checked before being saved here.

When behavior intentionally changes, review the expected intervals and policy
first. Do not refresh expected outputs just to make a failing test pass.
# MultiQC integration

`python3 tests/integrations/multiqc.py target/debug/gvcf-audit /path/to/multiqc`
runs four small audits through MultiQC in strict mode and checks its parsed General Statistics.
It covers numeric and quoted sample names, a leading `#`, failed/passed/unrequested gates,
and the explicit unset-filter override. CI uses MultiQC 1.35; this dependency is only for the test.

## Cohort checks

The Rust tests compare hand-calculated gene/exon matrices, check incompatible and damaged
inputs, and exercise 500 distinct saved samples. The source gVCF and FASTA are removed in
one test to check that summaries are self-contained. `oracles/cohort.py` independently
computes per-base callability and overlapping target unions for 12 randomized samples,
then checks both matrices and group statistics. It uses the same replayable seed mechanism
as the other suites.

The public-genome GIAB example has five additional small interval checks in
`examples/giab/test_overlap.py`. CI runs them without downloading genomes.
