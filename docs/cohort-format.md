# Cohort summaries

```sh
gvcf-audit summarize sample-a/audit sample-b/audit --out cohort
gvcf-audit summarize --audit-list audits.txt --level exon --out exon-cohort
```

Each input directory needs a complete `report.json`, `callable.bed` and `unresolved.bed`,
generated using `--gene-targets`. The original gVCFs and FASTAs are not needed.
One input is allowed; multiple inputs must have distinct sample names.

`--audit-list` accepts one directory path per line, relative to the list file's directory.
Blank lines and lines beginning with `#` are ignored. Paths are not shell-expanded or
unquoted. Positional inputs come first, followed by listed inputs; that order determines
matrix columns and sample summary rows.

Inputs must have the same target union, selected gene/exon identifiers and merged intervals,
depth/GQ cutoffs, filter settings, caller policy and policy definition. Annotation row order
and duplicate rows may differ when merged targets agree. BED partitions and their global,
contig, gene and exon counts are checked against each JSON report.

Reference paths, sizes and modification times must match. For relocated files,
`--assume-same-reference` explicitly asserts identical reference sequences and records
warnings about differing metadata. Neither matching metadata nor this option verifies
sequence identity; input reports contain no sequence checksums. Different reblocking
status also produces a warning.

## Files

| File | Contents |
|---|---|
| `callable_percent.tsv` | Gene (or exon) rows and sample columns, percentages from 0 to 100 |
| `callable_bases.tsv` | The same matrix with exact integer callable-base counts |
| `groups.tsv` | Group size, minimum/mean/maximum callability, samples below the threshold and fully callable samples |
| `samples.tsv` | Overall target size, callable bases/percentage and quality-gate status |
| `cohort.json` | Versioned metadata, input provenance, policies, warnings and per-sample state counts |

Matrix keys are `contig`, `gene`, and optionally `exon`, followed by `requested_bases`.
Sample columns use the prefix `sample:` to avoid collisions with these keys. Groups are
sorted by contig and annotation IDs. Fields containing quotes are quoted with doubled
internal quotes. Identifiers containing control characters are rejected.

Percentages are rounded to six decimals; use the integer matrix for exact fractions.
Means weight samples equally because every sample has the same group targets. Overlaps
are merged within a group. Different groups can overlap and must not be summed.

`--min-callable-percent` defaults to 95. Samples strictly below this unrounded group
percentage are counted. This summarizes groups; it does **not** set a cohort exit-code gate.
Samples that failed their own audit quality gates remain included and are flagged in TSV/JSON.
Exit 0 means a valid summary was written. Invalid/incompatible input returns 1 and leaves
no completed output directory. Existing outputs are not overwritten. Keep inputs unchanged
while the command runs.

## JSON v1

`schema` is `gvcf-audit-cohort-v1`; `status` is `complete`. `version`, `level`,
`sample_count`, `group_count` and `min_callable_percent` describe the summary.
`policy` contains shared audit options, caller and policy definition. `compatibility`
records checks and reference assertions; `reference_identity_verified` is always false.

`samples` follows input order. Each row records the sample, canonical audit directory,
generating version, requested/callable bases and percentage, quality checks, original
warnings, reference metadata, reblocking status and base counts by state.
Top-level `warnings` describe cohort limitations; `files` names the TSV exports.
Consumers should allow additional fields within v1.

## Scale and interpretation

Reports and BEDs are validated one sample at a time. Memory holds one input report,
target intervals, sample metadata and an integer matrix for the selected groups.
The values for 20,000 genes and 500 samples take about 80 MB, plus report and target
memory. Exon matrices can be much larger, so genes and exons use separate runs.

Tests cover 500 small saved audits and independent randomized gene/exon comparisons.
This checks accounting and column handling, not performance on 500 production genomes.
Callability describes usable evidence, not variant accuracy. Matching rules do not
establish comparable assays or caller calibration. Use the tables to investigate uneven
missingness before downstream analysis, not to infer variant absence in unresolved regions.
