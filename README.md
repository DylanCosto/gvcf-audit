# gvcf-audit

When a gVCF reports no variant, does that mean the sample matches the reference, or that there
wasn't enough evidence to tell?

**gvcf-audit checks where a sample's gVCF has usable genotype evidence and explains the gaps.**
It applies explicit depth, genotype-quality and filter rules across the whole reference or your chosen
regions, with summaries for genes and exons when you supply their coordinates.

For example, an exon with 200 bases might have 180 that meet your requirements, 15 with low depth,
and 5 with no genotype evidence. The report identifies those 20 uncertain bases and their coordinates.
It does not tell you whether they contain a missed variant. Even 100% callable does not prove every
variant was detected correctly. The reference FASTA is used to check reference agreement; it is not
an expected list of variants for the sample.

[**Explore the example reports**](https://DylanCosto.github.io/gvcf-audit/) without installing anything.
They use the small, fictional dataset included in this repository.

**Version 0.5.0.** Runs locally, with offline HTML reports and BED, TSV and JSON exports.
It does not call variants, inspect BAM reads, impute missing genotypes or provide clinical interpretation.
Use the reported coordinates to investigate the underlying reads in a genome viewer such as IGV.

## Install and run

Download the latest published Linux binary or source archive from the
[releases page](https://github.com/DylanCosto/gvcf-audit/releases).
A source checkout or unpacked source archive can be installed with:

```sh
cargo install --locked --path .
gvcf-audit --version
```

Alternatively, make the supplied Linux x86-64 binary executable and run it directly. It needs no Python,
Java, or caller installation. The supplied binary requires glibc 2.34 or newer. On older Linux systems,
build from source. Other platforms can build from source; only Linux has been validated here.

Build with Rust 1.88 or later:

```sh
cargo build --release --locked
./target/release/gvcf-audit \
    --gvcf sample.g.vcf.gz \
    --reference reference.fa \
    --bed targets.bed \
    --out audit
```

Try the included example without downloading anything:

```sh
./target/release/gvcf-audit --gvcf examples/sample.g.vcf \
    --reference examples/reference.fa --bed examples/targets.bed --out example-audit
```

The example has 600 callable bases out of 1,000. It includes low depth, missing calls, an ambiguous
reference stretch, a variant and a gap.

The output directory must be new. Open `audit/report.html` in a browser; it works offline.

- `--gvcf`: text VCF, gzip or BGZF. Multi-sample files require `--sample NAME`.
- `--reference`: the same uncompressed FASTA used to call the sample, with its `.fai` index.
  Create an index with `samtools faidx reference.fa` if needed.
- `--bed`: optional tab-separated BED3 or BED4, plain or gzip. Coordinates are 0-based, half-open.
  Without either target option, the scope is **every FASTA contig**, including unobserved decoys and alternate contigs.
- `--min-dp 10 --min-gq 20`: defaults; change them to suit the assay.

Overlapping BED rows get independent summaries. Overall totals count their union once. BED12 blocks are
not interpreted: split exons into separate BED intervals first. Extra BED columns are ignored.

The tool makes one sequential scan of the complete VCF, even when a BED is supplied.
No VCF index or checksum pass is required. It uses interval boundaries rather than expanding gVCF blocks
into per-base records. The FASTA layout is checked against its index once, then one contig at a time is memory-mapped
and scanned for ambiguous reference runs. Memory also holds BED intervals, reference ambiguity runs,
and simultaneously overlapping records. Report memory and file size also grow with the number of
input rows and gene/exon groups. Identical annotation rows can be deduplicated before an audit if
separate transcript-row summaries are not needed.

## Outputs

| File | Contents |
|---|---|
| `report.html` | Offline overview, primary reasons and searchable, sortable requested regions |
| `report.json` | Full region results, policy, caller adaptations, input paths/sizes/mtimes and warnings |
| `regions.tsv` | Every BED row, its callable fraction and base counts for each status |
| `callable.bed` | Merged intervals labelled `callable_reference` or `callable_variant` |
| `unresolved.bed` | Merged intervals labelled with their primary exclusion reason |

Both BEDs use reference contig names and 0-based, half-open coordinates. Together they partition the scope
exactly. Contigs are grouped in input encounter order, followed by unobserved reference contigs; they are
not necessarily lexicographically ordered. Outputs are staged and published together only after a
successful scan. Invalid/truncated input returns a nonzero exit code and leaves no completed report.

## Convert GTF or GFF3 annotations

Create the gene targets file directly from an annotation:

```sh
gvcf-audit targets --annotation genes.gtf.gz --out targets.tsv
gvcf-audit --gvcf sample.g.vcf.gz --reference reference.fa \
    --gene-targets targets.tsv --out gene-audit
```

GFF3 works the same way. Use `--contig chr22` to select a chromosome, or repeat `--gene ID` to
select exact annotation gene IDs. Both filters can be repeated and, when combined, must both match.
Gene symbols are not substituted for IDs. Unknown requested IDs or contigs produce an error.

The converter uses exon features, subtracts one from each start coordinate, and preserves the end
coordinate on both strands. It sorts the result and removes identical rows. It preserves distinct
exon IDs and overlapping exons; the audit merges their coverage within each gene.

GTF uses `gene_id` and `exon_id`. GFF3 follows `Parent` links to gene features and uses their `ID`,
with exon `ID` as the exon label. Shared exons can belong to multiple transcripts or genes, and
parents may appear after their children. Missing exon IDs get stable coordinate-and-strand labels.
Missing gene assignments are errors. No transcript is selected as canonical; all supplied exons
that pass the filters are included.

Input may be plain text, gzip or BGZF. The filename (`.gtf`, `.gff` or `.gff3`, optionally `.gz` or
`.bgz`) selects the format; a GFF3 version header also works. Use `--format gtf` or `--format gff3`
for other filenames. GFF1/GFF2 are not supported. `--out` is a new TSV **file**, not a directory.

Use the same genome assembly as your gVCF and reference. This command does not rename contigs,
lift coordinates between assemblies, or check sequence identity. The audit checks the resulting
intervals against its FASTA. See [annotation conversion details](docs/annotation-targets.md) for
identifier conventions, supported GFF3 gene types and limits.

## Gene and exon targets

Use `--gene-targets targets.tsv` instead of `--bed` to summarize named exons and genes. The file can be
plain or gzip-compressed and must have exactly these five tab-separated columns, including the header:

```text
contig	start	end	gene	exon
chr1	10000	10200	GENE_A	exon_1
chr1	11000	11200	GENE_A	exon_2
```

Coordinates are **0-based, half-open**. Supply targets from the same assembly as the reference.
This is an annotated TSV, not BED5 (whose fifth column normally means score). Blank lines and lines
starting with `#` are ignored; gene and exon identifiers must be present and cannot be `.`.

```sh
gvcf-audit --gvcf sample.g.vcf.gz --reference reference.fa \
    --gene-targets targets.tsv --out gene-audit
```

This adds `genes.tsv`, `exons.tsv`, and searchable gene/exon tables in HTML. JSON includes `genes` and
`exons`, each with callable and unresolved bases, the main exclusion reason, and all state counts.
Genes also report how many targeted exons have any unresolved bases. These are summaries of the
**supplied targets**, not claims about every exon or the full genomic span of a gene.

Repeated or overlapping intervals are merged within each gene and within each exon identifier. Genes
are keyed by reference contig and gene identifier; exons by contig, gene and exon identifier. Repeated
exon identifiers are treated as one exon, even when their intervals are disjoint. Use distinct identifiers
if you want transcript-specific exons kept separate. Different genes and exon identifiers may overlap,
so their counts must not be summed to obtain overall totals. Introns are not filled in between targets.

Try `--gene-targets examples/gene-targets.tsv` with the included example input and reference.

## Investigating exclusions

The HTML report shows record examples with coordinates, input line numbers, observed GT/DP/MIN_DP/GQ
and filters, and the rule that excluded them. Reference mismatches include the first differing base and
its 0-based coordinate. A separate set of examples shows final unresolved intervals, including gaps,
overlaps and ambiguous reference bases. The same information appears under `diagnostics` in JSON.

Examples are limited to five per reason and type by default. Set `--max-examples 0` to disable them or
choose a limit up to 100. Text fields in record examples are shortened to 80 characters; the first
reference mismatch is still located using the full REF allele. Examples are taken in encounter order,
not randomly sampled, and only records overlapping the requested scope are retained. Their full record
spans may extend beyond the targets. They show primary record reasons plus missing MIN_DP, not every
simultaneous failure. All records are still counted and checked when examples are disabled.

## Compare two audits

Compare saved audit directories without reading their original gVCFs or FASTAs again:

```sh
gvcf-audit compare --before old-audit --after new-audit --out comparison
```

Both directories must contain `report.json`, `callable.bed` and `unresolved.bed`. Complete reports from
0.2.0 onward are supported. Gene/exon comparisons require both audits to have been generated with the
same `--gene-targets` annotation. The output directory must be new.

The comparison counts positions that **gained** or **lost** callability, as well as the net change.
Ten gained bases and ten lost bases therefore remain visible even when the total is unchanged.
Changes between two callable states are recorded but are not gains or losses of callability.
The offline report includes searchable gene/exon tables and before/after reason totals.

- `comparison.json`: overall, contig, gene/exon and state-transition counts, conditions and provenance.
- `genes.tsv` and `exons.tsv`: before/after callability, gained/lost bases and net changes. Header-only
  files are written when the input reports have no gene/exon summaries.
- `changes.bed`: BED4 containing every changed primary state, named `before_state->after_state`.
- `report.html`: offline comparison report; keep it with the downloads.

The target union must match exactly. Group identifiers and each group's target union must also match;
input-row order and repeated rows can differ. Target differences are refused rather than silently
restricting the comparison to an intersection. Gene/exon totals can overlap and must not be summed.

Sample names, recorded reference metadata and quality rules must agree by default. Intentional
comparisons can use these explicit options, with the conditions retained in the report:

- `--allow-policy-change`: compare different depth/GQ thresholds, filter settings or caller policies.
  For example, add this flag when comparing audits made with `--min-gq 20` and `--min-gq 30`.
- `--allow-different-samples`: compare different sample names. Differences cannot be attributed to a
  pipeline change alone.
- `--assume-same-reference`: assert that identical reference sequences were used when their recorded
  paths, sizes or modification times differ, such as after relocating a FASTA.

Matching reference metadata does **not** prove sequence identity. No reference checksums are available
or computed, and the report says so. The reference override is an assertion, not a lift-over operation.
These comparisons describe reported genotype evidence, not biological accuracy.

The command validates both BED partitions and their global, contig and group counts against the reports.
It reads audit BEDs once to index contig byte ranges and once to compare them, allowing different contig
orders without keeping genome-wide intervals in memory. Memory still scales with report targets and
groups. Keep the input audit directories unchanged while comparing them.

Exit 0 means a valid comparison was written, even if callability fell or an input audit failed its
quality gate. Such gate failures are shown as warnings. Invalid/incompatible input returns 1 and
leaves no completed output directory. There is no automatic regression threshold in this version.
See [the comparison format](docs/comparison-format.md) for field definitions.

## Pipeline use

Set explicit quality limits when an automated workflow needs a pass/fail decision:

```sh
gvcf-audit --gvcf sample.g.vcf.gz --reference reference.fa --bed targets.bed \
    --out sample-audit --min-callable-percent 95 --max-reference-mismatch-bases 0 --quiet
```

- Exit `0`: processing completed and any requested limits passed.
- Exit `1`: invalid arguments/input or an I/O error; no completed output directory is published.
- Exit `2`: processing completed but a requested quality limit failed. The complete report is kept.

The limits apply to the union of requested regions. The reference-mismatch limit counts the primary
`reference_mismatch` status; it is not a count of every mismatching record under overlaps or ambiguous
reference bases. Mismatching records anywhere in the input also produce a report warning.
Use `--max-reference-mismatch-records 0` to fail on any mismatching record whose span overlaps the
requested scope, even when ambiguity or another record masks it in base counts. The report displays
full-input record counts, in-scope record counts, and primary base counts together. A mismatch outside
the scope remains a warning but does not fail either scope limit. The existing base-limit behavior is
unchanged.
No pass/fail limits are enabled by default. `--quiet` suppresses progress; errors and failed limits remain visible.

`report.json` identifies its format as `gvcf-audit-report-v1`. Its schema is in
[`docs/report.schema.json`](docs/report.schema.json). Consumers should allow additional fields within
v1 and check both `status` and `quality_checks`. Chromosome totals appear in `contigs`; overlapping target
rows remain independent in `regions`. See [`docs/report-format.md`](docs/report-format.md).

## What callable means

A callable base is covered by one accepted small-variant genotype record or reference block. Its
reference sequence agrees with the FASTA, its haploid or diploid genotype is complete, its filters are
accepted, and its reported depth and genotype quality meet the selected thresholds.

This measures the **evidence reported by the caller**. It does not prove that every possible variant type
was excluded, or that a caller's confident genotype is correct. In particular, reference-block confidence
does not establish sensitivity to arbitrary insertions, structural variants or repeat expansions.

The policy is deliberately conservative:

- Reference blocks use `MIN_DP`; variant records use `DP`. Both require `GQ`.
- Haploid and diploid calls are accepted without changing dosage. Other ploidies are unresolved.
- Missing genotype quality/depth is unresolved, never treated as zero or filled in.
- Two overlapping records make the overlap unresolved, including records from split multiallelic sites.
- Literal SNVs and equal-length substitutions can be callable. Indel and structural-record spans are
  unresolved; this version does not reconstruct haplotypes or infer callability inside a deletion.
- A gap in a plain VCF is `no_record`, never assumed to be homozygous reference.
- Non-ACGT reference bases are `reference_ambiguous`, even inside passing blocks.
- All reference alleles are checked. Header contig lengths are compared where available; matching lengths
  alone do not establish assembly identity. Use the correct reference yourself.
- Exact contig names take precedence. Unambiguous `chr` and mitochondrial aliases are recognized.

Each base has one primary status. Reference ambiguity takes precedence over overlapping records. For a
single record, the order is reference mismatch, ambiguous REF allele, missing/partial genotype,
unsupported ploidy, genotype/site filter, missing quality, low depth, low GQ, unsupported allele,
complex variant, then a callable status. Counts are not an inventory of every simultaneous problem.

## Caller conventions

`--caller auto` detects supported header signatures; ambiguous or unknown signatures use generic rules.
Override with `--caller deepvariant`, `gatk`, `dragen` or `generic` when the provenance is known.

For identified DeepVariant, GATK and DRAGEN inputs, homozygous-reference blocks with `FILTER=.` are
assessed using MIN_DP and GQ, since site variant filters are not applied to those blocks. This adaptation
is counted in JSON. DeepVariant's hom-ref `RefCall` convention is also recognized and counted.

Variant records otherwise require `FILTER=PASS`. A failing sample `FORMAT/FT` excludes a record.
Two explicit options relax the defaults and are recorded in every report:

- `--allow-unfiltered`: accept `FILTER=.` for other records too. For example, raw HaplotypeCaller variant
  records may have no site filtering yet; use this only when the quality policy is appropriate.
- `--allow-block-dp`: when MIN_DP is absent, use DP for a reference block. DP is not guaranteed to be the
  minimum depth across that block, so the report warns about this choice.

Caller labels are format adaptations, not claims of equivalent calibration across technologies.
Full public genomes were checked for DeepVariant 1.10.0, GATK HaplotypeCaller 4.beta.5-SNAPSHOT,
GATK ReblockGVCF 4.2.2.0 and DRAGEN 3.5.7b. These are the tested versions, not a guarantee for every
release or workflow. See [the validation notes](docs/validation.md) for results and input links.

ReblockGVCF headers are recognized and flagged in the report. Reblocking can remove MIN_DP and reduce
GQ resolution. Blocks without MIN_DP remain `quality_missing` by default. `--allow-block-dp` lets you
audit their retained DP, but cannot recover minimum depth. Unfiltered variant records still need
`--allow-unfiltered` if you choose to accept them. Neither option restores information discarded by
reblocking.

## Current limits

- No indexed target-only reading, cohort aggregation or automatic annotation download.
  Supply a target TSV directly or create one from a matching GTF/GFF3 with `targets`.
- No BCF or compressed FASTA support in this version.
- No automatic repair or normalization of the input.
- Every VCF record and requested BED interval must resolve in the FASTA. Unused header contigs absent
  from the reference are listed in JSON and generate a warning; they are outside the audit scope.
- Input contigs must be contiguous and positions nondecreasing; coordinate sorting is not performed.
- The reference and input files must not be modified while the audit is running.
- Reports contain sample identifiers, paths and genomic intervals. They remain local unless you share them.
- Platform builds beyond the supplied Linux x86-64 binary have not been tested.

## Source layout

`vcf.rs` handles record policy, `reference.rs` reference access, `intervals.rs` interval accounting,
`report.rs` report generation, `genes.rs` target grouping, `targets.rs` annotation conversion,
`diagnostics.rs` bounded examples, and
`main.rs` the CLI and streaming orchestration.

Licensed under MIT. No production genomes or scoring weights are included.
