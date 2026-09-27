# gvcf-audit

Explain where a sample has usable genotype evidence, and why other regions remain unresolved.

**Version 0.2.0.** A local command-line tool for auditing one selected sample using explicit depth,
genotype-quality and filter rules. It has no imputation, variant calling or clinical interpretation.

## Install and run

Download the Linux binary or source archive from the
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
  Without it, the scope is **every FASTA contig**, including unobserved decoys and alternate contigs.
- `--min-dp 10 --min-gq 20`: defaults; change them to suit the assay.

Overlapping BED rows get independent summaries. Overall totals count their union once. BED12 blocks are
not interpreted: split exons into separate BED intervals first. Extra BED columns are ignored.

The tool makes one sequential scan of the complete VCF, even when a BED is supplied.
No VCF index or checksum pass is required. It uses interval boundaries rather than expanding gVCF blocks
into per-base records. The FASTA layout is checked against its index once, then one contig at a time is memory-mapped
and scanned for ambiguous reference runs. Memory also holds BED intervals, reference ambiguity runs,
and simultaneously overlapping records.

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

- No indexed target-only reading, cohort aggregation, run-to-run comparison or automatic gene annotation.
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
`report.rs` report generation, and `main.rs` the CLI and streaming orchestration.

Licensed under MIT. No production genomes or scoring weights are included.
