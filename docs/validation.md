# Reproducible checks

The repository now includes 14 Rust integration tests and 362 synthetic checks in
Python, using only its standard library. See [the test instructions](../tests/README.md).
They cover the audit, comparisons, annotation conversion, gene/exon unions, reference
validation, errors and output protection. CI runs both suites on stable Rust and Rust
1.88, with clippy and formatting on stable. No downloaded genomes are needed.

The release notes below describe checks run at the time of each release. Their statements
about external validation scripts describe those historical releases. The latest synthetic
suites have since been brought into the repository. Full public-genome runs remain separate
historical evidence; they are not part of CI. Test counts are checks, not a coverage percentage.

## Validation for 0.5.1

The two added Rust cases cover soft-masked reference bases and a deletion followed by
an overlapping `*` spanning-deletion record. All 14 Rust integration tests pass.
The 362 synthetic checks pass with the existing fixed seeds and fresh seed
6988092914462208025; replaying that seed also passes on Rust 1.88. A failed suite
prints its seed and replay command. These checks use only small synthetic data.

Main-branch pushes and pull requests run the fixed suites. A separate weekly job
runs fresh seeded cases and records the seed in its log. It does not replace the
fixed suite, and has not yet run on its scheduled trigger at release time.

## Validation for 0.5.0

The annotation converter passed 136 checks, including 30 independent interval models represented in
both GTF and GFF3. Checks cover both strands, single-base exons, repeated rows, shared exons/genes,
parents after children, missing exon IDs, percent escapes, quoted GTF attributes and comments,
malformed links/coordinates, format and gene/contig filters, gzip input and truncated gzip errors.
Existing destinations and simultaneous output writers cannot overwrite a completed target file.
Validation programs remain outside the repository; no regression tests were added.

All 679,562 unique targets (78,733 genes, 25 contigs) from the complete GENCODE v50 basic GTF agree
with an independent conversion of 3,122,490 exon records. The chr22 subset also matches the 15,130
unique targets used in the prior real-genome audit, retaining its 3,690,789-base union. All 7,507
Ensembl 115 yeast GFF3 targets (7,127 genes, 17 contigs) agree with independent parent mapping.
That file exercises missing exon IDs and gene, ncRNA_gene, pseudogene and transposable_element_gene
parents on both strands. This checks coordinate and grouping consistency, not annotation correctness.

The full human conversion took 8.35 seconds and 186,456 KiB peak RSS. The chr22-filtered
conversion took 3.25 seconds and 8,864 KiB; the yeast conversion took 0.05 seconds and 19,004 KiB.
These are single local runs, not cross-machine performance guarantees. Conversion reads the full
annotation even when filtering; it does not read a genome sequence or gVCF.

The example GTF and GFF3 produce identical TSVs. Running an audit on those converted targets gives
identical counts and BED/TSV outputs to the published 0.4.0 binary with independently prepared targets.
The existing comparison example still matches. The 79 existing audit checks and 25 release checks pass.
Release builds pass on Rust 1.88.0 and 1.98.0 with identical example conversion output.

[Machine-readable results and source URLs](annotation-validation.json) describe the tested scope.

## Validation for 0.4.0

The audit command's 79 existing controlled checks and 25 release checks still pass. The new comparison
passed 59 checks, including 30 randomized comparisons against independent per-base calculations, plus
nine compatibility/export checks. These cover reversed contig order, repeated and overlapping groups,
equal totals hiding changed positions, missing/invalid BEDs, inconsistent counts, sample/reference/policy
guards, failed quality gates, empty callable/unresolved files, legacy 0.2.0 reports and BED4 output.
Validation programs remain outside the repository; no regression-test files were added.

The real comparison used the previous NA12878 GATK audit with GENCODE v50 basic chr22 exon targets.
A second complete gVCF scan changed only the minimum GQ from 20 to 30. All 3,690,789 target positions
were compared independently from the two audit BED partitions, including every group total and every
state transition. Exactly 7,612 bases lost callability, affecting 756 genes and 2,037 exons; none gained
callability. Counts for overlapping genes/exons must not be summed. A separate pysam/HTSlib parser
agreed at all 44,804 sampled bases across 200 exons under the new policy.

The final comparison took 0.68 seconds (0.61 user / 0.06 system CPU) and peaked at 215,612 KiB RSS.
This was one local cached run including exports, after the original audits had finished. It measures
the comparison only, not the full gVCF scans, and is not a cross-machine performance guarantee.
No original gVCF or reference was read by the comparison command. Reference identity remains an
explicit limitation because audits have metadata, not sequence checksums.

Release builds passed on Rust 1.88.0 and 1.98.0; their controlled comparison outputs match byte for byte.
46 comparison reports passed the new JSON schema. Browser checks covered search, sort, pagination and
empty results on the real 1,747-gene/15,130-exon report, with no logged JavaScript errors.
[Machine-readable results](comparison-validation.json) record the exact scope and limitations.
These checks establish calculation consistency, not biological accuracy or every caller/version's behavior.

## Validation for 0.3.0

The new diagnostics and gene/exon summaries were checked on Linux x86-64. No changes to the
callability policy were intended. The previous callability policy is retained.

- All 79 previous controlled checks, 25 release checks and three reblocking checks passed again.
- 61 new checks passed, including 45 randomized gene/exon union comparisons. These cover repeated
  and overlapping targets, shared gene coordinates, scope boundaries, malformed targets, bounded
  examples, reference mismatch details, HTML escaping and the separate mismatch-record limit.
- An additional overlapping-record example confirmed that the record limit still fails when the
  final interval is classified as an overlap instead of a reference mismatch.
- 60 JSON reports validated against the expanded schema, including a saved 0.2.0 report.
- Gene/exon search, sorting, pagination and empty results were exercised in Chrome with up to
  125 groups. These checks do not establish compatibility with every browser.

The complete reblocked NA12878 input described below was scanned again: 42,025,432 records across
3,217,346,917 requested bases. Both BED outputs and regions.tsv match the saved 0.2.0 outputs byte
for byte. Global, contig and record classifications are unchanged. The report now distinguishes
84 mismatching records in scope from zero primary mismatch bases; the separate record limit correctly
returned exit status 2 while retaining the completed report.

That run took 34.15 seconds (33.61 user and 0.51 system CPU seconds), with peak RSS of 253,536 KiB.
This is one local run, not evidence of a speedup. Checksums were not computed. The whole-genome comparison used the existing BED scope.

Gene/exon totals describe only supplied intervals. They do not establish completeness of an annotation
or assess unlisted exons. Validation scripts remain outside the source repository.

## Real gene/exon annotations

[GENCODE v50 basic chromosome annotations](https://www.gencodegenes.org/human/release_50.html)
were used to select every chr22 exon row. GTF starts were reduced by one to produce 0-based half-open
coordinates; ends were unchanged. Versioned gene_id and exon_id values supplied the group identifiers.
Transcript repeats were intentionally retained: 78,148 input rows describe 15,130 distinct gene/exon
pairs across 1,747 genes and 8,209 transcripts, on both strands. This is the basic annotation subset,
not every annotated transcript or an assertion that all possible exons are covered.

The complete public NA12878 HaplotypeCaller input was scanned (154,179,502 records), using its matching
Broad hg38 reference and the default policy. Of 3,690,789 distinct requested bases, 3,554,088 were
callable (96.2962%). Every input row, gene, exon, TSV state count and all 40,036 exported BED intervals
agreed with an independent calculation using sets of base positions and the saved 0.2.0 whole-genome
classifications. Repeated transcript rows did not inflate group totals. A separate pysam/HTSlib parser
also agreed at all 44,804 sampled bases across 200 exons. Independent classification was sampled,
not repeated at every requested base.

The run took 120.97 seconds, with 118.72 user and 0.94 system CPU seconds and 830,252 KiB peak RSS.
The many retained rows produced a 77 MiB JSON report and a 15 MiB HTML report. Memory and report size
scale with target rows and group counts; the smaller earlier BED audits are not a memory estimate for
large annotation sets. Search worked in Chrome on the real 1,747-gene/15,130-exon report without logged
JavaScript errors. This evidence covers one caller/sample and annotation release. Machine-readable
results and the annotation source URL are in [gene-validation.json](gene-validation.json).

The earlier release evidence and its limits follow.

## Validation for 0.2.0

This release was validated on Linux x86-64. These checks establish behavior under the documented policy;
they do not establish the biological accuracy of an input caller or clinical suitability.

## Full whole-genome run

Input: [public HG003 DeepVariant 1.10.0 WGS gVCF](https://storage.googleapis.com/deepvariant/case-study-outputs/1.10.0/wgs/wgs_HG003.deepvariant.g.vcf.gz),
1,023,866,857 compressed bytes. A local GRCh38 FASTA including EBV supplied the reference.
A BED covered all 195 declared contigs, totalling 3,099,922,541 bases. The gVCF itself has records on
25 contigs; declared contigs without records remain unresolved. This is an intentional whole-reference
scope, not a gene-panel coverage estimate.

| Measurement | Result |
|---|---:|
| Input records | 92,732,214 |
| Elapsed time | 69.33 seconds |
| User CPU | 68.58 seconds |
| System CPU | 0.70 seconds |
| Peak resident memory | 252,904 KiB (about 247 MiB) |
| Callable bases under default thresholds | 2,817,532,549 (90.8904%) |
| Reference-mismatch bases | 0 |
| Exported BED intervals | 32,913,514 |

This was one final run on an internal NVMe with cached inputs, including FASTA layout validation and
all exports. It is not a cross-machine performance guarantee. The uncompressed BED outputs total about
1.2 GiB; reserve output space for a full-genome audit. File checksums were not computed.

The BEDs were independently merged and checked to partition every requested base exactly once.
All global and contig totals agree with JSON. A separate implementation using pysam/HTSlib checked
458,903 sampled bases across 459 windows against the input records and reference. Every sampled status
matched, including gaps, overlapping records, complex alleles, ambiguous reference, low depth and low GQ.
Independent genotype classification was sampled; it was not repeated for every base in the genome.

Before limiting reference mapping to one contig, the same scan took 67.40 seconds and peaked at
3,204,888 KiB. BED and TSV outputs before and after that change match byte for byte. These two runs
support the memory improvement, not a claim of improved elapsed time.

## Full GATK and DRAGEN runs

All four runs below scanned each complete downloaded gVCF. BED scopes covered every declared contig
present in the matching reference, including contigs without records. The defaults were DP/MIN_DP >= 10
and GQ >= 20; the last row explicitly relaxed two policies as noted below.

| Input | Records | Elapsed | User / system CPU | Peak RSS | Callable / requested bases |
|---|---:|---:|---:|---:|---:|
| NA12878 HaplotypeCaller 4.beta.5-SNAPSHOT | 154,179,502 | 127.95 s | 125.39 / 1.04 s | 252,812 KiB | 2,803,398,314 / 3,105,715,063 |
| HG00096 DRAGEN 3.5.7b | 177,086,704 | 182.79 s | 179.04 / 1.43 s | 253,616 KiB | 2,689,763,407 / 3,215,250,450 |
| NA12878 ReblockGVCF 4.2.2.0, defaults | 42,025,432 | 36.64 s | 33.21 / 0.73 s | 253,084 KiB | 0 / 3,217,346,917 |
| Same reblocked file, explicit DP and unfiltered options | 42,025,432 | 34.34 s | 33.77 / 0.56 s | 253,148 KiB | 2,797,364,341 / 3,217,346,917 |

These are individual local NVMe runs, including reference validation and all exports. Cache state and
other machine activity were not controlled. Different samples, references and policies mean the
callable fractions and runtimes should not be used to rank callers. No file checksums were computed.

The original GATK file has unfiltered variant records, which remain `filter_not_assessed` under the
strict defaults. The reblocked file additionally lacks MIN_DP in all 37,003,552 reference blocks.
Its zero default callable count reflects missing required evidence, not a claim of zero sequencing
coverage. With `--allow-block-dp --allow-unfiltered`, retained DP and GQ can be evaluated, but DP does
not establish minimum depth. Reblocking can also reduce GQ resolution; see
[Broad's reblocking description](https://broadinstitute.github.io/warp/docs/Pipelines/Whole_Genome_Germline_Single_Sample_Pipeline/README#reblocking).
The retained header identifies ReblockGVCF 4.2.2.0; it does not establish the original caller version.

Every exported BED interval was checked for a complete partition, and every global and contig total
was compared with JSON. Independent pysam/HTSlib classification agreed at 550,000 sampled bases for
GATK, 975,903 for DRAGEN and 550,000 for each reblocked policy. Classification was sampled, not checked
independently at every base. The corresponding interval counts were 23,132,985, 29,663,085, 9,880,970
and 19,337,964. Machine-readable summaries are in `validation-results.json` beside this document.

The raw GATK and DRAGEN files each contain 92 REF mismatches, and the reblocked file contains 84.
Every mismatch was individually located: the VCF uses N where the FASTA uses another IUPAC ambiguity
code, such as R, Y or M. Their entire spans are ambiguous reference bases and correctly remain
`reference_ambiguous`. The primary `reference_mismatch` base count is therefore zero, while the record
counts and warnings retain the discrepancies. This is not a claim of exact REF agreement everywhere.

### Public input sources

- [GATK NA12878 gVCF](https://gatk-test-data.s3.amazonaws.com/wgs_gvcf/PlatinumGenomes_hg38/NA12878.g.vcf.gz)
  and [Broad hg38 reference](https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta).
- [Reblocked NA12878 gVCF](https://storage.googleapis.com/gatk-test-data/wgs_reblocked_gvcf/NA12878.rb.g.vcf.gz),
  using that same Broad reference.
- [DRAGEN HG00096 gVCF](https://1000genomes-dragen.s3.amazonaws.com/data/dragen-3.5.7b/hg38_altaware_nohla-cnv-anchored/HG00096/HG00096.hard-filtered.gvcf.gz)
  and [the published DRAGEN reference](https://1000genomes-dragen.s3.amazonaws.com/reference/hg38_alt_aware_nohla.fa).

The index files are available at these URLs with `.tbi` or `.fai` appended. Indexes were used by the
independent validator; the audit itself does not require a VCF index.

## Caller examples and controlled checks

- Three [DeepVariant r1.10 fixtures](https://github.com/google/deepvariant/tree/r1.10/deepvariant/testdata):
  standard, haploid and median-depth gVCFs, each containing 228 records over 10,001 bases.
- [GATK 4 fixture](https://github.com/broadinstitute/gatk/blob/master/src/test/resources/org/broadinstitute/hellbender/tools/haplotypecaller/expected.testGVCFMode.gatk4.g.vcf):
  1,291 records over a 100,001-base window.
- [GATK 3.5 fixture](https://github.com/broadinstitute/gatk/blob/master/src/test/resources/org/broadinstitute/hellbender/tools/haplotypecaller/expected.testGVCFMode.gatk3.5.g.vcf):
  1,572 records over that window. GATK examples used an explicit caller profile because their retained
  headers do not uniquely identify the producer.
- [DRAGEN fixture in GATK](https://github.com/broadinstitute/gatk/blob/master/src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/ReblockGVCF/dragen.g.vcf):
  five records. This checks format handling, including unused HLA header contigs, with an explicit
  DRAGEN profile. The complete DRAGEN 3.5.7b genome above provides additional evidence.
- All six positive caller fixtures had zero reference mismatches against their supplied references.
- 79 controlled checks passed, including 45 randomized comparisons against an independent per-base
  interval calculation. Another 25 release checks covered quality limits, exit codes, stale indexes,
  CRLF references, missing final newlines, malformed fields, missing contigs and HTML escaping.
- Three additional controlled ReblockGVCF checks passed: strict defaults, explicit DP fallback, and
  DP fallback plus unfiltered records. GQ=0 blocks remained unresolved.
- Plain VCF and BGZF forms of a controlled example produced matching counts and region results.
- The JSON schema and example, whole-genome and quality-limit reports were validated independently.
- HTML search, sorting and pagination were exercised in Chrome.
- Release builds and Clippy checks passed. Both Rust 1.88.0 and 1.98.0 built the release successfully.

Validation programs and downloaded data are kept outside the source repository. The included `examples`
directory is a runnable demonstration, not a regression-test suite.

Two additional downloaded examples were not counted as positive compatibility evidence: a DeepVariant
PacBio fixture had REF mismatches against the supplied reference, and a reblocked DRAGEN chrY fixture
contained a repeated column header and was rejected. Inputs were not silently repaired.

## Limits

Full-genome evidence covers the exact versions and inputs above. Other versions,
technologies, operating systems and malformed-input patterns may expose additional issues. Indels,
structural alleles and overlaps deliberately remain unresolved. The report describes a quality policy
for genotype evidence, not comprehensive variant sensitivity.
