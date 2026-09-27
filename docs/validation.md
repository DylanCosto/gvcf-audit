# Validation for 0.2.0

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
