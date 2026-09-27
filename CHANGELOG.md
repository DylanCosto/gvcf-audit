# Changes

## 0.3.0

- Add bounded record and interval examples with observed values and exclusion rules.
- Show full-input and in-scope reference mismatch records alongside primary base counts.
- Add an optional mismatch-record limit without changing the existing base limit.
- Add explicit gene/exon targets, summaries and searchable offline tables. Overlaps count once within each group.
- Validate gene/exon totals against GENCODE v50 chromosome 22 annotations and a public GATK genome.

## 0.2.0

- Validate FASTA names, lengths, offsets and line layout against its index.
- Keep only one reference contig mapped at a time.
- Add optional callability and reference-mismatch limits with distinct exit codes.
- Add chromosome totals and a versioned JSON schema.
- Make report regions searchable, sortable and paginated offline.
- Allow unused header contigs missing from the reference, with a warning and a recorded list.
- Reject malformed depth values and mixed missing ALT alleles.
- Recognize older HaplotypeCaller header signatures.
- Recognize ReblockGVCF headers and warn when reference blocks lack MIN_DP.
- Validate complete public GATK, reblocked GATK and DRAGEN genomes alongside DeepVariant.
- Include a small runnable example and installation instructions.

## 0.1.0

Initial prototype with streaming gVCF auditing, reason-labelled BEDs, target summaries and an offline report.
