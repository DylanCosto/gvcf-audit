# Changes

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
