# Changes

## 0.5.1

- Add small Rust fixtures, portable synthetic checks and CI on stable Rust and Rust 1.88.
- Warn in the report and terminal when unset site filters exclude records from callable totals.
- Cover soft-masked references and spanning deletions, and add weekly checks with replayable random seeds.
- Run push CI only on main to avoid duplicate pull-request checks.

## 0.5.0

- Convert GTF and GFF3 exon annotations into gene targets with the `targets` command.
- Follow GFF3 parent links, handle shared exons, and filter by contig or gene ID.
- Preserve annotation IDs, generate coordinate labels when exon IDs are absent, and remove duplicate targets.

## 0.4.0

- Compare saved audits with exact callability gains/losses, state transitions and gene/exon summaries.
- Check matching targets and flag intentional sample, reference-metadata and quality-rule differences.
- Add changed intervals, comparison JSON/TSVs and an offline report without rescanning gVCFs.
- Clarify what the audit measures and add browsable examples using fictional data.

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
