# Report format v1

The JSON schema is `report.schema.json`. The `schema` value identifies the report format;
`version` identifies the executable. Additive fields may appear without changing the schema identifier.
Consumers should ignore fields they do not use and reject an unfamiliar schema identifier.

`status: complete` means input processing finished successfully. It does not mean the sample passed
quality limits. Inspect `quality_checks.requested`, `quality_checks.passed`, and `quality_checks.failures`
for that decision. A quality-limit failure retains all outputs and returns exit status 2.

Coordinates in reports and BED files are 0-based, half-open. Fractions are between 0 and 1; the TSV
`callable_percent` column uses 0 to 100. Counts are integers. The JSON schema checks structure and types;
the executable also verifies base-count sums before writing the report.

`bases_by_state` partitions the union of requested regions. `contigs` partitions that same union by
reference contig. `regions` describes each input BED row independently, including overlaps; its counts
must not be summed to obtain overall totals. Without BED input, one region covers each FASTA contig.

`callable.bed` and `unresolved.bed` together cover the requested union without gaps or overlaps. Adjacent
intervals with the same status are merged. Contig names come from the FASTA. Contigs are grouped in VCF
encounter order, followed by unobserved FASTA contigs. Sort explicitly if a downstream tool needs a
different contig order. BEDs do not contain track or comment lines.

`record_states_before_overlap_and_reference_mask` counts records, not bases. A reference block can
span many bases, and overlapping records can change the final base-level status. These counts therefore
cannot be directly compared to `bases_by_state`.

`header.unmapped_header_contigs` lists header declarations absent from the supplied FASTA. They have
no actual records, otherwise the run would fail. They are not included in the requested scope.

`inputs` records local paths, sizes and modification times. `checksums_computed` is false: neither the
provenance fields nor a matching REF allele establish whole-file identity. Inputs must remain unchanged
during processing. `options`, `policy`, `adaptations` and `warnings` capture the rules used for the audit.

The HTML embeds its region data and display code. It works offline without a server or network access.
Sample and region names are escaped and displayed as text. The BED, TSV and JSON downloads remain
separate files, so keep the output directory together when sharing it.

`header.reblocked` records whether a ReblockGVCF command header was detected.
`adaptations.reference_blocks_missing_MIN_DP` counts reference blocks missing MIN_DP, including
blocks excluded for another reason. The warning does not enable DP fallback automatically.
