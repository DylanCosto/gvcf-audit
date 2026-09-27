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

## Additions in 0.3.0

`reference_mismatches` reports `all_input_records`, `records_overlapping_scope`, and
`primary_bases_in_scope`. The new record gate uses the middle count; the existing base gate uses the
last. Record overlap is based on the entire record span, not just the position of the differing REF
base. All counts remain available when examples are disabled.

`diagnostics.record_examples` stores bounded examples before overlap/reference masking. The reason
`reference_blocks_missing_MIN_DP` is an additional observation and can accompany another record reason.
Fields under `observed` contain the selected sample's values; numeric nulls mean missing. String values
are shortened to 80 characters with an ellipsis when truncated. `first_reference_mismatch.position`
is 0-based and is found using the full REF allele. `line` is the 1-based line number in the decompressed
VCF. Record examples must overlap the requested scope but retain their full span.

`diagnostics.interval_examples` contains bounded final unresolved intervals within the scope. These
can differ from a record's reason because overlap and ambiguous reference take precedence. The first
examples encountered are retained, with adjacent pieces of a retained interval merged. They are not
an unbiased sample. The limit applies separately to each record/interval reason.

`genes` and `exons` are empty unless `--gene-targets` is supplied. Group identity is contig + gene for
genes and contig + gene + exon for exons. Counts cover the union of that group's intervals. The JSON
`intervals` array contains the actual merged targets; gaps between them are not counted. Gene counts
are calculated from the union, never by adding overlapping exon counts. `targeted_exons` counts unique
exon identifiers, and `unresolved_exons` counts those with any unresolved base. `main_reason` is the
unresolved state with the most bases (null if all bases are callable); ties follow the state order in
`STATES` in the source. All state counts are retained so ties need not be inferred from the main reason.
Groups are ordered by reference contig, gene identifier and exon identifier. The HTML initially sorts
by lowest callability. `regions` still reports input rows independently, now with `gene` and `exon`
fields for annotated input; `regions.tsv` keeps its existing columns.

These fields are optional in the v1 schema so older reports remain valid. Quality rules and existing
BED/TSV outputs for unchanged inputs and options retain their previous meanings.
