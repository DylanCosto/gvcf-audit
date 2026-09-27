# Comparison format v1

`gvcf-audit compare` writes `comparison.json` with schema `gvcf-audit-comparison-v1`.
Its JSON schema is [comparison.schema.json](comparison.schema.json). The original audit JSON format
is unchanged. All comparisons have `status: complete`; invalid comparisons are not published.

`before` and `after` record each audit directory, version, sample/caller metadata, recorded reference
provenance, quality rules and quality-gate results. Source paths are informational; original gVCF and
FASTA files need not exist. `warnings` records inherited audit warnings and any comparison conditions.
An audit with a failed quality gate may still be compared because its processing completed.

`summary` counts each position in the shared target union once. `contigs` partitions that union.
`genes` and `exons` count each group's actual merged target intervals once; their totals can overlap.
Each row contains:

- `bases`: the shared target length.
- `before_callable_bases`, `after_callable_bases` and corresponding fractions in [0,1].
- `gained_callable_bases`: unresolved before, callable after.
- `lost_callable_bases`: callable before, unresolved after.
- `net_callable_bases`: after minus before, also gained minus lost.
- `callable_percentage_point_change`: 100 times the difference in callable fractions.
- `before_bases_by_state`, `after_bases_by_state` and `bases_by_state_delta` (after minus before).

Group identity is contig + gene, or contig + gene + exon. Group intervals must match after unioning;
renamed genes/exons are not inferred or remapped. Reports without gene/exon summaries produce empty
arrays and header-only TSVs. TSVs include percentages (0–100), not fractions.

`transitions` counts each before/after state pair across the entire target union, including unchanged
states. It includes only pairs with nonzero bases. Both callable states count as callable; a change
from callable reference to callable variant is a state change but not a callability gain or loss.

`changes.bed` is BED4: contig, start, end, name. The name is `before_state->after_state`. It contains
only primary state changes. Coordinates are 0-based half-open, contigs are lexicographically ordered,
and intervals are ascending within each contig. Do not assume adjacent equal transitions are merged.
Gene groups never duplicate rows in this BED.

`compatibility` records equal target/group unions, the three override flags, and
`reference_identity_verified: false`. Reference provenance uses paths, byte sizes and modification
times, not sequence checksums. `--assume-same-reference` records the user's assertion only. No reference
or sample identity is inferred from matching callable totals.

The reader requires complete supported v1 audit reports, known state names, consistent counts and
complete, non-overlapping BED partitions. Target unions and group targets must match. It accepts
additional report fields; quality-gate settings and example limits do not change per-base policy.
Different depth/GQ/filter settings, detected caller policies or policy definitions require
`--allow-policy-change`. Unknown policy IDs remain unsupported even with that flag.

Inputs must stay unchanged during the comparison. This format reports evidence changes, not their
cause, biological correctness or statistical significance. Numeric counts are limited to signed
64-bit target totals so signed deltas remain representable.
