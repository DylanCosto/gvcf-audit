# Annotation targets

`gvcf-audit targets` converts exon annotations into the five-column TSV accepted by `--gene-targets`.
It does not need a gVCF or a FASTA. Choose annotations from the same assembly as the subsequent audit.

```sh
gvcf-audit targets --annotation genes.gtf.gz --contig chr22 --out chr22.tsv
gvcf-audit targets --annotation genes.gff3.gz --gene gene:ENSG000001 --out selected.tsv
```

The gene ID in the second command is illustrative: use the exact ID found in your file, including
any prefix or version. `--gene` selects IDs, not symbols or display names. Repeat `--gene` and `--contig`
for multiple selections. Within each option selections are combined; when both options are present,
both must match. Each requested selector must have at least one output exon or the command fails.
Without filters, all exon features in the annotation are included. All transcript isoforms are included.

## Coordinates and identifiers

Only `exon` and its Sequence Ontology accession `SO:0000147` produce targets. CDS, UTR, gene and
transcript spans do not become targets. There is no intron filling, padding or inferred exon structure.

GTF/GFF3 coordinates are 1-based inclusive. The TSV uses `start - 1` and unchanged `end` on both strands.
A one-base exon at position 7 becomes `[6, 7)`. Strand is not a column in the audit target format;
negative-strand coordinates are not reversed.

- **GTF:** `gene_id` is required on each exon; `exon_id` is used when present. Attributes may have
  quoted values or single unquoted values such as `exon_number 1`. Repeated unrelated attributes
  such as GENCODE `tag` are allowed. End-of-line comments outside quoted values are ignored.
  Duplicate gene/exon identifier attributes are errors.
- **GFF3:** gene groups use the `ID` of an ancestor `gene`, `pseudogene`, `ncRNA_gene` or
  `transposable_element_gene`. The first three also accept `SO:0000704`, `SO:0000336` and `SO:0001263`.
  Exon labels use the exon feature's `ID`. `Name`, `gene_id` and `exon_id` attributes do not replace
  GFF3 feature IDs. This keeps GFF3's explicit parent relationships authoritative. Exons may point
  directly to a gene or reach it through transcripts or other intermediate features.
- **Missing exon IDs:** the label is `gvcf-audit:exon:START-END:STRAND`, using the converted coordinates.
  These labels are grouped within each contig and gene. Identical unlabeled exons shared by transcripts
  therefore collapse. This prefix is reserved; explicit exon IDs using it are rejected to avoid collisions.

IDs are preserved, including versions already in the ID. Separate version attributes are not appended.
Gene symbols can be duplicated, so they are not used for grouping. GTF and GFF3 exports of the same
annotation may use different exon IDs; identical coordinates alone do not make their exon labels equivalent.

GFF3 `Parent` lists are split before percent decoding. An encoded comma inside an ID remains part of
that ID, and `+` remains `+`. Shared exons produce one target per distinct gene ancestor. Parents may
appear later in the file. Repeated IDs for discontinuous features are accepted only when feature type,
contig, strand and parent IDs agree. Their separate intervals are retained. Missing parents, cycles or
incompatible contigs/known strands on the exon ancestry are errors, not guesses. Chains longer than
256 levels are refused. Unrecognized gene feature types are not inferred from their names.

The output is sorted by contig (lexically), numeric start, end, gene and exon. Identical five-column rows
are removed; distinct IDs, overlapping intervals and separate parts of a discontinuous exon remain.
The audit merges overlapping bases within each group. This can reduce the number of region rows compared
with a transcript-expanded annotation without changing the gene/exon interval unions.

## Input and output behavior

The format is selected by `.gtf`, `.gff3` or `.gff`, optionally followed by `.gz` or `.bgz`, or by a GFF3
version header. Other filenames require `--format gtf` or `--format gff3`. A GFF version header that
conflicts with GTF selection is an error. GFF1 and GFF2 are not supported. Compression is detected from
the file contents, including concatenated gzip/BGZF members. In GFF3, `##FASTA` ends the annotation;
embedded sequence is not used or validated.

The converter scans the annotation once. It retains unique target rows and, for GFF3, feature IDs,
parent links and exon records on selected contigs. Memory grows with these data. A contig filter reduces
retained data but does not avoid reading the rest of the annotation. Basic column, coordinate and strand
checks apply before contig filtering. Attribute and parent checks apply on selected contigs; GTF checks
exon attributes, while GFF3 retains links from other feature types too. Parent resolution precedes the
gene filter. This is a target converter, not a complete GTF/GFF3 validator.

If supplied, GFF3 `##sequence-region` bounds constrain output intervals. Circular-origin wrapping is
not performed. Without declared bounds, reference-length checks happen during the subsequent audit.
The converter does not infer assemblies, rename contigs, compare sequences or perform liftover.

`--out` is a new, plain TSV file. Its parent directory is created if needed. Existing files, directories
and symlinks are never replaced. Output is written to a temporary file, flushed, then published with a
hard link on the same filesystem; a filesystem without hard-link support returns an error. Errors leave
no completed target file. Exit 0 means success, exit 1 means failure. `--quiet` hides the completion
summary. The summary counts exon records on selected contigs before gene filtering, including repeats,
then reports the number of unique output rows and contig/gene groups.

## Format references

- [GENCODE GTF format](https://www.gencodegenes.org/pages/data_format.html)
- [Sequence Ontology GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
