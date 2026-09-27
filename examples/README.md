# Example

This synthetic 1,000-base example contains 600 callable bases. Run it with the command in the main README,
or download the example-report archive from the release and open `report.html` locally.

The saved report is in `report/`. Its JSON input paths were shortened for sharing; counts and policy
are unchanged. The report describes invented data, not a real person's genome.

For the gene/exon summaries, replace `--bed examples/targets.bed` with
`--gene-targets examples/gene-targets.tsv` and choose a new output directory.
The gene names are invented. Overlapping exon targets demonstrate that each gene's bases are counted once.
The saved `report/` directory is the original 0.2.0 example; run the current binary to see the new features.
