# Example

This fictional input spans 1,000 bases. The BED example in the main README contains 600 callable bases
out of 1,000. The gene-target example examines 750 unique bases, of which 450 are callable.

[Browse both reports online](https://DylanCosto.github.io/gvcf-audit/), or download the example-reports
archive from the release and open `index.html` locally.

The saved report is in `report/`. Its JSON input paths were shortened for sharing; counts and policy
are unchanged. The report describes invented data, not a real person's genome.

For the gene/exon summaries, replace `--bed examples/targets.bed` with
`--gene-targets examples/gene-targets.tsv` and choose a new output directory.
The gene names are invented. Overlapping exon targets demonstrate that each gene's bases are counted once.
The saved `report/` directory was generated with 0.4.0 using the gene targets.
The same report and a stricter-policy comparison are in `docs/demo/`.

To demonstrate a policy comparison with these fictional targets:

```sh
gvcf-audit --gvcf examples/sample.g.vcf --reference examples/reference.fa \
  --gene-targets examples/gene-targets.tsv --out before-audit
gvcf-audit --gvcf examples/sample.g.vcf --reference examples/reference.fa \
  --gene-targets examples/gene-targets.tsv --min-gq 40 --out after-audit
gvcf-audit compare --before before-audit --after after-audit \
  --allow-policy-change --out example-comparison
```

The explicit policy override keeps the changed GQ threshold visible in the report.

The comparison loses 200 callable bases when minimum GQ rises from 20 to 40. This is the expected
effect of stricter requirements, not evidence that the underlying sequence changed. All saved example
reports use shortened input paths for sharing; counts, input sizes, timestamps and quality rules are unchanged.
