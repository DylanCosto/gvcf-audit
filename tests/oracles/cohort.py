"""Check cohort matrices against independent per-base calls and target unions."""

import csv
import json
from pathlib import Path
import random
import subprocess

from support import random_seed, workspace

repo, study, exe, temporary = workspace()
rng = random.Random(random_seed(60271))
sequences = {"chr1": "A" * 50 + "N" + "A" * 49, "chr2": "C" * 100}
reference = study / "reference.fa"
with reference.open("wb") as fasta, Path(str(reference) + ".fai").open("w") as fai:
    for contig, sequence in sequences.items():
        fasta.write(f">{contig}\n".encode())
        offset = fasta.tell()
        fasta.write((sequence + "\n").encode())
        fai.write(f"{contig}\t100\t{offset}\t100\t101\n")

rows = []
gene_positions, exon_positions = {}, {}
for contig in sequences:
    for gene in range(6):
        for exon in range(2):
            start = rng.randrange(90)
            end = rng.randrange(start + 1, 101)
            rows.append((contig, start, end, f"G{gene}", f"E{exon}"))
            gene_positions.setdefault((contig, f"G{gene}"), set()).update(range(start, end))
            exon_positions[(contig, f"G{gene}", f"E{exon}")] = set(range(start, end))
rows.append(rows[0])  # Repeated annotation rows must not inflate the denominator.
targets = study / "targets.tsv"
samples = ["001", "#sample", 'sample "quoted"'] + [f"S{i}" for i in range(9)]
expected = {}
audits = []
for i, sample in enumerate(samples):
    rng.shuffle(rows)
    targets.write_text("contig\tstart\tend\tgene\texon\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))
    vcf = study / "sample.g.vcf"
    text = "##fileformat=VCFv4.2\n##source=DeepVariant\n"
    text += "".join(f"##contig=<ID={c},length=100>\n" for c in sequences)
    text += f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
    calls = {}
    for contig, sequence in sequences.items():
        calls[contig] = set()
        for pos, base in enumerate(sequence):
            if rng.random() < 0.15:
                continue
            dp, gq = rng.choice([2, 20, 30]), rng.choice([10, 20, 30])
            text += f"{contig}\t{pos + 1}\t.\t{base}\t<NON_REF>\t.\t.\tEND={pos + 1}\tGT:MIN_DP:GQ\t0/0:{dp}:{gq}\n"
            if base in "ACGT" and dp >= 10 and gq >= 20:
                calls[contig].add(pos)
    vcf.write_text(text)
    out = study / f"audit-{i}"
    result = subprocess.run([str(exe), "--gvcf", str(vcf), "--reference", str(reference),
                             "--gene-targets", str(targets), "--out", str(out), "--quiet"], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    audits.append(out)
    expected[sample] = calls

def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

for level, positions in [("gene", gene_positions), ("exon", exon_positions)]:
    out = study / level
    result = subprocess.run([str(exe), "summarize", *map(str, audits), "--level", level,
                             "--min-callable-percent", "37.5", "--out", str(out)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    counts = table(out / "callable_bases.tsv")
    percents = table(out / "callable_percent.tsv")
    summary = table(out / "groups.tsv")
    assert len(counts) == len(positions) == len(percents) == len(summary)
    for key, count, pct, group in zip(sorted(positions), counts, percents, summary):
        assert (count["contig"], count["gene"]) == key[:2]
        if level == "exon":
            assert count["exon"] == key[2]
        total = len(positions[key])
        assert int(count["requested_bases"]) == total
        values = []
        for sample in samples:
            called = len(positions[key] & expected[sample][key[0]])
            assert int(count[f"sample:{sample}"]) == called
            percent = 100 * called / total
            assert abs(float(pct[f"sample:{sample}"]) - percent) <= 0.00000051
            values.append(percent)
        for field, value in [("min_callable_percent", min(values)), ("mean_callable_percent", sum(values) / len(values)), ("max_callable_percent", max(values))]:
            assert abs(float(group[field]) - value) <= 0.00000051
        assert int(group["samples_below_threshold"]) == sum(v < 37.5 for v in values)
        assert int(group["samples_fully_callable"]) == sum(v == 100 for v in values)
    assert [r["sample"] for r in table(out / "samples.tsv")] == samples
    assert json.loads((out / "cohort.json").read_text())["sample_count"] == len(samples)
print("2 randomized cohort comparisons passed (12 samples, gene and exon matrices).")
