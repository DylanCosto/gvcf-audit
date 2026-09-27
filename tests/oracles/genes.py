import collections
import csv
import gzip
import json
import random
import subprocess
from pathlib import Path

from support import random_seed, workspace

repo, study, exe, temporary = workspace()
results = []
seqs = {"chr1": "A" * 50 + "R" + "A" * 149, "chr2": "C" * 100}
ref = study / "reference.fa"
with ref.open("wb") as f, Path(str(ref) + ".fai").open("w") as idx:
    for c, seq in seqs.items():
        f.write(f">{c}\n".encode())
        offset = f.tell()
        f.write((seq + "\n").encode())
        idx.write(f"{c}\t{len(seq)}\t{offset}\t{len(seq)}\t{len(seq) + 1}\n")
header = (
    "##fileformat=VCFv4.2\n##source=DeepVariant\n"
    + "".join(f"##contig=<ID={c},length={len(seq)}>\n" for c, seq in seqs.items())
    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
)


def record(
    c,
    s,
    e,
    gt="0/0",
    dp="20",
    gq="30",
    filt="PASS",
    base=None,
    fmt="GT:MIN_DP:GQ",
    alt="<NON_REF>",
):
    return f"{c}\t{s + 1}\t.\t{base or seqs[c][s]}\t{alt}\t.\t{filt}\tEND={e}\t{fmt}\t{gt}:{dp}:{gq}"


vcf = study / "input.vcf"
rows = [
    record("chr1", 0, 10),
    record("chr1", 10, 20, dp="2"),
    record("chr1", 20, 30, gq="1"),
    record("chr1", 30, 40, gt="./."),
    record("chr1", 40, 45),
    record("chr1", 42, 48),
    record("chr1", 50, 51, base="N"),
    record("chr1", 70, 80, fmt="GT:DP:GQ"),
    record("chr1", 100, 110, filt="."),
    record("chr1", 150, 160, base="T"),
    record("chr2", 0, 30),
]
vcf.write_text(header + "\n".join(rows) + "\n")


def run(name, flags=(), code=0, program=exe, input=vcf):
    out = study / name
    p = subprocess.run(
        [
            str(program),
            "--gvcf",
            str(input),
            "--reference",
            str(ref),
            "--out",
            str(out),
            *flags,
        ],
        capture_output=True,
        text=True,
    )
    assert p.returncode == code, (name, p.returncode, p.stderr)
    assert out.exists() == (code in (0, 2))
    assert not list(study.glob(".gvcf-audit-*"))
    results.append({"name": name, "passed": True})
    return json.loads((out / "report.json").read_text()) if out.exists() else None, out


new, newout = run("plain")
# Independent per-base labels for the short input above, including overlap precedence.
maps = {c: dict.fromkeys(range(len(seq)), "no_record") for c, seq in seqs.items()}
for start, end, state in [
    (0, 10, "callable_reference"),
    (10, 20, "low_depth"),
    (20, 30, "low_gq"),
    (30, 40, "no_call"),
    (40, 42, "callable_reference"),
    (42, 45, "overlapping_records"),
    (45, 48, "callable_reference"),
    (50, 51, "reference_ambiguous"),
    (70, 80, "quality_missing"),
    (100, 110, "callable_reference"),
    (150, 160, "reference_mismatch"),
]:
    maps["chr1"].update(dict.fromkeys(range(start, end), state))
maps["chr2"].update(dict.fromkeys(range(30), "callable_reference"))
observed = {c: {} for c in seqs}
for name in ["callable.bed", "unresolved.bed"]:
    for line in (newout / name).read_text().splitlines():
        c, start, end, state = line.split("\t")
        for pos in range(int(start), int(end)):
            assert pos not in observed[c]
            observed[c][pos] = state
assert observed == maps
assert new["reference_mismatches"] == {
    "all_input_records": 2,
    "records_overlapping_scope": 2,
    "primary_bases_in_scope": 10,
}
reasons = collections.Counter(
    r["reason"] for r in new["diagnostics"]["record_examples"]
)
assert reasons["reference_mismatch"] == 2
m = [
    e
    for e in new["diagnostics"]["record_examples"]
    if e["reason"] == "reference_mismatch"
][0]
assert m["observed"]["first_reference_mismatch"] == {
    "position": 50,
    "vcf_base": "N",
    "reference_base": "R",
}
for e in new["diagnostics"]["interval_examples"]:
    assert all(maps[e["contig"]][p] == e["reason"] for p in range(e["start"], e["end"]))
for label, bedtext, primary, records in [
    ("masked", "chr1\t50\t51\n", 0, 1),
    ("outside", "chr2\t50\t60\n", 0, 0),
    ("overlap-span", "chr1\t155\t157\n", 2, 1),
]:
    p = study / f"{label}.bed"
    p.write_text(bedtext)
    d, _ = run(
        label,
        ["--bed", str(p), "--max-reference-mismatch-records", "0"],
        2 if records else 0,
    )
    assert d["reference_mismatches"] == {
        "all_input_records": 2,
        "records_overlapping_scope": records,
        "primary_bases_in_scope": primary,
    }
    assert (
        all(e["contig"] == "chr1" for e in d["diagnostics"]["record_examples"])
        if records
        else not d["diagnostics"]["record_examples"]
    )
run(
    "masked-base-pass",
    ["--bed", str(study / "masked.bed"), "--max-reference-mismatch-bases", "0"],
)
d, _ = run(
    "no-examples", ["--max-examples", "0", "--max-reference-mismatch-records", "0"], 2
)
assert (
    not d["diagnostics"]["record_examples"]
    and not d["diagnostics"]["interval_examples"]
)
assert d["reference_mismatches"]["records_overlapping_scope"] == 2
d, _ = run("one-example", ["--max-examples", "1"])
assert (
    max(
        collections.Counter(
            x["reason"] for x in d["diagnostics"]["record_examples"]
        ).values()
    )
    == 1
)
run("bad-example-limit", ["--max-examples", "101"], 1)
seed = random_seed(71942)
rng = random.Random(seed)
for n in range(45):
    targets = []
    for i in range(rng.randint(3, 35)):
        c = rng.choice(list(seqs))
        s = rng.randrange(len(seqs[c]))
        e = rng.randint(s + 1, min(len(seqs[c]), s + 60))
        g = rng.choice(["GeneA", "GeneB", "GeneC"])
        ex = rng.choice(["exon1", "exon2", "exon3"])
        targets.append((c, s, e, g, ex))
    targets += targets[:2]
    if n == 0:
        targets += [
            ("chr1", 0, 70, "<script>bad</script>", 'E"&1'),
            ("chr2", 10, 20, "<script>bad</script>", 'E"&1'),
        ]
    p = study / f"targets-{n}.tsv"
    text = "contig\tstart\tend\tgene\texon\n" + "".join(
        "\t".join(map(str, t)) + "\n" for t in targets
    )
    if n % 2:
        p = p.with_suffix(".tsv.gz")
        p.write_bytes(gzip.compress(text.encode()))
    else:
        p.write_text(text)
    d, out = run(f"genes-{n}", ["--gene-targets", str(p)])
    gene_positions = collections.defaultdict(set)
    exon_positions = collections.defaultdict(set)
    union = collections.defaultdict(set)
    for c, s, e, g, ex in targets:
        gene_positions[c, g].update(range(s, e))
        exon_positions[c, g, ex].update(range(s, e))
        union[c].update(range(s, e))
    for name, expected, keyfields in [
        ("genes", gene_positions, ["contig", "gene"]),
        ("exons", exon_positions, ["contig", "gene", "exon"]),
    ]:
        assert len(d[name]) == len(expected)
        for row in d[name]:
            key = tuple(row[k] for k in keyfields)
            positions = expected[key]
            tally = collections.Counter(maps[key[0]][pos] for pos in positions)
            assert row["bases"] == len(positions) and row["unresolved_bases"] == sum(
                v for k, v in tally.items() if not k.startswith("callable_")
            )
            assert {k: v for k, v in row["bases_by_state"].items() if v} == dict(tally)
            assert set(p for a, b in row["intervals"] for p in range(a, b)) == positions
            if name == "genes":
                exons = [ps for k, ps in exon_positions.items() if k[:2] == key]
                assert row["targeted_exons"] == len(exons)
                assert row["unresolved_exons"] == sum(
                    any(not maps[key[0]][pos].startswith("callable_") for pos in ps)
                    for ps in exons
                )
        with (out / f"{name}.tsv").open() as f:
            tsv = list(csv.DictReader(f, delimiter="\t"))
            assert len(tsv) == len(d[name])
            for row, jrow in zip(tsv, d[name]):
                for k in ["bases", "callable_bases", "unresolved_bases"]:
                    assert int(row[k]) == jrow[k]
                for state, num in jrow["bases_by_state"].items():
                    assert int(row[state]) == num
    assert d["requested_bases"] == sum(map(len, union.values()))
    assert "<script>bad</script>" not in (out / "report.html").read_text()
    for example in d["diagnostics"]["interval_examples"]:
        assert all(
            p in union[example["contig"]]
            and maps[example["contig"]][p] == example["reason"]
            for p in range(example["start"], example["end"])
        )
# Error handling, CRLF and long reference alleles.
for name, text in [
    ("header", "chr1\t0\t1\tG\tE\n"),
    ("missing-name", "contig\tstart\tend\tgene\texon\nchr1\t0\t1\t.\tE\n"),
    ("empty", "contig\tstart\tend\tgene\texon\n"),
    ("outside-ref", "contig\tstart\tend\tgene\texon\nchr1\t0\t201\tG\tE\n"),
]:
    p = study / f"bad-{name}.tsv"
    p.write_text(text)
    run(f"bad-{name}", ["--gene-targets", str(p)], 1)
p = study / "crlf.tsv"
p.write_bytes(b"contig\tstart\tend\tgene\texon\r\n1\t0\t20\tG\tE\r\n")
d, _ = run("crlf", ["--gene-targets", str(p)])
assert d["genes"][0]["contig"] == "chr1"
run(
    "conflicting-targets",
    ["--bed", str(study / "masked.bed"), "--gene-targets", str(p)],
    1,
)
v = study / "long-ref.vcf"
v.write_text(
    header
    + f"chr1\t81\t.\t"
    + ("A" * 90 + "T")
    + "\tA\t.\tPASS\t.\tGT:DP:GQ\t0/1:20:30\n"
)
d, _ = run("long-mismatch", input=v)
m = d["diagnostics"]["record_examples"][0]["observed"]
assert m["first_reference_mismatch"]["position"] == 170 and len(m["REF"]) == 81
summary = {
    "checks_passed": len(results),
    "checks": results,
    "random_gene_union_comparisons": 45,
    "hand_calculated_base_labels_match": True,
    "seed": seed,
}
(study / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps({k: v for k, v in summary.items() if k != "checks"}, indent=2))
