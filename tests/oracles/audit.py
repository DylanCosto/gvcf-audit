import gzip
import json
import random
import subprocess

from support import random_seed, workspace

repo, study, exe, temporary = workspace()
run_dir = study / "cases"
if run_dir.exists():
    raise RuntimeError("Use a new validation directory")
run_dir.mkdir()
seqs = {"chr1": "A" * 100 + "NNR" + "A" * 97, "chr2": "C" * 80, "chrX": "G" * 60}
ref = study / "reference.fa"
with ref.open("wb") as f, (study / "reference.fa.fai").open("w") as idx:
    for name, seq in seqs.items():
        f.write(f">{name}\n".encode())
        offset = f.tell()
        for i in range(0, len(seq), 10):
            f.write(seq[i : i + 10].encode() + b"\n")
        idx.write(f"{name}\t{len(seq)}\t{offset}\t10\t11\n")
header = (
    "##fileformat=VCFv4.2\n"
    + "".join(f"##contig=<ID={c},length={len(s)}>\n" for c, s in seqs.items())
    + '##source=DeepVariant\n##FILTER=<ID=RefCall,Description="Genotyping model thinks this site is reference.">\n'
)
columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
results = []


def run(
    label,
    records,
    options=(),
    bed=None,
    compressed=False,
    expect_error=False,
    head=None,
):
    path = run_dir / f"{label}.vcf"
    text = (head if head is not None else header + columns) + "\n".join(records) + "\n"
    if compressed:
        path = path.with_suffix(".vcf.gz")
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)
    out = run_dir / label
    cmd = [
        str(exe),
        "--gvcf",
        str(path),
        "--reference",
        str(ref),
        "--out",
        str(out),
        *options,
    ]
    if bed is not None:
        b = run_dir / f"{label}.bed"
        b.write_text("".join(f"{c}\t{s}\t{e}\t{name}\n" for c, s, e, name in bed))
        cmd += ["--bed", str(b)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if expect_error:
        assert p.returncode != 0, (label, p.stderr)
        assert not out.exists(), label
        assert not list(run_dir.glob(".gvcf-audit-*")), label
        results.append(
            {"name": label, "passed": True, "expected_failure": p.stderr.strip()}
        )
        return
    assert p.returncode == 0, (label, p.stderr)
    data = json.loads((out / "report.json").read_text())
    mapped = {c: {} for c in seqs}
    for file in ["callable.bed", "unresolved.bed"]:
        for line in (out / file).read_text().splitlines():
            c, s, e, state = line.split("\t")
            for pos in range(int(s), int(e)):
                assert pos not in mapped[c], (label, "overlapping output", c, pos)
                mapped[c][pos] = state
    assert sum(map(len, mapped.values())) == data["requested_bases"], label
    results.append({"name": label, "passed": True})
    return data, mapped, out


def record(
    c,
    s,
    e,
    gt="0/0",
    dp="20",
    gq="30",
    filt="PASS",
    alt="<NON_REF>",
    fmt="GT:MIN_DP:GQ",
    extra="",
):
    return f"{c}\t{s + 1}\t.\t{seqs[c][s]}\t{alt}\t.\t{filt}\tEND={e}\t{fmt}\t{gt}:{dp}:{gq}{extra}"


# Independent per-base oracle for interval logic, quality precedence and overlapping BED rows.
seed = random_seed(892)
rng = random.Random(seed)
for case in range(45):
    rows = []
    oracle = {c: [[] for _ in seq] for c, seq in seqs.items()}
    for c, seq in seqs.items():
        blocks = []
        for i in range(rng.randint(0, 65)):
            s = rng.randrange(len(seq))
            e = rng.randint(s + 1, min(len(seq), s + 20))
            choice = rng.randrange(6)
            gt, dp, gq, filt, state = [
                ("0/0", "20", "30", "PASS", "callable_reference"),
                ("./.", "20", "30", "PASS", "no_call"),
                ("0/0", "2", "30", "PASS", "low_depth"),
                ("0/0", "20", "3", "PASS", "low_gq"),
                ("0/0", ".", "30", "PASS", "quality_missing"),
                ("0/0", "20", "30", "LowQual", "site_filtered"),
            ][choice]
            if seq[s] not in "ACGT":
                state = "reference_ambiguous"
            blocks.append((s, record(c, s, e, gt, dp, gq, filt)))
            for p in range(s, e):
                oracle[c][p].append(state)
        rows += [line for _, line in sorted(blocks, key=lambda x: x[0])]
    bed = (
        None
        if case % 3 == 0
        else [
            (c, s, min(len(seqs[c]), s + rng.randint(1, 40)), f"region_{i}")
            for i, (c, s) in enumerate(
                (lambda c: (c, rng.randrange(len(seqs[c]))))(rng.choice(list(seqs)))
                for _ in range(12)
            )
        ]
    )
    data, mapped, out = run(f"random-{case}", rows, bed=bed, compressed=case % 2 == 0)
    expected = {c: {} for c in seqs}
    universe = {
        c: set(range(len(seq))) if bed is None else set() for c, seq in seqs.items()
    }
    if bed:
        for c, s, e, _ in bed:
            universe[c].update(range(s, e))
    for c, positions in universe.items():
        for p in positions:
            labels = oracle[c][p]
            expected[c][p] = (
                "reference_ambiguous"
                if seqs[c][p] not in "ACGT"
                else (
                    "overlapping_records"
                    if len(labels) > 1
                    else labels[0]
                    if labels
                    else "no_record"
                )
            )
    assert expected == mapped, (case, "oracle mismatch")
    for reg in data["regions"]:
        tally = {state: 0 for state in data["bases_by_state"]}
        for pos in range(reg["start"], reg["end"]):
            tally[expected[reg["contig"]][pos]] += 1
        assert tally == reg["bases_by_state"], (case, "region mismatch")

cases = [
    ("dv-refcall", record("chr1", 0, 10, filt="RefCall"), [], "callable_reference"),
    (
        "generic-refcall",
        record("chr1", 0, 10, filt="RefCall"),
        ["--caller", "generic"],
        "site_filtered",
    ),
    (
        "dragen-haploid",
        record("chrX", 0, 10, gt="0"),
        ["--caller", "dragen"],
        "callable_reference",
    ),
    ("gatk-block", record("chr1", 0, 10), ["--caller", "gatk"], "callable_reference"),
    ("missing-min-dp", record("chr1", 0, 10, fmt="GT:DP:GQ"), [], "quality_missing"),
    (
        "dp-fallback",
        record("chr1", 0, 10, fmt="GT:DP:GQ"),
        ["--allow-block-dp"],
        "callable_reference",
    ),
    (
        "unfiltered",
        record("chr1", 0, 10, filt="."),
        ["--caller", "generic"],
        "filter_not_assessed",
    ),
    (
        "dv-block-unset-filter",
        record("chr1", 0, 10, filt="."),
        [],
        "callable_reference",
    ),
    (
        "gatk-block-unset-filter",
        record("chr1", 0, 10, filt="."),
        ["--caller", "gatk"],
        "callable_reference",
    ),
    (
        "dragen-block-unset-filter",
        record("chr1", 0, 10, filt="."),
        ["--caller", "dragen"],
        "callable_reference",
    ),
    (
        "dv-low-quality-unset-filter",
        record("chr1", 0, 10, filt=".", gq="1"),
        [],
        "low_gq",
    ),
    (
        "dv-variant-unfiltered",
        "chr1\t1\t.\tA\tG\t.\t.\t.\tGT:DP:GQ\t0/1:20:30",
        [],
        "filter_not_assessed",
    ),
    (
        "allow-unfiltered",
        record("chr1", 0, 10, filt="."),
        ["--allow-unfiltered"],
        "callable_reference",
    ),
    ("partial", record("chr1", 0, 10, gt="0/."), [], "partial_no_call"),
    ("polyploid", record("chr1", 0, 10, gt="0/0/0"), [], "unsupported_ploidy"),
    ("symbolic-called", record("chr1", 0, 10, gt="0/1"), [], "unsupported_allele"),
    (
        "ft-filter",
        record("chr1", 0, 10, fmt="GT:MIN_DP:GQ:FT", extra=":Fail"),
        [],
        "genotype_filtered",
    ),
    (
        "snv",
        "chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t0/1:20:30",
        [],
        "callable_variant",
    ),
    (
        "indel",
        "chr1\t1\t.\tAA\tA\t.\tPASS\t.\tGT:DP:GQ\t0/1:20:30",
        [],
        "complex_variant",
    ),
    (
        "mismatch",
        "chr1\t1\t.\tT\tG\t.\tPASS\t.\tGT:DP:GQ\t0/1:20:30",
        [],
        "reference_mismatch",
    ),
]
for name, line, opts, state in cases:
    _, m, _ = run(name, [line], options=opts)
    c = line.split("\t")[0]
    assert m[c][0] == state, (name, m[c][0])
run("unsorted", [record("chr1", 5, 7), record("chr1", 2, 3)], expect_error=True)
run(
    "revisited-contig",
    [record("chr1", 0, 1), record("chr2", 0, 1), record("chr1", 3, 4)],
    expect_error=True,
)
run("bad-end", [record("chr1", 3, 2)], expect_error=True)
run("bad-GT", [record("chr1", 0, 1, gt="2/2")], expect_error=True)
run("nonfinite-GQ", [record("chr1", 0, 1, gq="NaN")], expect_error=True)
run("empty-bed", [], bed=[], expect_error=True)
run("unknown-bed-contig", [], bed=[("missing", 0, 1, "bad")], expect_error=True)
run(
    "header-length",
    [],
    head=(header + columns).replace("length=200", "length=201"),
    expect_error=True,
)
run(
    "multiple-samples",
    [],
    head=(header + columns).replace("\tS\n", "\tS\tT\n"),
    expect_error=True,
)
_, m, out = run(
    "select-sample",
    [record("chr1", 0, 10) + "\t./.:20:30"],
    head=(header + columns).replace("\tS\n", "\tS\tT\n"),
    options=["--sample", "T"],
)
assert m["chr1"][0] == "no_call"
_, m, out = run(
    "alias",
    [record("chr1", 0, 10).replace("chr1\t", "1\t")],
    head=(header + columns).replace("ID=chr1,", "ID=1,"),
)
assert m["chr1"][0] == "callable_reference"
_, _, out = run(
    "html-escape",
    [record("chr1", 0, 10)],
    bed=[("chr1", 0, 10, "<script>alert(1)</script>")],
)
assert "<script>alert(1)</script>" not in (out / "report.html").read_text()
p = subprocess.run(
    [
        str(exe),
        "--gvcf",
        str(run_dir / "alias.vcf"),
        "--reference",
        str(ref),
        "--out",
        str(out),
    ],
    capture_output=True,
)
assert p.returncode != 0
results.append({"name": "existing-output-refused", "passed": True})
# A truncated gzip must not leave a report behind.
pth = run_dir / "truncated.vcf.gz"
pth.write_bytes(
    gzip.compress((header + columns + record("chr1", 0, 10) + "\n").encode())[:-5]
)
p = subprocess.run(
    [
        str(exe),
        "--gvcf",
        str(pth),
        "--reference",
        str(ref),
        "--out",
        str(run_dir / "truncated"),
    ],
    capture_output=True,
)
assert p.returncode != 0 and not (run_dir / "truncated").exists()
results.append({"name": "truncated-gzip-refused", "passed": True})
(study / "baseline-v050-results.json").write_text(json.dumps(results, indent=2) + "\n")
print(
    f"{len(results)} validation checks passed, including 45 random interval comparisons."
)
