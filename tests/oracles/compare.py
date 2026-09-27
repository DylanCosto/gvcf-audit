import csv
import json
import random
import shutil
import subprocess
from collections import Counter, defaultdict

from support import workspace

repo, study, exe, temporary = workspace()
cases = study / "cases"
cases.mkdir()
results = []


def command(args, code=0):
    p = subprocess.run([str(exe), *map(str, args)], capture_output=True, text=True)
    assert p.returncode == code, (args, p.returncode, p.stderr)
    return p


ref = cases / "reference.fa"
ref.write_text(">chr1\n" + "A" * 100 + "\n>chr2\n" + "A" * 80 + "\n")
ref.with_suffix(".fa.fai").write_text("chr1\t100\t6\t100\t101\nchr2\t80\t113\t80\t81\n")
header = '##fileformat=VCFv4.2\n##source=DeepVariant\n##contig=<ID=chr1,length=100>\n##contig=<ID=chr2,length=80>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
rng = random.Random(40521)


def audit(name, quality, targets, reverse=False, extra=()):
    vcf = cases / f"{name}.vcf"
    lines = []
    for c in ["chr2", "chr1"] if reverse else ["chr1", "chr2"]:
        for pos, values in enumerate(quality[c]):
            gt, dp, gq = values
            lines.append(
                f"{c}\t{pos + 1}\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t{gt}:{dp}:{gq}"
            )
    vcf.write_text(header + "\n".join(lines) + "\n")
    tsv = cases / f"{name}.tsv"
    tsv.write_text(
        "contig\tstart\tend\tgene\texon\n"
        + "".join("\t".join(map(str, r)) + "\n" for r in targets)
    )
    out = cases / name
    command(
        [
            "--gvcf",
            vcf,
            "--reference",
            ref,
            "--gene-targets",
            tsv,
            "--out",
            out,
            *extra,
        ],
        2 if "--min-callable-percent" in extra else 0,
    )
    return out


def base_map(folder):
    out = {}
    for name in ["callable.bed", "unresolved.bed"]:
        for line in (folder / name).read_text().splitlines():
            c, s, e, state = line.split("\t")
            for pos in range(int(s), int(e)):
                assert (c, pos) not in out
                out[c, pos] = state
    return out


def iscall(s):
    return s.startswith("callable_")


def compare(name, before, after, opts=(), code=0):
    out = cases / name
    p = command(
        ["compare", "--before", before, "--after", after, "--out", out, *opts], code
    )
    assert out.exists() == (code == 0)
    assert not list(cases.glob(".gvcf-audit-*"))
    results.append({"name": name, "passed": True, "exit_code": code})
    return (json.loads((out / "comparison.json").read_text()), out) if code == 0 else p


patterns = [
    ("0/0", 20, 30),
    ("0/1", 20, 30),
    ("0/0", 2, 30),
    ("0/0", 20, 2),
    ("./.", 20, 30),
    ("0/0", ".", 30),
]
first = None
for trial in range(30):
    targets = []
    for i in range(30):
        c = rng.choice(["chr1", "chr2"])
        s = rng.randrange(70)
        e = s + rng.randrange(1, 10)
        targets.append((c, s, e, "G" + str(i % 4), "E" + str(i % 8)))
    targets += targets[:3]
    old = {
        c: [rng.choice(patterns) for _ in range(n)]
        for c, n in [("chr1", 100), ("chr2", 80)]
    }
    new = {
        c: [rng.choice(patterns) for _ in range(n)]
        for c, n in [("chr1", 100), ("chr2", 80)]
    }
    a = audit(f"before-{trial}", old, targets)
    b = audit(f"after-{trial}", new, list(reversed(targets)), reverse=True)
    d, out = compare(f"random-{trial}", a, b)
    left = base_map(a)
    right = base_map(b)
    assert left.keys() == right.keys()
    gains = {p for p in left if not iscall(left[p]) and iscall(right[p])}
    losses = {p for p in left if iscall(left[p]) and not iscall(right[p])}
    assert d["summary"]["gained_callable_bases"] == len(gains) and d["summary"][
        "lost_callable_bases"
    ] == len(losses)
    expected = Counter((left[p], right[p]) for p in left)
    assert {(v["before"], v["after"]): v["bases"] for v in d["transitions"]} == expected
    changes = {}
    for line in (out / "changes.bed").read_text().splitlines():
        c, s, e, label = line.split("\t")
        x, y = label.split("->")
        assert x != y
        for pos in range(int(s), int(e)):
            assert (c, pos) not in changes
            changes[c, pos] = (x, y)
    assert changes == {p: (left[p], right[p]) for p in left if left[p] != right[p]}
    groups = defaultdict(set)
    for c, s, e, g, x in targets:
        for key in [(c, g, None), (c, g, x)]:
            groups[key].update((c, p) for p in range(s, e))
    for kind in ["genes", "exons"]:
        table = list(csv.DictReader((out / f"{kind}.tsv").open(), delimiter="\t"))
        for row, tsv in zip(d[kind], table, strict=True):
            positions = groups[row["contig"], row["gene"], row.get("exon")]
            assert row["bases"] == len(positions)
            assert row["gained_callable_bases"] == len(gains & positions) and row[
                "lost_callable_bases"
            ] == len(losses & positions)
            assert row["net_callable_bases"] == len(gains & positions) - len(
                losses & positions
            )
            for k in ["before", "after"]:
                counts = Counter(
                    (left if k == "before" else right)[p] for p in positions
                )
                assert {
                    s: n for s, n in row[k + "_bases_by_state"].items() if n
                } == counts
            for key in [
                "gained_callable_bases",
                "lost_callable_bases",
                "net_callable_bases",
            ]:
                assert int(tsv[key]) == row[key]
    if first is None:
        first = (a, b)
a, b = first
self_result, self_out = compare("self", a, a)
assert (
    self_result["summary"]["gained_callable_bases"] == 0
    and self_result["summary"]["lost_callable_bases"] == 0
    and not (self_out / "changes.bed").read_bytes()
)
reverse, _ = compare("reverse", b, a)
original = json.loads((cases / "random-0/comparison.json").read_text())
assert (
    reverse["summary"]["gained_callable_bases"]
    == original["summary"]["lost_callable_bases"]
)


def altered(name, edit, base=a):
    p = cases / ("input-" + name)
    shutil.copytree(base, p)
    j = json.loads((p / "report.json").read_text())
    edit(j, p)
    (p / "report.json").write_text(json.dumps(j))
    return p


def edit_at(path, value):
    def edit(j, p):
        obj = j
        for key in path[:-1]:
            obj = obj[key]
        obj[path[-1]] = value

    return edit


for name, path, value, flag in [
    ("policy", ["options", "min_gq"], 30, "--allow-policy-change"),
    ("sample", ["header", "sample"], "OTHER", "--allow-different-samples"),
    (
        "reference",
        ["inputs", "reference", "path"],
        "/moved/reference.fa",
        "--assume-same-reference",
    ),
    ("caller", ["header", "caller"], "gatk", "--allow-policy-change"),
]:
    changed = altered(name, edit_at(path, value))
    compare(name + "-refused", a, changed, code=1)
    d, _ = compare(name + "-allowed", a, changed, [flag])
    assert d["warnings"]
for name, path, value in [
    ("incomplete", ["status"], "partial"),
    ("schema", ["schema"], "unknown"),
    ("unknown-policy", ["policy", "id"], "unknown"),
    ("global-count", ["callable_bases"], 99999),
    ("scope", ["regions", 0, "end"], 99),
    ("group-id", ["genes", 0, "gene"], "different"),
    ("group-interval", ["exons", 0, "intervals", 0, 1], 99),
    ("contig-count", ["contigs", 0, "callable_bases"], 9999),
    ("missing-groups", ["exons"], []),
]:
    changed = altered(name, edit_at(path, value))
    compare(name, a, changed, code=1)
for name, edit in [
    (
        "bed-gap",
        lambda j, p: (p / "callable.bed").write_text(
            "\n".join((p / "callable.bed").read_text().splitlines()[1:]) + "\n"
        ),
    ),
    (
        "bed-overlap",
        lambda j, p: (p / "callable.bed").write_text(
            (p / "callable.bed").read_text() * 2
        ),
    ),
    (
        "bed-state",
        lambda j, p: (p / "callable.bed").write_text(
            (p / "callable.bed").read_text().replace("callable_reference", "low_depth")
        ),
    ),
    ("missing-bed", lambda j, p: (p / "unresolved.bed").unlink()),
    (
        "unknown-contig",
        lambda j, p: (p / "callable.bed").write_text(
            (p / "callable.bed").read_text().replace("chr1", "outside")
        ),
    ),
    ("truncated-json", None),
]:
    if edit is None:
        changed = cases / "input-truncated-json"
        shutil.copytree(a, changed)
        (changed / "report.json").write_text("{")
    else:
        changed = altered(name, edit)
    compare(name, a, changed, code=1)
# Reclassifying BED intervals must be rejected even when they remain structurally valid.
changed = altered(
    "bed-counts",
    lambda j, p: (p / "callable.bed").write_text(
        (p / "callable.bed")
        .read_text()
        .replace("callable_reference", "callable_variant")
    ),
)
compare("bed-counts", a, changed, code=1)
failed = altered("failed-gate", edit_at(["quality_checks", "passed"], False))
d, _ = compare("failed-gate", a, failed)
assert any("failed its quality" in w for w in d["warnings"])
# Equal totals can hide a complete replacement of callable positions.
targets = [("chr1", 0, 20, "SWAP", "EXON")]
q1 = {
    "chr1": [("0/0", 20 if i < 10 else 1, 30) for i in range(100)],
    "chr2": [patterns[0]] * 80,
}
q2 = {
    "chr1": [("0/0", 1 if i < 10 else 20, 30) for i in range(100)],
    "chr2": [patterns[0]] * 80,
}
a1 = audit("swap-before", q1, targets)
b1 = audit("swap-after", q2, targets)
d, _ = compare("equal-total-swap", a1, b1)
assert (
    d["summary"]["gained_callable_bases"] == 10
    and d["summary"]["lost_callable_bases"] == 10
    and d["summary"]["net_callable_bases"] == 0
)
# HTML must treat sample identifiers as text, even with an intentional sample override.
injected = altered(
    "escaping", edit_at(["header", "sample"], "</script><img src=x onerror=alert(1)>")
)
d, out = compare("escaping", a, injected, ["--allow-different-samples"])
assert "<img src=x" not in (out / "report.html").read_text()
command(["compare", "--before", a, "--after", b, "--out", cases / "self"], 1)
assert not list(cases.glob(".gvcf-audit-*"))
report = {
    "checks_passed": len(results),
    "random_per_base_comparisons": 30,
    "seed": 40521,
    "checks": results,
}
(study / "controlled-results.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps({k: v for k, v in report.items() if k != "checks"}, indent=2))
