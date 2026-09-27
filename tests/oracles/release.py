import json
import subprocess
from pathlib import Path

from support import workspace

repo, study, exe, temporary = workspace()
outroot = study / "cases"
outroot.mkdir()
results = []


def run(name, opts=(), code=0, ref=None, vcf=None):
    out = outroot / name
    p = subprocess.run(
        [
            str(exe),
            "--gvcf",
            str(vcf or repo / "examples/sample.g.vcf"),
            "--reference",
            str(ref or repo / "examples/reference.fa"),
            "--bed",
            str(repo / "examples/targets.bed"),
            "--out",
            str(out),
            *opts,
        ],
        capture_output=True,
        text=True,
    )
    assert p.returncode == code, (name, p.returncode, p.stderr)
    assert out.exists() == (code in (0, 2)), name
    assert not list(outroot.glob(".gvcf-audit-*")), name
    j = json.loads((out / "report.json").read_text()) if out.exists() else None
    results.append({"name": name, "passed": True, "exit_code": code})
    return j, p


j, _ = run("demo")
assert j["callable_bases"] == 600 and j["requested_bases"] == 1000
j, _ = run(
    "quality-pass",
    ["--min-callable-percent", "60", "--max-reference-mismatch-bases", "0"],
)
assert j["quality_checks"]["passed"]
j, _ = run("quality-fail", ["--min-callable-percent", "60.001"], 2)
assert not j["quality_checks"]["passed"]
j, p = run("quiet", ["--quiet"])
assert not p.stderr
run("bad-quality-limit", ["--min-callable-percent", "101"], 1)
run("nan-quality-limit", ["--min-callable-percent", "NaN"], 1)
run("bad-cli", ["--not-a-real-option"], 1)
original = (repo / "examples/sample.g.vcf").read_text()
for name, text, code, opts in [
    (
        "missing-unused-header",
        original.replace(
            "##source=DeepVariant",
            "##source=DeepVariant\n##contig=<ID=unused,length=25>",
        ),
        0,
        [],
    ),
    ("missing-used-header", original.replace("demo\t", "unused\t"), 1, []),
    (
        "duplicate-header",
        original.replace(
            "##source=DeepVariant",
            "##source=DeepVariant\n##contig=<ID=demo,length=1000>",
        ),
        1,
        [],
    ),
    ("decimal-depth", original.replace("0/0:30:40", "0/0:30.5:40"), 1, []),
    ("scientific-depth", original.replace("0/0:30:40", "0/0:3e1:40"), 1, []),
    ("mixed-missing-alt", original.replace("\tG\t.", "\tG,.\t."), 1, []),
    (
        "reference-limit",
        original.replace("demo\t701\t.\tA", "demo\t701\t.\tT"),
        2,
        ["--max-reference-mismatch-bases", "0"],
    ),
    (
        "html-injection",
        original.replace("\tExample\n", "\t</script><img src=x onerror=alert(1)>\n"),
        0,
        [],
    ),
]:
    p = outroot / (name + ".vcf")
    p.write_text(text)
    j, _ = run(name, opts, code, vcf=p)
    if name == "missing-unused-header":
        assert j["header"]["unmapped_header_contigs"] == ["unused"] and j["warnings"]
    if name == "html-injection":
        assert "<img src=x" not in (outroot / name / "report.html").read_text()
ref = (repo / "examples/reference.fa").read_bytes()
idx = (repo / "examples/reference.fa.fai").read_text()
for name, seq, index, code in [
    ("stale-name", ref, idx.replace("demo", "wrong"), 1),
    ("stale-length", ref, idx.replace("1000", "999"), 1),
    ("stale-offset", ref, idx.replace("\t6\t", "\t7\t"), 1),
    ("stale-width", ref, idx.replace("\t51\n", "\t52\n"), 1),
    ("truncated-reference", ref[:-2], idx, 1),
    (
        "blank-between-contigs",
        ref + b"\n>unused\nAAAA\n",
        idx + f"unused\t4\t{len(ref) + 9}\t4\t5\n",
        0,
    ),
    (
        "CRLF-reference",
        ref.replace(b"\n", b"\r\n"),
        idx.replace("\t6\t50\t51", "\t7\t50\t52"),
        0,
    ),
    ("no-terminal-newline", ref[:-1], idx, 0),
]:
    p = outroot / (name + ".fa")
    p.write_bytes(seq)
    Path(str(p) + ".fai").write_text(index)
    run(name, code=code, ref=p)
for arg in ["--version", "--help"]:
    p = subprocess.run([str(exe), arg], capture_output=True, text=True)
    assert p.returncode == 0
    results.append({"name": arg, "passed": True})
(study / "release-v050-results.json").write_text(json.dumps(results, indent=2) + "\n")
print(f"{len(results)} release checks passed")
