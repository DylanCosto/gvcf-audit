import json
import subprocess

from support import workspace

repo, study, exe, temporary = workspace()
s = study
(s / "ref.fa").write_text(">chr1\n" + "A" * 100 + "\n")
(s / "ref.fa.fai").write_text("chr1\t100\t6\t100\t101\n")
header = '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##GATKCommandLine=<ID=ReblockGVCF,Version="4.2.2.0">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
(s / "reblock.vcf").write_text(
    header
    + "chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=20\tGT:DP:GQ\t0/0:30:30\nchr1\t21\t.\tA\t<NON_REF>\t.\t.\tEND=40\tGT:DP:GQ\t0/0:30:0\nchr1\t41\t.\tA\t<NON_REF>\t.\t.\tEND=60\tGT:DP:GQ\t0/0:.:30\nchr1\t61\t.\tA\tC,<NON_REF>\t.\t.\t.\tGT:DP:GQ\t0/1:30:30\n"
)
cases = [
    ("default", [], 0, 60, 0, 1),
    ("dp", ["--allow-block-dp"], 20, 20, 20, 1),
    ("both", ["--allow-block-dp", "--allow-unfiltered"], 21, 20, 20, 0),
]
results = []
for name, flags, callable, missing, low, unfiltered in cases:
    out = s / name
    p = subprocess.run(
        [
            str(exe),
            "--gvcf",
            str(s / "reblock.vcf"),
            "--reference",
            str(s / "ref.fa"),
            "--out",
            str(out),
            *flags,
        ],
        capture_output=True,
        text=True,
    )
    assert p.returncode == 0, p.stderr
    d = json.loads((out / "report.json").read_text())
    b = d["bases_by_state"]
    assert d["header"]["caller"] == "gatk" and d["header"]["reblocked"]
    assert (
        d["callable_bases"] == callable
        and b["quality_missing"] == missing
        and b["low_gq"] == low
        and b["filter_not_assessed"] == unfiltered
    ), (name, b)
    assert d["adaptations"]["reference_blocks_missing_MIN_DP"] == 3
    assert any("ReblockGVCF" in w for w in d["warnings"])
    results.append(
        {
            "case": name,
            "passed": True,
            "callable_bases": callable,
            "missing_quality_bases": missing,
            "low_gq_bases": low,
        }
    )
(study / "results.json").write_text(json.dumps(results, indent=2) + "\n")
print("3 reblocking policy checks passed")
