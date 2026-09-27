"""Run with a gvcf-audit binary and the MultiQC executable (tested with 1.35)."""

import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


binary = Path(sys.argv[1]).resolve()
multiqc = Path(sys.argv[2]).resolve()
fixture = Path(__file__).resolve().parents[1] / "fixtures" / "core"
cases = [
    ("001", [], 0, 50, 20, 46, 2, "not requested"),
    ("#sample", ["--min-callable-percent", "21", "--quiet"], 2, 50, 20, 46, 2, "failed"),
    ('sample "quoted"', ["--min-callable-percent", "20"], 0, 50, 20, 46, 2, "passed"),
    ("allowed", ["--allow-unfiltered"], 0, 50, 22, 46, 0, "not requested"),
]

with tempfile.TemporaryDirectory(prefix="gvcf-audit-multiqc-") as tmp:
    root = Path(tmp)
    for i, (sample, args, code, *_) in enumerate(cases):
        work = root / "audits" / str(i)
        shutil.copytree(fixture, work)
        vcf = work / "sample.g.vcf"
        vcf.write_text(vcf.read_text().replace("\tS\n", f"\t{sample}\n"))
        result = subprocess.run(
            [str(binary), "--gvcf", str(vcf), "--reference", str(work / "reference.fa"),
             "--out", str(work / "audit"), *args],
            capture_output=True, text=True,
        )
        assert result.returncode == code, result.stderr
    subprocess.run(
        [str(multiqc), "--strict", "--no-version-check", "--no-ai", "--data-format", "json",
         str(root / "audits"), "--outdir", str(root / "multiqc")],
        check=True,
    )
    data = json.loads((root / "multiqc/multiqc_data/multiqc_data.json").read_text())
    rows = data["report_general_stats_data"]["custom_content"]
    assert set(rows) == {case[0] for case in cases}, rows
    for sample, args, _, total, callable_pct, missing_pct, filter_pct, gate in cases:
        assert rows[sample] == {
            "requested_bases": total,
            "callable_percent": callable_pct,
            "no_record_percent": missing_pct,
            "filter_not_assessed_percent": filter_pct,
            "quality_gate": gate,
            "min_dp": 10,
            "min_gq": 20,
            "caller": "gatk",
            "allow_unfiltered": "true" if "--allow-unfiltered" in args else "false",
            "allow_block_dp": "false",
        }, rows[sample]
    assert (root / "multiqc/multiqc_report.html").is_file()
print("MultiQC parsed all four audits with the expected names, counts and policies.")
