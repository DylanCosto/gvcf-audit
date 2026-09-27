"""Run the independent synthetic checks without downloaded genomes."""

import subprocess
import sys
from pathlib import Path

if len(sys.argv) != 2:
    raise SystemExit("Usage: python3 tests/oracles/run.py /path/to/gvcf-audit")
if not __debug__:
    raise SystemExit("Run the checks without Python's -O option")
binary = Path(sys.argv[1]).resolve(strict=True)
for name in ("audit", "release", "compare", "targets", "reblocking", "genes"):
    print(f"Running {name} checks", flush=True)
    subprocess.run(
        [
            sys.executable,
            "-B",
            str(Path(__file__).with_name(name + ".py")),
            str(binary),
        ],
        check=True,
    )
