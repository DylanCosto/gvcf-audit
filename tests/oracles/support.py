"""Shared setup for the portable, standard-library-only checks."""

import sys
import tempfile
from pathlib import Path


def workspace():
    if not __debug__:
        raise RuntimeError("Run the checks without Python's -O option")
    if len(sys.argv) != 2:
        raise SystemExit("Pass the path to the gvcf-audit binary")
    binary = Path(sys.argv[1]).resolve(strict=True)
    repo = Path(__file__).resolve().parents[2]
    temporary = tempfile.TemporaryDirectory(prefix="gvcf-audit-checks-")
    return repo, Path(temporary.name), binary, temporary
