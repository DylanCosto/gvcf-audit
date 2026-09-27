"""Run the independent synthetic checks without downloaded genomes."""

import argparse
import os
import secrets
import shlex
import subprocess
import sys
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("binary", type=Path)
seeds = parser.add_mutually_exclusive_group()
seeds.add_argument(
    "--seed", type=int, help="Replay this seed in every randomized suite"
)
seeds.add_argument(
    "--random-seed", action="store_true", help="Choose and print a fresh seed"
)
args = parser.parse_args()
if not __debug__:
    raise SystemExit("Run the checks without Python's -O option")
binary = args.binary.resolve(strict=True)
seed = secrets.randbits(64) if args.random_seed else args.seed
environment = os.environ.copy()
environment.pop("GVCF_AUDIT_SEED", None)
replay = [sys.executable, str(Path(__file__).resolve()), str(binary)]
if seed is not None:
    environment["GVCF_AUDIT_SEED"] = str(seed)
    replay.extend(["--seed", str(seed)])
print(f"Seed: {seed if seed is not None else 'fixed suite defaults'}", flush=True)
print(f"Replay: {shlex.join(replay)}", flush=True)
for name in ("audit", "release", "compare", "targets", "reblocking", "genes", "cohort"):
    print(f"Running {name} checks", flush=True)
    result = subprocess.run(
        [
            sys.executable,
            "-B",
            str(Path(__file__).with_name(name + ".py")),
            str(binary),
        ],
        env=environment,
    )
    if result.returncode:
        print(f"Failed suite: {name}. Replay: {shlex.join(replay)}", file=sys.stderr)
        raise SystemExit(1)
