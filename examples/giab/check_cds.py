"""Independent per-base check of the complete autosomal CDS overlap."""

from collections import defaultdict, Counter
import gzip
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("data_dir", type=Path)
parser.add_argument("audit_dir", type=Path)
parser.add_argument("overlap_json", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
root, audit = args.data_dir.resolve(), args.audit_dir.resolve()
report = json.loads((audit / "report.json").read_text())
lengths = {r["contig"]: r["requested_bases"] for r in report["contigs"]}
cds = defaultdict(list)
with gzip.open(root / "GRCh38_refseq_cds.bed.gz", "rt") as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        contig, start, end, *_ = line.split()
        if contig in {f"chr{i}" for i in range(1, 23)}:
            cds[contig].append((int(start), int(end)))
counts = Counter()
for filename in ["callable.bed", "unresolved.bed"]:
    current = None
    mask = bytearray()
    with (audit / filename).open() as handle:
        for line in handle:
            contig, start, end, state = line.split()
            if contig not in cds:
                continue
            if contig != current:
                current = contig
                mask = bytearray(lengths[contig])
                for a, b in cds[contig]:
                    mask[a:b] = b"\1" * (b-a)
            count = mask[int(start):int(end)].count(1)
            if count:
                counts[state] += count
expected = json.loads(args.overlap_json.read_text())
assert report["header"]["sample"] == expected["sample"]
assert dict(counts) == expected["results"]["refseq_cds"]["bases_by_state"], counts
result = {"method": "Independent per-base byte mask from raw CDS rows, without interval merging or binary searches",
       "scope": "All autosomal RefSeq CDS bases", "bases_by_state": dict(counts), "match": True}
with args.output.open("x") as handle:
    handle.write(json.dumps(result, indent=2) + "\n")
print("All autosomal CDS state counts match the independent per-base calculation.")
