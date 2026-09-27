"""Summarize an autosomal HG002/HG003 audit against the GIAB example regions."""

import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
import gzip
from heapq import merge as merge_sorted
import json
from pathlib import Path


def merge(rows):
    result = []
    for start, end in sorted(rows):
        if not 0 <= start < end:
            raise ValueError("Invalid BED coordinates")
        if result and start <= result[-1][1]:
            result[-1] = result[-1][0], max(end, result[-1][1])
        else:
            result.append((start, end))
    return result


def read_bed(path, contigs):
    rows = defaultdict(list)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            contig, start, end, *_ = line.split()
            if contig in contigs:
                rows[contig].append((int(start), int(end)))
    return {c: merge(v) for c, v in rows.items()}


def intersect(a, b):
    out = {}
    for contig, left in a.items():
        right = b.get(contig, [])
        rows = []
        i = j = 0
        while i < len(left) and j < len(right):
            lo = max(left[i][0], right[j][0])
            hi = min(left[i][1], right[j][1])
            if lo < hi:
                rows.append((lo, hi))
            if left[i][1] <= right[j][1]:
                i += 1
            else:
                j += 1
        out[contig] = rows
    return out


def audit_rows(path, ranks, allowed_states):
    previous = (-1, -1)
    with path.open() as handle:
        for line in handle:
            contig, start, end, state = line.split()
            start, end = int(start), int(end)
            key = ranks[contig], start
            if key < previous or start >= end or state not in allowed_states:
                raise ValueError(f"Invalid or unsorted audit BED: {path}")
            previous = key
            yield key[0], start, end, state, contig


def summarize(data_dir, audit_dir, sample):
    report = json.loads((audit_dir / "report.json").read_text())
    if report['status'] != 'complete' or report['header']['sample'] != sample:
        raise ValueError("Expected a complete audit for the selected sample")
    contigs = {f"chr{i}" for i in range(1, 23)}
    scope = read_bed(data_dir / "autosomes.bed", contigs)
    if set(scope) != contigs or any(len(v) != 1 or v[0][0] != 0 for v in scope.values()):
        raise ValueError("autosomes.bed must contain full chr1-22 spans")
    actual = defaultdict(list)
    for row in report['regions']:
        actual[row['contig']].append((row['start'], row['end']))
    if {c: merge(v) for c, v in actual.items()} != scope:
        raise ValueError("The audit must cover exactly the autosomal scope")
    sets = {"autosomes": scope}
    files = {
        "giab_high_confidence": f"{sample}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
        "refseq_cds": "GRCh38_refseq_cds.bed.gz",
        "low_mappability": "GRCh38_lowmappabilityall.bed.gz",
        "segmental_duplications": "GRCh38_segdups.bed.gz",
        "difficult_regions": "GRCh38_alldifficultregions.bed.gz",
    }
    for name, filename in files.items():
        sets[name] = read_bed(data_dir / filename, contigs)
        if intersect(sets[name], scope) != sets[name]:
            raise ValueError(f"Regions extend beyond the reference scope: {filename}")
    for name in files:
        if name != "refseq_cds":
            sets["cds_" + name] = intersect(sets["refseq_cds"], sets[name])
    counts = {name: Counter() for name in sets}
    indexed = defaultdict(list)
    for name, regions in sets.items():
        for contig, rows in regions.items():
            indexed[contig].append((name, rows, [end for _, end in rows]))
    ranks = {f"chr{i}": i for i in range(1, 23)}
    states = set(report['bases_by_state'])
    called_states = {'callable_reference', 'callable_variant'}
    streams = [audit_rows(audit_dir / 'callable.bed', ranks, called_states),
               audit_rows(audit_dir / 'unresolved.bed', ranks, states - called_states)]
    totals = Counter()
    ends = defaultdict(int)
    for n, (_, start, end, state, contig) in enumerate(merge_sorted(*streams), 1):
        if start != ends[contig] or end > scope[contig][0][1]:
            raise ValueError(f"Audit BED gap, overlap or out-of-scope interval on {contig}")
        ends[contig] = end
        totals[state] += end - start
        for name, rows, boundaries in indexed[contig]:
            i = bisect_right(boundaries, start)
            while i < len(rows) and rows[i][0] < end:
                a, b = rows[i]
                counts[name][state] += min(end, b) - max(start, a)
                i += 1
        if n % 5_000_000 == 0:
            print(f"Checked {n:,} audit intervals", flush=True)
    if dict(ends) != {c: v[0][1] for c, v in scope.items()}:
        raise ValueError("Audit BEDs do not cover the scope")
    if totals != Counter(report['bases_by_state']):
        raise ValueError("Audit BED counts disagree with the report")
    results = {}
    for name, regions in sets.items():
        total = sum(b-a for rows in regions.values() for a,b in rows)
        if not total or total != sum(counts[name].values()):
            raise ValueError(f"Incomplete region accounting: {name}")
        called = sum(counts[name][s] for s in called_states)
        results[name] = {"requested_bases": total, "callable_bases": called,
                         "callable_percent": 100 * called / total,
                         "bases_by_state": dict(sorted(counts[name].items()))}
    return {
        "sample": sample, "audit_version": report['version'],
        "caller": "DeepTrio 1.10.0 WGS" if sample == 'HG002' else "DeepVariant 1.10.0 WGS",
        "scope": "GRCh38 chr1-22", "stratifications_version": "3.6", "benchmark_version": "4.2.1",
        "policy": {k: report['options'][k] for k in ['min_dp', 'min_gq', 'allow_block_dp', 'allow_unfiltered']},
        "note": "Base-level overlap, not variant precision or recall. Region categories overlap.",
        "results": results,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data_dir', type=Path)
    parser.add_argument('audit_dir', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--sample', required=True, choices=['HG002', 'HG003'])
    args = parser.parse_args()
    result = summarize(args.data_dir, args.audit_dir, args.sample)
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
