"""Small independent checks for the published GIAB overlap example."""

import gzip
import json
from pathlib import Path
import tempfile
import subprocess
import sys
import unittest

from overlap import intersect, merge, summarize


class OverlapTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.data = self.root / 'data'
        self.audit = self.root / 'audit'
        self.data.mkdir()
        self.audit.mkdir()
        self.contigs = [f'chr{i}' for i in range(1, 23)]
        (self.data / 'autosomes.bed').write_text(''.join(f'{c}\t0\t100\n' for c in self.contigs))
        specs = {
            'GRCh38_refseq_cds.bed.gz': [(10, 20), (15, 25), (40, 70)],
            'GRCh38_lowmappabilityall.bed.gz': [(15, 55)],
            'GRCh38_segdups.bed.gz': [(40, 80)],
            'GRCh38_alldifficultregions.bed.gz': [(0, 15), (50, 80)],
            'HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed': [(0, 60)],
        }
        for name, rows in specs.items():
            opener = gzip.open if name.endswith('.gz') else open
            with opener(self.data / name, 'wt') as f:
                f.write(''.join(f'{c}\t{a}\t{b}\n' for c in self.contigs for a, b in rows))
        self.report = {'status': 'complete', 'header': {'sample': 'HG003'}, 'version': 'fixture',
                       'options': dict(min_dp=10, min_gq=20, allow_block_dp=False, allow_unfiltered=False),
                       'regions': [dict(contig=c, start=0, end=100) for c in self.contigs],
                       'contigs': [dict(contig=c, requested_bases=100) for c in self.contigs],
                       'bases_by_state': {'callable_reference': 1100, 'low_depth': 1100}}
        self.write_report()
        (self.audit / 'callable.bed').write_text(''.join(f'{c}\t0\t50\tcallable_reference\n' for c in self.contigs))
        (self.audit / 'unresolved.bed').write_text(''.join(f'{c}\t50\t100\tlow_depth\n' for c in self.contigs))

    def write_report(self):
        (self.audit / 'report.json').write_text(json.dumps(self.report))

    def test_hand_calculated_intersections(self):
        results = summarize(self.data, self.audit, 'HG003')['results']
        expected = {'refseq_cds': (45, 25), 'cds_giab_high_confidence': (35, 25),
                    'cds_low_mappability': (25, 20), 'cds_segmental_duplications': (30, 10),
                    'cds_difficult_regions': (25, 5)}
        for name, (total, called) in expected.items():
            self.assertEqual(results[name]['requested_bases'], 22 * total)
            self.assertEqual(results[name]['callable_bases'], 22 * called)

    def test_independent_checker(self):
        result = self.root / 'overlap.json'
        result.write_text(json.dumps(summarize(self.data, self.audit, 'HG003')))
        checked = self.root / 'checked.json'
        script = Path(__file__).with_name('check_cds.py')
        subprocess.run([sys.executable, str(script), str(self.data), str(self.audit),
                        str(result), str(checked)], check=True, capture_output=True)
        self.assertTrue(json.loads(checked.read_text())['match'])

    def test_gap_and_overlap_rejected(self):
        path = self.audit / 'unresolved.bed'
        original = path.read_text()
        for start in [49, 51]:
            path.write_text(original.replace('chr1\t50\t100', f'chr1\t{start}\t100'))
            with self.assertRaisesRegex(ValueError, 'gap, overlap'):
                summarize(self.data, self.audit, 'HG003')

    def test_sample_and_scope_rejected(self):
        with self.assertRaisesRegex(ValueError, 'selected sample'):
            summarize(self.data, self.audit, 'HG002')
        self.report['regions'][0]['start'] = 1
        self.write_report()
        with self.assertRaisesRegex(ValueError, 'exactly the autosomal scope'):
            summarize(self.data, self.audit, 'HG003')

    def test_interval_operations(self):
        self.assertEqual(merge([(4, 8), (1, 5), (8, 9), (11, 12)]), [(1, 9), (11, 12)])
        self.assertEqual(intersect({'c': [(1, 9), (11, 12)]}, {'c': [(0, 2), (3, 4), (6, 11)]}),
                         {'c': [(1, 2), (3, 4), (6, 9)]})


if __name__ == '__main__':
    unittest.main()
