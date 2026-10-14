"""Check that the independent benchmark oracle rejects plausible corrupted evidence.

Uses one tiny real native run, then changes its results while refreshing output hashes so
the arithmetic and attribution gates (rather than only checksums) must detect corruption.
"""

import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import robustness_scale as scale


class OracleEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.binary = Path(os.environ.get('PGSUM_BENCH_BINARY', 'target/release/pgsum')).resolve()
        if not cls.binary.is_file():
            raise unittest.SkipTest('Build the release binary, or set PGSUM_BENCH_BINARY')
        cls.temporary = tempfile.TemporaryDirectory(prefix='pgsum-scale-oracle-')
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.root = Path(cls.temporary.name)
        output = cls.root / 'source'
        subprocess.run([
            sys.executable, str(Path(scale.__file__).resolve()), '--binary', str(cls.binary),
            '--output', str(output), '--terms', '220', '--samples', '4', '--missing-rates',
            '0.005', '--threads', '1', '--keep-inputs',
        ], check=True, capture_output=True, text=True)
        cls.native = output / 'case-220-4-missing-0.005' / 'run-t1-r0'
        cls.expected = scale.expected_results(220, 4, 0.005, 42, 'site', False)
        cls.binary_hash = scale.sha256(cls.binary)

    def setUp(self):
        self.output = self.root / self.id().split('.')[-1]
        shutil.copytree(self.native, self.output)

    def verify(self):
        return scale.verify(self.output, self.expected, 220, 4, False, self.binary_hash)

    def report(self, change):
        path = self.output / 'report.json'
        value = json.loads(path.read_text())
        change(value)
        path.write_text(json.dumps(value))

    def refresh_hash(self, filename):
        self.report(lambda value: value['outputs'].update({filename: scale.sha256(self.output / filename)}))

    def test_real_evidence_passes(self):
        self.verify()

    def test_score_corruption_fails_with_refreshed_checksum(self):
        path = self.output / 'samples.tsv'
        with path.open() as source:
            reader = csv.DictReader(source, delimiter='\t')
            fields, rows = reader.fieldnames, list(reader)
        rows[0]['scenario_score'] = scale.decimal_units(scale.score_units(rows[0]['scenario_score']) + 1)
        with path.open('w') as output:
            writer = csv.DictWriter(output, fieldnames=fields, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        self.refresh_hash('samples.tsv')
        with self.assertRaises(AssertionError):
            self.verify()

    def test_mask_counts_are_not_inferred_from_nonzero_changes(self):
        def corrupt(value):
            arm = value['scores'][0]['scenarios'][1]
            arm['masked_observations'] = sum(self.expected['scenarios']['random20']['changed'])
        # Zero dosage/weight masks count as missing observations even though delta is zero.
        self.assertGreater(sum(self.expected['scenarios']['random20']['masked']),
                           sum(self.expected['scenarios']['random20']['changed']))
        self.report(corrupt)
        with self.assertRaises(AssertionError):
            self.verify()

    def test_top_ordinal_corruption_fails_with_refreshed_checksum(self):
        path = self.output / 'top-changes.jsonl'
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        next(row for row in rows if row['top'])['top'][0]['ordinal'] += 1
        path.write_text(''.join(json.dumps(row) + '\n' for row in rows))
        self.refresh_hash('top-changes.jsonl')
        with self.assertRaises(AssertionError):
            self.verify()

    def test_weight_coverage_corruption_fails(self):
        self.report(lambda value: value['scores'][0].update({'baseline_weight_fraction': 1.0}))
        self.assertLess(self.expected['baseline_weight'], self.expected['total_weight'])
        with self.assertRaises(AssertionError):
            self.verify()

    def test_metrics_corruption_fails(self):
        self.report(lambda value: value['scores'][0]['scenarios'][1]['metrics'].update({'mae': 0.0}))
        with self.assertRaises(AssertionError):
            self.verify()


if __name__ == '__main__':
    unittest.main()
