import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from check_group_summary import check


class GroupCheckerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.native = self.root/'native'; self.native.mkdir()
        self.out = self.root/'grouped'; self.out.mkdir()
        self.mapping = self.root/'groups.tsv'
        self.mapping.write_text('sample\tgroup\na\tfirst\nb\tfirst\nc\tsecond\nd\tsecond\n')
        (self.native/'samples.tsv').write_text('pgs_id\tscenario\tsample\tbaseline_score\tscenario_score\n'+
            ''.join(f'S\tx\t{name}\t{base}\t{base+10}\n' for name,base in zip('abcd',[0,10,20,30])))
        (self.native/'references.tsv').write_text('pgs_id\tscenario\tsample\treference\n')
        identity = dict(pgs_id='S',baseline_terms=2,baseline_term_sha256='fixture',total_model_terms=2,excluded_terms={})
        native = dict(experiment=dict(cutoffs=[50],reference_sets=[]),scores=[identity],outputs={})
        (self.native/'report.json').write_text(json.dumps(native))
        rows = []
        for group, shift, up in [('first',25.0,1),('second',18.75,0)]:
            rows.append(dict(scenario='x',group=group,metrics=dict(n=2,mae=10.0,rmse=10.0,max_abs_change=10.0,
                pearson=1.0,spearman=1.0,mean_abs_percentile_shift=shift,max_abs_percentile_shift=25.0,
                crossings=[dict(cutoff=50,up=up,down=0,unchanged=2-up,evaluable=2)],unavailable_reasons={}),reference_comparisons=[]))
        self.summary = dict(schema='pgsum-robustness-groups-v1',percentile_basis='full_target_baseline',correlation_basis='within_group',
            groups=[dict(id='first',samples=2,members=['a','b']),dict(id='second',samples=2,members=['c','d'])],
            scores=[dict(**identity,source_samples=4,scenarios=rows)],outputs={})

        sources = [self.mapping, self.native/'report.json', self.native/'samples.tsv', self.native/'references.tsv']
        self.summary['source_sha256'] = {str(path): 'sha256:'+hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}

    def run_check(self, report=None):
        (self.out/'groups.json').write_text(json.dumps(report or self.summary))
        return check(self.native,self.mapping,self.out)

    def test_hand_calculated_global_thresholds(self):
        result = self.run_check()
        self.assertEqual(result['score_metric_checks'],14)
        self.assertEqual(result['cutoff_checks'],2)

    def test_group_local_percentiles_are_rejected(self):
        changed=copy.deepcopy(self.summary)
        changed['scores'][0]['scenarios'][0]['metrics']['mean_abs_percentile_shift']=37.5
        with self.assertRaises(ValueError): self.run_check(changed)

    def test_literal_quoted_labels_are_preserved(self):
        self.mapping.write_text(self.mapping.read_text().replace('first', '"first"'))
        self.summary['source_sha256'][str(self.mapping)] = 'sha256:'+hashlib.sha256(self.mapping.read_bytes()).hexdigest()
        self.summary['groups'][0]['id'] = '"first"'
        self.summary['scores'][0]['scenarios'][0]['group'] = '"first"'
        self.assertTrue(self.run_check()['passed'])

    def test_changed_source_is_rejected(self):
        self.mapping.write_text(self.mapping.read_text().replace('first', 'renamed'))
        with self.assertRaisesRegex(ValueError, 'source fingerprint'):
            self.run_check()

    def test_membership_and_denominator_changes_are_rejected(self):
        for mutate in [lambda s:s['groups'][0].update(members=['a','c']),
                       lambda s:s['scores'][0]['scenarios'][0]['metrics'].update(n=3),
                       lambda s:s['scores'][0]['scenarios'][0]['metrics']['crossings'][0].update(up=0)]:
            changed=copy.deepcopy(self.summary);mutate(changed)
            with self.assertRaises(ValueError): self.run_check(changed)

    def test_changed_payload_is_rejected(self):
        file=self.out/'metrics.tsv';file.write_text('original')
        self.summary['outputs']['metrics.tsv']='sha256:'+hashlib.sha256(file.read_bytes()).hexdigest()
        file.write_text('changed')
        with self.assertRaises(ValueError): self.run_check()


if __name__ == '__main__':
    unittest.main()
