"""Hand-calculated and adversarial checks for the independent study verifier."""
import copy
import csv
import tempfile
import unittest
from pathlib import Path

from check_study_summary import percentile_evidence, check_percentile_summary


class StudyVerifierTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.grouped = {}
        self.reference_rows = []
        scenarios = []
        for name, scores in [('baseline', [0, 1, 2]), ('reverse', [2, 1, 0])]:
            self.grouped['S', name] = [dict(sample=str(i), baseline_score=str(i), scenario_score=str(v)) for i, v in enumerate(scores)]
            scenarios.append(dict(id=name, metrics=dict(crossings=[dict(cutoff=50, up=int(name=='reverse'), down=int(name=='reverse'), unchanged=1 if name=='reverse' else 3, evaluable=3)]),
                                  reference_comparisons=[dict(id='R', reference_samples=4, shared_sample_ids=0, evaluable=3, mean_baseline_percentile=50,
                                                             mean_scenario_percentile=50, mean_absolute_shift=160/3 if name=='reverse' else 0,
                                                             max_absolute_shift=80 if name=='reverse' else 0)]))
            for i, v in enumerate(scores):
                a, b = 10+40*i, 10+40*v
                self.reference_rows.append(['S', name, str(i), 'R', 4, 0, a, b, b-a])
        self.native = dict(experiment=dict(cutoffs=[50], reference_sets=[dict(id='R')]), scores=[dict(pgs_id='S', scenarios=scenarios)])

    def write_references(self):
        with (self.directory/'references.tsv').open('w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['pgs_id','scenario','sample','reference','reference_samples','shared_sample_ids','baseline_percentile','scenario_percentile','percentile_delta'])
            writer.writerows(self.reference_rows)

    def test_hand_calculated_ties_boundary_and_frozen_reference(self):
        self.write_references()
        expected, count = percentile_evidence(self.grouped, self.native, self.directory)
        self.assertEqual(count, 6)
        reverse = expected['S', 'reverse']
        # The middle participant is exactly at the 50th percentile and remains a member.
        self.assertEqual(reverse['cutoff',50,None,'changed_count'], (2,3))
        self.assertEqual(reverse['cutoff',50,None,'changed_fraction'], (2/3,3))
        self.assertEqual(reverse['reference',None,'R','mean_absolute_shift'], (160/3,3))

    def test_undefined_percentiles_are_unavailable(self):
        for rows in self.grouped.values():
            for r in rows:
                r['baseline_score'] = '1'
        for r in self.reference_rows:
            r[6:] = ['', '', '']
        for s in self.native['scores'][0]['scenarios']:
            s['metrics']['crossings'] = [dict(cutoff=50, up=0, down=0, unchanged=0, evaluable=0)]
            r = s['reference_comparisons'][0]
            r.update(evaluable=0, mean_baseline_percentile=None, mean_scenario_percentile=None, mean_absolute_shift=None, max_absolute_shift=None)
        self.write_references()
        expected, _ = percentile_evidence(self.grouped, self.native, self.directory)
        self.assertTrue(all(v==(None,0) for result in expected.values() for v in result.values()))

    def test_corrupted_reference_or_native_counts_rejected(self):
        saved = copy.deepcopy(self.reference_rows)
        for mode in ['duplicate', 'missing', 'delta', 'baseline', 'native_count']:
            with self.subTest(mode=mode):
                self.reference_rows = copy.deepcopy(saved)
                native = copy.deepcopy(self.native)
                if mode == 'duplicate': self.reference_rows.append(self.reference_rows[0])
                if mode == 'missing': self.reference_rows.pop()
                if mode == 'delta': self.reference_rows[-1][-1] = 99
                if mode == 'baseline': self.reference_rows[-1][-3] = 80; self.reference_rows[-1][-1] = -70
                if mode == 'native_count': native['scores'][0]['scenarios'][0]['metrics']['crossings'][0]['unchanged'] = 2
                self.write_references()
                with self.assertRaises(ValueError):
                    percentile_evidence(self.grouped, native, self.directory)

    def test_missing_seed_and_wrong_denominator_rejected(self):
        # One pair is unavailable: its missing member must not turn into a zero difference.
        field = ('reference',None,'R','mean_absolute_shift')
        expected = {('S','a1'):{field:(10,3)}, ('S','a2'):{field:(None,0)},
                    ('S','b1'):{field:(8,3)}, ('S','b2'):{field:(5,3)}}
        arms = {('a',0.5,1):'a1',('a',0.5,2):'a2',('b',0.5,1):'b1',('b',0.5,2):'b2'}
        def stats(mean, sd, lo, hi, available):
            return dict(requested=2, available=available, unavailable=2-available, mean=mean, median=mean, min=lo, max=hi, sample_sd=sd)
        common=dict(rate=0.5, scope='reference', cutoff=None, reference_set='R', metric='mean_absolute_shift')
        score=dict(pgs_id='S', percentile_summaries=[
            dict(common, policy='a', distribution=stats(10,None,10,10,1), seeds=[dict(seed=1,value=10,evaluable=3),dict(seed=2,value=None,evaluable=0)]),
            dict(common, policy='b', distribution=stats(6.5,3/2**0.5,5,8,2), seeds=[dict(seed=1,value=8,evaluable=3),dict(seed=2,value=5,evaluable=3)])],
            percentile_comparisons=[dict(common, policy='b', reference_policy='a', distribution=stats(-2,None,-2,-2,1), pairs=[
                dict(seed=1,policy_value=8,reference_value=10,difference=-2,policy_evaluable=3,reference_evaluable=3),
                dict(seed=2,policy_value=5,reference_value=None,difference=None,policy_evaluable=3,reference_evaluable=0)])])
        self.assertEqual(check_percentile_summary(score, expected, arms, 'a'),3)
        for mode in ['denominator','missing_seed','null_difference','mean']:
            with self.subTest(mode=mode):
                bad=copy.deepcopy(score)
                row=bad['percentile_comparisons'][0]
                if mode=='denominator': row['pairs'][0]['reference_evaluable']=2
                if mode=='missing_seed': row['pairs'].pop()
                if mode=='null_difference': row['pairs'][1]['difference']=0
                if mode=='mean': row['distribution']['mean']=1.5
                with self.assertRaises(ValueError):
                    check_percentile_summary(bad, expected, arms, 'a')



if __name__ == '__main__':
    unittest.main()
