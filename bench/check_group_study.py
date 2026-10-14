#!/usr/bin/env python3
"""Check paired group studies from native exact sums. Python 3.9+, stdlib.

Arguments: PLAN NATIVE_RESULTS GROUP_MAPPING GROUP_STUDY.
Descriptive sensitivity only; no genotype reconstruction or phenotype validation.
"""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import statistics
import sys

from check_group_summary import read_rows
from check_study_summary import close, distribution, metrics, percentiles, require, check_percentile_summary


def check(plan_dir, native_dir, mapping_path, output):
    plan = json.loads((plan_dir/'plan.json').read_text())
    native = json.loads((native_dir/'report.json').read_text())
    summary = json.loads((output/'summary.json').read_text())
    require(summary['schema'] == 'pgsum-robustness-group-study-v1', 'group study schema')
    require(summary['percentile_basis'] == 'full_target_baseline' and summary['correlation_basis'] == 'within_group', 'statistical basis')
    fingerprints = set(summary['source_sha256'].values())
    for path in [plan_dir/'plan.json', plan_dir/'experiment.json', native_dir/'report.json', mapping_path]:
        require('sha256:'+hashlib.sha256(path.read_bytes()).hexdigest() in fingerprints, 'source fingerprint')
    require(plan['experiment_sha256'] == 'sha256:'+hashlib.sha256((plan_dir/'experiment.json').read_bytes()).hexdigest(), 'plan fingerprint')
    for directory, record in [(native_dir, native), (output, summary)]:
        for name, digest in record['outputs'].items():
            require(Path(name).name == name, 'output path')
            require('sha256:'+hashlib.sha256((directory/name).read_bytes()).hexdigest() == digest, 'output fingerprint')
            if directory == native_dir:
                require(digest in fingerprints, 'native evidence fingerprint')
    mapping = read_rows(mapping_path)
    assignments = {r['sample']: r['group'] for r in mapping}
    require(len(assignments) == len(mapping) and assignments, 'mapping inventory')
    groups = {}
    for sample, group in assignments.items():
        groups.setdefault(group, set()).add(sample)
    require(len(summary['groups']) == len(groups), 'group inventory')
    require({g['id']: set(g['members']) for g in summary['groups']} == groups, 'group membership')
    for g in summary['groups']:
        require(g['samples'] == len(g['members']) == len(groups[g['id']]), 'group denominator')
    rows = read_rows(native_dir/'samples.tsv')
    grouped = {}
    identities = set()
    for row in rows:
        key = row['pgs_id'], row['scenario']
        identity = key+(row['sample'],)
        require(identity not in identities, 'duplicate native sample')
        identities.add(identity)
        grouped.setdefault(key, []).append(row)
    refs = {}
    for row in read_rows(native_dir/'references.tsv'):
        key = row['pgs_id'], row['scenario'], row['reference'], row['sample']
        require(key not in refs, 'duplicate reference row')
        refs[key] = row
    arms = {(a['policy'],a['rate'],a['seed']): a['scenario'] for a in plan['scenarios'] if a['seed'] is not None}
    require(len(arms) == sum(a['seed'] is not None for a in plan['scenarios']), 'duplicate plan cell')
    score_ids = {s['pgs_id'] for s in native['scores']}
    wanted = {(score,a['scenario']) for score in score_ids for a in plan['scenarios']}
    require(set(grouped) == wanted, 'native scenario inventory')
    actual = [(s['pgs_id'],s['group']) for s in summary['scores']]
    require(len(actual) == len(set(actual)) and set(actual) == {(s,g) for s in score_ids for g in groups}, 'summary inventory')
    checked = percentile_checked = 0
    for score in summary['scores']:
        pgs, group = score['pgs_id'], score['group']
        members = groups[group]
        original = next(s for s in native['scores'] if s['pgs_id'] == pgs)
        for field in ('baseline_terms','baseline_term_sha256','total_model_terms','baseline_term_fraction','baseline_weight_fraction','excluded_terms'):
            require(score[field] == original[field], 'baseline identity')
        require(score['samples'] == len(members) and score['source_samples'] == len(assignments), 'participant counts')
        expected, percentile_expected, reference_identity = {}, {}, {}
        for (sid, scenario), all_rows in grouped.items():
            if sid != pgs:
                continue
            require({r['sample'] for r in all_rows} == set(assignments), 'native participant inventory')
            selected = [r for r in all_rows if r['sample'] in members]
            baseline = [Fraction(r['baseline_score']) for r in all_rows]
            a = [Fraction(r['baseline_score']) for r in selected]
            b = [Fraction(r['scenario_score']) for r in selected]
            pa, pb = percentiles(a, baseline), percentiles(b, baseline)
            shifts = [abs(x-y) for x,y in zip(pa,pb)] if pa is not None else []
            values = metrics(selected)
            values.update(mean_abs_percentile_shift=float(sum(shifts)/len(shifts)) if shifts else None,
                          max_abs_percentile_shift=float(max(shifts)) if shifts else None)
            expected[scenario] = values
            evidence = {}
            for cutoff in native['experiment']['cutoffs']:
                pairs = list(zip(pa,pb)) if pa is not None else []
                up = sum(float(x) < cutoff <= float(y) for x,y in pairs)
                down = sum(float(y) < cutoff <= float(x) for x,y in pairs)
                n = len(pairs)
                for metric, numerator in dict(up_count=up, down_count=down, unchanged_count=n-up-down,
                                             changed_count=up+down,up_fraction=up,down_fraction=down,changed_fraction=up+down).items():
                    v = numerator/n if n and metric.endswith('_fraction') else numerator
                    evidence['cutoff',cutoff,None,metric] = (v if n else None,n)
            for ref in native['experiment']['reference_sets']:
                rr = [refs[pgs,scenario,ref['id'],sample] for sample in sorted(members)]
                pairs = [(float(r['baseline_percentile']),float(r['scenario_percentile'])) for r in rr if r['baseline_percentile'] and r['scenario_percentile']]
                shifts = [abs(x-y) for x,y in pairs]
                n = len(pairs)
                ref_values = dict(mean_baseline_percentile=statistics.mean(x for x,y in pairs) if n else None,
                                  mean_scenario_percentile=statistics.mean(y for x,y in pairs) if n else None,
                                  mean_absolute_shift=statistics.mean(shifts) if n else None,
                                  max_absolute_shift=max(shifts) if n else None)
                for metric,v in ref_values.items():
                    evidence['reference',None,ref['id'],metric] = (v,n)
                identity = dict(reference_samples=int(rr[0]['reference_samples']),source_cohort_shared_sample_ids=int(rr[0]['shared_sample_ids']),
                                evaluable=n,mean_baseline_percentile=ref_values['mean_baseline_percentile'])
                require(ref['id'] not in reference_identity or reference_identity[ref['id']] == identity, 'frozen reference identity')
                reference_identity[ref['id']] = identity
            percentile_expected[pgs,scenario] = evidence
        require(len(score['reference_sets']) == len(reference_identity) and {r['id'] for r in score['reference_sets']} == set(reference_identity), 'reference inventory')
        for ref in score['reference_sets']:
            require('shared_sample_ids' not in ref, 'ambiguous reference overlap')
            for field,v in reference_identity[ref['id']].items():
                close(ref[field],v)
        controls = {a['scenario'] for a in plan['scenarios'] if a['seed'] is None}
        require(len(score['controls']) == len(controls) and {c['id'] for c in score['controls']} == controls, 'control inventory')
        for control in score['controls']:
            for metric,v in expected[control['id']].items():
                close(control['metrics'][metric],v)
        for kind in ('summaries','comparisons'):
            paired = kind == 'comparisons'
            wanted = {(policy,rate,m) for policy,rate,seed in arms for m in next(iter(expected.values())) if not paired or policy != summary['reference_policy']}
            actual = [(r['policy'],r['rate'],r['metric']) for r in score[kind]]
            require(len(actual) == len(wanted) and set(actual) == wanted, 'metric grid')
            for row in score[kind]:
                seeds = {seed for policy,rate,seed in arms if (policy,rate) == (row['policy'],row['rate'])}
                values = []
                if paired:
                    require(row['reference_policy'] == summary['reference_policy'], 'reference policy')
                    require(len(row['pairs']) == len(seeds) and {p['seed'] for p in row['pairs']} == seeds, 'paired seed inventory')
                    for pair in row['pairs']:
                        v = expected[arms[row['policy'],row['rate'],pair['seed']]][row['metric']]
                        r = expected[arms[row['reference_policy'],row['rate'],pair['seed']]][row['metric']]
                        delta = v-r if v is not None and r is not None else None
                        close(pair['policy_value'],v); close(pair['reference_value'],r); close(pair['difference'],delta)
                        values.append(delta)
                else:
                    values = [expected[arms[row['policy'],row['rate'],seed]][row['metric']] for seed in seeds]
                distribution(row['distribution'],values)
                checked += 1
        percentile_checked += check_percentile_summary(score, percentile_expected, arms, summary['reference_policy'])
    return dict(passed=True,groups=len(groups),score_group_rows=len(summary['scores']),
                seed_distribution_checks=checked,percentile_distribution_checks=percentile_checked,
                verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                supporting_verifiers={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ('check_group_summary.py','check_study_summary.py')},
                scope='Within-group paired seed distributions from native exact sums and reference rows; full-target percentile basis. No genotype reconstruction or phenotype validation.')


if __name__ == '__main__':
    print(json.dumps(check(*map(Path,sys.argv[1:])),indent=2))
