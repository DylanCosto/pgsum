#!/usr/bin/env python3
"""Independently check study metrics/seed aggregation from exact native sample sums.

This validates descriptive summaries, not genotype-to-score arithmetic or prediction.
Python 3.9+, standard library only. Inputs: plan directory, results directory, summary directory.
"""
import csv
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys


def require(ok, message):
    if not ok:
        raise ValueError(message)


def close(actual, expected):
    require((actual is None and expected is None) or
            (actual is not None and expected is not None and math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-10)),
            f'{actual!r} != {expected!r}')


def percentiles(values, reference):
    if len(reference) < 2 or len(set(reference)) < 2:
        return None
    return [Fraction(100 * (2 * sum(r < v for r in reference) + sum(r == v for r in reference)), 2 * len(reference)) for v in values]


def correlation(a, b):
    if a is None or b is None or len(a) < 2:
        return None
    ac = [float(x - sum(a) / len(a)) for x in a]
    bc = [float(x - sum(b) / len(b)) for x in b]
    denominator = math.sqrt(math.fsum(x*x for x in ac) * math.fsum(x*x for x in bc))
    return math.fsum(x*y for x,y in zip(ac,bc)) / denominator if denominator else None


def metrics(rows):
    a = [Fraction(r['baseline_score']) for r in rows]
    b = [Fraction(r['scenario_score']) for r in rows]
    d = [abs(x-y) for x,y in zip(a,b)]
    pa, pb = percentiles(a,a), percentiles(b,a)
    shift = [abs(x-y) for x,y in zip(pa,pb)] if pa is not None else None
    return dict(mae=float(sum(d)/len(d)), rmse=math.sqrt(float(sum(x*x for x in d)/len(d))),
                max_abs_change=float(max(d)), pearson=correlation(a,b),
                spearman=correlation(percentiles(a,a),percentiles(b,b)),
                mean_abs_percentile_shift=float(sum(shift)/len(shift)) if shift is not None else None,
                max_abs_percentile_shift=float(max(shift)) if shift is not None else None)


def distribution(actual, values):
    xs = [v for v in values if v is not None]
    require((actual['requested'],actual['available'],actual['unavailable']) == (len(values),len(xs),len(values)-len(xs)), 'availability counts')
    for key, expected in dict(mean=statistics.mean(xs) if xs else None,
                              median=statistics.median(xs) if xs else None,
                              min=min(xs) if xs else None, max=max(xs) if xs else None,
                              sample_sd=statistics.stdev(xs) if len(xs)>1 else None).items():
        close(actual[key],expected)


def percentile_evidence(grouped, native, results_dir):
    """Cutoffs from exact sums; frozen-reference shifts from participant reference rows."""
    declared = native['experiment']
    refs = {r['id'] for r in declared['reference_sets']}
    with (results_dir/'references.tsv').open() as stream:
        reference_rows = list(csv.DictReader(stream, delimiter='\t'))
    reference_groups = {}
    identities = set()
    for r in reference_rows:
        identity = (r['pgs_id'], r['scenario'], r['reference'], r['sample'])
        require(identity not in identities, 'duplicate reference participant')
        identities.add(identity)
        reference_groups.setdefault(identity[:3], []).append(r)
    wanted = {(score, scenario, ref) for score, scenario in grouped for ref in refs}
    require(set(reference_groups) == wanted, 'reference row inventory')
    expected = {}
    frozen = {}
    for key, rows in grouped.items():
        out = {}
        a = [Fraction(r['baseline_score']) for r in rows]
        b = [Fraction(r['scenario_score']) for r in rows]
        pa, pb = percentiles(a, a), percentiles(b, a)
        for cutoff in declared['cutoffs']:
            pairs = list(zip(pa, pb)) if pa is not None else []
            up = sum(float(x) < cutoff <= float(y) for x, y in pairs)
            down = sum(float(y) < cutoff <= float(x) for x, y in pairs)
            n = len(pairs)
            for metric, numerator in dict(up_count=up, down_count=down, unchanged_count=n-up-down,
                                         changed_count=up+down, up_fraction=up, down_fraction=down,
                                         changed_fraction=up+down).items():
                v = numerator/n if n and metric.endswith('_fraction') else numerator
                out['cutoff', cutoff, None, metric] = (v if n else None, n)
        for ref in refs:
            rr = reference_groups[key + (ref,)]
            require({r['sample'] for r in rr} == {r['sample'] for r in rows}, 'reference sample identity')
            counts = {(int(r['reference_samples']), int(r['shared_sample_ids'])) for r in rr}
            require(len(counts) == 1, 'reference counts change within scenario')
            samples, shared = next(iter(counts))
            require(samples > 0 and 0 <= shared <= min(samples, len(rows)), 'reference counts invalid')
            pairs = []
            for r in rr:
                vals = [float(r[k]) if r[k] else None for k in ('baseline_percentile', 'scenario_percentile', 'percentile_delta')]
                x, y, delta = vals
                require(all(v is None for v in vals) or all(v is not None and math.isfinite(v) for v in vals), 'partial reference row')
                if x is not None:
                    require(0 <= x <= 100 and 0 <= y <= 100 and samples >= 2, 'reference percentile bounds')
                    close(delta, y-x)
                    pairs.append((x, y))
                identity = (key[0], ref, r['sample'])
                fingerprint = (samples, shared, x)
                require(identity not in frozen or frozen[identity] == fingerprint, 'frozen reference baseline changed')
                frozen[identity] = fingerprint
            n = len(pairs)
            shifts = [abs(y-x) for x,y in pairs]
            values = dict(mean_baseline_percentile=statistics.mean(x for x,y in pairs) if n else None,
                          mean_scenario_percentile=statistics.mean(y for x,y in pairs) if n else None,
                          mean_absolute_shift=statistics.mean(shifts) if n else None,
                          max_absolute_shift=max(shifts) if n else None)
            for metric, v in values.items():
                out['reference', None, ref, metric] = (v, n)
        expected[key] = out
    for score in native['scores']:
        for scenario in score['scenarios']:
            evidence = expected[score['pgs_id'], scenario['id']]
            require(len(scenario['metrics']['crossings']) == len(declared['cutoffs']), 'native cutoff count')
            require({c['cutoff'] for c in scenario['metrics']['crossings']} == set(declared['cutoffs']), 'native cutoff inventory')
            for c in scenario['metrics']['crossings']:
                for field in ('up', 'down', 'unchanged'):
                    v, n = evidence['cutoff', c['cutoff'], None, field+'_count']
                    close(c[field], v if n else 0)
                    require(c['evaluable'] == n, 'native cutoff denominator')
            require(len(scenario['reference_comparisons']) == len(refs), 'native reference count')
            require({r['id'] for r in scenario['reference_comparisons']} == refs, 'native reference inventory')
            for r in scenario['reference_comparisons']:
                samples, shared, _ = frozen[score['pgs_id'], r['id'], next(iter(grouped[score['pgs_id'],scenario['id']]))['sample']]
                require((r['reference_samples'], r['shared_sample_ids']) == (samples, shared), 'native reference counts')
                for dimension, (v, n) in evidence.items():
                    if dimension[0] == 'reference' and dimension[2] == r['id']:
                        close(r[dimension[3]], v)
                        require(r['evaluable'] == n, 'native reference denominator')
    return expected, len(reference_rows)


def check_percentile_summary(score, expected, arms, reference_policy):
    def dim(row):
        return row['scope'], row['cutoff'], row['reference_set'], row['metric']
    fields = set(next(v for (pgs, scenario), v in expected.items() if pgs == score['pgs_id']))
    checked = 0
    for kind in ('percentile_summaries', 'percentile_comparisons'):
        paired = kind.endswith('comparisons')
        wanted = {(policy, rate, field) for policy, rate, seed in arms for field in fields if not paired or policy != reference_policy}
        actual = [(r['policy'], r['rate'], dim(r)) for r in score[kind]]
        require(len(actual) == len(wanted) and set(actual) == wanted, 'percentile summary grid')
        for row in score[kind]:
            seeds = {seed for policy, rate, seed in arms if policy == row['policy'] and rate == row['rate']}
            observations = row['pairs' if paired else 'seeds']
            require(len(observations) == len(seeds) and {p['seed'] for p in observations} == seeds, 'percentile seed inventory')
            values = []
            for p in observations:
                key = row['policy'], row['rate'], p['seed']
                v, n = expected[score['pgs_id'], arms[key]][dim(row)]
                if paired:
                    require(row['reference_policy'] == reference_policy, 'wrong paired reference policy')
                    ref = reference_policy, row['rate'], p['seed']
                    r, rn = expected[score['pgs_id'], arms[ref]][dim(row)]
                    d = v-r if v is not None and r is not None else None
                    close(p['policy_value'], v); close(p['reference_value'], r); close(p['difference'], d)
                    require((p['policy_evaluable'], p['reference_evaluable']) == (n, rn), 'paired percentile denominators')
                    values.append(d)
                else:
                    close(p['value'], v)
                    require(p['evaluable'] == n, 'percentile denominator')
                    values.append(v)
            distribution(row['distribution'], values); checked += 1
    return checked


def main():
    plan_dir, results_dir, summary_dir = map(Path, sys.argv[1:])
    plan=json.loads((plan_dir/'plan.json').read_text())
    summary=json.loads((summary_dir/'summary.json').read_text())
    with (results_dir/'samples.tsv').open() as stream:
        rows=list(csv.DictReader(stream, delimiter='\t'))
    grouped={}
    for row in rows:
        grouped.setdefault((row['pgs_id'],row['scenario']),[]).append(row)
    identities=[(r['pgs_id'],r['scenario'],r['sample']) for r in rows]
    require(len(set(identities))==len(identities),'duplicate sample/scenario row')
    per_score={}
    for (score,scenario),values in grouped.items():
        baseline={r['sample']:Fraction(r['baseline_score']) for r in values}
        require(score not in per_score or per_score[score]==baseline,'baseline/sample identity changed')
        per_score[score]=baseline
    expected={key:metrics(values) for key,values in grouped.items()}
    arms={(a['policy'],a['rate'],a['seed']):a['scenario'] for a in plan['scenarios'] if a['seed'] is not None}
    require(len(arms)==sum(a['seed'] is not None for a in plan['scenarios']),'duplicate plan cell')
    require(set(grouped)=={(score,a['scenario']) for score in per_score for a in plan['scenarios']},'incomplete sample evidence')
    native=json.loads((results_dir/'report.json').read_text())
    native_keys=[(s['pgs_id'],r['id']) for s in native['scores'] for r in s['scenarios']]
    require(len(native_keys)==len(expected) and set(native_keys)==set(expected),'native scenario inventory')
    for score in native['scores']:
        for scenario in score['scenarios']:
            for metric,v in expected[score['pgs_id'],scenario['id']].items():
                close(scenario['metrics'][metric],v)
    percentile_expected, reference_rows = percentile_evidence(grouped, native, results_dir)
    checked=0
    percentile_checked=0
    require({s['pgs_id'] for s in summary['scores']}==set(per_score),'summary score inventory')
    require(len(summary['scores'])==len(per_score),'duplicate summary score')
    for score in summary['scores']:
        native_score = next(s for s in native['scores'] if s['pgs_id'] == score['pgs_id'])
        refs = native_score['scenarios'][0]['reference_comparisons']
        require(len(score['reference_sets']) == len(refs) and {r['id'] for r in score['reference_sets']} == {r['id'] for r in refs}, 'summary reference identity')
        for r in score['reference_sets']:
            expected_ref = next(v for v in refs if v['id'] == r['id'])
            for key in ('reference_samples','shared_sample_ids','evaluable','mean_baseline_percentile'):
                close(r[key], expected_ref[key])
        percentile_checked += check_percentile_summary(score, percentile_expected, arms, summary['reference_policy'])
        fields=set(next(iter(expected.values())))
        wanted={(policy,rate,metric) for policy,rate,seed in arms for metric in fields}
        actual=[(r['policy'],r['rate'],r['metric']) for r in score['summaries']]
        require(len(actual)==len(wanted) and set(actual)==wanted,'summary grid incomplete')
        wanted_pairs={k for k in wanted if k[0]!=summary['reference_policy']}
        actual=[(r['policy'],r['rate'],r['metric']) for r in score['comparisons']]
        require(len(actual)==len(wanted_pairs) and set(actual)==wanted_pairs,'comparison grid incomplete')
        for row in score['summaries']:
            values=[expected[score['pgs_id'],scenario][row['metric']] for (policy,rate,seed),scenario in arms.items() if policy==row['policy'] and rate==row['rate']]
            distribution(row['distribution'],values); checked+=1
        for row in score['comparisons']:
            require(row['reference_policy']==summary['reference_policy'],'wrong reference policy')
            deltas=[]
            seeds={seed for policy,rate,seed in arms if policy==row['policy'] and rate==row['rate']}
            require(len(row['pairs'])==len(seeds) and {p['seed'] for p in row['pairs']}==seeds,'paired seeds incomplete')
            for pair in row['pairs']:
                key=(row['policy'],row['rate'],pair['seed']);ref=(row['reference_policy'],row['rate'],pair['seed'])
                v=expected[score['pgs_id'],arms[key]][row['metric']];r=expected[score['pgs_id'],arms[ref]][row['metric']]
                d=v-r if v is not None and r is not None else None
                close(pair['policy_value'],v);close(pair['reference_value'],r);close(pair['difference'],d);deltas.append(d)
            distribution(row['distribution'],deltas);checked+=1
    print(json.dumps({'passed':True,'sample_scenario_rows':len(rows),'scenario_metric_checks':len(expected)*7,'seed_distribution_checks':checked,'percentile_distribution_checks':percentile_checked,'reference_participant_rows':reference_rows,
                      'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      'scope':'Score/cutoff metrics recalculated from native exact sums; frozen-reference shifts from participant reference rows. No new source-genotype or phenotype validation.'},indent=2))


if __name__=='__main__':
    main()
