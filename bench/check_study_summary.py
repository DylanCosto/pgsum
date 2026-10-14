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


def main():
    plan_dir, results_dir, summary_dir = map(Path, sys.argv[1:])
    plan=json.loads((plan_dir/'plan.json').read_text())
    summary=json.loads((summary_dir/'summary.json').read_text())
    rows=list(csv.DictReader((results_dir/'samples.tsv').open(), delimiter='\t'))
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
    checked=0
    require({s['pgs_id'] for s in summary['scores']}==set(per_score),'summary score inventory')
    require(len(summary['scores'])==len(per_score),'duplicate summary score')
    for score in summary['scores']:
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
            deltas=[]
            seeds={seed for policy,rate,seed in arms if policy==row['policy'] and rate==row['rate']}
            require(len(row['pairs'])==len(seeds) and {p['seed'] for p in row['pairs']}==seeds,'paired seeds incomplete')
            for pair in row['pairs']:
                key=(row['policy'],row['rate'],pair['seed']);ref=(row['reference_policy'],row['rate'],pair['seed'])
                v=expected[score['pgs_id'],arms[key]][row['metric']];r=expected[score['pgs_id'],arms[ref]][row['metric']]
                d=v-r if v is not None and r is not None else None
                close(pair['policy_value'],v);close(pair['reference_value'],r);close(pair['difference'],d);deltas.append(d)
            distribution(row['distribution'],deltas);checked+=1
    print(json.dumps({'passed':True,'sample_scenario_rows':len(rows),'scenario_metric_checks':len(expected)*7,'seed_distribution_checks':checked,
                      'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      'scope':'Summary metrics recalculated from native exact sample sums; no new source-genotype or phenotype validation.'},indent=2))


if __name__=='__main__':
    main()
