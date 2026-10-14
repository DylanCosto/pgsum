#!/usr/bin/env python3
"""Independently verify group diagnostics from native exact sums and reference rows.

Python 3.9+, standard library. This does not reconstruct scores from source genotypes.
Usage: check_group_summary.py NATIVE_RESULTS GROUP_MAPPING_TSV GROUP_RESULTS
"""
import csv
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import statistics
import sys

from check_study_summary import close, metrics, percentiles, require


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t', quoting=csv.QUOTE_NONE))


def check(native_dir, mapping_path, group_dir):
    native = json.loads((native_dir/'report.json').read_text())
    summary = json.loads((group_dir/'groups.json').read_text())
    require(summary['schema'] == 'pgsum-robustness-groups-v1', 'group schema')
    require(summary['percentile_basis'] == 'full_target_baseline', 'percentile basis')
    require(summary['correlation_basis'] == 'within_group', 'correlation basis')
    for directory, report in [(native_dir, native), (group_dir, summary)]:
        for name, digest in report['outputs'].items():
            require(Path(name).name == name, 'unexpected output path')
            require('sha256:'+hashlib.sha256((directory/name).read_bytes()).hexdigest() == digest, 'output fingerprint')
    required_sources = [mapping_path, native_dir/'report.json', native_dir/'samples.tsv', native_dir/'references.tsv']
    recorded_sources = set(summary['source_sha256'].values())
    for path in required_sources:
        require('sha256:'+hashlib.sha256(path.read_bytes()).hexdigest() in recorded_sources, 'source fingerprint missing or changed')
    mapping = read_rows(mapping_path)
    assignments = {r['sample']: r['group'] for r in mapping}
    require(len(assignments) == len(mapping) and assignments, 'duplicate/empty participant mapping')
    groups = {}
    for sample, group in assignments.items():
        groups.setdefault(group, set()).add(sample)
    require({g['id']: set(g['members']) for g in summary['groups']} == groups, 'group membership')
    require(len(summary['groups']) == len(groups), 'duplicate groups')
    for g in summary['groups']:
        require(g['samples'] == len(g['members']) == len(groups[g['id']]), 'group size')
    grouped = {}
    identities = set()
    rows = read_rows(native_dir/'samples.tsv')
    for r in rows:
        key = (r['pgs_id'], r['scenario'])
        identity = key+(r['sample'],)
        require(identity not in identities, 'duplicate native participant row')
        identities.add(identity)
        grouped.setdefault(key, []).append(r)
    ref_rows = read_rows(native_dir/'references.tsv')
    refs = {(r['pgs_id'], r['scenario'], r['reference'], r['sample']): r for r in ref_rows}
    require(len(refs) == len(ref_rows), 'duplicate reference rows')
    expected = {}
    for key, pairs in grouped.items():
        require({r['sample'] for r in pairs} == set(assignments), 'participant inventory')
        reference = [Fraction(r['baseline_score']) for r in pairs]
        for group, members in groups.items():
            selected = [r for r in pairs if r['sample'] in members]
            values = metrics(selected)
            before = percentiles([Fraction(r['baseline_score']) for r in selected], reference)
            after = percentiles([Fraction(r['scenario_score']) for r in selected], reference)
            shifts = [abs(a-b) for a,b in zip(before, after)] if before is not None else []
            values['mean_abs_percentile_shift'] = float(sum(shifts)/len(shifts)) if shifts else None
            values['max_abs_percentile_shift'] = float(max(shifts)) if shifts else None
            expected[key+(group,)] = (selected, values, before, after)
    actual = [(s['pgs_id'],r['scenario'],r['group']) for s in summary['scores'] for r in s['scenarios']]
    require(len(actual) == len(set(actual)) and set(actual) == set(expected), 'group scenario inventory')
    metric_checks = cutoff_checks = reference_checks = 0
    for score in summary['scores']:
        original = next(s for s in native['scores'] if s['pgs_id'] == score['pgs_id'])
        for key in ['baseline_terms', 'baseline_term_sha256', 'total_model_terms', 'excluded_terms']:
            require(score[key] == original[key], 'baseline identity: '+key)
        require(score['source_samples'] == len(assignments), 'original cohort size')
        for row in score['scenarios']:
            key = (score['pgs_id'], row['scenario'], row['group'])
            selected, values, before, after = expected[key]
            require(row['metrics']['n'] == len(selected), 'metric denominator')
            for metric, value in values.items():
                close(row['metrics'][metric], value)
                if value is None:
                    require(row['metrics']['unavailable_reasons'].get(metric), 'missing unavailability reason')
                metric_checks += 1
            crossings = row['metrics']['crossings']
            require([c['cutoff'] for c in crossings] == native['experiment']['cutoffs'], 'cutoff inventory')
            for crossing in crossings:
                cutoff = crossing['cutoff']
                pairs = list(zip(before, after)) if before is not None else []
                up = sum(float(a) < cutoff <= float(b) for a,b in pairs)
                down = sum(float(b) < cutoff <= float(a) for a,b in pairs)
                require((crossing['up'],crossing['down'],crossing['unchanged'],crossing['evaluable']) == (up,down,len(pairs)-up-down,len(pairs)), 'global cutoff transitions')
                cutoff_checks += 1
            declared_refs = {r['id'] for r in native['experiment']['reference_sets']}
            require({r['id'] for r in row['reference_comparisons']} == declared_refs and len(row['reference_comparisons']) == len(declared_refs), 'reference inventory')
            for ref in row['reference_comparisons']:
                rr = [refs[(key[0],key[1],ref['id'],r['sample'])] for r in selected]
                metadata = {(int(r['reference_samples']),int(r['shared_sample_ids'])) for r in rr}
                require(metadata == {(ref['reference_samples'],ref['source_cohort_shared_sample_ids'])}, 'reference metadata')
                pairs = [(float(r['baseline_percentile']),float(r['scenario_percentile'])) for r in rr if r['baseline_percentile'] and r['scenario_percentile']]
                require((ref['evaluable'],ref['unavailable']) == (len(pairs),len(selected)-len(pairs)), 'reference denominators')
                shifts = [abs(a-b) for a,b in pairs]
                for field, vals in [('mean_baseline_percentile',[a for a,b in pairs]), ('mean_scenario_percentile',[b for a,b in pairs]), ('mean_absolute_shift',shifts)]:
                    close(ref[field], statistics.mean(vals) if vals else None)
                close(ref['max_absolute_shift'], max(shifts) if shifts else None)
                reference_checks += 1
    return dict(passed=True, sample_scenario_rows=len(rows), groups=len(groups), group_scenario_rows=len(actual),
                score_metric_checks=metric_checks, cutoff_checks=cutoff_checks, reference_checks=reference_checks,
                verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                supporting_verifier_sha256=hashlib.sha256(Path(__file__).with_name('check_study_summary.py').read_bytes()).hexdigest(),
                scope='Group metrics derived from native exact sums with whole-cohort percentile thresholds; reference shifts derived from native participant reference rows. No source-genotype or phenotype validation.')


if __name__ == '__main__':
    print(json.dumps(check(*map(Path, sys.argv[1:])), indent=2))
