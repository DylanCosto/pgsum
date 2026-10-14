#!/usr/bin/env python3
"""Run and independently check a wholly local, synthetic robustness example."""

import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b''):
            result.update(chunk)
    return 'sha256:' + result.hexdigest()


def read_rows(path):
    with path.open(newline='') as source:
        return list(csv.DictReader(source, delimiter='\t'))


def percentile(value, reference):
    less = sum(item < value for item in reference)
    equal = sum(item == value for item in reference)
    return 100 * (less + equal / 2) / len(reference)


def run(binary, out, command, commands):
    full = [str(binary)] + list(map(str, command))
    commands.append(full)
    (out / 'commands.json').write_text(json.dumps(commands, indent=2) + '\n')
    with (out / ('command-%02d.log' % len(commands))).open('w') as log:
        subprocess.run(full, stdout=log, stderr=subprocess.STDOUT, check=True,
                       env={**os.environ, 'PGSUM_CACHE_DIR': str(out / 'cache')})


def verify(out, expected):
    results = out / 'results'
    report = json.loads((results / 'report.json').read_text())
    require(report['scores'][0]['baseline_terms'] == 3, 'all three terms must remain eligible')
    base = {sample: Decimal(score) for sample, score in expected['baseline'].items()}
    references = [Decimal(score) for score in expected['reference'].values()]
    values = {scenario: {sample: Decimal(score) for sample, score in samples.items()}
              for scenario, samples in expected['scenarios'].items()}
    sample_rows = read_rows(results / 'samples.tsv')
    require(len(sample_rows) == 16, 'expected four samples in four scenarios')
    seen = set()
    for row in sample_rows:
        scenario, sample = row['scenario'], row['sample']
        key = (scenario, sample)
        require(key not in seen, 'duplicate sample/scenario row')
        seen.add(key)
        score = values[scenario][sample]
        require(Decimal(row['baseline_score']) == base[sample], 'baseline arithmetic differs')
        require(Decimal(row['scenario_score']) == score, 'scenario arithmetic differs')
        require(Decimal(row['delta']) == score - base[sample], 'delta arithmetic differs')
        masked = expected['masked_terms'][scenario]
        filled = expected['filled_terms'][scenario]
        require(int(row['baseline_terms']) == 3, 'baseline term count differs')
        require(int(row['masked_terms']) == masked, 'masked term count differs')
        require(int(row['filled_terms']) == filled, 'filled term count differs')
        require(int(row['retained_terms']) == 3 - masked, 'retained term count differs')
        require(int(row['unfilled_terms']) == masked - filled, 'unfilled term count differs')
        require(int(row['unavailable_alternate_terms']) == 0, 'unexpected unavailable alternate call')
        require(float(row['baseline_percentile']) == percentile(base[sample], list(base.values())),
                'baseline percentile differs')
        require(float(row['fixed_baseline_percentile']) == percentile(score, list(base.values())),
                'fixed-baseline percentile differs')
        if scenario == 'all_missing':
            require(row['reranked_percentile'] == '', 'constant scenarios have no informative percentile')
        else:
            require(float(row['reranked_percentile']) == percentile(score, list(values[scenario].values())),
                    'reranked percentile differs')
        require(Decimal(row['top_terms_delta']) + Decimal(row['other_terms_delta']) == score - base[sample],
                'attribution does not reconcile')
    reference_rows = read_rows(results / 'references.tsv')
    require(len(reference_rows) == 16, 'expected 16 reference comparisons')
    for row in reference_rows:
        scenario, sample = row['scenario'], row['sample']
        before = percentile(base[sample], references)
        after = percentile(values[scenario][sample], references)
        require(int(row['shared_sample_ids']) == 0, 'target/reference samples overlap')
        require(int(row['reference_samples']) == 4, 'reference sample count differs')
        require(float(row['baseline_percentile']) == before, 'reference baseline percentile differs')
        require(float(row['scenario_percentile']) == after, 'reference scenario percentile differs')
        require(float(row['percentile_delta']) == after - before, 'reference percentile delta differs')
    region_rows = read_rows(results / 'regions.tsv')
    require(len(region_rows) == 16, 'expected 16 region rows')
    for row in region_rows:
        scenario, sample = row['scenario'], row['sample']
        original = Decimal(expected['last_site_contribution'][sample])
        delta = Decimal(0) if scenario == 'control' else (Decimal(3) if scenario == 'drop_last_filled' else Decimal(0)) - original
        require(Decimal(row['delta']) == delta, 'regional delta differs')
    for name, expected_hash in report['outputs'].items():
        require(digest(results / name) == expected_hash, 'output checksum differs: ' + name)
    return {'sample_scenario_checks': len(sample_rows), 'reference_checks': len(reference_rows),
            'region_checks': len(region_rows), 'output_hash_checks': len(report['outputs'])}


def study(binary, out, source, commands):
    """Use the same synthetic inputs to demonstrate paired policies and missing metrics."""
    template = json.loads((out / 'experiment.json').read_text())
    template['scenarios'] = [
        {'id': 'omit'},
        {'id': 'fill', 'fill_frequencies': 'reference.pgsf',
         'frequency_source_note': 'Frozen frequencies from synthetic R1-R4; no target participants contribute.'}
    ]
    policies = out / 'policies.json'
    policies.write_text(json.dumps(template, indent=2) + '\n')
    plan, results, summary = (out / name for name in ('study-plan', 'study-results', 'study-summary'))
    run(binary, out, ['robustness-plan', '--template', policies, '--rates', '0,0.1,0.2,1',
                     '--seeds', '11,22,33', '--out', plan], commands)
    run(binary, out, ['robustness', '--experiment', plan / 'experiment.json', '--out', results], commands)
    run(binary, out, ['robustness-summarize', '--plan', plan, '--results', results,
                     '--reference-policy', 'omit', '--out', summary], commands)
    checker = source.parents[1] / 'bench/check_study_summary.py'
    command = [sys.executable, str(checker), str(plan), str(results), str(summary)]
    commands.append(command)
    (out / 'commands.json').write_text(json.dumps(commands, indent=2) + '\n')
    checker_hash = digest(checker)
    verification = subprocess.check_output(command, text=True)
    require(digest(checker) == checker_hash, 'study checker changed during verification')
    checks = json.loads(verification)
    require(checks['passed'] and checks['sample_scenario_rows'] == 104, 'incomplete study example')
    result = json.loads((summary / 'summary.json').read_text())
    for name, expected_hash in result['outputs'].items():
        require(digest(summary / name) == expected_hash, 'study output checksum differs: ' + name)
    # Boundary rate 1 erases target variation; correlation must be unavailable, not zero.
    rows = result['scores'][0]['summaries']
    missing = [r for r in rows if r['policy'] == 'omit' and r['rate'] == 1 and r['metric'] in ('pearson', 'spearman')]
    require(len(missing) == 2 and all(r['distribution']['available'] == 0 for r in missing),
            'constant boundary scenario must have unavailable correlations')
    mapping = out / 'groups.tsv'
    mapping.write_text('sample\tgroup\nA\tsynthetic_singleton\nB\tsynthetic_group\nC\tsynthetic_group\nD\tsynthetic_group\n')
    grouped = out / 'study-groups'
    run(binary, out, ['robustness-groups', '--results', results, '--groups', mapping, '--out', grouped], commands)
    group_checker = source.parents[1] / 'bench/check_group_summary.py'
    command = [sys.executable, str(group_checker), str(results), str(mapping), str(grouped)]
    commands.append(command)
    (out / 'commands.json').write_text(json.dumps(commands, indent=2) + '\n')
    group_checker_hash = digest(group_checker)
    group_verification = subprocess.check_output(command, text=True)
    require(digest(group_checker) == group_checker_hash, 'group checker changed during verification')
    group_checks = json.loads(group_verification)
    require(group_checks['passed'] and group_checks['group_scenario_rows'] == 52, 'incomplete group example')
    checks['groups'] = group_checks
    group_study = out / 'group-study'
    run(binary, out, ['robustness-summarize', '--plan', plan, '--results', results,
                     '--groups', mapping, '--reference-policy', 'omit', '--out', group_study], commands)
    group_study_checker = source.parents[1] / 'bench/check_group_study.py'
    command = [sys.executable, str(group_study_checker), str(plan), str(results), str(mapping), str(group_study)]
    commands.append(command)
    (out / 'commands.json').write_text(json.dumps(commands, indent=2) + '\n')
    checker_hash = digest(group_study_checker)
    group_study_verification = subprocess.check_output(command, text=True)
    require(digest(group_study_checker) == checker_hash, 'group study checker changed during verification')
    checks['group_study'] = json.loads(group_study_verification)
    require(checks['group_study']['passed'] and checks['group_study']['score_group_rows'] == 2, 'group study inventory')
    # An aggregate-only checker could miss broken pairing. Corrupt one seed value and
    # require rejection while retaining the original verified result unmodified.
    tampered = out / 'group-study-tampered'
    shutil.copytree(group_study, tampered)
    altered = json.loads((tampered / 'summary.json').read_text())
    pair = altered['scores'][0]['comparisons'][0]['pairs'][0]
    pair['difference'] = 123456.0
    (tampered / 'summary.json').write_text(json.dumps(altered))
    rejection = subprocess.run(command[:-1] + [str(tampered)], capture_output=True, text=True)
    require(rejection.returncode != 0 and '123456' in rejection.stderr, 'checker accepted a changed seed pair')
    checks['group_study']['changed_pair_rejected'] = True
    (out / 'group-study-verification.json').write_text(json.dumps(checks['group_study'], indent=2) + '\n')
    shutil.rmtree(tampered)
    (out / 'group-verification.json').write_text(group_verification)
    (out / 'study-verification.json').write_text(json.dumps(checks, indent=2) + '\n')
    print('Offline group analysis passed: ' + str(grouped / 'groups.html'))
    print('Offline study passed: ' + str(out / 'study-verification.json'))
    print('Open study summary: ' + str(summary / 'summary.html'))
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True, help='new output directory')
    parser.add_argument('--study', action='store_true', help='also plan, run and independently check paired seed/rate comparisons (Python 3.9+)')
    args = parser.parse_args()
    require(not args.study or sys.version_info >= (3, 9), '--study requires Python 3.9 or newer')
    binary = args.binary.expanduser().resolve()
    require(binary.is_file(), 'binary does not exist: ' + str(binary))
    source = Path(__file__).resolve().parent
    out = args.out.expanduser().resolve()
    out.mkdir()  # Refuse to overwrite an earlier run.
    inputs = out / 'inputs'
    inputs.mkdir()
    for name in ['reference.fa', 'reference.fa.fai', 'score.tsv', 'targets.vcf', 'reference.vcf',
                 'frequencies.vcf', 'drop-last.bed']:
        shutil.copyfile(source / name, inputs / name)
    with (inputs / 'frequencies.vcf.gz').open('wb') as target:
        with gzip.GzipFile(filename='', mode='wb', fileobj=target, mtime=0) as compressed:
            compressed.write((inputs / 'frequencies.vcf').read_bytes())
    shutil.copyfile(source / 'experiment.json', out / 'experiment.json')
    identity = digest(binary)
    commands = []
    run(binary, out, ['compile', inputs / 'score.tsv', '--reference', inputs / 'reference.fa',
                      '--out', out / 'packs'], commands)
    for name in ['targets', 'reference']:
        run(binary, out, ['extract', '--gvcf', inputs / (name + '.vcf'), '--reference', inputs / 'reference.fa',
                          '--pack', out / 'packs', '--all-samples', '--accept-missing-quality',
                          '--out', out / (name + '.pgsc')], commands)
    run(binary, out, ['frequencies', '--vcf', inputs / 'frequencies.vcf.gz', '--field', 'AF',
                      '--out', out / 'reference.pgsf'], commands)
    run(binary, out, ['robustness', '--experiment', out / 'experiment.json', '--out', out / 'preflight', '--preflight', '--terms'], commands)
    run(binary, out, ['robustness', '--experiment', out / 'experiment.json', '--out', out / 'results', '--terms'], commands)
    checks = verify(out, json.loads((source / 'expected.json').read_text()))
    preflight = json.loads((out / 'preflight/preflight.json').read_text())
    native = json.loads((out / 'results/report.json').read_text())
    require(preflight['schema'] == 'pgsum-robustness-preflight-v1' and preflight['scoring_performed'] is False, 'preflight scope')
    require(preflight['input_sha256'] == native['input_sha256'], 'preflight/native input inventory differs')
    require(preflight['dimensions']['candidate_term_sample_scenario_cells'] == 48, 'preflight model dimensions')
    require(preflight['dimensions']['samples_tsv_rows_if_complete'] == 16, 'preflight score rows')
    require(preflight['dimensions']['reference_tsv_rows_if_complete'] == 16, 'preflight reference rows')
    require(preflight['reference_sets'][0]['shared_target_sample_count'] == 0, 'preflight target/reference overlap')
    require(not (out / 'preflight/report.json').exists() and not (out / 'preflight/samples.tsv').exists(), 'preflight emitted scoring evidence')
    require(preflight['outputs']['preflight.txt'] == digest(out / 'preflight/preflight.txt'), 'preflight text checksum')
    checks['preflight'] = dict(passed=True, input_inventory_matches_native=True, scoring_performed=False)
    if args.study:
        checks['study'] = study(binary, out, source, commands)
    require(digest(binary) == identity, 'binary changed during the example')
    record = {'schema': 'pgsum-robustness-example-verification-v1', 'status': 'passed',
              'binary_sha256': identity, 'version': subprocess.check_output([str(binary), '--version'], text=True).strip(),
              'checks': checks, 'scope': 'Synthetic offline arithmetic and package integration only.'}
    (out / 'verification.json').write_text(json.dumps(record, indent=2) + '\n')
    print('Offline example passed: ' + str(out / 'verification.json'))
    print('Open report: ' + str(out / 'results/report.html'))


if __name__ == '__main__':
    main()
