#!/usr/bin/env python3
"""Compare native robustness memory while requiring identical numerical evidence.

Uses an existing experiment, two prebuilt native binaries and GNU time. Standard library,
Python 3.9+. Runs sequentially; elapsed times on a shared host are not speedup evidence.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

PAYLOADS = ('samples.tsv', 'top-changes.jsonl', 'regions.tsv', 'references.tsv')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as source:
        for chunk in iter(lambda: source.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def save(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def measure(binary, experiment, out, label, threads):
    binary_sha = sha(binary)
    report_dir = out / label
    timing = out / (label + '.time.json')
    command = [str(binary), '--threads', str(threads), 'robustness', '--experiment', str(experiment), '--out', str(report_dir)]
    with (out / (label + '.log')).open('w') as log:
        subprocess.run(['/usr/bin/time', '-f', '{"elapsed_seconds":%e,"peak_rss_kib":%M,"user_seconds":%U,"system_seconds":%S}',
                        '-o', str(timing), *command], stdout=log, stderr=subprocess.STDOUT, check=True)
    require(sha(binary) == binary_sha, 'executable changed during measurement')
    report = json.loads((report_dir / 'report.json').read_text())
    require(report['schema'] == 'pgsum-robustness-report-v1', 'unexpected report schema')
    require(report['executable_sha256'] == 'sha256:' + binary_sha, 'native executable identity differs')
    require(report['worker_threads'] == threads and not report['full_term_audit'], 'unexpected measurement settings')
    for name in (*PAYLOADS, 'report.html'):
        require(report['outputs'][name] == 'sha256:' + sha(report_dir / name), 'output fingerprint differs: ' + name)
    result = dict(binary=str(binary), binary_sha256=binary_sha, command=command,
                  measurement=json.loads(timing.read_text()), report_sha256=sha(report_dir / 'report.json'))
    print(label + ' completed: ' + json.dumps(result['measurement']), flush=True)
    return result, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('before', 'after', 'experiment', 'out'):
        parser.add_argument('--' + name, required=True, type=Path)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    require(args.threads > 0, 'threads must be positive')
    before, after, experiment, out = (getattr(args, name).resolve() for name in ('before', 'after', 'experiment', 'out'))
    for path in (before, after, experiment):
        require(path.is_file(), 'missing input: ' + str(path))
    experiment_sha = sha(experiment)
    out.mkdir()  # No reuse of partial or earlier measurements.
    old, old_report = measure(before, experiment, out, 'before', args.threads)
    new, new_report = measure(after, experiment, out, 'after', args.threads)
    require(sha(experiment) == experiment_sha, 'experiment changed during comparison')
    for field in ('experiment', 'input_sha256', 'mask_algorithm', 'scores'):
        require(old_report[field] == new_report[field], 'numerical or input evidence differs: ' + field)
    payloads = {}
    for name in PAYLOADS:
        old_sha, new_sha = sha(out / 'before' / name), sha(out / 'after' / name)
        require(old_sha == new_sha, 'numerical payload differs: ' + name)
        payloads[name] = new_sha
    old_rss, new_rss = (r['measurement']['peak_rss_kib'] for r in (old, new))
    require(old_rss > 0 and new_rss > 0, 'invalid peak memory observation')
    record = dict(schema='pgsum-robustness-memory-comparison-v1', passed=True,
                  verifier_sha256=sha(Path(__file__)), experiment_sha256=experiment_sha,
                  before=old, after=new, numerical_payload_sha256=payloads,
                  peak_rss_reduction_fraction=1-new_rss/old_rss,
                  scenarios=len(old_report['experiment']['scenarios']),
                  scores=[{key: score[key] for key in ('pgs_id','baseline_terms','baseline_term_sha256')}
                          for score in old_report['scores']],
                  notes=['One sequential before/after observation on a shared host; no controlled timing or universal memory claim.',
                         'All native score summaries and numerical payload bytes must agree; this comparison is not an independent source-genotype oracle.'])
    save(out / 'comparison.json', record)
    print('PASS: identical numerical evidence; memory record at ' + str(out / 'comparison.json'), flush=True)


if __name__ == '__main__':
    main()
