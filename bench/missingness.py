#!/usr/bin/env python3
"""Compare missingness-report overhead and verify bounds with an independent additive-score oracle."""
import argparse
import hashlib
import json
import platform
import os
from pathlib import Path
import statistics
from decimal import Decimal
from performance import generate, measured, weight


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', required=True, type=Path)
    p.add_argument('--baseline', type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--terms', type=int, default=100000)
    p.add_argument('--threads', nargs='+', type=int, default=[1, 4])
    p.add_argument('--repeats', type=int, default=3)
    args = p.parse_args()
    if min([args.terms, args.repeats, *args.threads]) < 1:
        p.error('counts must be positive')
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    inputs = root / 'inputs'
    generate(inputs, args.terms, 1)
    header = '##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\n'
    expected = []
    for missing in [False, True]:
        name = 'half-missing' if missing else 'complete'
        sums = [0, 0]
        lower = upper = count = 0
        with (inputs / f'{name}.vcf').open('w') as f:
            f.write(header)
            for i in range(args.terms):
                absent = missing and i % 2 == 0
                genotype = './.' if absent else '0/1'
                f.write(f'chr1\t{i+1}\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t{genotype}:30:60\n')
                if absent:
                    count += 1
                    lower += min(0, 2*weight(i))
                    upper += max(0, 2*weight(i))
                else:
                    sums[0] += weight(i)
                    sums[1] += weight(i)
        expected.append((name, sums, count, lower, upper))
    binaries = [('candidate', args.binary.resolve())]
    if args.baseline:
        binaries.insert(0, ('baseline', args.baseline.resolve()))
    report = {'platform': platform.platform(), 'python': platform.python_version(), 'logical_cpus': os.cpu_count(),
              'terms_per_score': args.terms, 'scores': 2, 'repeats': args.repeats,
              'threads': args.threads, 'measurements': [], 'binaries': {},
              'oracle': 'Independent integer additive contributions and missing-term extrema; exact legacy result parity'}
    fingerprints = {}
    diagnostics = {}
    for label, binary in binaries:
        report['binaries'][label] = {'path': str(binary), 'sha256': hashlib.sha256(binary.read_bytes()).hexdigest()}
        work = root / label
        work.mkdir()
        packs = work / 'packs'
        measured([binary, 'compile', 'ALT.tsv', 'REF.tsv', '--reference', 'ref.fa', '--out', packs, '--threads', 1], inputs, work / 'compile.log')
        for threads in args.threads:
            for name, sums, count, lower, upper in expected:
                out = work / f'{name}-{threads}'
                for stage in ['run', 'score']:
                    command = [binary, stage, '--pack', packs, '--out', out, '--threads', threads]
                    command += (['--gvcf', f'{name}.vcf', '--reference', 'ref.fa'] if stage == 'run'
                                else ['--genotypes', out / 'genotypes.pgsg'])
                    measured(command, inputs, work / f'{name}-{threads}-{stage}-warmup.log')
                    runs = [measured(command, inputs, work / f'{name}-{threads}-{stage}-{i}.log') for i in range(args.repeats)]
                    row = dict(binary=label, workload=name, stage=stage, threads=threads, runs=runs,
                               median_seconds=statistics.median(r['wall_seconds'] for r in runs),
                               max_rss_mib=max(r['peak_rss_mib'] for r in runs))
                    report['measurements'].append(row)
                    print(f"{label} {name} {stage} {threads}: {row['median_seconds']:.4f}s {row['max_rss_mib']:.1f} MiB", flush=True)
                    for score, expected_sum in zip(['ALT', 'REF'], sums):
                        result = json.loads((out / f'{score}.score.json').read_text())
                        assert Decimal(result['partial']['raw_score']) == Decimal(expected_sum).scaleb(-5)
                        assert result['scorable_terms'] == args.terms - count
                        missingness = result.pop('missingness', None)
                        if label == 'candidate':
                            assert missingness['missing_terms'] == count
                            assert missingness['unbounded_terms'] == 0
                            for key, value in [('lower', lower), ('upper', upper)]:
                                assert Decimal(missingness['bounded_missing_contribution'][key]) == Decimal(value).scaleb(-5)
                                assert Decimal(missingness['completion_score_bounds'][key]) == Decimal(expected_sum+value).scaleb(-5)
                            key = (name, score)
                            digest = hashlib.sha256(json.dumps(missingness, sort_keys=True).encode()).hexdigest()
                            assert diagnostics.setdefault(key, digest) == digest
                        result.pop('pgsum_version', None)
                        digest = hashlib.sha256(json.dumps(result, sort_keys=True).encode()).hexdigest()
                        key = (name, score)
                        assert fingerprints.setdefault(key, digest) == digest
    report['verified_legacy_results'] = {':'.join(k): v for k, v in fingerprints.items()}
    report['verified_missingness_reports'] = {':'.join(k): v for k, v in diagnostics.items()}
    (root / 'results.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    import sys
    if sys.flags.optimize:
        raise RuntimeError('Run without -O: oracle assertions are required')
    main()
