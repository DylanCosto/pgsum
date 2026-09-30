#!/usr/bin/env python3
"""Deterministic CLI benchmark and independent integer scoring oracle (Linux/macOS).

Requires a release pgsum binary, zstd, and GNU time (gtime on macOS). No genome downloads or Python packages.
Uses synthetic workloads; these measurements are not real-genome throughput claims.
"""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import time
import sys
from decimal import Decimal


def weight(i):
    # Exact units of 10^-5, with positive, negative and zero weights.
    return ((i * 17) % 101 - 50) * (1 if i % 3 == 0 else 100)


def genotype(i, sample, dense):
    v = (i * 37 + sample * 13) % 100
    if v == 99:
        return 3
    return v % 3 if dense else (1 if v < 3 else 2 if v == 3 else 0)


def generate(root, terms, samples):
    root.mkdir(parents=True, exist_ok=True)
    seq = 'A' * terms
    (root / 'ref.fa').write_text('>chr1\n' + seq + '\n')
    (root / 'ref.fa.fai').write_text(f'chr1\t{terms}\t6\t{terms}\t{terms + 1}\n')
    for name, effect, other in [('ALT', 'G', 'A'), ('REF', 'A', 'G')]:
        with (root / f'{name}.tsv').open('w') as f:
            f.write(f'#pgs_id={name}\n#genome_build=GRCh38\n')
            f.write('chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n')
            for i in range(terms):
                f.write(f'1\t{i+1}\t{effect}\t{other}\t{Decimal(weight(i)).scaleb(-5)}\n')
    header = (f'##fileformat=VCFv4.2\n##contig=<ID=chr1,length={terms}>\n'
              '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
              '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Block depth">\n'
              '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
    expected = {}
    for dense in [False, True]:
        name = 'dense' if dense else 'sparse'
        sums = [[0, 0, 0, 0] for _ in range(samples)]
        with (root / f'{name}.vcf').open('w') as f:
            f.write(header + '\t'.join(f'S{s}' for s in range(samples)) + '\n')
            for i in range(terms):
                codes = [genotype(i, s, dense) for s in range(samples)]
                f.write(f'chr1\t{i+1}\t.\tA\tG\t.\tPASS\t.\tGT\t')
                f.write('\t'.join(['0|0', '0|1', '1|1', './.'][c] for c in codes) + '\n')
                for s, c in enumerate(codes):
                    if c != 3:
                        sums[s][0] += weight(i) * c
                        sums[s][1] += weight(i) * (2-c)
                        sums[s][2] += 1
                        sums[s][3] += abs(weight(i))
        expected[name] = sums
    sums = [0, 0, terms, sum(abs(weight(i)) for i in range(terms))]
    with (root / 'sample.g.vcf').open('w') as f:
        f.write(header + 'S0\n')
        for i in range(0, terms, 4):
            f.write(f'chr1\t{i+1}\t.\tA\tG,<NON_REF>\t.\tPASS\t.\tGT:DP:GQ\t0/1:30:60\n')
            if i + 1 < terms:
                f.write(f'chr1\t{i+2}\t.\tA\t<NON_REF>\t.\tPASS\tEND={min(i+4, terms)}\tGT:MIN_DP:GQ\t0/0:30:60\n')
        for i in range(terms):
            c = int(i % 4 == 0)
            sums[0] += weight(i) * c
            sums[1] += weight(i) * (2-c)
    expected['gvcf'] = [sums]
    return expected


def measured(command, cwd, log):
    timer = shutil.which('gtime') or '/usr/bin/time'
    stats = log.with_suffix('.time')
    start = time.perf_counter()
    with log.open('wb') as output:
        result = subprocess.run([timer, '-f', '%U %S %M', '-o', str(stats),
                                 *[str(x) for x in command]], cwd=cwd, stdout=output, stderr=output)
    if result.returncode:
        raise RuntimeError(f'Command failed ({result.returncode}): {command}\n{log.read_text()}')
    user, system, rss = map(float, stats.read_text().split())
    return dict(wall_seconds=time.perf_counter()-start, cpu_seconds=user+system, peak_rss_mib=rss/1024)



def verify(directory, expected, terms, total_weight, cohort):
    if cohort:
        raw = subprocess.check_output(['zstd', '-dc', str(directory / 'cohort-scores.tsv.zst')])
        rows = list(csv.DictReader(io.StringIO(raw.decode()), delimiter='\t'))
        assert len(rows) == len(expected) * 2
        seen = set()
        for r in rows:
            key = (int(r['sample'][1:]), r['pgs_id'])
            assert key not in seen
            seen.add(key)
            e = expected[key[0]]
            assert Decimal(r['partial_raw_score']) == Decimal(e[key[1] == 'REF']).scaleb(-5)
            assert int(r['scorable_terms']) == e[2] and int(r['total_terms']) == terms
            assert abs(float(r['weight_coverage']) - e[3]/total_weight) <= 0.00000051
            assert abs(float(r['term_coverage']) - e[2]/terms) <= 0.00000051
        return hashlib.sha256(raw).hexdigest()
    results = []
    for idx, name in enumerate(['ALT', 'REF']):
        data = json.loads((directory / f'{name}.score.json').read_text())
        assert Decimal(data['partial']['raw_score']) == Decimal(expected[0][idx]).scaleb(-5)
        assert data['scorable_terms'] == terms
        # Version changes are allowed; all other fields must agree across binaries and thread counts.
        data.pop('pgsum_version', None)
        results.append(data)
    return hashlib.sha256(json.dumps(results, sort_keys=True).encode()).hexdigest()


def main():
    if sys.flags.optimize:
        raise RuntimeError('Run without -O: arithmetic assertions are required')
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, default=Path('target/release/pgsum'))
    p.add_argument('--baseline', type=Path)
    p.add_argument('--output', type=Path, required=True, help='New directory for inputs, logs, and results.json')
    p.add_argument('--terms', type=int, default=20000)
    p.add_argument('--samples', type=int, default=256)
    p.add_argument('--threads', type=int, nargs='+', default=[1, 4])
    p.add_argument('--repeats', type=int, default=3)
    p.add_argument('--max-score-slowdown', type=float, help='Optional baseline gate, e.g. 1.15; use on a quiet dedicated host')
    args = p.parse_args()
    if not shutil.which('zstd'):
        p.error('Requires the zstd command')
    timer = shutil.which('gtime') or '/usr/bin/time'
    check = subprocess.run([timer, '--version'], capture_output=True, text=True)
    if check.returncode or 'GNU' not in check.stdout:
        p.error('Requires GNU time (/usr/bin/time on Linux, gtime on macOS)')
    if min(args.terms, args.samples, args.repeats, *args.threads) < 1:
        p.error('Sizes, repeats and thread counts must be positive')
    if args.max_score_slowdown is not None and (not args.baseline or args.max_score_slowdown <= 0):
        p.error('--max-score-slowdown requires --baseline and a positive ratio')
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binaries = [('candidate', args.binary.resolve())]
    if args.baseline:
        binaries.insert(0, ('baseline', args.baseline.resolve()))
    expected = generate(root / 'inputs', args.terms, args.samples)
    total_weight = sum(abs(weight(i)) for i in range(args.terms))
    if not total_weight:
        p.error('Workload must have nonzero total weight')
    report = dict(platform=platform.platform(), cpu_count=os.cpu_count(),
                  terms=args.terms, samples=args.samples, threads=args.threads, repeats=args.repeats,
                  cache_policy='One unmeasured warm-up per command; OS caches are not flushed',
                  binaries={}, measurements=[],
                  input_sha256={f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                                for f in sorted((root / 'inputs').iterdir())},
                  harness_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    fingerprints = {}
    for label, binary in binaries:
        report['binaries'][label] = dict(path=str(binary), sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                                        version=subprocess.check_output([binary, '--version'], text=True).strip())
        work = root / label
        work.mkdir()
        inputs = root / 'inputs'
        packs = work / 'packs'
        def run(stage, command):
            measured(command, inputs, work / 'warmup.log')
            times = [measured(command, inputs, work / f'{stage}-{r}.log') for r in range(args.repeats)]
            record = dict(binary=label, stage=stage, runs=times,
                          median_seconds=statistics.median(t['wall_seconds'] for t in times),
                          max_rss_mib=max(t['peak_rss_mib'] for t in times))
            report['measurements'].append(record)
            print(f"{label:9} {stage:22} {record['median_seconds']:.4f}s {record['max_rss_mib']:.1f} MiB", flush=True)
        run('compile', [binary, 'compile', 'ALT.tsv', 'REF.tsv', '--reference', 'ref.fa', '--out', packs, '--threads', 1])
        for threads in args.threads:
            for workload in ['gvcf', 'sparse', 'dense']:
                cohort = workload != 'gvcf'
                source = f'{workload}.vcf' if cohort else 'sample.g.vcf'
                table = work / (workload + ('.pgsc' if cohort else '.pgsg'))
                out = work / f'{workload}-t{threads}'
                common = ['--pack', packs, '--threads', threads]
                run(f'{workload}-extract-t{threads}', [binary, 'extract', '--gvcf', source, '--reference', 'ref.fa',
                    '--out', table, *common, *(['--all-samples', '--accept-missing-quality'] if cohort else [])])
                run(f'{workload}-score-t{threads}', [binary, 'score', '--genotypes', table, '--out', out, *common])
                digest = verify(out, expected[workload], args.terms, total_weight, cohort)
                previous = fingerprints.setdefault(workload, digest)
                assert digest == previous, f'{workload}: output changed across binaries/threads'
                if not cohort:
                    run(f'gvcf-run-t{threads}', [binary, 'run', '--gvcf', source, '--reference', 'ref.fa', '--out', out, *common])
                    assert verify(out, expected[workload], args.terms, total_weight, False) == digest
    report['verified_output_sha256'] = fingerprints
    report['oracle'] = 'All sample sums and cohort coverage checked against independent integer arithmetic'
    if args.baseline:
        by_key = {(r['binary'], r['stage']): r for r in report['measurements']}
        report['speedup_baseline_over_candidate'] = {
            stage: by_key[('baseline', stage)]['median_seconds'] / by_key[('candidate', stage)]['median_seconds']
            for label, stage in by_key if label == 'candidate'}
    (root / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
    print(f'Exact-output checks passed; results: {root / "results.json"}')
    if args.max_score_slowdown is not None:
        regressions = [stage for stage, speedup in report['speedup_baseline_over_candidate'].items()
                       if '-score-' in stage and 1 / speedup > args.max_score_slowdown]
        if regressions:
            raise SystemExit(f'Scoring slowdown exceeded {args.max_score_slowdown}: {regressions}')


if __name__ == '__main__':
    main()
