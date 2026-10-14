#!/usr/bin/env python3
"""Synthetic 22-autosome robustness scaling with an independent exact integer oracle.

No genome downloads or Python packages. Requires a prebuilt pgSum binary and GNU time.
Synthetic timings and numerical reconstruction are not predictive-performance evidence.
"""

import argparse
import bisect
import csv
from fractions import Fraction
import hashlib
import heapq
import json
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import struct
import subprocess
import sys
import time


MODEL = 'SYNTHETIC_DENSE'
SCHEMA = 'pgsum-robustness-scale-v1'
GENERATOR = 'pgsum-robustness-scale-generator-v1'
MASK_DOMAIN = b'pgsum-robustness-mask-sha256-v1\0'
WEIGHT_SCALE = 10000
MISSING_MODULUS = 1000003
TOP_TERMS = 5
CUTOFFS = [50.0, 90.0]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def decimal_units(units):
    """Render integer 10^-4 units without floating-point arithmetic."""
    return f'{"-" if units < 0 else ""}{abs(units) // WEIGHT_SCALE}.{abs(units) % WEIGHT_SCALE:04d}'


def weight(index):
    return (index * 17) % 101 - 50


def dosage(index, sample):
    return ((index * 37 + sample * 13) % 101) % 3


def missing(index, sample, threshold):
    """Deterministic synthetic cell pattern, not an empirical or independent-error model.

    A prime-field arithmetic grid avoids retaining random state; cohort prefixes are nested.
    The threshold is floor(rate * 1,000,003), so actual missingness is measured, not assumed.
    """
    return threshold > 0 and (
        (index + 1) * 2654435761 + (sample + 1) * 2246822519
        + (index + 1) * (sample + 1) * 3266489917
    ) % MISSING_MODULUS < threshold


def loci(terms):
    """Genome-order global ordinal, autosome, one-based position; balanced across 22 contigs."""
    index = 0
    for contig in range(1, 23):
        count = terms // 22 + int(contig <= terms % 22)
        for position in range(1, count + 1):
            yield index, contig, position
            index += 1


def random_draw(seed, contig, position, sample, unit):
    encoded = MASK_DOMAIN + struct.pack('<QBI', seed, contig, position)
    encoded += b'\0' if unit == 'site' else b'\1'
    if unit == 'sample_site':
        name = sample.encode('utf-8')
        encoded += struct.pack('<Q', len(name)) + name
    return (int.from_bytes(hashlib.sha256(encoded).digest()[:8], 'little') >> 11) / (1 << 53)


def generate_model(root, terms):
    """Stream the source model and a valid FASTA/FAI, retaining no genome-sized Python list."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=False)
    with (root / 'reference.fa').open('wb') as fasta, (root / 'reference.fa.fai').open('w') as fai:
        for contig in range(1, 23):
            count = terms // 22 + int(contig <= terms % 22)
            fasta.write(f'>chr{contig}\n'.encode())
            offset = fasta.tell()
            width = min(count, 80)
            fai.write(f'chr{contig}\t{count}\t{offset}\t{width}\t{width + 1}\n')
            for start in range(0, count, 80):
                fasta.write(b'A' * min(80, count - start) + b'\n')
    with (root / 'model.tsv').open('w') as score:
        score.write(f'#pgs_id={MODEL}\n#genome_build=GRCh38\n#variants_number={terms}\n')
        score.write('chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n')
        for index, contig, position in loci(terms):
            effect, other = ('G', 'A') if index % 2 == 0 else ('A', 'G')
            score.write(f'{contig}\t{position}\t{effect}\t{other}\t{decimal_units(weight(index))}\n')
    with (root / 'chr9.bed').open('w') as bed:
        bed.write(f'chr9\t0\t{terms // 22 + int(9 <= terms % 22)}\n')


def generate_cohort(root, terms, samples, missing_rate):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=False)
    threshold = int(Fraction(str(missing_rate)) * MISSING_MODULUS)
    missing_cells = 0
    with (root / 'cohort.vcf').open('w') as vcf:
        vcf.write('##fileformat=VCFv4.3\n')
        for contig in range(1, 23):
            count = terms // 22 + int(contig <= terms % 22)
            vcf.write(f'##contig=<ID=chr{contig},length={count}>\n')
        vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Synthetic diploid hard call">\n')
        vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
        vcf.write('\t'.join(f'S{s:06d}' for s in range(samples)) + '\n')
        for index, contig, position in loci(terms):
            calls = []
            for sample in range(samples):
                absent = missing(index, sample, threshold)
                missing_cells += absent
                calls.append('./.' if absent else ('0/0', '0/1', '1/1')[dosage(index, sample)])
            vcf.write(f'chr{contig}\t{position}\t.\tA\tG\t.\tPASS\t.\tGT\t' + '\t'.join(calls) + '\n')
    return {'missing_cells': missing_cells, 'genotype_cells': terms * samples,
            'actual_missing_fraction': missing_cells / (terms * samples), 'missing_threshold': threshold}


def make_manifest(model_root, case_root, seed, unit):
    model_root, case_root = Path(model_root), Path(case_root)
    return {
        'schema': 'pgsum-robustness-experiment-v1',
        'genotypes': str(case_root / 'cohort.pgsc'),
        'packs': [str(model_root / 'packs' / f'{MODEL}.pgsp')],
        'cutoffs': CUTOFFS, 'top_terms': TOP_TERMS,
        'scenarios': [
            {'id': 'control'},
            {'id': 'random20', 'masks': [{'kind': 'random', 'rate': 0.2, 'seed': seed, 'unit': unit}]},
            {'id': 'chr9-dropout', 'masks': [{'kind': 'bed', 'path': str(model_root / 'chr9.bed'), 'mode': 'drop'}]},
        ],
        'regions': [{'id': 'chr9', 'path': str(model_root / 'chr9.bed')}],
    }


def audit_update(digest, scenario, sample, ordinal, contig, position, base, changed):
    reason = 'masked'
    row = [MODEL, scenario, f'S{sample:06d}', ordinal, contig, position, reason, base, changed, changed - base]
    digest.update(('\t'.join(map(str, row)) + '\n').encode())


def expected_results(terms, samples, missing_rate, seed, unit, audit):
    """Formula-only integer oracle: never reads packs, extracted tables or native contributions.

    Only one row of calls and bounded top-K heaps are retained. A million-term run does not
    construct a sample-by-term or sample-by-scenario-by-term Python object graph.
    """
    threshold = int(Fraction(str(missing_rate)) * MISSING_MODULUS)
    baseline = [0] * samples
    scenarios = {
        name: {'delta': [0] * samples, 'masked': [0] * samples, 'changed': [0] * samples,
               'region': [0] * samples, 'top': [[] for _ in range(samples)]}
        for name in ['control', 'random20', 'chr9-dropout']
    }
    eligible = baseline_weight = total_weight = decisions = changed_terms = 0
    universe = hashlib.sha256()
    audit_digest = hashlib.sha256()
    by_contig = {str(contig): 0 for contig in range(1, 23)}
    for index, contig, position in loci(terms):
        w = weight(index)
        total_weight += abs(w)
        if threshold and any(missing(index, sample, threshold) for sample in range(samples)):
            continue
        eligible += 1
        baseline_weight += abs(w)
        by_contig[str(contig)] += 1
        universe.update(struct.pack('<QBI', index + 1, contig, position))
        values = [w * (dosage(index, sample) if index % 2 == 0 else 2 - dosage(index, sample))
                  for sample in range(samples)]
        for sample, value in enumerate(values):
            baseline[sample] += value
        site_hidden = random_draw(seed, contig, position, '', unit) < 0.2 if unit == 'site' else None
        for name in ['random20', 'chr9-dropout']:
            arm = scenarios[name]
            for sample, value in enumerate(values):
                hidden = (contig == 9 if name == 'chr9-dropout' else site_hidden if unit == 'site'
                          else random_draw(seed, contig, position, f'S{sample:06d}', unit) < 0.2)
                decisions += 1
                if not hidden:
                    continue
                arm['masked'][sample] += 1
                arm['delta'][sample] -= value
                if contig == 9:
                    arm['region'][sample] -= value
                if not value:
                    continue
                changed_terms += 1
                arm['changed'][sample] += 1
                # Largest absolute delta first, ties broken by the original model ordinal.
                entry = (abs(value), -(index + 1), -value, contig, position, value)
                heap = arm['top'][sample]
                if len(heap) < TOP_TERMS:
                    heapq.heappush(heap, entry)
                elif entry[:2] > heap[0][:2]:
                    heapq.heapreplace(heap, entry)
                if audit:
                    audit_update(audit_digest, name, sample, index + 1, contig, position, value, 0)
    if not eligible:
        raise RuntimeError('Synthetic case has no common-callable terms; reduce missingness/cohort size')
    for arm in scenarios.values():
        for sample in range(samples):
            arm['top'][sample] = sorted(arm['top'][sample], key=lambda item: (-item[0], -item[1]))
    return {
        'baseline': baseline, 'scenarios': scenarios, 'eligible': eligible,
        'baseline_weight': baseline_weight, 'total_weight': total_weight,
        'baseline_term_sha256': universe.hexdigest(), 'common_terms_by_contig': by_contig,
        'mask_decisions_checked': decisions, 'nonzero_changes': changed_terms,
        'full_audit_sha256': audit_digest.hexdigest() if audit else None,
    }


def score_units(value):
    exact = Fraction(value) * WEIGHT_SCALE
    assert exact.denominator == 1, f'Unexpected precision: {value}'
    return exact.numerator


def midpoint_percentiles(values, reference):
    if len(reference) < 2 or min(reference) == max(reference):
        return [None] * len(values)
    ordered = sorted(reference)
    return [Fraction(50 * (bisect.bisect_left(ordered, value) + bisect.bisect_right(ordered, value)), len(ordered))
            for value in values]


def pearson(x, y):
    if len(x) < 2:
        return None
    n = len(x)
    sx, sy = sum(x), sum(y)
    xx = n * sum(value * value for value in x) - sx * sx
    yy = n * sum(value * value for value in y) - sy * sy
    xy = n * sum(a * b for a, b in zip(x, y)) - sx * sy
    return float(xy) / math.sqrt(float(xx * yy)) if xx and yy else None


def doubled_ranks(values):
    ordered = sorted(values)
    return [bisect.bisect_left(ordered, value) + bisect.bisect_right(ordered, value) + 1 for value in values]


def close(actual, expected):
    if expected is None:
        assert actual is None, (actual, expected)
    else:
        assert actual is not None and math.isfinite(float(actual)), (actual, expected)
        assert math.isclose(float(actual), float(expected), rel_tol=1e-10, abs_tol=1e-12), (actual, expected)


def tsv_rows(path):
    with Path(path).open() as source:
        yield from csv.DictReader(source, delimiter='\t')


def validate_metrics(metrics, baseline, changed):
    n = len(baseline)
    differences = [b - a for a, b in zip(baseline, changed)]
    assert metrics['n'] == n
    close(metrics['mae'], Fraction(sum(abs(value) for value in differences), n * WEIGHT_SCALE))
    close(metrics['rmse'], math.sqrt(sum(value * value for value in differences) / n) / WEIGHT_SCALE)
    close(metrics['max_abs_change'], Fraction(max(map(abs, differences)), WEIGHT_SCALE))
    close(metrics['pearson'], pearson(baseline, changed))
    close(metrics['spearman'], pearson(doubled_ranks(baseline), doubled_ranks(changed)))
    before, after = midpoint_percentiles(baseline, baseline), midpoint_percentiles(changed, baseline)
    shifts = [abs(b - a) for a, b in zip(before, after) if a is not None and b is not None]
    close(metrics['mean_abs_percentile_shift'], sum(shifts) / n if shifts else None)
    close(metrics['max_abs_percentile_shift'], max(shifts) if shifts else None)
    assert [crossing['cutoff'] for crossing in metrics['crossings']] == CUTOFFS
    for crossing in metrics['crossings']:
        threshold = Fraction(str(crossing['cutoff']))
        pairs = [(a >= threshold, b >= threshold) for a, b in zip(before, after)
                 if a is not None and b is not None]
        assert crossing['evaluable'] == len(pairs)
        assert crossing['up'] == sum(not a and b for a, b in pairs)
        assert crossing['down'] == sum(a and not b for a, b in pairs)
        assert crossing['unchanged'] == sum(a == b for a, b in pairs)


def verify(root, expected, terms, samples, audit, binary_hash):
    root = Path(root)
    report = json.loads((root / 'report.json').read_text())
    assert report['executable_sha256'].removeprefix('sha256:') == binary_hash
    assert len(report['scores']) == 1
    score = report['scores'][0]
    assert score['pgs_id'] == MODEL and score['total_model_terms'] == terms
    assert score['baseline_terms'] == expected['eligible']
    assert score['baseline_term_sha256'] == expected['baseline_term_sha256']
    close(score['baseline_term_fraction'], Fraction(expected['eligible'], terms))
    close(score['baseline_weight_fraction'], Fraction(expected['baseline_weight'], expected['total_weight']))
    excluded = terms - expected['eligible']
    assert score['excluded_terms'] == ({'baseline_target_not_common_callable': excluded} if excluded else {})
    assert score['inventory_consistent']
    assert [arm['id'] for arm in score['scenarios']] == list(expected['scenarios'])
    for arm in score['scenarios']:
        oracle = expected['scenarios'][arm['id']]
        assert arm['masked_observations'] == sum(oracle['masked'])
        assert arm['unavailable_alternate_observations'] == arm['filled_observations'] == 0
        assert arm['unfilled_observations'] == sum(oracle['masked'])
        changed = [base + delta for base, delta in zip(expected['baseline'], oracle['delta'])]
        validate_metrics(arm['metrics'], expected['baseline'], changed)
    seen = set()
    for row in tsv_rows(root / 'samples.tsv'):
        name, sample = row['scenario'], int(row['sample'][1:])
        assert row['sample'] == f'S{sample:06d}' and 0 <= sample < samples
        assert (name, sample) not in seen
        seen.add((name, sample))
        arm = expected['scenarios'][name]
        assert row['pgs_id'] == MODEL
        assert score_units(row['baseline_score']) == expected['baseline'][sample]
        assert score_units(row['scenario_score']) == expected['baseline'][sample] + arm['delta'][sample]
        assert score_units(row['delta']) == arm['delta'][sample]
        assert int(row['baseline_terms']) == expected['eligible']
        assert int(row['masked_terms']) == int(row['unfilled_terms']) == arm['masked'][sample]
        assert int(row['retained_terms']) == expected['eligible'] - arm['masked'][sample]
        assert int(row['unavailable_alternate_terms']) == int(row['filled_terms']) == 0
        assert int(row['changed_terms']) == arm['changed'][sample]
        top = sum(item[2] for item in arm['top'][sample])
        assert score_units(row['top_terms_delta']) == top
        assert score_units(row['other_terms_delta']) == arm['delta'][sample] - top
        changed = [base + delta for base, delta in zip(expected['baseline'], arm['delta'])]
        for column, percentiles in [
            ('baseline_percentile', midpoint_percentiles(expected['baseline'], expected['baseline'])),
            ('fixed_baseline_percentile', midpoint_percentiles(changed, expected['baseline'])),
            ('reranked_percentile', midpoint_percentiles(changed, changed)),
        ]:
            close(float(row[column]) if row[column] else None, percentiles[sample])
    assert len(seen) == 3 * samples
    seen.clear()
    with (root / 'top-changes.jsonl').open() as source:
        for line in source:
            row = json.loads(line)
            name, sample = row['scenario'], int(row['sample'][1:])
            assert row['pgs_id'] == MODEL and (name, sample) not in seen
            seen.add((name, sample))
            arm = expected['scenarios'][name]
            assert row['changed_terms'] == arm['changed'][sample]
            assert len(row['top']) == len(arm['top'][sample])
            for actual, wanted in zip(row['top'], arm['top'][sample]):
                _, negative_ordinal, delta, contig, position, base = wanted
                assert actual['ordinal'] == -negative_ordinal
                assert actual['contig'] == f'chr{contig}'
                assert actual['pos1'] == actual['source_pos1'] == position
                assert actual['reason'] == 'masked'
                assert score_units(actual['baseline']) == base
                assert score_units(actual['scenario']) == 0
                assert score_units(actual['delta']) == delta
            top = sum(item[2] for item in arm['top'][sample])
            assert score_units(row['top_terms_delta']) == top
            assert score_units(row['other_terms_delta']) == arm['delta'][sample] - top
            assert score_units(row['delta']) == arm['delta'][sample]
    assert len(seen) == 3 * samples
    seen.clear()
    for row in tsv_rows(root / 'regions.tsv'):
        name, sample = row['scenario'], int(row['sample'][1:])
        assert row['pgs_id'] == MODEL and row['region'] == 'chr9' and (name, sample) not in seen
        seen.add((name, sample))
        assert score_units(row['delta']) == expected['scenarios'][name]['region'][sample]
    assert len(seen) == 3 * samples
    assert list(tsv_rows(root / 'references.tsv')) == []
    if audit:
        digest, count = hashlib.sha256(), 0
        process = subprocess.Popen(['zstd', '-dc', str(root / 'terms.tsv.zst')], stdout=subprocess.PIPE, text=True)
        try:
            for row in csv.DictReader(process.stdout, delimiter='\t'):
                assert row['pgs_id'] == MODEL and row['reason'] == 'masked'
                base, changed = score_units(row['baseline']), score_units(row['scenario_contribution'])
                assert score_units(row['delta']) == changed - base != 0
                audit_update(digest, row['scenario'], int(row['sample'][1:]), int(row['ordinal']),
                             int(row['contig'].removeprefix('chr')), int(row['pos1']), base, changed)
                count += 1
        except BaseException:
            process.stdout.close()
            if process.poll() is None:
                process.terminate()
            process.wait()
            raise
        process.stdout.close()
        status = process.wait()
        assert status == 0 and count == expected['nonzero_changes']
        assert digest.hexdigest() == expected['full_audit_sha256']
    for name, digest in report['outputs'].items():
        assert sha256(root / name) == digest.removeprefix('sha256:')
    for name, digest in report['input_sha256'].items():
        assert sha256(name) == digest.removeprefix('sha256:')
    # These files contain numerical evidence, not run paths, thread counts, timing or versions.
    evidence = {name: sha256(root / name) for name in report['outputs'] if name != 'report.html'}
    evidence['summaries'] = hashlib.sha256(json.dumps(score, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return evidence


def measured(command, cwd, log, timer):
    command = list(map(str, command))
    stats = log.with_suffix('.time')
    started = time.perf_counter()
    with log.open('wb') as output:
        result = subprocess.run([timer, '-f', '%e %U %S %M', '-o', str(stats), *command],
                                cwd=cwd, stdout=output, stderr=output)
    if result.returncode:
        raise RuntimeError(f'Command failed ({result.returncode}); inspect {log}: {command}')
    elapsed, user, system, rss = map(float, stats.read_text().split())
    return {'command': command, 'cwd': str(cwd), 'log': str(log), 'wall_seconds': time.perf_counter() - started,
            'gnu_time_elapsed_seconds': elapsed, 'cpu_seconds': user + system, 'peak_rss_mib': rss / 1024,
            'measurement_scope': 'native child process only; Python generation and oracle excluded'}


def protocol_self_test():
    # Independent published vectors pinned in the robustness mask specification/tests.
    assert random_draw(12345, 1, 1000, 'ignored', 'site') == 3076170996652624 / (1 << 53)
    assert random_draw(12345, 1, 1000, 'NA12878', 'sample_site') == 3112782897755573 / (1 << 53)
    assert random_draw(12345, 1, 1000, 'sample α', 'sample_site') == 2833699031287665 / (1 << 53)
    assert decimal_units(-51) == '-0.0051' and score_units('-0.0051000') == -51
    assert midpoint_percentiles([0, 1, 2, 3], [1, 1, 2]) == [Fraction(0), Fraction(100, 3), Fraction(250, 3), Fraction(100)]


def environment(timer_version):
    result = {'platform': platform.platform(), 'python': platform.python_version(),
              'logical_cpus': os.cpu_count(), 'machine': platform.machine(),
              'gnu_time': timer_version}
    cpuinfo = Path('/proc/cpuinfo')
    if cpuinfo.is_file():
        result['cpu_model'] = next((line.split(':', 1)[1].strip() for line in cpuinfo.read_text().splitlines()
                                    if line.startswith('model name')), None)
    for name, command in [('available_rustc', ['rustc', '--version', '--verbose']),
                          ('checkout_revision', ['git', 'rev-parse', 'HEAD'])]:
        if shutil.which(command[0]):
            output = subprocess.run(command, capture_output=True, text=True)
            if output.returncode == 0:
                result[name] = output.stdout.strip()
    return result


def main():
    if sys.flags.optimize:
        raise RuntimeError('Run without -O: exact-arithmetic and integrity assertions are required')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=Path('target/release/pgsum'))
    parser.add_argument('--output', type=Path, required=True, help='New benchmark directory; never overwritten')
    parser.add_argument('--terms', type=int, nargs='+', default=[10000, 100000], help='Terms per model; Cartesian product with samples/rates')
    parser.add_argument('--samples', type=int, nargs='+', default=[32, 128], help='Nested cohort sizes')
    parser.add_argument('--missing-rates', type=float, nargs='+', default=[0.0, 0.005])
    parser.add_argument('--threads', type=int, nargs='+', default=[1, 4])
    parser.add_argument('--repeats', type=int, default=1, help='Measured fresh-directory repetitions per thread count; no warm-up')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--mask-unit', choices=['site', 'sample_site'], default='site')
    parser.add_argument('--full-audit', action='store_true', help='Write and independently stream-check every changed contribution; requires zstd')
    parser.add_argument('--keep-inputs', action='store_true', help='Retain synthetic VCFs, packs, FASTA and extracted cohorts; default deletes each case after verification')
    args = parser.parse_args()
    for name in ['terms', 'samples', 'threads', 'missing_rates']:
        values = getattr(args, name)
        if len(values) != len(set(values)):
            parser.error(f'--{name.replace("_", "-")} must not contain duplicates')
    if min(args.terms) < 22 or min(args.samples + args.threads + [args.repeats]) < 1:
        parser.error('At least 22 terms and positive sample/thread/repeat counts are required')
    if any(not math.isfinite(rate) or not 0 <= rate < 1 for rate in args.missing_rates):
        parser.error('Missing rates must be finite and in [0,1)')
    if not 0 <= args.seed < 1 << 64:
        parser.error('Seed must be a u64')
    if args.full_audit and not shutil.which('zstd'):
        parser.error('--full-audit requires the zstd command')
    timer = shutil.which('gtime') or '/usr/bin/time'
    check = subprocess.run([timer, '--version'], capture_output=True, text=True)
    if check.returncode or 'GNU' not in check.stdout:
        parser.error('GNU time is required (/usr/bin/time on Linux, gtime on macOS)')
    binary = args.binary.resolve(strict=True)
    protocol_self_test()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary_hash = sha256(binary)
    report = {
        'schema': SCHEMA, 'generator': GENERATOR, 'passed': False,
        'scope': 'Synthetic mixed REF/ALT additive diploid GT model on 22 autosomes; score reconstruction only. No phenotype, portability or clinical-performance validation.',
        'oracle': 'Integer 10^-4-unit source formulas, original global ordinals, independently reproduced mask protocol. No engine pack/table/contribution as oracle input.',
        'command': [sys.executable, *sys.argv], 'cwd': str(Path.cwd()),
        'harness_sha256': sha256(Path(__file__)),
        'environment': environment(check.stdout.splitlines()[0]),
        'binary': {'path': str(binary), 'sha256': binary_hash,
                   'version': subprocess.check_output([binary, '--version'], text=True).strip()},
        'parameters': {'terms': args.terms, 'samples': args.samples, 'missing_rates': args.missing_rates,
                       'threads': args.threads, 'repeats': args.repeats, 'seed': args.seed,
                       'mask_unit': args.mask_unit, 'full_audit': args.full_audit, 'keep_inputs': args.keep_inputs},
        'case_count': len(args.terms) * len(args.samples) * len(args.missing_rates),
        'native_robustness_runs': len(args.terms) * len(args.samples) * len(args.missing_rates) * len(args.threads) * args.repeats,
        'measurement_policy': 'No warm-up; no OS-cache flush. Every native invocation uses GNU time in its own child process; generation, oracle construction and verification are timed separately. Shared-runner wall times are observations, not performance gates.',
        'retention_policy': 'Logs, manifests, summaries, output evidence and input hashes retained. Without --keep-inputs, generated VCF/cohort and model inputs are deleted only after verification; regenerate from pinned harness/parameters.',
        'models': [], 'cases': [],
    }
    for terms in sorted(args.terms):
        model_root = root / f'model-{terms}'
        started = time.perf_counter()
        generate_model(model_root, terms)
        model = {'terms': terms, 'generation_seconds': time.perf_counter() - started,
                 'input_sha256': {name: sha256(model_root / name) for name in ['model.tsv', 'reference.fa', 'reference.fa.fai', 'chr9.bed']}}
        model['compile'] = measured(
            [binary, 'compile', model_root / 'model.tsv', '--reference', model_root / 'reference.fa',
             '--out', model_root / 'packs', '--threads', 1], root, model_root / 'compile.log', timer)
        model['pack_sha256'] = sha256(model_root / 'packs' / f'{MODEL}.pgsp')
        report['models'].append(model)
        eligible_by_case = {}
        for samples in sorted(args.samples):
            for rate in sorted(args.missing_rates):
                case_root = root / f'case-{terms}-{samples}-missing-{rate:g}'
                started = time.perf_counter()
                input_stats = generate_cohort(case_root, terms, samples, rate)
                case = {'terms': terms, 'samples': samples, 'missing_rate': rate, 'generation_seconds': time.perf_counter() - started,
                        'source': input_stats, 'source_vcf_sha256': sha256(case_root / 'cohort.vcf'),
                        'source_vcf_bytes': (case_root / 'cohort.vcf').stat().st_size, 'runs': []}
                case['extract'] = measured(
                    [binary, 'extract', '--gvcf', case_root / 'cohort.vcf', '--reference', model_root / 'reference.fa',
                     '--pack', model_root / 'packs', '--out', case_root / 'cohort.pgsc', '--all-samples',
                     '--accept-missing-quality', '--threads', 1], root, case_root / 'extract.log', timer)
                case['cohort_sha256'] = sha256(case_root / 'cohort.pgsc')
                case['cohort_bytes'] = (case_root / 'cohort.pgsc').stat().st_size
                manifest_path = case_root / 'experiment.json'
                manifest_path.write_text(json.dumps(make_manifest(model_root, case_root, args.seed, args.mask_unit), indent=2) + '\n')
                started = time.perf_counter()
                expected = expected_results(terms, samples, rate, args.seed, args.mask_unit, args.full_audit)
                case['oracle_seconds'] = time.perf_counter() - started
                case['oracle_checks'] = {key: expected[key] for key in ['eligible', 'baseline_term_sha256', 'common_terms_by_contig', 'mask_decisions_checked', 'nonzero_changes', 'full_audit_sha256']}
                case['oracle_checks']['paired_scores'] = 3 * samples
                case['shared_term_fraction'] = expected['eligible'] / terms
                eligible_by_case[(samples, rate)] = expected['eligible']
                fingerprint = None
                for threads in args.threads:
                    for repeat in range(args.repeats):
                        out = case_root / f'run-t{threads}-r{repeat}'
                        command = [binary, 'robustness', '--experiment', manifest_path, '--out', out, '--threads', threads]
                        if args.full_audit:
                            command.append('--terms')
                        measurement = measured(command, root, case_root / f'run-t{threads}-r{repeat}.log', timer)
                        measurement.update(threads=threads, repeat=repeat)
                        started = time.perf_counter()
                        evidence = verify(out, expected, terms, samples, args.full_audit, binary_hash)
                        measurement['verification_seconds'] = time.perf_counter() - started
                        measurement['evidence_sha256'] = evidence
                        if fingerprint is None:
                            fingerprint = evidence
                        assert evidence == fingerprint, 'Numerical evidence changed across repeats/thread counts'
                        case['runs'].append(measurement)
                        print(f'{terms} terms / {samples} samples / missing={rate:g} / t{threads} r{repeat}: '
                              f'{measurement["wall_seconds"]:.3f}s {measurement["peak_rss_mib"]:.1f} MiB; '
                              f'{expected["eligible"]}/{terms} common terms; oracle passed', flush=True)
                case['summary_by_threads'] = [
                    {'threads': threads,
                     'median_wall_seconds': statistics.median(run['wall_seconds'] for run in case['runs'] if run['threads'] == threads),
                     'max_peak_rss_mib': max(run['peak_rss_mib'] for run in case['runs'] if run['threads'] == threads)}
                    for threads in args.threads
                ]
                report['cases'].append(case)
                assert sha256(binary) == binary_hash, 'Native binary changed during the benchmark'
                if not args.keep_inputs:
                    (case_root / 'cohort.vcf').unlink()
                    (case_root / 'cohort.pgsc').unlink()
                (root / 'progress.json').write_text(json.dumps(report, indent=2) + '\n')
        # Nested sample prefixes and missingness thresholds must never increase the common universe.
        for rate in sorted(args.missing_rates):
            counts = [eligible_by_case[(samples, rate)] for samples in sorted(args.samples)]
            assert counts == sorted(counts, reverse=True)
        for samples in sorted(args.samples):
            counts = [eligible_by_case[(samples, rate)] for rate in sorted(args.missing_rates)]
            assert counts == sorted(counts, reverse=True)
        if not args.keep_inputs:
            for name in ['reference.fa', 'reference.fa.fai', 'model.tsv']:
                (model_root / name).unlink()
            shutil.rmtree(model_root / 'packs')
    assert sha256(binary) == binary_hash, 'Native binary changed during the benchmark'
    assert sha256(Path(__file__)) == report['harness_sha256'], 'Benchmark harness changed during the run'
    report['passed'] = True
    temporary = root / 'results.json.tmp'
    temporary.write_text(json.dumps(report, indent=2) + '\n')
    temporary.replace(root / 'results.json')
    (root / 'progress.json').unlink(missing_ok=True)
    print(f'Verified {report["native_robustness_runs"]} native runs; {root / "results.json"}', flush=True)


if __name__ == '__main__':
    main()
