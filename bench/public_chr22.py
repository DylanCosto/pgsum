#!/usr/bin/env python3
"""Reproducible chromosome-22 comparison on public 1000 Genomes data.

Preparation preserves all genotype records for 32 deterministically selected samples,
and every chr22 term from two published harmonized scoring files. Outputs are partial
chromosome contributions, NOT complete published PGS values.
"""
import argparse
import csv
import collections
import bisect
import io
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import urllib.request
from decimal import Decimal, getcontext

import pysam
import zstandard

SCORES = ['PGS001229', 'PGS000018']
PIPELINE_COMMIT = '72ee54f2119d2fe65ac671327342dfc412bd0a3b'
SOURCE_SHA256 = {'chr22.vcf.gz': 'dc2b718c00d1a8220f19561bc647eea08fc25ce641c8f2c12e540fb31f073f88', 'chr22.fa.gz': '05f9d97d6fbfd08a44ca45b50837ca2ae9c471f35ba79dffec04d2cb5eaaf695', 'PGS001229.txt.gz': '35620a9b1c6c640698f85c84a183d130e6c56ef61d697c31a7c9590ff133ef58', 'PGS000018.txt.gz': 'ed5dcd5359213364a182dd7972edf524cc3e161f8332ca54ca5515c95cbf8d2a'}

URLS = {
    'chr22.vcf.gz': 'https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz',
    'chr22.fa.gz': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz',
    **{s + '.txt.gz': f'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{s}/ScoringFiles/Harmonized/{s}_hmPOS_GRCh38.txt.gz' for s in SCORES},
}


def sha(path):
    with open(path, 'rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def save(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def fetch(root):
    root.mkdir(exist_ok=True)
    for name, url in URLS.items():
        path = root / name
        if path.exists():
            if sha(path) != SOURCE_SHA256[name]:
                raise ValueError(f'cached source has unexpected content: {path}')
            continue
        temporary = root / (name + f'.{os.getpid()}.tmp')
        try:
            with urllib.request.urlopen(url, timeout=60) as response, temporary.open('xb') as out:
                shutil.copyfileobj(response, out)
            if sha(temporary) != SOURCE_SHA256[name]:
                raise ValueError(f'public source changed: {url}; review before updating fingerprints')
            os.link(temporary, path)  # Atomic no-clobber publication.
        finally:
            temporary.unlink(missing_ok=True)
    save(root / 'sources.json', URLS)


def prepare(data, root):
    sources = {name: {'url': url, 'sha256': sha(data / name), 'bytes': (data / name).stat().st_size}
               for name, url in URLS.items()}
    assert {k:v['sha256'] for k,v in sources.items()} == SOURCE_SHA256, 'public source content changed; inspect before updating the benchmark'
    root.mkdir()
    with gzip.open(data / 'chr22.fa.gz', 'rb') as src, (root / 'chr22.fa').open('wb') as out:
        shutil.copyfileobj(src, out)
    pysam.faidx(str(root / 'chr22.fa'))
    score_counts = {}
    for score in SCORES:
        with gzip.open(data / f'{score}.txt.gz', 'rt') as f:
            metadata = {}
            for line in f:
                if not line.startswith('#'):
                    header = line.rstrip('\n').split('\t')
                    break
                key, sep, value = line[1:].rstrip('\n').partition('=')
                if sep:
                    metadata[key] = value
            assert metadata['HmPOS_build'] == 'GRCh38', metadata
            rows = []
            total = 0
            for ordinal, row in enumerate(csv.DictReader(f, fieldnames=header, delimiter='\t'), 1):
                total += 1
                if row['hm_chr'] == '22':
                    row['source_ordinal'] = str(ordinal)
                    row['chr_name'] = row['hm_chr']
                    row['chr_position'] = row['hm_pos']
                    rows.append(row)
        # Keep scoring/model annotations; remove harmonization columns after explicit coordinate mapping.
        fields = [k for k in header if not k.startswith('hm_')]
        score_id = f'CHR22_{score}'
        with (root / f'{score_id}.txt').open('w') as out:
            out.write(f"#pgs_id={score_id}\n#genome_build=GRCh38\n#variants_number={len(rows)}\n#trait_reported=Chromosome 22 contribution: {metadata['trait_reported']}\n")
            w = csv.DictWriter(out, fieldnames=fields, delimiter='\t', extrasaction='ignore', lineterminator='\n')
            w.writeheader()
            w.writerows(rows)
        (root / f'{score_id}.source-ordinals.tsv').write_text('ordinal\tsource_ordinal\n' + ''.join(f"{i}\t{row['source_ordinal']}\n" for i, row in enumerate(rows, 1)))
        score_counts[score] = {'original_terms': total, 'chr22_terms': len(rows), 'source_metadata': metadata}
    with pysam.VariantFile(str(data / 'chr22.vcf.gz')) as source:
        names = list(source.header.samples)
        selected = [names[i * (len(names) - 1) // 31] for i in range(32)]
        source.subset_samples(selected)
        with pysam.VariantFile(str(root / 'public32.vcf.gz'), 'wz', header=source.header) as out:
            count = 0
            for record in source:
                assert record.contig == '22'
                out.write(record)
                count += 1
    pysam.tabix_index(str(root / 'public32.vcf.gz'), preset='vcf')
    (root / 'samplesheet.csv').write_text(f'sampleset,path_prefix,chrom,format\npublic32,{root}/public32,22,vcf\n')
    save(root / 'preparation.json', {'schema': 'pgsum-public-chr22-preparation-v1',
         'sources': sources, 'scores': score_counts, 'source_samples': len(names), 'samples': selected,
         'sample_selection': '32 evenly spaced source-header indices, floor(i*(n-1)/31), i=0..31',
         'records': count, 'pysam_version': pysam.__version__,
         'outputs_sha256': {p.name: sha(p) for p in root.iterdir() if p.is_file()},
         'scope': 'All chr22 records retained; all published chr22 terms retained. No complete PGS or population-performance claim.'})
    print(json.dumps({'records': count, 'samples': len(selected), 'terms': {s: c['chr22_terms'] for s, c in score_counts.items()}}, indent=2))


def run_checked(root, args, env=None):
    root.mkdir(exist_ok=True)
    commands_path = root / 'commands.json'
    commands = json.loads(commands_path.read_text()) if commands_path.exists() else []
    args = list(map(str, args))
    commands.append(args)
    save(commands_path, commands)
    with (root / f'command-{len(commands):02d}.log').open('w') as log:
        result = subprocess.run(args, cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'command failed; see {root}/command-{len(commands):02d}.log')


def verify_prepared(root):
    manifest = json.loads((root / 'preparation.json').read_text())
    assert {k:v['sha256'] for k,v in manifest['sources'].items()} == SOURCE_SHA256
    for name, digest in manifest['outputs_sha256'].items():
        assert sha(root / name) == digest, f'changed prepared input: {name}'
    return manifest


def native(prepared, out, pgsum):
    verify_prepared(prepared)
    executable_digest = sha(pgsum)
    out.mkdir()
    run_checked(out, [pgsum, '--threads', '2', 'compile', *[prepared / f'CHR22_{s}.txt' for s in SCORES],
                     '--reference', prepared / 'chr22.fa', '--out', out / 'packs'])
    run_checked(out, [pgsum, '--threads', '2', 'extract', '--gvcf', prepared / 'public32.vcf.gz',
                     '--reference', prepared / 'chr22.fa', '--pack', out / 'packs', '--all-samples',
                     '--accept-missing-quality', '--out', out / 'cohort.pgsc'])
    run_checked(out, [pgsum, '--threads', '2', 'score', '--genotypes', out / 'cohort.pgsc', '--pack', out / 'packs',
                     '--out', out / 'scores'])
    # A representative single-sample export supports direct per-term diagnosis.
    sample = json.loads((prepared / 'preparation.json').read_text())['samples'][0]
    run_checked(out, [pgsum, '--threads', '2', 'run', '--gvcf', prepared / 'public32.vcf.gz',
                     '--reference', prepared / 'chr22.fa', '--pack', out / 'packs', '--sample', sample,
                     '--accept-missing-quality', '--terms', '--out', out / 'single'])
    assert sha(pgsum) == executable_digest, 'native executable changed during scoring'
    save(out / 'identity.json', {'executable_sha256': executable_digest, 'sample': sample, 'version': subprocess.check_output([pgsum, '--version'], text=True).strip()})


def pgsc(prepared, out, nextflow, pipeline, plink):
    verify_prepared(prepared)
    commit = subprocess.check_output(['git', '-C', pipeline, 'rev-parse', 'HEAD'], text=True).strip()
    assert commit == PIPELINE_COMMIT, 'use the pinned pgsc_calc v2.3.0 source commit'
    subprocess.run(['git', '-C', pipeline, 'diff', '--exit-code', 'HEAD'], check=True)
    plink_digest = sha(plink)
    out.mkdir()
    # Local runtime dependencies are explicit; the official pipeline source remains unmodified.
    bindir = out / 'bin'
    bindir.mkdir()
    (bindir / 'plink2').symlink_to(plink)
    env = dict(os.environ)
    env['PATH'] = str(bindir) + os.pathsep + str(Path(sys.executable).parent) + os.pathsep + env['PATH']
    env['NXF_HOME'] = str(out.parent / 'nextflow-home')
    config = out / 'local.config'
    config.write_text("process.executor = 'local'\nprocess.cpus = 2\nprocess.memory = '4 GB'\n"
                      "process.maxForks = 2\nexecutor.queueSize = 2\nconda.enabled = false\n"
                      "docker.enabled = false\nsingularity.enabled = false\n")
    run_checked(out, [nextflow, 'run', pipeline / 'main.nf', '-c', config,
                     '--input', prepared / 'samplesheet.csv', '--scorefile', str(prepared / 'CHR22_*.txt'),
                     '--target_build', 'GRCh38', '--outdir', out / 'results', '--only_score',
                     '--max_cpus', '2', '--max_memory', '4.GB', '--max_time', '2.h',
                     '-ansi-log', 'false'], env)
    aggregates = list((out / 'work').glob('*/*/aggregated_scores.txt.gz'))
    assert len(aggregates) == 1, aggregates
    shutil.copyfile(aggregates[0], out / 'aggregated_scores.txt.gz')
    commit = subprocess.check_output(['git', '-C', pipeline, 'rev-parse', 'HEAD'], text=True).strip()
    assert sha(plink) == plink_digest, 'PLINK executable changed during scoring'
    save(out / 'identity.json', {'pipeline_commit': commit, 'plink_sha256': plink_digest,
         'nextflow_launcher_sha256': sha(nextflow), 'plink_version': subprocess.check_output([plink, '--version'], text=True).strip(),
         'packages': subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'], text=True).splitlines()})


def read_tsv(path):
    if path.suffix == '.gz':
        with gzip.open(path, 'rt') as f:
            return list(csv.DictReader(f, delimiter='\t'))
    if path.suffix == '.zst':
        with path.open('rb') as raw, zstandard.ZstdDecompressor().stream_reader(raw) as decoded:
            return list(csv.DictReader(io.TextIOWrapper(decoded), delimiter='\t'))
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def read_score(path):
    with path.open() as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.rstrip('\n').split('\t')
                return fields, list(csv.DictReader(f, fieldnames=fields, delimiter='\t'))
    raise ValueError('empty score')


def close_enough(a, b):
    return abs(a - b) <= Decimal('0.00001') + Decimal('0.000005') * max(abs(a), abs(b))


def analyze(prepared, native_root, pgsc_root, out, pgsum, plink):
    getcontext().prec = 100
    executable_digests = {'pgsum': sha(pgsum), 'plink2': sha(plink)}
    assert executable_digests['pgsum'] == json.loads((native_root / 'identity.json').read_text())['executable_sha256']
    assert executable_digests['plink2'] == json.loads((pgsc_root / 'identity.json').read_text())['plink_sha256']
    manifest = verify_prepared(prepared)
    out.mkdir()
    samples = manifest['samples']
    single = samples[0]
    aggregates = list((pgsc_root / 'work').glob('*/*/aggregated_scores.txt.gz'))
    assert len(aggregates) == 1
    aggregate_path = aggregates[0]
    external = read_tsv(aggregate_path)
    external_values = {(r['IID'], r['PGS']): Decimal(r['SUM']) for r in external}
    assert len(external_values) == len(external) == len(SCORES) * len(samples)
    native_values = {(r['sample'], r['pgs_id']): Decimal(r['partial_raw_score']) for r in
                     read_tsv(native_root / 'scores/cohort-scores.tsv.zst')}
    assert set(native_values) == set(external_values)
    compare_flags = ['--absolute-tolerance', '0.00001', '--relative-tolerance', '0.000005']
    run_checked(out, [pgsum, 'compare', '--left', native_root / 'scores/cohort-scores.tsv.zst',
                     '--right', aggregate_path, '--right-format', 'pgsc-calc', '--right-sampleset', 'public32',
                     '--out', out / 'default-comparison', *compare_flags])
    match_root = pgsc_root / 'results/public32/match'
    with gzip.open(match_root / 'public32_log.csv.gz', 'rt') as f:
        matches = list(csv.DictReader(f))
    matched = {}
    for r in matches:
        if r['match_status'] == 'matched':
            key = (r['accession'], int(r['row_nr']) + 1)
            assert key not in matched, key
            matched[key] = r
    need_ids = {r['ID'] for r in matched.values()}
    positions = sorted({int(r['chr_position']) for score in SCORES for r in read_score(prepared / f'CHR22_{score}.txt')[1] if r['chr_position'].isdigit()})
    overlaps = collections.defaultdict(list)
    variants = {}
    with pysam.VariantFile(str(prepared / 'public32.vcf.gz')) as source:
        assert list(source.header.samples) == samples
        for r in source:
            key = f'{r.contig}:{r.pos}:{r.ref}:{",".join(r.alts)}'
            lo, hi = bisect.bisect_left(positions, r.pos), bisect.bisect_right(positions, r.stop)
            if lo < hi:
                record = (r.pos, r.ref, r.alts, [r.samples[s]['GT'] for s in samples], tuple(r.filter), tuple(r.format))
                for pos in positions[lo:hi]:
                    overlaps[pos].append(record)
            if key in need_ids:
                assert key not in variants, f'nonunique scored variant {key}'
                assert len(r.alts) == 1, key
                variants[key] = (r.ref, r.alts[0], [r.samples[s]['GT'] for s in samples])
    assert set(variants) == need_ids
    # Verify actual matched scorefile weights against source rows, including reference effect alleles.
    weights = {}
    for path in sorted(match_root.glob('*.scorefile.gz')):
        for row in read_tsv(path):
            for s in SCORES:
                accession = f'CHR22_{s}'
                key = (accession, row['ID'], row['effect_allele'])
                assert key not in weights, key
                weights[key] = Decimal(row[accession])
    required_fields = ['ordinal', 'contig', 'pos', 'model', 'ref', 'alt', 'effect_is_alt', 'weights',
                       'status', 'call_state', 'effect_dosage', 'contribution']
    with (out / 'external-single.tsv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(external[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(r for r in external if r['IID'] == single)
    results = {}
    oracle_rows = []
    for score in SCORES:
        accession = f'CHR22_{score}'
        fields, source_rows = read_score(prepared / f'{accession}.txt')
        native_terms = read_tsv(native_root / 'single' / f'{accession}.terms.tsv')
        assert len(source_rows) == len(native_terms)
        expected = [Decimal(0) for _ in samples]
        matched_weights = {}
        common = []
        categories = collections.Counter()
        delta_categories = collections.defaultdict(Decimal)
        external_terms = []
        common_expected = [Decimal(0) for _ in samples]
        for ordinal, (source_row, term) in enumerate(zip(source_rows, native_terms), 1):
            assert int(term['ordinal']) == ordinal
            source_weight = Decimal(source_row['effect_weight'])
            contribution = None
            right = {'ordinal': ordinal, 'contig': 'chr22', 'pos': source_row['chr_position'],
                     'model': 'additive', 'ref': '', 'alt': '', 'effect_is_alt': '',
                     'weights': str(source_weight), 'status': 'not_matched_by_pgsc_calc',
                     'call_state': '', 'effect_dosage': '', 'contribution': ''}
            match = matched.get((accession, ordinal))
            if match is not None:
                key = (accession, match['ID'], match['matched_effect_allele'])
                assert key not in matched_weights, key
                matched_weights[key] = source_weight
                assert weights[key] == source_weight, (key, weights[key], source_weight)
                assert Decimal(match['effect_weight']) == source_weight
                ref, alt, calls = variants[match['ID']]
                allele = match['matched_effect_allele']
                assert allele in (ref, alt)
                copies = []
                for i, call in enumerate(calls):
                    assert len(call) == 2 and all(x in (0, 1, None) for x in call), call
                    d = None if None in call else sum((ref, alt)[x] == allele for x in call)
                    copies.append(d)
                    if d is not None:
                        expected[i] += source_weight * d
                contribution = None if copies[0] is None else source_weight * copies[0]
                right.update(ref=ref, alt=alt, effect_is_alt=int(allele == alt),
                             status='scorable_observation' if contribution is not None else 'missing_genotype',
                             call_state='missing' if contribution is None else 'observed_reference' if sum(calls[0]) == 0 else 'observed_variant',
                             effect_dosage='' if copies[0] is None else str(copies[0]),
                             contribution='' if contribution is None else str(contribution))
                # Shared-set arithmetic check, explicitly separate from default policy comparison.
                if term['contribution'] != '' and all(d is not None for d in copies) and term['ref'] == ref and term['alt'] == alt and term['effect_is_alt'] == str(int(allele == alt)):
                    assert Decimal(term['contribution']) == contribution, (score, ordinal)
                    common.append((ordinal, source_row, match))
                    for i, d in enumerate(copies):
                        common_expected[i] += source_weight * d
            left_value = None if term['contribution'] == '' else Decimal(term['contribution'])
            if left_value is not None and contribution is not None:
                category = 'both_scored_equal' if left_value == contribution else 'both_scored_different'
            elif left_value is None and contribution is not None:
                category = 'pgsum_excluded:' + term['status']
            elif left_value is not None:
                category = 'pgsc_calc_not_scored'
            else:
                category = 'neither_scored'
            categories[category] += 1
            delta_categories[category] += (contribution or Decimal(0)) - (left_value or Decimal(0))
            external_terms.append(right)
        # Zero coefficient cells added by scorefile pivoting are the only unmatched cells allowed.
        assert all(weight == matched_weights.get(key, Decimal(0)) for key, weight in weights.items() if key[0] == accession)
        for i, sample in enumerate(samples):
            assert close_enough(expected[i], external_values[(sample, accession)]), (sample, score, expected[i], external_values[(sample, accession)])
            oracle_rows.append({'sample': sample, 'score': accession, 'exact_external_contribution_sum': str(expected[i]),
                                'external_reported_sum': str(external_values[(sample, accession)]),
                                'native_reported_sum': str(native_values[(sample, accession)]),
                                'reported_difference': str(external_values[(sample, accession)] - native_values[(sample, accession)])})
        native_sum = sum((Decimal(t['contribution']) for t in native_terms if t['contribution']), Decimal(0))
        assert native_sum == native_values[(single, accession)]
        assert sum(delta_categories.values(), Decimal(0)) == expected[0] - native_sum
        right_terms = out / f'{accession}.external-terms.tsv'
        with right_terms.open('w') as f:
            w = csv.DictWriter(f, fieldnames=required_fields, delimiter='\t', lineterminator='\n')
            w.writeheader()
            w.writerows(external_terms)
        run_checked(out, [pgsum, 'compare', '--left', native_root / 'single/scores.tsv', '--left-sample', single,
                         '--right', out / 'external-single.tsv', '--right-format', 'pgsc-calc',
                         '--left-terms', native_root / 'single' / f'{accession}.terms.tsv', '--right-terms', right_terms,
                         '--term-sample', single, '--term-score', accession, '--out', out / f'terms-{score}', *compare_flags])
        term_report = json.loads((out / f'terms-{score}/comparison.json').read_text())['terms']
        assert term_report['left']['score_reconciliation'] == 'exact'
        assert term_report['right']['score_reconciliation'] in ('exact', 'within_tolerance')
        # Build a labelled common-term model preserving original weights/annotations and ordinal mapping.
        common_id = 'COMMON_' + score
        common_file = out / f'{common_id}.txt'
        with common_file.open('w') as f:
            f.write(f'#pgs_id={common_id}\n#genome_build=GRCh38\n#variants_number={len(common)}\n')
            w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
            w.writeheader()
            w.writerows(row for _, row, _ in common)
        with (out / f'{common_id}.weights').open('w') as f:
            f.write(f'ID\tALLELE\t{common_id}\n')
            for _, row, match in common:
                f.write(f"{match['ID']}\t{match['matched_effect_allele']}\t{row['effect_weight']}\n")
        (out / f'{common_id}.ordinals.txt').write_text(''.join(str(ordinal) + '\n' for ordinal, _, _ in common))
        run_checked(out, [pgsum, 'compile', common_file, '--reference', prepared / 'chr22.fa', '--out', out / 'common-packs'])
        run_checked(out, [pgsum, '--threads', '2', 'extract', '--gvcf', prepared / 'public32.vcf.gz',
                         '--reference', prepared / 'chr22.fa', '--pack', out / 'common-packs' / f'{common_id}.pgsp',
                         '--all-samples', '--accept-missing-quality', '--out', out / f'{common_id}.pgsc'])
        run_checked(out, [pgsum, '--threads', '2', 'score', '--genotypes', out / f'{common_id}.pgsc',
                         '--pack', out / 'common-packs' / f'{common_id}.pgsp', '--out', out / common_id])
        run_checked(out, [plink, '--threads', '2', '--vcf', prepared / 'public32.vcf.gz', '--set-all-var-ids', '@:#:$r:$a',
                         '--new-id-max-allele-len', '1000', '--max-alleles', '2', '--extract', out / f'{common_id}.weights',
                         '--score', out / f'{common_id}.weights', '1', '2', '3', 'header-read', 'no-mean-imputation',
                         'cols=+scoresums', '--out', out / f'plink-{score}'])
        run_checked(out, [pgsum, 'compare', '--left', out / common_id / 'cohort-scores.tsv.zst',
                         '--right', out / f'plink-{score}.sscore', '--right-format', 'plink2',
                         '--right-score-column', f'{common_id}_SUM', '--right-score-id', common_id,
                         '--out', out / f'common-{score}', '--fail-on-difference', *compare_flags])
        comparison_rows = read_tsv(out / f'common-{score}/scores.tsv')
        assert len(comparison_rows) == len(samples)
        for row in comparison_rows:
            i = samples.index(row['sample'])
            assert Decimal(row['left']) == common_expected[i], (score, row, common_expected[i])
            assert close_enough(Decimal(row['right']), common_expected[i])
        results[score] = {'source_terms': len(source_rows), 'pgsc_matched_terms': len(matched_weights),
                         'single_sample': single, 'categories': dict(categories),
                         'exact_delta_by_category': {k: str(v) for k, v in delta_categories.items()},
                         'native_single_sum': str(native_sum), 'external_exact_single_sum': str(expected[0]),
                         'common_terms': len(common), 'common_comparison': json.loads((out / f'common-{score}/comparison.json').read_text())['counts']}
    full_oracle = verify_default_samples(prepared, native_root, out, pgsum, samples, matched, variants, overlaps, native_values, external_values)
    save(out / 'oracle-scores.json', oracle_rows)
    assert executable_digests == {'pgsum': sha(pgsum), 'plink2': sha(plink)}, 'executable changed during validation'
    save(out / 'results.json', {'schema': 'pgsum-public-chr22-validation-v1', 'preparation': manifest,
         'native_identity': json.loads((native_root / 'identity.json').read_text()),
         'pgsc_identity': json.loads((pgsc_root / 'identity.json').read_text()), 'scores': results,
         'default_policy_oracle': full_oracle,
         'default_comparison': json.loads((out / 'default-comparison/comparison.json').read_text())['counts'],
         'oracle': 'Python Decimal over original retained VCF GT and verified source weights; checks all external sums and all common-set native/PLINK sums; default call contributions independently checked and policy differences reconciled per term for all 32 samples.',
         'scope': 'Chromosome-22 partial contributions only. 32 selected public samples. Default policies differ; common-set agreement is conditional on identical retained terms. No full-genome, ancestry, DS/GP or gVCF validation.',
         'runner_sha256': sha(Path(__file__)),
         'outputs_sha256': {p.name: sha(p) for p in out.iterdir() if p.is_file()}})
    print(json.dumps(results, indent=2))



def verify_default_samples(prepared, native_root, out, pgsum, samples, matched, variants, overlaps, native_values, external_values):
    """Check each retained native contribution against source weights, FASTA and original GT.

    Eligibility exclusions remain tool policy decisions, not assertions of biological truth.
    """
    complement = str.maketrans('ACGT', 'TGCA')
    source = {score: read_score(prepared / f'CHR22_{score}.txt')[1] for score in SCORES}
    provenance = {score: {int(r['ordinal']): int(r['source_ordinal']) for r in read_tsv(prepared / f'CHR22_{score}.source-ordinals.tsv')} for score in SCORES}
    results = []
    files = {}
    passing = 0
    (out / 'default-samples').mkdir()
    diff_path = out / 'default-term-differences.tsv'
    with pysam.FastaFile(str(prepared / 'chr22.fa')) as fasta, diff_path.open('w') as differences:
        w = csv.writer(differences, delimiter='\t', lineterminator='\n')
        w.writerow(['sample', 'score', 'ordinal', 'published_source_ordinal', 'position', 'category', 'native_contribution', 'external_contribution', 'exact_delta', 'source_records_at_target_json'])
        for i, sample in enumerate(samples):
            directory = native_root / 'single' if i == 0 else out / 'default-samples' / sample
            if i:
                run_checked(out, [pgsum, '--threads', '2', 'run', '--gvcf', prepared / 'public32.vcf.gz',
                                 '--reference', prepared / 'chr22.fa', '--pack', native_root / 'packs',
                                 '--sample', sample, '--accept-missing-quality', '--terms', '--out', directory])
            for score in SCORES:
                accession = f'CHR22_{score}'
                path = directory / f'{accession}.terms.tsv'
                terms = read_tsv(path)
                files[f'{sample}/{accession}.terms.tsv'] = sha(path)
                assert len(terms) == len(source[score])
                ns, es = Decimal(0), Decimal(0)
                categories = collections.Counter()
                deltas = collections.defaultdict(Decimal)
                for ordinal, (row, term) in enumerate(zip(source[score], terms), 1):
                    weight = Decimal(row['effect_weight'])
                    assert term['model'] == 'additive'
                    left = None
                    if term['contribution']:
                        pos = int(row['chr_position'])
                        base = fasta.fetch('chr22', pos-1, pos).upper()
                        ea, oa = row['effect_allele'], row['other_allele']
                        assert len(ea) == len(oa) == 1 and ea in 'ACGT' and oa in 'ACGT' and ea != oa
                        pairs = [pair for pair in [(ea, oa), (ea.translate(complement), oa.translate(complement))] if base in pair]
                        assert len(pairs) == 1, (score, ordinal, pairs)
                        forward_effect, forward_other = pairs[0]
                        alternate = forward_other if forward_effect == base else forward_effect
                        assert (term['ref'], term['alt'], term['effect_is_alt']) == (base, alternate, str(int(forward_effect == alternate)))
                        assert len(overlaps[pos]) == 1, (score, ordinal, term, overlaps[pos])
                        start, ref, alts, calls, filters, formats = overlaps[pos][0]
                        assert start == pos and ref == base and filters == ('PASS',) and formats == ('GT',)
                        call = calls[i]
                        assert len(call) == 2 and None not in call
                        alleles = (ref,) + alts
                        actual_alleles = [alleles[x] for x in call]
                        assert all(a in (base, alternate) for a in actual_alleles), (sample, score, ordinal)
                        copies = actual_alleles.count(forward_effect)
                        left = weight * copies
                        assert Decimal(term['weights']) == weight
                        assert Decimal(term['effect_dosage']) == copies
                        assert Decimal(term['contribution']) == left
                        ns += left
                        passing += 1
                    right = None
                    match = matched.get((accession, ordinal))
                    if match is not None:
                        ref, alt, calls = variants[match['ID']]
                        call = calls[i]
                        if None not in call:
                            right = weight * sum((ref, alt)[x] == match['matched_effect_allele'] for x in call)
                            es += right
                    if left is not None and right is not None:
                        category = 'both_scored_equal' if left == right else 'both_scored_different'
                    elif right is not None:
                        category = 'pgsum_excluded:' + term['status']
                    elif left is not None:
                        category = 'pgsc_calc_not_scored'
                    else:
                        category = 'neither_scored'
                    delta = (right or Decimal(0)) - (left or Decimal(0))
                    categories[category] += 1
                    deltas[category] += delta
                    if category not in ('both_scored_equal', 'neither_scored'):
                        w.writerow([sample, accession, ordinal, provenance[score][ordinal], row['chr_position'], category,
                                    '' if left is None else str(left), '' if right is None else str(right), str(delta),
                                    json.dumps([{'position': rec[0], 'ref': rec[1], 'alt': rec[2], 'gt': rec[3][i]} for rec in overlaps.get(int(row['chr_position']), [])], separators=(',', ':'))])
                assert ns == native_values[(sample, accession)]
                assert close_enough(es, external_values[(sample, accession)])
                assert sum(deltas.values(), Decimal(0)) == es - ns
                assert categories['both_scored_different'] == 0, (sample, score, categories)
                results.append({'sample': sample, 'score': accession, 'native_exact': str(ns), 'external_exact': str(es),
                                'categories': dict(categories), 'delta_by_category': {k: str(v) for k,v in deltas.items()}})
    save(out / 'default-policy-oracle.json', results)
    save(out / 'default-term-file-hashes.json', files)
    return {'sample_score_pairs': len(results), 'native_contributions_checked': passing,
            'both_scored_different': 0, 'results_sha256': sha(out / 'default-policy-oracle.json'),
            'differences_sha256': sha(diff_path), 'term_exports_sha256': sha(out / 'default-term-file-hashes.json')}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=['fetch', 'prepare', 'native', 'pgsc', 'analyze'])
    p.add_argument('--data', type=Path)
    p.add_argument('--prepared', type=Path)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--pgsum', type=Path)
    p.add_argument('--plink2', type=Path)
    p.add_argument('--nextflow', type=Path)
    p.add_argument('--pipeline', type=Path)
    p.add_argument('--native', type=Path)
    p.add_argument('--pgsc', type=Path)
    a = p.parse_args()
    required = {'fetch': [], 'prepare': ['data'], 'native': ['prepared','pgsum'], 'pgsc': ['prepared','nextflow','pipeline','plink2'], 'analyze': ['prepared','native','pgsc','pgsum','plink2']}
    for key in required[a.stage]:
        if getattr(a, key) is None:
            p.error(f'--{key} is required for {a.stage}')
    if not __debug__:
        p.error('Do not use python -O: validation assertions must be enabled')
    for name in ['data', 'prepared', 'out', 'pgsum', 'plink2', 'nextflow', 'pipeline', 'native', 'pgsc']:
        if getattr(a, name) is not None:
            setattr(a, name, getattr(a, name).resolve())
    if a.stage == 'fetch': fetch(a.out)
    elif a.stage == 'prepare': prepare(a.data, a.out)
    elif a.stage == 'native': native(a.prepared, a.out, a.pgsum)
    elif a.stage == 'pgsc': pgsc(a.prepared, a.out, a.nextflow, a.pipeline, a.plink2)
    elif a.stage == 'analyze': analyze(a.prepared, a.native, a.pgsc, a.out, a.pgsum, a.plink2)


if __name__ == '__main__': main()
