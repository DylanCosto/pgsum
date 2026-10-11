#!/usr/bin/env python3
"""Public CAD-score missingness experiment with independent source-GT arithmetic.

All published rows of four sparse models are retained. Robustness is measured against
an explicitly conditional common-callable baseline, never as disease prediction.
Requires pysam (the pinned public_chr22.requirements.txt environment also works).
"""
import argparse
import collections
import concurrent.futures
import csv
from decimal import Decimal, getcontext
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import urllib.request

import pysam

SCORES = {
    'PGS000010': ('2d8d217d407d50805681323ce4264a1b905bc17fb1b1a7b4e5d873cf7fe08549', 27),
    'PGS000011': ('c258514c3c908b8f4bf6e6d882d3f7ffb3d9b7bba529994e39072497f25900ea', 50),
    'PGS000057': ('fb61fbbec7381bce2be871de1c5d87ccd8c62289569279bb7d759e1f46357744', 57),
    'PGS000349': ('03451f003bd05ddb82532b28aac35b8f31fd095db1ce09e33a76eb23d5bbfefe', 70),
}
METADATA_SHA256 = {
    'PGS000010':'cfb5ccb9fcb61d5e03cc037c3eb8ab5a0ff4cfff462942f08eb9fa1a38060987',
    'PGS000011':'bea497c43f3acb137286d0caef9a1fe2105fd4a18c7e7262e778df9667ddc073',
    'PGS000057':'38ac150b6c7e7fde109d30d4b0dc35bc3e778f7f7eaca629428c2b59a577f846',
    'PGS000349':'4a86cbdde67246a60d7e4c062427f1c466b6d4658018151876e3bc2f301ebedb',
}
VCF_TEMPLATE = ('https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/'
                'ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz')
SCHEMA = 'pgsum-public-robustness-preparation-v1'
PINNED_SELECTION_SHA256 = 'e40d931358126755cb626e51dda49bd8974fe57574326366d263fc2e604f9a1c'
REFERENCE_URL = ('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/'
                 'GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/'
                 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz')


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def genotype_fingerprint(path):
    """Hash sample order, coordinate/allele/filter records and GT arrays independent of gzip."""
    digest = hashlib.sha256()
    with pysam.VariantFile(str(path)) as source:
        samples = list(source.header.samples)
        digest.update((json.dumps(samples, separators=(',', ':')) + '\n').encode())
        for record in source:
            row = [record.contig, record.pos, record.ref, list(record.alts), list(record.filter),
                   [list(record.samples[sample]['GT']) for sample in samples]]
            digest.update((json.dumps(row, separators=(',', ':')) + '\n').encode())
    return digest.hexdigest()


def save(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def read_score(path):
    metadata = {}
    with gzip.open(path, 'rt') as source:
        for line in source:
            if not line.startswith('#'):
                rows = list(csv.DictReader(source, fieldnames=line.rstrip('\n').split('\t'), delimiter='\t'))
                return metadata, rows
            key, separator, value = line[1:].rstrip('\n').partition('=')
            if separator:
                metadata[key] = value
    raise ValueError(f'empty score: {path}')


def download(url, path, expected=None):
    if path.exists():
        if expected is not None and sha(path) != expected:
            raise ValueError(f'cached content changed: {path}')
        return
    temporary = path.with_name(path.name + '.part')
    try:
        with urllib.request.urlopen(url, timeout=60) as source, temporary.open('xb') as target:
            shutil.copyfileobj(source, target)
        if expected is not None and sha(temporary) != expected:
            raise ValueError(f'upstream content changed: {url}; inspect before changing fingerprints')
        temporary.rename(path)
    finally:
        temporary.unlink(missing_ok=True)


def fetch(out):
    out.mkdir(parents=True, exist_ok=True)
    for score, (digest, _) in SCORES.items():
        url = f'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{score}/ScoringFiles/Harmonized/{score}_hmPOS_GRCh38.txt.gz'
        download(url, out / f'{score}.txt.gz', digest)
        download(f'https://www.pgscatalog.org/rest/score/{score}', out / f'{score}.metadata.json', METADATA_SHA256[score])
    print(f'Verified {len(SCORES)} published scoring files', flush=True)


def intervals(positions, padding=100):
    result = []
    for pos in sorted(set(positions)):
        start, end = max(0, pos - 1 - padding), pos + padding
        if result and start <= result[-1][1]:
            result[-1][1] = max(result[-1][1], end)
        else:
            result.append([start, end])
    return result


def prepare(data, out, reference):
    out.mkdir()
    sources, by_chrom, model_info = {}, collections.defaultdict(set), {}
    for score, (digest, count) in SCORES.items():
        source = data / f'{score}.txt.gz'
        assert sha(source) == digest
        metadata, rows = read_score(source)
        assert metadata['HmPOS_build'] == 'GRCh38'
        assert len(rows) == count == int(metadata['variants_number'])
        assert metadata['trait_efo'] == 'MONDO_0005010'
        assert all(row['hm_chr'].isdigit() and 1 <= int(row['hm_chr']) <= 22 and row['hm_pos'].isdigit() for row in rows)
        for row in rows:
            by_chrom[row['hm_chr']].add(int(row['hm_pos']))
        shutil.copyfile(source, out / source.name)
        metadata_path = data / f'{score}.metadata.json'
        assert sha(metadata_path) == METADATA_SHA256[score]
        shutil.copyfile(metadata_path, out / metadata_path.name)
        sources[metadata_path.name] = {'sha256':METADATA_SHA256[score], 'bytes':metadata_path.stat().st_size, 'url':f'https://www.pgscatalog.org/rest/score/{score}'}
        sources[source.name] = {'sha256': digest, 'bytes': source.stat().st_size,
            'url': f'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{score}/ScoringFiles/Harmonized/{score}_hmPOS_GRCh38.txt.gz'}
        model_info[score] = {'published_rows_retained': count, 'metadata': metadata}
    query = {chrom: intervals(positions) for chrom, positions in by_chrom.items()}
    # One header determines selection; subsequent chromosome headers must have identical samples.
    first_chrom = min(by_chrom, key=int)
    with pysam.VariantFile(VCF_TEMPLATE.format(chrom=first_chrom)) as source:
        all_samples = list(source.header.samples)
    selected = [all_samples[i * (len(all_samples) - 1) // 127] for i in range(128)]
    targets, refs = selected[::4], [sample for i, sample in enumerate(selected) if i % 4]
    assert len(set(targets) & set(refs)) == 0 and len(targets) == 32 and len(refs) == 96
    chromosome_dir = out / 'source-regions'
    chromosome_dir.mkdir()
    def chromosome(chrom):
        url = VCF_TEMPLATE.format(chrom=chrom)
        index_path = chromosome_dir / f'chr{chrom}.source.tbi'
        download(url + '.tbi', index_path)
        path = chromosome_dir / f'chr{chrom}.vcf.gz'
        with pysam.VariantFile(url, index_filename=str(index_path)) as source:
            assert list(source.header.samples) == all_samples
            source.subset_samples(selected)
            with pysam.VariantFile(str(path), 'wz', header=source.header) as dest:
                previous, records = None, 0
                for start, end in query[chrom]:
                    for record in source.fetch(chrom, start, end):
                        identity = (record.pos, record.ref, record.alts)
                        if identity != previous:
                            dest.write(record)
                            records += 1
                        previous = identity
        pysam.tabix_index(str(path), preset='vcf')
        print(f'chr{chrom}: retained {records} source records across {len(query[chrom])} intervals', flush=True)
        return chrom, {'url': url, 'query_intervals_0based_halfopen': query[chrom], 'records': records,
                       'subset_sha256': sha(path), 'source_index_sha256': sha(index_path),
                       'source_samples': len(all_samples)}
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        genotype_sources = dict(pool.map(chromosome, sorted(by_chrom, key=int)))
    # A fresh header avoids carrying source population AFs into reference-frequency estimation.
    header = pysam.VariantHeader()
    header.add_line('##fileformat=VCFv4.2')
    with pysam.FastaFile(str(reference)) as fasta:
        for chrom in sorted(by_chrom, key=int):
            name = 'chr' + chrom if 'chr' + chrom in fasta.references else chrom
            header.contigs.add(chrom, length=fasta.get_reference_length(name))
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Public source GT">')
    for sample in selected:
        header.add_sample(sample)
    combined = out / 'public128.vcf.gz'
    with pysam.VariantFile(str(combined), 'wz', header=header) as dest:
        for chrom in sorted(by_chrom, key=int):
            with pysam.VariantFile(str(chromosome_dir / f'chr{chrom}.vcf.gz')) as source:
                for record in source:
                    new = dest.new_record(contig=chrom, start=record.start, stop=record.stop,
                                          alleles=record.alleles, filter=tuple(record.filter))
                    for sample in selected:
                        new.samples[sample]['GT'] = record.samples[sample]['GT']
                        new.samples[sample].phased = record.samples[sample].phased
                    dest.write(new)
    pysam.tabix_index(str(combined), preset='vcf')
    assert genotype_fingerprint(combined) == PINNED_SELECTION_SHA256, 'public genotype selection changed; inspect before updating fingerprint'
    for label, samples in [('targets', targets), ('reference', refs)]:
        with pysam.VariantFile(str(combined)) as source:
            source.subset_samples(samples)
            with pysam.VariantFile(str(out / f'{label}.vcf.gz'), 'wz', header=source.header) as dest:
                for record in source:
                    dest.write(record)
        pysam.tabix_index(str(out / f'{label}.vcf.gz'), preset='vcf')
    # This literal INFO decimal, not a rounded float from an in-memory VCF API, is the fill oracle input.
    frequency_rows, common_sites = [], set()
    with pysam.VariantFile(str(combined)) as source:
        for record in source:
            if len(record.alts) != 1:
                continue
            calls = [allele for sample in refs for allele in record.samples[sample]['GT'] if allele is not None]
            if not calls:
                continue
            assert all(allele in (0, 1) for allele in calls)
            alt_count = sum(calls)
            af = Decimal(alt_count) / Decimal(len(calls))
            frequency_rows.append((record.contig, record.pos, record.ref, record.alts[0], str(af), alt_count, len(calls)))
            if min(af, 1 - af) >= Decimal('0.05'):
                common_sites.add((record.contig, record.pos))
    with gzip.open(out / 'reference-frequencies.vcf.gz', 'wt') as dest:
        dest.write('##fileformat=VCFv4.2\n##INFO=<ID=REF_AF,Number=A,Type=Float,Description="ALT AF in 96 disjoint reference samples">\n')
        dest.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for chrom, pos, ref, alt, af, _, _ in frequency_rows:
            dest.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tREF_AF={af}\n')
    with (out / 'reference-frequencies.tsv').open('w') as dest:
        writer = csv.writer(dest, delimiter='\t', lineterminator='\n')
        writer.writerow(['chrom', 'pos', 'ref', 'alt', 'af', 'alt_count', 'allele_count'])
        writer.writerows(frequency_rows)
    (out / 'reference-common.bed').write_text(''.join(f'{chrom}\t{pos - 1}\t{pos}\n' for chrom, pos in sorted(common_sites, key=lambda x:(int(x[0]),x[1]))))
    (out / 'regional-dropout.bed').write_text('9\t21500000\t22500000\n')
    reference = reference.resolve()
    manifest = {'schema': SCHEMA, 'scores': model_info, 'sources': sources,
                'genotype_sources': genotype_sources, 'source_samples': len(all_samples),
                'selection_sha256': PINNED_SELECTION_SHA256,
                'sample_selection': '128 evenly spaced original VCF header indices floor(i*(n-1)/127); every fourth selected sample is a target, other 96 are reference',
                'targets': targets, 'reference_samples': refs, 'reference_panels': {'panel_a': refs[::2], 'panel_b': refs[1::2]},
                'reference': {'path': str(reference), 'url': REFERENCE_URL, 'sha256': sha(reference), 'fai_sha256': sha(str(reference) + '.fai')},
                'pysam_version': pysam.__version__,
                'scope': 'All four model inventories retained. Evaluation is conditional on common-callable terms; no phenotype or predictive accuracy assessment.',
                'coverage_mask': 'Keep sites with MAF >= 0.05 in 96 disjoint reference samples. Approximation of common-variant coverage, not a commercial array manifest.',
                'regional_mask': 'Controlled contiguous dropout of GRCh38 chr9:21,500,000-22,500,000 (zero-based half-open), not an empirical sequencing failure map.',
                'outputs_sha256': {str(p.relative_to(out)): sha(p) for p in out.rglob('*') if p.is_file()}}
    save(out / 'preparation.json', manifest)
    print(json.dumps({'models': {s:n for s,(_,n) in SCORES.items()}, 'targets':len(targets), 'reference':len(refs),
                      'source_records':sum(x['records'] for x in genotype_sources.values())}), flush=True)


def verify_prepared(prepared):
    manifest = json.loads((prepared / 'preparation.json').read_text())
    assert manifest['schema'] == SCHEMA
    assert genotype_fingerprint(prepared / 'public128.vcf.gz') == PINNED_SELECTION_SHA256
    for name, digest in manifest['outputs_sha256'].items():
        assert sha(prepared / name) == digest, f'changed prepared input: {name}'
    for score, (digest, count) in SCORES.items():
        assert sha(prepared / f'{score}.txt.gz') == digest
        assert sha(prepared / f'{score}.metadata.json') == METADATA_SHA256[score]
        assert manifest['scores'][score]['published_rows_retained'] == count
    reference = Path(manifest['reference']['path'])
    assert sha(reference) == manifest['reference']['sha256']
    assert sha(str(reference) + '.fai') == manifest['reference']['fai_sha256']
    return manifest


def run_checked(out, command):
    path = out / 'commands.json'
    commands = json.loads(path.read_text()) if path.exists() else []
    command = list(map(str, command))
    commands.append(command)
    save(path, commands)
    log_path = out / f'command-{len(commands):02d}.log'
    with log_path.open('w') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'command failed: {log_path}')


def native(prepared, out, pgsum):
    manifest = verify_prepared(prepared)
    out.mkdir()
    identity = {'pgsum_sha256': sha(pgsum), 'version': subprocess.check_output([pgsum, '--version'], text=True).strip()}
    run_checked(out, [pgsum, 'compile', *[prepared / f'{score}.txt.gz' for score in SCORES],
                      '--reference', manifest['reference']['path'], '--out', out / 'packs'])
    for label in ['targets', 'reference']:
        run_checked(out, [pgsum, 'extract', '--gvcf', prepared / f'{label}.vcf.gz', '--reference', manifest['reference']['path'],
                         '--pack', out / 'packs', '--all-samples', '--accept-missing-quality', '--out', out / f'{label}.pgsc'])
    run_checked(out, [pgsum, 'frequencies', '--vcf', prepared / 'reference-frequencies.vcf.gz', '--field', 'REF_AF', '--out', out / 'reference.pgsf'])
    scenarios = []
    for rate in [0.05, 0.2]:
        for seed in [42, 271828]:
            scenarios.append({'id': f'random_{int(rate*100):02d}_s{seed}', 'masks':[{'kind':'random','rate':rate,'seed':seed,'unit':'sample_site'}]})
    scenarios.extend([
        {'id':'regional_dropout', 'masks':[{'kind':'bed','path':str(prepared / 'regional-dropout.bed'),'mode':'drop'}]},
        {'id':'reference_common_coverage', 'masks':[{'kind':'bed','path':str(prepared / 'reference-common.bed'),'mode':'keep'}]},
    ])
    scenarios = [variant for scenario in scenarios for variant in [scenario,dict(scenario,id=scenario['id']+'_fill',fill_frequencies='reference.pgsf', frequency_source_note='ALT frequencies frozen from the 96 named reference samples in preparation.json; all 32 evaluation sample IDs are disjoint; no evaluation GT contributes.')]]
    experiment = {'schema':'pgsum-robustness-experiment-v1', 'genotypes':'targets.pgsc', 'packs':[str(path.relative_to(out)) for path in sorted((out / 'packs').glob('*.pgsp'))],
                  'scoring':{'allow_inferred_other_allele':True}, 'cutoffs':[90,95], 'scenarios':scenarios,
                  'reference_sets':[{'id':name,'genotypes':'reference.pgsc','samples':samples} for name,samples in manifest['reference_panels'].items()],
                  'regions':[{'id':'chr9_dropout_region','path':str(prepared / 'regional-dropout.bed')}], 'top_terms':10}
    save(out / 'experiment.json', experiment)
    run_checked(out, [pgsum, 'robustness', '--experiment', out / 'experiment.json', '--out', out / 'robustness', '--terms'])
    assert sha(pgsum) == identity['pgsum_sha256'], 'executable changed during experiment'
    save(out / 'identity.json', identity)
    print(f'Experiment complete: {out / "robustness"}', flush=True)


def read_tsv(path):
    if path.suffix == '.zst':
        import zstandard
        with path.open('rb') as raw, zstandard.ZstdDecompressor().stream_reader(raw) as decoded:
            return list(csv.DictReader(io.TextIOWrapper(decoded), delimiter='\t'))
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def independent_terms(prepared, manifest):
    """Orient literal SNVs from FASTA and source allele text; derive GT dosage independently.

    This intentionally supports only the four pinned simple additive inventories. It
    does not read engine packs, extracted tables, or engine contributions as an oracle.
    """
    from fractions import Fraction
    selected = manifest['targets'] + manifest['reference_samples']
    positions = collections.defaultdict(set)
    model_rows = {}
    for score in SCORES:
        _, model_rows[score] = read_score(prepared / f'{score}.txt.gz')
        for row in model_rows[score]:
            positions[row['hm_chr']].add(int(row['hm_pos']))
    overlaps = collections.defaultdict(list)
    with pysam.VariantFile(str(prepared / 'public128.vcf.gz')) as source:
        for record in source:
            for position in positions[record.contig]:
                if record.pos <= position <= record.stop:
                    overlaps[(record.contig, position)].append(record.copy())
    frequencies = {(row['chrom'], int(row['pos']), row['ref'], row['alt']): Fraction(row['af'])
                   for row in read_tsv(prepared / 'reference-frequencies.tsv')}
    # Audit disjointness and recompute the frequency counts directly from source calls.
    assert set(manifest['targets']).isdisjoint(manifest['reference_samples'])
    for row in read_tsv(prepared / 'reference-frequencies.tsv'):
        key = (row['chrom'], int(row['pos']))
        records = [record for record in overlaps[key] if record.pos == key[1] and record.ref == row['ref'] and record.alts == (row['alt'],)]
        if not records:  # Context records outside model anchors do not affect this experiment.
            continue
        assert len(records) == 1
        calls = [a for s in manifest['reference_samples'] for a in records[0].samples[s]['GT'] if a is not None]
        assert (sum(calls), len(calls)) == (int(row['alt_count']), int(row['allele_count']))
        assert Decimal(row['af']) == Decimal(sum(calls)) / Decimal(len(calls))
    complement = str.maketrans('ACGT', 'TGCA')
    result = {}
    with pysam.FastaFile(manifest['reference']['path']) as fasta:
        for score, rows in model_rows.items():
            included, excluded = [], collections.Counter()
            for ordinal, row in enumerate(rows, 1):
                chrom, pos = row['hm_chr'], int(row['hm_pos'])
                name = 'chr' + chrom if 'chr' + chrom in fasta.references else chrom
                ref = fasta.fetch(name, pos - 1, pos).upper()
                author_other = row.get('other_allele', '')
                other = author_other or row.get('hm_inferOtherAllele', '')
                effect = row['effect_allele']
                assert all(not row.get(key) for key in ['is_haplotype','is_diplotype','is_interaction','is_dominant','is_recessive','variant_description'])
                if len(effect) != 1 or len(other) != 1 or effect not in 'ACGT' or other not in 'ACGT' or effect == other:
                    excluded['model_or_orientation_unscorable'] += 1
                    continue
                orientations = [(e,o) for e,o in [(effect,other),(effect.translate(complement),other.translate(complement))] if ref in (e,o)]
                if len(orientations) != 1:
                    excluded['model_or_orientation_unscorable'] += 1
                    continue
                effect, other = orientations[0]
                records = overlaps[(chrom,pos)]
                if len(records) != 1:
                    excluded['baseline_target_not_common_callable'] += 1
                    continue
                record = records[0]
                assert record.ref == ref and record.pos == pos and len(record.alts) == 1, (score,ordinal)
                calls = {sample:record.samples[sample]['GT'] for sample in selected}
                def callable_sample(sample):
                    call = calls[sample]
                    return (tuple(record.filter) in [('PASS',), ()] and len(call) == 2 and all(x in (0,1) for x in call)
                            and all(record.alleles[x] in (effect,other) for x in call))
                if not all(callable_sample(sample) for sample in manifest['targets']):
                    excluded['baseline_target_not_common_callable'] += 1
                    continue
                if not all(callable_sample(sample) for sample in manifest['reference_samples']):
                    excluded['reference_not_common_callable'] += 1
                    continue
                weight = Fraction(row['effect_weight'])
                values = {sample:weight * sum(record.alleles[x] == effect for x in calls[sample]) for sample in selected}
                af = frequencies[(chrom,pos,record.ref,record.alts[0])]
                # Production fill requires authored alleles; opt-in genotype orientation does not
                # remove the missing-author-allele review reason for hypothetical frequency filling.
                fill = 2 * weight * (af if effect == record.alts[0] else 1-af) if author_other else None
                included.append({'ordinal':ordinal,'chrom':chrom,'pos':pos,'values':values,'fill':fill,
                                 'source_weight':row['effect_weight'],'effect':effect,'other':other})
            result[score] = {'terms':included,'excluded':dict(excluded)}
    return result


def hidden(mask, term, sample, beds):
    if mask['kind'] == 'random':
        key = (b'pgsum-robustness-mask-sha256-v1\0' + int(mask['seed']).to_bytes(8,'little')
               + int(term['chrom']).to_bytes(1,'little') + term['pos'].to_bytes(4,'little'))
        per_sample = mask['unit'] == 'sample_site'
        key += bytes([per_sample])
        if per_sample:
            encoded = sample.encode()
            key += len(encoded).to_bytes(8,'little') + encoded
        draw = (int.from_bytes(hashlib.sha256(key).digest()[:8],'little') >> 11) / (1 << 53)
        return draw < mask['rate']
    present = any(chrom.removeprefix('chr') == term['chrom'] and start <= term['pos']-1 < end for chrom,start,end in beds[mask['path']])
    return present if mask['mode'] == 'drop' else not present


def percentile(value, reference):
    if len(reference) < 2 or len(set(reference)) == 1:
        return None
    return 100 * (sum(x < value for x in reference) + 0.5 * sum(x == value for x in reference)) / len(reference)


def verify(prepared, native_root, out, pgsum):
    from fractions import Fraction
    import math
    import statistics
    manifest = verify_prepared(prepared)
    identity = json.loads((native_root / 'identity.json').read_text())
    assert sha(pgsum) == identity['pgsum_sha256'], 'use the same executable as the native experiment'
    report_root = native_root / 'robustness'
    report = json.loads((report_root / 'report.json').read_text())
    for filename,digest in report['outputs'].items():
        assert sha(report_root / filename) == digest.removeprefix('sha256:')
    for filename,digest in report['input_sha256'].items():
        assert sha(filename) == digest.removeprefix('sha256:')
    experiment = report['experiment']
    terms = independent_terms(prepared, manifest)
    scenarios = {scenario['id']:scenario for scenario in experiment['scenarios']}
    beds = {}
    for scenario in scenarios.values():
        for mask in scenario['masks']:
            if mask['kind'] == 'bed' and mask['path'] not in beds:
                beds[mask['path']] = [(fields[0],int(fields[1]),int(fields[2])) for fields in
                    (line.split() for line in Path(mask['path']).read_text().splitlines()) if fields and not fields[0].startswith('#')]
    actual = {(r['pgs_id'],r['scenario'],r['sample']):r for r in read_tsv(report_root / 'samples.tsv')}
    assert len(actual) == len(SCORES)*len(scenarios)*len(manifest['targets'])
    expected_changes, expected_regions, baseline, expected_scenario = {}, {}, {}, {}
    mask_checks = 0
    counts = {}
    near = lambda x,y: math.isclose(float(x),float(y),rel_tol=1e-12,abs_tol=1e-12)
    for score,item in terms.items():
        engine = next(row for row in report['scores'] if row['pgs_id'] == score)
        assert engine['total_model_terms'] == SCORES[score][1]
        assert engine['baseline_terms'] == len(item['terms']), (score,engine['baseline_terms'],len(item['terms']),item['excluded'])
        assert engine['excluded_terms'] == item['excluded'], (score,engine['excluded_terms'],item['excluded'])
        universe = b''.join(term['ordinal'].to_bytes(8,'little') + int(term['chrom']).to_bytes(1,'little') + term['pos'].to_bytes(4,'little') for term in item['terms'])
        assert engine['baseline_term_sha256'] == hashlib.sha256(universe).hexdigest()
        values = {sample:sum((term['values'][sample] for term in item['terms']),Fraction(0)) for sample in manifest['targets']+manifest['reference_samples']}
        baseline[score] = values
        base_vector = [values[sample] for sample in manifest['targets']]
        counts[score] = {'published_terms':SCORES[score][1],'independently_verified_baseline_terms':len(item['terms']), 'excluded':item['excluded']}
        for scenario_id,scenario in scenarios.items():
            scenario_values, scenario_counts = {}, collections.Counter()
            for sample in manifest['targets']:
                key = score,scenario_id,sample
                value = Fraction(0)
                masked = filled = unfilled = changed = 0
                regional_delta = Fraction(0)
                for term in item['terms']:
                    omitted = any(hidden(mask,term,sample,beds) for mask in scenario['masks'])
                    mask_checks += 1
                    before = term['values'][sample]
                    after = before
                    if omitted:
                        masked += 1
                        if scenario.get('fill_frequencies') and term['fill'] is not None:
                            filled += 1; after = term['fill']
                        else:
                            unfilled += 1; after = Fraction(0)
                    value += after
                    if after != before:
                        changed += 1
                        reason = 'masked_frequency_fill' if omitted and scenario.get('fill_frequencies') and term['fill'] is not None else 'masked'
                        expected_changes[(score,scenario_id,sample,term['ordinal'])] = (term,before,after,after-before,reason)
                        if term['chrom'] == '9' and 21500000 <= term['pos']-1 < 22500000:
                            regional_delta += after-before
                row = actual[key]
                assert Fraction(row['baseline_score']) == values[sample], ('baseline',key)
                assert Fraction(row['scenario_score']) == value, ('scenario',key,Fraction(row['scenario_score']),value)
                assert Fraction(row['delta']) == value-values[sample], ('delta',key)
                assert Fraction(row['top_terms_delta'])+Fraction(row['other_terms_delta']) == value-values[sample]
                assert [int(row[x]) for x in ['baseline_terms','retained_terms','masked_terms','filled_terms','unfilled_terms','changed_terms']] == [len(item['terms']),len(item['terms'])-masked,masked,filled,unfilled,changed], key
                assert int(row['unavailable_alternate_terms']) == 0
                assert near(row['baseline_percentile'],percentile(values[sample],base_vector))
                assert near(row['fixed_baseline_percentile'],percentile(value,base_vector))
                expected_regions[key] = regional_delta
                scenario_values[sample] = value
                scenario_counts.update(masked=masked,filled=filled,unfilled=unfilled)
            expected_scenario[(score,scenario_id)] = scenario_values
            vector = [scenario_values[sample] for sample in manifest['targets']]
            for sample in manifest['targets']:
                assert near(actual[score,scenario_id,sample]['reranked_percentile'],percentile(scenario_values[sample],vector))
            engine_scenario = next(row for row in engine['scenarios'] if row['id'] == scenario_id)
            for count in ['masked','filled','unfilled']:
                assert engine_scenario[count+'_observations'] == scenario_counts[count]
            m = engine_scenario['metrics']
            errors = [float(a-b) for a,b in zip(vector,base_vector)]
            assert near(m['mae'],sum(abs(e) for e in errors)/len(errors))
            assert near(m['rmse'],math.sqrt(sum(e*e for e in errors)/len(errors)))
            assert near(m['max_abs_change'],max(abs(e) for e in errors))
            rank_a = [percentile(x,base_vector) for x in base_vector]
            rank_b = [percentile(x,vector) for x in vector]
            assert near(m['spearman'],statistics.correlation(rank_a,rank_b))
            assert near(m['pearson'],statistics.correlation(list(map(float,base_vector)),list(map(float,vector))))
            fixed = [percentile(x,base_vector) for x in vector]
            shifts = [abs(a-b) for a,b in zip(rank_a,fixed)]
            assert near(m['mean_abs_percentile_shift'],sum(shifts)/len(shifts))
            assert near(m['max_abs_percentile_shift'],max(shifts))
            for crossing in m['crossings']:
                cutoff = crossing['cutoff']
                up = sum(a<cutoff<=b for a,b in zip(rank_a,fixed))
                down = sum(b<cutoff<=a for a,b in zip(rank_a,fixed))
                assert [crossing[x] for x in ['up','down','unchanged','evaluable']] == [up,down,len(vector)-up-down,len(vector)]
    top_rows = [json.loads(line) for line in (report_root / 'top-changes.jsonl').read_text().splitlines()]
    assert len(top_rows) == len(actual)
    grouped_changes = collections.defaultdict(list)
    for (score,scenario_id,sample,ordinal),change in expected_changes.items():
        grouped_changes[(score,scenario_id,sample)].append((ordinal,change))
    for row in top_rows:
        key = row['pgs_id'],row['scenario'],row['sample']
        expected = sorted(grouped_changes[key], key=lambda pair:(-abs(pair[1][3]),pair[0]))
        top = expected[:experiment['top_terms']]
        assert row['changed_terms'] == len(expected) and len(row['top']) == len(top)
        for actual_change,(ordinal,(term,before,after,delta,reason)) in zip(row['top'],top):
            assert actual_change['ordinal'] == ordinal
            assert (Fraction(actual_change['baseline']),Fraction(actual_change['scenario']),Fraction(actual_change['delta']),actual_change['reason']) == (before,after,delta,reason)
        top_delta = sum((change[3] for _,change in top),Fraction(0))
        all_delta = sum((change[3] for _,change in expected),Fraction(0))
        assert (Fraction(row['top_terms_delta']),Fraction(row['other_terms_delta']),Fraction(row['delta'])) == (top_delta,all_delta-top_delta,all_delta)
    changed_rows = read_tsv(report_root / 'terms.tsv.zst')
    assert len(changed_rows) == len(expected_changes)
    for row in changed_rows:
        key = row['pgs_id'],row['scenario'],row['sample'],int(row['ordinal'])
        term,before,after,delta,reason = expected_changes.pop(key)
        assert (row['contig'].removeprefix('chr'),int(row['pos1'])) == (term['chrom'],term['pos'])
        assert (Fraction(row['baseline']),Fraction(row['scenario_contribution']),Fraction(row['delta']),row['reason']) == (before,after,delta,reason)
    assert not expected_changes
    region_rows = read_tsv(report_root / 'regions.tsv')
    assert len(region_rows) == len(actual)
    for row in region_rows:
        assert row['region'] == 'chr9_dropout_region'
        assert Fraction(row['delta']) == expected_regions[(row['pgs_id'],row['scenario'],row['sample'])]
    reference_rows = read_tsv(report_root / 'references.tsv')
    assert len(reference_rows) == len(actual)*2
    for row in reference_rows:
        score,scenario_id,sample,panel = row['pgs_id'],row['scenario'],row['sample'],row['reference']
        reference = [baseline[score][name] for name in manifest['reference_panels'][panel]]
        before,after = percentile(baseline[score][sample],reference),percentile(expected_scenario[(score,scenario_id)][sample],reference)
        assert int(row['shared_sample_ids']) == 0 and int(row['reference_samples']) == 48
        assert near(row['baseline_percentile'],before) and near(row['scenario_percentile'],after)
        assert near(row['percentile_delta'],after-before)
    assert sha(pgsum) == identity['pgsum_sha256']
    out.mkdir()
    result = {'schema':'pgsum-public-robustness-validation-v1','passed':True,'models':counts,
              'target_samples':32,'disjoint_reference_samples':96,'reference_panels':2,'scenarios':len(scenarios),
              'sample_scenario_scores_checked_exactly':len(actual),'independent_mask_decisions':mask_checks,
              'nonzero_term_changes_checked_exactly':len(changed_rows),'top_attribution_rows_checked_exactly':len(top_rows),'regional_deltas_checked_exactly':len(region_rows),
              'external_reference_percentile_rows_checked':len(reference_rows),
              'oracle':'Python fractions.Fraction on original published weight strings, literal reference AF strings, FASTA orientation, and retained public VCF GTs. No pack, extracted genotype table, or engine contribution used as arithmetic input.',
              'aggregate_tolerance':'Absolute or relative 1e-12 for floating-point summary statistics and displayed percentiles; exact equality for scores and contributions.',
              'scope':'Four full published inventories, conditional common-callable autosomal baselines. Real genotype masking and score reconstruction only; no phenotype prediction, clinical calibration, ancestry normalization, or imputation-accuracy claim.',
              'preparation':manifest,'identity':identity,'native_report_sha256':sha(report_root/'report.json'),
              'report':report}
    save(out / 'results.json',result)
    print(json.dumps({key:value for key,value in result.items() if key not in ['preparation','report']},indent=2),flush=True)


def main():
    if not __debug__:
        raise SystemExit('Do not use -O: assertions are validation gates')
    getcontext().prec = 100
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('fetch'); p.add_argument('--out', required=True, type=Path)
    p = sub.add_parser('prepare'); p.add_argument('--data', required=True, type=Path); p.add_argument('--out', required=True, type=Path); p.add_argument('--reference', required=True, type=Path)
    p = sub.add_parser('native'); p.add_argument('--prepared', required=True, type=Path); p.add_argument('--out', required=True, type=Path); p.add_argument('--pgsum', required=True, type=Path)
    p = sub.add_parser('verify'); p.add_argument('--prepared', required=True, type=Path); p.add_argument('--native-root', required=True, type=Path); p.add_argument('--out', required=True, type=Path); p.add_argument('--pgsum', required=True, type=Path)
    args = parser.parse_args()
    kwargs = {key:(value.resolve() if isinstance(value,Path) else value) for key,value in vars(args).items() if key != 'command'}
    globals()[args.command](**kwargs)


if __name__ == '__main__':
    main()
