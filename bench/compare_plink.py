#!/usr/bin/env python3
"""Independent synthetic GT oracle and actual PLINK 2 / pgsum migration check.

No downloads. Supply both executables; a fresh output directory retains every input,
command, tool version, hash and comparison. This is not a public real-data benchmark.
"""
import argparse
import csv
from decimal import Decimal, getcontext
import hashlib
import json
from pathlib import Path
import subprocess

getcontext().prec = 100


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pgsum', type=Path, required=True)
    p.add_argument('--plink2', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    pgsum, plink = a.pgsum.resolve(), a.plink2.resolve()
    root = a.out.resolve()
    root.mkdir()
    commands = []
    binary_hashes = {'pgsum': sha(pgsum), 'plink2': sha(plink)}

    def run(args):
        args = list(map(str, args))
        commands.append(args)
        result = subprocess.run(args, cwd=root, capture_output=True, text=True)
        (root / f'command-{len(commands):02d}.log').write_text(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(f'{args}:\n{result.stdout}\n{result.stderr}')
        return result.stdout.strip()

    versions = {'pgsum': run([pgsum, '--version']), 'plink2': run([plink, '--version'])}
    samples, variants = 32, 48
    reference = root / 'reference.fa'
    reference.write_text('>1\n' + 'A' * variants + '\n')
    Path(str(reference) + '.fai').write_text(f'1\t{variants}\t3\t{variants}\t{variants + 1}\n')
    calls = [[None if s == samples - 1 or (v * 7 + s * 11) % 17 == 0 else (v + s) % 3
              for s in range(samples)] for v in range(variants)]
    vcf = root / 'cohort.vcf'
    with vcf.open('w') as f:
        f.write('##fileformat=VCFv4.2\n##contig=<ID=1,length=48>\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' +
                '\t'.join(f'S{s}' for s in range(samples)) + '\n')
        for v, row in enumerate(calls):
            f.write(f'1\t{v+1}\tv{v}\tA\tG\t.\tPASS\t.\tGT\t' +
                    '\t'.join({None: './.', 0: '0/0', 1: '0/1', 2: '1/1'}[x] for x in row) + '\n')
    models = ['ADDITIVE', 'DOMINANT', 'RECESSIVE', 'PRECISE']
    expected = {}
    for model in models:
        values = [Decimal(0)] * samples
        with (root / f'{model}.txt').open('w') as f, (root / f'{model}.weights').open('w') as w:
            f.write(f'#pgs_id={model}\n#genome_build=GRCh38\n'
                    'chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\tis_dominant\tis_recessive\n')
            w.write(f'ID\tALLELE\t{model}\n')
            for v, row in enumerate(calls):
                ref_effect = v % 3 == 0
                effect, other = ('A', 'G') if ref_effect else ('G', 'A')
                weight = (Decimal('0.123456789123456789') if model == 'PRECISE'
                          else Decimal((v % 7) - 3) / Decimal(8))
                f.write(f'1\t{v+1}\t{effect}\t{other}\t{weight}\t{model == "DOMINANT"}\t{model == "RECESSIVE"}\n')
                w.write(f'v{v}\t{effect}\t{weight}\n')
                for s, call in enumerate(row):
                    if call is None:
                        continue  # Neither tool imputes absent calls in this check.
                    dosage = 2 - call if ref_effect else call
                    multiplier = int(dosage > 0) if model == 'DOMINANT' else int(dosage == 2) if model == 'RECESSIVE' else dosage
                    values[s] += multiplier * weight
        expected[model] = values
    run([pgsum, '--threads', '2', 'compile', *[root / f'{m}.txt' for m in models],
         '--reference', reference, '--out', root / 'packs'])
    run([pgsum, '--threads', '2', 'extract', '--gvcf', vcf, '--reference', reference,
         '--pack', root / 'packs', '--all-samples', '--accept-missing-quality', '--out', root / 'cohort.pgsc'])
    results = {}
    for model in models:
        run([pgsum, '--threads', '2', 'score', '--genotypes', root / 'cohort.pgsc',
             '--pack', root / 'packs' / f'{model}.pgsp', '--out', root / model])
        modifiers = [model.lower()] if model in ['DOMINANT', 'RECESSIVE'] else []
        run([plink, '--vcf', vcf, '--score', root / f'{model}.weights', '1', '2', '3',
             'header-read', 'no-mean-imputation', 'cols=+scoresums', *modifiers,
             '--threads', '2', '--out', root / f'plink-{model}'])
        tolerance = ['--absolute-tolerance', '0.00001', '--relative-tolerance', '0.000005'] if model == 'PRECISE' else []
        comparison = root / f'compare-{model}'
        run([pgsum, 'compare', '--left', root / model / 'cohort-scores.tsv.zst',
             '--right', root / f'plink-{model}.sscore', '--right-format', 'plink2',
             '--right-score-column', f'{model}_SUM', '--right-score-id', model,
             '--out', comparison, '--fail-on-difference', *tolerance])
        rows = list(csv.DictReader((comparison / 'scores.tsv').open(), delimiter='\t'))
        assert len(rows) == samples
        for row in rows:
            s = int(row['sample'][1:])
            assert Decimal(row['left']) == expected[model][s], (model, s, row, expected[model][s])
            # Independent check of PLINK's rounded output, separate from pgsum compare.
            actual = Decimal(row['right'])
            limit = Decimal('0.00001') + Decimal('0.000005') * max(abs(actual), abs(expected[model][s])) if model == 'PRECISE' else Decimal(0)
            assert abs(actual - expected[model][s]) <= limit, (model, s, row)
        results[model] = json.loads((comparison / 'comparison.json').read_text())
    assert binary_hashes == {'pgsum': sha(pgsum), 'plink2': sha(plink)}, 'executable changed during validation'
    summary = {
        'schema': 'pgsum-plink-synthetic-migration-v1', 'versions': versions,
        'binaries_sha256': binary_hashes,
        'runner_sha256': sha(__file__), 'samples': samples, 'variants': variants,
        'oracle': 'Independent Python Decimal: 128 sample/score sums, ref and alt effect alleles, missing calls, three GT models, high precision weights.',
        'limitations': 'Synthetic biallelic diploid autosomes only; no public real-data, DS/GP, gVCF, indel, multiallelic, sex-chromosome or pgsc_calc execution validation.',
        'inputs_sha256': {path.name: sha(path) for path in root.iterdir() if path.suffix in ['.txt', '.vcf', '.fa', '.fai', '.weights']},
        'comparisons': results, 'commands': commands,
    }
    (root / 'results.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({m: r['counts'] for m, r in results.items()}, indent=2))


if __name__ == '__main__':
    main()
