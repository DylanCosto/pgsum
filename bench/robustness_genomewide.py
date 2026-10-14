#!/usr/bin/env python3
"""Whole-inventory PGS000018 public-genotype robustness and scaling validation.

Preparation streams chromosome records and uses sorted position indexes. The independent
oracle reads original allele/weight strings and public GTs, never engine contributions.
"""
import argparse
import bisect
import collections
import concurrent.futures
import csv
from decimal import Decimal
from fractions import Fraction
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import urllib.request

import numpy as np
import pysam
import pysam.bcftools

SCORE = 'PGS000018'
SCORE_SHA256 = 'ed5dcd5359213364a182dd7972edf524cc3e161f8332ca54ca5515c95cbf8d2a'
SCORE_ROWS = 1745179
VCF_URL = ('https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/'
           'ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz')
SCORE_URL = 'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000018/ScoringFiles/Harmonized/PGS000018_hmPOS_GRCh38.txt.gz'
METADATA_URL = 'https://www.pgscatalog.org/rest/score/PGS000018'
METADATA_SHA256 = 'eac884834abaeb0bf07d630ccfdcb6b9f6b46fa69fa6703741efe4061620edd7'
SOURCE_SHA256 = {'1': ('5598cad22ca866457eb9ff5baeeade6e47e2ebd2b34d1dc2a180cf6f399f5044',
       '057d36e5edfb8abe26051c00bdbcef4c1d5ef89c4dda77edd1dd3eabd2ac58be'),
 '2': ('8c2a40e796e6a47855d4b825754ba01180f855286f5f2703f2ea1c849e1500c3',
       '64b4d71c02db97eed35363942df0fc42362ad122304d62657a5c07cb55dfa4cd'),
 '3': ('f8819c823e66edc4f5f00bbafbd176dfa65ef5be74d2afac85db88707bb1899b',
       'c5dc443576576ff08c36bbc6aaf0c5f23ece8c5cd262278c86657c6d9e169be5'),
 '4': ('81b0d81476d750dfb62d368062b4a0f60ccea80d197f3aa1990fd4c3be3fe55a',
       '5dcec0451c2d5f185a160f5d6a2828d5c71ac4e784970776621150982985fdb4'),
 '5': ('f36cb8ae8e47e297fa60ac5814bf145dede1b84f985e5e60a3dd8935777d0fd3',
       '7f78818c4a290451769df0e702aef7e53ffc581da609510b9f7bce86b947ec23'),
 '6': ('5432a12a8d63b9fba217156292a56257f2afdc65bf61ff7d60d15797a7c609ca',
       '6bf30b4d84c4879829a878a5e6298cb60550274ccce816b700c751d5dac5b36f'),
 '7': ('9b3812d5f7db9d674785d9607b87e9bb55ca5da149f171b57cbf1ee649bea673',
       '3996821957644808feb23608b09082ece91be986b2c947565253fd4574a6225f'),
 '8': ('0856971ddbc278c3e4c3209a90aa17112139bb7591f042aa5aba100c57d5edda',
       '599d7fb91e1b6d39e3133a75a4736ebe1566b33fc1a59e3b7daf364e5a79de54'),
 '9': ('ed72c73985aba82582959d31075a355905e69cac341366b85d7acf7f2205e262',
       '23b91c1f4e838387262e0d08a7b03f91fff665e003cc32f03ddf454f2520a10c'),
 '10': ('17e067b902ca9a66505c1fd6ce88f170f80cf8c6572c1f8ed874d26c4a4b8a57',
        '4e2984374546a95b703197ef83e74f6378fa4df6958d2f9ed6ccfa516809efcb'),
 '11': ('f7a273a73035ee422f243419ae1169cba97a8c1e40932b78dc2bc7ce01e6f081',
        'bb51925586e69cdc7594bb2891db9c251cbe9e17451f854120560920c955a426'),
 '12': ('bf0120b7ca8705f747dec6b5aafe2ed6b1d858185fa14763f01cf403145a5e5e',
        'c7016abd6121a8a676c424539ff12124c89788b99b7fc87b4d724e182f5f2a3a'),
 '13': ('7f404c013e70497d2722ad2a54dfb23bbad4c495473c8de7e8a0a5edaff0928e',
        '75216d24801dbadbfcbe9683eec55edf70a35e603f3d5c614c04930d4b225f05'),
 '14': ('bdd4f4b911c945971d4927a5890613fadfbdf90560929738934a9d75cf80ed76',
        '7314f2341fe96fa1ed31c456cf1a4626593f6d215769c00032d52ad5f0141508'),
 '15': ('0b694817f4c6fd1e780c54a50a4c48fff81e2236968ccf2cc5b753802c56a98c',
        '7354c23d8df53fe56fa5128f370e668734abc53577a57dd92635553392bc2b5f'),
 '16': ('74a6ab9f91ec15edb2bf5a3cb5aa630240adaef575df62515c2a52009193367c',
        'd233a28c97ebfe8f6987cc6649236d6d6d37a8b559e207342159ef2ddda913bd'),
 '17': ('27542c440b24e5ab39e429dd897b0ced89dd22041051e1e1830c753fd260355c',
        'c41e59d65a5615421c5a7ef15d829143d1e22a34d4807576d7fb7053f079ce12'),
 '18': ('5446ec34682f707f41a388e8444c5f3b0c3023f363918cea9d20e2b8a39e2ac0',
        'b4dde6835e53d9f841226ca70c6d26c2b9f9e56e6c0c516a909e41eedc8fa211'),
 '19': ('2cb80e80fc6f843d384abbe1f150ee9a0bcfd5a4a4b348dcc2b52a3633de0b62',
        '193c904b39182a3f040806f3c6c871f0289513271eada67683007c4ffc2e3144'),
 '20': ('0629e7895185e5082d6baba1a3965156d17635dd252de05ced350e2075cf4f48',
        '2f970da60a6e23905027e9f2c468a9c80702d295e24ccafcf242f57a31689c15'),
 '21': ('44d9f3aab459edf80ed51f069ad0acd77f5859f3d1101eabd435b4f73626ee1b',
        '3888183c963b6782a171a0d01a92a37f2350e21b6348ff0c96879c189e05026c'),
 '22': ('dc2b718c00d1a8220f19561bc647eea08fc25ce641c8f2c12e540fb31f073f88',
        '643106c3ddcd0fe95538aa0be8255873339f034b8f40137213e92130283b8b67')}
REFERENCE_SHA256 = '9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7'
COHORT_SIZES = (8, 32, 128)
REFERENCE_SIZE = 64


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def download(url, path, expected=None):
    if path.exists():
        if expected and sha(path) != expected:
            raise ValueError(f'cached source changed: {path}')
        return
    temporary = path.with_suffix(path.suffix + '.part')
    digest = hashlib.sha256()
    started = time.monotonic()
    with urllib.request.urlopen(url, timeout=120) as source, temporary.open('xb') as target:
        expected_size = int(source.headers.get('Content-Length', 0))
        written = 0
        while chunk := source.read(1 << 20):
            target.write(chunk)
            digest.update(chunk)
            written += len(chunk)
    assert not expected_size or written == expected_size, (url, written, expected_size)
    assert expected is None or digest.hexdigest() == expected, f'upstream content changed: {url}'
    os.replace(temporary, path)
    print(f'downloaded {path.name}: {written:,} bytes in {time.monotonic()-started:.1f}s', flush=True)


def fetch(out, reuse=None):
    out.mkdir(parents=True, exist_ok=True)
    for name,url,expected in [(SCORE+'.txt.gz',SCORE_URL,SCORE_SHA256), (SCORE+'.metadata.json',METADATA_URL,METADATA_SHA256)]:
        path = out/name
        if reuse and (reuse/name).exists() and not path.exists():
            shutil.copyfile(reuse/name,path)
        download(url,path,expected)
    if reuse and (reuse/'chr22.vcf.gz').exists() and not (out/'chr22.vcf.gz').exists():
        (out/'chr22.vcf.gz').symlink_to((reuse/'chr22.vcf.gz').resolve())
    def chromosome(chrom):
        url = VCF_URL.format(chrom=chrom)
        path = out/f'chr{chrom}.vcf.gz'
        download(url,path,SOURCE_SHA256[str(chrom)][0])
        download(url+'.tbi',out/f'chr{chrom}.vcf.gz.tbi',SOURCE_SHA256[str(chrom)][1])
        return str(chrom),{'url':url,'bytes':path.stat().st_size,'sha256':sha(path),
                           'index_sha256':sha(out/f'chr{chrom}.vcf.gz.tbi')}
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        chromosomes = dict(pool.map(chromosome, range(1,23)))
    save(out/'sources.json', {'chromosomes':chromosomes,'score':{'url':SCORE_URL,'sha256':sha(out/(SCORE+'.txt.gz'))},
                             'catalog_metadata':{'url':METADATA_URL,'sha256':sha(out/(SCORE+'.metadata.json'))}})


def timed(out, label, command):
    command = list(map(str, command))
    timing = out / (label + '.time.json')
    with (out / (label + '.log')).open('w') as log:
        result = subprocess.run(['/usr/bin/time','-f','{"elapsed_seconds":%e,"peak_rss_kib":%M,"user_seconds":%U,"system_seconds":%S}',
                                 '-o',str(timing),*command],stdout=log,stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'{label} failed; inspect {out / (label + ".log")}')
    measurement = json.loads(timing.read_text())
    measurement['command'] = command
    return measurement


def score_rows(path):
    with gzip.open(path,'rt') as stream:
        yield from csv.DictReader((line for line in stream if not line.startswith('#')),delimiter='\t')


def _subset(source, samples, bed, out):
    options = ['-S',str(samples),'--no-update','--no-version','-Ob','-o',str(out)]
    if bed:
        options += ['-T',str(bed),'--targets-overlap','1']
    pysam.bcftools.view(*options,str(source),catch_stdout=False)
    pysam.bcftools.index('-f',str(out),catch_stdout=False)


def _oracle_chrom(model, genotypes, reference, selection, out):
    """Sweep original records against sorted source loci; store exact source-derived dosages."""
    data = np.load(model)
    pos, ordinals, effects, others = (data[k] for k in ('pos','ordinal','effect','other'))
    chrom = int(data['chrom'][0])
    names = json.loads(selection.read_text())
    samples = names['selected']
    with pysam.VariantFile(str(genotypes)) as source:
        assert list(source.header.samples) == samples, 'source genotype column order changed'
    ref_indices = np.array([samples.index(s) for s in names['references']],dtype=np.int64)
    query = out / f'chr{chrom}.gt.tsv'
    pysam.bcftools.query('-f','%POS\t%END\t%REF\t%ALT\t%FILTER[\t%GT]\n','-o',str(query),str(genotypes),catch_stdout=False)
    with pysam.FastaFile(str(reference)) as fasta:
        name = f'chr{chrom}' if f'chr{chrom}' in fasta.references else str(chrom)
        bases = np.array([ord(fasta.fetch(name,int(p)-1,int(p)).upper()) for p in pos],dtype=np.uint8)
    complement = np.arange(256,dtype=np.uint8)
    for before,after in zip(b'ACGT',b'TGCA'):
        complement[before] = after
    direct = (bases==effects)|(bases==others)
    reverse = (bases==complement[effects])|(bases==complement[others])
    oriented = np.logical_xor(direct,reverse) & np.isin(effects,list(b'ACGT')) & np.isin(others,list(b'ACGT')) & (effects != others)
    effect = np.where(direct,effects,complement[effects])
    alternate = np.where(bases==effect,np.where(direct,others,complement[others]),effect)
    calls = np.full((len(pos),len(samples)),255,dtype=np.uint8)
    fill_count = np.full(len(pos),-1,dtype=np.int16)
    frequency_path = out / f'chr{chrom}.frequencies.vcf.gz'
    records_seen = 0
    with gzip.open(frequency_path,'wt') as frequency:
        frequency.write('##fileformat=VCFv4.2\n##INFO=<ID=REF_AF,Number=A,Type=Float,Description="ALT frequency from 64 disjoint reference people">\n')
        frequency.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        def records():
            nonlocal records_seen
            with query.open('rb') as stream:
                for line in stream:
                    p,end,ref,alt,flt,gt = line.rstrip(b'\n').split(b'\t',5)
                    # This pinned diploid cohort uses single-digit allele indexes. Fail closed on other encodings.
                    assert len(gt)+1 == 4*len(samples), 'oracle requires diploid single-digit GT allele indexes'
                    chars = np.frombuffer(gt+b'\t',dtype=np.uint8).reshape(-1,4)
                    assert len(chars)==len(samples) and np.all(np.isin(chars[:,1],[ord('|'),ord('/')])) and np.all(chars[:,3]==9)
                    a,b = chars[:,0].astype(np.int16)-48,chars[:,2].astype(np.int16)-48
                    alts = alt.split(b',')
                    assert np.all((a==-2)|((a>=0)&(a<=len(alts)))) and np.all((b==-2)|((b>=0)&(b<=len(alts))))
                    reference_calls = np.concatenate([a[ref_indices],b[ref_indices]])
                    reference_calls = reference_calls[reference_calls>=0]
                    counts = [int(np.count_nonzero(reference_calls==i)) for i in range(1,len(alts)+1)]
                    af = [str(Decimal(count)/Decimal(len(reference_calls))) if len(reference_calls) else '.' for count in counts]
                    frequency.write(f'{chrom}\t{int(p)}\t.\t{ref.decode()}\t{alt.decode()}\t.\tPASS\tREF_AF={",".join(af)}\n')
                    records_seen += 1
                    yield (int(p),int(end),ref,alts,flt,a,b,reference_calls)
        iterator = iter(records())
        pending = next(iterator,None)
        active = []
        for i,p in enumerate(pos):
            p = int(p)
            active = [r for r in active if r[1]>=p]
            while pending is not None and pending[0]<=p:
                if pending[1]>=p:
                    active.append(pending)
                pending = next(iterator,None)
            if not oriented[i] or len(active)!=1:
                continue
            start,end,ref,alts,flt,a,b,reference_calls = active[0]
            if start!=p or len(ref)!=1 or ref[0]!=int(bases[i]) or flt not in (b'PASS',b'.'):
                continue
            matching = [j+1 for j,allele in enumerate(alts) if allele==bytes([int(alternate[i])])]
            assert len(matching)<=1
            target_alt = matching[0] if matching else -1
            good = ((a==0)|(a==target_alt)) & ((b==0)|(b==target_alt))
            dosage = (a==target_alt).astype(np.uint8)+(b==target_alt).astype(np.uint8)
            if int(effect[i])==int(bases[i]):
                dosage = 2-dosage
            calls[i,good] = dosage[good]
            if matching and len(alts)==1 and len(reference_calls)==2*REFERENCE_SIZE:
                count = int(np.count_nonzero(reference_calls==target_alt))
                fill_count[i] = count if int(effect[i])!=int(bases[i]) else 2*REFERENCE_SIZE-count
        # Includes any context record yielded after the final source locus, if present.
        if pending is not None:
            for _ in iterator:
                pass
    np.savez_compressed(out/f'chr{chrom}.oracle.npz',ordinal=ordinals,pos=pos,weights=data['weight'],
                        oriented=oriented,calls=calls,fill_effect_alleles=fill_count)
    query.unlink()
    save(out/f'chr{chrom}.oracle.json',{'chromosome':chrom,'source_loci':len(pos),'records':records_seen,
        'oriented_terms':int(oriented.sum()),'common_callable_192':int((oriented & np.all(calls!=255,axis=1)).sum()),
        'oracle_sha256':sha(out/f'chr{chrom}.oracle.npz'),'frequency_sha256':sha(frequency_path)})


def prepare(sources, reference, out):
    out.mkdir()
    assert sha(sources/(SCORE+'.txt.gz'))==SCORE_SHA256
    assert sha(sources/(SCORE+'.metadata.json'))==METADATA_SHA256
    shutil.copyfile(sources/(SCORE+'.txt.gz'),out/(SCORE+'.txt.gz'))
    shutil.copyfile(sources/(SCORE+'.metadata.json'),out/(SCORE+'.metadata.json'))
    # Reference and evaluation IDs are selected before reading any score values or phenotypes.
    with pysam.VariantFile(str(sources/'chr22.vcf.gz'),index_filename=str(sources/'chr22.vcf.gz.tbi')) as source:
        all_samples = list(source.header.samples)
    selected = [all_samples[i*(len(all_samples)-1)//191] for i in range(192)]
    refs = selected[2::3]
    targets = [sample for i,sample in enumerate(selected) if i%3 != 2]
    assert len(targets)==128 and len(refs)==64 and set(targets).isdisjoint(refs)
    selection = {'source_samples':len(all_samples),'selected':selected,'references':refs,'targets':targets,
                 'method':'192 evenly spaced VCF-header indices floor(i*(n-1)/191); index%3==2 supplies 64 fixed references; other 128 are targets; cohorts use first8/32/128 targets'}
    save(out/'selection.json',selection)
    (out/'samples192.txt').write_text('\n'.join(selected)+'\n')
    for size in COHORT_SIZES:
        (out/f'targets{size}.txt').write_text('\n'.join(targets[:size])+'\n')
    (out/'reference64.txt').write_text('\n'.join(refs)+'\n')
    buckets = collections.defaultdict(list)
    total_abs_weight = 0
    unplaced = []
    for index,row in enumerate(score_rows(out/(SCORE+'.txt.gz'))):
        value = Fraction(row['effect_weight'])*10**12
        assert value.denominator==1
        weight = value.numerator
        total_abs_weight += abs(weight)
        assert len(row['effect_allele'])==len(row['other_allele'])==1
        chrom = row['hm_chr']
        if chrom not in {str(i) for i in range(1,23)} or not row['hm_pos'].isdigit():
            unplaced.append(index+1)
            continue
        buckets[int(chrom)].append((index+1,int(row['hm_pos']),ord(row['effect_allele']),ord(row['other_allele']),weight))
    assert index+1==SCORE_ROWS
    assert total_abs_weight*2*REFERENCE_SIZE < np.iinfo(np.int64).max, 'fixed-point oracle overflow bound exceeded'
    for chrom,rows in buckets.items():
        rows.sort(key=lambda row:(row[1],row[0]))
        values = np.array(rows,dtype=np.int64)
        np.savez_compressed(out/f'chr{chrom}.model.npz',chrom=np.array([chrom],dtype=np.uint8),ordinal=values[:,0].astype(np.uint32),
                            pos=values[:,1].astype(np.uint32),effect=values[:,2].astype(np.uint8),other=values[:,3].astype(np.uint8),weight=values[:,4])
        positions = sorted(set(row[1] for row in rows))
        (out/f'chr{chrom}.bed').write_text(''.join(f'{chrom}\t{p-1}\t{p}\n' for p in positions))
    script = Path(__file__).resolve()
    def chromosome(chrom):
        source = sources/f'chr{chrom}.vcf.gz'
        deadline = time.monotonic()+3600
        while not source.exists() or not (sources/f'chr{chrom}.vcf.gz.tbi').exists():
            if time.monotonic()>deadline:
                raise RuntimeError(f'public download not ready after one hour: {source}')
            time.sleep(1)
        subset = out/f'chr{chrom}.bcf'
        measurement = timed(out,f'chr{chrom}-subset',[sys.executable,script,'_subset','--source',source,'--samples',out/'samples192.txt','--bed',out/f'chr{chrom}.bed','--out',subset])
        oracle_measurement = timed(out,f'chr{chrom}-oracle',[sys.executable,script,'_oracle_chrom','--model',out/f'chr{chrom}.model.npz','--genotypes',subset,
            '--reference',reference,'--selection',out/'selection.json','--out',out])
        print(f'prepared chr{chrom}: subset {measurement["elapsed_seconds"]:.1f}s; oracle {oracle_measurement["elapsed_seconds"]:.1f}s',flush=True)
        return str(chrom),{'subset':measurement,'oracle':oracle_measurement,**json.loads((out/f'chr{chrom}.oracle.json').read_text())}
    # Four C-level streaming readers, each with bounded chromosome-sized oracle arrays.
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        measurements = dict(pool.map(chromosome,range(1,23)))
    combined = out/'public192.bcf'
    pysam.bcftools.concat('--no-version','-Ob','-o',str(combined),*[str(out/f'chr{c}.bcf') for c in range(1,23)],catch_stdout=False)
    pysam.bcftools.index('-f',str(combined),catch_stdout=False)
    for label in [f'targets{size}' for size in COHORT_SIZES]+['reference64']:
        timed(out,label+'-subset',[sys.executable,script,'_subset','--source',combined,'--samples',out/(label+'.txt'),'--out',out/(label+'.bcf')])
    manifest = {'schema':'pgsum-genomewide-preparation-v1','score':SCORE,'published_file_rows':SCORE_ROWS,'file_declared_rows':1745180,
                'unplaced_source_ordinals':unplaced,'selection':selection,'weight_scale':10**12,'exact_oracle_scale':10**12*REFERENCE_SIZE,
                'absolute_weight_sum_scaled':total_abs_weight,'fixed_point_overflow_bound':total_abs_weight*2*REFERENCE_SIZE,
                'reference':{'path':str(reference),'sha256':sha(reference),'fai_sha256':sha(str(reference)+'.fai')},
                'sources':json.loads((sources/'sources.json').read_text()),'measurements':measurements,
                'versions':{'python':sys.version,'pysam':pysam.__version__,'numpy':np.__version__,'bcftools':pysam.__samtools_version__},
                'scope':'Full published scoring-file inventory. Common-callable baseline re-evaluated separately at each target cohort size with fixed64 disjoint references; no phenotype or clinical calibration.',
                'outputs_sha256':{p.name:sha(p) for p in out.iterdir() if p.is_file()}}
    save(out/'preparation.json',manifest)
    print('whole-genome preparation complete',flush=True)


def verify_prepared(prepared, reference_files=False):
    manifest = json.loads((prepared/'preparation.json').read_text())
    assert manifest['schema']=='pgsum-genomewide-preparation-v1'
    assert manifest['reference']['sha256']==REFERENCE_SHA256
    assert manifest['sources']['catalog_metadata']['sha256']==METADATA_SHA256
    assert manifest['sources']['score']['sha256']==SCORE_SHA256
    for chrom,source in manifest['sources']['chromosomes'].items():
        assert (source['sha256'],source['index_sha256'])==SOURCE_SHA256[chrom]
    for name,digest in manifest['outputs_sha256'].items():
        assert sha(prepared/name)==digest, f'prepared input changed: {name}'
    if reference_files:
        assert sha(manifest['reference']['path'])==manifest['reference']['sha256']
        assert sha(manifest['reference']['path']+'.fai')==manifest['reference']['fai_sha256']
    return manifest


def native_inputs(prepared, out, pgsum):
    manifest = verify_prepared(prepared,reference_files=True)
    out.mkdir()
    identity = {'pgsum_sha256':sha(pgsum),'version':subprocess.check_output([pgsum,'--version'],text=True).strip(),
                'prepared_manifest_sha256':sha(prepared/'preparation.json')}
    measurements = {}
    # Validate source-specific assumptions for every chromosome, independently of cached oracle arrays.
    for chrom in range(1,23):
        with pysam.VariantFile(str(prepared/f'chr{chrom}.bcf')) as source:
            assert list(source.header.samples)==manifest['selection']['selected']
            assert set(source.header.formats)=={'GT'}, 'oracle is scoped to this GT-only public cohort'
        assert not pysam.bcftools.query('-i','N_ALT!=1','-f','%CHROM:%POS\n',str(prepared/f'chr{chrom}.bcf')), 'pinned source is expected to be biallelic'
        assert not pysam.bcftools.query('-i','ALT~"[^ACGT]" || REF~"[^ACGT]"','-f','%CHROM:%POS\n',str(prepared/f'chr{chrom}.bcf')), 'oracle requires literal DNA alleles, without symbolic or spanning-deletion records'
    measurements['compile'] = timed(out,'compile',[pgsum,'--threads','4','compile',prepared/(SCORE+'.txt.gz'),'--reference',manifest['reference']['path'],'--out',out/'packs'])
    for label in [f'targets{s}' for s in COHORT_SIZES]+['reference64']:
        measurements['extract_'+label] = timed(out,'extract-'+label,[pgsum,'--threads','4','extract','--gvcf',prepared/(label+'.bcf'),
            '--reference',manifest['reference']['path'],'--pack',out/'packs','--all-samples','--accept-missing-quality','--out',out/(label+'.pgsc')])
    measurements['frequencies'] = timed(out,'frequencies',[pgsum,'--threads','4','frequencies','--vcf',*[prepared/f'chr{c}.frequencies.vcf.gz' for c in range(1,23)],'--field','REF_AF','--out',out/'reference.pgsf'])
    (out/'9p21-dropout.bed').write_text('9\t21500000\t22500000\n')
    assert sha(pgsum)==identity['pgsum_sha256'], 'executable changed during input construction'
    assert sha(manifest['reference']['path'])==manifest['reference']['sha256']
    assert sha(manifest['reference']['path']+'.fai')==manifest['reference']['fai_sha256']
    identity['measurements'] = measurements
    identity['public_input_gates'] = {'all_chromosome_sample_orders_verified':True,'biallelic_literal_DNA_records_only':True,
                                    'diploid_single_digit_GT_parser_fail_closed':True}
    identity['generated_inputs_sha256'] = {str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file() and p.suffix in ('.pgsp','.pgsc','.pgsf','.bed')}
    save(out/'identity.json',identity)
    print('native genotype and frequency inputs complete',flush=True)


def native(prepared, out, pgsum, cohort_size):
    manifest = json.loads((prepared/'preparation.json').read_text())
    identity = json.loads((out/'identity.json').read_text())
    assert identity['prepared_manifest_sha256']==sha(prepared/'preparation.json')
    assert sha(pgsum)==identity['pgsum_sha256'], 'use one executable for all measurements'
    for name,digest in identity['generated_inputs_sha256'].items():
        assert sha(out/name)==digest, f'native input changed: {name}'
    scenarios = [{'id':'no_mask','masks':[]}]
    for scenario in [
        {'id':'random20','masks':[{'kind':'random','rate':0.2,'seed':42,'unit':'site'}]},
        {'id':'9p21_dropout','masks':[{'kind':'bed','path':'../9p21-dropout.bed','mode':'drop'}]},
    ]:
        scenarios += [scenario,dict(scenario,id=scenario['id']+'_fill',fill_frequencies='../reference.pgsf',
            frequency_source_note='ALT frequencies computed only from the 64 named disjoint reference people in preparation.json. All evaluation target IDs are excluded; eligible baseline sites have 128 called reference alleles.')]
    run = out/f'cohort{cohort_size}'
    run.mkdir()
    experiment = {'schema':'pgsum-robustness-experiment-v1','genotypes':f'../targets{cohort_size}.pgsc',
        'packs':['../'+str(p.relative_to(out)) for p in sorted((out/'packs').glob('*.pgsp'))],
        'scoring':{'allow_inferred_other_allele':False,'allow_inferred_palindromes':False,'allow_inferred_indels':False,'accept_informational_descriptions':False},
        'cutoffs':[90,95],'scenarios':scenarios,'reference_sets':[{'id':'fixed_reference64','genotypes':'../reference64.pgsc'}],
        'regions':[{'id':'9p21','path':'../9p21-dropout.bed'}],'top_terms':10}
    save(run/'experiment.json',experiment)
    measured = timed(run,'robustness',[pgsum,'--threads','4','robustness','--experiment',run/'experiment.json','--out',run/'robustness'])
    assert sha(pgsum)==identity['pgsum_sha256']
    save(run/'measurement.json',measured)
    print(f'native cohort{cohort_size}: {measured["elapsed_seconds"]}s, {measured["peak_rss_kib"]} KiB RSS',flush=True)


def read_tsv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream,delimiter='\t'))


def site_mask(chrom, positions):
    prefix = b'pgsum-robustness-mask-sha256-v1\0'+(42).to_bytes(8,'little')+bytes([chrom])
    return np.fromiter((int.from_bytes(hashlib.sha256(prefix+int(p).to_bytes(4,'little')+b'\0').digest()[:8],'little')>>11 < 0.2*(1<<53) for p in positions),dtype=bool,count=len(positions))


def top_candidates(delta, ordinals, limit=10):
    nonzero = np.flatnonzero(delta)
    if len(nonzero)>limit:
        magnitude = np.abs(delta[nonzero])
        threshold = np.partition(magnitude,-limit)[-limit]
        nonzero = nonzero[magnitude>=threshold]
    return sorted((int(i) for i in nonzero),key=lambda i:(-abs(int(delta[i])),int(ordinals[i])))[:limit]


def source_arithmetic(prepared, cohort_size):
    """Source-derived fixed-point oracle: no pack, PGSC, or native contribution input."""
    manifest = json.loads((prepared/'preparation.json').read_text())
    selection = manifest['selection']
    samples = selection['selected']
    target_cols = np.array([samples.index(s) for s in selection['targets'][:cohort_size]])
    ref_cols = np.array([samples.index(s) for s in selection['references']])
    assert set(target_cols).isdisjoint(ref_cols)
    scenarios = ['no_mask','random20','random20_fill','9p21_dropout','9p21_dropout_fill']
    baseline = np.zeros(cohort_size,dtype=np.int64)
    refs = np.zeros(REFERENCE_SIZE,dtype=np.int64)
    results = {s:{'delta':np.zeros(cohort_size,dtype=np.int64),'region':np.zeros(cohort_size,dtype=np.int64),
        'masked':0,'filled':0,'changed':np.zeros(cohort_size,dtype=np.int64),'top':[[] for _ in range(cohort_size)],
        'fill_omissions':collections.Counter()} for s in scenarios}
    excluded = collections.Counter({'unresolved_position':len(manifest['unplaced_source_ordinals'])})
    universe = np.zeros(SCORE_ROWS,dtype=[('ordinal','<u8'),('contig','u1'),('pos','<u4')])
    present = np.zeros(SCORE_ROWS,dtype=bool)
    chromosome_terms = {}
    eligible_count,weight_sum = 0,0
    for chrom in range(1,23):
        data = np.load(prepared/f'chr{chrom}.oracle.npz')
        model = np.load(prepared/f'chr{chrom}.model.npz')
        oriented,calls = data['oriented'],data['calls']
        target_called = np.all(calls[:,target_cols]!=255,axis=1)
        reference_called = np.all(calls[:,ref_cols]!=255,axis=1)
        excluded['model_or_orientation_unscorable'] += int((~oriented).sum())
        excluded['baseline_target_not_common_callable'] += int((oriented & ~target_called).sum())
        excluded['reference_not_common_callable'] += int((oriented & target_called & ~reference_called).sum())
        eligible = oriented & target_called & reference_called
        ordinal,pos,weights = (data[k][eligible] for k in ('ordinal','pos','weights'))
        effects,others = model['effect'][eligible],model['other'][eligible]
        fill_count = data['fill_effect_alleles'][eligible]
        fillable = fill_count>=0
        target = calls[eligible][:,target_cols].astype(np.int64)
        reference = calls[eligible][:,ref_cols].astype(np.int64)
        universe['ordinal'][ordinal-1]=ordinal
        universe['contig'][ordinal-1]=chrom
        universe['pos'][ordinal-1]=pos
        present[ordinal-1]=True
        chromosome_terms[str(chrom)]=len(ordinal)
        eligible_count += len(ordinal)
        weight_sum += sum(map(int,np.abs(weights)))
        # Bound checked from all original absolute weights in prepare(), before any int64 arithmetic.
        contributions = weights[:,None]*target*REFERENCE_SIZE
        baseline += contributions.sum(axis=0)
        refs += (weights[:,None]*reference*REFERENCE_SIZE).sum(axis=0)
        random = site_mask(chrom,pos)
        region = (chrom==9)&(pos>21500000)&(pos<=22500000)
        masks = {'no_mask':np.zeros(len(pos),dtype=bool),'random20':random,'9p21_dropout':region}
        for scenario,result in results.items():
            with_fill = scenario.endswith('_fill')
            mask = masks[scenario.removesuffix('_fill')]
            result['masked'] += int(mask.sum())
            filled = mask & fillable if with_fill else np.zeros(len(pos),dtype=bool)
            result['filled'] += int(filled.sum())
            # A callable unique SNV without a matching biallelic frequency is an allele mismatch.
            # Multi-ALT records are not present in the pinned public source (asserted in prepare).
            if with_fill:
                result['fill_omissions']['alleles_differ'] += int((mask & ~fillable).sum())*cohort_size
            delta = np.where(mask[:,None],-contributions,0)
            if with_fill:
                delta += np.where(filled,weights*fill_count,0)[:,None]
            result['delta'] += delta.sum(axis=0)
            result['region'] += delta[region].sum(axis=0)
            result['changed'] += np.count_nonzero(delta,axis=0)
            for j in range(cohort_size):
                top = result['top'][j]
                for i in top_candidates(delta[:,j],ordinal):
                    before = int(contributions[i,j])
                    d = int(delta[i,j])
                    top.append({'ordinal':int(ordinal[i]),'contig':f'chr{chrom}','pos1':int(pos[i]),'source_pos1':int(pos[i]),
                        'effect_allele':chr(effects[i]),'other_allele':chr(others[i]),
                        'reason':'masked_frequency_fill' if filled[i] else 'masked','baseline':before,'scenario':before+d,'delta':d})
                top.sort(key=lambda r:(-abs(r['delta']),r['ordinal']))
                del top[10:]
        print(f'oracle cohort{cohort_size}: chr{chrom}, {len(ordinal):,} eligible terms',flush=True)
    assert sum(excluded.values())+eligible_count==SCORE_ROWS
    assert int(present.sum())==eligible_count
    assert universe.dtype.itemsize==13  # unpadded u64 ordinal + u8 chromosome + u32 coordinate
    digest = hashlib.sha256(universe[present].tobytes()).hexdigest()
    for result in results.values():
        result['fill_omissions']={k:v for k,v in result['fill_omissions'].items() if v}
    return {'baseline':baseline,'references':refs,'results':results,'baseline_terms':eligible_count,
        'baseline_weight_fraction':weight_sum/manifest['absolute_weight_sum_scaled'],
        'excluded_terms':{k:v for k,v in excluded.items() if v},'baseline_term_sha256':digest,
        'chromosome_baseline_terms':chromosome_terms,'sample_ids':selection['targets'][:cohort_size]}


def percentile(values, reference):
    reference = sorted(map(int,reference))
    assert len(reference)>1 and reference[0]!=reference[-1], 'this benchmark requires a varying reference'
    return np.array([100*(bisect.bisect_left(reference,int(v))+bisect.bisect_right(reference,int(v)))/(2*len(reference)) for v in values])


def assert_close(actual, expected, name):
    assert actual is not None and math.isclose(float(actual),float(expected),rel_tol=1e-10,abs_tol=1e-10), (name,actual,expected)


def exact_number(value):
    result = Fraction(value)*(10**12*REFERENCE_SIZE)
    assert result.denominator==1, ('unexpected precision',value)
    return result.numerator


def verify(prepared, out, cohort_size):
    verifier_sha256 = sha(Path(__file__))
    manifest = verify_prepared(prepared)
    identity = json.loads((out/'identity.json').read_text())
    assert identity['prepared_manifest_sha256']==sha(prepared/'preparation.json')
    for name,digest in identity['generated_inputs_sha256'].items():
        assert sha(out/name)==digest
    started = time.monotonic()
    truth = source_arithmetic(prepared,cohort_size)
    directory = out/f'cohort{cohort_size}'/'robustness'
    report = json.loads((directory/'report.json').read_text())
    assert report['executable_sha256'].removeprefix('sha256:')==identity['pgsum_sha256']
    assert report['full_term_audit'] is False
    assert report['worker_threads']==4
    for name,digest in report['outputs'].items():
        assert sha(directory/name)==digest.removeprefix('sha256:')
    for name,digest in report['input_sha256'].items():
        assert sha(name)==digest.removeprefix('sha256:')
    score = report['scores'][0]
    assert len(report['scores'])==1 and score['pgs_id']==SCORE and score['total_model_terms']==SCORE_ROWS
    assert score['inventory_consistent'] is False, 'source declares 1,745,180 but contains 1,745,179 rows'
    assert [s['id'] for s in score['scenarios']]==list(truth['results'])
    assert [s['id'] for s in report['experiment']['scenarios']]==list(truth['results'])
    assert report['experiment']['cutoffs']==[90,95]
    for name in ['baseline_terms','excluded_terms','baseline_term_sha256']:
        actual = score[name].removeprefix('sha256:') if isinstance(score[name],str) else score[name]
        assert actual==truth[name], (name,actual,truth[name])
    assert_close(score['baseline_term_fraction'],truth['baseline_terms']/SCORE_ROWS,'baseline term fraction')
    assert_close(score['baseline_weight_fraction'],truth['baseline_weight_fraction'],'baseline weight fraction')
    baseline,refs = truth['baseline'],truth['references']
    base_pct,reference_pct = percentile(baseline,baseline),percentile(baseline,refs)
    rows = read_tsv(directory/'samples.tsv')
    top_rows = [json.loads(line) for line in (directory/'top-changes.jsonl').read_text().splitlines()]
    region_rows = read_tsv(directory/'regions.tsv')
    reference_rows = read_tsv(directory/'references.tsv')
    expected_keys = {(s,n) for s in truth['results'] for n in truth['sample_ids']}
    for table in [rows,top_rows,region_rows,reference_rows]:
        assert len(table)==len(expected_keys)
        assert {(r['scenario'],r['sample']) for r in table}==expected_keys
        assert all(r['pgs_id']==SCORE for r in table)
    tops = {(r['scenario'],r['sample']):r for r in top_rows}
    regions = {(r['scenario'],r['sample']):r for r in region_rows}
    references = {(r['scenario'],r['sample']):r for r in reference_rows}
    samples = {(r['scenario'],r['sample']):r for r in rows}
    for scenario,result in truth['results'].items():
        changed = baseline+result['delta']
        fixed,rerank = percentile(changed,baseline),percentile(changed,changed)
        refpct = percentile(changed,refs)
        for j,sample in enumerate(truth['sample_ids']):
            key = (scenario,sample)
            row,top,region,ref = samples[key],tops[key],regions[key],references[key]
            top_sum = sum(r['delta'] for r in result['top'][j])
            values = {'baseline_score':int(baseline[j]),'scenario_score':int(changed[j]),'delta':int(result['delta'][j]),
                      'top_terms_delta':top_sum,'other_terms_delta':int(result['delta'][j])-top_sum}
            for field,expected in values.items():
                assert exact_number(row[field])==expected,(key,field,row[field],expected)
            counts = {'baseline_terms':truth['baseline_terms'],'retained_terms':truth['baseline_terms']-result['masked'],
                'masked_terms':result['masked'],'unavailable_alternate_terms':0,'filled_terms':result['filled'],
                'unfilled_terms':result['masked']-result['filled'],'changed_terms':int(result['changed'][j])}
            for field,expected in counts.items():
                assert int(row[field])==expected,(key,field,row[field],expected)
            for field,expected in [('baseline_percentile',base_pct[j]),('fixed_baseline_percentile',fixed[j]),('reranked_percentile',rerank[j])]:
                assert_close(row[field],expected,(key,field))
            assert top['changed_terms']==counts['changed_terms']
            for field in ['delta','top_terms_delta','other_terms_delta']:
                assert exact_number(top[field])==values[field]
            assert len(top['top'])==len(result['top'][j])
            for actual,expected in zip(top['top'],result['top'][j]):
                for field,value in expected.items():
                    observed = exact_number(actual[field]) if field in ('baseline','scenario','delta') else actual[field]
                    assert observed==value,(key,'top',field,observed,value)
            assert region['region']=='9p21' and exact_number(region['delta'])==int(result['region'][j])
            assert ref['reference']=='fixed_reference64' and int(ref['reference_samples'])==REFERENCE_SIZE and int(ref['shared_sample_ids'])==0
            for field,value in [('baseline_percentile',reference_pct[j]),('scenario_percentile',refpct[j]),('percentile_delta',refpct[j]-reference_pct[j])]:
                assert_close(ref[field],value,(key,field))
        summary = next(r for r in score['scenarios'] if r['id']==scenario)
        for field,value in [('masked_observations',result['masked']*cohort_size),('unavailable_alternate_observations',0),
                            ('filled_observations',result['filled']*cohort_size),('unfilled_observations',(result['masked']-result['filled'])*cohort_size)]:
            assert summary[field]==value,(scenario,field,summary[field],value)
        assert summary['fill_omissions']==result['fill_omissions'],(scenario,summary['fill_omissions'],result['fill_omissions'])
        delta = result['delta'].astype(float)/(10**12*REFERENCE_SIZE)
        metrics = summary['metrics']
        assert metrics['n']==cohort_size and metrics['unavailable_reasons']=={}
        expected_metrics = {'mae':float(np.mean(np.abs(delta))),'rmse':float(np.sqrt(np.mean(delta**2))),
            'max_abs_change':float(np.max(np.abs(delta))),'pearson':float(np.corrcoef(baseline-baseline[0],changed-changed[0])[0,1]),
            'spearman':float(np.corrcoef(base_pct,rerank)[0,1]),'mean_abs_percentile_shift':float(np.mean(np.abs(fixed-base_pct))),
            'max_abs_percentile_shift':float(np.max(np.abs(fixed-base_pct)))}
        for field,value in expected_metrics.items():
            assert_close(metrics[field],value,(scenario,field))
        assert [r['cutoff'] for r in metrics['crossings']]==[90,95]
        for crossing in metrics['crossings']:
            before,after = base_pct>=crossing['cutoff'],fixed>=crossing['cutoff']
            assert crossing['up']==int((~before & after).sum()) and crossing['down']==int((before & ~after).sum())
            assert crossing['unchanged']==int((before==after).sum()) and crossing['evaluable']==cohort_size
        refsum = summary['reference_comparisons']
        assert len(refsum)==1 and refsum[0]['id']=='fixed_reference64' and refsum[0]['reference_samples']==64 and refsum[0]['shared_sample_ids']==0 and refsum[0]['evaluable']==cohort_size
        for field,value in [('mean_baseline_percentile',reference_pct.mean()),('mean_scenario_percentile',refpct.mean()),
                            ('mean_absolute_shift',np.abs(refpct-reference_pct).mean()),('max_absolute_shift',np.abs(refpct-reference_pct).max())]:
            assert_close(refsum[0][field],value,(scenario,field))
    evidence = {'cohort_size':cohort_size,'reference_samples':64,'exact_sample_scenario_checks':len(rows),
        'verifier_sha256':verifier_sha256,
        'baseline_terms':truth['baseline_terms'],'total_model_terms':SCORE_ROWS,'baseline_term_fraction':score['baseline_term_fraction'],
        'inventory_consistent':score['inventory_consistent'],
        'baseline_weight_fraction':score['baseline_weight_fraction'],'excluded_terms':truth['excluded_terms'],
        'baseline_term_sha256':truth['baseline_term_sha256'],'chromosome_baseline_terms':truth['chromosome_baseline_terms'],
        'native_measurement':json.loads((out/f'cohort{cohort_size}'/'measurement.json').read_text()),
        'source_oracle_seconds':round(time.monotonic()-started,3),'scenario_summaries':score['scenarios'],
        'source_reference_score_sha256':hashlib.sha256(refs.astype('<i8').tobytes()).hexdigest(),
        'report_sha256':sha(directory/'report.json'),'report_path':str(directory/'report.html'),
        'max_exact_raw_score_error':0,'max_exact_delta_error':0,'all_top_terms_and_regional_deltas_exact':True,
        'all_fixed_and_reranked_percentiles_verified':True,'masked_counts_include_zero_delta':True}
    assert sha(Path(__file__))==verifier_sha256, 'verifier changed during validation'
    save(out/f'cohort{cohort_size}'/'verification.json',evidence)
    print(f'PASS cohort{cohort_size}: {len(rows)} source-verified sample/scenario rows, {truth["baseline_terms"]:,} common-callable terms',flush=True)


def record(prepared, out, result):
    manifest = verify_prepared(prepared)
    identity = json.loads((out/'identity.json').read_text())
    cohorts = [json.loads((out/f'cohort{n}'/'verification.json').read_text()) for n in COHORT_SIZES]
    assert all(c['verifier_sha256']==sha(Path(__file__)) for c in cohorts)
    assert all(b['baseline_terms']<=a['baseline_terms'] for a,b in zip(cohorts,cohorts[1:]))
    summary = {'schema':'pgsum-genomewide-robustness-validation-v1','date_utc':time.strftime('%Y-%m-%d',time.gmtime()),
        'score':SCORE,'actual_published_rows':SCORE_ROWS,'file_declared_rows':manifest['file_declared_rows'],
        'score_selection':'Full unmodified PGS000018 harmonized GRCh38 file; all rows compiled and retained in model inventory.',
        'scope':'Numerical omission robustness conditional on each cohort-specific common-callable baseline. No phenotype, predictive-accuracy, population calibration, or complete published-score claim.',
        'oracle':'Independent source-allele orientation against GRCh38 FASTA and a sorted overlap sweep of original public diploid GTs. NumPy signed int64 arithmetic at exact 64e12 units, after proving absolute-sum overflow bound from all original rational weight strings. Native packs, genotype caches, and contribution outputs are not oracle inputs.',
        'checks':['all original model rows retained','all22 autosomes queried','each cohort intersection recalculated independently',
                  'every baseline/scenario raw sum and delta exact','every top10 term identity, sign and residual exact',
                  'every regional delta exact','mask/retained/fill/change counts','baseline ordinal/coordinate digest',
                  'frozen separate-reference percentiles','fixed-cohort and reranked percentiles','summary correlations/error metrics/cutoff crossings'],
        'scenarios':'no mask; 20% deterministic site-level random omissions (seed42); chr9:21500001-22500000 regional omission; each omission without and with disjoint64-reference allele-frequency fill',
        'frequency_independence':'Reference64 IDs are disjoint from all nested target8/32/128 IDs; no target genotype enters AF estimation.',
        'sampling':manifest['selection'],'sources':manifest['sources'],'reference':manifest['reference'],
        'versions':manifest['versions'],'native_identity':identity,'cohorts':cohorts,
        'preparation_manifest_sha256':sha(prepared/'preparation.json'),'harness_sha256':sha(Path(__file__)),
        'compiler':subprocess.check_output(['rustc','--version'],text=True).strip(),
        'hardware':{'platform':sys.platform,'cpu':next((line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')),None),
                    'logical_cpus':os.cpu_count(),'memory_total_kib':int(Path('/proc/meminfo').read_text().splitlines()[0].split()[1])},
        'timing_note':'One measured run per cohort, same machine/executable, no cold-cache control. Source acquisition/preparation and compile/extract/frequency build are reported separately from robustness computation. No cross-tool speed claim.',
        'preparation_measurements':manifest['measurements']}
    save(result,summary)
    print(f'wrote {result}',flush=True)


def main():
    if not __debug__:
        raise SystemExit('Do not use -O: assertions are validation gates')
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command',required=True)
    command = commands.add_parser('fetch')
    command.add_argument('--out',type=Path,required=True)
    command.add_argument('--reuse',type=Path)
    command = commands.add_parser('prepare')
    for name in ('sources','reference','out'):
        command.add_argument('--'+name,type=Path,required=True)
    command = commands.add_parser('_subset')
    for name in ('source','samples','out'):
        command.add_argument('--'+name,type=Path,required=True)
    command.add_argument('--bed',type=Path)
    command = commands.add_parser('_oracle_chrom')
    for name in ('model','genotypes','reference','selection','out'):
        command.add_argument('--'+name,type=Path,required=True)
    for name in ('native_inputs','native'):
        command = commands.add_parser(name)
        for field in ('prepared','out','pgsum'):
            command.add_argument('--'+field,type=Path,required=True)
        if name=='native':
            command.add_argument('--cohort-size',type=int,choices=COHORT_SIZES,required=True)
    command = commands.add_parser('verify')
    for field in ('prepared','out'):
        command.add_argument('--'+field,type=Path,required=True)
    command.add_argument('--cohort-size',type=int,choices=COHORT_SIZES,required=True)
    command = commands.add_parser('record')
    for field in ('prepared','out','result'):
        command.add_argument('--'+field,type=Path,required=True)
    args = parser.parse_args()
    kwargs = {k:(v.resolve() if isinstance(v,Path) else v) for k,v in vars(args).items() if k != 'command'}
    globals()[args.command](**kwargs)


if __name__ == '__main__':
    main()
