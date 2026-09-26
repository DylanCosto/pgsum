# pgsum

Polygenic score calculation from gVCFs, exact and without imputation.

**Status: v0.1, early development.** The full pipeline works and has been validated on one genome, GIAB
HG002 (see [Validation](#validation)); broader validation is under way.

## What's different

Most PGS tools score a VCF and treat a site with no variant record as missing, then fill it with a mean
dosage. A gVCF says more than that: a passing reference block is a confident homozygous-reference call.
pgsum reads the gVCF directly and uses those blocks, and it never imputes.

Every score gets two answers:

- **strict**: a raw score only when every term in the model has a usable call (otherwise withheld, with
  the reasons);
- **partial**: the exact sum over usable terms, labelled as partial, with the share of terms and of total
  effect it covers.

Weights are summed exactly (no floating point), and every output records the digests of its inputs.

By default pgsum follows a conservative reference implementation exactly. Opt-ins widen what can be scored,
each validated on real data and labelled in every result that uses it (see DESIGN.md):

| Option | Scores terms that… |
|---|---|
| `--allow-inferred-other-allele` | publish only an effect allele (about 30% of Catalog terms) |
| `--allow-inferred-palindromes` | are A/T or C/G SNVs, in scores whose other SNVs are on the forward strand |
| `--allow-inferred-indels` | are indels or multi-base variants (needs `compile --public-variants`) |
| `--accept-informational-descriptions` | carry a `variant_description` that is only an annotation |
| `extract --haploid-xy-as-homozygous` | sit on haploid chrX/chrY calls from callers that write them |

## Usage

```sh
# Once per PGS Catalog release: compile harmonized GRCh38 scoring files (each next to its
# <PGS_ID>.metadata.json from the Catalog REST API) into packs.
pgsum compile PGS000001_hmPOS_GRCh38.txt.gz --reference GRCh38.fa --out packs/

# Your own score: a TSV with #pgs_id= and #genome_build=GRCh38 (see DESIGN.md, "Custom scores").
pgsum compile my_score.tsv --reference GRCh38.fa --out packs/

# Per sample: read genotypes at every pack site and score, in one step ...
pgsum run --gvcf sample.g.vcf.gz --reference GRCh38.fa --pack packs/PGS000001.pgsp --out results/

# ... or in two, reusing the genotype table for more packs later.
pgsum extract --gvcf sample.g.vcf.gz --reference GRCh38.fa --pack packs/PGS000001.pgsp --out sample.pgsg
pgsum score --genotypes sample.pgsg --pack packs/PGS000001.pgsp --out results/ --terms

# Look inside packs and genotype tables.
pgsum inspect packs/PGS000001.pgsp [--header] [--genotypes sample.pgsg]
```

Many scores at once: `--pack` takes files or directories and can be repeated, `--pack-list` reads paths from
a file, and `--ids PGS000001,PGS000013` keeps only those scores. `extract` reads the gVCF once for all
selected packs. `score` runs packs in parallel and splits large packs into chunks across cores, so a
13-million-term score uses the whole machine. `--threads` sets the worker count for any command.

```sh
pgsum run --gvcf sample.g.vcf.gz --reference GRCh38.fa --pack packs/ --ids PGS000013,PGS000018 --out results/
```

`results/` gets `<PGS_ID>.score.json` per score, `scores.tsv` across scores, and with `--terms` a per-term
TSV. With `--bundle`, every result goes into one `results.jsonl.zst` (one JSON object per line) instead of a
file per score: for the whole Catalog that is 4 MB rather than 6,990 small files, which on an exFAT drive with
1 MB allocation blocks take 14 GB.

### Inputs

pgsum reads gVCFs from DeepVariant (long- and short-read), DRAGEN and GATK HaplotypeCaller, and plain VCFs,
on GRCh38:

- bgzipped (`.vcf.gz`), plain gzip or uncompressed; BCF is not read (convert with `bcftools view -Oz`);
- contigs named `chr1` or `1` (and `chrM` or `MT`), in the VCF and in the reference FASTA alike;
- one sample, or several with `--sample NAME`;
- a `.tbi` or `.csi` index next to a bgzipped file lets `extract` read only the parts the scores need, so a
  run with a few scores takes a fraction of a second (`tabix -p vcf sample.g.vcf.gz` creates one). pgsum
  decides per run (`--scan auto`, the default) and the table is identical either way.

A VCF whose `##contig` lengths differ from the reference (for example GRCh37) is refused. A plain VCF has no
reference blocks, so sites without a record are `unknown_no_record` rather than assumed reference. DRAGEN
writes male chrX and chrY as haploid; pass `--haploid-xy-as-homozygous` to `extract` or `run`.

pgsum records the SHA-256 of every input. Hashing a 3 GB reference takes a second, so digests are remembered
in `~/.cache/pgsum/digests.tsv` (or `$XDG_CACHE_HOME/pgsum`, or `$PGSUM_CACHE_DIR`), keyed by path, size,
modification time and inode; `PGSUM_NO_DIGEST_CACHE=1` turns this off.

The whole Catalog: `pgsum fetch --all --reference GRCh38.fa --out packs/` downloads and compiles all 6,990
scores (4.5 billion terms, 31.6 GB of packs); an interrupted run resumes. For many samples against the same
packs, pass `--targets-cache packs.pgst` to `extract` or `run`: the first run records which sites the packs
need, later runs skip reading every pack to find out.

### Speed

On a 12-core Mac with the packs on a USB SSD and a 36-million-record HG002 gVCF:

| Run | Time |
|---|---|
| `run`, one score (77 terms), indexed gVCF | 0.4 s |
| `run`, four small scores (588 targets), indexed gVCF | 0.8 s |
| `extract`, all 6,990 scores (43.3M targets, with the target index) | 12 s |
| `score`, all 6,990 scores (4.5 billion terms) | 141 s |
| `score`, three scores against a saved whole-Catalog genotype table | 0.8 s |

Without an index, `extract` reads the whole file on all cores: 5 s for a 154-million-record GATK gVCF.

## Validation

What has been checked, and against what:

- **Genotype calls against a truth set.** At every score site on chr1–22 inside the GIAB v4.2.1 HG002
  benchmark regions (40.6 million sites; `bench/giab_concordance.py`):

  | HG002 input | Sites called | SNV agreement | Indel agreement | Wrong hom-ref calls |
  |---|---|---|---|---|
  | DeepVariant, PacBio Revio | 98.9% | 99.9989% | 99.958% | 390 of 34.0M |
  | DeepVariant, Illumina | 99.0% | 99.9966% | 99.986% | 1,050 of 34.1M |
  | DRAGEN 3.7.6, Illumina | 98.7% | 99.9974% | 99.994% | 498 of 34.0M |

  "Wrong hom-ref calls" are sites read as homozygous reference (mostly from gVCF reference blocks) where
  GIAB has a variant. For the long-read gVCF, 886 of the 887 disagreements are the genotype the caller wrote
  in the gVCF; the other is a rule noted in DESIGN.md (Open questions).
- **Arithmetic against plink2.** On the same HG002 genotypes, pgsum's sums match plink2 `--score` for 108
  scores to plink2's printed precision (`bench/README.md`). This does not test genotype calling: plink2 was
  given pgsum's calls.
- **Against pgsc_calc.** Where both score the same variants the sums agree (7 of 8 development scores; the
  rest differ in which variants each includes, by design).
- **Against the reference implementation of the same rules.** pgsum was written to reproduce the PGS scorer
  in an existing Python implementation by the same author, which is not public. Compile output is identical
  for 52 scores, and every term's call, dosage and contribution, and every sum, for 8 scores (11.7 million
  terms) on HG002. This is a check that the two implementations agree, not an independent check that either
  is right.
- **Across technologies.** HG002 scores complete from both the long-read and an Illumina gVCF are identical
  for 258 of 259 scores (213 of 215 against DRAGEN); each difference is one genotype call.

## Limitations

- GRCh38 only; VCFs on another assembly are refused.
- One sample per run. Many-sample performance has not been benchmarked.
- Raw scores only: no ancestry adjustment, percentiles or absolute risk.
- A plain VCF (no reference blocks) scores poorly: pgsum does not assume the reference where a VCF is silent.
- Validated end to end on one genome (HG002); calls from other callers and samples follow the same rules
  but have not been checked against a truth set beyond it.

## Scoring files and licences

pgsum does not ship any PGS weights. Download scoring files from the
[PGS Catalog](https://www.pgscatalog.org/) yourself. Individual scores carry their own licence terms, and
pgsum keeps each score's licence field in its output.

## Design

See [DESIGN.md](DESIGN.md) for the pipeline, genotype rules, completeness rules and output format.

## Licence

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or [MIT license](LICENSE-MIT) at
your option.
