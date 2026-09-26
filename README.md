# pgsum

Polygenic score calculation from gVCFs, exact and without imputation.

**Status: v0.3, early development.** The full pipeline works; genotype calls are validated against the seven
GIAB truth sets (see [Validation](#validation)).

## What's different

Most PGS tools score a VCF and treat a site with no variant record as missing, then fill it with a mean
dosage. A gVCF says more than that: a passing reference block is a confident homozygous-reference call.
pgsum reads the gVCF directly and uses those blocks, and it never imputes.

Every score gets two answers:

- **strict**: a raw score only when every term in the model has a usable call (otherwise withheld, with
  the reasons);
- **partial**: the exact sum over usable terms, labelled as partial, with the share of terms and of total
  effect it covers.

### Which number to use

The strict score is rarely available: for HG002, 189 of 6,990 Catalog scores have every term usable by
default (329 with every opt-in), because most large scores include a few sites a genome can't call. In
practice, use the **partial score when it covers at least 99% of the terms and 99% of the total weight**
(`partial.meets_coverage_guideline`, and the `meets_coverage_guideline` column of `scores.tsv`). For HG002
that is 1,950 scores by default and 4,120 with every opt-in. Below that, the missing terms can move the
score noticeably.

Both are raw sums on the author's scale, not percentiles. A partial score can only be compared with other
scores computed over the same terms: a reference population scored on a different set of sites is not a
valid comparison. `score --reference-panel` does that comparison for you (below).

Sums are exact decimals, so the same inputs give byte-identical output on any machine and thread count,
and every output records the digests of its inputs. The exactness is for reproducibility: a result like
`-2461760.133929721328533131` is exact, but the weights behind it carry only a few significant figures.

By default pgsum follows a conservative reference implementation exactly. Opt-ins widen what can be scored,
each validated on real data and labelled in every result that uses it (see DESIGN.md):

| Option | Scores terms that… |
|---|---|
| `--allow-inferred-other-allele` | publish only an effect allele (about 30% of Catalog terms) |
| `--allow-inferred-palindromes` | are A/T or C/G SNVs, in scores whose other SNVs are on the forward strand |
| `--allow-inferred-indels` | are indels or multi-base variants (needs `compile --public-variants`) |
| `--accept-informational-descriptions` | carry a `variant_description` that is only an annotation |
| `extract --haploid-xy-as-homozygous` | sit on haploid chrX/chrY calls from callers that write them |
| `extract --accept-missing-quality` | come from genotype-only VCFs (imputed, array, joint-called: `GT` without depth or GQ) |

## Install

Prebuilt binaries for macOS (Apple Silicon, Intel) and Linux (x86-64, ARM) are attached to each
[release](https://github.com/DylanCosto/pgsum/releases). Or build from source with Rust 1.88 or later:

```sh
cargo install --git https://github.com/DylanCosto/pgsum
```

You also need a GRCh38 reference FASTA with its `.fai` (`samtools faidx`).

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

### Many samples and reference percentiles

A multi-sample VCF (a joint-called cohort, or a panel such as 1000 Genomes) is read once for every sample,
and `score` then scores every sample, with the same rules and exact sums as a single-sample run:

```sh
pgsum extract --gvcf cohort.vcf.gz --all-samples --reference GRCh38.fa --pack packs/ --out cohort.pgsc
pgsum score --genotypes cohort.pgsc --pack packs/ --out cohort-results/   # cohort-scores.tsv.zst
```

A partial score only means something next to scores over the same terms. Given a panel extracted with the
same packs and its sample groups, `score` places each score among the panel's scores computed over exactly
the terms scorable in your genome and in every panel sample, and assigns your genome the nearest group:

```sh
pgsum extract --gvcf 1kgp.vcf.gz --all-samples --accept-missing-quality --reference GRCh38.fa --pack packs/ --out 1kgp.pgsc
pgsum score --genotypes sample.pgsg --pack packs/ --out results/ \
    --reference-panel 1kgp.pgsc --reference-groups integrated_call_samples_v3.20130502.ALL.panel
```

Each result then has `reference` (matched terms and coverage, percentile among all panel samples and within
each group, group mean, SD and z-score) and `ancestry` (the nearest group), and `scores.tsv` gains the
percentile in the nearest group. For HG002 against 1000 Genomes this takes 33 s for 100 scores and assigns
EUR. Percentiles are uncalibrated: no ancestry adjustment beyond the choice of group, and no absolute risk.

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

- **Genotype calls against truth sets, on seven genomes.** At every score site on chr1–22 inside each GIAB
  v4.2.1 benchmark region (about 40 million sites per genome; `bench/giab_concordance.py`):

  | Genome (ancestry) | Input | Sites called | SNV agreement | Indel agreement | Wrong hom-ref calls |
  |---|---|---|---|---|---|
  | HG001 (European) | DRAGEN 3.7.6 gVCF, Illumina | 98.9% | 99.9968% | 99.995% | 698 of 33.6M |
  | HG002 (Ashkenazi) | DeepVariant gVCF, PacBio Revio | 98.9% | 99.9989% | 99.958% | 390 of 34.0M |
  | HG002 | DeepVariant gVCF, Illumina | 99.0% | 99.9966% | 99.986% | 1,050 of 34.1M |
  | HG002 | DRAGEN 3.7.6 gVCF, Illumina | 98.7% | 99.9974% | 99.994% | 498 of 34.0M |
  | HG003 (Ashkenazi) | DeepVariant 1.5 gVCF, PacBio Revio | 98.8% | 99.9990% | 99.983% | 97 of 33.8M |
  | HG004 (Ashkenazi) | DRAGEN 3.7.6 gVCF, Illumina | 98.7% | 99.9974% | 99.994% | 595 of 33.7M |
  | HG005 (Han Chinese) | DRAGEN 4.2.4 VCF (no gVCF public) | 14.8%¹ | 99.9988% | 99.981% | — |
  | HG006 (Han Chinese) | DRAGEN 4.2.4 VCF | 14.7%¹ | 99.9980% | 99.980% | — |
  | HG007 (Han Chinese) | DRAGEN 4.2.4 VCF | 14.8%¹ | 99.9978% | 99.968% | — |

  ¹ A plain VCF has records only where the caller saw a variant, so only those sites are called; for these
  three genomes the check covers variant calls, not reference blocks.

  "Wrong hom-ref calls" are sites read as homozygous reference (mostly from gVCF reference blocks) where
  GIAB has a variant. For the long-read HG002 gVCF, 886 of the 887 disagreements are the genotype the caller
  wrote in the gVCF; the other is a rule measured in DESIGN.md (Open questions).
- **Arithmetic against plink2.** On the same HG002 genotypes, pgsum's sums match plink2 `--score` for 108
  scores to plink2's printed precision (`bench/README.md`). This does not test genotype calling: plink2 was
  given pgsum's calls.
- **Against pgsc_calc.** Where both score the same variants the sums agree (7 of 8 development scores; the
  rest differ in which variants each includes, by design).
- **Against the reference implementation of the same rules.** pgsum was written to reproduce the PGS scorer
  in an existing Python implementation by the same author, which is not public. On HG002, 133 scores (91.5
  million terms: 8 development scores, the 27 report models and 100 random Catalog scores) give
  identical compile output, identical calls, dosages and contributions for every term, and identical exact
  sums; the case's saved production evidence for its 30 installed models is identical too. This shows the two
  implementations agree, not that either is right: that is what the GIAB check above is for.
- **Across technologies.** HG002 scores complete from both the long-read and an Illumina gVCF are identical
  for 258 of 259 scores (213 of 215 against DRAGEN); each difference is one genotype call.

## Limitations

- GRCh38 only; VCFs on another assembly are refused.
- Percentiles are relative to a reference panel you supply (for example 1000 Genomes), within the nearest
  group; there is no ancestry adjustment beyond that choice and no absolute risk. Placing the whole Catalog
  needs a panel extracted for all of it, which is large; panels are meant for a chosen set of scores.
- Cohort scoring costs grow with scores × samples: 100 scores for 2,504 samples take about 2 minutes to
  score after a 3.5-minute read, and the whole Catalog for a large cohort is not practical on one machine.
- A plain VCF (no reference blocks) scores poorly: pgsum does not assume the reference where a VCF is silent.
- chrX dosage: terms on chrX count 0, 1 or 2 copies as the gVCF's diploid calls give them, so a male's
  hemizygous ALT counts as 2. Scores differ in whether their authors coded males 0/1 or 0/2; results for
  scores with chrX terms say so (`sex_chromosomes`).
- Genotype calls are checked against truth sets for the seven GIAB genomes (above); no East Asian gVCF is
  public, so reference-block reading is checked on European and Ashkenazi genomes only.

## Scoring files and licences

pgsum does not ship any PGS weights. Download scoring files from the
[PGS Catalog](https://www.pgscatalog.org/) yourself. Individual scores carry their own licence terms, and
pgsum keeps each score's licence field in its output.

## Design

See [DESIGN.md](DESIGN.md) for the pipeline, genotype rules, completeness rules and output format.

## Licence

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or [MIT license](LICENSE-MIT) at
your option.
