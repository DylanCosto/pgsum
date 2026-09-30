# pgsum

Polygenic score calculation from gVCFs, with exact observed scores and optional population-frequency filling.

**Status: early development; latest tagged release v0.5.1.** The full pipeline works; genotype calls are
validated against GIAB truth sets on seven genomes, reference blocks on four of them (see [Validation](#validation)).

## What's different

Most PGS tools score a VCF and treat a site with no variant record as missing, then fill it with a mean
dosage. A gVCF says more than that: a passing reference block is a confident homozygous-reference call.
pgsum reads the gVCF directly and uses those blocks. Strict and partial scores use observed hard calls by
default, or explicitly selected DS/GP expected contributions. Optional missing-term frequency filling and
reference-panel comparisons are reported separately.

Every score gets two answers:

- **strict**: a raw score only when every term in the model has a usable call (otherwise withheld, with
  the reasons);
- **partial**: the exact sum over usable terms, labelled as partial, with the share of terms and of total
  effect it covers.

### Which number to use

The strict score is rarely available: for HG002, 189 of 6,990 Catalog scores have every term usable by
default (329 with every opt-in), because most large scores include a few sites a genome can't call. In
practice, use the **partial score when it covers at least 99% of the terms and 99% of the total weight**
(`partial.meets_coverage_guideline`, and the `meets_coverage_guideline` column of `scores.tsv`). This is a
completeness guideline, not a guarantee of percentile stability or predictive validity. For HG002
that is 1,950 scores by default and 4,120 with every opt-in. Below that, the missing terms can move the
score noticeably.

Both are raw sums on the author's scale, not percentiles. A partial score can only be compared with other
scores computed over the same terms: a reference population scored on a different set of sites is not a
valid comparison. `score --reference-panel` does that comparison for you (below).

Sums are exact decimals, so the same scoring inputs and policies give identical score text across machines
and thread counts. Outputs record input digests; execution manifests also identify the software and settings. The exactness is for reproducibility: a result like
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
| `extract --skip-structural-alleles` | sit under structural-variant records (`<DEL>`, `<INV>`, …) in panels that carry them |
| `extract --merge-split-records` | sit at multi-allelic sites a panel splits into one record per ALT |

## Compare with an existing workflow

`pgsum compare` joins pgsum, PLINK 2 and pgsc_calc reported sums by sample and score, retaining exact
decimal differences and explicit tolerance settings. It rejects ambiguous identities and averages,
reports unmatched or unavailable values, and can compare per-term evidence from two runs. A validation
gate checks term sums too, so matching totals cannot hide different contributions. See
[comparison and migration instructions](docs/compare.md) for formats, commands and policy limitations.
The [public-data comparison](docs/public-validation.md) verifies shared-term arithmetic against PLINK
and traces default-policy differences from the actual pgsc_calc workflow on chromosome 22.

## Install

Prebuilt binaries for macOS (Apple Silicon, Intel) and Linux (x86-64, ARM) are attached to each
[release](https://github.com/DylanCosto/pgsum/releases). Or build from source with Rust 1.88 or later:

```sh
cargo install --git https://github.com/DylanCosto/pgsum
```

You also need a GRCh38 reference FASTA with its `.fai` (`samtools faidx`).

## Usage

```sh
# Download, compile, cache and score selected PGS Catalog IDs in one command.
pgsum run --gvcf sample.g.vcf.gz --reference GRCh38.fa --ids PGS000001,PGS000013 --out results/

# Check the input header and packs before scoring (may download missing scores).
pgsum run --gvcf sample.g.vcf.gz --reference GRCh38.fa --ids PGS000001 --preflight --out checks/

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

# Per-term evidence rows: status, call, dosage, contribution and source records (extract with
# --term-positions; see DESIGN.md).
pgsum evidence packs/PGS000001.pgsp --genotypes sample.pgsg --reference GRCh38.fa
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

### What is missing from a score?

Single-sample results (including each sample in `batch`) include a `missingness` report. It counts
unscored terms by reason, ranks the 20 terms with the largest possible absolute contribution, and
shows their alleles, published weights, call states, review reasons and retained source records.
Up to three records per term are shown, with truncation and omission counts; their hashes identify the
full retained text. Model-excluded positions need `extract --term-positions` to retain source evidence.
An absent record is distinguished from evidence that was never extracted.

For a supported contribution model, the report evaluates every permitted coded dosage and sums the
per-term minima and maxima exactly. `completion_score_bounds` adds those missing-contribution bounds
to the unfilled partial score. It is a conservative envelope, **not a confidence interval or a risk
estimate**: linkage and shared positions may make its endpoints impossible to attain simultaneously.
Scored DS/GP expectations stay fixed. Explicit XY chrX dosage conventions are respected.

Unsupported ploidy, ambiguous model definitions and invalid weights have no bounds. If any missing
term cannot be bounded—or inventory/publication metadata is uncertain—the whole-score envelope is
withheld. The report still gives the known subset's bounds and examples of the unbounded terms.
Frequency filling and reference placement do not change this diagnostic.

With `--terms`, the TSV includes `missing_contribution_lower`, `missing_contribution_upper` and
`missing_bounds_unavailable_because`, before any optional frequency-fill columns. The report and
columns do not change the raw or partial score.

For joint cohorts, request the detailed report explicitly:

```sh
pgsum score --genotypes cohort.pgsc --pack packs/ --out cohort-results/ --missingness
```

`cohort-missingness.jsonl.zst` contains one report per score and sample, with the same bounds and
rankings. Source references identify the original VCF/BCF by its hash, one-based record ordinal
(headers excluded), and zero-based sample index. Full cohort records are not copied into reports.
New v3 cohort files retain exclusion reasons; older v1/v2 files remain readable, but missing reasons
are labeled `missing_reason_not_retained` and cannot support complete-score bounds. Re-extract from
the original input to recover those reasons. Model-excluded positions do not retain source evidence
in cohort files. Ordinary scoring does not load the diagnostic frames.

`cohort-scores.metadata.json` identifies the report and its hash, or records `missingness_report: null`
when it was not requested. Treat that metadata as authoritative if an output directory contains a
report from an earlier invocation. Detailed reporting adds a second term/sample pass and maintains
bounded rankings for each sample; it is optional because its cost grows with cohort size.

### Catalog cache and run provenance

With `run` or `batch`, `--ids` without `--pack`/`--pack-list` downloads missing Catalog scores and
compiles them against the supplied reference. With explicit packs, `--ids` remains a selection filter.
`--catalog-cache DIR` chooses the cache root; otherwise pgsum uses `$PGSUM_CACHE_DIR`,
`$XDG_CACHE_HOME/pgsum`, or `~/.cache/pgsum`. Cache entries are separated by reference identity and
compilation rules, and verified against their file hashes before reuse. New downloads require a
consistent scoring-file/Catalog term inventory. Score licences remain in pack metadata and results.

`--offline` requires verified cached packs and makes no Catalog requests. Normal runs reuse cached
scores; `--refresh-catalog` explicitly downloads current records and scoring files again. Immutable
pack generations preserve previous versions for reproducibility. These three cache flags apply only
to ID-based selection, not explicit packs. `fetch` remains a separate bulk-download command.

`run --preflight` checks sample selection, stated contig lengths against the reference, and pack
integrity, then writes `preflight.json` without extracting or scoring. Missing contig lengths are
reported as absent assembly evidence; matching lengths alone cannot prove assembly identity. Full
record validation happens during extraction. After a successful `run`, `execution-manifest.json`
records the executable hash, input/reference identities, selected packs and source hashes, policies,
settings, and output hashes. A rerun removes the previous completion marker before writing results.
Only listed outputs belong to that invocation; unrelated older files can remain in the directory.

### Separate gVCFs in a resumable batch

Create a tab-separated sample sheet (`sample_id`, `input`, and optional `sample` for a selected VCF
column). Relative input paths are resolved against the sheet's directory:

```tsv
sample_id	input	sample
participant_1	inputs/one.g.vcf.gz	VCF_SAMPLE_1
participant_2	inputs/cohort.vcf.gz	VCF_SAMPLE_2
```

```sh
pgsum batch --samplesheet samples.tsv --reference GRCh38.fa --ids PGS000001,PGS000013 \
    --out batch-results/ --threads 8 --jobs 2 --resume
```

`--jobs` bounds concurrent samples; the total `--threads` budget is divided among them. Every input
is hashed before resume decisions. Successful samples are published atomically under
`samples/<sample_id>/<signature>/`, with a `run-manifest.json`. `--resume` reuses a generation only
when inputs, reference, packs, executable, policies, and output hashes match. New sheet rows run
independently; failed samples retain logs and do not discard successful ones. The aggregate
`batch-scores.tsv` and `batch-results.json` report completed/reused/failed samples. Concurrent writers
to the same batch directory are refused. Old generations and failed staging directories are retained.

Workers share the target index, but separate input files are decoded independently. Concurrency is
bounded; a total memory-budget scheduler and large-cohort scaling validation are still in development.

### Many samples and reference percentiles

A multi-sample VCF (a joint-called cohort, or a panel such as 1000 Genomes) is read once for every sample,
and `score` then scores every sample, with the same rules and exact sums as a single-sample run:

```sh
pgsum extract --gvcf cohort.vcf.gz --all-samples --reference GRCh38.fa --pack packs/ --out cohort.pgsc
pgsum score --genotypes cohort.pgsc --pack packs/ --out cohort-results/   # cohort-scores.tsv.zst
```

A partial score only means something next to scores over the same terms. Given a panel extracted with the
same packs and its sample groups, `score` places each score among the panel's scores computed over exactly
the terms scorable in your genome and called in at least 99% of listed panel samples, and assigns your genome
the nearest group. Remaining missing panel calls are filled from group frequencies. Terms below that call
rate are excluded from both sides and count against matched coverage:

```sh
pgsum extract --gvcf 1kgp.vcf.gz --all-samples --accept-missing-quality --reference GRCh38.fa --pack packs/ --out 1kgp.pgsc
pgsum score --genotypes sample.pgsg --pack packs/ --out results/ \
    --reference-panel 1kgp.pgsc --reference-groups integrated_call_samples_v3.20130502.ALL.panel
```

Panel inputs include 1000 Genomes phase 3 lifted to GRCh38 and the NYGC 30× release (native GRCh38, called
like a modern WGS genome). Extract the 30× release with `--skip-structural-alleles --merge-split-records`, since
it carries structural variants and splits multi-allelic sites. The comparisons in `bench/README.md` used
earlier panel rules; coverage and percentile agreement need revalidation under the current exclusion rule.

Each result then has `reference` (matched terms and coverage, `excluded_low_call_rate_terms`, percentile
among all panel samples and within each group, group mean, SD and z-score) and `ancestry` (the nearest group), and `scores.tsv` gains the
percentile in the nearest group. Read a percentile only when `reference.meets_coverage_guideline` is set (the
matched terms hold at least 99% of the terms and of the weight); below that it places a subset of the score.
If no terms remain, no percentile is reported. Z-scores are `null` when the reference group has fewer than
two samples or zero score variance. Percentiles are uncalibrated: no ancestry adjustment beyond the choice
of group, and no absolute risk.

### Filling missing terms and placing within one group

The steps a report usually needs after scoring (a completeness gate, filling the few missing terms, a
percentile within one reference population, and where an extreme result comes from) run in seconds:

```sh
# Once: an allele-frequency table from population VCFs (here 1000 Genomes phase 3 on GRCh38, EUR_AF).
pgsum frequencies --vcf phase3.chr{1..22}.GRCh38.GT.crossmap.vcf.gz --field EUR_AF \
    --manifest phase3.crossmap.GRCh38.07302021.manifest.tsv --out 1kg-eur.pgsf

# Fill missing terms (w × 2 × effect-allele frequency), and list the positions a panel needs.
pgsum score --genotypes sample.pgsg --pack packs/ --fill-frequencies 1kg-eur.pgsf \
    --placement-positions positions.txt --out results/

# Cut a PLINK panel to those positions (for example the PGS Catalog's pgsc_1000G_v1), then place within EUR.
plink2 --pfile GRCh38_1000G_ALL vzs --set-all-var-ids '@:#:$r:$a' --extract range positions.txt \
    --keep phase3-samples.txt --make-bed --out panel
pgsum score --genotypes sample.pgsg --pack packs/ --fill-frequencies 1kg-eur.pgsf \
    --reference-panel panel.bed --reference-groups integrated_call_samples_v3.20130502.ALL.panel \
    --reference-group EUR --contribution-region GRCh38_MHC.bed.gz --out results/
```

Each result gains `fill` (terms filled and left out, by method and reason, exact sums, and a `filled_score`
when at least 99% of terms and of |weight| are scorable) and `placement` (percentile, Z, the group's mean
and SD, terms without a panel line, a homozygous-reference sensitivity result, and the share of the deviation
from the group mean in the top 1 Mb window and each region given). `<PGS_ID>.reference-scores.tsv` has every
panel sample's exact score. For HG002 and HG003, 28 scores (27 million terms) take about 20 s including the
fill, placement and regions. `--karyotype XX|XY --x-dosage-model` reads chrX by sex (DESIGN.md, "Sex
chromosomes"). See DESIGN.md, "Filling missing terms" and "Placing within a reference group".

JSON output uses `pgsum-score-v4`. `dosage_field` identifies GT, DS or GP, and
`expected_genotype_terms` counts scored DS/GP contributions. `imputation.observed_scores` is true when
raw/partial scores use any such expected contribution; it stays false for GT-only scoring.
`imputation.sample_fill` reports frequency-based contributions in `fill`, including when `filled_score`
is withheld, and `imputation.reference_panel` reports expected contributions used in a panel comparison.
`imputation_performed` is the OR of these three flags, including zero-valued contributions.
Missing-term frequency fills remain separate from `raw_score` and `partial.raw_score`.

For single-sample imputed VCFs (or a selected sample), explicitly choose the input field:

```sh
pgsum run --gvcf imputed.vcf.gz --reference GRCh38.fa --pack packs/ \
  --dosage-field ds --accept-missing-quality --out results/
```

`--dosage-field gt` is the default and preserves hard-call behavior. `ds` preserves decimal dosage text
exactly for additive scores. `gp` uses the expected contribution under the genotype probabilities,
including dominant, recessive and dosage-weight models. A DS value alone cannot identify those
nonadditive expectations; those terms are withheld with `genotype_probabilities_required`.
Missing or invalid selected fields never fall back to GT. Passing gVCF reference blocks still supply
confident reference calls. FILTER, FT and any reported quality values remain enforced.

The DS/GP path accepts biallelic diploid records in single samples and cohorts. GP values must be finite, nonnegative, at most
one, and sum exactly to one; no silent renormalization occurs. Multiallelic mapping and hemizygous dosage conventions are under development.
For cohorts, use `extract --all-samples --dosage-field ds` (or `gp`) followed by `score`.
Cohort TSVs include the selected field and each sample’s expected-genotype term count;
`cohort-scores.metadata.json` records input, reference, pack and policy provenance.
Cohort scoring releases decoded genotype blocks instead of preloading the whole file. GT and
unordered quantitative scores use an evicting cache. Quantitative scores ordered by cohort block
share a pass, decoding small groups of blocks in parallel and releasing them after use. A score's
term order and arithmetic are preserved in both paths.

`score --cohort-cache-mib 64` changes the decoded-payload target (default 128 MiB).
`cohort-scores.metadata.json` records cache activity, peak decoded-wave payload and which scores shared
a pass. This is **not a total RAM limit**: packs, active handles, decoder buffers, diagnostic frames,
score accumulators and mapped files use additional memory. One block can exceed the target; a parallel
wave can also exceed it if later blocks are larger than the observed sizes used to plan that wave.
Smaller caches can increase redecoding for scores whose targets jump between blocks. Quantitative
packs are opened in groups of two to eight, according to the thread count, rather than all at once.

Quantitative cohorts are not yet supported as reference panels. See
[the implementation ledger](docs/implementation-roadmap.md) for the full expansion scope.

### Inputs

pgsum reads gVCFs from DeepVariant (long- and short-read), DRAGEN and GATK HaplotypeCaller, and plain VCFs,
on GRCh38:

- bgzipped (`.vcf.gz`), plain gzip or uncompressed VCF; native streaming BCF (compressed or uncompressed),
  detected by file contents rather than extension;
- contigs named `chr1` or `1` (and `chrM` or `MT`), in the VCF and in the reference FASTA alike;
- one sample, or several with `--sample NAME`;
- a `.tbi` or `.csi` index next to a bgzipped file lets `extract` read only the parts the scores need, so a
  run with a few scores takes a fraction of a second (`tabix -p vcf sample.g.vcf.gz` creates one). pgsum
  decides per run (`--scan auto`, the default) and the table is identical either way.

BCF uses the same genotype rules and retains the original input hash without writing a converted file.
BCF currently scans the full input; `--scan indexed` is unsupported for BCF. The decoder rejects
unsupported uncommon field encodings. DS/GP values have the precision stored in BCF's float fields;
original text precision cannot be recovered.

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

For repeatable local performance checks with synthetic gVCFs and cohort VCFs, see
[the performance runner](bench/README.md#reproducible-performance-checks). Upcoming changes and the score
JSON migration are described in [release notes](RELEASE_NOTES.md).

## Validation

What has been checked, and against what:

- **Public, hand-calculated reference checks.** `cargo test --locked --test reference_validation` creates
  tiny synthetic inputs and checks the 99% panel call-rate boundary, missingness sensitivity, withheld
  percentiles, zero-variance Z-scores, and filling metadata. No external data or private implementation is
  required. See `bench/README.md`, "Reproducible reference checks". These regression checks do not establish
  calibration on real genomes.

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
- **Against the reference implementation of the same rules.** pgsum was written to reproduce an existing
  Python implementation of these rules by the same author, which is not public. On HG002, 133 scores (91.5
  million terms: 8 development scores, 27 report scores and 100 random Catalog scores) give identical compile
  output, identical calls, dosages and contributions for every term, and identical exact sums; saved
  production results for 30 installed scores are identical too. This shows the two
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

pgsum does not ship any PGS weights. Retrieve scoring files from the
[PGS Catalog](https://www.pgscatalog.org/) with `run --ids`, `batch --ids`, `fetch`, or manually. Individual scores carry their own licence terms, and
pgsum keeps each score's licence field in its output.

## Design

See [DESIGN.md](DESIGN.md) for the pipeline, genotype rules, completeness rules and output format.

## Licence

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or [MIT license](LICENSE-MIT) at
your option.
