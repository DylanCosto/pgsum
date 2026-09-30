# Benchmarks

Summary (12-core Mac, 26 GB RAM, GIAB HG002 DeepVariant gVCF; pgsc_calc and plink2 measured 2026-09-25,
pgsum re-measured 2026-09-26 at commit `d6999ad`):

| Per sample | pgsum | pgsc_calc v2.3.0 steps | plink2 `--score` alone |
|---|---|---|---|
| 8 scores (11.7M terms) | 4.4–6.3 s | 432 s | 5.4 s |
| 100 random scores (56.1M terms) | 6.9–7.3 s (3.4 GB peak) | 2,224 s (29.5 GB peak) | 11.1 s |
| Whole Catalog (6,990 scores, 4.50B terms) | about 2.5 min | not run (see below) | not run |

pgsum ranges are a first and a second run (the second with the gVCF in the file cache). On 2026-09-25 the
same runs took 10.1 s and 11.8 s; the difference is the parallel gVCF scan, reuse of parsed records and
cached input digests (see DESIGN.md).

pgsum's times start from the gVCF (`pgsum run`: extract + score). pgsc_calc and plink2 start from a VCF of
genotypes at the score sites that pgsum produced for them; turning a gVCF into that is not counted for them.
plink2 alone is only the scoring engine: it also needs its weight tables built (pgsc_calc's format and match
steps).


## Many samples (2026-09-26)

Same Mac, pgsum at the commit that added `--accept-missing-quality`.

**Per-sample gVCFs**, pgsum's intended input: 20 1000 Genomes samples (4 per superpopulation) called with
DeepVariant 0.8 on Illumina, about 97 million records each, from Google's public
`brain-genomics-public/research/cohort/1KGP/dv_vcf/v1` (no index, so each gVCF is read whole):

| Scores | Per sample | Peak memory | Throughput |
|---|---|---|---|
| 100 random (56.1M terms) | median 6.6 s (5.9–9.1 s) | 3.7 GB | about 545 samples/hour |
| 100 random, two samples at a time | 5.3 s per sample | 3.8 GB each | about 685 samples/hour |
| Whole Catalog (6,990 scores, 4.50B terms) | median 145 s (139–154 s) | up to 14.5 GB | about 25 samples/hour |

plink2 and pgsc_calc cannot read per-sample gVCFs; scoring these with them needs joint genotyping first,
which was not measured.

**One joint-called cohort VCF**, plink2's intended input: 1000 Genomes phase 3 lifted to GRCh38 (2,504
samples, 81.6 million records, 15.5 GB), the same 100 scores:

| | Time |
|---|---|
| plink2: import to pgen at the 7.2M score sites | 340 s |
| plink2: `--score` for all 2,504 samples (alt and ref tables) | 14 s |
| **plink2 total, 2,504 samples** | **354 s (0.14 s per sample)** |
| pgsum `extract --all-samples --accept-missing-quality`, all 2,504 samples in one pass | 209 s (0.8 GB peak) |
| pgsum `score` on the cohort file, 100 scores × 2,504 samples (exact sums) | 115 s |
| **pgsum total, 2,504 samples** | **324 s (0.13 s per sample)** |
| pgsum `run --sample S`, one sample at a time (before the cohort mode) | 103–113 s per sample |

The cohort mode reads the file once for every sample. Every one of the 200 scores of HG00096 and NA12878 is
identical, as text, to their single-sample runs (`tests/cohort_fixture.rs` checks the same on a fixture).
Scoring is slower than plink2's because every sum is exact; the read is faster than plink2's import.

Placing HG002 among these 2,504 samples for the same 100 scores (`score --reference-panel`), over the terms
scorable in both, takes 33 s and assigns HG002 to EUR.

The file has only `GT`, so by default pgsum calls nothing (`quality_missing`); `--accept-missing-quality`
accepts genotype-only calls. With it, pgsum's sums for HG00096 match plink2's to plink2's precision on the
sites plink2 scores; the tools differ in which sites they use: pgsum does not score the 1.5% of target sites
where the file has several records at one position (split multi-allelic sites), and it scores multi-allelic
sites written as one record, which plink2 cannot match to `chr:pos:ref:alt` IDs.

## Reference panels (2026-09-26)

**Historical results:** these measurements predate the `pgsum-score-v3` correction. Current `.pgsc` panel
comparisons exclude terms called in fewer than 99% of listed panel samples from both sides, count them
against coverage, and report undefined Z-scores as null. The coverage and percentile comparisons below
have not been rerun under this rule; they are not validation of the current behavior.

HG002 placed among the 2,504 unrelated 1000 Genomes samples for the same 100 scores, with two versions of
the panel: phase 3 lifted to GRCh38 (low coverage, imputed) and the NYGC 30× release (native GRCh38,
`20220422_3202_phased_SNV_INDEL_SV`). A term is used only if it is scorable in HG002 and in every panel
sample.

| Panel | Placed | Meet the 99% guideline | Matched weight, median (min) | Extract, 2,504/3,202 samples |
|---|---|---|---|---|
| Phase 3, lifted | 64 | 26 | 93.7% (41%) | 209 s |
| NYGC 30× | 64 | 0 | 75.4% (34%) | 322 s |
| NYGC 30×, `--skip-structural-alleles` | 64 | 2 | 92.6% (44%) | 335 s |
| Phase 3, `--merge-split-records` | 64 | 26 | 93.7% (41%) | 178 s |
| NYGC 30×, both options | 64 | 2 | 92.8% (44%) | 279 s |
| Phase 3, with panel fill (2026-09-26 rule) | 64 | 28 | 93.8% | — |
| NYGC 30×, both options, with panel fill | 64 | 26 | 93.6% | — |

The 30× release carries structural variants as symbolic records spanning many kilobases; every score site
under one has two overlapping records, which pgsum's genotype rules call ambiguous, so the site is unusable
for all samples (23% of target calls). `--skip-structural-alleles` leaves those records out and recovers most
of the matched weight, but fewer scores reach the guideline than with phase 3: the 30× release also splits
multi-allelic sites into separate records five times as often (37,292 positions on chr22 against 7,169) and
calls more indels over SNVs, and a single ambiguous panel sample at a site removes the term. For the two
scores that meet the guideline on both panels the EUR percentiles agree within 0.2 points.

`--merge-split-records` reads records split from one multi-allelic site as that site; it cuts ambiguous
calls by 73% on the 30× panel and 71% on phase 3, but not the number of scores meeting the guideline. What
limits both panels now is the rule that a term is used only when every one of the 2,504 panel samples has a
passing call there: one sample carrying a third allele, or a site absent from the panel, removes the term
(the 30× release has more of both: 639 against 500 million `other_called_allele` calls and 360 against 202
million target calls without a record).

Since then a term is used when it is called in the panel at all: panel calls missing where at least 99% of the
panel is called are filled with the group's expected contribution, and a term called in fewer gives every
panel sample that expectation (DESIGN.md, "Reference panels"). With this rule both panels place the same scores
about equally well (28 and 26 meet the guideline), filling 566,198 and 1,090,870 panel genotypes over the 100
scores, at most 388 terms for any one panel sample. For the 25 scores meeting the guideline on both, the EUR
percentiles differ by a median of 0.8 points (90th percentile 2.8, largest 10.3).

## plink2 `--score` (2026-09-25)

Same machine (12-core Mac, 26 GB RAM), same sample (GIAB HG002, DeepVariant gVCF), same scores, same
genotypes. plink2 cannot read gVCF reference blocks, so it is given pgsum's calls as a VCF
(`plink2-export`), i.e. plink2 gets pgsum's gVCF handling for free. Both score exactly the same terms (no
review reasons, resolved orientation, additive); plink2 runs with `no-mean-imputation cols=+scoresums`, and
its SCORE_SUM (ALT-effect table + REF-effect table) matches pgsum's partial raw score for every score to
plink2's printed precision (max relative difference 6.4e-6).

plink2 v2.0.0-a.6.39 M1 (19 Sep 2026, macOS arm64), `--threads 12`. HG002 needs `--psam` (sex 2, so chrX
calls stay diploid as pgsum scores them) and `--split-par hg38`.

| Scores | Terms | Sites | plink2: import + score | pgsum `score` (2026-09-25) | pgsum `score` (2026-09-26) |
|---|---|---|---|---|---|
| 8 development | 11.7M | 6.8M | 2.1 s + 3.3 s = 5.4 s | 3.7 s | 1.9 s |
| 100 random | 56.1M | 7.2M | 3.9 s + 7.2 s = 11.1 s | 5.5 s | 3.0 s |

pgsum times are against a genotype table for the whole Catalog (39.3M targets then, 43.3M now). Since
2026-09-26 tables are stored in blocks and `score` decompresses only those its packs need.

Not included on the plink2 side: building its dense weight tables (`plink2-export` took 13 s and 30 s; in
pgsc_calc this is the combine and match steps), and turning a gVCF into genotypes at the score sites (pgsum
`extract`: 12–15 s for the whole Catalog with a target index). At the whole Catalog (39M sites x 6,990 scores)
a single dense table is impractical, so plink2 has to run in batches.

To reproduce:

```sh
cargo build --release --manifest-path bench/plink2-export/Cargo.toml
bench/plink2-export/target/release/plink2-export genotypes.pgsg packs.txt in/
plink2 --vcf in/genotypes.vcf --psam sample.psam --split-par hg38 --make-pgen --out g --threads 12
plink2 --pfile g --score in/alt.scores.tsv 1 2 header-read no-mean-imputation cols=+scoresums \
  --score-col-nums 3-<N+2> --threads 12 --out alt     # and the same for in/ref.scores.tsv
pgsum score --genotypes genotypes.pgsg --pack packs/ --ids <the PGS IDs in in/columns.txt> --out results/
```

## pgsc_calc (2026-09-25)

pgsc_calc v2.3.0 (June 2026) is a Nextflow pipeline. Its heavy steps are the `pgscatalog-utils` tools and
plink2; they were run natively (arm64) with the versions pgsc_calc pins (`pgscatalog.core` 1.1.1,
`pgscatalog.match` 0.4.0, `pgscatalog.calc` 0.3.1) and plink2 a6.39, with the command lines from its
modules, and without Nextflow, so no workflow or container overhead is counted. (On this Mac, pgsc_calc's own
conda and Docker profiles would run x86 builds under emulation: bioconda has no arm64 plink2 and there is no
Java.) Target: the same genotype VCF as for plink2 (chromosomes renamed `1`, `2`, … to match scoring files),
`--min_overlap 0` so no score is dropped. Scripts: `run.sh` in the benchmark folder.

| Step | 8 scores | 100 scores |
|---|---|---|
| `pgscatalog-format` | 198.1 s | 564.9 s |
| plink2 import | 2.3 s | 2.6 s |
| `pgscatalog-match` | 16.5 s | 153.6 s |
| `pgscatalog-matchmerge` | 201.9 s | 1,448.9 s |
| plink2 `--score` | 4.0 s | 24.2 s |
| `pgscatalog-aggregate` | 9.0 s | 29.8 s |
| **Total** | **431.8 s** | **2,224 s** (peak 29.5 GB, above this Mac's RAM) |
| pgsum `run` from the gVCF (2026-09-25) | 10.1 s | 11.8 s (peak 4.6 GB) |
| pgsum `run` from the gVCF (2026-09-26) | 4.4–6.3 s | 6.9–7.3 s (peak 3.4 GB) |

pgsum's format-and-match equivalent (`compile`) runs once per Catalog release (5.3 s for the 8 scores), not
per sample. For the whole Catalog, pgsc_calc's per-term cost here (about 10 µs to format, 26 µs to
match-merge) would be on the order of 40 hours per run on this machine, and its memory grew past RAM at 100
scores, so it would have to run in batches; this was not measured.

Agreement: where both tools score the same rows the sums agree to pgsc_calc's printed precision (7 of the 8
development scores; 55 of 99 scored random ones with `--allow-inferred-other-allele`). The others differ in
which rows are included, not in arithmetic: pgsum's review rules exclude rows pgsc_calc scores
(`variant_description` in PGS000667, missing other alleles pgsum cannot infer, ambiguous overlapping gVCF
records), and pgsc_calc only sees sites present in the VCF it was given, which for effect-allele-only scores
misses sites pgsum infers (e.g. PGS000054, whose effect allele is the reference at every site).

## Reproducible reference checks

From a checkout with Rust 1.88 or later:

```sh
cargo test --locked --test reference_validation
```

This public validation suite needs no downloads beyond the pinned Cargo dependencies, no private Python
implementation, and no real genome. Its complete synthetic input generator and expected answers are
versioned in `tests/reference_validation.rs`. It exercises the compiled CLI, including JSON and bundled
output. Temporary inputs and outputs are removed after each test.

The reference sequence is `ACGT`. Two additive terms have weight 1: chr1:1 A>G and chr1:4 T>C. The sample
has dosages 1 and 2. There are 100 panel samples, all in one group. At the first term, 34 have dosage 0,
33 have dosage 1 and 33 have dosage 2. At the second term every called sample has dosage 1. The suite
varies the number called at that second term:

| Called at term 2 | Matched terms | Sample comparison score | Mid-rank percentile | Coverage | Filled panel calls |
|---|---|---|---|---|---|
| 100/100 | 2 | 3 | 83.5 | 100% | 0 |
| 99/100 | 2 | 3 | 83.5 | 100% | 1 |
| 98/100 | 1 | 1 | 50.5 | 50% | 0 |
| 0/100 | 1 | 1 | 50.5 | 50% | 0 |

With both terms, 67 panel scores are below 3 and 33 equal it: `100 × (67 + 33/2) / 100 = 83.5`.
After excluding term 2, 34 are below 1 and 33 equal it: `100 × (34 + 33/2) / 100 = 50.5`. Changing the
sample's excluded dosage from 2 to 0 must leave the entire reference result unchanged while changing its
observed score. These answers come from counting dosages, independently of the scoring implementation.

Additional cases cover exclusion of the only term (no percentile), a constant retained panel (`z: null`),
and a zero-frequency sample fill (reported as filling even though its contribution is zero and the final
filled score is withheld). The full suite also checks PLINK zero-variance Z-scores and panel-fill metadata:

```sh
cargo test --locked
```

These are regression and arithmetic checks. They do not validate ancestry assignment, real-data
calibration, or broad percentile stability. The next real-data validation should rerun the two panel
comparisons above with pinned source manifests and per-score outputs, including missingness sensitivity.

## Reproducible performance checks

Build a release binary and run the deterministic synthetic workload generator. It needs Python 3,
`zstd`, and GNU `time` (`/usr/bin/time` on Linux; `gtime` on macOS). It downloads no genomes or weights.
Use a new output directory for each run; inputs and logs are retained with `results.json`.

```sh
cargo build --release --locked
python3 bench/performance.py --output bench/results/local \
  --terms 100000 --samples 257 --threads 1 4 --repeats 5
```

The runner measures compilation separately, then extraction and scoring for a single gVCF with reference
blocks and for sparse and dense joint VCFs. It also measures the combined gVCF `run` command. Two scores
use opposite effect alleles, signed weights, and zero weights. Cohorts include missing calls; 257 samples
exercise a partial packed byte. All sums and cohort coverage are checked against an independent integer
calculation, and full scoring outputs are compared across thread counts (excluding the version field).
The existing Rust tests separately cover indexed reads, non-additive models, very wide decimal weights,
indels, and other input policies. Synthetic throughput does not predict whole-genome throughput.

For a before/after comparison, retain the previous release binary outside `target/`, then pass it explicitly:

```sh
python3 bench/performance.py --baseline /path/to/previous/pgsum \
  --binary target/release/pgsum --output bench/results/comparison \
  --terms 100000 --samples 257 --threads 1 4 --repeats 5
```

Both binaries process the same generated inputs. Output differences fail the run; use comparable schema
versions. Each command gets one unmeasured warm-up and the requested number of timed runs. The report
records each run's wall time, CPU time, peak resident memory, binary/input hashes, environment, and output
hashes. OS caches are not flushed; compilation and extraction are not included in isolated scoring times.
Baseline runs precede candidate runs, so repeat on an idle machine if temperature/load drift matters.

On a dedicated host, `--max-score-slowdown 1.15` also fails if any scoring median exceeds 115% of its baseline.
Choose sufficiently large workloads and multiple repetitions for this gate. CI uses a small release-mode
run to gate correctness and uploads the measurements, without treating noisy shared-host timing as a gate.

### Cohort optimization measurement

Measured on Linux x86-64, Intel Core i5-12600K, using 100,000 terms per score, two scores, 257 samples,
five measured runs after warm-up. Baseline: PR #1 commit `1f41987`, before the cohort-loop optimization.
Candidate: packed-byte skipping and precomputed exact contribution differences. These are synthetic
measurements and are separate from the earlier real-genome benchmarks above.

| Scoring workload | Threads | Baseline median | Optimized median | Time reduction |
|---|---:|---:|---:|---:|
| Sparse cohort | 1 | 0.0913 s | 0.0692 s | 24% |
| Dense cohort | 1 | 0.3795 s | 0.1960 s | 48% |
| Sparse cohort | 4 | 0.0553 s | 0.0431 s | 22% |
| Dense cohort | 4 | 0.2045 s | 0.1072 s | 48% |

All cohort output bytes and single-sample score JSON matched across the two binaries and thread counts.
The independent integer oracle also passed. gVCF end-to-end times were approximately unchanged
(0.151 to 0.153 s at one thread); this optimization targets cohort scoring. Cohort extraction still takes
0.26–0.57 s here, so complete-workflow gains are smaller than scoring-only gains. Peak cohort scoring RSS
remained below 33 MiB for both binaries; no memory reduction is claimed.
