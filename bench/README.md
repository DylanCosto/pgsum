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
