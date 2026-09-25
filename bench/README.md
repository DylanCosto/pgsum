# Benchmarks

## plink2 `--score` (2026-09-25)

Same machine (12-core Mac, 26 GB RAM), same sample (GIAB HG002, DeepVariant gVCF), same scores, same
genotypes. plink2 cannot read gVCF reference blocks, so it is given pgsum's calls as a VCF
(`plink2-export`), i.e. plink2 gets pgsum's gVCF handling for free. Both score exactly the same terms (no
review reasons, resolved orientation, additive); plink2 runs with `no-mean-imputation cols=+scoresums`, and
its SCORE_SUM (ALT-effect table + REF-effect table) matches pgsum's partial raw score for every score to
plink2's printed precision (max relative difference 6.4e-6).

plink2 v2.0.0-a.6.39 M1 (19 Sep 2026, macOS arm64), `--threads 12`. HG002 needs `--psam` (sex 2, so chrX
calls stay diploid as pgsum scores them) and `--split-par hg38`.

| Scores | Terms | Sites | plink2: import + score | pgsum `score` |
|---|---|---|---|---|
| 8 development | 11.7M | 6.8M | 2.1 s + 3.3 s = 5.4 s | 3.7 s |
| 100 random | 56.1M | 7.2M | 3.9 s + 7.2 s = 11.1 s | 5.5 s |

pgsum times include loading a genotype table for the whole Catalog (39.3M targets, about 2.9 s).

Not included on the plink2 side: building its dense weight tables (`plink2-export` took 13 s and 30 s; in
pgsc_calc this is the combine and match steps), and turning a gVCF into genotypes at the score sites (pgsum
`extract`: 17 s for the whole Catalog with a target index). At the whole Catalog (39M sites x 6,990 scores)
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
