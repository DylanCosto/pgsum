# Public-data scoring concordance

The chromosome 22 validation runs pgsum, the unmodified pgsc_calc v2.3.0 scoring workflow, and a separate
PLINK 2 scoring run on real public genotypes. It preserves default-policy differences and then checks
arithmetic on explicitly shared terms. It is a correctness check, not a runtime benchmark or validation
of complete genome-wide PGS values.

## Inputs and scope

The runner selects 32 samples at evenly spaced indices in the original VCF header, retaining **all
1,059,079 chromosome 22 records** for those samples. The source contains 2,548 samples. Selection is
reproducible and is not intended to represent ancestry groups or disease outcomes.

| Published score | Original rows present in harmonized file | All chromosome 22 rows retained |
|---|---:|---:|
| PGS001229 | 51,209 | 849 |
| PGS000018 | 1,745,179 | 26,529 |

Source coordinates are explicitly mapped from `hm_chr`/`hm_pos` to GRCh38. The two chromosome-restricted
models are labelled `CHR22_PGS001229` and `CHR22_PGS000018`, never as the complete published scores.
Weights and model annotations are preserved. Separate ordinal maps retain each original scoring-file
row number. The original PGS000018 header declares 1,745,180 terms while its harmonized file contains
1,745,179 rows; the preparation manifest retains that discrepancy rather than treating the subset as a
repair of the published inventory.

Sources are the [UCSC-hosted GRCh38 1000 Genomes files](https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/),
[GRCh38 chromosome FASTAs](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/), and the Catalog's
harmonized files for [PGS001229](https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS001229/ScoringFiles/Harmonized/PGS001229_hmPOS_GRCh38.txt.gz)
and [PGS000018](https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000018/ScoringFiles/Harmonized/PGS000018_hmPOS_GRCh38.txt.gz).
The runner pins and verifies every source file's SHA-256; an upstream change requires inspection before
updating those fingerprints.

## What passed, and what differed

The checked [measurement record](../bench/results/public-chr22-validation.json) records the actual
pgsc_calc source commit, package versions, binary identities, source hashes and all validation results.
The pipeline ran with `--only_score`, covering input formatting, VCF conversion, variant matching,
scoring and verified aggregation. Ancestry adjustment and the HTML report were not run. Its source was
unmodified. This local run used PLINK v2.0.0-a.7.10LM in place of the workflow environment recipe's
2.00a5.10. The scoring packages matched that recipe: pgscatalog.core 1.1.1, pgscatalog.match 0.4.0 and
pgscatalog.calc 0.3.1. Nextflow was 24.10.0; the complete actual Python dependency set is recorded.

For these 32 samples, pgsc_calc's scoring stage uses no mean imputation. pgsum explicitly accepts
GT-only calls with `--accept-missing-quality`. No other pgsum eligibility opt-ins are enabled.

- **Default-policy comparison:** 62 of 64 reported sums differ beyond the explicit tolerance; two
  agree within it. This is a recorded result, not a failed check that is hidden by filtering.
- **Shared-term comparison:** 769 terms for PGS001229 and 25,896 for PGS000018 are selected from terms
  scored by both tools in the representative sample, with identical oriented alleles and complete
  diploid calls in all 32 samples. Separate pgsum and PLINK runs agree for all 64 shared-set scores
  within the rendered-output tolerance. These labelled `COMMON_*` models are conditional arithmetic
  checks, not substitutes for the original or chromosome-restricted scores.
- **Independent arithmetic:** Python Decimal reconstructs every external sum from retained VCF GT
  calls and the actual matched weights, verified against the original published weight text. It also
  checks all 853,312 retained default pgsum contributions against source weights, FASTA-based allele
  orientation and original VCF genotypes. All 64 default pgsum sums match these reconstructions exactly;
  all external reported sums match within tolerance. All shared-set pgsum sums match the independent
  calculation exactly too.
- **Term differences:** all 853,280 sample/term contributions scored by both tools are exactly equal.
  Every default score gap is reconciled from term eligibility differences, with exact differences and
  original overlapping-record evidence saved for all 32 samples.

The explicit comparison tolerance is `0.00001 + 0.000005 * max(abs(left), abs(right))`, applied only to
reported cross-tool values. The source-weight/genotype arithmetic checks against pgsum use exact
Decimal equality. This tolerance is specific to this validation and is not a general acceptance rule.

The representative sample HG00096 illustrates the differences:

| Chromosome-restricted model | pgsum exact sum | Reconstructed pgsc_calc sum before output rounding | External minus pgsum |
|---|---:|---:|---:|
| CHR22_PGS001229 | 0.464297034716 | 0.462057342716 | −0.002239692000 |
| CHR22_PGS000018 | 0.016987548083 | 0.023783132969 | +0.006795584886 |

The largest eligibility difference is overlapping records: pgsum withholds five PGS001229 terms and
151 PGS000018 terms that pgsc_calc scores by matching a specific variant. These produce 4,992 differing
sample/term inclusion decisions across the 32 samples. An additional PGS001229 indel is withheld by
pgsum's default orientation policy. In the other direction, pgsum scores a PGS000018 reference-effect
term where a passing reference genotype establishes the dosage, while pgsc_calc has no matching target
variant for that score allele pair. These decisions remain explicit; this benchmark does not establish
which policy is biologically preferable for every representation.

`default-term-differences.tsv` includes source ordinals, contributions, exact deltas and overlapping VCF
records with the selected sample's GT. `default-policy-oracle.json` reconciles each sample/score pair.
The HG00096 exports also exercise `pgsum compare`'s per-term interface and sum-reconciliation reports.

## Reproduce

Use Python 3.12, Java 21, Nextflow 24.10.0 and the tested PLINK binary. Install Python dependencies in an
isolated environment; the pinned requirements capture the actual validation environment. The fetch
step downloads approximately 236 MB. Prepared genomes, workflow intermediates and per-sample term
exports require additional disk space. Each stage output directory must be new, except that `fetch`
can reuse files whose pinned hashes match.

```sh
python3.12 -m venv .venv-public
.venv-public/bin/pip install -r bench/public_chr22.requirements.txt
cargo build --release

git clone --branch v2.3.0 --depth 1 https://github.com/PGScatalog/pgsc_calc.git pgsc-calc-v2.3.0
# Required pipeline commit: 72ee54f2119d2fe65ac671327342dfc412bd0a3b

.venv-public/bin/python bench/public_chr22.py fetch --out public-sources
.venv-public/bin/python bench/public_chr22.py prepare --data public-sources --out public-inputs
.venv-public/bin/python bench/public_chr22.py native --prepared public-inputs \
  --pgsum target/release/pgsum --out public-native
.venv-public/bin/python bench/public_chr22.py pgsc --prepared public-inputs \
  --nextflow /path/to/nextflow --pipeline pgsc-calc-v2.3.0 \
  --plink2 /path/to/plink2 --out public-pgsc
.venv-public/bin/python bench/public_chr22.py analyze --prepared public-inputs \
  --native public-native --pgsc public-pgsc --pgsum target/release/pgsum \
  --plink2 /path/to/plink2 --out public-analysis
```

[Nextflow 24.10.0](https://github.com/nextflow-io/nextflow/releases/tag/v24.10.0) and the
[tested Linux PLINK binary](https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260929.zip)
are external dependencies. The runner uses the official pipeline's local execution configuration;
container engines are disabled for this isolated dependency environment. It verifies the pipeline's
commit and checks for tracked source changes before executing it. Its `--only_score` aggregate is
copied from the completed workflow's working directory because the report-publishing stage is skipped.

Every command and log is retained. `public-analysis/results.json` is written only after all checks
pass. Do not use Python's `-O` option: the runner rejects it because assertions are validation gates.
The executable hashes must match between scoring and analysis and remain stable during each run.

This validation covers GT-only, diploid chromosome 22 data and partial contributions from two models.
It does not establish whole-genome model coverage, gVCF reference-block concordance, DS/GP cross-tool
agreement, sex-chromosome behavior, ancestry adjustment, or biobank-scale performance. Those remain
separate validation requirements in the roadmap.
