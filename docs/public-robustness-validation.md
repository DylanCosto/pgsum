# Public missingness experiment

`bench/robustness_public.py` runs the robustness command on real public genotypes and
independently reconstructs its results from the original scoring-file weights, FASTA
bases and VCF genotypes. The checked [measurement record](../bench/results/public-robustness-validation.json)
contains source fingerprints, sample IDs, binary identity, the full report and verification counts.
This tests score reconstruction after masking; no phenotype data are used.

## Models and the available-data baseline

All published rows of four coronary artery disease/heart disease models are compiled,
without chromosome restriction, weight changes or removal of inconvenient rows.
They share the Catalog trait mapping `MONDO_0005010`.

| Catalog score | Publication | Published rows retained | Common-callable baseline | Excluded: orientation/model | Excluded: original calls |
|---|---|---:|---:|---:|---:|
| [PGS000010](https://www.pgscatalog.org/score/PGS000010/) | Mega et al., Lancet (2015) | 27 | 18 | 9 | 0 |
| [PGS000011](https://www.pgscatalog.org/score/PGS000011/) | Tada et al., European Heart Journal (2015) | 50 | 43 | 6 | 1 |
| [PGS000057](https://www.pgscatalog.org/score/PGS000057/) | Natarajan et al., Circulation (2017) | 57 | 31 | 21 | 5 |
| [PGS000349](https://www.pgscatalog.org/score/PGS000349/) | Pechlivanis et al., BMC Medical Genetics (2020) | 70 | 59 | 9 | 2 |

These baselines are **conditional partial scores, not complete published PGS values**.
The command requires one frozen term set callable in every target and selected reference
sample. Ambiguous orientation, unavailable genotypes and mismatched source alleles remain
in the model inventory and exclusion totals. Each retained ordinal/coordinate is included
in the baseline fingerprint. No additional terms are removed to make a scenario look better.

PGS000010 and PGS000057 lack authored other alleles. This experiment explicitly opts into
Catalog-inferred other alleles for baseline genotype scoring. It does not enable inferred
palindromes, indels or informational descriptions. Frequency filling retains the production
eligibility rules: these two models' missing-author-allele review reasons prevent filling,
and every withheld fill is reported. PGS000011 and PGS000349 exercise successful fills.

## Public cohort, reference samples and masks

The source is the [UCSC-hosted GRCh38 1000 Genomes VCF release](https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/).
From its 2,548 VCF-header samples the runner selects 128 evenly spaced indices
`floor(i * (2548 - 1) / 127)`, for `i = 0..127`. Every fourth selected sample is a target
(32 people); the other 96 supply the frequency reference. The latter are divided into
two interleaved 48-person comparison panels. The identities are disjoint; ancestry balance,
relatedness and clinical representativeness are not established by this selection.

Indexed requests retain 636 source VCF records across 91 merged intervals surrounding
every model locus. Each locus receives 100 bases of context on either side, with overlapping
VCF records retained. Both reference and target data use the same unchanged source GTs;
GT-only extraction explicitly uses `--accept-missing-quality`. The reference FASTA used in
the checked run is the NCBI GRCh38 no-alt analysis set; its exact digest and index digest
are recorded. An existing local reference avoids downloading it again.

The experiment has 12 scenarios, comprising six masks with and without hypothetical
frequency filling:

- Random 5% and 20% sample/site masking, each with seeds 42 and 271828. The documented
  SHA-256 mask algorithm is reproduced separately in Python; zero-contribution masked
  genotypes still count as masked observations.
- A contiguous dropout of GRCh38 chromosome 9 `[21,500,000, 22,500,000)`. This models a
  regional loss, with the interval explicitly chosen; it is not an empirical failure map.
- A common-variant coverage mask, retaining sites with minor-allele frequency at least
  5% in the 96 reference samples. This is a reference-derived coverage approximation,
  not a particular genotyping array or an imputed dataset.

Allele frequencies are computed only from the 96 reference samples, frozen before masking
and written as decimal text with 100-digit Decimal precision. Original allele counts and
counts of called chromosomes are retained. The oracle verifies those counts against the
source GTs. It uses the exact retained decimal strings for arithmetic, so division rounding
is explicit and is never confused with an exact rational population frequency.

Each scenario keeps the model's same common-callable term set. External-panel comparisons
also use this set, so a change of comparison samples cannot silently change variant coverage.
Percentiles use empirical mid-CDF ranks against the frozen available-data baseline; within-cohort
reranking is reported separately. These values describe numerical movement, not calibrated
population risk. Cutoffs of 90 and 95 are descriptive experiment settings.

## What the independent checks establish

The separate Python oracle does not take engine packs, extracted tables or reported
contributions as its arithmetic inputs. It independently resolves literal SNV alleles
against the source FASTA, checks original VCF calls for all 128 samples, determines the
same common-callable term sets, and verifies their ordinal/coordinate fingerprints.
Python `fractions.Fraction` then reconstructs scores exactly from the original decimal
weights and GT dosage, and from the frozen reference-frequency strings when filling.

The run checks:

- 1,536 sample/scenario scores and their baseline sums, exact differences, observation
  counts and attribution reconciliation.
- 57,984 independently generated sample/term masking decisions.
- All 4,063 nonzero changed term contributions and 1,536 top-attribution rows, including
  signed ordering and the remainder outside the top ten.
- All 1,536 regional deltas and 3,072 external-reference percentile rows.
- Mean absolute error, RMSE, maximum change, Pearson and Spearman correlations,
  fixed-baseline percentile movement, within-cohort reranking and cutoff crossings.

Scores, fills, contributions and deltas require exact rational equality. Floating-point
summary statistics and displayed percentiles use absolute or relative tolerance `1e-12`.
Input/output hashes are rechecked. The runner refuses Python `-O`, which would disable
its assertions, and writes `results.json` only after all gates pass.

A useful result is that mean filling need not improve every metric. For PGS000011 with
20% random masking and seed 42, mean absolute score error decreases from **0.6782 to 0.1660**.
For its common-variant coverage mask, mean absolute error instead increases from
**0.04447 to 0.06067**, even while RMSE decreases from **0.11691 to 0.10925**. This is an
observed property of these conditional scores and selected samples, not a general
recommendation about missing-data methods. The complete results include every scenario.

## Reproduce

Use the Python environment pinned for the existing public chromosome-22 benchmark
(Python 3.12, `pysam` and `zstandard`), and an indexed GRCh38 no-alt analysis FASTA.
The runner downloads small scoring files and metadata, then uses indexed remote VCF
queries; it does not download entire chromosomes. The reference FASTA is a separate
input of roughly 3 GB uncompressed. Do not place public data or scoring files in Git.

```sh
python3.12 -m venv .venv-public
.venv-public/bin/pip install -r bench/public_chr22.requirements.txt
cargo build --locked --release

.venv-public/bin/python bench/robustness_public.py fetch --out data/robustness-sources
.venv-public/bin/python bench/robustness_public.py prepare \
  --data data/robustness-sources --reference /path/to/GRCh38_no_alt.fasta \
  --out data/robustness-inputs
.venv-public/bin/python bench/robustness_public.py native \
  --prepared data/robustness-inputs --pgsum target/release/pgsum \
  --out data/robustness-native
.venv-public/bin/python bench/robustness_public.py verify \
  --prepared data/robustness-inputs --native-root data/robustness-native \
  --pgsum target/release/pgsum --out data/robustness-verified
```

Every stage output directory must be new except the source cache used by `fetch`.
Score files and Catalog metadata have pinned SHA-256 hashes; the selected public
records have a pinned semantic GT fingerprint independent of gzip encoding. An upstream
change causes a strict failure and needs inspection before updating pins. Preparation
records all interval queries, retained source files, sample identities, source-index
hashes, reference identity and output hashes. Every native command and its log are retained.
The native output includes the editable experiment manifest, HTML report, full per-sample,
reference, regional and changed-term tables. The scoring and verification executable
must have the same fingerprint.

This validation covers these four sparse additive models, diploid GTs and controlled
missingness. It does not establish genome-wide dense-model performance, phenotype
prediction, clinical calibration, ancestry normalization, empirical sequencing-quality
thresholds, array imputation performance or performance in an independent population.
