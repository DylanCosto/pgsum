# 0.6.0 candidate (unreleased)

This is release preparation, not an announcement of a published release. The latest
previous tag is v0.5.1. The original dense genome-wide validation is complete. The additional
paired-study checks below accompany this follow-on candidate; no version has been tagged.

## Robustness experiments

`pgsum robustness --experiment experiment.json --out results` runs controlled sensitivity
experiments from a JSON manifest. It supports deterministic random sample/site or whole-site
masking, BED coverage/dropout masks, frozen population-frequency filling, alternative
genotype extractions, and comparisons against fixed reference distributions.

Every arm uses the same explicit common-callable autosomal term set. Existing missing calls,
model/orientation exclusions, and incomplete coverage remain visible. Results include exact
scores and deltas, independent reranking and fixed-baseline percentiles, cutoff crossings,
variant/region attribution, a self-contained HTML summary, and input/output/software hashes.
Reported top contributions plus their exact remainder must reconcile to every score delta.
Normalized indel anchors receive consistent masks across source representations.

Filling is an experimental reconstruction, including below production coverage gates. It
never changes ordinary scoring rules, replaces observed calls, or uses hidden evaluation
calls to estimate frequencies. The source's independence is a recorded user declaration.
See the [experiment guide](docs/robustness.md) for the manifest and interpretation.

## Paired robustness studies

`robustness-plan` expands authored policies across omission rates and seeds while retaining
controls. `robustness-summarize` verifies the completed experiment and compares score errors,
correlations, cutoff transitions and frozen-reference percentile shifts. Comparisons pair
policies within each seed before calculating descriptive distributions. Undefined metrics and
participant denominators remain explicit. See the [study guide](docs/robustness.md).

Repeated scenarios share decoded frequency tables and reuse deterministic random draws at
each locus. These changes preserve the public pilot's numerical outputs. The recorded pilot
memory comparison fell from 2,060.48 to 501.14 MiB peak RSS; observed timings are shared-host
measurements, not controlled speed claims. The follow-on passes all 179 tests on stable Rust and
Rust 1.88, formatting, Clippy and twelve synthetic full-attribution runs.

The completed [128-person study](docs/robustness.md#completed-128-person-public-study)
checks 42 scenarios across four omission rates and five seeds. All 5,376 participant/scenario
rows, 294 score-metric checks and 300 distribution checks pass. Its exact baselines match the
source-verified genome-wide run; new masked scores are checked descriptively, not independently
reconstructed from source genotypes. Frozen-reference filling lowers mean absolute error at
each tested rate in this cohort, without establishing predictive or clinical validity.

The optional offline paired study uses Python 3.9+ and independently checks its summaries:

```sh
python3 examples/robustness/run.py --binary ./pgsum --out example-output --study
```

## Progress and participant-group diagnostics

Long robustness runs now report stages and processed term counts, with human-readable, JSON
or quiet modes. Completion is emitted only after final input verification and publication of
the report. Progress settings leave scientific outputs unchanged.

`robustness-groups` summarizes completed experiments for an explicit, exhaustive participant
partition. Score errors and correlations use within-group pairs; percentile shifts and cutoff
transitions retain the same full-cohort baseline for every group. Reports include group and
scenario filters, unavailable reasons, membership and input/output fingerprints. The
[checked 128-person technical partition](bench/results/group-summary-public128.json) passes
588 score-metric, 168 cutoff and 84 reference comparisons. These artificial groups validate
reporting only; they do not establish demographic or ancestry differences. The extracted
offline study also checks a singleton group and a three-person group.

## Scoring and workflow changes since v0.5.1

- Dense model compilation now accepts legitimate IUPAC ambiguity in the reference FASTA.
  Ambiguous bases remain unknown, affected terms stay unresolved, and indel normalization
  cannot invent an ambiguous anchor. Invalid non-IUPAC characters still fail. This fixes
  compilation of the full PGS000018 file against the official GRCh38 reference; compiled
  Catalog caches are invalidated by the updated rule identity.
- `compare` joins identities explicitly and compares pgsum partial sums, selected PLINK 2
  sum columns, and pgsc_calc `SUM` results with exact decimal differences and tolerances.
  Paired term exports support attribution and reconciliation. See [comparison guidance](docs/compare.md).
- GT, additive DS and model-specific GP expectations preserve decimal arithmetic. A
  missing selected field does not silently fall back to GT. Native streaming BCF uses
  the same rules as VCF; indexed BCF remains unsupported.
- `batch` adds sample sheets, bounded workers, verified resume, incremental sample runs,
  aggregate tables and per-sample failure logs. Catalog-backed `run --ids` and `batch --ids`
  support explicit offline reuse/refresh, preflight checks and execution manifests.
- Missingness reports retain exact contribution envelopes, reasons, source evidence and
  term rankings. Cohort reports are optional; ordinary scoring skips diagnostic frames.
- Cohort readers use an evicting decoded-block cache and share ordered quantitative passes.
  Packed homozygous-reference calls and exact contribution differences reduce repeated work.
  The cache setting is not a total process-memory limit.

## Compatibility and interpretation

Current score output is `pgsum-score-v4`, genotype tables are v5, and new cohort extractions
are v3. Older genotype/cohort formats remain readable as described in [DESIGN](DESIGN.md);
legacy cohorts cannot recover diagnostic reasons that were never stored. Consumers upgrading
from score v2/v3 should review the migration notes, selected dosage fields, imputation fields
and optional missingness outputs before deployment.

Reference-panel terms below 99% panel call rate are excluded from both comparison sides.
Undefined Z-scores are null. Observed raw/partial sums remain separate from frequency fills.
Historical reference-panel measurements have not all been rerun under the current exclusion
rule; do not interpret historical and current percentiles as interchangeable.

Robustness describes numerical reconstruction and rank sensitivity. It does not establish
phenotype prediction, clinical calibration, ancestry portability or valid risk percentiles
for incomplete models. Continuous ancestry normalization and independent outcome validation
are outside this release.

## Validation evidence and candidate gates

The checked [public robustness experiment](docs/public-robustness-validation.md) retains the
published inventories of four sparse coronary-artery-disease models and uses 32 target plus
96 disjoint reference samples. A separate source-genotype oracle verifies exact scores,
mask decisions, fills, attribution, ranks and reference comparisons. Its conditional
common-callable baselines are explicitly partial scores, not complete published PGS values.

The [genome-wide public benchmark](docs/genomewide-robustness-validation.md) compiles all
1,745,179 actual rows of PGS000018 and checks all 22 autosomes in nested 8/32/128-person
cohorts with 64 disjoint reference people. All 840 sample/scenario checks pass independent
exact reconstruction, attribution, count and percentile checks. The 128-person common-callable
baseline retains 1,712,276 terms; five scenarios took 119.46 seconds and 585.40 MiB peak RSS
on the recorded host. These are conditional scores and observational timings, not evidence
of predictive accuracy or a controlled speed comparison. The upstream header's one-row
inventory discrepancy remains explicit.

The separate [public chromosome-22 comparison](docs/public-validation.md) has already run
the actual pgsc_calc workflow and PLINK 2 alongside pgsum. Its shared-term arithmetic checks
cover chromosome contributions and do not establish whole-genome or phenotype performance.

The original Linux candidate passed formatting, Clippy, and all 169 tests on both stable Rust
and the Rust 1.88 minimum; one optional local-data parity test is ignored. Six deliberate
output-corruption checks and 16 smaller full-attribution benchmark runs also pass. Hosted
CI covers macOS/Linux tests and extracted-package checks; the release workflow adds native
package checks for all four supported targets. The checked records and their limits accompany
the documentation.

The [synthetic genome-wide scale run](bench/README.md#checked-scale-run-2026-10-11)
passed 24 runs across 10,000, 100,000 and 1,000,000 model terms, 32/128 samples,
zero/0.5% requested missingness and one/four threads. The million-term, 128-sample
complete-data case took 29–30 seconds and about 106 MiB for three scenarios on the
recorded shared host with uncontrolled background load. Exact evidence matched across thread counts. Coverage attrition
is reported separately; these synthetic observations do not imply real-data performance.

## Complete offline packages

Linux x86-64/ARM64 and macOS x86-64/ARM64 targets remain supported by the release workflow.
Archives now include the binary, documentation, needed benchmark scripts and checked result
records, licences, a payload hash manifest, and an entirely synthetic runnable example.
Real Catalog scoring files and participant genotypes are not redistributed in the package.

After unpacking an archive, run the example with Python 3.8+ and no additional packages:

```sh
python3 examples/robustness/run.py --binary ./pgsum --out example-output
```

It compiles the miniature model, extracts target/reference cohorts, masks/fills genotypes,
and checks hand-calculated scores, ranks, counts, regions and output hashes. Outputs go to a
new directory; inputs remain unchanged. Each of the four release targets has a native build
runner configured to execute this example using the extracted binary, including paths with
spaces and a different working directory. macOS uses the supported macOS 15 ARM and Intel
images following the [macOS 14 runner retirement notice](https://github.com/actions/runner-images/issues/13518).

Each archive has an adjacent SHA-256 checksum and an internal `package-manifest.json` with
every payload file's hash, size and mode. A local archive can be built and checked with:

```sh
python3 scripts/package_release.py build --binary target/release/pgsum \
  --target x86_64-unknown-linux-gnu --output dist --run-example
python3 scripts/package_release.py verify dist/pgsum-v0.6.0-x86_64-unknown-linux-gnu.tar.gz \
  --run-example
```

The helper fixes archive ordering, permissions, owner fields and timestamps. Identical
payload bytes, version, target, `SOURCE_DATE_EPOCH`, and packaging/compression runtime
produce identical archives; this does not assert reproducible compilation. The default
epoch is zero. Local Markdown
links/anchors are checked against the packaged files, including linked benchmark records.

Manual release-workflow dispatch builds and validates artifacts without publishing.
Publication is limited to a version-tag push that matches Cargo's package version and
passes the release jobs. No tag or release is created by package preparation.
