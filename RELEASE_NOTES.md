# Scoring workflow development (unreleased)

Added `compare` for pgsum partial scores, explicitly selected PLINK 2 sum columns and pgsc_calc `SUM`
outputs, including gzip/zstd input, exact decimal differences/tolerances, complete identity joins and
an optional failure gate. Paired term exports stream through evidence comparison and reconciliation
against each reported sum. Reports record source/output hashes and keep unavailable values separate
from zero. See `docs/compare.md`; public real-data cross-tool validation remains outstanding.

Single-sample/cohort extraction and `run` accept `--dosage-field gt|ds|gp`. DS additive contributions and GP
model-specific expectations retain exact decimal arithmetic. Missing selected fields do not fall back
to GT, and confident reference blocks retain their existing rules. See README for current biallelic,
ploidy, reference-panel and strict probability-mass limitations.

Score JSON is now v4, with selected field and expected-contribution counts. Genotype files are v5;
v1–v4 remain readable. See DESIGN for migration details. This is one part of the broader workflow
expansion, not a completed release of all planned features.

Native streaming BCF now uses the same GT/DS/GP rules for single samples and cohorts. Indexed BCF and
uncommon codec encodings remain unsupported. BCF dosage precision follows its stored float values.

`batch` accepts sample sheets with bounded worker concurrency, atomic sample generations, verified
resume, incremental additions, aggregate output and per-sample failure logs. It does not yet enforce a
total memory budget. `run --ids` and `batch --ids` download and compile missing Catalog scores into a
verified reference-specific cache; `--offline` and `--refresh-catalog` make reuse/freshness explicit.
`run --preflight` validates headers/packs before extraction. Successful runs write a versioned execution
manifest with software, source, policy and output identities. See README for commands and limitations.

Single-sample score JSON now includes a missingness report with exact contribution envelopes, ranked
missing terms, exclusion reasons and retained source evidence. Ambiguous models and unsupported ploidy
withhold whole-score bounds. These diagnostics preserve existing raw/partial results and remain separate
from filling and panel placement. Per-term TSVs add bound columns before optional frequency-fill columns.
New cohort extractions write v3 with compact exclusion reasons and original source-record ordinals.
`score --missingness` optionally writes cohort reports with the same bounds/rankings; numeric scoring
continues to skip diagnostic frames. Readers preserve v1/v2 scoring and label unrecoverable legacy
missing reasons explicitly. Cohort extraction now also checks declared contig lengths against the
reference, matching the single-sample safeguard. See README for report provenance and memory limits.

Cohort scoring now loads decoded blocks through an evicting cache, defaulting to 128 MiB of retained
payload. `score --cohort-cache-mib` adjusts it; score metadata records cache activity and scope.
Workers reuse their current block, and scoring no longer eagerly decodes the whole cohort. Ordered
quantitative scores also share decoded waves, with source-order fallback for unordered targets.
Source-reference reporting gathers evidence in key order to avoid repeated probability-block decoding. The
cache is not a total process RAM cap; see README for active-block/decoder overhead and ordering tradeoffs.

# Next release (unreleased)

pgsum scores individual gVCFs and cohort VCFs with exact decimal sums. This release improves cohort
scoring and makes the performance and correctness checks reproducible without downloading genomes.

## Scoring performance

- Skip four homozygous-reference cohort samples at a time from the packed genotype representation.
- Precompute exact differences from the reference-dosage contribution once per term, instead of
  subtracting and adding full contributions for every non-reference sample. Wide weights and different
  dosage-weight exponents retain the general arithmetic path.
- Preserve score text, missing-call handling, coverage, and decimal precision. No file-format change is
  needed for existing packs or genotype tables.
- `bench/performance.py` measures compilation, extraction, scoring, and gVCF end-to-end runs, verifies sums
  against an independent integer oracle, and compares outputs across thread counts and optional baseline
  binaries. See `bench/README.md` for measurements and reproduction commands.

## Score JSON migration

Score results now use `pgsum-score-v3`. Consumers of v2 should update before upgrading:

- Reference-panel terms below 99% panel call rate are excluded from both the sample and reference sums
  and from matched coverage. `reference.excluded_low_call_rate_terms` replaces
  `reference.fills.terms_constant`. Percentiles and coverage can change.
- Undefined Z-scores are `null`, including groups with fewer than two samples or zero variance.
- `imputation_performed` accounts for frequency filling. The new `imputation` fields distinguish sample
  filling, reference-panel filling, and observed scores. Observed strict and partial sums remain unfilled.

Full details are in `DESIGN.md`. Historical real-genome reference-panel results have not been rerun under
this exclusion rule; the public synthetic tests verify the arithmetic and boundary behavior.

## Release checks

CI tests Linux on stable Rust and the minimum supported Rust 1.88, plus macOS on stable Rust. A release-mode
performance smoke test checks exact outputs and saves timing/memory results; shared-runner timings are
informational. Release packaging checks that the tag matches the package version, tests native targets
with locked dependencies, and smoke-tests their binaries. Cross-compiled binaries are built but not run.

Before tagging, review the checks, set the package version and lockfile to the intended release, and finalize
these notes; release packaging includes them automatically. This file does not announce a published release.
