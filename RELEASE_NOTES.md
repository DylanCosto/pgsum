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
