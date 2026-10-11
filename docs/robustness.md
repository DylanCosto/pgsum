# Reproducible PGS robustness experiments

`pgsum robustness` measures how genotype availability, replacement policies, alternative
extractions, and reference distributions change a score. It produces exact score differences,
rank diagnostics, cutoff crossings, and variant/region explanations. It does not estimate
phenotype prediction, ancestry portability, or absolute risk.

## Run an experiment

First compile score packs and extract the baseline cohort with the normal pipeline. A
single-sample `.pgsg` also works; meaningful population ranks then require a reference set.

```sh
pgsum extract --gvcf cohort.vcf.gz --all-samples --reference GRCh38.fa \
  --pack packs/ --out cohort.pgsc
pgsum robustness --experiment experiment.json --out experiment-results --terms
```

Use a **new** output directory for each invocation. Relative paths in the manifest resolve
against the manifest's directory, not the shell's working directory. Packs must be explicit
file paths, not directories. Example `experiment.json`:

```json
{
  "schema": "pgsum-robustness-experiment-v1",
  "genotypes": "cohort.pgsc",
  "packs": ["packs/PGS000011.pgsp"],
  "cutoffs": [90, 95],
  "top_terms": 10,
  "scenarios": [
    {"id": "control"},
    {
      "id": "missing_5_percent",
      "masks": [{"kind": "random", "rate": 0.05, "seed": 42, "unit": "sample_site"}]
    },
    {
      "id": "missing_5_percent_filled",
      "masks": [{"kind": "random", "rate": 0.05, "seed": 42, "unit": "sample_site"}],
      "fill_frequencies": "independent-reference.pgsf",
      "frequency_source_note": "Frequencies computed from the separate reference participants listed in the study manifest; no evaluation participant was included."
    },
    {
      "id": "restricted_marker_coverage",
      "masks": [{"kind": "bed", "path": "assayed-markers.bed", "mode": "keep"}]
    },
    {
      "id": "regional_dropout",
      "masks": [{"kind": "bed", "path": "difficult-regions.bed", "mode": "drop"}]
    }
  ],
  "reference_sets": [
    {"id": "reference_A", "genotypes": "reference-A.pgsc"},
    {"id": "reference_B", "genotypes": "reference-B.pgsc"}
  ],
  "regions": [{"id": "region_of_interest", "path": "regions.bed"}]
}
```

Frequency files, reference sets, region annotations and the complete term audit are optional.
Remove their entries if those files do not exist. A reference set can select specific
`"samples": ["S1", "S2"]`; omitting its `genotypes` selects from the baseline cohort.
The report records overlap in sample IDs rather than representing an overlapping set as
independent. Changed IDs cannot prove independence; study design remains the researcher's
responsibility. A frequency source note is required whenever filling is requested. Its contents
are recorded as a declaration, not verified as a fact.

Optional `scoring` has the four boolean production opt-ins:
`allow_inferred_other_allele`, `allow_inferred_palindromes`, `allow_inferred_indels`,
and `accept_informational_descriptions`. Defaults are false. They apply identically across
every arm; masks never override term review or orientation rules.

## What is held fixed

For each score, the baseline includes only autosomal terms that are scorable in **every
baseline sample and every selected reference sample**, using the declared scoring rules.
This conservative intersection gives all comparisons an identical, explicit variant set.
Pre-existing missing genotypes are not invented and are never treated as hidden truth.
The report retains the original model inventory, exclusion reasons, term/weight coverage,
and a digest of the included source ordinals and normalized target positions. Common-callable
coverage can shrink substantially in large or heterogeneous cohorts.

This is an available-data, shared-term baseline. Even when the published scoring file is
complete, the experimental baseline may be only a subset. A no-mask control equals the
ordinary cohort score only when all its model terms pass these eligibility rules.
Inventory inconsistency and publication metadata are retained, not turned into validation.
No eligible baseline terms is an error, with no completed report.

Baseline genotypes and reference distributions never change when a target is masked.
Reference comparisons are the **frozen-common-terms** analysis: adding another reference
set can reduce the common universe, so comparisons across separate experiments must also
check the baseline-term digest. The first release does not silently rerun the production
panel-matching pipeline with a different term set in each arm.

## Masks and alternative extractions

- `random` / `sample_site` independently hides positions for each sample identity.
- `random` / `site` hides the same positions in every target sample.
- `bed` / `keep` hides positions outside the intervals; `drop` hides positions inside them.
- Multiple masks combine as a union of hidden positions.
- Rates range from 0 to 1. Use explicit scenarios with different seeds for replicates.
  Identical seeds reuse identical patterns; increasing the rate nests the hidden sets.

BED coordinates are zero-based, half-open. A one-base marker at position 101 uses `chr1 100 101`.
Contig aliases such as `1` and `chr1` are equivalent. Sequence variants use their normalized
genotype-target anchor, not their original scoring-file coordinate or full allele span.
Terms or models using that same position receive the same mask. JSON variant explanations
retain both normalized `pos1` and scoring-file `source_pos1`.

The versioned SHA-256 mask protocol is specified in `src/robustness/config.rs` and recorded
in every report. It depends on seed, canonical chromosome/position, mask unit, and (for
`sample_site`) the UTF-8 sample identity. It is independent of sample order, pack order,
scenario ID and thread count.

A scenario's optional `genotypes` path compares an independently extracted table with the
same sample identities (ordering may differ), the same packs, and the same reference identity.
This supports GQ/DP policy or GT/DS/GP comparisons after re-extraction. Cached genotype outcomes
do not retain all raw quality measurements; the command therefore does not pretend to change
quality thresholds on an existing cache. Alternate input cannot expand the frozen baseline.
Native GT, additive DS and supported GP expectations retain the existing scoring semantics.

## Replacement policies

Without `fill_frequencies`, hidden or unavailable alternate contributions are omitted (zero
added). With a frozen frequency table, eligible missing additive terms receive `weight × 2 ×
effect-allele frequency`, using the existing strict allele/frequency matching. Nonadditive,
ambiguous, unresolved or unsupported terms remain unfilled, with reasons. Observed contributions
are never overwritten by frequency fills.

These are explicitly **experimental reconstructions**, including below the production 99%
coverage gate. The command does not modify ordinary `score` behavior or emit a production
`filled_score`. Estimate frequencies from disjoint reference participants or an appropriate
external source, never the hidden evaluation calls. Genotype masking does not simulate linkage-
based imputation or every source of sequencing/array error.

## Outputs and interpretation

| File | Contents |
|---|---|
| `report.html` | Self-contained visual summary, scope, cutoff crossings, and declared settings |
| `report.json` | Completion marker, experiment, source metadata, software/input/output hashes, per-score summaries |
| `samples.tsv` | Exact baseline/scenario sums and deltas, coverage counts, fixed-distribution and reranked percentiles |
| `top-changes.jsonl` | Largest absolute variant deltas per sample/arm, with signed contributions and exact remainder |
| `regions.tsv` | Signed deltas inside each supplied region; overlapping regions must not be summed |
| `references.tsv` | Baseline/scenario placement against each unchanged reference distribution and sample-ID overlap |
| `terms.tsv.zst` | Optional complete **nonzero change** audit, in source-term order; potentially large |

Score differences, term attribution, ties, and rank ordering use exact decimals. Summary error
and correlation metrics use floating point only after exact differences/centering; unrepresentable
metrics are null with reasons. Raw errors cannot be pooled across different PGS scales.

Percentiles use empirical mid-CDF: `100 × (number below + 0.5 × number equal) / N`.
Fewer than two values or a constant reference distribution produces no percentile. Values beyond
the reference range can produce 0 or 100. Fixed-baseline percentiles measure numerical movement
relative to the unchanged baseline; `reranked_percentile` uses the scenario cohort's distribution.
Spearman correlation ranks baseline and scenario separately. Cutoff membership is `percentile >=
cutoff`; crossings count paired changes using the fixed baseline distribution.

Partial or perturbed scores placed against an unmodified distribution are **diagnostics**, not
valid risk percentiles. A sensitivity range is not a confidence interval. Stable scores and high
coverage do not establish predictive accuracy. Predictive discrimination, calibration, and
population portability require measured outcomes in an appropriate independent cohort.

Each sample's score delta is checked against the sum of **all** changed contributions, then
split into the reported top terms and an exact remainder. This explains raw-score changes;
percentile shifts have no additive per-variant decomposition.

The final report is published atomically after output hashes and a second input-hash check.
An interrupted or failed run may leave diagnostic files but no `report.json`; it is incomplete.
The command refuses to overwrite an existing output directory. Inputs must remain immutable
while mapped by the existing genotype readers.

## Validation

Run `cargo test --locked` for arithmetic, masks, schema validation, input-policy parity,
attribution and CLI regression checks. The public coronary-artery-disease experiment is
documented in [public-robustness-validation.md](public-robustness-validation.md), including
the complete score inventories, actual shared-term coverage, independent oracle and its limits.
