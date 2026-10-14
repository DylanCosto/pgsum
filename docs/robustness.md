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

Try the full workflow without downloading data using the included
[offline paired study](../examples/robustness/README.md#paired-missingness-study):

```sh
python3 examples/robustness/run.py --binary ./pgsum --study --out study-example-output
```

## Monitor a long experiment

By default, `robustness` writes readable stage changes and scoring updates to stderr.
Updates show the score's index, participant/scenario counts, processed terms, retained
common-callable baseline terms, exclusions and total elapsed time. The default scoring
interval is 30 seconds; use a positive `--progress-interval` to change it:

```sh
pgsum robustness --experiment experiment.json --out results \
  --progress json --progress-interval 10 2> progress.jsonl
```

`--progress json` emits one `pgsum-robustness-progress-v1` JSON object per update, with
`phase`, `elapsed_seconds`, a `score` object when processing a score, and a report path
only on `complete`. `--progress quiet` suppresses these updates. Errors still appear on
stderr as ordinary diagnostics; workflow tools should handle a nonzero exit status and
separate JSON progress objects from error lines.

Stages are `preparing`, `fingerprinting`, `scoring`, `summarizing`, `writing`,
`verifying_inputs` and `complete`. Per-score counts reset at the next score. Processed
terms include exclusions and count completed terms only; they are not a percentage of
elapsed runtime. Interval checks occur between blocks of 256 terms, and stage boundaries
emit immediately. This is not a fixed-frequency heartbeat during input loading, hashing,
summarization or a slow term. No completion-time estimate is made.

Finishing a term scan does not mean the experiment succeeded. A `complete` event is emitted
only after the final input verification and publication of `report.json`, the authoritative
completion marker. An interruption can leave partial files without that marker. Progress
settings and timestamps never enter scientific output files or their fingerprints.

The [checked public progress run](../bench/results/progress-public-pilot.json) scans the full
1,745,179-row PGS000018 model with eight targets and 14 scenarios. Its structured updates
reconcile 1,712,283 baseline terms and 32,896 exclusions, and all numerical output bytes and
scenario summaries match the earlier checked public pilot.

## Compare participant groups

`robustness-groups` reuses a completed experiment to show how numerical sensitivity differs
across an explicit participant partition, without reloading genotypes or rerunning masks:

```sh
pgsum robustness-groups --results study-results --groups participants.tsv --out group-results
python3 bench/check_group_summary.py study-results participants.tsv group-results
```

The mapping is an unquoted TSV with exactly `sample` and `group` columns:

```text
sample	group
participant_A	assay_batch_1
participant_B	assay_batch_2
participant_C	assay_batch_1
```

Use literal tabs. Include every target participant exactly once. Missing, extra, duplicate or
blank identities are rejected; row order does not matter. Labels are supplied by the researcher
and are never inferred from names, genotype data or score values. Analyze another grouping
variable with a separate mapping and output directory. This command does not normalize ancestry
or validate the truth of demographic labels.

Each group keeps the original whole-cohort common-callable term set. Score-error metrics and
correlations use that group's aligned pairs, but **percentile shifts and cutoff membership use
the same frozen full-target baseline distribution for every group**. Recalculating percentile
cutoffs separately inside each group would change their meaning and hide some differences.
Spearman correlation still ranks the selected pairs within each group. A singleton can have a
valid full-cohort percentile shift while its correlation is unavailable.

External-reference summaries retain the original reference distribution. The field
`source_cohort_shared_sample_ids` reports overlap with the entire source cohort; it is not a
count of overlap within an individual group. Evaluable/unavailable counts and null metric
reasons remain explicit. Different group sizes, label selection and conditional baselines
limit interpretation; these are descriptive summaries, not group-effect significance tests,
predictive performance or causal comparisons.

Outputs are `metrics.tsv`, `cutoffs.tsv`, `references.tsv`, `groups.html` and the `groups.json`
completion marker. The HTML report filters by participant group and scenario. JSON retains
group memberships, baseline identity, source/output fingerprints
and unavailable reasons. Native output hashes are verified before analysis; inputs are checked
again before completion. Partial or changed evidence cannot publish a successful group report.
Each scenario remains separate; this command does not pool seeds or select a best policy.

The offline `--study` example also demonstrates an explicit singleton/three-person partition
and runs the independent group checker. The checker derives score/cutoff statistics from exact
native sum strings and frozen-reference summaries from participant reference rows. It does not
reconstruct these scores from original genotypes.

The [checked 128-person group report](../bench/results/group-summary-public128.json) splits the
public study into two artificial 64-person groups using target selection order. Across 42
scenarios, the independent checker passes 588 score-metric, 168 cutoff and 84 reference
comparisons. This validates the reporting path; the technical split supports no ancestry or
demographic conclusion.

## Paired policy comparisons within groups

Pass the same exhaustive mapping to `robustness-summarize` to compare repeated policies
within each group without rerunning the experiment:

```sh
pgsum robustness-summarize --plan planned-study --results study-results \
  --groups participants.tsv --reference-policy omit --out group-study
python3 bench/check_group_study.py planned-study study-results participants.tsv group-study
```

The command verifies native results and calculates each group's scenario metrics before
pairing policies by rate and seed. Every group keeps the original full-target percentile
thresholds and common-callable terms. Error and correlation calculations use selected pairs
within the group. Groups never become replicates, and participants are not pooled across seeds.

This opt-in output uses `pgsum-robustness-group-study-v1`. Each JSON score entry has `group`,
`samples` and `source_samples`; the four TSVs add those dimensions after `pgs_id`. JSON records
explicit memberships, per-seed comparisons, availability counts, controls and fingerprints.
`reference_sets` retains `source_cohort_shared_sample_ids`, which describes the entire source
cohort even when it exceeds a group's size. The HTML tables identify groups and both sample counts.
Without `--groups`, the existing whole-cohort schema and table columns stay unchanged.

The independent checker reconstructs group errors and whole-cohort cutoff membership from
exact native score strings, then checks paired distributions, missing metrics and reference
shifts. These summaries describe masking variability within an authored partition. They do not
estimate uncertainty about populations, validate group labels, or establish causal group effects.

The [checked 128-person group study](../bench/results/group-study-public128.json) passes 168
score-metric and 432 cutoff/reference distribution checks across two artificial 64-person groups,
four omission rates and five seeds. Whole-cohort outputs remain byte-identical without grouping.

## Plan repeated sensitivity experiments

`robustness-plan` expands policies over a Cartesian product of missingness rates and
explicit seeds. Planning reads only the template; it does not load genotypes or run scores.
For example, save this as `policies.json` beside your extracted inputs:

```json
{
  "schema": "pgsum-robustness-experiment-v1",
  "genotypes": "cohort.pgsc",
  "packs": ["packs/PGS000011.pgsp"],
  "scenarios": [
    {"id": "omit"},
    {
      "id": "fill",
      "fill_frequencies": "independent-reference.pgsf",
      "frequency_source_note": "Frozen frequencies from disjoint reference participants."
    }
  ]
}
```

```sh
pgsum robustness-plan --template policies.json --rates 0.01,0.05,0.1,0.2 \
  --seeds 41,42,43,44,45 --mask-unit sample-site --out planned-study
pgsum robustness --experiment planned-study/experiment.json --out study-results
```

This produces 42 scenarios: the two original policies as controls plus two policies × four
rates × five seeds. `plan.json` maps each generated scenario to its policy, rate and seed,
and records the template and generated-manifest hashes. `experiment.json` is an ordinary,
explicit experiment: inspect it before running, and retain `plan.json` beside it. Relative
input paths are resolved to absolute paths against the template's directory so the plan can
run from another working directory. Relocating the underlying inputs requires a new plan.
The output directory must be new; existing files are never overwritten. A completed plan
has `experiment.json`; a failed write may leave diagnostic files without that completion marker.

Matching seeds reuse the same random draws across policies, allowing paired comparisons of
omission and replacement. The same seed across rates nests the randomly hidden positions.
Template BED masks and alternative genotype inputs remain in each arm, so total missingness
may additionally reflect those policies. Template scenarios must not already contain random
masks. Rates and seeds must be distinct; rate zero and one are accepted as boundary checks.
The default expansion limit is 10,000 scenarios, including controls; `--max-scenarios` changes
that limit before expansion. Plan size does not estimate the memory required to run it.

Use `--mask-unit site` to hide the same positions across participants, or the default
`sample-site` for participant-specific omissions. Site-level filling adds the same offset to
each participant when the same terms are filled, so it cannot repair within-cohort ordering.
Replicate variation describes sensitivity to the chosen masking design, not population
uncertainty, predictive validity or clinical confidence intervals. Reports retain individual arms; the study summarizer below describes seed variation without
selecting a best policy.

## Summarize a completed study

```sh
pgsum robustness-summarize --plan planned-study --results study-results \
  --reference-policy omit --out study-summary
```

The command verifies the generated manifest's fingerprint, checks that the completed report
used that experiment and mask algorithm, and verifies its output-file hashes. Each planned
arm must match its actual policy, rate, seed and mask unit; policies must share a complete
rate-by-seed grid. Missing, duplicated or incompatible arms are errors. The command accepts
one completed experiment, keeping every model's common-callable baseline and scale separate.
The original genotype files need not be reopened for this descriptive summary.

`summary.html` presents four tables. `summaries.tsv` describes each policy/rate/metric across
seeds, with the requested, available and unavailable counts; mean, median, minimum, maximum
and sample standard deviation. The seven metrics are MAE, RMSE, maximum absolute score
change, Pearson and Spearman correlation, and mean/maximum absolute fixed-percentile shift.
`comparisons.tsv` describes **policy minus reference policy** differences, calculated within
matching seeds before aggregation. Lower differences favor the policy for error and shift
metrics; higher differences favor it for correlations. These are descriptive comparisons,
without significance tests or automatic policy selection.

A missing metric is never treated as zero. A paired difference requires both values and a
representable difference. Sample SD is unavailable with fewer than two available seeds;
unrepresentable aggregate values are also null. Identical boundary masks can give zero
variation even with distinct seeds. This does not establish certainty about real data.

`summary.json` includes per-seed paired values, controls, model coverage and baseline digests,
source fingerprints and output hashes. Controls remain separate from random-mask replicates.
`percentile_summaries.tsv` and `percentile_comparisons.tsv` add two structured scopes:

- **Cutoff:** upward, downward, unchanged and total changed participant counts; upward,
  downward and total changed fractions. Membership uses `percentile >= cutoff` against the
  frozen target baseline distribution. Each declared cutoff stays separate. Fractions range
  from 0 to 1; paired fraction differences are in fraction units, not percent units.
- **Reference:** mean baseline/scenario percentile and mean/maximum absolute shift against
  each authored frozen reference distribution, in percentile points. Each named reference
  set stays separate, with reference sample count and target overlap in `summary.json`.

The JSON retains each seed's evaluable participant count and both denominators for each
paired difference. With no evaluable participants, cutoff counts **and** fractions are null;
this must not be read as evidence of no crossings. Reference metrics are likewise unavailable
for a constant or singleton reference distribution. Frozen reference identity, baseline
percentile and eligibility must agree across scenarios; duplicate/missing cutoff or reference
entries, inconsistent counts and out-of-range values are rejected. Participant counts are
summarized across seeds, never pooled as if they represented distinct people.

The new output directory is never overwritten; `summary.json` marks successful completion. Preserve
the plan and original results together with the summary for inspection and reproduction.

### Public-data study pilot

The checked [study pilot record](../bench/results/study-summary-pilot.json) uses the same
full PGS000018 model and nested eight-target/64-reference public inputs described in the
[genome-wide validation](genomewide-robustness-validation.md). It tests sample-site omissions
at 5% and 20%, seeds 41/42/43, and omission versus frozen-reference filling: 14 arms including
two controls, producing 112 sample/scenario rows. Its baseline contains 1,712,283 terms.

The independent [summary checker](../bench/check_study_summary.py) recalculates all seven
scenario metrics from the exact native baseline/scenario sum strings, then recalculates
all seed distributions and paired comparisons with Python's standard library. It also derives
cutoff membership from those exact sums and checks frozen-reference shifts against individual
rows in `references.tsv`. All 98 score-metric checks, 42 score-distribution checks and 108
cutoff/reference-distribution checks pass, including 112 participant/reference rows. This
checks descriptive summaries; it does not independently reconstruct these new masked scores
or reference scores from source genotypes.

| Missingness | Mean MAE, omission | Mean MAE, filling | Mean paired MAE difference (fill minus omit) |
|---:|---:|---:|---:|
| 5% | 0.21847 | 0.03561 | -0.18286 |
| 20% | 0.34825 | 0.13162 | -0.21663 |

The added summaries illustrate why score error alone is insufficient: at 5% missingness,
filling reduces mean absolute score error in this pilot but increases the mean fraction
crossing the 90th-percentile cutoff from 0.04167 to 0.08333 (about 4.17% to 8.33%). The mean
absolute shift against the frozen 64-person reference distribution falls from 7.8125 to
1.5625 percentile points. These are different descriptive measures, not a universal policy
ranking or clinical decision thresholds.

These are observations in one small, fixed public cohort and one conditional score universe.
Three masking seeds do not support population or clinical conclusions. Unlike site-level
omission, sample-site omission can hide different terms in different people, so frequency
filling can change their relative ordering. No phenotypes enter this pilot.

After running a planned experiment and summarizing it, verify the same evidence with:

```sh
python3 bench/check_study_summary.py planned-study study-results study-summary
```

### Completed 128-person public study

The [checked larger-study record](../bench/results/public-study128.json) contains 128 targets,
64 disjoint references and the full PGS000018 input inventory across all 22 autosomes.
Omission and frequency filling are paired across 1%, 5%, 10% and 20% sample-site omissions
and seeds 11/22/33/44/55, with two controls: 42 scenarios and **5,376 participant/scenario rows**.
The frozen baseline retains 1,712,276 common-callable terms.

All 128 exact baseline scores, the baseline term digest and the source-input fingerprints
match the earlier source-verified genome-wide run. The independent summary checker passes
294 score-metric checks, 84 score-distribution checks and 216 cutoff/reference-distribution
checks, including 5,376 participant/reference rows. The new masks are not independently
reconstructed from source genotypes; this adds descriptive seed/policy validation to the
separately source-verified baseline, not phenotype validation.

| Missingness | Mean MAE, omission | Mean MAE, filling |
|---:|---:|---:|
| 1% | 0.08205 | 0.01897 |
| 5% | 0.19317 | 0.04999 |
| 10% | 0.26757 | 0.07836 |
| 20% | 0.37458 | 0.12975 |

In this fixed cohort, filling lowers mean absolute score error at each tested rate. At 20%
omission, the mean fraction crossing the frozen target 90th-percentile cutoff falls from
9.84375% to 6.09375%; the mean absolute shift against the 64-person reference distribution
falls from 17.41699 to 5.35400 percentile points. Five masking seeds characterize this
experiment's sensitivity, not population uncertainty or a general recommendation.

The native run took 2,070.04 seconds with 6,527.38 MiB peak RSS using the original frozen
release binary, before frequency sharing and random-draw caching. Development checks ran
concurrently. These measurements do not describe the optimized candidate's 128-person
performance and are not controlled scaling measurements.

### Reproduce the larger public study

First complete the input preparation and 128-person run in the
[genome-wide guide](genomewide-robustness-validation.md#reproduce). From the same checkout,
create two policies beside that run's manifest, preserving its relative input paths:

```sh
python3 - "$EVIDENCE/native/cohort128/experiment.json" <<'PYTHON'
import json
from pathlib import Path
import sys
source = Path(sys.argv[1])
experiment = json.loads(source.read_text())
policies = []
for original, name in [('no_mask', 'omit'), ('random20_fill', 'fill')]:
    policy = dict(next(s for s in experiment['scenarios'] if s['id'] == original))
    policy.update(id=name, masks=[])
    policies.append(policy)
experiment['scenarios'] = policies
with source.with_name('study-policies.json').open('x') as output:
    json.dump(experiment, output, indent=2)
PYTHON
pgsum robustness-plan --template "$EVIDENCE/native/cohort128/study-policies.json" \
  --rates 0.01,0.05,0.1,0.2 --seeds 11,22,33,44,55 --mask-unit sample-site \
  --out "$EVIDENCE/planned-study128"
pgsum --threads 4 robustness --experiment "$EVIDENCE/planned-study128/experiment.json" \
  --out "$EVIDENCE/study128-results"
pgsum robustness-summarize --plan "$EVIDENCE/planned-study128" \
  --results "$EVIDENCE/study128-results" --reference-policy omit \
  --out "$EVIDENCE/study128-summary"
python3 bench/check_study_summary.py "$EVIDENCE/planned-study128" \
  "$EVIDENCE/study128-results" "$EVIDENCE/study128-summary"
```

Use the candidate's `pgsum` executable on your `PATH`. This creates 42 scenarios, including
two controls, and retains the same 90/95 percentile cutoffs, 9p21 region and disjoint reference
set. New output directories are required. Relocating inputs changes path-bearing manifest
fingerprints; retain each run's own plan, native results and summary together.

### Memory use when reusing frequencies

Scenarios that name the same resolved frequency-table path share its decoded blocks for
the duration of the experiment. This avoids retaining one large frequency table for every
rate, seed and policy. Different paths keep separate tables and their own input fingerprints.

The [checked memory comparison](../bench/results/shared-frequency-memory.json) reruns the
14-scenario public pilot above using the original release-candidate binary and a build with
shared frequency tables. Peak RSS fell from **2,060.48 MiB to 501.14 MiB** (75.7%). All score
summaries and all four numerical output files were identical, and those files also match
the independently checked public pilot. This is one sequential observation on a shared host,
with the 128-person study running concurrently; it does not establish a universal memory
reduction or controlled runtime improvement.

To compare two native builds on your own experiment (Linux with GNU `time`, Python 3.9+):

```sh
python3 bench/compare_robustness_memory.py --before /path/to/previous/pgsum \
  --after /path/to/current/pgsum --experiment planned-study/experiment.json \
  --out memory-comparison --threads 4
```

The [comparison harness](../bench/compare_robustness_memory.py) requires identical inputs,
score summaries and numerical payload bytes before writing `comparison.json`. It retains
individual logs and peak-memory measurements. It does not rerun a source-genotype oracle.

### Reusing paired random draws

Within each genomic locus, scenarios with the same seed and mask unit share the same
random draw. Each rate still applies its own threshold, and sample-site draws remain tied
to the original sample identity. The cache is discarded as the locus changes, so its size
depends on seeds and participants rather than the number of genome-wide terms.

The [checked draw-cache comparison](../bench/results/mask-cache-comparison.json) preserves
all numerical output bytes from the public pilot. Its observed elapsed time changed from
45.81 to 35.25 seconds, with approximately 501 MiB peak RSS in both runs. This single shared-host
measurement is not a controlled speedup estimate. Twelve additional synthetic runs check both
mask units, repeated runs, one/four threads and full contribution attribution against an
independent oracle. Boundary rates, Unicode sample identities and revisited loci are covered
by the cache regression test.

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

The [dense genome-wide benchmark](genomewide-robustness-validation.md) runs the downloaded
PGS000018 inventory against public genotypes on all 22 autosomes, with nested target cohorts
and a fixed, disjoint reference group. Each experiment freezes its own common-callable
baseline across scenarios. Coverage and the baseline digest are compared across cohort sizes,
so exclusions caused by adding samples remain visible.

The [synthetic scale benchmark](../bench/README.md#synthetic-genome-wide-robustness-scaling)
separately varies model size, sample count and missingness. Its exact source-formula oracle
checks scores, masks, ranks, cutoff crossings, attribution and thread-count repeatability.
Recorded time and memory cover native commands separately from data generation and checking.

For a quick local check, the [offline example](../examples/robustness/README.md) includes
synthetic genotypes, a miniature model, a reference FASTA and hand-calculated expectations.
It requires only the `pgsum` binary and Python 3.8+; no downloads are needed. Candidate
release archives include the example and verify it after extraction.
