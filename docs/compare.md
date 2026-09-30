# Compare existing scores with pgsum

`pgsum compare` checks reported sums without converting averages, filling missing scores or changing
allele matching. It writes a full outer join of sample/score identities, exact decimal differences,
input/output hashes and an optional comparison of per-term evidence. Matching numbers establish
arithmetic agreement for those reported values; they do not prove matching scoring policies or inputs.

## Score tables

For a single-sample pgsum run and one PLINK 2 score:

```sh
pgsum compare \
  --left results/scores.tsv --left-sample HG002 \
  --right plink.sscore --right-format plink2 \
  --right-score-column PGS000001_SUM --right-score-id PGS000001 \
  --absolute-tolerance 0.00001 --relative-tolerance 0.000005 \
  --out comparison/ --fail-on-difference
```

The example assumes `results/scores.tsv` contains only PGS000001. A PLINK input selects **one** named
sum column; additional scores on the other side are reported as unmatched. Compare separate score
outputs or prepare an explicitly selected subset to compare one score from a larger run. Tolerances
above are examples for rounded output, not defaults or universal acceptance thresholds.

For a cohort or batch against pgsc_calc:

```sh
pgsum compare \
  --left results/cohort-scores.tsv.zst \
  --right aggregated_scores.txt.gz --right-format pgsc-calc \
  --right-sampleset target \
  --out comparison/
```

Native batch `batch-scores.tsv` is also supported. Gzip and zstd are detected by file contents; filenames
do not control decoding. External tables may be tab- or whitespace-delimited. Native tables with empty
fields must retain their tab separators. Both `--left-format` and `--right-format` accept `pgsum`,
`plink2`, or `pgsc-calc` and have corresponding sample, score-column, score-id and sampleset options.

| Format | Sample identity | Score identity | Numeric value |
|---|---|---|---|
| pgsum | `sample_id` (batch), `sample` (cohort), or explicit `--left-sample` / `--right-sample` for single-sample output | `pgs_id` | `partial_raw_score` |
| PLINK 2 `.sscore` | `IID` | explicit `--left-score-id` / `--right-score-id` | explicitly named score column ending in `_SUM`; dosage-sum and average columns are rejected |
| pgsc_calc aggregated scores | `IID` | `PGS` | `SUM` |

A sample override cannot replace an existing sample column. Multiple pgsc_calc samplesets require an
explicit selection. Duplicate sample/score rows or repeated IIDs with different FID/SID identities
within an input are errors. Conflicting FID/SID values on joined rows are reported as identity
mismatches. When one file lacks these identifiers, their consistency cannot be verified. Renaming
samples, mapping score IDs and remapping families are deliberate preprocessing steps, never guesses
based on filenames or row order. pgsc_calc versions with different column schemas need an explicit
conversion to a supported schema; the reader does not guess columns.

## Numeric agreement and output

The output directory must not already exist. A successful command writes:

- `scores.tsv`: every sample/score pair from either input, original numeric text, status,
  `delta_right_minus_left`, and all other input columns as JSON evidence.
- `terms.jsonl`: every source ordinal from paired term files, if requested.
- `comparison.json`: versioned completion marker with counts, selected columns, input hashes,
  tolerances, scope, term reconciliation and output hashes. An interrupted directory without this
  marker is incomplete. Ordinary failures remove files created by that comparison.

The default tolerances are zero. Arithmetic uses arbitrary-precision decimal integers, including the
comparison itself:

```text
abs(right - left) <= absolute_tolerance + relative_tolerance * max(abs(left), abs(right))
```

Numerically equal decimals are `exact` even if their textual precision differs. Other paired numbers
are `within_tolerance` or `different`. Differences use an exact compact coefficient/exponent notation.
Missing rows are `left_only` or `right_only`. Empty values, `.`, `NA`, `N/A`, `nan` and `NaN` are
`unavailable`; even two unavailable values do not count as numeric agreement. Infinity, malformed
numbers, negative tolerances and ambiguous headers are errors. The comparison decimal parser supports
up to 32,768 input bytes and exponents from −16,384 to +16,384 after decimal-point adjustment, covering
pgsum's generated sums without converting to floating point.

By default a completed comparison exits successfully even when values differ. `--fail-on-difference`
returns failure **after saving the report** for unmatched, unavailable, identity-conflicting or
out-of-tolerance scores. Supplied term evidence must also agree and reconcile for this gate to pass.
The report's `agrees` field describes this same scoped gate. It does not certify biological equivalence.
Score summaries are retained in memory, proportional to the number of sample/score pairs and their
metadata; term comparisons stream with only the current rows and running sums retained.

## Per-term differences

For two runs using the same original score term order:

```sh
pgsum compare \
  --left before/scores.tsv --left-sample HG002 \
  --right after/scores.tsv --right-sample HG002 \
  --left-terms before/PGS000001.terms.tsv \
  --right-terms after/PGS000001.terms.tsv \
  --term-sample HG002 --term-score PGS000001 \
  --out term-comparison/ --fail-on-difference
```

Generate these files with `run --terms` or single-sample `score --terms`. One invocation can attach
one sample/score pair's term comparison. Both selected score rows must be present. Each term export
requires the following tab-separated columns (additional evidence columns are preserved):

```text
ordinal contig pos model ref alt effect_is_alt weights status call_state effect_dosage contribution
```

This line lists column names only; actual files must use tabs. Ordinals must be positive, unique,
strictly increasing and correspond to the **same original scoring-file rows** on both sides. External
per-variant outputs require an explicit conversion to this schema. PLINK `.sscore` and pgsc_calc
aggregate files do not contain per-variant contributions, so the command cannot reconstruct them.
The comparison does not perform strand alignment, variant normalization or row-number remapping.

Every ordinal is reported, including missing rows, changed positions/alleles, changed models/weights,
call states, dosages and inclusion decisions. Explanations identify observed differences in matching,
model, eligibility/policy acceptance, or dosage evidence. These are evidence categories, not inferred
causes from a numeric discrepancy alone. Extra columns such as fill diagnostics and missingness bounds
are compared too. Textual changes in evidence remain visible even when numeric contributions agree.

The exact sum of available contributions on **each** side is checked against its selected aggregate
score using the same tolerance. This exposes incomplete exports and compensating term differences
that a total-only check would miss. Blank term contributions are omitted from these partial sums, as
in pgsum scoring; they remain unavailable in the term report. Identical unavailable terms can satisfy
the evidence gate if the partial sums reconcile, but do not establish complete score coverage.
No whole-score contribution attribution is claimed when reconciliation fails.

## Interpreting a discrepancy

pgsum's comparison deliberately selects the partial sum and keeps strict-score withholding, coverage,
status and other columns as evidence. Frequency-filled scores and ancestry-adjusted values are not
substituted. Before interpreting cross-tool agreement, align the selected scoring file/build, retained
variants, allele orientation, GT/DS/GP choice, quality filters, missing-call handling, ploidy and genetic
model. A score match alone cannot establish these conditions.

PLINK's default scoring can impute missing genotypes and reports averages separately from sums. Use
its explicit sum output and choose the missingness policy you intend to compare; the independent
synthetic check uses `no-mean-imputation`. See [PLINK scoring documentation](https://www.cog-genomics.org/plink/2.0/score).
pgsc_calc's aggregate schema defines `SUM` separately from `AVG` and ancestry adjustments; precision
of rendered output can limit a comparison. See [pgsc_calc output documentation](https://pgsc-calc.readthedocs.io/en/latest/explanation/output.html).

`bench/compare_plink.py` validates real PLINK output against pgsum and an independent Python Decimal
oracle on synthetic diploid GT data. See [the benchmark instructions](../bench/README.md#migration-comparison-check).
A separate [public chromosome 22 validation](public-validation.md) exercises the actual pgsc_calc
scoring workflow and PLINK on real genotypes, including default-policy differences and per-term
reconciliation. Whole-genome and additional format/ploidy/model comparisons remain separate gates.
