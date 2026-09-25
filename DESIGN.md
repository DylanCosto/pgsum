# pgsum design

Status: draft for v0. Everything here is open to change until the parity test passes.

## Goal

Calculate polygenic scores from **one sample's gVCF**, for many scores at once, fast, and with genotype
rules that treat a passing reference block as a confident homozygous-reference call rather than a missing
genotype.

## Non-goals (v0)

- No imputation. A term without a usable call is never filled in with a mean or reference dosage.
- No silent partial scores. The strict raw score exists only when every term is scorable; a partial sum is
  reported separately, labelled, with its coverage (see [Completeness](#completeness)).
- No percentiles, ancestry projection or absolute risk. Output is a raw, uncalibrated score.
- No clinical interpretation.
- GRCh38 only.

## Pipeline

```
PGS Catalog harmonized scoring files ──compile──▶ packs (binary, one per score)
                                                      │
single-sample gVCF + GRCh38 FASTA ──────────extract──▶ genotype table at the union of pack sites
                                                      │
                                         packs + genotype table ──score──▶ one result per score
```

`pgsum run` does `extract` then `score` in one command.

### 1. `compile`: scoring file → pack

Done once per PGS Catalog release, independent of any sample.

- Input: a PGS Catalog *harmonized* scoring file (`hm_chr`, `hm_pos`, `effect_allele`, `other_allele`,
  `effect_weight` or `dosage_{0,1,2}_weight`, model flags) and its header metadata.
- Each term gets a **description** (below). Terms needing review are kept in the pack with their reasons,
  never dropped, so completeness can be judged later.
- SNV orientation against the reference FASTA is resolved at compile time (it does not depend on the
  sample), so `score` only needs genotypes.
- Output: terms in source order, the source file's SHA-256, the score's metadata (weight type, licence,
  publication), and per-term review reasons. Weights are stored exactly (see [Pack format](#pack-format-pgsum-pack-v1)).
- **Packs and scoring files are never committed to this repo.** Users download scoring files themselves;
  individual scores carry their own licence terms, and the pack keeps each score's licence field.

### 2. `extract`: gVCF → genotype table

- Input: a bgzipped single-sample gVCF (no index needed), the reference FASTA and `.fai`, and packs compiled
  against that same reference (checked by digest).
- Targets: every distinct oriented SNV `(contig, pos, REF, ALT)` from a pack term with no review reasons and a
  resolved orientation. Packs are read one at a time, so memory holds one pack's terms at most.
- One pass over the file: BGZF blocks are decompressed on all cores (`noodles-bgzf` with libdeflate) and one
  thread scans the text, reading only `CHROM`, `POS`, `REF` and `INFO/END` from each record. Records on
  contigs other than chr1–22, X, Y and M are skipped. Records must be sorted and each contig contiguous.
- Every record whose interval `[POS, INFO/END or POS+len(REF)-1]` overlaps a target is kept verbatim. After
  the scan, targets are assessed in parallel from their records with the genotype rules below.
- Output: a genotype table (`pgsum-genotypes-v1`, `.pgsg`) with each target's state, ALT dosage and flags, the
  kept gVCF lines each call was made from, and a header with the gVCF, reference and pack digests, the
  sample ID, the DeepVariant version, whether the header defines `RefCall`, and counts per state.

Why one stream rather than per-contig index queries: on HG002 (437 MB, 36.0M records) decompression takes
0.3 s on 12 cores and scanning every line 0.9 s, so reading the whole file is not the bottleneck. For
comparison, full record parsing with htslib took 18.2 s.

### 3. `score`: packs × genotype table → results

Each term's outcome is `model_term_requires_review` (any review reason), `unresolved_orientation`, the call
state, or `scorable_observation` with its effect-allele dosage (the ALT dosage when the effect allele is the
ALT, otherwise 2 − ALT dosage) and exact contribution. Output per score: a JSON result
(`<pgs_id>.score.json`), optionally a per-term TSV (`--terms`), and a `scores.tsv` summary across scores.
Packs are scored in parallel, and each pack is split into fixed chunks of 200,000 terms scored in parallel
and combined in order (`rayon` work-stealing across both), so one large pack uses every core and results do
not depend on the thread count. Weights with 64-bit coefficients are summed in 128-bit integers; only wider
ones use arbitrary precision. All 8 development packs (11.7M terms) score against HG002 in 0.4 s after the
genotype table loads (0.7 s), about 29M terms/s on 12 cores.

## Term description (compile time)

A term needs review, and cannot contribute, when any of these apply:

| Reason | Condition |
|---|---|
| `requires_special_model:*` | `is_haplotype`, `is_diplotype` or `is_interaction` is true |
| `invalid_model_flag:*` | a model flag is not `TRUE`/`FALSE`/empty/`.`/`NA` |
| `conditional_inclusion_requires_review` | `inclusion_criteria` is set |
| `variant_description_requires_review` | `variant_description` is set and not a recognised form |
| `specific_imputation_method_requires_review` | `imputation_method` is set |
| `conflicting_weight_models` | dominant and recessive, or explicit dosage weights plus either flag |
| `invalid_dosage_weights` / `invalid_effect_weight` | a weight is not a finite decimal (≤128 chars, exponent −1000..100) |
| `unresolved_harmonized_position` | `hm_chr`/`hm_pos` missing or invalid |
| `harmonization_mismatch:*` | `hm_match_chr` or `hm_match_pos` is `False` |
| `author_other_allele_missing` | `other_allele` is empty |
| `nonliteral_alleles` | an allele is not `[ACGT]+` |
| `identical_effect_other_alleles` | effect allele equals other allele |

Weight model: `dosage_weights` if any `dosage_N_weight` is set, else `dominant`, `recessive` or `additive`
from the flags. Contribution for effect-allele dosage `d ∈ {0,1,2}`:

- additive: `w × d`
- dominant: `w × [d > 0]`
- recessive: `w × [d = 2]`
- dosage_weights: `dosage_d_weight`

## Orientation (compile time)

- **SNVs:** compare the reference base at `hm_pos` with the effect and other alleles, directly and
  complemented. Exactly one match resolves the term; otherwise it is `palindromic_orientation_unresolved`
  (A/T, C/G) or `allele_representation_unresolved`.
- **Longer alleles (indels, MNVs):** `sequence_reference_evidence_required` in v0. The effect/other strings
  alone don't say which allele is the reference sequence.
- Frequency-based resolution of palindromic SNVs is out of scope for v0.

## Genotype rules (extract time)

Default policy `pgsum-diploid-dp10-gq20-pass-v1`. Checks run in this order; the first failure is the state.

Before the site checks, for every overlapping record:

1. `FORMAT/FT` present and not `PASS`/`.` → `genotype_filtered`.
2. Reference block whose REF base differs from the FASTA → `reference_anchor_mismatch`.
3. DeepVariant convention: when the header defines the `RefCall` filter and a record is `0/0` with
   `FILTER=RefCall` only, the filter is treated as `PASS` and the adaptation is recorded.

Site checks:

| # | State | Condition |
|---|---|---|
| 1 | `reference_mismatch` | the term's REF differs from the FASTA at the target span |
| 2 | `unknown_no_record` | no record overlaps the target |
| 3 | `ambiguous_overlapping_records` | more than one record overlaps |
| 4 | `unsupported_ploidy` | GT is not diploid |
| 5 | `no_call` / `partial_no_call` | both / one allele is `.` |
| 6 | `filtered` | site FILTER not `PASS`/`.` (or FT not passing) |
| 7 | `quality_missing` | depth or GQ missing (depth is `MIN_DP` for reference blocks, `DP` otherwise) |
| 8 | `low_quality` | depth < 10 or GQ < 20 |
| 9 | `unsupported_symbolic_genotype` | reference block with a non-`0/0` GT |
| 10 | `incomplete_reference_span` | reference block ends before the target span ends |
| 11 | `unsupported_allele_representation` | variant record's POS or REF differs from the target |
| 12 | `other_called_allele` | a called ALT is not the term's ALT |
| — | `observed_reference` / `observed_variant` | passing call; dosage = count of the term's ALT |

A record is a reference block when every ALT is `<NON_REF>`, `<*>` or `.`.

Phase (`|` with a numeric `PS`) is preserved in the output but does not affect additive scoring.

## Completeness

Two sums are reported for every score.

**Strict.** `status` is `complete_uncalibrated_score` and `raw_score` is set only when:

- every term is `scorable_observation`,
- the inventory is consistent (scoring-file header, Catalog record and actual term count agree), and
- the Catalog records that the scoring file matches its publication.

Otherwise `status` is `score_withheld`, `raw_score` is null, and `withheld_because` lists each unmet
condition.

**Partial.** `partial.raw_score` is always the exact sum over scorable terms. Unscorable terms contribute
nothing; no genotype is imputed. It carries its coverage: the share of terms, and the share of summed
absolute per-term effect (`|w|`, or `max(|w1 − w0|, |w2 − w0|)` for dosage weights) held by scorable terms.
A partial sum is only comparable with another score or a reference distribution computed over the same
terms. When the strict score exists, the two are equal.

Why both: on HG002, none of the 8 development scores meets the strict rule, although the four genome-wide
ones have 99.5–99.8% of terms scorable. The remaining terms are low-quality calls, overlapping records,
deletions spanning the site, no-calls and unresolved palindromic SNVs, which any real genome will have.

## Arithmetic

Weights are exact decimals and the sum is exact (no floating point). Contributions are grouped by exponent
in 128-bit accumulators and combined with arbitrary-precision integers at the end. The sum starts from
decimal zero, so its exponent is the smallest exponent among the scorable contributions (and at most 0),
and it is written in Python `decimal.Decimal` text form. Raw scores are therefore identical as strings to
the reference implementation's, trailing zeros included (e.g. `1332.64750200000000000000000`).

Coverage fractions are floating point; they describe the sum and never enter it.

## Output (`pgsum-score-v1`)

One JSON object per score (`src/score.rs`, `ScoreResult`):

- `pgs_id`, `sample_id`, `policy`, `status`, `raw_score` (string or null), `withheld_because`
- `partial`: `raw_score`, `coverage` (`scorable_terms`, `total_terms`, `term_fraction`, `weight_fraction`,
  `terms_without_weight`) and a note on comparability
- `required_terms`, `scorable_terms`, `states` (term status → count)
- `weight_type` (Catalog), `license`, `matches_publication`, `inventory_consistent`
- `imputation_performed: false`, `calibration: "uncalibrated: …"`
- `inputs`: SHA-256 of the pack records, the scoring file, the gVCF, the genotype table body and the FASTA

Per-term TSV (`--terms`): the `inspect` columns, then status, call state, effect dosage and contribution.

## Provenance

Every output records the digests of all inputs and the pgsum version and policy ID. The same inputs always
give byte-identical output: sorted keys, deterministic ordering, and no timestamps in the scored content.

## Parity

The v0 bar: on GIAB HG002 (GRCh38, public data), for a fixed set of scores, every term's state and dosage
and every score's value match an existing reference implementation of these rules. See
`tests/parity_hg002.rs`.

## Pack format (`pgsum-pack-v1`)

One file per score, `<pgs_id>.pgsp`: magic bytes, a JSON header, then one zstd frame of fixed-layout term
records in source order (layout in `src/pack.rs`). The header carries the scoring file's name, SHA-256 and
size, its header lines and columns, the Catalog REST record and its digest, the licence, both weight types
(scoring file and Catalog, kept side by side rather than reconciled), the reference FASTA and `.fai`
digests, the inventory (declared, Catalog and actual term counts) and per-state counts. Weights are stored
exactly: a 64-bit coefficient and exponent, or the original text when the coefficient is wider.

`pgsum inspect <pack>` prints every term as TSV; `--header` prints the header.

## Parity status

**M1 (compile), 2026-09-25:** term description, review reasons, orientation (method, REF, ALT, effect
direction) and exact weight text are identical to the reference implementation for:

- the 8 development scores, 11,689,907 terms: PGS000001, PGS000004, PGS000013, PGS000018, PGS000027,
  PGS000662, PGS000667, PGS002724;
- the synthetic `tests/fixtures/PGS999999` file (47 terms, at least one per rule), checked in CI by
  `tests/compile_fixture.rs`.

**M2 (extract), 2026-09-25:** on GIAB HG002 (DeepVariant 1.10.0 gVCF, 36.0M records), every term's status,
call state and effect-allele dosage is identical to the reference implementation for all 11,689,907 terms of
the 8 development scores, and for the synthetic `tests/fixtures/PGS999998` + `synthetic.g.vcf.gz` pair (one
scenario per genotype rule), checked in CI by `tests/extract_fixture.rs`. `extract` over all 8 packs
(6,825,089 distinct targets) takes 10.2 s on 12 cores with a 3.9 GB peak.

**M3 (score), 2026-09-25:** on HG002, every term's status, effect dosage and contribution text, and every
score's exact partial sum, are identical to the reference implementation for all 8 development scores
(checked by `tests/parity_hg002.rs` with `PGSUM_PARITY_DIR`). For the small scores and the synthetic fixtures
the strict status and raw score also match the reference implementation's own complete-score reduction,
including `PGS999997`, which is complete: `1332.64750200000000000000000`.

| Score | Terms | Scorable | Weight covered | Partial raw score |
|---|---|---|---|---|
| PGS000001 | 77 | 90.91% | 91.17% | 1.150607011830277419 |
| PGS000004 | 313 | 72.84% | 73.26% | 1.2048 |
| PGS000013 | 6,630,150 | 99.47% | 99.54% | 17.98628904251893594 |
| PGS000018 | 1,745,179 | 99.62% | 99.19% | -1.914331603599 |
| PGS000027 | 2,100,302 | 99.79% | 99.78% | 38.733425392830412602414284095 |
| PGS000662 | 269 | 81.04% | 82.06% | 20.236805567 |
| PGS000667 | 43 | 0% | 0% | 0 |
| PGS002724 | 1,213,574 | 99.77% | 99.80% | -2461760.133929721328533131 |

All are `score_withheld` under the strict rule.

PGS000018's scoring file declares 1,745,180 variants and contains 1,745,179; its pack records the inventory
as inconsistent, so its score will be withheld.

Known differences from the reference implementation, none of which occur in current Catalog files:

- Digits are ASCII only; the reference implementation's regular expressions also accept other Unicode
  digits.
- The term limit is 50 million per score (the reference implementation stops at 10 million; the largest
  Catalog score has 13.1 million).
- Scoring-file and Catalog weight types are both recorded; pgsum does not reject a disagreement.

## Open questions

- Compile memory: about 4 GB peak for the development set with four files in parallel. Records could be
  written as they are read instead of held in memory.
- chrX/chrY ploidy for male samples: v0 requires diploid calls, which withholds scores with X terms for
  those samples.
- Whether to add frequency-supported orientation for palindromic SNVs, and from which reference panel.
