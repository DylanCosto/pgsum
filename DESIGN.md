# pgsum design

Status: draft for v0. Everything here is open to change until the parity test passes.

## Goal

Calculate polygenic scores from **one sample's gVCF**, for many scores at once, fast, and with genotype
rules that treat a passing reference block as a confident homozygous-reference call rather than a missing
genotype.

## Non-goals (v0)

- No imputation. A term without a usable call is never filled in with a mean or reference dosage.
- No partial scores by default. A score is emitted only when every term in the model is scorable (see
  [Completeness](#completeness)).
- No percentiles, ancestry projection or absolute risk. Output is a raw, uncalibrated score.
- No clinical interpretation.
- GRCh38 only.

## Pipeline

```
PGS Catalog harmonized scoring files ──compile──▶ packs (sorted, binary, one per score)
                                                      │
single-sample gVCF + index + GRCh38 FASTA ──extract──▶ genotype table at the union of pack sites
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
- Output: terms sorted by `(contig, pos)`, the source file's SHA-256, the score's metadata (weight type,
  licence, publication), and per-term review reasons. The weights stay as the exact decimal strings from the
  source file.
- **Packs and scoring files are never committed to this repo.** Users download scoring files themselves;
  individual scores carry their own licence terms, and the pack keeps each score's licence field.

### 2. `extract`: gVCF → genotype table

- Input: bgzipped single-sample gVCF with a tabix or CSI index, the reference FASTA and `.fai`, and the
  union of target sites across the packs being scored.
- Parallel by region: each thread seeks its own contig (or contig slice) through the index. No single
  process reads the whole file.
- For each target site, keep every record whose interval `[POS, INFO/END or POS+len(REF)-1]` overlaps it.
- Output: one row per target site with the assessed call (below), plus the overlapping records' digests.

### 3. `score`: packs × genotype table → results

Exact decimal sum of each scorable term's contribution, per score.

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

Per score, the result records the count of terms in each state. The raw score is emitted only when:

- every term in the source file was processed,
- every term is `scorable_observation`, and
- the pack's inventory matches the source file.

Otherwise the status is `score_withheld` with the state counts, so the caller can see exactly why. A
partial-score mode may be added later as an explicit opt-in that labels its output as partial.

## Arithmetic

Weights are exact decimals and the sum is exact (arbitrary precision, no floating point). The raw score is
written as a decimal string. Parity checks compare numeric value, not string formatting.

## Output (`pgsum-score-v1`)

One JSON object per score:

- `pgs_id`, `status`, `raw_score` (string or null), `weight_type`
- `processed_terms`, `required_terms`, `scorable_terms`, `states` (state → count)
- `policy`, `imputation_performed: false`, `partial_score_emitted: false`
- `inputs`: SHA-256 of the gVCF, its index, the FASTA, the `.fai` and the pack; the pgsum version

Optional per-term TSV (`--terms`): ordinal, contig, pos, effect/other allele, state, effect dosage,
contribution.

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
- gVCF reading: `noodles` (pure Rust) or `rust-htslib`? noodles avoids a C dependency; htslib is faster
  on BGZF.
- chrX/chrY ploidy for male samples: v0 requires diploid calls, which withholds scores with X terms for
  those samples.
- Whether to add frequency-supported orientation for palindromic SNVs, and from which reference panel.
