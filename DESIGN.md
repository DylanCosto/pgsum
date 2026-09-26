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
  resolved orientation. Packs are read in parallel, decoding only the columns that identify targets (not the
  weights), and the union is deduplicated whenever it doubles.
- Target index (`--targets-cache`, `.pgst`): the union of targets with the identity (PGS ID and records
  digest) of every pack it came from and the reference digest. A later run reads only the packs' headers and
  reuses the index if all three match; otherwise it collects targets again and rewrites the index.
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

## Custom scores

Any score, not only the Catalog's, can be compiled from a tab-separated file (plain or gzipped) laid out like
a Catalog scoring file, with the author's own GRCh38 positions:

```
#pgs_id=MY_SCORE
#genome_build=GRCh38
#weight_type=beta              optional, as are #pgs_name, #license and #variants_number
chr_name	chr_position	effect_allele	other_allele	effect_weight
chr1	1005806	T	C	0.0112
2	21263900	A	G	-0.0231
```

- `#pgs_id=` is required: 1–64 letters, digits, `_`, `.` or `-`, starting with a letter or digit, and not a
  Catalog ID (`PGS` + six digits), so a custom score can never be taken for a Catalog one.
- `#genome_build=GRCh38` (or `hg38`) is required. Other builds are rejected; lift positions over first.
- Required columns: `chr_name`, `chr_position`, `effect_allele`, and `effect_weight` or
  `dosage_0_weight`…`dosage_2_weight`. Optional: `other_allele`, `is_dominant`, `is_recessive`,
  `is_haplotype`, `is_diplotype`, `is_interaction`, `inclusion_criteria`, `variant_description`,
  `imputation_method`, `hm_inferOtherAllele`, with their Catalog meanings. `chr_name` may be `1` or `chr1`
  (`X`, `Y`, `MT`/`M` likewise).
- Every term rule is the same as for Catalog files. A `<pgs_id>.metadata.json` next to the file is used if
  present; otherwise the pack's metadata is built from the header. The strict score does not require a
  Catalog publication match, and `#variants_number`, if given, must equal the number of terms.

A file is treated as a Catalog file when its header has `#HmPOS_build=` or a Catalog `#pgs_id=`, and as
custom otherwise.

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

**Sex chromosomes.** DeepVariant writes a male's chrX and chrY as diploid calls (`1/1` for a hemizygous ALT,
`0/0` for the reference; HG002 has 122,899 `1/1` and 4,846,880 `0/0` block calls outside the X
pseudoautosomal regions), so these are scored like autosomes: a hemizygous ALT counts as two copies
("dosage compensation"). Callers that write haploid calls (`1`) are rejected as `unsupported_ploidy` by
default; `extract --haploid-xy-as-homozygous` reads a haploid chrX/chrY call as homozygous, and the genotype
table's policy ID becomes `pgsum-dp10-gq20-pass-haploid-xy-homozygous-v1`. pgsum does not infer sex.

## Inferred other alleles (opt-in)

2,214 of the 6,990 Catalog scores (September 2026) give only an effect allele. Their terms carry the review
reason `author_other_allele_missing` and are unscorable by default, as in the reference implementation.
`score --allow-inferred-other-allele` scores such a term (when that is its only review reason) with an
orientation its pack inferred at compile time, in this order:

1. **Reference convention of the score.** Among the score's terms without an author other allele whose
   effect allele is one base at a position with an A, C, G or T reference base: if at least 99% have an
   effect allele different from the reference, the effect allele is the ALT and the other allele the
   reference (`reference_anchored_effect_is_alt`); if at least 99% have the reference as effect allele, the
   effect dosage is the number of reference alleles, so any non-reference allele counts against it
   (`reference_anchored_effect_is_ref`). Terms not following the convention fall through to 2. The
   score-level threshold also guards strand: a score coded on the opposite strand has an effect allele equal
   to the reference at palindromic sites and cannot reach 99%.
2. **The Catalog's `hm_inferOtherAllele`**, when it names exactly one base and that pair orients like an
   author pair (not palindromic, exactly one strand matches the reference)
   (`catalog_inferred_other_allele`).

`extract` always includes inferred targets, so one genotype table serves both modes. Results report
`inferred_other_allele`: whether it was allowed, the scorable terms per method, and the pack's convention;
`--terms` adds the method per term. Packs record the evidence (`inference` in the header: eligible terms,
counts on each side of the convention, the convention chosen, terms per method).

Evidence and effect (September 2026 Catalog, GIAB HG002):

- **The rule's premise, checked on author data.** In 52 scores that do publish author alleles (8 development
  + 44 sampled), 8,970,169 SNV terms have an effect allele that is not the reference; in all but 8 of them
  the author's other allele is the reference.
- **What the packs inferred.** 2,160 packs have terms without an author other allele (1.34B terms). 1,539
  scores follow the effect-is-ALT convention and 13 the effect-is-REF one; 608 follow neither (effect allele
  equal to the reference about half the time, typical of risk-allele coding) and rely on the Catalog's
  allele. Inferred orientations cover 94.6% of those terms: 73.0% effect-is-ALT, 20.0% Catalog-inferred,
  1.6% effect-is-REF.
- **HG002 across all 6,990 scores:**

  | | Default | `--allow-inferred-other-allele` |
  |---|---|---|
  | Complete (strict) | 189 | 289 |
  | Term coverage ≥ 99% | 1,973 | 3,512 |
  | Term coverage ≥ 90% | 2,946 | 4,502 |
  | No scorable term | 2,214 | 114 |

  Scores without inferred terms give byte-identical results in both modes.

## Inferred palindromes (opt-in)

Palindromic SNVs (A/T, C/G) are 3.0% of Catalog terms and unresolvable from the alleles alone: both strands
fit the reference. They block completeness for 2,591 scores. Most Catalog scores are reported on the forward
strand (across 40 checked scores, 0.05% of resolved SNVs needed the complement), but not all: small scores
in particular mix strands.

Compile records each score's strand evidence (resolved non-palindromic SNVs on the forward strand and on the
complement) in the pack header (`palindromes`). When at least 100 such SNVs exist and at least 99.9% are on
the forward strand, each palindromic SNV gets its forward-strand reading as an inferred orientation
(`strand_consistent_palindrome`). `score --allow-inferred-palindromes` uses it; results report
`inferred_palindromes` (allowed, scorable terms, whether the score met the rule).

Validation against 1000 Genomes phase 3 (GRCh38 lift-over), independent of pgsum: in 30 scores that publish
`allelefrequency_effect`, at palindromic SNVs where both the author frequency and the 1000 Genomes frequency
are outside 0.35–0.65, the forward-strand reading puts the effect allele on the same side of 0.5 as the
author in 99.90% of 238,518 sites across the 19 scores that meet the rule, and in 88.45% of 1,403 sites in
the 11 that do not (three of those are at 50–60%, i.e. genuinely mixed strands).

## Inferred indels (opt-in)

Indels and multi-base terms are 0.74% of Catalog terms, in 1,738 scores. Their alleles do not say which is
the reference: with anchored notation (`AT`/`A`) the shorter allele is a prefix of the longer, so "insertion"
and "deletion" readings both fit whenever the next reference base matches. In 12 indel-rich scores (392,866
indel terms), 79% fitted both readings, 14% neither, and only 7% exactly one.

Compile orients such a term (inference kind in the pack, variant in a sparse `sequences` column; pack v4):

1. `reference_fit_sequence`: exactly one reading (effect allele as REF, or other allele as REF) matches the
   reference at `hm_pos`;
2. `public_sequence_pair`: otherwise, exactly one record in a public variant set (`--public-variants`) at that
   position has exactly the term's two alleles; its REF is the reference allele. This is the reference
   implementation's "exact public sequence pair" rule.

The chosen `(pos, REF, ALT)` is normalized (left-aligned and trimmed, `src/alleles.rs`). `extract` finds every
record overlapping the target's span; variant records match when one of their ALT alleles normalizes to the
same triple (so `GTT→GT` matches `GT→G`); a REF call must span the target; several passing `0/0` blocks that
cover the span without a gap count as homozygous reference. SNV rules are unchanged. `score
--allow-inferred-indels` uses these orientations; results report `inferred_indels` by method.

Public variant set used here: 1000 Genomes phase 3 on GRCh38 (a CrossMap lift-over),
FILTER=PASS non-SNV records, 3,245,341 records (`CHROM POS REF ALT`, 18 MB gzipped).

Evidence:

- Of the 12 scores' indel terms, public pairs orient 74.1% (and 1000 Genomes decides 273,950 of the "both
  fit" cases: effect allele is REF in 28%, other allele in 72%). Where both methods apply they never
  disagree (17,132 terms).
- On HG002, 78.0% of public-pair terms and 92.5% of reference-fit terms are scorable; unmatched allele
  representations are about 0%. The rest: overlapping records in repeats (13.4% of public-pair terms), a
  different called allele (5%), no-calls and low quality.

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

## Pack format (`pgsum-pack-v3`, v4 with indel orientations)

One file per score, `<pgs_id>.pgsp`: magic bytes, a JSON header, then one zstd frame (level 3) holding the
terms in source order as columns (layout in `src/pack.rs`): contig, delta-encoded position, model, allele
kind, flags, orientation, REF/ALT, review reasons, and weights split into tag, coefficient, exponent and text
columns, plus four columns for inferred orientations (v3; v2 packs, without them, are still read). The header
carries the scoring file's name, SHA-256 and size, its header lines and columns, the
Catalog REST record and its digest, the licence, both weight types (scoring file and Catalog, kept side by
side rather than reconciled), the reference FASTA and `.fai` digests, the inventory (declared, Catalog and
actual term counts) and per-state counts. Weights are stored exactly: a 64-bit coefficient and exponent, or
the original text when the coefficient is wider.

Why columns: on PGS000013 (6.6M terms) the v1 row layout took 45 MB and columns take 32.6 MB. What remains is
mostly the weights' significant digits, which do not compress; zstd level 19 saves only another 4% at a
hundredth of the speed. The 8 development packs take 64 MB (5.5 bytes per term), which puts the whole
Catalog (4.50 billion terms, September 2026) at about 25 GB.

`pgsum inspect <pack>` prints every term as TSV; `--header` prints the header.

### Fetching the Catalog

`pgsum fetch --ids …` or `--all` reads each score's Catalog REST record, downloads its harmonized GRCh38
file, compiles it, and deletes the download (`--keep-downloads` keeps it). A score whose pack exists with an
identical Catalog record and reference is skipped, so an interrupted run resumes. Transient HTTP failures
(network, 429, 5xx) are retried with exponential backoff. `--jobs` sets how many scores are in flight
(default 4) and `--max-variants` skips very large scores. Every score's outcome goes to `<out>/fetch.tsv`.

## Parity status

**M1 (compile), 2026-09-25:** term description, review reasons, orientation (method, REF, ALT, effect
direction) and exact weight text are identical to the reference implementation for:

- the 8 development scores, 11,689,907 terms: PGS000001, PGS000004, PGS000013, PGS000018, PGS000027,
  PGS000662, PGS000667, PGS002724;
- the synthetic `tests/fixtures/PGS999999` file (47 terms, at least one per rule), checked in CI by
  `tests/compile_fixture.rs`;
- 44 more Catalog scores fetched with `pgsum fetch` (4.2M terms), chosen for variety: every score with
  interaction terms and up to three of each common weight type (beta, OR, HR, log2 OR, dosage, unweighted,
  MetaPRS and others).

**M2 (extract), 2026-09-25:** on GIAB HG002 (DeepVariant 1.10.0 gVCF, 36.0M records), every term's status,
call state and effect-allele dosage is identical to the reference implementation for all 11,689,907 terms of
the 8 development scores, and for the synthetic `tests/fixtures/PGS999998` + `synthetic.g.vcf.gz` pair (one
scenario per genotype rule), checked in CI by `tests/extract_fixture.rs`. `extract` over all 8 packs
(6,825,089 distinct targets) takes 10.2 s on 12 cores with a 3.9 GB peak.

**M4 (the whole Catalog), 2026-09-25:** `pgsum fetch --all` compiled 6,990 of the 6,991 Catalog scores (4.50
billion terms) into 31.6 GB of packs in about 1 h 50 min, limited by download speed from EBI. PGS005164 fails
because its scoring file has a quoted allele containing a line break, which splits one row in two.
Compile output is identical to the reference implementation for the 52 scores checked (8 development + 44
sampled).

HG002 against all 6,990 packs (12-core Mac, packs on a USB SSD):

| Step | Time | Peak memory |
|---|---|---|
| `extract`, first run (collects targets from every pack, writes the index) | 120 s | 6.6 GB |
| `extract` with the target index | 17.1 s | 5.4 GB |
| `score`, all 6,990 packs (4.50B terms, about 48M terms/s) | 95 s | 8.1 GB |

After the 2,160 packs with inferred orientations were recompiled, targets rose to 39.3M; rebuilding the
target index took 137 s, and scoring took 163–171 s with the packs no longer in the file cache (35 GB of
packs, 26 GB of RAM), so full-Catalog scoring is bound by reading packs from the drive.

`extract` found 30.5M distinct targets and kept 9.8M gVCF records. Reading packs one at a time, the first
run took 378 s (362 s of it collecting targets); the table is byte-identical either way.

Under the strict rule 189 scores are complete for HG002. 2,214 have no scorable term at all: their scoring
files give only an effect allele, so every term is `author_other_allele_missing` (see Open questions).

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
- A repeated descriptive header key (e.g. two `#citation=` lines in PGS003767–PGS003769) is accepted: every
  line is kept in the pack's `scoring_file_header` and the key is listed in `duplicate_header_keys`. The
  reference implementation rejects any repeated key; pgsum still rejects repeats of `format_version`,
  `pgs_id`, `variants_number`, `weight_type`, `genome_build` and `HmPOS_build`.

## Open questions

- Compile memory: the resident size peaks near 3.8 GB for the development set, most of it pages of the
  memory-mapped reference touched during orientation (reclaimable file cache). Terms are held as columns
  while a file is compiled.
- chrX/chrY ploidy for male samples: v0 requires diploid calls, which withholds scores with X terms for
  those samples.
- Whether to add frequency-supported orientation for palindromic SNVs, and from which reference panel.
