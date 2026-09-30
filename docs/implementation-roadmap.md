# Scoring workflow expansion

Authorized scope: implement all six additions discussed after PR #1. This checklist is an implementation ledger, not a substitute for the requested features. An item is complete only when its behavior and validation are present.

- [ ] Input compatibility: explicit GT/DS/GP selection, exact decimal dosage arithmetic, model-specific probability expectations, single-sample and cohort persistence, native BCF, PGEN/BGEN inputs, GRCh37 and GRCh38 with assembly safeguards.
- [ ] Batch workflow: sample sheets for separate gVCFs, bounded concurrent work and memory, resumable atomic outputs, incremental samples, shared variant matching and decoding, end-to-end scale benchmarks.
- [ ] Migration: `compare` command, score and per-variant differences with policy explanations, reproducible public real-data PLINK/pgsc_calc comparisons, independent edge-case oracles.
- [ ] Unified entry point: PGS IDs to cached fetch/compile/extract/score, preflight diagnostics, versioned run manifests with hashes and policies, Bioconda packaging and container, documented installation.
- [x] Missingness: rank consequential missing terms, source/reason reporting, mathematically valid contribution bounds with assumptions, single/batch/cohort report integration.
- [ ] Optional ancestry module: reference PCA fit/projection, continuous mean and variance normalization, reusable reference artifacts, validation against established methods, separate raw-scoring path.

## Validation gates

Preserve hard-call numerical results and old table readers. Test missing, malformed, filtered and multi-allelic inputs, allele reversals, sex/ploidy, nonadditive models, decimal precision, file round trips and CLI behavior. Publish measured end-to-end time and memory rather than extrapolating from summation benchmarks. Test resume invalidation when inputs, packs, parameters or software change. Ancestry adjustments require independent numerical and reference-workflow comparisons, not only internally consistent tests.

## In progress

Starting from merged main c8acca6. The sections below record local development checkpoints; the category checklist tracks the full authorized scope.

### Quantitative input foundation (working branch)

Implemented explicit GT/DS/GP for single samples, exact expected contributions, v5 genotype tables with v1–v4 reading, and v4 score metadata. Cohort v2 storage and scoring now retain quantitative inputs, with bounded extraction blocks and single-sample parity tests. Remaining quantitative requirements include multiallelic marginalization, probability rounding policy, sex-chromosome dosage conventions, and quantitative reference panels. Current non-diploid calls and multi-allelic DS/GP records are explicitly withheld; no fallback or rounding is used. These limitations must be addressed before input compatibility is checked off.

Validation so far: 97 tests pass on stable Rust and Rust 1.88, one external HG002 test remains ignored. Formatting and clippy pass. The first-stage hard-call performance smoke (10,000 terms, 65 samples, 1/4 threads) passed independent-oracle and exact-output checks. The cohort changes also passed a fresh hard-call smoke and the extended DS/GP oracle (5,000 terms, 65 samples, 1/4 threads), including multiple quantitative blocks. This is not evidence for the real-data cross-tool benchmark or full-scale memory requirement.

Quantitative performance remains an optimization target: in the 5,000-term/65-sample smoke, DS extraction used about 93 MiB and GP extraction about 168 MiB; DS/GP scoring used roughly 74–129 MiB. Current exact-decimal storage is much heavier than packed GT storage. Improve representation and bounded scoring caches before claiming broad cohort efficiency. Results are in the local work/quantitative-oracle-smoke-verified/results.json.

### BCF and resumable batch foundation (local, unpublished)

Added streaming BCF decoding for single-sample and cohort extraction, with compressed and uncompressed parity fixtures for GT/DS/GP. Original input hashes are preserved. Indexed BCF scans and uncommon codec encodings are not yet supported; unsupported codec branches return an input error. BCF dosage precision is limited by the stored float representation.

Added sample-sheet batch execution with bounded worker concurrency, a shared target cache, per-sample atomic publication, content-hashed manifests, verified resume, incremental samples, and aggregated scores. Tests cover changed bytes with preserved timestamps, changed packs/reference, corrupt or incomplete outputs, failed samples, and output locking. A general memory budget and scale benchmarks remain outstanding.

Current local validation: 105 tests pass on stable Rust, one external HG002 test remains ignored; formatting and clippy pass. These additions are not in merged PR #1. The six scope categories above remain incomplete.

### Catalog-driven workflow (local, unpublished)

`run --ids` and `batch --ids` now fetch and compile requested Catalog scores when no explicit packs
are supplied. A reference/rules-specific cache verifies pack hashes, preserves immutable generations,
and supports explicit offline reuse and refresh. Downloads use the existing Catalog client; malformed
responses/downloads and inconsistent term inventories cannot publish a valid cache receipt.

`run --preflight` checks input headers, selected sample, available assembly evidence and pack integrity
before extraction. Successful runs publish an execution manifest with executable/input/reference/pack
identities, policies, settings and output hashes. Failed reruns remove their prior completion marker.
Documentation now covers these workflows, batch resume, native BCF and their remaining limits.

Validation: 109 tests pass on stable Rust and Rust 1.88; one external HG002 test remains ignored.
The Catalog tests use a real local HTTP server to exercise download, compile, cache refresh/integrity,
reference invalidation, offline CLI scoring and batch execution. BCF preflight is checked for compressed
and uncompressed inputs. Formatting and clippy pass. Real public-data cross-tool validation, packaging,
full scale/memory work and the other unchecked requirements are still outstanding.

### Missingness diagnostics (local development)

Single-sample scoring, including batch children, now reports exact missing-contribution envelopes,
ranked consequential terms, exclusion reasons and bounded retained-record evidence. Unsupported model
semantics/ploidy cannot receive invented bounds; unknown terms or uncertain inventory/publication
metadata withhold a whole-score envelope. The report remains independent of filling and placement.
Per-term TSVs include bounds and their unavailable reasons. Rankings retain only a fixed number of
candidates per chunk and use exact decimal comparisons with deterministic ties.

At that checkpoint, joint-cohort `.pgsc` diagnostics remained outstanding. The new benchmark verifies additive envelopes independently and checks legacy-result parity;
hand-calculated fixtures enumerate all completions across additive, dominant, recessive and arbitrary
dosage-weight models, with separate extreme-precision, karyotype and metadata/legacy-JSON checks.

Missingness validation: the full regression suite and added tests cover 115 passing tests on stable and
Rust 1.88, with the external HG002 test still ignored; formatting and clippy pass. The larger synthetic
comparison (two 500,000-term scores, five repeats, 1/4 threads) preserved all legacy result fields and
passed an independent bound oracle. With 50% missing terms, reporting added about 4–5% end-to-end time
and 8–16% scoring time; peak scoring memory was nearly unchanged. The hard-call/DS/GP regression smoke
also passed against `1521f9f`. Exact benchmark binary hashes and runs are stored locally; the reproducible
runner and measured table are in `bench/`.

### Joint-cohort diagnostics (local development)

New `.pgsc` files now use v3, retaining per-sample exclusion reasons in compact diagnostic frames and
original VCF/BCF source-record ordinals. Numeric scoring does not read those frames. The optional
`score --missingness` report uses the same bound/ranking logic as single samples and records its output
hash in cohort score metadata. Legacy v1/v2 files still score; their missing reasons are explicitly
unavailable and cannot justify complete-score bounds. Source ordinals count all original records,
including ignored contigs, while excluding headers. Model-excluded cohort positions still do not
retain source evidence; reports state that limitation explicitly.

Validation covers exact GT/DS/GP parity with independently extracted selected samples, report stability
across thread counts, genuine old-writer fixtures, chunk-spanning ordinals, VCF/BCF source references,
corrupt frame digests and inconsistent diagnostic contents. Cohort extraction now applies the existing
single-sample contig-length safeguard. The full suite passes 120 tests on stable and Rust 1.88, with
one external HG002 test ignored; formatting and clippy pass.

At that checkpoint, the eight-frame diagnostic cache was bounded, but the packed/measurement cache
still retained all loaded blocks. Optional reporting keeps a bounded candidate list for every sample and adds a
term-by-sample pass. This does not complete the separate memory-budget or scale-benchmark requirement.


The independent benchmark oracle now checks cohort bounds, ranking, source ordinals, exact numerical
outputs and thread stability. The GT comparison used two 200,000-term scores and 65 samples, five
repeats at 1/4 threads; DS/GP used 10,000 terms and 65 samples, three repeats. Ordinary GT scoring
changed roughly −2% to +3%, while retaining diagnostics added 12–18% extraction time. Optional reports
add a substantial separate pass. See `bench/README.md` for measured times/memory and exact binary
identities. This completes the scoped missingness feature, with the documented source-evidence and
legacy-file limitations; the other five categories remain open.


### Cohort memory and shared quantitative scoring (local development)

Decoded genotype/probability frames now use an evicting payload cache with a configurable target;
workers borrow rows from their current shared block. Ordinary scoring no longer preloads every
frame. Quantitative scores with block-ordered targets share a streaming pass, decoding bounded waves
in parallel and processing each score in its original term order. The CLI opens at most two to eight
quantitative packs together; unordered packs use the independent cache path. GT retains its parallel
pack schedule. Metadata exposes the chosen paths, loads, evictions and peak wave payload.

Source evidence is gathered once per retained target in key order. Previously validated diagnostic
lookups use a small genotype-key prefix rather than reparsing probabilities. Numeric reloads verify
frame digests before reusing successful probability validation; invalid field types and failed reads
cannot be marked valid. Report generation skips redundant expectation arithmetic for known passing
measurements. Tests cover live handles after eviction, concurrent readers, exact DS/GP model results,
non-monotone fallback, excluded terms before/after scored terms, parallel waves, wrong field types and
bounded report decoding over 24 blocks.

This does not complete the batch/memory category: the payload target is not a total process RAM cap,
wave sizing uses observed block sizes, pack data/accumulators and decoding buffers use extra memory,
and separate-file batch scheduling still needs a global memory policy. General shared scheduling for
unordered packs and broader real-data scaling validation also remain outstanding.


Cohort memory validation: 122 tests pass on stable Rust and Rust 1.88, with one external HG002 test
still ignored; formatting and clippy pass. Final synthetic benchmarks preserve every tested score
and report. For two 50,000-term GP scores across 65 samples, one-thread scoring used 132.5 MiB versus
731.5 MiB previously, with essentially unchanged runtime (6.585 versus 6.614 s). At four threads it
used 456.5 MiB versus 912.1 MiB, taking 3.084 versus 2.502 s. The 200,000-term GT comparison and a
16 MiB cache-target GT/DS/GP smoke also passed. Full measurement/source fingerprints are committed in
`bench/results/cohort-memory-synthetic.json`; see the adjacent README for scope and tradeoffs. The
batch/memory category remains unchecked for the outstanding global-budget and unordered-sharing work.

### Comparison and migration gate (local development)

Added `compare` for native single/batch/cohort summaries, explicit PLINK 2 sum columns and pgsc_calc
IID/PGS/SUM aggregate files. Exact arbitrary-precision decimal arithmetic covers differences and
absolute/relative tolerances. Gzip/zstd readers preserve empty TSV cells; duplicate rows, ambiguous
identities, mixed samplesets, missing values and output overwrite attempts cannot silently pass.
Cross-file family conflicts remain visible and fail the gate.

Paired same-source-ordinal term exports stream through per-term differences, observed policy/matching/
model/dosage evidence changes and independent sum reconciliation on both sides. The command retains
input hashes and output hashes in a completion marker; optional failure exits leave complete reports
available for inspection. Aggregate score tables are retained in memory. External per-variant exports
need explicit conversion; the command does not invent contributions or causes from aggregate scores.

Validation: 130 tests pass on stable Rust and Rust 1.88, one external HG002 test remains ignored;
formatting and clippy pass. The actual PLINK 2 a7.10 executable and an independent Python Decimal oracle
agree with pgsum for 128 synthetic sample/score pairs (32 samples, 48 variants, additive/dominant/
recessive and high-precision additive weights, both effect-allele directions, missing calls). 97 pairs
match exactly; 31 match within the explicitly configured rendering tolerance. The reproducible runner
and measurement identities are in `bench/compare_plink.py` and `bench/results/compare-plink-synthetic.json`.
Public real-data PLINK/pgsc_calc execution comparisons remain outstanding; the migration category is
not checked off by these synthetic and format-fixture checks.

### Public-data migration validation (local development)

Added a pinned-source, reproducible runner on 32 real public chromosome 22 samples (all 1,059,079
records retained) and all chr22 terms of PGS001229 and PGS000018 (849 and 26,529). It executes the
unmodified pgsc_calc v2.3.0 Nextflow scoring workflow, pgsum and a separate PLINK 2 shared-term run.
Custom `CHR22_*` and `COMMON_*` model names preserve the distinction from complete published PGS values;
all original term ordinals and source/model metadata are retained. Package/tool/source identities and
per-sample evidence are recorded. Source hashes, prepared-input hashes and executable identities are
validated before reuse; all stage outputs are fresh.

The full default comparison retains 62/64 differences beyond the explicit output-rounding tolerance.
All 64 shared-term scores agree with PLINK within tolerance and with an independent Python Decimal
oracle exactly on the pgsum side. For the default policies, all 853,312 retained native contributions
are independently checked against original genotype calls, reference-based allele orientation and
published weights. All 853,280 contributions scored by both tools are equal. Each default score gap
is reconciled from eligibility choices across all 32 samples; the main difference is overlapping-record
exclusion, with one unresolved indel and one reference-effect term not matched by pgsc_calc.

The final fresh prepare/native/pgsc/analyze run passes every gate. The checked measurement record is
`bench/results/public-chr22-validation.json`; `docs/public-validation.md` explains the findings and
reproduction. This is public real-data correctness evidence for GT-only chr22 partial contributions,
not whole-genome, gVCF, DS/GP, sex-chromosome, ancestry or performance validation. Broader migration
validation and the remaining input/batch/packaging/ancestry requirements remain open. Rust scoring code
was unchanged in this checkpoint.
