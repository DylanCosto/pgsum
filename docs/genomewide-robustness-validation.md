# Genome-wide public-data robustness validation

This experiment compiles the full published **PGS000018 metaGRS_CAD** scoring file and tests numerical sensitivity to omitted genotype calls on public 1000 Genomes data across all 22 autosomes. It uses nested cohorts of 8, 32 and 128 target people and a fixed, disjoint set of 64 reference people. Each cohort has its own common-callable baseline, held fixed across its five scenarios. The experiment does not use phenotypes or establish predictive accuracy, disease risk, or calibrated population percentiles.

All three cohorts passed independent verification. Measurements and source fingerprints are in the checked [validation record](../bench/results/genomewide-robustness-validation.json). The whole published model is the input inventory; the numerical scores are conditional on the common-callable subset. They must not be presented as complete published scores.

## Checked results (2026-10-13)

| Targets | References | Common-callable terms | Exact sample/scenario checks | Observed time | Peak RSS |
|---:|---:|---:|---:|---:|---:|
| 8 | 64 | 1,712,283 | 40 | 10.30 s | 512.66 MiB |
| 32 | 64 | 1,712,281 | 160 | 34.22 s | 526.69 MiB |
| 128 | 64 | 1,712,276 | 640 | 119.46 s | 585.40 MiB |

All 840 sample/scenario checks agree exactly with the independent source calculation. The 8- and 32-person runs were completed on October 11 under shared-host load; the final 128-person run was completed on October 13 after the busy application stopped. These observations are not controlled throughput comparisons.

## Sources and selection

- [PGS Catalog PGS000018](https://www.pgscatalog.org/score/PGS000018/) and its [harmonized GRCh38 scoring file](https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000018/ScoringFiles/Harmonized/PGS000018_hmPOS_GRCh38.txt.gz). The pinned file contains **1,745,179 actual rows**, while its header declares **1,745,180**. Every actual row is compiled, and the report retains this inventory discrepancy. No score rows are removed to improve coverage.
- [1000 Genomes GRCh38 phased VCFs hosted by UCSC](https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/): `ALL.chr{1..22}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz`, plus tabix indexes. The full source VCFs total 13,322,146,206 bytes. They contain 2,548 people; this validation selects 192 evenly spaced header indexes before looking at scores, using `floor(i*(2548-1)/191)` for `i=0..191`.
- Selected indexes congruent to 2 modulo 3 supply the fixed 64 references. The other 128 supply the targets; the first 8 and 32 form the smaller nested cohorts. Sample IDs and selection order are recorded in the evidence. This is a deterministic public test cohort, not a population-representative or ancestry-stratified sample.
- The reference sequence is [NCBI GRCh38 no-alt analysis set](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz), supplied as an indexed, uncompressed FASTA. Its uncompressed SHA-256 is `9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7`.

The downloader pins the score, catalog metadata, every chromosome VCF and every index. A changed upstream file fails the fingerprint check. The run records the FASTA and index fingerprints, selected source BCFs, independent oracle arrays, native genotype tables, frequency table, executable, manifests and report outputs. Genotypes and published weights are not committed to the repository or included in the release archive.

## Experimental design

The native genotype extractor reads original diploid GT calls with `--accept-missing-quality`, because these public population VCFs do not supply per-person GQ/DP. No reference calls are invented for absent VCF records. For each target cohort separately, eligible terms must be autosomal, have a supported literal allele orientation and have valid calls in **every target and all 64 reference people**. The report counts all excluded terms by reason and fingerprints the eligible original source ordinals and coordinates. One original locus, chr3:16,816,464, lands on the official reference's IUPAC ambiguity code `Y`; it remains unresolved and contributes no guessed allele. A larger cohort may lose terms that the smaller cohort could retain; both the extent and digest of that intersection are reported.

Five scenarios are run against each frozen baseline:

1. No experimental omission.
2. Deterministic 20% random **site-level** omission with seed 42; a selected site is hidden in all target people.
3. The same random omissions with allele-frequency fills.
4. Dropout of chr9:21,500,001–22,500,000, represented by the half-open BED interval `9 21500000 22500000`.
5. The same regional dropout with allele-frequency fills.

These omissions are deliberately controlled sensitivity tests, not a measured consumer-array manifest or a model of clinical missingness. Site-level random masks complement the separate synthetic benchmark's sample-site masks. The regional scenario tests a contiguous locus loss. All frequency estimates use only the 64 named reference people; none of the target calls enters those estimates. Fill remains subject to allele matching and conservative eligibility. Fills below production scoring coverage gates are hypothetical robustness reconstructions.

Because a site-level mask removes the same terms from every target, its fixed-frequency fill adds the same score offset to every person. It can change raw-score error and frozen-reference percentiles, but cannot change the perturbed within-cohort ordering. The filled and unfilled versions therefore have identical Spearman correlations in this experiment.

Target scores are compared both with their frozen target-cohort baseline distribution and with the **unmasked, frozen 64-reference score distribution**. Separately reranked target percentiles are labeled as such. Cutoffs of 90 and 95 are descriptive numerical thresholds, not clinical decisions. Reference raw sums are recalculated for each cohort's changed term universe, then remain fixed across that cohort's scenarios.

## Independent verification

[`bench/robustness_genomewide.py`](../bench/robustness_genomewide.py) performs a sorted overlap sweep of retained source VCF records and source score positions. It derives allele orientation directly from the reference FASTA and authored effect/other alleles, then obtains each person's dosage from the original GTs. This is bounded by chromosome-sized arrays and avoids rescanning all source records for every score term. The oracle never reads native packs, native genotype caches, or native contribution output as scoring truth.

The pinned score weights are integral after multiplying by `10^12`. Complete 64-person reference calls give 128 alleles per eligible site; consequently, a filled expected diploid contribution is exact at a common integer scale of `64 × 10^12`. The preparer checks an absolute-sum overflow bound over the entire original inventory before using signed 64-bit array arithmetic. Native decimal strings are checked using exact rational arithmetic, including a denominator gate that rejects extra tiny tails. The source-specific parser fails closed unless calls are diploid with single-digit allele indexes, source sample order matches, and retained records are biallelic.

The verifier checks every sample/scenario baseline sum, perturbed sum and delta exactly. It also checks all retained/masked/filled/unfilled/changed counts; all ten largest changes and source-ordinal tie ordering; signed regional changes; top-plus-residual reconciliation; baseline term counts and source-order fingerprints; fixed and reranked percentiles; reference percentiles; correlations, error summaries and cutoff transitions. Mask counts include terms whose contribution does not change because their dosage or weight is zero. The dense run does not request the potentially very large per-term TSV; exact independent aggregate and top-term checks still cover every eligible contribution.

The final verifier hash identifies the code used for validation. Preparation ran while assertion and orchestration checks were being completed; its original input/output fingerprints and timing records are retained, rather than attributing every earlier preparation subprocess to that final script hash.

## Reproduce

This harness is tested on Linux and requires GNU `time` at `/usr/bin/time`; its machine metadata reads Linux `/proc` files. These harness instructions do not claim macOS or Windows support.

From a source checkout, use Python 3.12 (the tested version) and the small pinned dependency file. The script itself requires Python 3.11 or newer, but other interpreter/dependency combinations were not tested:

```sh
python3.12 -m venv .venv-genomewide
.venv-genomewide/bin/pip install -r bench/robustness_genomewide.requirements.txt
cargo build --release --locked
```

Provide the indexed GRCh38 FASTA above as `REFERENCE`. The source download requires about 13.4 GB; allow additional space for prepared source subsets, oracle arrays and native caches. All evidence stays outside the checkout in this example.

```sh
PY=.venv-genomewide/bin/python
EVIDENCE=../pgsum-genomewide-evidence
REFERENCE=/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
$PY bench/robustness_genomewide.py fetch --out "$EVIDENCE/sources"
$PY bench/robustness_genomewide.py prepare --sources "$EVIDENCE/sources" \
  --reference "$REFERENCE" --out "$EVIDENCE/prepared"
$PY bench/robustness_genomewide.py native_inputs --prepared "$EVIDENCE/prepared" \
  --pgsum target/release/pgsum --out "$EVIDENCE/native"
for N in 8 32 128; do
  $PY bench/robustness_genomewide.py native --prepared "$EVIDENCE/prepared" \
    --pgsum target/release/pgsum --out "$EVIDENCE/native" --cohort-size "$N"
  $PY bench/robustness_genomewide.py verify --prepared "$EVIDENCE/prepared" \
    --out "$EVIDENCE/native" --cohort-size "$N"
done
$PY bench/robustness_genomewide.py record --prepared "$EVIDENCE/prepared" \
  --out "$EVIDENCE/native" --result "$EVIDENCE/validation.json"
```

Each `native/cohortN/robustness/` directory contains the HTML report, sample scores, largest changes, regional deltas, reference comparisons, experiment identity and output fingerprints. Preparation and native input building have separate measurements from the robustness command. The native commands request four workers. The 8- and 32-person elapsed times are shared-host observations with background application activity, without cold-cache control; they do not establish throughput or comparisons with another tool. An earlier 128-person attempt was interrupted and cancelled; its partial outputs are excluded. The checked 128-person result comes from a fresh, unthrottled run using the same executable and preserved inputs. Input construction overlapped automated verification; the robustness commands began after the separate synthetic benchmark finished.
