//! Score packs against a genotype table.
//!
//! Each term's outcome follows `DESIGN.md`: `model_term_requires_review`, `unresolved_orientation`, the call
//! state, or `scorable_observation` with an effect-allele dosage and an exact contribution. Sums are exact
//! and formatted like Python's `decimal.Decimal`.
//!
//! Two sums are reported. The strict raw score exists only when every term is scorable, the inventory is
//! consistent and the Catalog says the scoring file matches its publication. The partial raw score sums
//! the scorable terms only; unusable terms contribute nothing and nothing is imputed.

use std::collections::BTreeMap;
use std::io::Write;

use num_bigint::BigInt;
use num_traits::{Signed, Zero};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::decimal::Decimal;
use crate::genotypes::{Entry, GenotypeTable, target_key};
use crate::orient::Status;
use crate::pack::{Inference, Oriented, Pack, TermRecord, Waivers, Weight};
use crate::term::{CONTIGS, Model, Reason};
use crate::{Result, invalid};

pub const SCHEMA: &str = "pgsum-score-v3";

/// A contribution's coefficient: a weight's 64-bit coefficient times a dosage multiplier fits in 128 bits;
/// wider weights use arbitrary precision.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Coefficient {
    Small(u128),
    Big(BigInt),
}

/// An exact contribution `(-1)^negative × coefficient × 10^exponent`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Contribution {
    pub negative: bool,
    pub coefficient: Coefficient,
    pub exponent: i64,
}

impl Contribution {
    pub fn to_decimal(&self) -> Decimal {
        Decimal {
            negative: self.negative,
            coefficient: match &self.coefficient {
                Coefficient::Small(c) => c.to_string(),
                Coefficient::Big(c) => c.to_string(),
            },
            exponent: self.exponent,
        }
    }

    /// The nearest `f64` (for placing scores in a distribution, never for the reported sums).
    pub fn to_f64(&self) -> f64 {
        let magnitude = match &self.coefficient {
            Coefficient::Small(c) => *c as f64,
            Coefficient::Big(c) => c.to_string().parse::<f64>().unwrap_or(f64::NAN),
        };
        let v = magnitude * 10f64.powi(self.exponent as i32);
        if self.negative { -v } else { v }
    }

    /// The same amount with the opposite sign.
    pub fn negated(&self) -> Contribution {
        Contribution {
            negative: !self.negative,
            ..self.clone()
        }
    }

    /// A cheap exact difference when both coefficients fit and use the same exponent.
    /// Keeping that exponent (even for zero) preserves the reported decimal precision.
    pub(crate) fn small_difference(&self, base: &Contribution) -> Option<Contribution> {
        if self.exponent != base.exponent {
            return None;
        }
        let signed = |c: &Contribution| match c.coefficient {
            Coefficient::Small(v) => i128::try_from(v).ok().map(|v| if c.negative { -v } else { v }),
            Coefficient::Big(_) => None,
        };
        let difference = signed(self)?.checked_sub(signed(base)?)?;
        Some(Contribution {
            negative: difference < 0,
            coefficient: Coefficient::Small(difference.checked_abs()? as u128),
            exponent: self.exponent,
        })
    }

    fn signed_big(&self) -> BigInt {
        let magnitude = match &self.coefficient {
            Coefficient::Small(c) => BigInt::from(*c),
            Coefficient::Big(c) => c.clone(),
        };
        if self.negative { -magnitude } else { magnitude }
    }
}

/// One term's outcome against a genotype table.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Outcome {
    pub status: &'static str,
    pub call_state: &'static str,
    pub effect_dosage: Option<u8>,
    pub contribution: Option<Contribution>,
    /// Set when the term was oriented with an inferred other allele.
    pub inferred: Option<Inference>,
    /// Set when the term was scored despite an informational `variant_description`.
    pub informational_description_accepted: bool,
}

impl Outcome {
    fn without_call(status: &'static str) -> Outcome {
        Outcome {
            status,
            call_state: "",
            effect_dosage: None,
            contribution: None,
            inferred: None,
            informational_description_accepted: false,
        }
    }
}

/// The weight applied at an effect-allele dosage: `w × d`, `w × [d > 0]`, `w × [d = 2]` or `dosage_d_weight`.
pub(crate) fn contribution(model: Model, weights: &[Weight], dosage: u8) -> Option<Contribution> {
    let (weight, multiplier) = match model {
        Model::Additive => (weights.first()?, dosage as u32),
        Model::Dominant => (weights.first()?, (dosage > 0) as u32),
        Model::Recessive => (weights.first()?, (dosage == 2) as u32),
        Model::DosageWeights => (weights.get(dosage as usize)?, 1),
    };
    match weight {
        Weight::Invalid => None,
        Weight::Small {
            negative,
            coefficient,
            exponent,
        } => Some(Contribution {
            negative: *negative,
            coefficient: Coefficient::Small(*coefficient as u128 * multiplier as u128),
            exponent: *exponent as i64,
        }),
        Weight::Text(text) => {
            let d = Decimal::parse(text)?;
            Some(Contribution {
                negative: d.negative,
                coefficient: Coefficient::Big(d.coefficient.parse::<BigInt>().ok()? * multiplier),
                exponent: d.exponent,
            })
        }
    }
}

/// Scoring choices beyond the default rules.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Options {
    /// Score terms whose only review reason is a missing author `other_allele` with the orientation their
    /// pack inferred for them (see `pack::Inference`).
    pub allow_inferred_other_allele: bool,
    /// Score terms whose `variant_description` is informational (fine-mapping statistics, the author's variant
    /// IDs, known notes) despite it; see `term::informational_description`.
    pub accept_informational_descriptions: bool,
    /// Score palindromic SNVs on the forward strand when the score's other SNVs are (almost) all there (see
    /// `pack::Inference::StrandConsistentPalindrome`).
    pub allow_inferred_palindromes: bool,
    /// Score indels and multi-base terms with the orientation their pack inferred (reference fit, or a public
    /// variant set); see `pack::Inference::ReferenceFitSequence` and `PublicSequencePair`.
    pub allow_inferred_indels: bool,
}

impl Options {
    fn waivers(&self) -> Waivers {
        Waivers {
            informational_descriptions: self.accept_informational_descriptions,
            inferred_other_allele: self.allow_inferred_other_allele,
            inferred_palindromes: self.allow_inferred_palindromes,
            inferred_indels: self.allow_inferred_indels,
        }
    }
}

/// Where a term's genotype comes from, before any sample is looked at.
pub(crate) enum Plan {
    /// Not scorable for any sample: `model_term_requires_review` or `unresolved_orientation`.
    Unscorable(&'static str),
    Snv {
        key: u64,
        effect_is_alt: bool,
    },
    Sequence {
        variant: crate::alleles::Variant,
        effect_is_alt: bool,
    },
}

/// A term's plan, the inference its orientation used, and whether an informational description was accepted.
pub(crate) fn plan(t: &TermRecord, options: &Options) -> (Plan, Option<Inference>, bool) {
    let Some((oriented, inferred)) = t.waived(options.waivers()) else {
        return (Plan::Unscorable("model_term_requires_review"), None, false);
    };
    let informational = t.reasons.0 & Reason::VariantDescription as u32 != 0;
    let plan = match oriented {
        Oriented::Snv(o) if o.status != Status::Resolved => Plan::Unscorable("unresolved_orientation"),
        Oriented::Snv(o) => Plan::Snv {
            key: target_key(t.contig, t.pos, o.ref_base, o.alt_base),
            effect_is_alt: o.effect_is_alt,
        },
        Oriented::Sequence { variant, effect_is_alt } => Plan::Sequence { variant, effect_is_alt },
    };
    (plan, inferred, informational)
}

/// A term's contribution at an effect-allele dosage, and its absolute effect (see `Coverage`).
pub(crate) fn term_contribution(t: &TermRecord, dosage: u8) -> Option<Contribution> {
    contribution(t.model, &t.weights, dosage)
}

pub(crate) fn term_effect(t: &TermRecord) -> Option<f64> {
    effect_size(t)
}

/// Classify one term against the genotype table.
pub fn outcome(t: &TermRecord, genotypes: &GenotypeTable, options: &Options) -> Result<Outcome> {
    let (plan, inferred, informational) = plan(t, options);
    let mut outcome = match plan {
        Plan::Unscorable(status) => return Ok(Outcome::without_call(status)),
        Plan::Snv { key, effect_is_alt } => called(t, genotypes.get(key)?, effect_is_alt, inferred)?,
        Plan::Sequence { variant, effect_is_alt } => {
            called(t, genotypes.get_sequence(t.contig, &variant)?, effect_is_alt, inferred)?
        }
    };
    outcome.informational_description_accepted = informational;
    Ok(outcome)
}

/// The outcome of a term from its target's entry in the genotype table.
fn called(t: &TermRecord, entry: Option<&Entry>, effect_is_alt: bool, inferred: Option<Inference>) -> Result<Outcome> {
    let Some(entry) = entry else {
        return invalid!(
            "the genotype table has no call for {}:{}; extract with this pack",
            CONTIGS[t.contig as usize - 1],
            t.pos
        );
    };
    let state = entry.state.as_str();
    let outcome = |status, effect_dosage, contribution| Outcome {
        status,
        call_state: state,
        effect_dosage,
        contribution,
        inferred,
        informational_description_accepted: false,
    };
    match entry.alt_dosage {
        Some(alt) if entry.state.is_passing() => {
            let dosage = if effect_is_alt { alt } else { 2 - alt };
            match contribution(t.model, &t.weights, dosage) {
                Some(c) => Ok(outcome("scorable_observation", Some(dosage), Some(c))),
                None => Ok(outcome("model_term_requires_review", None, None)),
            }
        }
        _ => Ok(outcome(state, None, None)),
    }
}

/// Overflowed amounts are kept at exponent −2000, below any accepted weight exponent (≥ −1127).
const FLOOR: i64 = -2000;

/// Exact running sum. Contributions are grouped by exponent in 128-bit accumulators and combined once at the
/// end; the result's exponent is the smallest exponent added (and at most 0), as for `Decimal(0) + …`.
#[derive(Clone, Default)]
pub struct ExactSum {
    buckets: Vec<(i64, i128)>,
    overflow: BigInt,
    min_exponent: i64,
}

impl ExactSum {
    pub fn add(&mut self, c: &Contribution) {
        self.min_exponent = self.min_exponent.min(c.exponent);
        let small = match c.coefficient {
            // At most (2^64 - 1) × 2, so it fits in i128 either sign.
            Coefficient::Small(v) => Some(if c.negative { -(v as i128) } else { v as i128 }),
            Coefficient::Big(_) => None,
        };
        match small {
            Some(v) => self.add_small(c.exponent, v),
            None => self.add_big(c.exponent, c.signed_big()),
        }
    }

    fn add_small(&mut self, exponent: i64, value: i128) {
        let slot = match self.buckets.iter().position(|(e, _)| *e == exponent) {
            Some(i) => i,
            None => {
                self.buckets.push((exponent, 0));
                self.buckets.len() - 1
            }
        };
        match self.buckets[slot].1.checked_add(value) {
            Some(total) => self.buckets[slot].1 = total,
            None => self.add_big(exponent, BigInt::from(value)),
        }
    }

    fn add_big(&mut self, exponent: i64, value: BigInt) {
        self.overflow += value * BigInt::from(10u8).pow((exponent - FLOOR) as u32);
    }

    /// Add another sum into this one (for sums computed in parallel).
    pub fn merge(&mut self, other: ExactSum) {
        self.min_exponent = self.min_exponent.min(other.min_exponent);
        for (e, v) in other.buckets {
            self.add_small(e, v);
        }
        self.overflow += other.overflow;
    }

    pub fn finish(&self) -> Decimal {
        let floor = FLOOR;
        let ten = BigInt::from(10u8);
        let mut total = self.overflow.clone();
        for &(e, acc) in &self.buckets {
            total += BigInt::from(acc) * ten.pow((e - floor) as u32);
        }
        // Rescale from the floor to the result exponent; the division is exact.
        total /= ten.pow((self.min_exponent - floor) as u32);
        Decimal {
            negative: total.is_negative(),
            coefficient: if total.is_zero() {
                "0".into()
            } else {
                total.abs().to_string()
            },
            exponent: self.min_exponent,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Coverage {
    pub scorable_terms: u64,
    pub total_terms: u64,
    /// `scorable_terms / total_terms`.
    pub term_fraction: f64,
    /// Share of the summed absolute per-term effect held by scorable terms. The effect of a term is `|w|` for
    /// additive, dominant and recessive terms and `max(|w1 − w0|, |w2 − w0|)` for dosage weights; terms without
    /// a valid weight are left out of both sides and counted in `terms_without_weight`.
    pub weight_fraction: f64,
    pub terms_without_weight: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Partial {
    pub status: String,
    pub raw_score: String,
    pub coverage: Coverage,
    /// Whether coverage reaches `COVERAGE_GUIDELINE` of both terms and weight, the point from which pgsum
    /// suggests using the partial score (see README, "Which number to use").
    pub meets_coverage_guideline: bool,
    pub note: String,
}

/// Share of terms and of weight a partial score should cover before it is used (the 99% rule).
pub const COVERAGE_GUIDELINE: f64 = 0.99;

/// A score's terms on the sex chromosomes and how pgsum counted them.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SexChromosomes {
    pub chrx_terms: u64,
    pub chrx_scorable: u64,
    pub chry_terms: u64,
    pub chry_scorable: u64,
    /// The karyotype and chrX coding the score was read with (`score --karyotype`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub karyotype: Option<Sex>,
    /// Present when the score has chrX or chrY terms.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Inputs {
    pub pack_records_sha256: String,
    pub scoring_file_sha256: String,
    pub gvcf_sha256: String,
    pub genotypes_body_sha256: String,
    pub reference_fasta_sha256: String,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ScoreResult {
    pub schema: String,
    pub pgsum_version: String,
    pub pgs_id: String,
    pub sample_id: String,
    pub policy: String,
    /// `complete_uncalibrated_score` or `score_withheld`.
    pub status: String,
    /// Present only with `complete_uncalibrated_score`.
    pub raw_score: Option<String>,
    pub withheld_because: Vec<String>,
    pub partial: Partial,
    pub required_terms: u64,
    pub scorable_terms: u64,
    pub states: BTreeMap<String, u64>,
    pub weight_type: Option<String>,
    pub license: Option<String>,
    pub matches_publication: Option<bool>,
    pub inventory_consistent: bool,
    pub formula: String,
    /// Any expected contribution was computed for sample filling or used in a reference comparison.
    /// The observed strict and partial scores never include imputation; see `imputation` for scope.
    pub imputation_performed: bool,
    pub imputation: Imputation,
    pub calibration: String,
    /// Whether scoring was allowed to use inferred other alleles, and how many scorable terms did.
    pub inferred_other_allele: InferredUse,
    /// Whether informational `variant_description`s were accepted, and how many scorable terms had one.
    pub informational_descriptions: InformationalUse,
    /// Whether palindromic SNVs could be read on the forward strand, and how many scorable terms were.
    pub inferred_palindromes: PalindromeUse,
    /// Whether indels and multi-base terms could be scored with inferred orientations, and how many were, by
    /// method.
    pub inferred_indels: IndelUse,
    /// Terms on chrX and chrY, with a note on the dosage convention when there are any.
    pub sex_chromosomes: SexChromosomes,
    /// The score placed among a reference panel's scores over the same terms (`score --reference`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference: Option<crate::panel::Placement>,
    /// The panel group nearest to the sample (`score --reference`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ancestry: Option<crate::panel::Ancestry>,
    /// The filled score placed within one group of a PLINK reference panel (`score --reference-group`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub placement: Option<crate::placement::GroupPlacement>,
    /// Missing terms filled from population allele frequencies (`score --fill-frequencies`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fill: Option<crate::fill::Fill>,
    pub inputs: Inputs,
}

/// Which outputs include expected contributions rather than only observed genotypes.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct Imputation {
    /// Always false: `raw_score` and `partial.raw_score` sum observed calls only.
    pub observed_scores: bool,
    /// `fill` contains at least one frequency-based contribution, even if its final score is withheld.
    pub sample_fill: bool,
    /// A reference comparison uses expected contributions (missing panel calls or shared filled terms).
    pub reference_panel: bool,
}

impl ScoreResult {
    /// Refresh the summary after attaching `reference` or `placement` results.
    /// Counts, rather than sums, detect fills even when their contributions are zero or cancel out.
    pub fn refresh_imputation(&mut self) {
        let sample_fill = self.fill.as_ref().is_some_and(|f| f.filled.part.terms > 0);
        let reference_panel = self.reference.as_ref().is_some_and(|p| p.fills.genotypes_filled > 0)
            || self
                .placement
                .as_ref()
                .is_some_and(|p| p.missing_genotypes.total > 0 || p.terms.absent_filled > 0 || sample_fill);
        self.imputation = Imputation {
            observed_scores: false,
            sample_fill,
            reference_panel,
        };
        self.imputation_performed = sample_fill || reference_panel;
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct IndelUse {
    pub allowed: bool,
    pub scorable_terms: BTreeMap<String, u64>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct PalindromeUse {
    pub allowed: bool,
    pub scorable_terms: u64,
    /// Whether the pack's strand evidence met the rule (see `pack::PalindromeSummary`).
    pub strand_consistent: Option<bool>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct InformationalUse {
    pub accepted: bool,
    pub scorable_terms: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct InferredUse {
    pub allowed: bool,
    /// Scorable terms oriented with an inferred other allele, by method.
    pub scorable_terms: BTreeMap<String, u64>,
    /// The pack's reference convention (`effect_is_alt`, `effect_is_ref`), if any.
    pub reference_convention: Option<String>,
}

/// A weight as floating point, for coverage only.
fn approximate(w: &Weight) -> Option<f64> {
    match w {
        Weight::Invalid => None,
        Weight::Small {
            negative,
            coefficient,
            exponent,
        } => {
            let v = *coefficient as f64 * 10f64.powi(*exponent);
            Some(if *negative { -v } else { v })
        }
        Weight::Text(_) => w.to_decimal()?.to_python_string().parse().ok(),
    }
}

fn effect_size(t: &TermRecord) -> Option<f64> {
    match t.model {
        Model::DosageWeights => {
            let v: Vec<f64> = t.weights.iter().map(approximate).collect::<Option<_>>()?;
            Some((v[1] - v[0]).abs().max((v[2] - v[0]).abs()))
        }
        _ => approximate(t.weights.first()?).map(f64::abs),
    }
}

/// Terms per parallel chunk. Fixed, so results do not depend on the thread count.
pub const CHUNK_TERMS: usize = 200_000;

/// The totals of one chunk of terms.
#[derive(Default)]
struct Tally {
    sum: ExactSum,
    states: BTreeMap<&'static str, u64>,
    total: u64,
    scorable: u64,
    without_weight: u64,
    effect_all: f64,
    effect_scorable: f64,
    inferred: BTreeMap<&'static str, u64>,
    informational: u64,
    /// Terms on chrX and chrY, and how many of them were scorable.
    sex_chromosome_terms: [u64; 2],
    sex_chromosome_scorable: [u64; 2],
    fill: crate::fill::Tally,
    tsv: Vec<u8>,
}

fn tally(
    pack: &Pack,
    genotypes: &GenotypeTable,
    options: &Options,
    extras: &Extras,
    start: (usize, usize),
    end: usize,
    terms: bool,
) -> Result<Tally> {
    let io = |e| crate::Error::Io {
        path: "<terms>".into(),
        source: e,
    };
    let mut t = Tally::default();
    let (offset, first) = start;
    for (i, term) in pack.terms_from(offset).take(end - first).enumerate() {
        let term = term?;
        let o = by_karyotype(&term, outcome(&term, genotypes, options)?, extras.sex);
        t.total += 1;
        let sex_chromosome = match term.contig {
            23 => Some(0),
            24 => Some(1),
            _ => None,
        };
        if let Some(x) = sex_chromosome {
            t.sex_chromosome_terms[x] += 1;
            t.sex_chromosome_scorable[x] += o.contribution.is_some() as u64;
        }
        *t.states.entry(o.status).or_default() += 1;
        let effect = effect_size(&term);
        match effect {
            Some(e) => t.effect_all += e,
            None => t.without_weight += 1,
        }
        let fill = match extras.fill {
            Some(table) => fill_term(&mut t.fill, pack, table, first + i, &term, &o)?,
            None => None,
        };
        if let Some(c) = &o.contribution {
            t.scorable += 1;
            t.effect_scorable += effect.unwrap_or(0.0);
            t.sum.add(c);
            if let Some(kind) = o.inferred {
                *t.inferred.entry(kind.as_str()).or_default() += 1;
            }
            t.informational += o.informational_description_accepted as u64;
        }
        if terms {
            crate::pack::write_term_line(&mut t.tsv, first + i + 1, &term)?;
            writeln!(
                t.tsv,
                "\t{}\t{}\t{}\t{}\t{}\t{}",
                o.status,
                o.call_state,
                o.effect_dosage.map(|d| d.to_string()).unwrap_or_default(),
                o.contribution
                    .map(|c| c.to_decimal().to_python_string())
                    .unwrap_or_default(),
                o.inferred.map_or("", Inference::as_str),
                if o.informational_description_accepted {
                    "accepted"
                } else {
                    ""
                }
            )
            .map_err(io)?;
            if extras.fill.is_some() {
                t.tsv.pop();
                let (outcome, frequency, value) = match &fill {
                    Some(crate::fill::Outcome::Filled {
                        method,
                        frequency,
                        contribution,
                    }) => (
                        format!("filled:{method}"),
                        frequency.as_str(),
                        contribution.to_decimal().to_python_string(),
                    ),
                    Some(crate::fill::Outcome::Omitted(reason)) => (format!("omitted:{reason}"), "", String::new()),
                    None => (String::new(), "", String::new()),
                };
                writeln!(t.tsv, "\t{outcome}\t{frequency}\t{value}").map_err(io)?;
            }
        }
    }
    Ok(t)
}

/// Inputs beyond the pack and genotypes that add to a score's result.
#[derive(Clone, Copy, Default)]
pub struct Extras<'a> {
    /// Fill missing terms from this frequency table (see `fill`).
    pub fill: Option<&'a crate::frequencies::FrequencyTable>,
    /// Score chrX by the sample's karyotype (see `Sex`).
    pub sex: Option<Sex>,
}

/// A sample's karyotype and how the score's authors coded chrX, for `score --karyotype`.
///
/// For an XX sample chrX terms stay diploid. For an XY sample, a chrX term outside the pseudoautosomal regions
/// is hemizygous: the gVCF's diploid call `0/0` or `1/1` is read as 0 or 1 copy, the contribution is
/// `w × k × copies` with `k` = 2 for authors who coded males 0/2 (dosage compensation) and 1 for 0/1, and a
/// heterozygous call there is missing (`heterozygous_call_in_hemizygous_region`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Sex {
    pub xy: bool,
    pub k: u8,
}

/// GRCh38 pseudoautosomal regions on chrX: `start < pos <= end`.
pub const PAR: [(u32, u32); 2] = [(10000, 2781479), (155701382, 156030895)];

pub fn in_par(pos: u32) -> bool {
    PAR.iter().any(|&(start, end)| start < pos && pos <= end)
}

/// A term's outcome for a sample of this karyotype (see `Sex`).
fn by_karyotype(term: &TermRecord, o: Outcome, sex: Option<Sex>) -> Outcome {
    let Some(sex) = sex else { return o };
    if !sex.xy || term.contig != 23 || in_par(term.pos) {
        return o;
    }
    match o.effect_dosage {
        Some(1) if o.contribution.is_some() => Outcome {
            status: "heterozygous_call_in_hemizygous_region",
            effect_dosage: None,
            contribution: None,
            ..o
        },
        Some(d) if o.contribution.is_some() => Outcome {
            contribution: contribution(term.model, &term.weights, sex.k * (d / 2)),
            ..o
        },
        _ => o,
    }
}

/// A term's `|effect weight|` and, if it is missing, its fill.
fn fill_term(
    t: &mut crate::fill::Tally,
    pack: &Pack,
    table: &crate::frequencies::FrequencyTable,
    index: usize,
    term: &TermRecord,
    o: &Outcome,
) -> Result<Option<crate::fill::Outcome>> {
    let abs = match (term.model, term.weights.first()) {
        (Model::DosageWeights, _) | (_, None) => None,
        (_, Some(w)) => {
            contribution(Model::Additive, std::slice::from_ref(w), 1).map(|c| Contribution { negative: false, ..c })
        }
    };
    match &abs {
        Some(a) => t.abs_total.add(a),
        None => t.without_effect_weight += 1,
    }
    let add = |slot: &mut (u64, ExactSum)| {
        slot.0 += 1;
        if let Some(a) = &abs {
            slot.1.add(a);
        }
    };
    if o.contribution.is_some() {
        if let Some(a) = &abs {
            t.abs_scorable.add(a);
        }
        return Ok(None);
    }
    add(t.by_state.entry(o.status).or_default());
    let w2 = contribution(Model::Additive, &term.weights, 2);
    let outcome = crate::fill::resolve(term, pack.alleles(index), table, w2)?;
    match &outcome {
        crate::fill::Outcome::Filled {
            method, contribution, ..
        } => {
            add(&mut t.filled);
            t.fill_sum.add(contribution);
            t.filled_palindromic += term.palindromic as u64;
            *t.by_method.entry(method.clone()).or_default() += 1;
        }
        crate::fill::Outcome::Omitted(reason) => {
            add(&mut t.omitted);
            add(t.by_reason.entry(reason).or_default());
        }
    }
    Ok(Some(outcome))
}

/// Score one pack. With `terms`, also write one TSV line per term.
pub fn score(
    pack: &Pack,
    genotypes: &GenotypeTable,
    options: &Options,
    terms: Option<&mut dyn Write>,
) -> Result<ScoreResult> {
    score_with(pack, genotypes, options, &Extras::default(), terms)
}

/// `score` with extra inputs (see `Extras`).
pub fn score_with(
    pack: &Pack,
    genotypes: &GenotypeTable,
    options: &Options,
    extras: &Extras,
    mut terms: Option<&mut dyn Write>,
) -> Result<ScoreResult> {
    let h = &pack.header;
    if !genotypes
        .header
        .packs
        .iter()
        .any(|p| p.pgs_id == h.pgs_id && p.records_sha256 == h.records_sha256)
    {
        return invalid!("the genotype table was not extracted with this {} pack", h.pgs_id);
    }
    let io = |e| crate::Error::Io {
        path: "<terms>".into(),
        source: e,
    };
    if let Some(out) = terms.as_deref_mut() {
        writeln!(
            out,
            "{}\tstatus\tcall_state\teffect_dosage\tcontribution\tinferred_other_allele\tinformational_description{}",
            crate::pack::TSV_HEADER.trim_end(),
            if extras.fill.is_some() {
                "\tfill\tfill_frequency\tfill_contribution"
            } else {
                ""
            }
        )
        .map_err(io)?;
    }
    // Chunks are scored in parallel and combined in order.
    let starts = pack.chunk_starts(CHUNK_TERMS)?;
    let total_terms = pack.header.inventory.actual_terms as usize;
    let want_tsv = terms.is_some();
    let tallies: Vec<Result<Tally>> = starts
        .par_iter()
        .enumerate()
        .map(|(i, &start)| {
            let end = starts.get(i + 1).map_or(total_terms, |s| s.1);
            tally(pack, genotypes, options, extras, start, end, want_tsv)
        })
        .collect();
    let mut sum = ExactSum::default();
    let mut states: BTreeMap<String, u64> = BTreeMap::new();
    let (mut total, mut scorable, mut without_weight) = (0u64, 0u64, 0u64);
    let (mut effect_all, mut effect_scorable) = (0f64, 0f64);
    let mut inferred: BTreeMap<String, u64> = BTreeMap::new();
    let mut informational = 0u64;
    let (mut sex_terms, mut sex_scorable) = ([0u64; 2], [0u64; 2]);
    let mut fill = crate::fill::Tally::default();
    for t in tallies {
        let t = t?;
        fill.merge(t.fill);
        informational += t.informational;
        for x in 0..2 {
            sex_terms[x] += t.sex_chromosome_terms[x];
            sex_scorable[x] += t.sex_chromosome_scorable[x];
        }
        for (k, v) in t.inferred {
            *inferred.entry(k.to_owned()).or_default() += v;
        }
        sum.merge(t.sum);
        for (k, v) in t.states {
            *states.entry(k.to_owned()).or_default() += v;
        }
        total += t.total;
        scorable += t.scorable;
        without_weight += t.without_weight;
        effect_all += t.effect_all;
        effect_scorable += t.effect_scorable;
        if let Some(out) = terms.as_deref_mut() {
            out.write_all(&t.tsv).map_err(io)?;
        }
    }
    if total as usize != total_terms {
        return invalid!("{}: pack holds {total} terms, header says {total_terms}", h.pgs_id);
    }
    let mut withheld = Vec::new();
    if scorable != total {
        withheld.push(format!("{} of {} terms are not scorable", total - scorable, total));
    }
    if !h.inventory.consistent {
        withheld.push(format!(
            "inventory inconsistent (scoring file declares {}, Catalog {}, file has {})",
            h.inventory.declared_terms.map_or("none".into(), |n| n.to_string()),
            h.inventory.catalog_terms.map_or("unknown".into(), |n| n.to_string()),
            h.inventory.actual_terms
        ));
    }
    if h.origin == crate::scoring_file::Origin::PgsCatalog && h.matches_publication != Some(true) {
        withheld.push("the Catalog does not record that the scoring file matches its publication".into());
    }
    let raw = sum.finish().to_python_string();
    let fill = extras
        .fill
        .map(|table| crate::fill::Fill::from_tally(&fill, table, total, scorable, &sum));
    let palindrome_terms = inferred
        .remove(Inference::StrandConsistentPalindrome.as_str())
        .unwrap_or(0);
    let indel_terms: BTreeMap<String, u64> = [Inference::ReferenceFitSequence, Inference::PublicSequencePair]
        .into_iter()
        .filter_map(|k| inferred.remove(k.as_str()).map(|n| (k.as_str().to_owned(), n)))
        .collect();
    let term_fraction = if total == 0 {
        0.0
    } else {
        scorable as f64 / total as f64
    };
    let weight_fraction = if effect_all == 0.0 {
        0.0
    } else {
        effect_scorable / effect_all
    };
    let mut result = ScoreResult {
        schema: SCHEMA.into(),
        pgsum_version: env!("CARGO_PKG_VERSION").into(),
        pgs_id: h.pgs_id.clone(),
        sample_id: genotypes.header.sample.sample_id.clone(),
        policy: genotypes.header.policy.id.clone(),
        status: if withheld.is_empty() {
            "complete_uncalibrated_score"
        } else {
            "score_withheld"
        }
        .into(),
        raw_score: withheld.is_empty().then(|| raw.clone()),
        withheld_because: withheld,
        partial: Partial {
            status: "partial_uncalibrated_score".into(),
            raw_score: raw,
            coverage: Coverage {
                scorable_terms: scorable,
                total_terms: total,
                term_fraction,
                weight_fraction,
                terms_without_weight: without_weight,
            },
            meets_coverage_guideline: term_fraction >= COVERAGE_GUIDELINE && weight_fraction >= COVERAGE_GUIDELINE,
            note: "Sum of scorable terms only. Unscorable terms contribute nothing; no genotype is imputed. Not \
                   comparable with complete scores or reference distributions unless their coverage matches."
                .into(),
        },
        required_terms: total,
        scorable_terms: scorable,
        states,
        weight_type: h.weight_type.catalog.clone(),
        license: h.license.clone(),
        matches_publication: h.matches_publication,
        inventory_consistent: h.inventory.consistent,
        formula: "Sum of published per-term contributions at the observed effect-allele dosage".into(),
        imputation_performed: false,
        imputation: Imputation::default(),
        calibration: "uncalibrated: raw scores and optional reference-panel percentiles; no absolute risk".into(),
        inferred_other_allele: InferredUse {
            allowed: options.allow_inferred_other_allele,
            scorable_terms: inferred,
            reference_convention: h.inference.as_ref().and_then(|i| i.reference_convention.clone()),
        },
        inferred_indels: IndelUse {
            allowed: options.allow_inferred_indels,
            scorable_terms: indel_terms,
        },
        inferred_palindromes: PalindromeUse {
            allowed: options.allow_inferred_palindromes,
            scorable_terms: palindrome_terms,
            strand_consistent: h.palindromes.as_ref().map(|p| p.applied),
        },
        informational_descriptions: InformationalUse {
            accepted: options.accept_informational_descriptions,
            scorable_terms: informational,
        },
        sex_chromosomes: SexChromosomes {
            chrx_terms: sex_terms[0],
            chrx_scorable: sex_scorable[0],
            chry_terms: sex_terms[1],
            chry_scorable: sex_scorable[1],
            karyotype: extras.sex,
            note: (sex_terms != [0, 0]).then(|| match extras.sex {
                None => "chrX/chrY terms are counted as 0, 1 or 2 copies, as the gVCF's diploid calls give them: a \
                         male's hemizygous ALT counts as 2. Authors differ in whether they coded males 0/1 or 0/2 on \
                         chrX, so compare with the score's publication before comparing male and female samples."
                    .into(),
                Some(Sex { xy: false, .. }) => "chrX terms are diploid for this XX sample.".into(),
                Some(Sex { xy: true, k }) => format!(
                    "chrX terms outside the pseudoautosomal regions are hemizygous for this XY sample: a 0/0 or 1/1 \
                     call is 0 or 1 copy, weighted x{k}; a heterozygous call there is missing."
                ),
            }),
        },
        reference: None,
        ancestry: None,
        placement: None,
        fill,
        inputs: Inputs {
            pack_records_sha256: h.records_sha256.clone(),
            scoring_file_sha256: h.source.sha256.clone(),
            gvcf_sha256: genotypes.header.gvcf.sha256.clone(),
            genotypes_body_sha256: genotypes.header.body_sha256.clone(),
            reference_fasta_sha256: genotypes.header.reference.fasta_sha256.clone(),
        },
    };
    result.refresh_imputation();
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(text: &str, multiplier: u32) -> Contribution {
        let d = Decimal::parse(text).unwrap();
        // As in packs: coefficients that fit in 64 bits take the fast path, wider ones arbitrary precision.
        let coefficient = match d.coefficient.parse::<u64>() {
            Ok(v) => Coefficient::Small(v as u128 * multiplier as u128),
            Err(_) => Coefficient::Big(d.coefficient.parse::<BigInt>().unwrap() * multiplier),
        };
        Contribution {
            negative: d.negative,
            coefficient,
            exponent: d.exponent,
        }
    }

    fn sum(parts: &[(&str, u32)]) -> String {
        let mut s = ExactSum::default();
        for (t, m) in parts {
            s.add(&c(t, *m));
        }
        s.finish().to_python_string()
    }

    #[test]
    fn sums_match_python_decimal() {
        // Expected values are str(Decimal(0) + Decimal(a) * m + ...).
        assert_eq!(sum(&[]), "0");
        assert_eq!(sum(&[("-0.5", 0)]), "0.0");
        assert_eq!(sum(&[("0.5", 1), ("-0.5", 1)]), "0.0");
        assert_eq!(sum(&[("1E+2", 1)]), "100");
        assert_eq!(sum(&[("1E+2", 2)]), "200");
        assert_eq!(sum(&[("0E-17", 1)]), "0E-17");
        assert_eq!(sum(&[("1.5", 1), ("-2.25E-3", 1)]), "1.49775");
        assert_eq!(
            sum(&[("0.12345678901234567890123", 2), ("1", 1)]),
            "1.24691357802469135780246"
        );
        assert_eq!(sum(&[("1e-1000", 1), ("1e100", 1)]).len(), "1".len() + 1100 + 1);
    }

    #[test]
    fn merged_sums_equal_one_sum() {
        let parts = [
            ("0.1", 2u32),
            ("-2.5E-3", 1),
            ("1234", 1),
            ("0.12345678901234567890123", 2),
            ("1E+2", 1),
            ("-0.75", 0),
        ];
        let mut whole = ExactSum::default();
        let (mut a, mut b) = (ExactSum::default(), ExactSum::default());
        for (i, (t, m)) in parts.iter().enumerate() {
            whole.add(&c(t, *m));
            if i % 2 == 0 { a.add(&c(t, *m)) } else { b.add(&c(t, *m)) }
        }
        a.merge(b);
        assert_eq!(a.finish(), whole.finish());
        let mut overflowing = ExactSum::default();
        for _ in 0..4 {
            overflowing.add_small(0, i128::MAX / 3);
        }
        assert_eq!(
            overflowing.finish().coefficient,
            (BigInt::from(i128::MAX / 3) * BigInt::from(4u8)).to_string()
        );
    }

    #[test]
    fn xy_chrx_terms_are_hemizygous_outside_the_pars() {
        use crate::orient::{Method, Orientation, Status};
        let term = |pos| TermRecord {
            contig: 23,
            pos,
            model: Model::Additive,
            allele_kind: crate::term::AlleleKind::LiteralSnv,
            palindromic: false,
            orientation: Orientation {
                status: Status::Resolved,
                method: Method::Direct,
                ref_base: b'A',
                alt_base: b'G',
                effect_is_alt: true,
            },
            reasons: crate::term::Reasons(0),
            weights: vec![Weight::Small {
                negative: false,
                coefficient: 3,
                exponent: -1,
            }],
            inferred: None,
            informational_description: false,
            inferred_sequence: None,
        };
        let called = |t: &TermRecord, d: u8| Outcome {
            status: "scorable_observation",
            call_state: "observed_variant",
            effect_dosage: Some(d),
            contribution: contribution(t.model, &t.weights, d),
            inferred: None,
            informational_description_accepted: false,
        };
        let text = |o: Outcome| o.contribution.map(|c| c.to_decimal().to_python_string());
        let (x, par) = (term(5_000_000), term(20_000));
        let xy = |k| Some(Sex { xy: true, k });
        // 1/1 is one copy: w × k × 1.
        assert_eq!(text(by_karyotype(&x, called(&x, 2), xy(2))).as_deref(), Some("0.6"));
        assert_eq!(text(by_karyotype(&x, called(&x, 2), xy(1))).as_deref(), Some("0.3"));
        assert_eq!(text(by_karyotype(&x, called(&x, 0), xy(2))).as_deref(), Some("0.0"));
        // A heterozygous call in a hemizygous region is missing.
        let het = by_karyotype(&x, called(&x, 1), xy(2));
        assert_eq!(
            (het.status, het.contribution),
            ("heterozygous_call_in_hemizygous_region", None)
        );
        // Pseudoautosomal and XX calls are unchanged.
        assert_eq!(by_karyotype(&par, called(&par, 1), xy(2)), called(&par, 1));
        assert_eq!(
            by_karyotype(&x, called(&x, 1), Some(Sex { xy: false, k: 2 })),
            called(&x, 1)
        );
        assert!(in_par(10_001) && in_par(2_781_479) && !in_par(10_000) && !in_par(2_781_480));
    }

    #[test]
    fn contribution_signs() {
        let w = [Weight::Small {
            negative: true,
            coefficient: 5,
            exponent: -1,
        }];
        let zero = contribution(Model::Additive, &w, 0).unwrap();
        assert_eq!(zero.to_decimal().to_python_string(), "-0.0");
        assert_eq!(
            contribution(Model::Dominant, &w, 2)
                .unwrap()
                .to_decimal()
                .to_python_string(),
            "-0.5"
        );
        assert_eq!(
            contribution(Model::Recessive, &w, 1)
                .unwrap()
                .to_decimal()
                .to_python_string(),
            "-0.0"
        );
    }
}
