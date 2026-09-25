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
use serde::{Deserialize, Serialize};

use crate::decimal::Decimal;
use crate::genotypes::{GenotypeTable, target_key};
use crate::orient::Status;
use crate::pack::{Pack, TermRecord, Weight};
use crate::term::{CONTIGS, Model};
use crate::{Result, invalid};

pub const SCHEMA: &str = "pgsum-score-v1";

/// An exact contribution `(-1)^negative × coefficient × 10^exponent`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Contribution {
    pub negative: bool,
    pub coefficient: BigInt,
    pub exponent: i64,
}

impl Contribution {
    pub fn to_decimal(&self) -> Decimal {
        Decimal {
            negative: self.negative,
            coefficient: self.coefficient.to_string(),
            exponent: self.exponent,
        }
    }
}

/// One term's outcome against a genotype table.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Outcome {
    pub status: &'static str,
    pub call_state: &'static str,
    pub effect_dosage: Option<u8>,
    pub contribution: Option<Contribution>,
}

impl Outcome {
    fn without_call(status: &'static str) -> Outcome {
        Outcome {
            status,
            call_state: "",
            effect_dosage: None,
            contribution: None,
        }
    }
}

/// The weight applied at an effect-allele dosage: `w × d`, `w × [d > 0]`, `w × [d = 2]` or `dosage_d_weight`.
fn contribution(model: Model, weights: &[Weight], dosage: u8) -> Option<Contribution> {
    let (weight, multiplier) = match model {
        Model::Additive => (weights.first()?, dosage as u32),
        Model::Dominant => (weights.first()?, (dosage > 0) as u32),
        Model::Recessive => (weights.first()?, (dosage == 2) as u32),
        Model::DosageWeights => (weights.get(dosage as usize)?, 1),
    };
    let (negative, coefficient, exponent) = match weight {
        Weight::Invalid => return None,
        Weight::Small {
            negative,
            coefficient,
            exponent,
        } => (*negative, BigInt::from(*coefficient), *exponent as i64),
        Weight::Text(text) => {
            let d = Decimal::parse(text)?;
            (d.negative, d.coefficient.parse().ok()?, d.exponent)
        }
    };
    Some(Contribution {
        negative,
        coefficient: coefficient * multiplier,
        exponent,
    })
}

/// Classify one term against the genotype table.
pub fn outcome(t: &TermRecord, genotypes: &GenotypeTable) -> Result<Outcome> {
    if !t.reasons.is_empty() {
        return Ok(Outcome::without_call("model_term_requires_review"));
    }
    let o = t.orientation;
    if o.status != Status::Resolved {
        return Ok(Outcome::without_call("unresolved_orientation"));
    }
    let key = target_key(t.contig, t.pos, o.ref_base, o.alt_base);
    let Some(entry) = genotypes.get(key) else {
        return invalid!(
            "the genotype table has no call for {}:{}; extract with this pack",
            CONTIGS[t.contig as usize - 1],
            t.pos
        );
    };
    let state = entry.state.as_str();
    match entry.alt_dosage {
        Some(alt) if entry.state.is_passing() => {
            let dosage = if o.effect_is_alt { alt } else { 2 - alt };
            match contribution(t.model, &t.weights, dosage) {
                Some(c) => Ok(Outcome {
                    status: "scorable_observation",
                    call_state: state,
                    effect_dosage: Some(dosage),
                    contribution: Some(c),
                }),
                None => Ok(Outcome {
                    status: "model_term_requires_review",
                    call_state: state,
                    effect_dosage: None,
                    contribution: None,
                }),
            }
        }
        _ => Ok(Outcome {
            status: state,
            call_state: state,
            effect_dosage: None,
            contribution: None,
        }),
    }
}

/// Exact running sum. Contributions are grouped by exponent in 128-bit accumulators and combined once at the
/// end; the result's exponent is the smallest exponent added (and at most 0), as for `Decimal(0) + …`.
#[derive(Default)]
pub struct ExactSum {
    buckets: Vec<(i64, i128)>,
    overflow: BigInt,
    min_exponent: i64,
}

impl ExactSum {
    pub fn add(&mut self, c: &Contribution) {
        self.min_exponent = self.min_exponent.min(c.exponent);
        let signed = if c.negative {
            -c.coefficient.clone()
        } else {
            c.coefficient.clone()
        };
        let small: Option<i128> = (&signed).try_into().ok();
        let slot = match self.buckets.iter().position(|(e, _)| *e == c.exponent) {
            Some(i) => i,
            None => {
                self.buckets.push((c.exponent, 0));
                self.buckets.len() - 1
            }
        };
        match small.and_then(|v| self.buckets[slot].1.checked_add(v)) {
            Some(total) => self.buckets[slot].1 = total,
            None => self.overflow += signed * BigInt::from(10u8).pow((c.exponent - self.min_exponent_floor()) as u32),
        }
    }

    /// Overflowed amounts are kept at exponent −2000, below any accepted weight exponent.
    fn min_exponent_floor(&self) -> i64 {
        -2000
    }

    pub fn finish(&self) -> Decimal {
        let floor = self.min_exponent_floor();
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
    pub note: String,
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
    pub imputation_performed: bool,
    pub calibration: String,
    pub inputs: Inputs,
}

fn magnitude(w: &Weight) -> Option<f64> {
    w.to_decimal()
        .and_then(|d| d.to_python_string().parse::<f64>().ok())
        .map(f64::abs)
}

fn effect_size(t: &TermRecord) -> Option<f64> {
    match t.model {
        Model::DosageWeights => {
            let v: Vec<f64> = t
                .weights
                .iter()
                .map(|w| w.to_decimal()?.to_python_string().parse().ok())
                .collect::<Option<_>>()?;
            Some((v[1] - v[0]).abs().max((v[2] - v[0]).abs()))
        }
        _ => magnitude(t.weights.first()?),
    }
}

/// Score one pack. With `terms`, also write one TSV line per term.
pub fn score(pack: &Pack, genotypes: &GenotypeTable, mut terms: Option<&mut dyn Write>) -> Result<ScoreResult> {
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
            "{}\tstatus\tcall_state\teffect_dosage\tcontribution",
            crate::pack::TSV_HEADER.trim_end()
        )
        .map_err(io)?;
    }
    let mut sum = ExactSum::default();
    let mut states: BTreeMap<String, u64> = BTreeMap::new();
    let (mut total, mut scorable, mut without_weight) = (0u64, 0u64, 0u64);
    let (mut effect_all, mut effect_scorable) = (0f64, 0f64);
    for (i, term) in pack.terms().enumerate() {
        let t = term?;
        let o = outcome(&t, genotypes)?;
        total += 1;
        *states.entry(o.status.to_owned()).or_default() += 1;
        let effect = effect_size(&t);
        match effect {
            Some(e) => effect_all += e,
            None => without_weight += 1,
        }
        if let Some(c) = &o.contribution {
            scorable += 1;
            effect_scorable += effect.unwrap_or(0.0);
            sum.add(c);
        }
        if let Some(out) = terms.as_deref_mut() {
            crate::pack::write_term_line(out, i + 1, &t)?;
            writeln!(
                out,
                "\t{}\t{}\t{}\t{}",
                o.status,
                o.call_state,
                o.effect_dosage.map(|d| d.to_string()).unwrap_or_default(),
                o.contribution
                    .map(|c| c.to_decimal().to_python_string())
                    .unwrap_or_default()
            )
            .map_err(io)?;
        }
    }
    let mut withheld = Vec::new();
    if scorable != total {
        withheld.push(format!("{} of {} terms are not scorable", total - scorable, total));
    }
    if !h.inventory.consistent {
        withheld.push(format!(
            "inventory inconsistent (scoring file declares {}, Catalog {}, file has {})",
            h.inventory.declared_terms,
            h.inventory.catalog_terms.map_or("unknown".into(), |n| n.to_string()),
            h.inventory.actual_terms
        ));
    }
    if h.matches_publication != Some(true) {
        withheld.push("the Catalog does not record that the scoring file matches its publication".into());
    }
    let raw = sum.finish().to_python_string();
    Ok(ScoreResult {
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
                term_fraction: if total == 0 {
                    0.0
                } else {
                    scorable as f64 / total as f64
                },
                weight_fraction: if effect_all == 0.0 {
                    0.0
                } else {
                    effect_scorable / effect_all
                },
                terms_without_weight: without_weight,
            },
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
        calibration: "uncalibrated: no percentile or absolute risk".into(),
        inputs: Inputs {
            pack_records_sha256: h.records_sha256.clone(),
            scoring_file_sha256: h.source.sha256.clone(),
            gvcf_sha256: genotypes.header.gvcf.sha256.clone(),
            genotypes_body_sha256: genotypes.header.body_sha256.clone(),
            reference_fasta_sha256: genotypes.header.reference.fasta_sha256.clone(),
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(text: &str, multiplier: u32) -> Contribution {
        let d = Decimal::parse(text).unwrap();
        Contribution {
            negative: d.negative,
            coefficient: d.coefficient.parse::<BigInt>().unwrap() * multiplier,
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
