//! Filling a score's missing terms from population allele frequencies (`score --fill-frequencies`).
//!
//! A term that is not scorable for the sample (no call, a failed call, another called allele, no orientation,
//! a review reason) is *missing*. With a frequency table, each missing additive term is filled with its expected
//! contribution `w × 2 × f`, where `f` is the population frequency of its effect allele, or omitted with a reason:
//!
//! * the term needs a GRCh38 position, no review reason other than an unresolved harmonized position, an
//!   additive model, and an autosomal position (frequencies of sex chromosomes depend on sex);
//! * its alleles are the resolved orientation's (reference strand); without one, a palindromic term is omitted,
//!   a term with unresolved alleles too, and any other term uses its published alleles as written;
//! * exactly one record at the position must carry the same two alleles; its frequency must be a number in
//!   0..=1. An effect allele equal to the record's ALT takes `f`; one equal to its REF takes `1 − f`, but only
//!   when the record is the position's only record and is not flagged as split from a multi-allelic site,
//!   since there `1 − f` is not the REF frequency.
//!
//! Arithmetic is exact. Fills are applied only when at least 99% of terms and 99% of summed `|effect weight|`
//! are scorable, measured before filling; the rest is reported either way.

use std::cmp::Ordering;
use std::collections::BTreeMap;

use num_bigint::BigInt;
use num_traits::{Signed, Zero};
use serde::{Deserialize, Serialize};

use crate::decimal::Decimal;
use crate::frequencies::{FrequencyTable, Record};
use crate::orient::Status;
use crate::pack::TermRecord;
use crate::score::{Coefficient, Contribution, ExactSum};
use crate::term::{AlleleKind, Model, Reason};
use crate::{Result, invalid};

/// Why a missing term was not filled.
pub const OMIT_REASONS: [&str; 11] = [
    "no_position",
    "requires_review",
    "non_additive_model",
    "no_frequency_source_for_contig",
    "palindrome_orientation_unresolved",
    "source_alleles_unresolved",
    "site_absent",
    "alleles_differ",
    "multiple_matching_records",
    "frequency_missing_or_invalid",
    "multiallelic_reference_effect_allele",
];

/// A missing term's fill.
#[derive(Clone, Debug, PartialEq)]
pub enum Outcome {
    Filled {
        /// `orientation_direct`, `orientation_complement` or `published_alleles`, then `:effect_is_alt` or
        /// `:effect_is_ref`.
        method: String,
        /// The record's frequency as written.
        frequency: String,
        contribution: Contribution,
    },
    Omitted(&'static str),
}

fn big(c: &Contribution) -> BigInt {
    let m = match &c.coefficient {
        Coefficient::Small(v) => BigInt::from(*v),
        Coefficient::Big(v) => v.clone(),
    };
    if c.negative { -m } else { m }
}

fn decimal_big(d: &Decimal) -> BigInt {
    let m: BigInt = d.coefficient.parse().expect("decimal digits");
    if d.negative { -m } else { m }
}

/// `1 − d`, with Python's exponent (the smaller of 0 and `d`'s).
fn one_minus(d: &Decimal) -> Decimal {
    let e = d.exponent.min(0);
    let value = BigInt::from(10u8).pow((-e) as u32) - decimal_big(d) * BigInt::from(10u8).pow((d.exponent - e) as u32);
    Decimal {
        negative: value.is_negative(),
        coefficient: value.abs().to_string(),
        exponent: e,
    }
}

/// `w2 × f` where `w2` is the weight times two: coefficients multiply, exponents add, signs combine.
fn times(w2: &Contribution, f: &Decimal) -> Contribution {
    let product = big(w2).abs() * decimal_big(f).abs();
    Contribution {
        negative: w2.negative != f.negative,
        coefficient: Coefficient::Big(product),
        exponent: w2.exponent + f.exponent,
    }
}

/// `0 <= d <= 1`.
fn is_probability(d: &Decimal) -> bool {
    let v = decimal_big(d);
    if v.is_negative() {
        return false;
    }
    if v.is_zero() || d.exponent >= 0 {
        return v.is_zero() || (d.exponent == 0 && v == BigInt::from(1u8));
    }
    v <= BigInt::from(10u8).pow((-d.exponent) as u32)
}

/// The fill of one missing term. `alleles` are its published effect and other allele (packs from v5).
pub fn resolve(
    term: &TermRecord,
    alleles: Option<(&str, &str)>,
    table: &FrequencyTable,
    w2: Option<Contribution>,
) -> Result<Outcome> {
    use Outcome::Omitted;
    if term.contig == 0 || term.pos == 0 {
        return Ok(Omitted("no_position"));
    }
    if term.reasons.0 & !(Reason::UnresolvedPosition as u32) != 0 {
        return Ok(Omitted("requires_review"));
    }
    if term.model != Model::Additive {
        return Ok(Omitted("non_additive_model"));
    }
    if term.contig > 22 {
        return Ok(Omitted("no_frequency_source_for_contig"));
    }
    let o = &term.orientation;
    let (effect, other, method): (String, String, String) = if o.status == Status::Resolved {
        let (r, a) = ((o.ref_base as char).to_string(), (o.alt_base as char).to_string());
        let (e, x) = if o.effect_is_alt { (a, r) } else { (r, a) };
        (e, x, format!("orientation_{}", o.method.as_str()))
    } else if term.palindromic {
        return Ok(Omitted("palindrome_orientation_unresolved"));
    } else if term.allele_kind == AlleleKind::Unresolved {
        return Ok(Omitted("source_alleles_unresolved"));
    } else {
        let Some((e, x)) = alleles else {
            return invalid!("the pack does not record published alleles; compile it with pgsum 0.5 or later");
        };
        (e.to_owned(), x.to_owned(), "published_alleles".to_owned())
    };
    resolve_alleles(term.contig, term.pos, &effect, &other, &method, table, w2)
}

/// The fill of an additive term at `contig`:`pos` with these effect and other alleles (reference strand, or as
/// published); `method` names where the alleles came from. `w2` is the term's weight times two.
pub fn resolve_alleles(
    contig: u8,
    pos: u32,
    effect: &str,
    other: &str,
    method: &str,
    table: &FrequencyTable,
    w2: Option<Contribution>,
) -> Result<Outcome> {
    use Outcome::Omitted;
    if contig == 0 || pos == 0 {
        return Ok(Omitted("no_position"));
    }
    if contig > 22 {
        return Ok(Omitted("no_frequency_source_for_contig"));
    }
    let records: Vec<&Record> = table.at(contig, pos)?;
    if let Some(r) = records.iter().find(|r| r.alt.contains(',')) {
        return invalid!(
            "the frequency table has a multi-ALT record at {}:{} ({}); split multi-allelic sites first",
            crate::term::CONTIGS[contig as usize - 1],
            pos,
            r.alt
        );
    }
    if records.is_empty() {
        return Ok(Omitted("site_absent"));
    }
    let same = |r: &&&Record| {
        let (mut a, mut b) = ([effect, other], [r.ref_allele.as_str(), r.alt.as_str()]);
        a.sort_unstable();
        b.sort_unstable();
        a == b
    };
    let matches: Vec<&&Record> = records.iter().filter(same).collect();
    let record = match matches[..] {
        [] => return Ok(Omitted("alleles_differ")),
        [one] => *one,
        _ => return Ok(Omitted("multiple_matching_records")),
    };
    let Some(f) = record
        .frequency
        .as_deref()
        .and_then(Decimal::parse)
        .filter(is_probability)
    else {
        return Ok(Omitted("frequency_missing_or_invalid"));
    };
    let Some(w2) = w2 else {
        return Ok(Omitted("non_additive_model"));
    };
    let text = record.frequency.clone().expect("parsed above");
    if effect == record.alt {
        return Ok(Outcome::Filled {
            method: format!("{method}:effect_is_alt"),
            frequency: text,
            contribution: times(&w2, &f),
        });
    }
    if record.multiallelic || records.len() > 1 {
        return Ok(Omitted("multiallelic_reference_effect_allele"));
    }
    Ok(Outcome::Filled {
        method: format!("{method}:effect_is_ref"),
        frequency: text,
        contribution: times(&w2, &one_minus(&f)),
    })
}

/// A count of terms with their summed `|effect weight|`.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct Part {
    pub terms: u64,
    pub abs_weight: String,
}

/// Totals per chunk of terms, merged in order.
#[derive(Clone, Default)]
pub struct Tally {
    pub abs_total: ExactSum,
    pub abs_scorable: ExactSum,
    /// Terms without a numeric effect weight (their `|w|` is undefined).
    pub without_effect_weight: u64,
    pub fill_sum: ExactSum,
    pub filled: (u64, ExactSum),
    pub filled_palindromic: u64,
    pub by_method: BTreeMap<String, u64>,
    pub omitted: (u64, ExactSum),
    pub by_reason: BTreeMap<&'static str, (u64, ExactSum)>,
    pub by_state: BTreeMap<&'static str, (u64, ExactSum)>,
}

impl Tally {
    pub fn merge(&mut self, other: Tally) {
        self.abs_total.merge(other.abs_total);
        self.abs_scorable.merge(other.abs_scorable);
        self.without_effect_weight += other.without_effect_weight;
        self.fill_sum.merge(other.fill_sum);
        self.filled.0 += other.filled.0;
        self.filled.1.merge(other.filled.1);
        self.filled_palindromic += other.filled_palindromic;
        for (k, v) in other.by_method {
            *self.by_method.entry(k).or_default() += v;
        }
        self.omitted.0 += other.omitted.0;
        self.omitted.1.merge(other.omitted.1);
        for (k, (n, s)) in other.by_reason {
            let e = self.by_reason.entry(k).or_default();
            e.0 += n;
            e.1.merge(s);
        }
        for (k, (n, s)) in other.by_state {
            let e = self.by_state.entry(k).or_default();
            e.0 += n;
            e.1.merge(s);
        }
    }
}

/// `sum × 100 ≥ whole × 99`, exactly.
fn at_least_99_percent(part: &Decimal, whole: &Decimal) -> bool {
    let e = part.exponent.min(whole.exponent);
    let scale = |d: &Decimal| decimal_big(d) * BigInt::from(10u8).pow((d.exponent - e) as u32);
    (scale(part) * BigInt::from(100u8)).cmp(&(scale(whole) * BigInt::from(99u8))) != Ordering::Less
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Filled {
    #[serde(flatten)]
    pub part: Part,
    pub sum: String,
    pub palindromic: u64,
    pub by_method: BTreeMap<String, u64>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Omitted {
    #[serde(flatten)]
    pub part: Part,
    pub by_reason: BTreeMap<String, Part>,
}

/// A score's fill (`ScoreResult::fill`).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Fill {
    /// The frequency table: its field and sources (with their digests).
    pub field: String,
    pub sources: Vec<crate::frequencies::Source>,
    pub terms: u64,
    pub scorable_terms: u64,
    /// Summed `|effect weight|` of all terms and of scorable ones; `None` when a term has no numeric effect
    /// weight (then the weight share, and so the fill, is undefined).
    pub abs_weight_total: Option<String>,
    pub abs_weight_scorable: Option<String>,
    /// At least 99% of terms, and of `|effect weight|`, are scorable (before filling).
    pub terms_at_least_99_percent: bool,
    pub weight_at_least_99_percent: bool,
    /// The sum of scorable contributions (the partial score).
    pub observed_sum: String,
    pub filled: Filled,
    pub omitted: Omitted,
    /// Missing terms by their status.
    pub missing_by_state: BTreeMap<String, Part>,
    /// `observed_sum + filled.sum`, when both 99% conditions hold.
    pub filled_score: Option<String>,
}

fn part((n, s): &(u64, ExactSum)) -> Part {
    Part {
        terms: *n,
        abs_weight: s.finish().to_python_string(),
    }
}

impl Fill {
    pub fn from_tally(t: &Tally, table: &FrequencyTable, terms: u64, scorable: u64, observed: &ExactSum) -> Fill {
        let defined = t.without_effect_weight == 0;
        let (total, scored) = (t.abs_total.finish(), t.abs_scorable.finish());
        let terms_ok = scorable * 100 >= terms * 99;
        let weight_ok = defined && !decimal_big(&total).is_zero() && at_least_99_percent(&scored, &total);
        let mut filled_score = observed.clone();
        filled_score.merge(t.fill_sum.clone());
        Fill {
            field: table.header.field.clone(),
            sources: table.header.sources.clone(),
            terms,
            scorable_terms: scorable,
            abs_weight_total: defined.then(|| total.to_python_string()),
            abs_weight_scorable: defined.then(|| scored.to_python_string()),
            terms_at_least_99_percent: terms_ok,
            weight_at_least_99_percent: weight_ok,
            observed_sum: observed.finish().to_python_string(),
            filled: Filled {
                part: part(&t.filled),
                sum: t.fill_sum.finish().to_python_string(),
                palindromic: t.filled_palindromic,
                by_method: t.by_method.clone(),
            },
            omitted: Omitted {
                part: part(&t.omitted),
                by_reason: t.by_reason.iter().map(|(k, v)| (k.to_string(), part(v))).collect(),
            },
            missing_by_state: t.by_state.iter().map(|(k, v)| (k.to_string(), part(v))).collect(),
            filled_score: (terms_ok && weight_ok).then(|| filled_score.finish().to_python_string()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn d(t: &str) -> Decimal {
        Decimal::parse(t).unwrap()
    }

    #[test]
    fn python_decimal_arithmetic() {
        assert_eq!(one_minus(&d("0.1233")).to_python_string(), "0.8767");
        assert_eq!(one_minus(&d("1")).to_python_string(), "0");
        assert_eq!(one_minus(&d("0")).to_python_string(), "1");
        let w2 = Contribution {
            negative: true,
            coefficient: Coefficient::Small(24),
            exponent: -2,
        };
        // Decimal("-0.12") * 2 * Decimal("0.1233")
        assert_eq!(times(&w2, &d("0.1233")).to_decimal().to_python_string(), "-0.029592");
        assert_eq!(times(&w2, &d("0")).to_decimal().to_python_string(), "-0.00");
        assert!(is_probability(&d("1.000")) && is_probability(&d("0")) && is_probability(&d("0.5")));
        assert!(!is_probability(&d("1.0001")) && !is_probability(&d("-0.1")) && !is_probability(&d("2")));
        assert!(at_least_99_percent(&d("99"), &d("100")) && !at_least_99_percent(&d("98.99"), &d("100")));
    }
}
