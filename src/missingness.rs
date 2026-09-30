//! Exact missing-contribution envelopes, without imputing genotypes or interpreting model-ambiguous terms.
use crate::Result;
use crate::decimal::Decimal;
use crate::genotypes::{GenotypeTable, position_key};
use crate::pack::{Pack, TermRecord};
use crate::score::{Contribution, ExactSum, Options, Outcome, Plan, Sex};
use crate::term::{CONTIGS, Reason};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::BTreeMap;

pub const TOP_TERMS: usize = 20;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Interval {
    pub lower: String,
    pub upper: String,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Term {
    pub ordinal: usize,
    pub contig: Option<String>,
    pub pos1: Option<u32>,
    pub model: String,
    pub effect_allele: Option<String>,
    pub other_allele: Option<String>,
    pub weights: Vec<Option<String>>,
    pub status: String,
    pub call_state: String,
    pub review_reasons: Vec<String>,
    pub contribution_bounds: Option<Interval>,
    pub max_absolute_contribution: Option<String>,
    pub allowed_effect_dosages: Vec<u8>,
    pub bounds_unavailable_because: Option<String>,
    pub evidence: Evidence,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Evidence {
    pub availability: String,
    pub retained_record_count: usize,
    pub records: Vec<Record>,
    pub omitted_record_count: usize,
    /// Cohort v3: original VCF/BCF record ordinals (one-based, excluding headers).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source_record_ordinals: Vec<u64>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Record {
    /// Digest of retained selected-sample VCF text, not a byte offset into the original file.
    pub sha256: String,
    pub text: String,
    pub truncated: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Report {
    pub schema: String,
    pub missing_terms: u64,
    pub bounded_terms: u64,
    pub unbounded_terms: u64,
    pub by_status: BTreeMap<String, u64>,
    /// Only the subset whose published contribution function is known.
    pub bounded_missing_contribution: Interval,
    /// Absent if any missing term is unbounded or score inventory/publication metadata is uncertain.
    pub completion_score_bounds: Option<Interval>,
    pub completion_bounds_withheld_because: Vec<String>,
    pub ranked_missing_terms: Vec<Term>,
    pub unbounded_examples: Vec<Term>,
    pub ranking: String,
    pub top_terms_limit: usize,
    pub source_input_sha256: String,
    pub assumptions: Vec<String>,
}

#[derive(Clone)]
struct Bounds {
    lower: Contribution,
    upper: Contribution,
    impact: Decimal,
    dosages: Vec<u8>,
}

fn bounds(term: &TermRecord, outcome: &Outcome, sex: Option<Sex>) -> std::result::Result<Bounds, &'static str> {
    // These reasons concern matching, not the definition of the contribution function. Other reasons
    // can encode interactions, centring, special inclusion, duplicates, or unknown effect conventions.
    let matching_reasons = Reason::UnresolvedPosition as u32
        | Reason::MismatchChr as u32
        | Reason::MismatchPos as u32
        | Reason::InvalidMatchChr as u32
        | Reason::InvalidMatchPos as u32
        | Reason::OtherAlleleMissing as u32;
    let allowed_reasons = matching_reasons
        | if term.informational_description {
            Reason::VariantDescription as u32
        } else {
            0
        };
    // Reject unknown future reason bits too, rather than treating them as model-free annotations.
    if term.reasons.0 & !allowed_reasons != 0 {
        return Err("ambiguous_or_unsupported_contribution_model");
    }
    if outcome.call_state == "missing_reason_not_retained" {
        return Err("missing_reason_not_retained");
    }
    if outcome.call_state == "unsupported_ploidy" {
        return Err("unsupported_ploidy");
    }
    let dosages = match sex {
        Some(Sex { xy: true, k }) if term.contig == 23 && term.pos > 0 && !crate::score::in_par(term.pos) => {
            if !matches!(k, 1 | 2) {
                return Err("unsupported_sex_dosage_convention");
            }
            vec![0, k]
        }
        _ => vec![0, 1, 2],
    };
    // Additive/dominant/recessive models are monotone up to the weight's sign. Their extrema
    // occur at the endpoints; avoid allocating/comparing decimal strings for intermediate dosages.
    if term.model != crate::term::Model::DosageWeights {
        let zero = crate::score::term_contribution(term, 0).ok_or("invalid_or_missing_weights")?;
        let end = crate::score::term_contribution(term, *dosages.last().expect("nonempty domain"))
            .ok_or("invalid_or_missing_weights")?;
        let impact = Decimal {
            negative: false,
            ..end.to_decimal()
        };
        let (lower, upper) = if end.negative { (end, zero) } else { (zero, end) };
        return Ok(Bounds {
            lower,
            upper,
            impact,
            dosages,
        });
    }
    let mut values = dosages
        .iter()
        .map(|d| crate::score::term_contribution(term, *d))
        .collect::<Option<Vec<_>>>()
        .ok_or("invalid_or_missing_weights")?;
    values.sort_by(|a, b| a.to_decimal().numeric_cmp(&b.to_decimal()));
    let lower = values.first().expect("at least two dosages").clone();
    let upper = values.last().expect("at least two dosages").clone();
    let absolute = |c: &Contribution| Decimal {
        negative: false,
        ..c.to_decimal()
    };
    let a = absolute(&lower);
    let b = absolute(&upper);
    let impact = if a.numeric_cmp(&b) == Ordering::Greater { a } else { b };
    Ok(Bounds {
        lower,
        upper,
        impact,
        dosages,
    })
}

#[derive(Clone)]
struct Candidate {
    ordinal: usize,
    term: TermRecord,
    outcome: Outcome,
    bounds: std::result::Result<Bounds, &'static str>,
}

fn priority(a: &Candidate, b: &Candidate) -> Ordering {
    match (&a.bounds, &b.bounds) {
        (Ok(a_bounds), Ok(b_bounds)) => a_bounds
            .impact
            .numeric_cmp(&b_bounds.impact)
            .then_with(|| b.ordinal.cmp(&a.ordinal)),
        _ => b.ordinal.cmp(&a.ordinal),
    }
}

fn keep(top: &mut Vec<Candidate>, candidate: Candidate) {
    if top.len() == TOP_TERMS && priority(&candidate, top.last().expect("nonempty")) != Ordering::Greater {
        return;
    }
    top.push(candidate);
    top.sort_by(|a, b| priority(b, a));
    top.truncate(TOP_TERMS);
}

#[derive(Default)]
pub(crate) struct Tally {
    missing: u64,
    bounded: u64,
    lower: ExactSum,
    upper: ExactSum,
    by_status: BTreeMap<String, u64>,
    ranked: Vec<Candidate>,
    unbounded: Vec<Candidate>,
}

impl Tally {
    pub fn add(&mut self, ordinal: usize, term: &TermRecord, outcome: &Outcome, sex: Option<Sex>) {
        if outcome.contribution.is_some() {
            return;
        }
        self.missing += 1;
        *self.by_status.entry(outcome.status.to_owned()).or_default() += 1;
        let bounds = bounds(term, outcome, sex);
        let top = match &bounds {
            Ok(b) => {
                self.bounded += 1;
                self.lower.add(&b.lower);
                self.upper.add(&b.upper);
                &mut self.ranked
            }
            Err(_) => &mut self.unbounded,
        };
        // Avoid cloning allele/model data for the great majority of terms below the retained cutoff.
        if top.len() == TOP_TERMS {
            let worse = match (&bounds, &top.last().expect("nonempty").bounds) {
                (Ok(a), Ok(b)) => {
                    a.impact
                        .numeric_cmp(&b.impact)
                        .then_with(|| top.last().unwrap().ordinal.cmp(&ordinal))
                        != Ordering::Greater
                }
                _ => ordinal >= top.last().unwrap().ordinal,
            };
            if worse {
                return;
            }
        }
        keep(
            top,
            Candidate {
                ordinal,
                term: term.clone(),
                outcome: outcome.clone(),
                bounds,
            },
        );
    }

    pub fn merge(&mut self, other: Self) {
        self.missing += other.missing;
        self.bounded += other.bounded;
        self.lower.merge(other.lower);
        self.upper.merge(other.upper);
        for (status, count) in other.by_status {
            *self.by_status.entry(status).or_default() += count;
        }
        for candidate in other.ranked {
            keep(&mut self.ranked, candidate);
        }
        for candidate in other.unbounded {
            keep(&mut self.unbounded, candidate);
        }
    }

    pub fn finish(self, pack: &Pack, genotypes: &GenotypeTable, options: &Options, sum: &ExactSum) -> Result<Report> {
        self.finish_with(pack, &genotypes.header.gvcf.sha256, sum, |c| {
            let evidence = single_evidence(&c, genotypes, options)?;
            describe(c, pack, evidence)
        })
    }

    fn finish_with(
        self,
        pack: &Pack,
        source_hash: &str,
        sum: &ExactSum,
        describe: impl Fn(Candidate) -> Result<Term>,
    ) -> Result<Report> {
        let interval = |lower: &ExactSum, upper: &ExactSum| Interval {
            lower: lower.finish().to_python_string(),
            upper: upper.finish().to_python_string(),
        };
        let bounded_missing_contribution = interval(&self.lower, &self.upper);
        let mut withheld = Vec::new();
        if self.missing != self.bounded {
            withheld.push(format!(
                "{} missing terms have no valid contribution bounds",
                self.missing - self.bounded
            ));
        }
        if !pack.header.inventory.consistent {
            withheld.push("inconsistent_term_inventory".into());
        }
        if pack.header.origin == crate::scoring_file::Origin::PgsCatalog
            && pack.header.matches_publication != Some(true)
        {
            withheld.push("scoring_file_not_confirmed_to_match_publication".into());
        }
        let completion = if withheld.is_empty() {
            let mut lower = sum.clone();
            lower.merge(self.lower);
            let mut upper = sum.clone();
            upper.merge(self.upper);
            Some(interval(&lower, &upper))
        } else {
            None
        };
        Ok(Report {
            schema: "pgsum-missingness-v1".into(), missing_terms: self.missing, bounded_terms: self.bounded,
            unbounded_terms: self.missing - self.bounded, by_status: self.by_status, bounded_missing_contribution,
            completion_score_bounds: completion, completion_bounds_withheld_because: withheld,
            ranked_missing_terms: self.ranked.into_iter().map(&describe).collect::<Result<_>>()?,
            unbounded_examples: self.unbounded.into_iter().map(&describe).collect::<Result<_>>()?,
            ranking: "descending maximum absolute possible contribution; source ordinal breaks ties".into(),
            top_terms_limit: TOP_TERMS, source_input_sha256: source_hash.to_owned(),
            assumptions: vec![
                "Bounds concern the published contribution functions for terms in this pack; they are not confidence intervals or risk estimates.".into(),
                "Scored contributions remain fixed, including selected DS/GP expectations. Missing contributions range over coded dosages 0,1,2, or 0,k for explicitly configured XY non-PAR chrX.".into(),
                "Per-term extrema are summed independently. Linkage and shared sites can make the envelope conservative; its endpoints need not be jointly attainable.".into(),
                "Ambiguous models and unsupported ploidy are not assigned bounds. Any such missing term withholds a whole-score envelope.".into(),
                "Frequency filling and reference-panel placement do not change this report; it describes the unfilled partial score.".into(),
            ],
        })
    }
}

fn single_evidence(c: &Candidate, genotypes: &GenotypeTable, options: &Options) -> Result<Evidence> {
    let t = &c.term;
    let key = match crate::score::plan(t, options).0 {
        Plan::Snv { key, .. } => Some(key),
        Plan::Sequence { variant, .. } => genotypes.get_sequence(t.contig, &variant)?.map(|e| e.key),
        Plan::Unscorable(_) if t.contig > 0 && t.pos > 0 => Some(position_key(t.contig, t.pos)),
        _ => None,
    };
    let retained = key.map(|key| genotypes.get(key)).transpose()?.flatten().is_some();
    let lines = key.map(|key| genotypes.records(key)).transpose()?.unwrap_or_default();
    let records = lines
        .iter()
        .take(3)
        .map(|text| {
            let truncated = text.chars().count() > 1024;
            Record {
                sha256: crate::pack::records_sha256(text.as_bytes()),
                text: text.chars().take(1024).collect(),
                truncated,
            }
        })
        .collect();
    let evidence = Evidence {
        availability: if !retained {
            "position_evidence_not_extracted"
        } else if lines.is_empty() {
            "no_overlapping_record"
        } else {
            "retained_records"
        }
        .into(),
        retained_record_count: lines.len(),
        records,
        omitted_record_count: lines.len().saturating_sub(3),
        source_record_ordinals: Vec::new(),
    };
    Ok(evidence)
}

fn describe(c: Candidate, pack: &Pack, evidence: Evidence) -> Result<Term> {
    let t = &c.term;
    let alleles = pack.alleles(c.ordinal - 1);
    let (contribution_bounds, impact, dosages, why) = match c.bounds {
        Ok(b) => (
            Some(Interval {
                lower: b.lower.to_decimal().to_python_string(),
                upper: b.upper.to_decimal().to_python_string(),
            }),
            Some(b.impact.to_python_string()),
            b.dosages,
            None,
        ),
        Err(why) => (None, None, Vec::new(), Some(why.to_owned())),
    };
    Ok(Term {
        ordinal: c.ordinal,
        contig: t.contig.checked_sub(1).map(|i| CONTIGS[i as usize].to_owned()),
        pos1: (t.pos != 0).then_some(t.pos),
        model: t.model.as_str().into(),
        effect_allele: alleles.map(|a| a.0.to_owned()),
        other_allele: alleles.map(|a| a.1.to_owned()),
        weights: t
            .weights
            .iter()
            .map(|w| w.to_decimal().map(|d| d.to_python_string()))
            .collect(),
        status: c.outcome.status.into(),
        call_state: c.outcome.call_state.into(),
        review_reasons: t.reasons.names().into_iter().map(str::to_owned).collect(),
        contribution_bounds,
        max_absolute_contribution: impact,
        allowed_effect_dosages: dosages,
        bounds_unavailable_because: why,
        evidence,
    })
}

/// Detailed cohort diagnostics, streamed one sample per JSON line. `scored` must use the same pack/options.
/// Scoring remains independent: diagnostic blocks are decompressed only when this report is requested.
pub fn write_cohort(
    pack: &Pack,
    cohort: &crate::cohort::CohortTable,
    options: &Options,
    scored: &crate::cohort::CohortScore,
    out: &mut dyn std::io::Write,
) -> Result<()> {
    use crate::cohort::{MISSING, code};
    let n = cohort.header.samples.len();
    if scored.pgs_id != pack.header.pgs_id || scored.sums.len() != n {
        return invalid!("cohort report and scores do not match");
    }
    let mut tallies: Vec<_> = (0..n).map(|_| Tally::default()).collect();
    let mut reader = cohort.reader();
    for (index, term) in pack.terms_from(0).enumerate() {
        let term = term?;
        let key = match crate::score::plan(&term, options).0 {
            Plan::Unscorable(status) => {
                let outcome = Outcome {
                    status,
                    call_state: "",
                    effect_dosage: None,
                    expected_effect_dosage: None,
                    contribution: None,
                    inferred: None,
                    informational_description_accepted: false,
                };
                for tally in &mut tallies {
                    tally.add(index + 1, &term, &outcome, None);
                }
                continue;
            }
            Plan::Snv { key, effect_is_alt } => Some((key, effect_is_alt)),
            Plan::Sequence { variant, effect_is_alt } => cohort
                .sequence_key(term.contig, &variant)
                .map(|key| (key, effect_is_alt)),
        };
        let Some((key, effect_is_alt)) = key else {
            return invalid!("cohort report target is missing");
        };
        let row = reader
            .row(key)?
            .ok_or_else(|| crate::Error::Invalid("cohort report target is missing".into()))?;
        let measurements: std::collections::HashMap<_, _> =
            row.measurements().iter().map(|(sample, m)| (*sample, m)).collect();
        let contributions =
            [0, 1, 2].map(|d| crate::score::term_contribution(&term, if effect_is_alt { d } else { 2 - d }));
        for (sample, tally) in tallies.iter_mut().enumerate() {
            let value = code(&row, sample);
            if value != MISSING && contributions[value as usize].is_some() {
                continue;
            }
            // Valid quantitative measurements with defined contributions are already counted by the
            // numeric scorer. Reporting only needs excluded terms, not a second expectation calculation.
            if contributions.iter().all(Option::is_some)
                && measurements.get(&sample).is_some_and(|m| {
                    matches!(m, crate::dosage::Measurement::GP(_)) || term.model == crate::term::Model::Additive
                })
            {
                continue;
            }
            // The measurement map already indexes this row. Avoid a linear sample search for every
            // quantitative cell, which would make wide-cohort diagnostics quadratic in sample count.
            let state = if let Some(m) = measurements.get(&sample) {
                Some(if m.alt_dosage().coefficient == "0" {
                    crate::genotype::State::ObservedReference
                } else {
                    crate::genotype::State::ObservedVariant
                })
            } else if value == MISSING {
                cohort.missing_state(key, sample)?
            } else {
                cohort.call_state(key, sample)?
            };
            let outcome = if let Some(state) = state {
                let entry = crate::genotypes::Entry {
                    key,
                    state,
                    alt_dosage: (value != MISSING).then_some(value),
                    measurement: measurements.get(&sample).map(|m| Box::new((*m).clone())),
                    refcall_adapted: false,
                    phased: false,
                    first_ref: 0,
                    ref_count: 0,
                };
                crate::score::called(&term, Some(&entry), effect_is_alt, None)?
            } else {
                Outcome {
                    status: "missing_reason_not_retained",
                    call_state: "missing_reason_not_retained",
                    effect_dosage: None,
                    expected_effect_dosage: None,
                    contribution: None,
                    inferred: None,
                    informational_description_accepted: false,
                }
            };
            tally.add(index + 1, &term, &outcome, None);
        }
    }
    let candidate_key = |c: &Candidate| match crate::score::plan(&c.term, options).0 {
        Plan::Snv { key, .. } => Some(key),
        Plan::Sequence { variant, .. } => cohort.sequence_key(c.term.contig, &variant),
        _ => None,
    };
    let keys: std::collections::BTreeSet<_> = tallies
        .iter()
        .flat_map(|t| t.ranked.iter().chain(&t.unbounded))
        .filter_map(candidate_key)
        .collect();
    let unavailable = if cohort.header.schema != crate::cohort::DIAGNOSTIC_SCHEMA {
        "cohort_diagnostics_not_retained"
    } else {
        "position_evidence_not_extracted"
    };
    let mut evidence = std::collections::BTreeMap::new();
    for key in keys {
        let source = cohort.source_records(key)?;
        let availability = match &source {
            Some(v) if v.is_empty() => "no_overlapping_record",
            Some(_) => "source_record_ordinals",
            None => unavailable,
        }
        .to_owned();
        let source = source.unwrap_or_default();
        evidence.insert(
            key,
            Evidence {
                availability,
                retained_record_count: source.len(),
                records: Vec::new(),
                omitted_record_count: source.len().saturating_sub(3),
                source_record_ordinals: source.into_iter().take(3).collect(),
            },
        );
    }
    for (sample, tally) in tallies.into_iter().enumerate() {
        if tally.missing != scored.total_terms - scored.scorable[sample] {
            return invalid!("cohort diagnostics disagree with scored missing-term count");
        }
        let report = tally.finish_with(pack, &cohort.header.gvcf.sha256, &scored.sums[sample], |c| {
            let records = candidate_key(&c)
                .and_then(|k| evidence.get(&k))
                .cloned()
                .unwrap_or_else(|| Evidence {
                    availability: unavailable.into(),
                    retained_record_count: 0,
                    records: Vec::new(),
                    omitted_record_count: 0,
                    source_record_ordinals: Vec::new(),
                });
            describe(c, pack, records)
        })?;
        let value = serde_json::json!({ "schema": "pgsum-cohort-missingness-v1", "pgs_id": pack.header.pgs_id,
            "sample": cohort.header.samples[sample], "sample_index": sample,
            "partial_raw_score": scored.text(sample), "dosage_field": cohort.header.policy.dosage_field,
            "policy": cohort.header.policy, "missingness": report,
            "source_reference": "Input hash plus one-based record ordinal (headers excluded) and zero-based sample_index identify the original VCF/BCF fields; raw cohort record text is not duplicated." });
        serde_json::to_writer(&mut *out, &value).map_err(|e| crate::Error::Invalid(e.to_string()))?;
        out.write_all(b"\n").map_err(crate::Error::io("<cohort missingness>"))?;
    }
    Ok(())
}

/// Extra TSV fields only for unscored terms; calculations use the same explicit assumptions as the report.
pub(crate) fn tsv_fields(term: &TermRecord, outcome: &Outcome, sex: Option<Sex>) -> [String; 3] {
    if outcome.contribution.is_some() {
        return Default::default();
    }
    match bounds(term, outcome, sex) {
        Ok(b) => [
            b.lower.to_decimal().to_python_string(),
            b.upper.to_decimal().to_python_string(),
            String::new(),
        ],
        Err(why) => [String::new(), String::new(), why.to_owned()],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::orient::{Method, Orientation, Status};
    use crate::pack::Weight;
    use crate::term::{AlleleKind, Model, Reasons};
    fn term(model: Model, weights: &[&str]) -> TermRecord {
        TermRecord {
            contig: 1,
            pos: 1,
            model,
            allele_kind: AlleleKind::LiteralSnv,
            palindromic: false,
            orientation: Orientation {
                status: Status::Resolved,
                method: Method::Direct,
                ref_base: b'A',
                alt_base: b'G',
                effect_is_alt: true,
            },
            reasons: Reasons::default(),
            weights: weights
                .iter()
                .map(|w| Weight::from_decimal(Decimal::parse(w).as_ref(), w))
                .collect(),
            inferred: None,
            informational_description: false,
            inferred_sequence: None,
        }
    }
    fn missing() -> Outcome {
        Outcome {
            status: "unknown_no_record",
            call_state: "unknown_no_record",
            effect_dosage: None,
            expected_effect_dosage: None,
            contribution: None,
            inferred: None,
            informational_description_accepted: false,
        }
    }
    #[test]
    fn exact_ordering_and_bounds_preserve_extreme_weights() {
        let values = [
            "-1e100",
            "-1.00000000000000000000000001",
            "-1",
            "-1e-1000",
            "0",
            "1e-1000",
            "1",
            "1.00000000000000000000000001",
            "1e100",
        ];
        for (i, a) in values.iter().enumerate() {
            for (j, b) in values.iter().enumerate() {
                assert_eq!(
                    Decimal::parse(a).unwrap().numeric_cmp(&Decimal::parse(b).unwrap()),
                    i.cmp(&j)
                );
            }
        }
        assert_eq!(
            Decimal::parse("-0.00")
                .unwrap()
                .numeric_cmp(&Decimal::parse("0").unwrap()),
            Ordering::Equal
        );
        assert_eq!(
            Decimal::parse("1.0")
                .unwrap()
                .numeric_cmp(&Decimal::parse("1.000").unwrap()),
            Ordering::Equal
        );
        let b = bounds(
            &term(Model::DosageWeights, &["1.00000000000000000000000001", "-1e-1000", "1"]),
            &missing(),
            None,
        )
        .unwrap();
        assert_eq!(b.lower.to_decimal().to_python_string(), "-1E-1000");
        assert_eq!(b.upper.to_decimal().to_python_string(), "1.00000000000000000000000001");
    }
    #[test]
    fn sex_chromosome_bounds_follow_coded_dosages_and_par_boundaries() {
        let mut t = term(Model::DosageWeights, &["1", "8", "-2"]);
        t.contig = 23;
        t.pos = 3_000_000;
        for (k, lower, upper) in [(1, "1", "8"), (2, "-2", "1")] {
            let b = bounds(&t, &missing(), Some(Sex { xy: true, k })).unwrap();
            assert_eq!(b.dosages, [0, k]);
            assert_eq!(b.lower.to_decimal().to_python_string(), lower);
            assert_eq!(b.upper.to_decimal().to_python_string(), upper);
        }
        t.pos = 10001;
        assert_eq!(
            bounds(&t, &missing(), Some(Sex { xy: true, k: 1 })).unwrap().dosages,
            [0, 1, 2]
        );
        t.pos = 10000;
        assert_eq!(
            bounds(&t, &missing(), Some(Sex { xy: true, k: 1 })).unwrap().dosages,
            [0, 1]
        );
        let mut o = missing();
        o.call_state = "unsupported_ploidy";
        assert!(bounds(&t, &o, None).is_err());
        t.reasons.insert(Reason::SpecialInteraction);
        assert!(bounds(&t, &missing(), None).is_err());
        t.reasons = Reasons(1 << 31);
        assert!(bounds(&t, &missing(), None).is_err());
    }
    #[test]
    fn ranking_is_bounded_and_chunk_merge_is_deterministic() {
        let mut all = Tally::default();
        let mut a = Tally::default();
        let mut b = Tally::default();
        for ordinal in (1..=300).rev() {
            let t = term(Model::Additive, &[if ordinal % 2 == 0 { "-2" } else { "1" }]);
            all.add(ordinal, &t, &missing(), None);
            if ordinal % 3 == 0 {
                a.add(ordinal, &t, &missing(), None);
            } else {
                b.add(ordinal, &t, &missing(), None);
            }
        }
        a.merge(b);
        let ids = |v: &Tally| v.ranked.iter().map(|c| c.ordinal).collect::<Vec<_>>();
        assert_eq!(ids(&all), (1..=20).map(|n| n * 2).collect::<Vec<_>>());
        assert_eq!(ids(&all), ids(&a));
        assert_eq!(all.lower.finish(), a.lower.finish());
        assert_eq!(all.upper.finish(), a.upper.finish());
    }
}
