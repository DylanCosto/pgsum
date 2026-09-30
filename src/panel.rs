//! Reference distributions: a sample's score placed among a panel's scores over the same terms.
//!
//! A partial score is only comparable with scores computed over the same terms. For each score, the terms used
//! are those scorable in the sample and called in at least 99% of the panel (remaining missing panel calls
//! filled; see `place`). The panel is a cohort file from `extract --all-samples` extracted with the same packs.
//! The sample's sum is exact; panel sums, percentiles and moments use floating point.
//! Panel samples belong to groups (for 1000 Genomes, superpopulations), and the sample is assigned the
//! nearest group by the likelihood of its genotypes under each group's allele frequencies at a subset of the
//! panel's SNV targets.

use std::collections::BTreeMap;
use std::path::Path;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::cohort::{CohortTable, MISSING, code};
use crate::genotypes::{ANY_ALT, GenotypeTable, is_sequence_key, unpack_key};
use crate::pack::Pack;
use crate::score::{ExactSum, Options, Plan, outcome, plan, term_contribution, term_effect};
use crate::{Error, Result, invalid};

/// Panel samples' groups, read from a TSV with a header naming a `sample` column and the group column.
/// Panel samples the file does not list (for example relatives of other panel samples) are left out of every
/// comparison.
pub struct Groups {
    pub column: String,
    pub names: Vec<String>,
    /// Per panel sample (in the cohort file's order): the index of its group in `names`, if it is listed.
    pub of_sample: Vec<Option<usize>>,
    /// The panel samples that are listed, in the cohort file's order.
    pub members: Vec<usize>,
}

pub fn read_groups(path: &Path, cohort: &CohortTable, column: &str) -> Result<Groups> {
    let text = std::fs::read_to_string(path).map_err(Error::io(path))?;
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap_or_default().split('\t').map(str::trim).collect();
    let (Some(s), Some(g)) = (
        header.iter().position(|h| *h == "sample"),
        header.iter().position(|h| *h == column),
    ) else {
        return invalid!(
            "{}: needs a header with `sample` and `{column}` columns",
            path.display()
        );
    };
    let mut group_of = std::collections::HashMap::new();
    for line in lines.filter(|l| !l.trim().is_empty()) {
        let f: Vec<&str> = line.split('\t').map(str::trim).collect();
        if let (Some(sample), Some(group)) = (f.get(s), f.get(g)) {
            group_of.insert((*sample).to_owned(), (*group).to_owned());
        }
    }
    let mut names: Vec<String> = Vec::new();
    let mut of_sample = Vec::with_capacity(cohort.header.samples.len());
    let mut members = Vec::new();
    for (index, sample) in cohort.header.samples.iter().enumerate() {
        let Some(group) = group_of.get(sample) else {
            of_sample.push(None);
            continue;
        };
        let i = match names.iter().position(|n| n == group) {
            Some(i) => i,
            None => {
                names.push(group.clone());
                names.len() - 1
            }
        };
        of_sample.push(Some(i));
        members.push(index);
    }
    if members.is_empty() {
        return invalid!("{}: lists none of the panel's samples", path.display());
    }
    Ok(Groups {
        column: column.to_owned(),
        names,
        of_sample,
        members,
    })
}

/// The panel group nearest to a sample.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Ancestry {
    pub nearest_group: String,
    /// Log-likelihood of the sample's genotypes under each group's allele frequencies, less the best.
    pub log_likelihood_vs_nearest: BTreeMap<String, f64>,
    pub sites: u64,
    pub note: String,
}

/// Every `SITE_STRIDE`-th panel SNV target is used to assign the nearest group.
const SITE_STRIDE: usize = 20;

/// Assign the sample to the panel group under whose allele frequencies its genotypes are most likely, at
/// every `SITE_STRIDE`-th SNV target where the sample has a passing call.
pub fn nearest_group(table: &GenotypeTable, cohort: &CohortTable, groups: &Groups) -> Result<Ancestry> {
    if cohort.header.policy.dosage_field != crate::dosage::Field::Gt {
        return invalid!(
            "quantitative reference panels require a dosage-aware ancestry model; use a GT reference panel"
        );
    }
    let k = groups.names.len();
    let per_block: Vec<Result<(Vec<f64>, u64)>> = (0..cohort.header.blocks.len())
        .into_par_iter()
        .map(|b| {
            let block = cohort.block_rows(b)?;
            let (keys, rows) = (block.keys(), block.rows());
            let width = rows.len() / keys.len().max(1);
            let mut ll = vec![0f64; k];
            let mut sites = 0u64;
            for (i, &key) in keys.iter().enumerate().step_by(SITE_STRIDE) {
                if is_sequence_key(key) || unpack_key(key).3 == ANY_ALT {
                    continue;
                }
                let Some(entry) = table.get(key)? else { continue };
                let Some(d) = entry.alt_dosage.filter(|_| entry.state.is_passing()) else {
                    continue;
                };
                let row = &rows[i * width..(i + 1) * width];
                let (mut alt, mut n) = (vec![0u32; k], vec![0u32; k]);
                for &s in &groups.members {
                    let g = groups.of_sample[s].expect("members are listed");
                    let c = code(row, s);
                    if c != MISSING {
                        alt[g] += c as u32;
                        n[g] += 2;
                    }
                }
                if n.contains(&0) {
                    continue;
                }
                sites += 1;
                for g in 0..k {
                    let p = (alt[g] as f64 / n[g] as f64).clamp(1e-3, 1.0 - 1e-3);
                    ll[g] += match d {
                        0 => 2.0 * (1.0 - p).ln(),
                        1 => (2.0 * p * (1.0 - p)).ln(),
                        _ => 2.0 * p.ln(),
                    };
                }
            }
            Ok((ll, sites))
        })
        .collect();
    let mut ll = vec![0f64; k];
    let mut sites = 0;
    for r in per_block {
        let (l, n) = r?;
        sites += n;
        for g in 0..k {
            ll[g] += l[g];
        }
    }
    if sites == 0 {
        return invalid!("no panel SNV target has a passing call in the sample");
    }
    let best = (0..k)
        .max_by(|&a, &b| ll[a].total_cmp(&ll[b]))
        .expect("at least one group");
    Ok(Ancestry {
        nearest_group: groups.names[best].clone(),
        log_likelihood_vs_nearest: (0..k).map(|g| (groups.names[g].clone(), ll[g] - ll[best])).collect(),
        sites,
        note: format!(
            "Nearest {} by genotype likelihood at {sites} score sites; sites are not independent (linkage), so \
             differences are larger than a formal test would give. Admixed samples may sit between groups.",
            groups.column
        ),
    })
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct GroupPlacement {
    pub group: String,
    pub samples: u64,
    pub mean: f64,
    pub sd: f64,
    /// Mid-rank percentile of the sample's score among the group's scores (0–100).
    pub percentile: f64,
    /// Undefined (JSON null) when the group has fewer than two samples or zero score variance.
    pub z: Option<f64>,
}

/// A score placed among the panel's scores over the same terms.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Placement {
    pub panel: String,
    pub panel_samples: u64,
    /// Terms scorable in the sample and called in at least 99% of panel samples, and their share of all terms and
    /// weight (panel calls missing at them are filled; see `fills`).
    pub matched_terms: u64,
    /// Scorable sample terms excluded from both sides because fewer than 99% of panel samples have calls.
    pub excluded_low_call_rate_terms: u64,
    pub matched_term_fraction: f64,
    pub matched_weight_fraction: f64,
    /// Whether the matched terms reach `score::COVERAGE_GUIDELINE` of both terms and weight. Below it the
    /// comparison rests on a subset of the score and the percentile should not be read as the score's.
    pub meets_coverage_guideline: bool,
    /// The sample's exact sum over the matched terms.
    pub score: String,
    /// Absent when no term is scorable in the sample and meets the panel call-rate threshold.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub percentile_all: Option<f64>,
    pub groups: Vec<GroupPlacement>,
    /// The group whose percentile to read, when the nearest group was assigned.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub nearest_group: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub nearest_group_percentile: Option<f64>,
    /// Panel calls filled in for the comparison (the sample's own calls never are).
    pub fills: PanelFills,
    pub note: String,
}

/// A matched term needs a passing call in this share of panel samples; the rest are filled (see `place`).
pub const PANEL_CALL_RATE: f64 = 0.99;

/// Panel genotypes filled in for one score.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct PanelFills {
    /// Terms called in at least `PANEL_CALL_RATE` of panel samples where some were missing, and the missing
    /// sample-term pairs, each filled with its group's expected contribution.
    pub terms_filled: u64,
    pub genotypes_filled: u64,
    pub max_per_sample: u32,
    pub samples_filled: u64,
    /// Filled terms per panel sample, in the order of the groups file's listed samples in the panel.
    pub per_sample: Vec<u32>,
}

fn percentile(value: f64, others: &[f64]) -> f64 {
    let below = others.iter().filter(|&&v| v < value).count() as f64;
    let equal = others.iter().filter(|&&v| v == value).count() as f64;
    100.0 * (below + equal / 2.0) / others.len().max(1) as f64
}

fn moments(v: &[f64]) -> (f64, f64) {
    let n = v.len().max(1) as f64;
    let mean = v.iter().sum::<f64>() / n;
    let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0).max(1.0);
    (mean, var.sqrt())
}

/// Place the sample's score for one pack among the panel's scores over the same terms.
pub fn place(
    pack: &Pack,
    table: &GenotypeTable,
    cohort: &CohortTable,
    panel_name: &str,
    groups: &Groups,
    ancestry: Option<&Ancestry>,
    options: &Options,
) -> Result<Placement> {
    if cohort.header.policy.dosage_field != crate::dosage::Field::Gt {
        return invalid!("quantitative reference-panel placement is not supported yet; use a GT reference panel");
    }
    let h = &pack.header;
    if !cohort
        .header
        .packs
        .iter()
        .any(|p| p.pgs_id == h.pgs_id && p.records_sha256 == h.records_sha256)
    {
        return invalid!("the reference panel was not extracted with this {} pack", h.pgs_id);
    }
    let members = &groups.members;
    let n = members.len();
    // The sample's sum is exact; the panel's sums only place it, so they are accumulated in floating point
    // (always in the same order, so reproducibly). `sums[i]` is for panel sample `members[i]`.
    let mut ours = ExactSum::default();
    let mut sums = vec![0f64; n];
    let as_f64 = |c: &crate::score::Contribution| c.to_f64();
    let k = groups.names.len();
    let (mut terms_filled, mut genotypes_filled, mut excluded_low_call_rate_terms) = (0u64, 0u64, 0u64);
    let mut filled = vec![0u32; n];
    let (mut matched, mut total, mut effect_all, mut effect_matched) = (0u64, 0u64, 0f64, 0f64);
    let mut reader = cohort.reader();
    for term in pack.terms_from(0) {
        let term = term?;
        total += 1;
        let effect = term_effect(&term);
        effect_all += effect.unwrap_or(0.0);
        let Some(c) = outcome(&term, table, options)?.contribution else {
            continue;
        };
        let (key, effect_is_alt) = match plan(&term, options).0 {
            Plan::Unscorable(_) => continue,
            Plan::Snv { key, effect_is_alt } => (Some(key), effect_is_alt),
            Plan::Sequence { variant, effect_is_alt } => (cohort.sequence_key(term.contig, &variant), effect_is_alt),
        };
        let Some(row) = key.map(|k| reader.row(k)).transpose()?.flatten() else {
            continue;
        };
        // The contribution at each effect-allele dosage; a term whose model leaves one undefined is not used.
        let Some(by_dosage) = [0u8, 1, 2]
            .map(|d| term_contribution(&term, d).as_ref().map(as_f64))
            .into_iter()
            .collect::<Option<Vec<f64>>>()
        else {
            continue;
        };
        let effect_dosage = |c: u8| if effect_is_alt { c } else { 2 - c };
        // Called panel samples and their effect-allele count, overall and per group.
        let mut called = 0usize;
        let (mut effect_alleles, mut group_called) = (vec![0u64; k], vec![0u64; k]);
        for &s in members {
            let c = code(&row, s);
            if c != MISSING {
                let g = groups.of_sample[s].expect("members are listed");
                called += 1;
                group_called[g] += 1;
                effect_alleles[g] += effect_dosage(c) as u64;
            }
        }
        if (called as f64) < PANEL_CALL_RATE * n as f64 {
            // Keep the same informative terms on both sides. Replacing only the panel's contributions
            // with constants would retain the sample's deviation while erasing the panel's variation.
            excluded_low_call_rate_terms += 1;
            continue;
        }
        let all_frequency = effect_alleles.iter().sum::<u64>() as f64 / (2 * called) as f64;
        // Expected contribution in each group under Hardy-Weinberg at the group's effect-allele frequency among
        // its called samples (the whole panel's, for a group with none called).
        let expected: Vec<f64> = (0..k)
            .map(|g| {
                let f = if group_called[g] > 0 {
                    effect_alleles[g] as f64 / (2 * group_called[g]) as f64
                } else {
                    all_frequency
                };
                (1.0 - f) * (1.0 - f) * by_dosage[0] + 2.0 * f * (1.0 - f) * by_dosage[1] + f * f * by_dosage[2]
            })
            .collect();
        matched += 1;
        effect_matched += effect.unwrap_or(0.0);
        ours.add(&c);
        let mut filled_here = 0u64;
        for (i, (&s, sum)) in members.iter().zip(sums.iter_mut()).enumerate() {
            let c = code(&row, s);
            if c == MISSING {
                *sum += expected[groups.of_sample[s].expect("members are listed")];
                filled[i] += 1;
                filled_here += 1;
            } else {
                *sum += by_dosage[effect_dosage(c) as usize];
            }
        }
        if filled_here > 0 {
            terms_filled += 1;
            genotypes_filled += filled_here;
        }
    }
    let scores = sums;
    let fills = PanelFills {
        terms_filled,
        genotypes_filled,
        max_per_sample: filled.iter().copied().max().unwrap_or(0),
        samples_filled: filled.iter().filter(|&&f| f > 0).count() as u64,
        per_sample: filled,
    };
    let ours_text = ours.finish().to_python_string();
    let v = ours_text.parse::<f64>().unwrap_or(f64::NAN);
    if matched == 0 {
        return Ok(Placement {
            panel: panel_name.to_owned(),
            panel_samples: n as u64,
            matched_terms: 0,
            excluded_low_call_rate_terms,
            matched_term_fraction: 0.0,
            matched_weight_fraction: 0.0,
            meets_coverage_guideline: false,
            score: ours_text,
            percentile_all: None,
            groups: Vec::new(),
            nearest_group: ancestry.map(|a| a.nearest_group.clone()),
            nearest_group_percentile: None,
            fills,
            note: "No scorable sample term reaches the 99% panel call-rate threshold, so the score cannot be placed. Low-call-rate terms are excluded from both sides.".into(),
        });
    }
    let mut placements = Vec::new();
    for (g, name) in groups.names.iter().enumerate() {
        let in_group: Vec<f64> = scores
            .iter()
            .zip(members)
            .filter(|(_, s)| groups.of_sample[**s] == Some(g))
            .map(|(v, _)| *v)
            .collect();
        let (mean, sd) = moments(&in_group);
        placements.push(GroupPlacement {
            group: name.clone(),
            samples: in_group.len() as u64,
            mean,
            sd,
            percentile: percentile(v, &in_group),
            z: (in_group.len() >= 2 && sd > 0.0).then(|| (v - mean) / sd),
        });
    }
    placements.sort_by(|a, b| a.group.cmp(&b.group));
    let nearest = ancestry.map(|a| a.nearest_group.clone());
    let nearest_percentile = nearest
        .as_ref()
        .and_then(|g| placements.iter().find(|p| &p.group == g))
        .map(|p| p.percentile);
    let term_fraction = if total == 0 { 0.0 } else { matched as f64 / total as f64 };
    let weight_fraction = if effect_all == 0.0 {
        0.0
    } else {
        effect_matched / effect_all
    };
    let meets =
        term_fraction >= crate::score::COVERAGE_GUIDELINE && weight_fraction >= crate::score::COVERAGE_GUIDELINE;
    let mut note = String::from(
        "Percentiles compare sums over the same terms: scorable in the sample and called in at least 99% of \
         panel samples. They are uncalibrated: no ancestry adjustment beyond choosing the group, and no absolute \
         risk. Coverage is a completeness check, not a guarantee of percentile stability. Z is null for a group \
         with fewer than two samples or zero score variance.",
    );
    if fills.terms_filled > 0 {
        note.push_str(&format!(
            " Panel calls were filled for the comparison only: {} genotypes at {} terms took their group's expected \
             contribution under Hardy–Weinberg at the effect-allele frequency among called samples; at most {} terms were \
             filled for one panel sample.",
            fills.genotypes_filled, fills.terms_filled, fills.max_per_sample
        ));
    }
    if excluded_low_call_rate_terms > 0 {
        note.push_str(&format!(
            " {excluded_low_call_rate_terms} low-call-rate terms were excluded from both the sample and panel sums and count against matched coverage."
        ));
    }
    if !meets {
        note.push_str(&format!(
            " The matched terms hold only {:.1}% of the terms and {:.1}% of the weight, below the 99% guideline: \
             this places a subset of the score, not the score.",
            100.0 * term_fraction,
            100.0 * weight_fraction
        ));
    }
    Ok(Placement {
        panel: panel_name.to_owned(),
        panel_samples: n as u64,
        matched_terms: matched,
        excluded_low_call_rate_terms,
        matched_term_fraction: term_fraction,
        matched_weight_fraction: weight_fraction,
        meets_coverage_guideline: meets,
        score: ours_text,
        percentile_all: Some(percentile(v, &scores)),
        groups: placements,
        nearest_group: nearest,
        nearest_group_percentile: nearest_percentile,
        fills,
        note,
    })
}
