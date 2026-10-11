//! Descriptive sensitivity statistics, not estimates of predictive accuracy.
//!
//! Percentiles always compare against the fixed baseline distribution. Spearman correlation instead
//! ranks each cohort independently, with average ranks for ties. Decimal ordering and subtraction are
//! exact; aggregate magnitudes and correlations are floating-point approximations.

use std::cmp::Ordering;
use std::collections::BTreeMap;

use num_bigint::BigInt;
use serde::Serialize;

use crate::decimal::Decimal;
use crate::score::{Coefficient, Contribution, ExactSum};
use crate::{Result, invalid};

/// Descriptive percentile-cutoff transitions. A percentile equal to the cutoff is on its upper side.
#[derive(Clone, Debug, Serialize)]
pub struct Crossing {
    pub cutoff: f64,
    pub up: usize,
    pub down: usize,
    pub unchanged: usize,
    pub evaluable: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct Summary {
    pub n: usize,
    pub mae: Option<f64>,
    pub rmse: Option<f64>,
    pub max_abs_change: Option<f64>,
    pub pearson: Option<f64>,
    pub spearman: Option<f64>,
    pub mean_abs_percentile_shift: Option<f64>,
    pub max_abs_percentile_shift: Option<f64>,
    pub crossings: Vec<Crossing>,
    /// Reasons are keyed by the corresponding unavailable numeric field, or `crossings`.
    pub unavailable_reasons: BTreeMap<String, String>,
}

/// Empirical mid-CDF percentiles: `100 * (less + 0.5 * equal) / reference.len()`.
///
/// Comparisons are exact, including ties with different decimal precision and signed zeros. Empty,
/// single-observation, constant, or invalid reference distributions yield `None` for every value:
/// these references do not support informative percentile movement. An invalid individual value
/// yields `None` for that value. Values below/above the reference range have percentiles 0/100.
pub fn percentiles(values: &[Decimal], reference: &[Decimal]) -> Vec<Option<f64>> {
    if reference_issue(reference).is_some() {
        return vec![None; values.len()];
    }
    let mut sorted: Vec<&Decimal> = reference.iter().collect();
    sorted.sort_by(|a, b| a.numeric_cmp(b));
    values
        .iter()
        .map(|value| {
            if !valid_decimal(value) {
                return None;
            }
            let less = sorted.partition_point(|candidate| candidate.numeric_cmp(value) == Ordering::Less);
            let through_equal = sorted.partition_point(|candidate| candidate.numeric_cmp(value) != Ordering::Greater);
            Some(100.0 * (less as f64 + 0.5 * (through_equal - less) as f64) / sorted.len() as f64)
        })
        .collect()
}

/// Summarize aligned baseline/scenario pairs. Length mismatch, malformed decimal representations, or
/// non-finite/out-of-range percentiles are input errors. Valid inputs with statistically undefined or
/// unrepresentable metrics return `None` and a reason, never NaN/infinity or a misleading zero.
pub fn summarize(baseline: &[Decimal], scenario: &[Decimal], cutoffs: &[f64]) -> Result<Summary> {
    if baseline.len() != scenario.len() {
        return invalid!(
            "baseline and scenario lengths differ: {} versus {}",
            baseline.len(),
            scenario.len()
        );
    }
    if baseline.iter().chain(scenario).any(|value| !valid_decimal(value)) {
        return invalid!("sensitivity statistics require valid finite decimal representations");
    }
    if cutoffs
        .iter()
        .any(|cutoff| !cutoff.is_finite() || !(0.0..=100.0).contains(cutoff))
    {
        return invalid!("percentile cutoffs must be finite numbers between 0 and 100");
    }

    let mut summary = Summary {
        n: baseline.len(),
        mae: None,
        rmse: None,
        max_abs_change: None,
        pearson: None,
        spearman: None,
        mean_abs_percentile_shift: None,
        max_abs_percentile_shift: None,
        crossings: cutoffs
            .iter()
            .map(|&cutoff| Crossing {
                cutoff,
                up: 0,
                down: 0,
                unchanged: 0,
                evaluable: 0,
            })
            .collect(),
        unavailable_reasons: BTreeMap::new(),
    };

    if baseline.is_empty() {
        mark_unavailable(
            &mut summary,
            &["mae", "rmse", "max_abs_change", "pearson", "spearman"],
            "no paired scores",
        );
    } else {
        let changes = baseline
            .iter()
            .zip(scenario)
            .map(|(base, changed)| exact_difference(changed, base).and_then(|d| decimal_float(&d)))
            .collect::<Option<Vec<_>>>();
        if let Some(changes) = changes {
            let max = changes.iter().map(|value| value.abs()).fold(0.0, f64::max);
            summary.max_abs_change = Some(max);
            if max == 0.0 {
                summary.mae = Some(0.0);
                summary.rmse = Some(0.0);
            } else {
                // Normalize before squaring/summing so representable large changes do not overflow.
                let mean = compensated_sum(changes.iter().map(|value| value.abs() / max)) / changes.len() as f64;
                let mean_square =
                    compensated_sum(changes.iter().map(|value| (value / max).powi(2))) / changes.len() as f64;
                summary.mae = positive_finite(max * mean);
                summary.rmse = positive_finite(max * mean_square.sqrt());
                for metric in ["mae", "rmse"] {
                    let missing = if metric == "mae" {
                        summary.mae.is_none()
                    } else {
                        summary.rmse.is_none()
                    };
                    if missing {
                        mark_unavailable(&mut summary, &[metric], "metric exceeds floating-point range");
                    }
                }
            }
        } else {
            mark_unavailable(
                &mut summary,
                &["mae", "rmse", "max_abs_change"],
                "an exact score difference exceeds supported arithmetic or floating-point range",
            );
        }

        if baseline.len() < 2 || is_constant(baseline) || is_constant(scenario) {
            mark_unavailable(
                &mut summary,
                &["pearson", "spearman"],
                "correlation requires at least two pairs and variation in both score vectors",
            );
        } else {
            // Subtract a shared origin exactly BEFORE conversion. This preserves differences such as
            // 9007199254740993 - 9007199254740992 that ordinary f64 subtraction would lose.
            summary.pearson = scaled_offsets(baseline)
                .zip(scaled_offsets(scenario))
                .and_then(|(x, y)| correlation(&x, &y));
            summary.spearman = correlation(&ranks(baseline), &ranks(scenario));
            if summary.pearson.is_none() {
                mark_unavailable(
                    &mut summary,
                    &["pearson"],
                    "exact score offsets cannot be represented with sufficient floating-point range",
                );
            }
            if summary.spearman.is_none() {
                mark_unavailable(&mut summary, &["spearman"], "rank correlation is numerically undefined");
            }
        }
    }

    if let Some(reason) = reference_issue(baseline) {
        mark_unavailable(
            &mut summary,
            &["mean_abs_percentile_shift", "max_abs_percentile_shift", "crossings"],
            reason,
        );
    } else {
        let before = percentiles(baseline, baseline);
        let after = percentiles(scenario, baseline);
        let shifts: Vec<f64> = before
            .iter()
            .zip(&after)
            .filter_map(|(&base, &changed)| base.zip(changed).map(|(b, s)| (s - b).abs()))
            .collect();
        summary.mean_abs_percentile_shift = Some(compensated_sum(shifts.iter().copied()) / shifts.len() as f64);
        summary.max_abs_percentile_shift = Some(shifts.iter().copied().fold(0.0, f64::max));
        for crossing in &mut summary.crossings {
            for (&base, &changed) in before.iter().zip(&after) {
                if let Some((base, changed)) = base.zip(changed) {
                    crossing.evaluable += 1;
                    match (base >= crossing.cutoff, changed >= crossing.cutoff) {
                        (false, true) => crossing.up += 1,
                        (true, false) => crossing.down += 1,
                        _ => crossing.unchanged += 1,
                    }
                }
            }
        }
    }
    Ok(summary)
}

fn mark_unavailable(summary: &mut Summary, metrics: &[&str], reason: &str) {
    for metric in metrics {
        summary.unavailable_reasons.insert((*metric).into(), reason.into());
    }
}

fn valid_decimal(value: &Decimal) -> bool {
    !value.coefficient.is_empty()
        && value.coefficient.bytes().all(|digit| digit.is_ascii_digit())
        && (value.coefficient == "0" || !value.coefficient.starts_with('0'))
        && i64::try_from(value.coefficient.len() - 1)
            .ok()
            .and_then(|digits| value.exponent.checked_add(digits))
            .is_some()
}

fn is_constant(values: &[Decimal]) -> bool {
    values
        .first()
        .is_none_or(|first| values.iter().all(|value| value.numeric_cmp(first) == Ordering::Equal))
}

fn reference_issue(reference: &[Decimal]) -> Option<&'static str> {
    if reference.iter().any(|value| !valid_decimal(value)) {
        Some("reference contains an invalid decimal")
    } else if reference.len() < 2 {
        Some("percentile reference requires at least two scores")
    } else if is_constant(reference) {
        Some("percentile reference is constant")
    } else {
        None
    }
}

fn exact_difference(value: &Decimal, origin: &Decimal) -> Option<Decimal> {
    // Score decimals ordinarily span only a few thousand digits. Bound arbitrary public Decimal
    // inputs before ExactSum allocates powers of ten; never silently wrap its u32 exponent conversion.
    const MAX_EXACT_DIGITS: i64 = 100_000;
    let floor = value.exponent.min(origin.exponent).min(0);
    for decimal in [value, origin] {
        let span = decimal.exponent.checked_sub(floor)?;
        let width = span.checked_add(i64::try_from(decimal.coefficient.len()).ok()?)?;
        if width > MAX_EXACT_DIGITS {
            return None;
        }
    }
    let mut sum = ExactSum::default();
    for (decimal, negate) in [(value, false), (origin, true)] {
        sum.add(&Contribution {
            negative: decimal.negative ^ negate,
            coefficient: Coefficient::Big(decimal.coefficient.parse::<BigInt>().ok()?),
            exponent: decimal.exponent,
        });
    }
    Some(sum.finish())
}

fn decimal_float(value: &Decimal) -> Option<f64> {
    // Parse the entire scientific decimal, avoiding `huge coefficient * underflowed power of ten`.
    let converted: f64 = format!(
        "{}{}e{}",
        if value.negative { "-" } else { "" },
        value.coefficient,
        value.exponent
    )
    .parse()
    .ok()?;
    (converted.is_finite() && (converted != 0.0 || value.coefficient == "0")).then_some(converted)
}

fn positive_finite(value: f64) -> Option<f64> {
    (value.is_finite() && value > 0.0).then_some(value)
}

fn scaled_offsets(values: &[Decimal]) -> Option<Vec<f64>> {
    let first = values.first()?;
    let offsets = values
        .iter()
        .map(|value| exact_difference(value, first))
        .collect::<Option<Vec<_>>>()?;
    let scale = offsets
        .iter()
        .filter(|value| value.coefficient != "0")
        .map(Decimal::adjusted)
        .max()?;
    offsets
        .into_iter()
        .map(|mut value| {
            value.exponent = value.exponent.checked_sub(scale)?;
            decimal_float(&value)
        })
        .collect()
}

/// Average ordinal ranks, computed with exact decimal ordering rather than converted score values.
fn ranks(values: &[Decimal]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|&a, &b| values[a].numeric_cmp(&values[b]));
    let mut ranks = vec![0.0; values.len()];
    let mut start = 0;
    while start < order.len() {
        let mut end = start + 1;
        while end < order.len() && values[order[start]].numeric_cmp(&values[order[end]]) == Ordering::Equal {
            end += 1;
        }
        let rank = 0.5 * ((start + 1) as f64 + end as f64);
        for &index in &order[start..end] {
            ranks[index] = rank;
        }
        start = end;
    }
    ranks
}

fn correlation(x: &[f64], y: &[f64]) -> Option<f64> {
    if x.len() < 2 || x.len() != y.len() {
        return None;
    }
    let mean_x = compensated_sum(x.iter().copied()) / x.len() as f64;
    let mean_y = compensated_sum(y.iter().copied()) / y.len() as f64;
    let xx = compensated_sum(x.iter().map(|value| (value - mean_x).powi(2)));
    let yy = compensated_sum(y.iter().map(|value| (value - mean_y).powi(2)));
    let xy = compensated_sum(x.iter().zip(y).map(|(a, b)| (a - mean_x) * (b - mean_y)));
    let denominator = xx.sqrt() * yy.sqrt();
    let r = xy / denominator;
    (denominator > 0.0 && r.is_finite()).then(|| r.clamp(-1.0, 1.0))
}

fn compensated_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut sum, mut correction) = (0.0, 0.0);
    for value in values {
        let corrected = value - correction;
        let next = sum + corrected;
        correction = (next - sum) - corrected;
        sum = next;
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;

    fn decimals(values: &[&str]) -> Vec<Decimal> {
        values.iter().map(|value| Decimal::parse(value).unwrap()).collect()
    }

    fn close(actual: Option<f64>, expected: f64) {
        assert!((actual.unwrap() - expected).abs() < 1e-12, "{actual:?} != {expected}");
    }

    #[test]
    fn hand_calculated_summary_uses_fixed_baseline_percentiles() {
        let base = decimals(&["1", "2", "3"]);
        let changed = decimals(&["2", "2", "1"]);
        let summary = summarize(&base, &changed, &[50.0]).unwrap();
        assert_eq!(summary.n, 3);
        close(summary.mae, 1.0);
        close(summary.rmse, (5.0_f64 / 3.0).sqrt());
        close(summary.max_abs_change, 2.0);
        close(summary.pearson, -3.0_f64.sqrt() / 2.0);
        close(summary.spearman, -3.0_f64.sqrt() / 2.0);
        close(summary.mean_abs_percentile_shift, 100.0 / 3.0);
        close(summary.max_abs_percentile_shift, 200.0 / 3.0);
        let crossing = &summary.crossings[0];
        assert_eq!(
            (crossing.up, crossing.down, crossing.unchanged, crossing.evaluable),
            (1, 1, 1, 3)
        );
        assert!(summary.unavailable_reasons.is_empty());
    }

    #[test]
    fn ties_negative_values_and_signed_zero_are_exact() {
        let reference = decimals(&["-2", "0", "-0.00", "2.0", "2.00"]);
        let values = decimals(&["-3", "-2", "-1", "-0", "1", "2", "3"]);
        assert_eq!(
            percentiles(&values, &reference),
            [0.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0].map(Some)
        );
        assert_eq!(ranks(&reference), vec![1.0, 2.5, 2.5, 4.5, 4.5]);
    }

    #[test]
    fn preserves_differences_above_float_integer_precision() {
        let base = decimals(&["9007199254740992", "9007199254740993", "9007199254740994"]);
        let changed = decimals(&["9007199254740993", "9007199254740994", "9007199254740995"]);
        let summary = summarize(&base, &changed, &[]).unwrap();
        close(summary.mae, 1.0);
        close(summary.rmse, 1.0);
        close(summary.max_abs_change, 1.0);
        close(summary.pearson, 1.0);
        close(summary.spearman, 1.0);
        close(summary.mean_abs_percentile_shift, 250.0 / 9.0);
        assert_eq!(
            percentiles(&base, &base),
            vec![Some(100.0 / 6.0), Some(50.0), Some(500.0 / 6.0)]
        );
    }

    #[test]
    fn big_coefficients_do_not_pass_through_small_integer_accumulators() {
        let base = decimals(&[
            "10000000000000000000000000000000000000000",
            "10000000000000000000000000000000000000001",
        ]);
        let changed = decimals(&[
            "10000000000000000000000000000000000000001",
            "10000000000000000000000000000000000000002",
        ]);
        let summary = summarize(&base, &changed, &[]).unwrap();
        close(summary.mae, 1.0);
        close(summary.pearson, 1.0);
        close(summary.spearman, 1.0);
    }

    #[test]
    fn tiny_values_keep_ranks_and_scaled_correlation_but_not_false_zero_errors() {
        let base = decimals(&["1e-1000", "2e-1000", "3e-1000"]);
        let changed = decimals(&["2e-1000", "4e-1000", "6e-1000"]);
        let summary = summarize(&base, &changed, &[]).unwrap();
        assert!(summary.mae.is_none());
        assert!(summary.rmse.is_none());
        assert!(summary.max_abs_change.is_none());
        assert!(summary.unavailable_reasons.contains_key("mae"));
        close(summary.pearson, 1.0);
        close(summary.spearman, 1.0);
        assert!(summary.mean_abs_percentile_shift.is_some());
    }

    #[test]
    fn empty_singleton_and_constant_references_are_unavailable() {
        for reference in [vec![], decimals(&["2"]), decimals(&["2", "2.00"])] {
            assert_eq!(percentiles(&decimals(&["1", "2", "3"]), &reference), vec![None; 3]);
            let summary = summarize(&reference, &reference, &[50.0]).unwrap();
            assert!(summary.pearson.is_none());
            assert!(summary.spearman.is_none());
            assert!(summary.mean_abs_percentile_shift.is_none());
            assert!(summary.max_abs_percentile_shift.is_none());
            assert_eq!(summary.crossings[0].evaluable, 0);
            assert!(summary.unavailable_reasons.contains_key("crossings"));
            if !reference.is_empty() {
                close(summary.mae, 0.0);
                close(summary.rmse, 0.0);
            }
        }
        let summary = summarize(&decimals(&["1"]), &decimals(&["3"]), &[]).unwrap();
        close(summary.mae, 2.0);
        close(summary.rmse, 2.0);
    }

    #[test]
    fn constant_scenario_still_has_descriptive_percentile_movement() {
        let summary = summarize(&decimals(&["-1", "0", "1"]), &decimals(&["0", "0", "0"]), &[50.0]).unwrap();
        assert!(summary.pearson.is_none());
        assert!(summary.spearman.is_none());
        close(summary.mean_abs_percentile_shift, 200.0 / 9.0);
        assert_eq!(summary.crossings[0].evaluable, 3);
    }

    #[test]
    fn ranks_are_recomputed_for_spearman() {
        let summary = summarize(&decimals(&["1", "2", "3"]), &decimals(&["101", "102", "103"]), &[90.0]).unwrap();
        close(summary.spearman, 1.0);
        close(summary.mean_abs_percentile_shift, 50.0);
        assert_eq!(summary.crossings[0].up, 3);
    }

    #[test]
    fn large_finite_errors_do_not_overflow_the_squared_error_sum() {
        let base = vec![
            Decimal {
                negative: false,
                coefficient: "0".into(),
                exponent: 0
            };
            2
        ];
        let changed = vec![
            Decimal {
                negative: false,
                coefficient: "1".into(),
                exponent: 200
            };
            2
        ];
        let summary = summarize(&base, &changed, &[]).unwrap();
        assert_eq!(summary.mae, Some(1e200));
        assert_eq!(summary.rmse, Some(1e200));
        assert_eq!(summary.max_abs_change, Some(1e200));
        assert!(serde_json::to_string(&summary).is_ok());
    }

    #[test]
    fn invalid_inputs_fail_without_nan_json() {
        assert!(summarize(&decimals(&["1"]), &[], &[]).is_err());
        for cutoff in [f64::NAN, f64::INFINITY, -1.0, 101.0] {
            assert!(summarize(&[], &[], &[cutoff]).is_err());
        }
        let malformed = Decimal {
            negative: false,
            coefficient: "NaN".into(),
            exponent: 0,
        };
        assert!(summarize(std::slice::from_ref(&malformed), std::slice::from_ref(&malformed), &[]).is_err());
        assert_eq!(percentiles(&[malformed], &decimals(&["1", "2"])), vec![None]);
        let huge = Decimal {
            negative: false,
            coefficient: "1".into(),
            exponent: i64::MAX,
        };
        let summary = summarize(&decimals(&["0"]), &[huge], &[]).unwrap();
        assert!(summary.mae.is_none());
        assert!(serde_json::to_string(&summary).is_ok());
    }
}
