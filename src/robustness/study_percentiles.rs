//! Seed-level cutoff transitions and comparisons to frozen reference distributions.
use super::*;

const CROSSING_METRICS: [&str; 7] = [
    "up_count",
    "down_count",
    "unchanged_count",
    "changed_count",
    "up_fraction",
    "down_fraction",
    "changed_fraction",
];
const REFERENCE_METRICS: [&str; 4] = [
    "mean_baseline_percentile",
    "mean_scenario_percentile",
    "mean_absolute_shift",
    "max_absolute_shift",
];

// Structured dimensions prevent authored reference IDs from colliding with cutoff/metric names.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Dimension {
    scope: &'static str,
    cutoff: Option<String>,
    reference_set: Option<String>,
    metric: &'static str,
}
impl Dimension {
    fn row(&self, policy: &str, rate: &str) -> Value {
        json!({"policy":policy,"rate":rate.parse::<f64>().unwrap(),"scope":self.scope,
            "cutoff":self.cutoff.as_ref().map(|r|r.parse::<f64>().unwrap()),
            "reference_set":self.reference_set,"metric":self.metric})
    }
}
#[derive(Clone, Copy, Serialize)]
struct Observation {
    value: Option<f64>,
    evaluable: u64,
}
fn count(v: &Value, key: &str) -> Result<u64> {
    v[key]
        .as_u64()
        .ok_or_else(|| Error::Invalid(format!("invalid percentile count {key}")))
}
fn percentile(v: &Value, key: &str, evaluable: u64) -> Result<Option<f64>> {
    let v = v
        .get(key)
        .ok_or_else(|| Error::Invalid(format!("missing reference metric {key}")))?;
    if evaluable == 0 && v.is_null() {
        return Ok(None);
    }
    if evaluable > 0
        && let Some(x) = v.as_f64().filter(|x| x.is_finite() && (0.0..=100.0).contains(x))
    {
        return Ok(Some(x));
    }
    invalid!("reference metric {key} disagrees with evaluable count or percentile bounds")
}

pub(super) fn summarize(
    experiment: &Experiment,
    groups: &BTreeMap<(String, String), BTreeMap<u64, &Arm>>,
    scenarios: &BTreeMap<&str, &Value>,
    reference_policy: &str,
    source_samples: Option<u64>,
) -> Result<(Vec<Value>, Vec<Value>, Vec<Value>)> {
    let expected_cutoffs: BTreeSet<_> = experiment.cutoffs.iter().map(|r| rate_key(*r)).collect();
    let expected_refs: BTreeSet<_> = experiment.reference_sets.iter().map(|r| r.id.as_str()).collect();
    let mut reference_identity = BTreeMap::new();
    let mut observations = BTreeMap::new();
    for (id, scenario) in scenarios {
        let n = count(&scenario["metrics"], "n")?;
        let mut values = BTreeMap::new();
        let mut cutoffs = BTreeSet::new();
        for c in array(&scenario["metrics"], "crossings")? {
            let cutoff = c["cutoff"]
                .as_f64()
                .filter(|x| x.is_finite() && (0.0..=100.0).contains(x))
                .ok_or_else(|| Error::Invalid("invalid crossing cutoff".into()))?;
            let key = rate_key(cutoff);
            if !cutoffs.insert(key.clone()) {
                return invalid!("duplicate cutoff crossing");
            }
            let (up, down, unchanged, evaluable) = (
                count(c, "up")?,
                count(c, "down")?,
                count(c, "unchanged")?,
                count(c, "evaluable")?,
            );
            if evaluable > n || up.checked_add(down).and_then(|v| v.checked_add(unchanged)) != Some(evaluable) {
                return invalid!("crossing counts do not reconcile with evaluable participants");
            }
            for (metric, numerator) in
                CROSSING_METRICS
                    .into_iter()
                    .zip([up, down, unchanged, up + down, up, down, up + down])
            {
                let value = (evaluable > 0).then(|| {
                    if metric.ends_with("_fraction") {
                        numerator as f64 / evaluable as f64
                    } else {
                        numerator as f64
                    }
                });
                values.insert(
                    Dimension {
                        scope: "cutoff",
                        cutoff: Some(key.clone()),
                        reference_set: None,
                        metric,
                    },
                    Observation { value, evaluable },
                );
            }
        }
        if cutoffs != expected_cutoffs {
            return invalid!("cutoff inventory differs from experiment");
        }
        let mut refs = BTreeSet::new();
        for r in array(scenario, "reference_comparisons")? {
            let reference = text(r, "id")?;
            if !refs.insert(reference) {
                return invalid!("duplicate reference comparison");
            }
            let (samples, shared, evaluable) = (
                count(r, "reference_samples")?,
                count(
                    r,
                    if source_samples.is_some() {
                        "source_cohort_shared_sample_ids"
                    } else {
                        "shared_sample_ids"
                    },
                )?,
                count(r, "evaluable")?,
            );
            if samples == 0
                || shared > source_samples.unwrap_or(n).min(samples)
                || evaluable > n
                || (evaluable > 0 && samples < 2)
            {
                return invalid!("invalid reference participant counts");
            }
            let baseline = percentile(r, "mean_baseline_percentile", evaluable)?;
            let mut identity = json!({"id":reference,"reference_samples":samples,
                "mean_baseline_percentile":baseline,"evaluable":evaluable});
            identity[if source_samples.is_some() {
                "source_cohort_shared_sample_ids"
            } else {
                "shared_sample_ids"
            }] = json!(shared);
            // The current engine uses one immutable reference distribution and complete target scores.
            // Its eligibility and baseline percentile must therefore be identical across scenarios.
            if let Some(previous) = reference_identity.insert(reference, identity.clone())
                && previous != identity
            {
                return invalid!("frozen reference identity or baseline changes across scenarios");
            }
            for metric in REFERENCE_METRICS {
                values.insert(
                    Dimension {
                        scope: "reference",
                        cutoff: None,
                        reference_set: Some(reference.to_owned()),
                        metric,
                    },
                    Observation {
                        value: percentile(r, metric, evaluable)?,
                        evaluable,
                    },
                );
            }
            if let (Some(mean), Some(max)) = (r["mean_absolute_shift"].as_f64(), r["max_absolute_shift"].as_f64())
                // Native means sum participant floats; allow accumulation roundoff for identical shifts.
                && mean > max + f64::EPSILON * evaluable as f64 * max.max(1.0)
            {
                return invalid!("mean reference shift exceeds maximum shift");
            }
        }
        if refs != expected_refs {
            return invalid!("reference inventory differs from experiment");
        }
        observations.insert(*id, values);
    }
    let mut summaries = Vec::new();
    let mut comparisons = Vec::new();
    for ((policy, rate), arms) in groups {
        let dimensions = &observations[arms.values().next().unwrap().scenario.as_str()];
        for dimension in dimensions.keys() {
            let seeds: Vec<_> = arms
                .iter()
                .map(|(seed, arm)| {
                    let obs = observations[arm.scenario.as_str()][dimension];
                    json!({"seed":seed,"value":obs.value,"evaluable":obs.evaluable})
                })
                .collect();
            let mut row = dimension.row(policy, rate);
            row["distribution"] = value(&distribution(
                &seeds.iter().map(|s| s["value"].as_f64()).collect::<Vec<_>>(),
            ))?;
            row["seeds"] = json!(seeds);
            summaries.push(row);
            if policy != reference_policy {
                let refs = &groups[&(reference_policy.to_owned(), rate.clone())];
                let pairs: Vec<_> = arms
                    .iter()
                    .map(|(seed, arm)| {
                        let obs = observations[arm.scenario.as_str()][dimension];
                        let reference = observations[refs[seed].scenario.as_str()][dimension];
                        json!({"seed":seed,"policy_value":obs.value,"reference_value":reference.value,
                        "policy_evaluable":obs.evaluable,"reference_evaluable":reference.evaluable,
                        "difference":obs.value.zip(reference.value).map(|(a,b)|a-b)})
                    })
                    .collect();
                let mut row = dimension.row(policy, rate);
                row["reference_policy"] = json!(reference_policy);
                row["distribution"] = value(&distribution(
                    &pairs.iter().map(|p| p["difference"].as_f64()).collect::<Vec<_>>(),
                ))?;
                row["pairs"] = json!(pairs);
                comparisons.push(row);
            }
        }
    }
    Ok((summaries, comparisons, reference_identity.into_values().collect()))
}
