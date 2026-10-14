//! Descriptive variation across masking seeds, within one verified experiment.
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use clap::Args;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use super::config::{self, Experiment, Mask};
use crate::{Error, Result, invalid};

#[path = "study_percentiles.rs"]
mod percentiles;

const METRICS: [&str; 7] = [
    "mae",
    "rmse",
    "max_abs_change",
    "pearson",
    "spearman",
    "mean_abs_percentile_shift",
    "max_abs_percentile_shift",
];

#[derive(Args, Debug)]
pub struct StudyArgs {
    /// Planner output directory containing plan.json and experiment.json.
    #[arg(long)]
    pub plan: PathBuf,
    /// Completed robustness output directory.
    #[arg(long)]
    pub results: PathBuf,
    /// Authored policy ID against which other policies are paired by seed.
    #[arg(long)]
    pub reference_policy: String,
    /// New directory for JSON, TSV and HTML summaries.
    #[arg(long)]
    pub out: PathBuf,
}

#[derive(Deserialize)]
struct Plan {
    schema: String,
    experiment_sha256: String,
    mask_algorithm: String,
    scenarios: Vec<Arm>,
}
#[derive(Clone, Deserialize)]
#[serde(deny_unknown_fields)]
struct Arm {
    scenario: String,
    policy: String,
    rate: Option<f64>,
    seed: Option<u64>,
}
#[derive(Debug, Serialize)]
struct Distribution {
    requested: usize,
    available: usize,
    unavailable: usize,
    mean: Option<f64>,
    median: Option<f64>,
    min: Option<f64>,
    max: Option<f64>,
    sample_sd: Option<f64>,
}
fn distribution(values: &[Option<f64>]) -> Distribution {
    let mut xs: Vec<_> = values.iter().flatten().copied().filter(|v| v.is_finite()).collect();
    xs.sort_by(f64::total_cmp);
    let n = xs.len();
    let finite = |x: f64| x.is_finite().then_some(x);
    let mean = (n > 0)
        .then(|| xs.iter().map(|x| x / n as f64).sum::<f64>())
        .and_then(finite);
    let median = (n > 0).then(|| {
        if n % 2 == 1 {
            xs[n / 2]
        } else {
            xs[n / 2 - 1] / 2.0 + xs[n / 2] / 2.0
        }
    });
    // Scale before squaring to avoid overflow for ordinary large-valued scores.
    let sample_sd = mean.and_then(|m| {
        if n < 2 {
            return None;
        }
        let scale = xs.iter().map(|x| x.abs()).fold(m.abs(), f64::max);
        if scale == 0.0 {
            return Some(0.0);
        }
        finite(scale * (xs.iter().map(|x| (x / scale - m / scale).powi(2)).sum::<f64>() / (n - 1) as f64).sqrt())
    });
    Distribution {
        requested: values.len(),
        available: n,
        unavailable: values.len() - n,
        mean,
        median,
        min: xs.first().copied(),
        max: xs.last().copied(),
        sample_sd,
    }
}
fn value<T: Serialize>(v: &T) -> Result<Value> {
    serde_json::to_value(v).map_err(|e| Error::Invalid(e.to_string()))
}
fn read<T: for<'a> Deserialize<'a>>(path: &Path) -> Result<T> {
    serde_json::from_slice(&std::fs::read(path).map_err(Error::io(path))?)
        .map_err(|e| Error::Invalid(format!("{}: {e}", path.display())))
}
fn array<'a>(v: &'a Value, key: &str) -> Result<&'a Vec<Value>> {
    v[key]
        .as_array()
        .ok_or_else(|| Error::Invalid(format!("missing array {key}")))
}
fn text<'a>(v: &'a Value, key: &str) -> Result<&'a str> {
    v[key]
        .as_str()
        .ok_or_else(|| Error::Invalid(format!("missing string {key}")))
}
fn rate_key(rate: f64) -> String {
    if rate == 0.0 { "0".into() } else { rate.to_string() }
}

fn validate_plan(plan: &Plan, experiment: &Experiment) -> Result<()> {
    if plan.schema != "pgsum-robustness-plan-v1" || plan.mask_algorithm != config::MASK_ALGORITHM {
        return invalid!("unsupported study plan or mask algorithm");
    }
    let scenarios: BTreeMap<_, _> = experiment.scenarios.iter().map(|s| (s.id.as_str(), s)).collect();
    let controls: BTreeMap<_, _> = plan
        .scenarios
        .iter()
        .filter(|a| a.rate.is_none() && a.seed.is_none())
        .map(|a| (a.policy.as_str(), a))
        .collect();
    let mut seen = BTreeSet::new();
    let mut cells = BTreeSet::new();
    let mut grids: BTreeMap<&str, BTreeSet<(String, u64)>> = BTreeMap::new();
    for arm in &plan.scenarios {
        if !seen.insert(arm.scenario.as_str()) {
            return invalid!("duplicate scenario in plan");
        }
        let scenario = scenarios
            .get(arm.scenario.as_str())
            .ok_or_else(|| Error::Invalid("plan scenario absent from experiment".into()))?;
        match (arm.rate, arm.seed) {
            (None, None) => {
                if arm.scenario != arm.policy || scenario.masks.iter().any(|m| matches!(m, Mask::Random { .. })) {
                    return invalid!("plan control differs from its policy");
                }
            }
            (Some(rate), Some(seed)) if rate.is_finite() && (0.0..=1.0).contains(&rate) => {
                let control = controls
                    .get(arm.policy.as_str())
                    .and_then(|a| scenarios.get(a.scenario.as_str()))
                    .ok_or_else(|| Error::Invalid("missing policy control".into()))?;
                let mut expected = (*scenario).clone();
                match expected.masks.pop() {
                    Some(Mask::Random { rate: r, seed: s, .. }) if r == rate && s == seed => (),
                    _ => return invalid!("plan rate/seed does not match the actual scenario mask"),
                }
                expected.id = control.id.clone();
                if value(&expected)? != value(control)? {
                    return invalid!("generated scenario changes more than its planned random mask");
                }
                let key = rate_key(rate);
                if !cells.insert((arm.policy.clone(), key.clone(), seed)) {
                    return invalid!("duplicate policy/rate/seed cell");
                }
                grids.entry(&arm.policy).or_default().insert((key, seed));
            }
            _ => return invalid!("plan arms need both rate and seed or neither"),
        }
    }
    if seen.len() != scenarios.len() || controls.is_empty() || grids.len() != controls.len() {
        return invalid!("incomplete study plan");
    }
    let grid = grids.values().next().unwrap();
    let rates: BTreeSet<_> = grid.iter().map(|(r, _)| r).collect();
    let seeds: BTreeSet<_> = grid.iter().map(|(_, s)| s).collect();
    if grid.len() != rates.len() * seeds.len() || grids.values().any(|g| g != grid) {
        return invalid!("policies must share a complete rate/seed grid");
    }
    // Pairing requires the same random-mask unit for every generated arm.
    let mut unit = None;
    for a in plan.scenarios.iter().filter(|a| a.seed.is_some()) {
        if let Some(Mask::Random { unit: u, .. }) = scenarios[a.scenario.as_str()].masks.last() {
            if unit.is_some_and(|v| v != *u) {
                return invalid!("paired policies use different random-mask units");
            }
            unit = Some(*u);
        }
    }
    Ok(())
}

fn summarize(plan: &Plan, experiment: &Experiment, report: &Value, reference: &str) -> Result<Value> {
    validate_plan(plan, experiment)?;
    if report["schema"] != "pgsum-robustness-report-v1"
        || report["mask_algorithm"] != config::MASK_ALGORITHM
        || report["experiment"] != value(experiment)?
    {
        return invalid!("results do not belong to the planned experiment");
    }
    if !plan.scenarios.iter().any(|a| a.seed.is_none() && a.policy == reference) {
        return invalid!("reference policy is absent from plan");
    }
    let mut groups: BTreeMap<(String, String), BTreeMap<u64, &Arm>> = BTreeMap::new();
    for arm in &plan.scenarios {
        if let (Some(rate), Some(seed)) = (arm.rate, arm.seed) {
            groups
                .entry((arm.policy.clone(), rate_key(rate)))
                .or_default()
                .insert(seed, arm);
        }
    }
    let mut scores = Vec::new();
    let mut ids = BTreeSet::new();
    for score in array(report, "scores")? {
        let id = text(score, "pgs_id")?;
        if !ids.insert(id) {
            return invalid!("duplicate score ID");
        }
        if score["baseline_terms"].as_u64().is_none_or(|n| n == 0) || text(score, "baseline_term_sha256")?.is_empty() {
            return invalid!("missing common-baseline identity");
        }
        let mut scenarios = BTreeMap::new();
        let mut n = None;
        for s in array(score, "scenarios")? {
            if scenarios.insert(text(s, "id")?, s).is_some() {
                return invalid!("duplicate result scenario");
            }
            let count = s["metrics"]["n"]
                .as_u64()
                .filter(|n| *n > 0)
                .ok_or_else(|| Error::Invalid("invalid sample count".into()))?;
            if n.is_some_and(|v| v != count) {
                return invalid!("scenario sample counts differ");
            }
            n = Some(count);
            for metric in METRICS {
                let v = s["metrics"]
                    .get(metric)
                    .ok_or_else(|| Error::Invalid(format!("missing metric {metric}")))?;
                if !v.is_null() && v.as_f64().is_none_or(|v| !v.is_finite()) {
                    return invalid!("invalid metric {metric}");
                }
            }
        }
        if scenarios.len() != plan.scenarios.len()
            || plan
                .scenarios
                .iter()
                .any(|a| !scenarios.contains_key(a.scenario.as_str()))
        {
            return invalid!("missing or unexpected scenario results");
        }
        let (percentile_summaries, percentile_comparisons, reference_sets) =
            percentiles::summarize(experiment, &groups, &scenarios, reference)?;
        let mut summaries = Vec::new();
        let mut comparisons = Vec::new();
        for ((policy, rate), arms) in &groups {
            for metric in METRICS {
                let values: Vec<_> = arms
                    .values()
                    .map(|a| scenarios[a.scenario.as_str()]["metrics"][metric].as_f64())
                    .collect();
                summaries.push(json!({"policy":policy,"rate":rate.parse::<f64>().unwrap(),"metric":metric,"distribution":distribution(&values)}));
                if policy != reference {
                    let refs = &groups[&(reference.to_owned(), rate.clone())];
                    let pairs: Vec<_> = arms
                        .iter()
                        .map(|(seed, a)| {
                            let v = scenarios[a.scenario.as_str()]["metrics"][metric].as_f64();
                            let r = scenarios[refs[seed].scenario.as_str()]["metrics"][metric].as_f64();
                            let delta = v.zip(r).and_then(|(v, r)| (v - r).is_finite().then_some(v - r));
                            json!({"seed":seed,"policy_value":v,"reference_value":r,"difference":delta})
                        })
                        .collect();
                    let deltas: Vec<_> = pairs.iter().map(|p| p["difference"].as_f64()).collect();
                    comparisons.push(json!({"policy":policy,"reference_policy":reference,"rate":rate.parse::<f64>().unwrap(),"metric":metric,"distribution":distribution(&deltas),"pairs":pairs}));
                }
            }
        }
        let controls: Vec<_> = plan
            .scenarios
            .iter()
            .filter(|a| a.seed.is_none())
            .map(|a| scenarios[a.scenario.as_str()].clone())
            .collect();
        scores.push(json!({"pgs_id":id,"samples":n,"baseline_terms":score["baseline_terms"],"total_model_terms":score["total_model_terms"],"baseline_term_sha256":score["baseline_term_sha256"],"baseline_term_fraction":score["baseline_term_fraction"],"baseline_weight_fraction":score["baseline_weight_fraction"],"excluded_terms":score["excluded_terms"],"controls":controls,"summaries":summaries,"comparisons":comparisons,"percentile_summaries":percentile_summaries,"percentile_comparisons":percentile_comparisons,"reference_sets":reference_sets}));
    }
    if scores.len() != experiment.packs.len() {
        return invalid!("result score inventory differs from experiment packs");
    }
    Ok(
        json!({"schema":"pgsum-robustness-study-v1","reference_policy":reference,"scores":scores,"notes":[
            "Descriptive variation across masking seeds on one frozen cohort and common-callable baseline; not confidence intervals, phenotype validation or population uncertainty.",
            "Paired differences are policy minus reference policy within the same score, rate and seed. Error metrics favor smaller values; correlation metrics favor larger values. No automatic best-policy selection is performed.",
            "Undefined metrics are excluded with requested/available/unavailable counts. Paired differences require both metrics to be available and their difference representable.",
            "Mean, median, range and sample standard deviation describe the available seed-level metrics, not pooled participant observations. Sample SD is unavailable with fewer than two values; unrepresentable aggregates are null.",
            "Cutoff counts and fractions describe membership transitions at percentile >= cutoff in the frozen target baseline distribution. Fractions range from 0 to 1; zero evaluable participants yields unavailable values, including counts. Per-seed denominators remain explicit.",
            "Reference metrics use the authored frozen reference score distribution. Percentiles and absolute shifts are in percentile points; reference identities, participant counts and target overlap are retained. Neither membership nor percentile shifts imply clinical risk.",
            "Controls are retained separately. Scores with different scales, baselines, cutoffs or reference sets are never pooled."
        ]}),
    )
}

pub fn run(args: &StudyArgs) -> Result<()> {
    let plan_path = args.plan.join("plan.json");
    let experiment_path = args.plan.join("experiment.json");
    let report_path = args.results.join("report.json");
    let mut provenance = BTreeMap::new();
    for p in [&plan_path, &experiment_path, &report_path] {
        provenance.insert(p.clone(), crate::digest::file_sha256(p)?);
    }
    let plan: Plan = read(&plan_path)?;
    if plan.experiment_sha256 != provenance[&experiment_path] {
        return invalid!("planned experiment hash differs");
    }
    let experiment = config::load(&experiment_path)?;
    let report: Value = read(&report_path)?;
    let inputs = report["input_sha256"]
        .as_object()
        .ok_or_else(|| Error::Invalid("missing report input hashes".into()))?;
    if !inputs
        .values()
        .any(|v| v.as_str() == Some(plan.experiment_sha256.as_str()))
    {
        return invalid!("report does not record the planned manifest hash");
    }
    let outputs = report["outputs"]
        .as_object()
        .ok_or_else(|| Error::Invalid("missing report output hashes".into()))?;
    for required in [
        "samples.tsv",
        "top-changes.jsonl",
        "regions.tsv",
        "references.tsv",
        "report.html",
    ] {
        if !outputs.contains_key(required) {
            return invalid!("incomplete report output inventory");
        }
    }
    for (name, hash) in outputs {
        if ![
            "samples.tsv",
            "top-changes.jsonl",
            "regions.tsv",
            "references.tsv",
            "report.html",
            "terms.tsv.zst",
        ]
        .contains(&name.as_str())
        {
            return invalid!("unsupported report output path {name:?}");
        }
        let p = args.results.join(name);
        let actual = crate::digest::file_sha256(&p)?;
        if hash.as_str() != Some(actual.as_str()) {
            return invalid!("report output hash differs: {name}");
        }
        provenance.insert(p, actual);
    }
    let mut summary = summarize(&plan, &experiment, &report, &args.reference_policy)?;
    summary["source_sha256"] = value(&provenance)?;
    summary["pgsum_version"] = json!(env!("CARGO_PKG_VERSION"));
    summary["executable_sha256"] = json!(crate::digest::file_sha256(
        &std::env::current_exe().map_err(Error::io("<executable>"))?
    )?);
    for (p, hash) in &provenance {
        if crate::digest::file_sha256(p)? != *hash {
            return invalid!("input changed during study summary");
        }
    }
    std::fs::create_dir(&args.out).map_err(Error::io(&args.out))?;
    let mut html = String::from(
        "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>PGS robustness study</title><style>body{font:16px system-ui;max-width:1200px;margin:2rem auto;padding:1rem}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:.4rem;border-bottom:1px solid #ddd;text-align:left}th{background:#eee}</style><h1>PGS robustness study</h1><p>Descriptive masking sensitivity, not predictive validation or confidence intervals. Paired differences are policy minus reference.</p>",
    );
    for note in array(&summary, "notes")? {
        html.push_str(&format!("<p>{}</p>", note.as_str().unwrap()));
    }
    for kind in [
        "summaries",
        "comparisons",
        "percentile_summaries",
        "percentile_comparisons",
    ] {
        let extended = kind.starts_with("percentile_");
        let title = match kind {
            "summaries" => "Score sensitivity across seeds",
            "comparisons" => "Paired policy differences",
            "percentile_summaries" => "Cutoff and reference sensitivity across seeds",
            _ => "Paired cutoff and reference differences",
        };
        let mut tsv = String::from(
            "pgs_id\tpolicy\treference_policy\trate\tmetric\trequested\tavailable\tunavailable\tmean\tmedian\tmin\tmax\tsample_sd\n",
        );
        if extended {
            tsv = tsv.replacen("metric\t", "scope\tcutoff\treference_set\tmetric\t", 1);
        }
        let dimensions = if extended {
            "<th>Scope</th><th>Cutoff</th><th>Reference set</th>"
        } else {
            ""
        };
        html.push_str(&format!("<h2>{title}</h2><div style=\"overflow-x:auto\"><table><tr><th>Score</th><th>Policy</th><th>Reference</th><th>Rate</th>{dimensions}<th>Metric</th><th>Seeds</th><th>Available</th><th>Unavailable</th><th>Mean</th><th>Median</th><th>Min</th><th>Max</th><th>SD</th></tr>"));
        for score in array(&summary, "scores")? {
            for row in array(score, kind)? {
                let mut cells = vec![
                    score["pgs_id"].clone(),
                    row["policy"].clone(),
                    row["reference_policy"].clone(),
                    row["rate"].clone(),
                ];
                if extended {
                    cells.extend([
                        row["scope"].clone(),
                        row["cutoff"].clone(),
                        row["reference_set"].clone(),
                    ]);
                }
                cells.push(row["metric"].clone());
                for field in [
                    "requested",
                    "available",
                    "unavailable",
                    "mean",
                    "median",
                    "min",
                    "max",
                    "sample_sd",
                ] {
                    cells.push(row["distribution"][field].clone());
                }
                let cells: Vec<String> = cells
                    .iter()
                    .map(|v| {
                        if v.is_null() {
                            String::new()
                        } else if let Some(s) = v.as_str() {
                            s.to_owned()
                        } else {
                            v.to_string()
                        }
                    })
                    .collect();
                if cells.iter().any(|s| s.contains(['\t', '\n', '\r'])) {
                    return invalid!("invalid tabular identity");
                }
                tsv.push_str(&cells.join("\t"));
                tsv.push('\n');
                html.push_str("<tr>");
                for (index, mut c) in cells.into_iter().enumerate() {
                    if index >= if extended { 11 } else { 8 }
                        && let Ok(v) = c.parse::<f64>()
                    {
                        c = if v != 0.0 && (v.abs() < 1e-4 || v.abs() >= 1e6) {
                            format!("{v:.5e}")
                        } else {
                            format!("{v:.5}")
                        };
                    }
                    html.push_str(&format!(
                        "<td>{}</td>",
                        c.replace('&', "&amp;")
                            .replace('<', "&lt;")
                            .replace('>', "&gt;")
                            .replace('"', "&quot;")
                            .replace('\'', "&#39;")
                    ));
                }
                html.push_str("</tr>");
            }
        }
        html.push_str("</table></div>");
        let p = args.out.join(format!("{kind}.tsv"));
        std::fs::write(&p, tsv).map_err(Error::io(&p))?;
    }
    html.push_str("<p>Blank values are unavailable. HTML numbers are rounded for display; JSON and TSV retain full precision. See summary.json for provenance, per-seed paired values, baseline coverage and controls.</p></html>");
    let p = args.out.join("summary.html");
    std::fs::write(&p, html).map_err(Error::io(&p))?;
    let mut hashes = BTreeMap::new();
    for name in [
        "summaries.tsv",
        "comparisons.tsv",
        "percentile_summaries.tsv",
        "percentile_comparisons.tsv",
        "summary.html",
    ] {
        hashes.insert(name, crate::digest::file_sha256(&args.out.join(name))?);
    }
    summary["outputs"] = value(&hashes)?;
    let tmp = args.out.join(".summary.json.pending");
    super::write_json(&tmp, &summary)?;
    std::fs::rename(&tmp, args.out.join("summary.json")).map_err(Error::io(&tmp))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn descriptive_statistics_handle_missing_singleton_and_large_values() {
        let d = distribution(&[Some(2.0), None, Some(4.0), Some(6.0)]);
        assert_eq!((d.requested, d.available, d.unavailable), (4, 3, 1));
        assert_eq!(
            (d.mean, d.median, d.min, d.max, d.sample_sd),
            (Some(4.0), Some(4.0), Some(2.0), Some(6.0), Some(2.0))
        );
        assert!(distribution(&[Some(2.0)]).sample_sd.is_none());
        let d = distribution(&[None, None]);
        assert!(d.mean.is_none() && d.median.is_none());
        let d = distribution(&[Some(1e300), Some(1e300)]);
        assert_eq!(d.mean, Some(1e300));
        assert_eq!(d.sample_sd, Some(0.0));
        assert!(distribution(&[Some(-f64::MAX), Some(f64::MAX)]).sample_sd.is_none());
    }
    fn fixture() -> (Plan, Experiment, Value) {
        let mut arms = Vec::new();
        let mut scenarios = Vec::new();
        let mut results = Vec::new();
        for policy in ["a", "b"] {
            for seed in [None, Some(1), Some(2)] {
                let id = format!("{policy}-{seed:?}");
                // Control IDs are required to equal their policy names.
                let id = if seed.is_none() { policy.to_owned() } else { id };
                arms.push(json!({"scenario":id,"policy":policy,"rate":seed.map(|_|0.5),"seed":seed}));
                scenarios.push(json!({"id":id,"masks":seed.map(|seed| vec![json!({"kind":"random","rate":0.5,"seed":seed,"unit":"sample_site"})]).unwrap_or_default()}));
                let mae = match (policy, seed) {
                    ("a", Some(1)) => Some(1.0),
                    ("a", Some(2)) => Some(100.0),
                    ("b", Some(1)) => Some(3.0),
                    ("b", Some(2)) => None,
                    _ => Some(0.0),
                };
                let mut metrics = json!({"n":3});
                for m in METRICS {
                    metrics[m] = if m == "mae" { json!(mae) } else { json!(0.5) };
                }
                metrics["crossings"] = json!([
                    {"cutoff":90,"up":0,"down":1,"unchanged":2,"evaluable":3},
                    {"cutoff":95,"up":0,"down":0,"unchanged":3,"evaluable":3}
                ]);
                results.push(json!({"id":id,"metrics":metrics,"reference_comparisons":[]}));
            }
        }
        let experiment: Experiment = serde_json::from_value(
            json!({"schema":config::EXPERIMENT_SCHEMA,"genotypes":"x.pgsc","packs":["x.pgsp"],"scenarios":scenarios}),
        )
        .unwrap();
        let plan:Plan=serde_json::from_value(json!({"schema":"pgsum-robustness-plan-v1","experiment_sha256":"x","mask_algorithm":config::MASK_ALGORITHM,"scenarios":arms})).unwrap();
        let report = json!({"schema":"pgsum-robustness-report-v1","mask_algorithm":config::MASK_ALGORITHM,"experiment":experiment,"scores":[{"pgs_id":"S","baseline_terms":1,"baseline_term_sha256":"x","scenarios":results}]});
        (plan, experiment, report)
    }
    #[test]
    fn differences_are_paired_before_missing_metrics_are_removed() {
        let (p, e, r) = fixture();
        let s = summarize(&p, &e, &r, "a").unwrap();
        let row = s["scores"][0]["comparisons"]
            .as_array()
            .unwrap()
            .iter()
            .find(|r| r["metric"] == "mae")
            .unwrap();
        assert_eq!(row["distribution"]["mean"], 2.0);
        assert_eq!(row["distribution"]["available"], 1);
        assert_eq!(row["distribution"]["unavailable"], 1);
        assert_eq!(row["pairs"][0]["seed"], 1);
        assert!(row["pairs"][1]["difference"].is_null());
    }
    #[test]
    fn percentile_summaries_preserve_dimensions_denominators_and_pairing() {
        let (p, mut e, mut r) = fixture();
        e.reference_sets = serde_json::from_value(json!([
            {"id":"cutoff/90/changed_fraction", "samples":["A","B"]},
            {"id":"singleton", "samples":["A"]}
        ]))
        .unwrap();
        r["experiment"] = json!(e);
        for (i, s) in r["scores"][0]["scenarios"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .enumerate()
        {
            let shift = [0.0, 10.0, 30.0, 0.0, 8.0, 20.0][i];
            s["reference_comparisons"] = json!([
                {"id":"cutoff/90/changed_fraction","reference_samples":2,"shared_sample_ids":2,"evaluable":3,
                 "mean_baseline_percentile":50,"mean_scenario_percentile":50,"mean_absolute_shift":shift,"max_absolute_shift":shift},
                {"id":"singleton","reference_samples":1,"shared_sample_ids":1,"evaluable":0,
                 "mean_baseline_percentile":null,"mean_scenario_percentile":null,"mean_absolute_shift":null,"max_absolute_shift":null}
            ]);
        }
        r["scores"][0]["scenarios"][1]["metrics"]["crossings"][0] =
            json!({"cutoff":90,"up":1,"down":1,"unchanged":1,"evaluable":3});
        let s = summarize(&p, &e, &r, "a").unwrap();
        let rows = s["scores"][0]["percentile_comparisons"].as_array().unwrap();
        let row = rows
            .iter()
            .find(|x| x["metric"] == "changed_fraction" && x["cutoff"] == 90.0)
            .unwrap();
        assert!((row["distribution"]["mean"].as_f64().unwrap() + 1.0 / 6.0).abs() < 1e-12);
        assert_eq!(row["pairs"][0]["policy_evaluable"], 3);
        let row = rows
            .iter()
            .find(|x| x["metric"] == "mean_absolute_shift" && x["reference_set"] == "cutoff/90/changed_fraction")
            .unwrap();
        assert_eq!(row["distribution"]["mean"], -6.0);
        assert_eq!(row["pairs"][0]["difference"], -2.0);
        let row = rows
            .iter()
            .find(|x| x["metric"] == "mean_absolute_shift" && x["reference_set"] == "singleton")
            .unwrap();
        assert_eq!(row["distribution"]["unavailable"], 2);
        assert!(row["distribution"]["mean"].is_null());
        assert_eq!(row["pairs"][0]["reference_evaluable"], 0);
        assert_eq!(s["scores"][0]["reference_sets"].as_array().unwrap().len(), 2);
        let saved = r.clone();
        // Three identical 0.1-point shifts have a floating mean slightly above their maximum.
        let mut rounded_experiment = e.clone();
        rounded_experiment.reference_sets[0].samples = None;
        rounded_experiment.reference_sets[0].genotypes = Some("reference500.pgsc".into());
        r["experiment"] = json!(rounded_experiment);
        for scenario in r["scores"][0]["scenarios"].as_array_mut().unwrap() {
            let reference = &mut scenario["reference_comparisons"][0];
            reference["reference_samples"] = json!(500);
            reference["mean_scenario_percentile"] = json!(50.1);
            reference["mean_absolute_shift"] = json!((0.1 + 0.1 + 0.1) / 3.0);
            reference["max_absolute_shift"] = json!(0.1);
        }
        assert!(summarize(&p, &rounded_experiment, &r, "a").is_ok());
        r["scores"][0]["scenarios"][1]["reference_comparisons"][0]["mean_absolute_shift"] = json!(0.2);
        assert!(summarize(&p, &rounded_experiment, &r, "a").is_err());
        // Reference metadata must agree across arms, and values must honor availability and bounds.
        for (key, bad) in [
            ("reference_samples", json!(3)),
            ("shared_sample_ids", json!(4)),
            ("mean_baseline_percentile", json!(49)),
            ("evaluable", json!(4)),
            ("mean_absolute_shift", json!(101)),
            ("max_absolute_shift", json!(null)),
        ] {
            let mut r = saved.clone();
            r["scores"][0]["scenarios"][1]["reference_comparisons"][0][key] = bad;
            assert!(summarize(&p, &e, &r, "a").is_err(), "{key}");
        }
        let mut r = saved.clone();
        r["scores"][0]["scenarios"][1]["reference_comparisons"]
            .as_array_mut()
            .unwrap()
            .pop();
        assert!(summarize(&p, &e, &r, "a").is_err());
        let mut r = saved;
        r["scores"][0]["scenarios"][1]["reference_comparisons"][1]["id"] = json!("cutoff/90/changed_fraction");
        assert!(summarize(&p, &e, &r, "a").is_err());
    }

    #[test]
    fn undefined_cutoffs_are_unavailable_and_malformed_crossings_fail() {
        let (p, e, mut r) = fixture();
        for scenario in r["scores"][0]["scenarios"].as_array_mut().unwrap() {
            for c in scenario["metrics"]["crossings"].as_array_mut().unwrap() {
                for k in ["up", "down", "unchanged", "evaluable"] {
                    c[k] = json!(0);
                }
            }
        }
        let s = summarize(&p, &e, &r, "a").unwrap();
        for row in s["scores"][0]["percentile_summaries"].as_array().unwrap() {
            assert_eq!(row["distribution"]["unavailable"], 2);
            assert!(row["distribution"]["mean"].is_null());
            assert_eq!(row["seeds"][0]["evaluable"], 0);
        }
        for (key, bad) in [
            ("up", json!(u64::MAX)),
            ("down", json!(-1)),
            ("evaluable", json!(4)),
            ("cutoff", json!(89)),
        ] {
            let (p, e, mut r) = fixture();
            r["scores"][0]["scenarios"][1]["metrics"]["crossings"][0][key] = bad;
            assert!(summarize(&p, &e, &r, "a").is_err(), "{key}");
        }
        let (p, e, mut r) = fixture();
        r["scores"][0]["scenarios"][1]["metrics"]["crossings"][1]["cutoff"] = json!(90);
        assert!(summarize(&p, &e, &r, "a").is_err());
        let (p, e, mut r) = fixture();
        r["scores"][0]["scenarios"][1]["metrics"]["crossings"]
            .as_array_mut()
            .unwrap()
            .pop();
        assert!(summarize(&p, &e, &r, "a").is_err());
    }
    #[test]
    fn rejects_incomplete_mislabeled_or_unpaired_evidence() {
        let (mut p, e, r) = fixture();
        p.scenarios.pop();
        assert!(summarize(&p, &e, &r, "a").is_err());
        let (mut p, e, r) = fixture();
        p.scenarios[1].seed = Some(999);
        assert!(summarize(&p, &e, &r, "a").is_err());
        let (p, mut e, mut r) = fixture();
        if let Mask::Random { unit, .. } = &mut e.scenarios[4].masks[0] {
            *unit = config::RandomUnit::Site;
        }
        r["experiment"] = json!(e);
        assert!(summarize(&p, &e, &r, "a").is_err());
        let (p, e, mut r) = fixture();
        r["scores"][0]["scenarios"].as_array_mut().unwrap().pop();
        assert!(summarize(&p, &e, &r, "a").is_err());
        let (p, e, mut r) = fixture();
        r["scores"][0]["scenarios"][1]["metrics"]["n"] = 4.into();
        assert!(summarize(&p, &e, &r, "a").is_err());
        let (p, e, mut r) = fixture();
        r["experiment"]["genotypes"] = "unrelated".into();
        assert!(summarize(&p, &e, &r, "a").is_err());
        let (p, e, r) = fixture();
        assert!(summarize(&p, &e, &r, "unknown").is_err());
    }
}
