//! Descriptive sensitivity within an explicit, exhaustive partition of the original cohort.
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use clap::Args;
use serde_json::{Value, json};

use super::stats;
use crate::decimal::Decimal;
use crate::{Error, Result, invalid};

#[derive(Args, Debug)]
pub struct GroupsArgs {
    /// Completed native robustness directory containing report.json and its output files.
    #[arg(long)]
    pub results: PathBuf,
    /// Tab-separated sample/group mapping: exactly one row for every original participant.
    #[arg(long)]
    pub groups: PathBuf,
    /// New output directory for group metrics, cutoff/reference tables and an HTML report.
    #[arg(long)]
    pub out: PathBuf,
}

#[derive(Clone)]
struct Pair {
    baseline: Decimal,
    scenario: Decimal,
    terms: usize,
}
type Rows = BTreeMap<(String, String), BTreeMap<String, Pair>>;
type ReferenceKey = (String, String, String, String);
#[derive(Clone, PartialEq)]
struct ReferencePair {
    before: Option<f64>,
    after: Option<f64>,
    samples: usize,
    shared: usize,
}

fn identity(s: &str) -> Result<String> {
    if s.trim().is_empty() || s.chars().any(char::is_control) {
        return invalid!("empty or invalid group/sample/score identity");
    }
    Ok(s.to_owned())
}
fn string<'a>(value: &'a Value, key: &str) -> Result<&'a str> {
    value[key]
        .as_str()
        .ok_or_else(|| Error::Invalid(format!("missing string {key}")))
}
fn array<'a>(value: &'a Value, key: &str) -> Result<&'a Vec<Value>> {
    value[key]
        .as_array()
        .ok_or_else(|| Error::Invalid(format!("missing array {key}")))
}
fn count(value: &str) -> Result<usize> {
    value
        .parse()
        .map_err(|_| Error::Invalid("invalid nonnegative count".into()))
}
fn score(value: &str) -> Result<Decimal> {
    Decimal::parse_score(value).ok_or_else(|| Error::Invalid("invalid exact score decimal".into()))
}
fn percentile(value: &str) -> Result<Option<f64>> {
    if value.is_empty() {
        return Ok(None);
    }
    let x: f64 = value
        .parse()
        .map_err(|_| Error::Invalid("invalid reference percentile".into()))?;
    if !x.is_finite() || !(0.0..=100.0).contains(&x) {
        return invalid!("reference percentile outside 0..100");
    }
    Ok(Some(x))
}

/// Native TSVs are unquoted. Preserve identities literally; never trim or infer membership.
fn table(
    path: &Path,
    required: &[&str],
    mut visit: impl FnMut(&[&str], &BTreeMap<String, usize>) -> Result<()>,
) -> Result<()> {
    let mut lines = BufReader::new(File::open(path).map_err(Error::io(path))?).lines();
    let header = lines
        .next()
        .ok_or_else(|| Error::Invalid("empty TSV".into()))?
        .map_err(Error::io(path))?;
    let fields: Vec<_> = header.split('\t').collect();
    let columns: BTreeMap<_, _> = fields.iter().enumerate().map(|(i, s)| (s.to_string(), i)).collect();
    if columns.len() != fields.len() || fields.contains(&"") || required.iter().any(|key| !columns.contains_key(*key)) {
        return invalid!("{}: missing or duplicate TSV columns", path.display());
    }
    for (i, line) in lines.enumerate() {
        let line = line.map_err(Error::io(path))?;
        let row: Vec<_> = line.split('\t').collect();
        if row.len() != fields.len() {
            return invalid!("{}: row {} has wrong column count", path.display(), i + 2);
        }
        visit(&row, &columns)?;
    }
    Ok(())
}

fn load_groups(path: &Path) -> Result<BTreeMap<String, String>> {
    let mut assignments = BTreeMap::new();
    table(path, &["sample", "group"], |row, columns| {
        if columns.len() != 2 {
            return invalid!("group mapping must have exactly sample and group columns");
        }
        let sample = identity(row[columns["sample"]])?;
        let group = identity(row[columns["group"]])?;
        if assignments.insert(sample, group).is_some() {
            return invalid!("duplicate participant in group mapping");
        }
        Ok(())
    })?;
    if assignments.is_empty() {
        return invalid!("empty group mapping");
    }
    Ok(assignments)
}

fn load_scores(path: &Path) -> Result<Rows> {
    let mut rows: Rows = BTreeMap::new();
    table(
        path,
        &[
            "pgs_id",
            "scenario",
            "sample",
            "baseline_terms",
            "baseline_score",
            "scenario_score",
        ],
        |row, c| {
            let key = (identity(row[c["pgs_id"]])?, identity(row[c["scenario"]])?);
            let sample = identity(row[c["sample"]])?;
            let pair = Pair {
                baseline: score(row[c["baseline_score"]])?,
                scenario: score(row[c["scenario_score"]])?,
                terms: count(row[c["baseline_terms"]])?,
            };
            if rows.entry(key).or_default().insert(sample, pair).is_some() {
                return invalid!("duplicate sample/scenario score row");
            }
            Ok(())
        },
    )?;
    Ok(rows)
}

fn load_references(path: &Path) -> Result<BTreeMap<ReferenceKey, ReferencePair>> {
    let mut rows = BTreeMap::new();
    table(
        path,
        &[
            "pgs_id",
            "scenario",
            "sample",
            "reference",
            "reference_samples",
            "shared_sample_ids",
            "baseline_percentile",
            "scenario_percentile",
            "percentile_delta",
        ],
        |row, c| {
            let key = (
                identity(row[c["pgs_id"]])?,
                identity(row[c["scenario"]])?,
                identity(row[c["reference"]])?,
                identity(row[c["sample"]])?,
            );
            let before = percentile(row[c["baseline_percentile"]])?;
            let after = percentile(row[c["scenario_percentile"]])?;
            let delta = row[c["percentile_delta"]];
            match (before, after) {
                (None, None) if delta.is_empty() => (),
                (Some(a), Some(b)) => {
                    let d: f64 = delta
                        .parse()
                        .map_err(|_| Error::Invalid("invalid reference delta".into()))?;
                    if !d.is_finite() || (b - a - d).abs() > 1e-10 {
                        return invalid!("reference delta does not reconcile");
                    }
                }
                _ => return invalid!("partial reference percentile pair"),
            }
            let value = ReferencePair {
                before,
                after,
                samples: count(row[c["reference_samples"]])?,
                shared: count(row[c["shared_sample_ids"]])?,
            };
            if rows.insert(key, value).is_some() {
                return invalid!("duplicate participant/reference row");
            }
            Ok(())
        },
    )?;
    Ok(rows)
}

fn mean(values: &[f64]) -> Option<f64> {
    (!values.is_empty()).then(|| values.iter().map(|x| x / values.len() as f64).sum())
}

fn summarize(
    report: &Value,
    rows: &Rows,
    references: &BTreeMap<ReferenceKey, ReferencePair>,
    assignments: &BTreeMap<String, String>,
) -> Result<Value> {
    if report["schema"] != "pgsum-robustness-report-v1" || report["mask_algorithm"] != super::config::MASK_ALGORITHM {
        return invalid!("unsupported native robustness report");
    }
    let experiment = &report["experiment"];
    let cutoffs: Vec<f64> =
        serde_json::from_value(experiment["cutoffs"].clone()).map_err(|e| Error::Invalid(e.to_string()))?;
    let ids = |items: &Vec<Value>| -> Result<BTreeSet<String>> {
        let values: Vec<_> = items
            .iter()
            .map(|v| identity(string(v, "id")?))
            .collect::<Result<_>>()?;
        let set: BTreeSet<_> = values.iter().cloned().collect();
        if set.len() != values.len() {
            return invalid!("duplicate declared identity");
        }
        Ok(set)
    };
    let scenarios = ids(array(experiment, "scenarios")?)?;
    let reference_ids = ids(array(experiment, "reference_sets")?)?;
    if scenarios.is_empty() {
        return invalid!("no scenarios");
    }
    let participants: BTreeSet<_> = assignments.keys().cloned().collect();
    let mut groups: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for (sample, group) in assignments {
        groups.entry(group.clone()).or_default().push(sample.clone());
    }
    let mut seen = BTreeSet::new();
    let mut used_rows = BTreeSet::new();
    let mut used_references = BTreeSet::new();
    let mut scores = Vec::new();
    for native in array(report, "scores")? {
        let id = identity(string(native, "pgs_id")?)?;
        if !seen.insert(id.clone()) {
            return invalid!("duplicate result score");
        }
        let terms = native["baseline_terms"]
            .as_u64()
            .filter(|x| *x > 0)
            .ok_or_else(|| Error::Invalid("missing baseline term count".into()))?;
        if string(native, "baseline_term_sha256")?.is_empty() {
            return invalid!("missing baseline digest");
        }
        if ids(array(native, "scenarios")?)? != scenarios {
            return invalid!("scenario inventory differs from experiment");
        }
        let mut baseline: Option<BTreeMap<String, Decimal>> = None;
        let mut fixed_reference: BTreeMap<(String, String), (Option<f64>, usize, usize)> = BTreeMap::new();
        let mut results = Vec::new();
        for scenario in array(native, "scenarios")? {
            let sid = string(scenario, "id")?;
            let key = (id.clone(), sid.to_owned());
            let pairs = rows
                .get(&key)
                .ok_or_else(|| Error::Invalid("missing score/scenario rows".into()))?;
            used_rows.insert(key);
            if pairs.keys().cloned().collect::<BTreeSet<_>>() != participants
                || scenario["metrics"]["n"].as_u64() != Some(participants.len() as u64)
            {
                return invalid!(
                    "group mapping must account for every original participant exactly once; unknown or missing identities"
                );
            }
            let current: BTreeMap<_, _> = pairs
                .iter()
                .map(|(sample, pair)| (sample.clone(), pair.baseline.clone()))
                .collect();
            if let Some(frozen) = &baseline {
                if current
                    .iter()
                    .any(|(sample, value)| value.numeric_cmp(&frozen[sample]) != std::cmp::Ordering::Equal)
                {
                    return invalid!("baseline scores change across scenarios");
                }
            } else {
                baseline = Some(current);
            }
            if pairs.values().any(|p| p.terms as u64 != terms) {
                return invalid!("baseline term count differs in sample rows");
            }
            let full_baseline: Vec<_> = baseline.as_ref().unwrap().values().cloned().collect();
            let native_references = array(scenario, "reference_comparisons")?;
            if ids(native_references)? != reference_ids {
                return invalid!("reference inventory differs from experiment");
            }
            for reference in native_references {
                let rid = string(reference, "id")?;
                for sample in &participants {
                    let key = (id.clone(), sid.to_owned(), rid.to_owned(), sample.clone());
                    let row = references
                        .get(&key)
                        .ok_or_else(|| Error::Invalid("missing participant/reference row".into()))?;
                    if row.samples == 0
                        || row.shared > row.samples.min(participants.len())
                        || reference["reference_samples"].as_u64() != Some(row.samples as u64)
                        || reference["shared_sample_ids"].as_u64() != Some(row.shared as u64)
                    {
                        return invalid!("reference size or source-cohort overlap differs");
                    }
                    let baseline_key = (rid.to_owned(), sample.clone());
                    let value = (row.before, row.samples, row.shared);
                    if fixed_reference
                        .insert(baseline_key, value)
                        .is_some_and(|old| old != value)
                    {
                        return invalid!("frozen reference baseline changes across scenarios");
                    }
                    used_references.insert(key);
                }
            }
            for (group, members) in &groups {
                let a: Vec<_> = members.iter().map(|s| pairs[s].baseline.clone()).collect();
                let b: Vec<_> = members.iter().map(|s| pairs[s].scenario.clone()).collect();
                let metrics = stats::summarize_against(&a, &b, &cutoffs, &full_baseline)?;
                let mut ref_metrics = Vec::new();
                for reference in native_references {
                    let rid = string(reference, "id")?;
                    let selected: Vec<_> = members
                        .iter()
                        .map(|sample| &references[&(id.clone(), sid.to_owned(), rid.to_owned(), sample.clone())])
                        .collect();
                    let before: Vec<_> = selected
                        .iter()
                        .filter_map(|r| r.before.zip(r.after).map(|(a, _)| a))
                        .collect();
                    let after: Vec<_> = selected
                        .iter()
                        .filter_map(|r| r.before.zip(r.after).map(|(_, b)| b))
                        .collect();
                    let shifts: Vec<_> = before.iter().zip(&after).map(|(a, b)| (b - a).abs()).collect();
                    ref_metrics.push(json!({"id":rid,"reference_samples":selected[0].samples,"source_cohort_shared_sample_ids":selected[0].shared,
                        "evaluable":shifts.len(),"unavailable":members.len()-shifts.len(),"mean_baseline_percentile":mean(&before),"mean_scenario_percentile":mean(&after),
                        "mean_absolute_shift":mean(&shifts),"max_absolute_shift":shifts.iter().copied().reduce(f64::max)}));
                }
                results
                    .push(json!({"scenario":sid,"group":group,"metrics":metrics,"reference_comparisons":ref_metrics}));
            }
        }
        scores.push(json!({"pgs_id":id,"source_samples":participants.len(),"baseline_terms":terms,"total_model_terms":native["total_model_terms"],
            "baseline_term_sha256":native["baseline_term_sha256"],"baseline_term_fraction":native["baseline_term_fraction"],"baseline_weight_fraction":native["baseline_weight_fraction"],
            "excluded_terms":native["excluded_terms"],"scenarios":results}));
    }
    if scores.is_empty()
        || scores.len() != array(experiment, "packs")?.len()
        || used_rows.len() != rows.len()
        || used_references.len() != references.len()
    {
        return invalid!("unexpected or missing score/reference evidence");
    }
    Ok(
        json!({"schema":"pgsum-robustness-groups-v1","percentile_basis":"full_target_baseline","correlation_basis":"within_group",
        "groups":groups.iter().map(|(id,members)|json!({"id":id,"samples":members.len(),"members":members})).collect::<Vec<_>>(),"scores":scores,
        "notes":["Groups are an explicit exhaustive partition supplied by the researcher; labels are not inferred or verified demographic or ancestry assignments.",
        "Every group retains the original whole-cohort common-callable term universe. Grouping never rescans genotypes or recovers excluded terms.",
        "Percentile movement and cutoff membership use the same full target baseline distribution for every group. Pearson and Spearman correlations use pairs within each group.",
        "External-reference shifts retain the authored frozen reference distribution. Source-cohort shared sample counts describe the whole input cohort, not overlap within a particular group.",
        "These are descriptive masking-sensitivity statistics, not predictive accuracy, calibration, ancestry normalization, population uncertainty or a causal comparison between groups.",
        "Small or constant groups may have unavailable correlations. Null metrics retain their reasons; group sizes and percentile denominators remain explicit."]}),
    )
}

fn escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}
fn cell(value: &Value) -> String {
    if value.is_null() {
        String::new()
    } else {
        value.as_str().map(str::to_owned).unwrap_or_else(|| value.to_string())
    }
}
fn write_tables(out: &Path, report: &Value) -> Result<BTreeMap<String, String>> {
    let mut tables: BTreeMap<&str, (Vec<&str>, Vec<Vec<String>>)> = BTreeMap::from([
        (
            "metrics",
            (
                vec![
                    "pgs_id",
                    "scenario",
                    "group",
                    "n",
                    "metric",
                    "value",
                    "unavailable_reason",
                ],
                Vec::new(),
            ),
        ),
        (
            "cutoffs",
            (
                vec![
                    "pgs_id",
                    "scenario",
                    "group",
                    "n",
                    "cutoff",
                    "up",
                    "down",
                    "unchanged",
                    "evaluable",
                    "changed_fraction",
                    "unavailable_reason",
                ],
                Vec::new(),
            ),
        ),
        (
            "references",
            (
                vec![
                    "pgs_id",
                    "scenario",
                    "group",
                    "n",
                    "reference",
                    "reference_samples",
                    "source_cohort_shared_sample_ids",
                    "evaluable",
                    "unavailable",
                    "mean_baseline_percentile",
                    "mean_scenario_percentile",
                    "mean_absolute_shift",
                    "max_absolute_shift",
                ],
                Vec::new(),
            ),
        ),
    ]);
    for score in array(report, "scores")? {
        for row in array(score, "scenarios")? {
            let base = vec![
                cell(&score["pgs_id"]),
                cell(&row["scenario"]),
                cell(&row["group"]),
                cell(&row["metrics"]["n"]),
            ];
            for metric in [
                "mae",
                "rmse",
                "max_abs_change",
                "pearson",
                "spearman",
                "mean_abs_percentile_shift",
                "max_abs_percentile_shift",
            ] {
                let mut fields = base.clone();
                fields.extend([
                    metric.into(),
                    cell(&row["metrics"][metric]),
                    cell(&row["metrics"]["unavailable_reasons"][metric]),
                ]);
                tables.get_mut("metrics").unwrap().1.push(fields);
            }
            for crossing in array(&row["metrics"], "crossings")? {
                let mut fields = base.clone();
                fields.extend(["cutoff", "up", "down", "unchanged", "evaluable"].map(|k| cell(&crossing[k])));
                let n = crossing["evaluable"].as_u64().unwrap();
                fields.push(if n == 0 {
                    String::new()
                } else {
                    ((crossing["up"].as_u64().unwrap() + crossing["down"].as_u64().unwrap()) as f64 / n as f64)
                        .to_string()
                });
                if n == 0 {
                    for value in &mut fields[5..8] {
                        value.clear();
                    }
                }
                fields.push(cell(&row["metrics"]["unavailable_reasons"]["crossings"]));
                tables.get_mut("cutoffs").unwrap().1.push(fields);
            }
            for reference in array(row, "reference_comparisons")? {
                let mut fields = base.clone();
                fields.extend(
                    [
                        "id",
                        "reference_samples",
                        "source_cohort_shared_sample_ids",
                        "evaluable",
                        "unavailable",
                        "mean_baseline_percentile",
                        "mean_scenario_percentile",
                        "mean_absolute_shift",
                        "max_absolute_shift",
                    ]
                    .map(|k| cell(&reference[k])),
                );
                tables.get_mut("references").unwrap().1.push(fields);
            }
        }
    }
    let mut html = String::from(
        "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>PGS robustness by group</title><style>body{font:16px system-ui;max-width:1200px;margin:2rem auto;padding:1rem;color:#182c3b}table{border-collapse:collapse;font-size:14px;min-width:100%}th,td{padding:.4rem;border-bottom:1px solid #ddd;text-align:left}th{background:#edf2f5}.scroll{overflow-x:auto}a{color:#12619b}.filters{display:flex;gap:1rem;flex-wrap:wrap;padding:1rem;background:#edf2f5}select{font:inherit;padding:.35rem}label{display:grid;gap:.3rem}</style><h1>PGS robustness by group</h1>",
    );
    for note in array(report, "notes")? {
        html.push_str(&format!("<p>{}</p>", escape(note.as_str().unwrap())));
    }
    html.push_str("<h2>Participant groups</h2><ul>");
    for group in array(report, "groups")? {
        html.push_str(&format!(
            "<li>{}: {} participants</li>",
            escape(string(group, "id")?),
            group["samples"]
        ));
    }
    html.push_str("</ul><p>Download <a href=\"metrics.tsv\">score metrics</a>, <a href=\"cutoffs.tsv\">cutoff transitions</a>, <a href=\"references.tsv\">reference shifts</a> or <a href=\"groups.json\">complete JSON with provenance</a>. Blank values are unavailable.</p>");
    html.push_str("<div class=\"filters\"><label>Participant group<select id=\"group-filter\"><option value=\"\">All groups</option>");
    for group in array(report, "groups")? {
        let id = escape(string(group, "id")?);
        html.push_str(&format!("<option value=\"{id}\">{id}</option>"));
    }
    html.push_str(
        "</select></label><label>Scenario<select id=\"scenario-filter\"><option value=\"\">All scenarios</option>",
    );
    let mut scenario_ids = BTreeSet::new();
    for score in array(report, "scores")? {
        for row in array(score, "scenarios")? {
            scenario_ids.insert(string(row, "scenario")?);
        }
    }
    for id in scenario_ids {
        let id = escape(id);
        html.push_str(&format!("<option value=\"{id}\">{id}</option>"));
    }
    html.push_str("</select></label></div><p id=\"visible-count\" aria-live=\"polite\"></p>");
    let mut hashes = BTreeMap::new();
    for name in ["metrics", "cutoffs", "references"] {
        let (headers, rows) = tables.remove(name).unwrap();
        let title = match name {
            "metrics" => "Score metrics",
            "cutoffs" => "Cutoff transitions",
            _ => "Frozen-reference shifts",
        };
        let mut tsv = headers.join("\t") + "\n";
        html.push_str(&format!("<h2>{title}</h2><div class=\"scroll\"><table><thead><tr>"));
        for header in &headers {
            html.push_str(&format!("<th>{header}</th>"));
        }
        html.push_str("</tr></thead><tbody>");
        for fields in rows {
            tsv.push_str(&fields.join("\t"));
            tsv.push('\n');
            html.push_str(&format!(
                "<tr data-group=\"{}\" data-scenario=\"{}\">",
                escape(&fields[2]),
                escape(&fields[1])
            ));
            for (header, field) in headers.iter().zip(fields) {
                let rounded = if [
                    "value",
                    "changed_fraction",
                    "mean_baseline_percentile",
                    "mean_scenario_percentile",
                    "mean_absolute_shift",
                    "max_absolute_shift",
                ]
                .contains(header)
                {
                    field.parse::<f64>().ok().map(|x| {
                        if x != 0.0 && (x.abs() < 1e-4 || x.abs() >= 1e6) {
                            format!("{x:.5e}")
                        } else {
                            format!("{x:.5}")
                        }
                    })
                } else {
                    None
                };
                html.push_str(&format!(
                    "<td title=\"{}\">{}</td>",
                    escape(&field),
                    escape(rounded.as_deref().unwrap_or(&field))
                ));
            }
            html.push_str("</tr>");
        }
        html.push_str("</tbody></table></div>");
        let filename = format!("{name}.tsv");
        let path = out.join(&filename);
        std::fs::write(&path, tsv).map_err(Error::io(&path))?;
        hashes.insert(filename, crate::digest::file_sha256(&path)?);
    }
    html.push_str(r#"<p>Displayed decimal metrics are rounded; TSV and JSON retain full precision. All rows remain available without JavaScript.</p><script>
const group=document.getElementById('group-filter'),scenario=document.getElementById('scenario-filter');
function filter(){let shown=0;const rows=document.querySelectorAll('tr[data-group]');for(const row of rows){row.hidden=Boolean((group.value&&row.dataset.group!==group.value)||(scenario.value&&row.dataset.scenario!==scenario.value));if(!row.hidden)shown++;}document.getElementById('visible-count').textContent=`Showing ${shown} of ${rows.length} metric, cutoff and reference rows.`;}
group.addEventListener('change',filter);scenario.addEventListener('change',filter);filter();
</script></html>"#);
    let path = out.join("groups.html");
    std::fs::write(&path, html).map_err(Error::io(&path))?;
    hashes.insert("groups.html".into(), crate::digest::file_sha256(&path)?);
    Ok(hashes)
}

pub(super) fn analyze(report: &Value, results: &Path, mapping: &Path) -> Result<Value> {
    let assignments = load_groups(mapping)?;
    let rows = load_scores(&results.join("samples.tsv"))?;
    let references = load_references(&results.join("references.tsv"))?;
    summarize(report, &rows, &references, &assignments)
}

pub fn run(args: &GroupsArgs) -> Result<()> {
    let report_path = args.results.join("report.json");
    let mut provenance = BTreeMap::new();
    for path in [&report_path, &args.groups] {
        provenance.insert(path.clone(), crate::digest::file_sha256(path)?);
    }
    let report: Value = serde_json::from_slice(&std::fs::read(&report_path).map_err(Error::io(&report_path))?)
        .map_err(|e| Error::Invalid(e.to_string()))?;
    super::verify_report_outputs(&args.results, &report, &mut provenance)?;
    let mut grouped = analyze(&report, &args.results, &args.groups)?;
    grouped["source_sha256"] = serde_json::to_value(&provenance).map_err(|e| Error::Invalid(e.to_string()))?;
    grouped["pgsum_version"] = json!(env!("CARGO_PKG_VERSION"));
    grouped["executable_sha256"] = json!(crate::digest::file_sha256(
        &std::env::current_exe().map_err(Error::io("<executable>"))?
    )?);
    std::fs::create_dir(&args.out).map_err(Error::io(&args.out))?;
    grouped["outputs"] = json!(write_tables(&args.out, &grouped)?);
    for (path, digest) in &provenance {
        if crate::digest::file_sha256(path)? != *digest {
            return invalid!("input changed during group analysis");
        }
    }
    let pending = args.out.join(".groups.json.pending");
    super::write_json(&pending, &grouped)?;
    std::fs::rename(&pending, args.out.join("groups.json")).map_err(Error::io(&pending))?;
    Ok(())
}
