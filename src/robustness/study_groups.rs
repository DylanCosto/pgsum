//! Reuse the verified paired study engine for each explicit participant partition.
use super::*;

pub(super) fn summarize(
    plan: &Plan,
    experiment: &Experiment,
    native: &Value,
    grouped: &Value,
    reference: &str,
) -> Result<Value> {
    // Validate native plan identity before adapting any per-scenario metrics.
    super::summarize(plan, experiment, native, reference)?;
    let groups = array(grouped, "groups")?;
    let source_samples: u64 = groups.iter().map(|g| g["samples"].as_u64().unwrap()).sum();
    let mut scores = Vec::new();
    let mut notes = Vec::new();
    for group in groups {
        let id = text(group, "id")?;
        let mut report = native.clone();
        let mut selected = Vec::new();
        for score in array(grouped, "scores")? {
            let mut selected_score = score.clone();
            selected_score["scenarios"] = json!(array(score, "scenarios")?
                .iter()
                .filter(|s| s["group"] == id)
                .map(|s| json!({"id":s["scenario"],"metrics":s["metrics"],"reference_comparisons":s["reference_comparisons"]}))
                .collect::<Vec<_>>());
            selected.push(selected_score);
        }
        report["scores"] = json!(selected);
        let summary = summarize_scoped(plan, experiment, &report, reference, Some(source_samples))?;
        notes = array(&summary, "notes")?.clone();
        for score in array(&summary, "scores")? {
            let mut score = score.clone();
            score["group"] = json!(id);
            score["source_samples"] = json!(source_samples);
            scores.push(score);
        }
    }
    notes.extend(array(grouped, "notes")?.clone());
    notes.push(json!("Policy differences are paired within each participant group, rate and seed. Groups are never pooled or compared as causal effects. Reference overlap counts remain source-cohort counts. The full-target baseline defines all percentile thresholds."));
    Ok(
        json!({"schema":"pgsum-robustness-group-study-v1","reference_policy":reference,
        "percentile_basis":"full_target_baseline","correlation_basis":"within_group",
        "groups":groups,"scores":scores,"notes":notes}),
    )
}
