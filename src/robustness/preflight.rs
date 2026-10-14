//! Configuration/input checks and explicit workload dimensions, without a score scan.
use super::*;
use serde_json::json;

pub(super) fn input_paths(manifest: &Path, experiment: &Experiment) -> BTreeSet<PathBuf> {
    let mut paths = BTreeSet::from([manifest.to_owned(), experiment.genotypes.clone()]);
    paths.extend(experiment.packs.iter().cloned());
    paths.extend(experiment.regions.iter().map(|r| r.path.clone()));
    paths.extend(experiment.reference_sets.iter().filter_map(|r| r.genotypes.clone()));
    for scenario in &experiment.scenarios {
        paths.extend(scenario.genotypes.iter().cloned());
        paths.extend(scenario.fill_frequencies.iter().cloned());
        paths.extend(scenario.masks.iter().filter_map(|m| match m {
            config::Mask::Bed { path, .. } => Some(path.clone()),
            _ => None,
        }));
    }
    paths
}

fn product(values: &[u64]) -> Option<u64> {
    if values.contains(&0) {
        return Some(0);
    }
    values.iter().try_fold(1u64, |a, b| a.checked_mul(*b))
}

pub(super) fn write(args: &RobustnessArgs, prepared: &Prepared) -> Result<()> {
    let n = prepared.samples.len() as u64;
    let arms = prepared.scenarios.len() as u64;
    let models = prepared.models.len() as u64;
    let terms = prepared
        .models
        .iter()
        .try_fold(0u64, |a, model| a.checked_add(model["terms"].as_u64().unwrap()));
    let candidate_cells = terms.and_then(|terms| product(&[terms, n, arms]));
    let target_names: BTreeSet<_> = prepared.samples.iter().collect();
    let references: Vec<_> = prepared
        .references
        .iter()
        .map(|r| {
            let shared: Vec<_> = r.samples.iter().filter(|s| target_names.contains(s)).collect();
            json!({"id":r.id,"samples":r.samples,"sample_count":r.samples.len(),
            "shared_target_sample_ids":shared,"shared_target_sample_count":shared.len()})
        })
        .collect();
    let mut observations = Vec::new();
    if n < 2 {
        observations.push("The target cohort has fewer than two participants; target-distribution percentiles and correlations will be unavailable.".to_owned());
    }
    for reference in &references {
        if reference["sample_count"].as_u64().unwrap() < 2 {
            observations.push(format!(
                "Reference {} has fewer than two participants; its percentiles will be unavailable.",
                reference["id"]
            ));
        }
        if reference["shared_target_sample_count"].as_u64().unwrap() > 0 {
            observations.push(format!(
                "Reference {} shares {} sample IDs with the targets; it is not disjoint by sample identity.",
                reference["id"], reference["shared_target_sample_count"]
            ));
        }
    }
    for model in &prepared.models {
        if model["inventory"]["consistent"] == false {
            observations.push(format!("Model {} has inconsistent declared and actual inventories; inspect its provenance before interpreting scores.", model["pgs_id"]));
        }
    }
    if prepared.scenarios.iter().any(|s| s.fill.is_some()) {
        observations.push("Frequency-source independence is a user declaration. Input checks cannot establish independence from the evaluation participants.".into());
    }
    let limits = vec![
        "Preflight validates the manifest, sample alignment, model payloads, extraction/reference identities and BED syntax. It performs no genotype-to-score scan.",
        "Cohort and frequency data blocks are decoded lazily by scoring; preflight does not validate every block or prove that scoring will succeed. Some older single-sample formats decode during opening.",
        "Common-callable terms, baseline coverage, missingness, model eligibility and score distributions are unknown until scoring. No score, percentile or clinical validity is established.",
        "Candidate term/sample/scenario cells are a structural upper bound before exclusions, not an operation count, runtime estimate, RAM estimate or output byte estimate. Baseline/reference work adds further processing.",
        "Input fingerprints identify file bytes. Inputs must remain immutable during preparation and scoring.",
        "Preflight is a separate check, not a reusable approval or resume checkpoint. Run the original manifest again into a new output directory to score it.",
    ];
    let mut report = json!({"schema":"pgsum-robustness-preflight-v1","status":"input_checks_passed","scoring_performed":false,
        "pgsum_version":env!("CARGO_PKG_VERSION"),"executable_sha256":crate::digest::file_sha256(&std::env::current_exe().map_err(Error::io("<executable>"))?)?,
        "experiment":prepared.experiment,"input_sha256":prepared.input_sha256,"source_metadata":prepared.metadata,
        "target_samples":prepared.samples,"models":prepared.models,"reference_sets":references,
        "dimensions":{"target_samples":n,"models":models,"scenarios":arms,"model_terms":terms,
            "candidate_term_sample_scenario_cells":candidate_cells,
            "samples_tsv_rows_if_complete":product(&[models,n,arms]),
            "reference_tsv_rows_if_complete":product(&[models,n,arms,prepared.references.len() as u64]),
            "region_tsv_rows_if_complete":product(&[models,n,arms,prepared.regions.len() as u64]),
            "full_term_audit_requested":args.terms,"full_term_audit_max_rows":if args.terms {candidate_cells} else {Some(0)},
            "count_overflow_behavior":"Counts that exceed u64 are null, never wrapped or truncated."},
        "observations":observations,"limitations":limits});
    let mut text =
        format!("PGS robustness preflight: input checks passed\n\nTargets: {n}\nModels: {models}\nScenarios: {arms}\n");
    for (label, key) in [
        ("Model terms", "model_terms"),
        (
            "Candidate term/sample/scenario cells",
            "candidate_term_sample_scenario_cells",
        ),
        ("Score rows if complete", "samples_tsv_rows_if_complete"),
        ("Reference rows if complete", "reference_tsv_rows_if_complete"),
        ("Region rows if complete", "region_tsv_rows_if_complete"),
        ("Maximum full-audit rows", "full_term_audit_max_rows"),
    ] {
        text.push_str(&format!("{label}: {}\n", report["dimensions"][key]));
    }
    for model in &prepared.models {
        text.push_str(&format!("\nModel {}: {} terms\n", model["pgs_id"], model["terms"]));
    }
    for reference in &references {
        text.push_str(&format!(
            "Reference {}: {} participants, {} shared target IDs\n",
            reference["id"], reference["sample_count"], reference["shared_target_sample_count"]
        ));
    }
    text.push_str("\nObservations\n");
    if observations.is_empty() {
        text.push_str("No observations from the checks performed. Scoring remains unchecked.\n");
    }
    for note in observations.iter().map(String::as_str).chain(limits.iter().copied()) {
        text.push_str(&format!("- {note}\n"));
    }
    std::fs::create_dir(&args.out).map_err(Error::io(&args.out))?;
    let text_path = args.out.join("preflight.txt");
    std::fs::write(&text_path, text).map_err(Error::io(&text_path))?;
    report["outputs"] = json!({"preflight.txt":crate::digest::file_sha256(&text_path)?});
    verify_inputs(&prepared.input_sha256)?;
    let pending = args.out.join(".preflight.json.pending");
    write_json(&pending, &report)?;
    std::fs::rename(&pending, args.out.join("preflight.json")).map_err(Error::io(&pending))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn dimension_overflow_is_unavailable() {
        assert_eq!(product(&[3, 128, 42]), Some(16128));
        assert_eq!(product(&[u64::MAX, 2]), None);
        assert_eq!(product(&[u64::MAX, 2, 0]), Some(0));
        assert_eq!(product(&[100, 0]), Some(0));
    }
}
