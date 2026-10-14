//! Expand paired rate/seed experiments without opening genotype data.
use std::path::PathBuf;

use clap::{Args, ValueEnum};
use serde::Serialize;

use super::config::{self, Experiment, Mask, RandomUnit};
use crate::{Error, Result, invalid};

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum Unit {
    Site,
    SampleSite,
}

#[derive(Args, Debug)]
pub struct PlannerArgs {
    /// Valid experiment whose scenarios define the policies to compare (no random masks).
    #[arg(long)]
    pub template: PathBuf,
    /// Missingness rates, comma-separated. The same seed nests masks across rates.
    #[arg(long, value_delimiter = ',', required = true)]
    pub rates: Vec<f64>,
    /// Distinct reproducible seeds, comma-separated; reused across every policy.
    #[arg(long, value_delimiter = ',', required = true)]
    pub seeds: Vec<u64>,
    #[arg(long, value_enum, default_value = "sample-site")]
    pub mask_unit: Unit,
    /// New directory for plan.json and the runnable experiment.json.
    #[arg(long)]
    pub out: PathBuf,
    /// Reject larger expansions before allocating scenarios (includes original controls).
    #[arg(long, default_value_t = 10000)]
    pub max_scenarios: usize,
}

#[derive(Serialize)]
struct Arm {
    scenario: String,
    policy: String,
    rate: Option<f64>,
    seed: Option<u64>,
}

fn expand(mut experiment: Experiment, args: &PlannerArgs) -> Result<(Experiment, Vec<Arm>)> {
    if args.rates.is_empty() || args.seeds.is_empty() {
        return invalid!("provide at least one rate and one seed");
    }
    for (i, rate) in args.rates.iter().enumerate() {
        if !rate.is_finite() || !(0.0..=1.0).contains(rate) || args.rates[..i].contains(rate) {
            return invalid!("rates must be distinct finite numbers between 0 and 1");
        }
    }
    for (i, seed) in args.seeds.iter().enumerate() {
        if args.seeds[..i].contains(seed) {
            return invalid!("duplicate seed {seed}: repeated seeds are not independent replicates");
        }
    }
    let count = args
        .rates
        .len()
        .checked_mul(args.seeds.len())
        .and_then(|n| n.checked_add(1))
        .and_then(|n| n.checked_mul(experiment.scenarios.len()));
    if count.is_none_or(|n| n > args.max_scenarios) {
        return invalid!("expanded scenario count exceeds max_scenarios {}", args.max_scenarios);
    }
    if experiment
        .scenarios
        .iter()
        .any(|s| s.masks.iter().any(|m| matches!(m, Mask::Random { .. })))
    {
        return invalid!(
            "template policies must not contain random masks; use BED masks or alternative inputs, then supply rates/seeds here"
        );
    }
    let policies = experiment.scenarios.clone();
    let mut arms: Vec<_> = policies
        .iter()
        .map(|s| Arm {
            scenario: s.id.clone(),
            policy: s.id.clone(),
            rate: None,
            seed: None,
        })
        .collect();
    for rate in &args.rates {
        for seed in &args.seeds {
            for policy in &policies {
                let mut scenario = policy.clone();
                scenario.id = format!("{}/rate={rate}/seed={seed}", policy.id);
                scenario.masks.push(Mask::Random {
                    rate: *rate,
                    seed: *seed,
                    unit: match args.mask_unit {
                        Unit::Site => RandomUnit::Site,
                        Unit::SampleSite => RandomUnit::SampleSite,
                    },
                });
                arms.push(Arm {
                    scenario: scenario.id.clone(),
                    policy: policy.id.clone(),
                    rate: Some(*rate),
                    seed: Some(*seed),
                });
                experiment.scenarios.push(scenario);
            }
        }
    }
    experiment.validate()?; // Includes collisions between generated IDs and authored policy IDs.
    Ok((experiment, arms))
}

pub fn run(args: &PlannerArgs) -> Result<()> {
    let before = crate::digest::file_sha256(&args.template)?;
    let (experiment, arms) = expand(config::load(&args.template)?, args)?;
    let mut bytes = serde_json::to_vec_pretty(&experiment).map_err(|e| Error::Invalid(e.to_string()))?;
    bytes.push(b'\n');
    if crate::digest::file_sha256(&args.template)? != before {
        return invalid!("template changed while planning");
    }
    std::fs::create_dir(&args.out).map_err(Error::io(&args.out))?;
    let temporary = args.out.join("experiment.json.pending");
    std::fs::write(&temporary, bytes).map_err(Error::io(&temporary))?;
    let plan = serde_json::json!({
        "schema": "pgsum-robustness-plan-v1", "pgsum_version": env!("CARGO_PKG_VERSION"),
        "template_sha256": before, "experiment_sha256": crate::digest::file_sha256(&temporary)?,
        "mask_algorithm": config::MASK_ALGORITHM, "scenarios": arms,
        "notes": ["Original policies are retained as controls without additional random omissions.",
            "Matching seeds pair masks across policies and nest them across rates; policy BED masks can add further omissions.",
            "Replicates describe masking sensitivity, not population uncertainty or clinical confidence intervals.",
            "Referenced input paths are absolute; inputs are validated when the experiment runs."]
    });
    let path = args.out.join("plan.json");
    std::fs::write(
        &path,
        serde_json::to_vec_pretty(&plan).map_err(|e| Error::Invalid(e.to_string()))?,
    )
    .map_err(Error::io(&path))?;
    let path = args.out.join("experiment.json");
    std::fs::rename(&temporary, &path).map_err(Error::io(&path))?;
    eprintln!("Planned {} scenarios: {}", experiment.scenarios.len(), path.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> (Experiment, PlannerArgs) {
        let experiment = serde_json::from_value(serde_json::json!({
            "schema":config::EXPERIMENT_SCHEMA,"genotypes":"calls.pgsc","packs":["score.pgsp"],
            "scenarios":[{"id":"omit"},{"id":"fill","fill_frequencies":"independent.pgsf","frequency_source_note":"disjoint reference"}]
        })).unwrap();
        (
            experiment,
            PlannerArgs {
                template: "template.json".into(),
                rates: vec![0.1, 0.5],
                seeds: vec![7, 9],
                mask_unit: Unit::Site,
                out: "out".into(),
                max_scenarios: 10,
            },
        )
    }
    #[test]
    fn expansion_preserves_fill_policy_and_records_every_arm() {
        let (experiment, args) = fixture();
        let (expanded, arms) = expand(experiment, &args).unwrap();
        assert_eq!(expanded.scenarios.len(), 10);
        assert_eq!(arms.len(), 10);
        for (s, a) in expanded.scenarios.iter().zip(&arms) {
            assert_eq!(s.id, a.scenario);
            assert_eq!(s.fill_frequencies.is_some(), a.policy == "fill");
            if a.seed.is_some() {
                assert!(matches!(
                    s.masks[0],
                    Mask::Random {
                        unit: RandomUnit::Site,
                        ..
                    }
                ));
            }
        }
    }
    #[test]
    fn rejects_duplicate_nonfinite_or_oversized_designs() {
        for rates in [
            vec![],
            vec![0.1, 0.1],
            vec![-0.0, 0.0],
            vec![f64::NAN],
            vec![f64::INFINITY],
            vec![-0.1],
            vec![1.1],
        ] {
            let (experiment, mut args) = fixture();
            args.rates = rates;
            assert!(expand(experiment, &args).is_err());
        }
        for seeds in [vec![], vec![7, 7]] {
            let (experiment, mut args) = fixture();
            args.seeds = seeds;
            assert!(expand(experiment, &args).is_err());
        }
        let (experiment, mut args) = fixture();
        args.max_scenarios = 9;
        assert!(expand(experiment, &args).is_err());
    }
    #[test]
    fn rejects_confounded_random_templates_and_generated_id_collisions() {
        let (mut experiment, args) = fixture();
        experiment.scenarios[0].masks.push(Mask::Random {
            rate: 0.2,
            seed: 8,
            unit: RandomUnit::Site,
        });
        assert!(expand(experiment, &args).is_err());
        let (mut experiment, args) = fixture();
        experiment.scenarios[1].id = "omit/rate=0.1/seed=7".into();
        assert!(expand(experiment, &args).is_err());
    }
}
