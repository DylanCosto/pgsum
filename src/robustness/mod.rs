//! Controlled score-reconstruction experiments. These diagnostics do not estimate phenotype prediction.
pub mod config;
pub mod groups;
pub mod planner;
mod preflight;
mod progress;
pub use progress::ProgressFormat;
mod report;
pub mod stats;
pub mod study;

use std::collections::{BTreeMap, BTreeSet};
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::rc::Rc;

use clap::Args;
use num_bigint::BigInt;
use serde::Serialize;
use sha2::{Digest, Sha256};

use crate::cohort::{CohortReader, CohortTable};
use crate::decimal::Decimal;
use crate::frequencies::FrequencyTable;
use crate::genotypes::GenotypeTable;
use crate::pack::{Pack, TermRecord};
use crate::score::{Contribution, ExactSum, Options, Plan};
use crate::{Error, Result, invalid};
use config::{Bed, Experiment, PreparedMasks, RandomDrawCache};

#[derive(Args, Debug)]
pub struct RobustnessArgs {
    /// JSON experiment manifest. Relative paths resolve beside this file.
    #[arg(long)]
    pub experiment: PathBuf,
    /// New directory for the report and auditable tables; existing directories are never overwritten.
    #[arg(long)]
    pub out: PathBuf,
    /// Validate experiment inputs and report planned dimensions without scoring.
    #[arg(long)]
    pub preflight: bool,
    /// Also write every changed term/sample/scenario to terms.tsv.zst (can be large).
    #[arg(long)]
    pub terms: bool,
    /// Progress on stderr: readable updates, JSON Lines for workflow tools, or none.
    #[arg(long, value_enum, default_value = "human")]
    pub progress: ProgressFormat,
    /// Minimum seconds between scoring updates (checked between blocks of 256 terms).
    #[arg(long, default_value_t = 30, value_parser = clap::value_parser!(u64).range(1..))]
    pub progress_interval: u64,
}

enum Source {
    Cohort(Box<CohortTable>),
    Single(Box<GenotypeTable>),
}
enum Reader<'a> {
    Cohort(&'a CohortTable, CohortReader<'a>),
    Single(&'a GenotypeTable),
}
impl Source {
    fn open(path: &Path) -> Result<Self> {
        let mut magic = [0; 8];
        File::open(path)
            .and_then(|mut f| f.read_exact(&mut magic))
            .map_err(Error::io(path))?;
        let source = if &magic == crate::cohort::MAGIC {
            Self::Cohort(Box::new(CohortTable::open(path)?))
        } else if &magic == crate::genotypes::MAGIC {
            Self::Single(Box::new(GenotypeTable::open(path)?))
        } else {
            return invalid!(
                "{}: expected an extracted .pgsc or .pgsg genotype table",
                path.display()
            );
        };
        let names = source.samples();
        let mut seen = BTreeSet::new();
        if names
            .iter()
            .any(|s| s.is_empty() || s.contains(['\t', '\n', '\r']) || !seen.insert(s))
        {
            return invalid!("{}: empty, duplicate or invalid sample identities", path.display());
        }
        Ok(source)
    }
    fn samples(&self) -> Vec<String> {
        match self {
            Self::Cohort(t) => t.header.samples.clone(),
            Self::Single(t) => vec![t.header.sample.sample_id.clone()],
        }
    }
    fn check(&self, pack: &Pack) -> Result<()> {
        let (reference, packs) = match self {
            Self::Cohort(t) => (&t.header.reference, &t.header.packs),
            Self::Single(t) => (&t.header.reference, &t.header.packs),
        };
        if reference.fasta_sha256 != pack.header.reference.fasta_sha256
            || reference.fai_sha256 != pack.header.reference.fai_sha256
        {
            return invalid!(
                "{}: reference identity differs from genotype extraction",
                pack.header.pgs_id
            );
        }
        if !packs
            .iter()
            .any(|p| p.pgs_id == pack.header.pgs_id && p.records_sha256 == pack.header.records_sha256)
        {
            return invalid!("genotype table was not extracted with this {} pack", pack.header.pgs_id);
        }
        Ok(())
    }
    fn reader(&self) -> Reader<'_> {
        match self {
            Self::Cohort(t) => Reader::Cohort(t, t.reader()),
            Self::Single(t) => Reader::Single(t),
        }
    }
    fn metadata(&self) -> serde_json::Value {
        match self {
            Self::Cohort(t) => serde_json::to_value(&t.header).expect("serializable header"),
            Self::Single(t) => serde_json::to_value(&t.header).expect("serializable header"),
        }
    }
}
impl Reader<'_> {
    fn contributions(&mut self, term: &TermRecord, options: &Options) -> Result<Vec<Option<Contribution>>> {
        match self {
            Self::Single(t) => Ok(vec![crate::score::outcome(term, t, options)?.contribution]),
            Self::Cohort(table, reader) => {
                let n = table.header.samples.len();
                let (key, alt) = match crate::score::plan(term, options).0 {
                    Plan::Unscorable(_) => return Ok(vec![None; n]),
                    Plan::Snv { key, effect_is_alt } => (key, effect_is_alt),
                    Plan::Sequence { variant, effect_is_alt } => (
                        table
                            .sequence_key(term.contig, &variant)
                            .ok_or_else(|| Error::Invalid("missing extracted sequence target".into()))?,
                        effect_is_alt,
                    ),
                };
                let row = reader
                    .row(key)?
                    .ok_or_else(|| Error::Invalid("missing extracted genotype target".into()))?;
                let mut result: Vec<_> = (0..n)
                    .map(|s| {
                        let dosage = crate::cohort::code(&row, s);
                        (dosage != crate::cohort::MISSING)
                            .then(|| crate::score::term_contribution(term, if alt { dosage } else { 2 - dosage }))
                            .flatten()
                    })
                    .collect();
                for (s, m) in row.measurements() {
                    result[*s] = m.contribution(term.model, &term.weights, alt);
                }
                Ok(result)
            }
        }
    }
}

struct ScenarioRuntime {
    masks: PreparedMasks,
    source: Option<Source>,
    indices: Vec<usize>,
    fill: Option<Rc<FrequencyTable>>,
}
struct ReferenceRuntime {
    id: String,
    source: Option<Source>,
    indices: Vec<usize>,
    samples: Vec<String>,
}

#[derive(Serialize)]
pub struct Report {
    pub schema: &'static str,
    pub pgsum_version: &'static str,
    pub executable_sha256: String,
    pub worker_threads: usize,
    pub full_term_audit: bool,
    pub experiment: Experiment,
    pub input_sha256: BTreeMap<String, String>,
    pub source_metadata: BTreeMap<String, serde_json::Value>,
    pub mask_algorithm: &'static str,
    pub notes: Vec<String>,
    pub scores: Vec<ScoreReport>,
    pub outputs: BTreeMap<String, String>,
}
#[derive(Serialize)]
pub struct ScoreReport {
    pub pgs_id: String,
    pub total_model_terms: usize,
    pub baseline_terms: usize,
    pub baseline_term_fraction: f64,
    pub baseline_weight_fraction: Option<f64>,
    pub excluded_terms: BTreeMap<String, usize>,
    pub baseline_term_sha256: String,
    pub inventory_consistent: bool,
    pub matches_publication: Option<bool>,
    pub license: Option<String>,
    pub catalog_metadata: serde_json::Value,
    pub scenarios: Vec<ScenarioReport>,
}
#[derive(Serialize)]
pub struct ScenarioReport {
    pub id: String,
    pub masked_observations: u64,
    pub unavailable_alternate_observations: u64,
    pub filled_observations: u64,
    pub unfilled_observations: u64,
    pub fill_omissions: BTreeMap<String, u64>,
    pub metrics: stats::Summary,
    pub reference_comparisons: Vec<ReferenceSummary>,
}
#[derive(Serialize)]
pub struct ReferenceSummary {
    pub id: String,
    pub reference_samples: usize,
    pub shared_sample_ids: usize,
    pub evaluable: usize,
    pub mean_baseline_percentile: Option<f64>,
    pub mean_scenario_percentile: Option<f64>,
    pub mean_absolute_shift: Option<f64>,
    pub max_absolute_shift: Option<f64>,
}
#[derive(Clone, Serialize)]
struct Change {
    ordinal: usize,
    contig: String,
    pos1: u32,
    source_pos1: u32,
    effect_allele: Option<String>,
    other_allele: Option<String>,
    reason: String,
    baseline: String,
    scenario: String,
    delta: String,
    #[serde(skip)]
    delta_exact: Decimal,
    #[serde(skip)]
    magnitude: Option<Decimal>,
}
struct Accumulator {
    sums: Vec<ExactSum>,
    deltas: Vec<ExactSum>,
    retained: Vec<u64>,
    masked: Vec<u64>,
    unavailable: Vec<u64>,
    filled: Vec<u64>,
    unfilled: Vec<u64>,
    fill_omissions: BTreeMap<String, u64>,
    changed: Vec<u64>,
    top: Vec<Vec<Change>>,
    regions: Vec<Vec<ExactSum>>,
}
impl Accumulator {
    fn new(n: usize, regions: usize) -> Self {
        Self {
            sums: vec![ExactSum::default(); n],
            deltas: vec![ExactSum::default(); n],
            retained: vec![0; n],
            masked: vec![0; n],
            unavailable: vec![0; n],
            filled: vec![0; n],
            unfilled: vec![0; n],
            fill_omissions: BTreeMap::new(),
            changed: vec![0; n],
            top: vec![Vec::new(); n],
            regions: vec![vec![ExactSum::default(); n]; regions],
        }
    }
}

fn decimal_contribution(d: &Decimal) -> Contribution {
    Contribution {
        negative: d.negative,
        coefficient: crate::score::Coefficient::Big(d.coefficient.parse::<BigInt>().expect("decimal")),
        exponent: d.exponent,
    }
}
fn difference(a: &Decimal, b: &Decimal) -> Decimal {
    let mut sum = ExactSum::default();
    sum.add(&decimal_contribution(a));
    sum.add(&decimal_contribution(b).negated());
    sum.finish()
}
fn zero() -> Contribution {
    Contribution {
        negative: false,
        coefficient: crate::score::Coefficient::Small(0),
        exponent: 0,
    }
}
fn index_samples(names: &[String], selected: Option<&[String]>) -> Result<Vec<usize>> {
    let Some(selected) = selected else {
        return Ok((0..names.len()).collect());
    };
    let mut seen = BTreeSet::new();
    selected
        .iter()
        .map(|s| {
            if !seen.insert(s) {
                return invalid!("duplicate selected sample {s}");
            }
            names
                .iter()
                .position(|x| x == s)
                .ok_or_else(|| Error::Invalid(format!("sample {s} is absent from genotype table")))
        })
        .collect()
}
fn write_json(path: &Path, value: &impl Serialize) -> Result<()> {
    let mut out = BufWriter::new(File::create(path).map_err(Error::io(path))?);
    serde_json::to_writer_pretty(&mut out, value).map_err(|e| Error::Invalid(e.to_string()))?;
    out.write_all(b"\n").and_then(|_| out.flush()).map_err(Error::io(path))
}
fn checked_text(value: &str) -> Result<()> {
    if value.is_empty() || value.contains(['\t', '\n', '\r']) {
        return invalid!("invalid tabular identity {value:?}");
    }
    Ok(())
}
fn optional(v: Option<f64>) -> String {
    v.map(|x| x.to_string()).unwrap_or_default()
}

struct Prepared {
    experiment: Experiment,
    baseline: Source,
    samples: Vec<String>,
    scenarios: Vec<ScenarioRuntime>,
    references: Vec<ReferenceRuntime>,
    regions: Vec<Bed>,
    metadata: BTreeMap<String, serde_json::Value>,
    input_sha256: BTreeMap<String, String>,
    models: Vec<serde_json::Value>,
}

fn prepare(path: &Path) -> Result<Prepared> {
    let manifest_hash = crate::digest::file_sha256(path)?;
    let experiment = config::load(path)?;
    let paths = preflight::input_paths(path, &experiment);
    let input_sha256: BTreeMap<String, String> = paths
        .iter()
        .map(|p| Ok((p.display().to_string(), crate::digest::file_sha256(p)?)))
        .collect::<Result<_>>()?;
    if input_sha256[&path.display().to_string()] != manifest_hash {
        return invalid!("manifest changed during preparation");
    }
    let baseline = Source::open(&experiment.genotypes)?;
    let samples = baseline.samples();
    let n = samples.len();
    let mut metadata = BTreeMap::from([(experiment.genotypes.display().to_string(), baseline.metadata())]);
    let mut scenarios = Vec::new();
    // Decoded frequency blocks are immutable but potentially large. A rate/seed sweep
    // reuses each authored input across many scenarios; keep one table per resolved path.
    let mut frequency_tables: BTreeMap<PathBuf, Rc<FrequencyTable>> = BTreeMap::new();
    for scenario in &experiment.scenarios {
        let source = scenario.genotypes.as_deref().map(Source::open).transpose()?;
        let indices = if let Some(t) = &source {
            if t.samples().len() != n {
                return invalid!(
                    "scenario {} must contain exactly the baseline sample identities",
                    scenario.id
                );
            }
            index_samples(&t.samples(), Some(&samples))?
        } else {
            (0..n).collect()
        };
        if let (Some(path), Some(t)) = (&scenario.genotypes, &source) {
            metadata.insert(path.display().to_string(), t.metadata());
        }
        let fill = scenario
            .fill_frequencies
            .as_deref()
            .map(|path| -> Result<Rc<FrequencyTable>> {
                match frequency_tables.entry(path.to_owned()) {
                    std::collections::btree_map::Entry::Occupied(entry) => Ok(Rc::clone(entry.get())),
                    std::collections::btree_map::Entry::Vacant(entry) => {
                        Ok(Rc::clone(entry.insert(Rc::new(FrequencyTable::open(path)?))))
                    }
                }
            })
            .transpose()?;
        if let (Some(path), Some(table)) = (&scenario.fill_frequencies, &fill) {
            metadata.insert(
                path.display().to_string(),
                serde_json::to_value(&table.header).expect("serializable header"),
            );
        }
        scenarios.push(ScenarioRuntime {
            masks: PreparedMasks::load(&scenario.masks)?,
            source,
            indices,
            fill,
        });
    }
    let mut references = Vec::new();
    for reference in &experiment.reference_sets {
        let source = reference.genotypes.as_deref().map(Source::open).transpose()?;
        let names = source.as_ref().map_or_else(|| samples.clone(), Source::samples);
        let indices = index_samples(&names, reference.samples.as_deref())?;
        if indices.is_empty() {
            return invalid!("reference {} selects no samples", reference.id);
        }
        if let (Some(path), Some(t)) = (&reference.genotypes, &source) {
            metadata.insert(path.display().to_string(), t.metadata());
        }
        references.push(ReferenceRuntime {
            id: reference.id.clone(),
            samples: indices.iter().map(|i| names[*i].clone()).collect(),
            source,
            indices,
        });
    }
    let regions: Vec<_> = experiment
        .regions
        .iter()
        .map(|r| Bed::load(&r.path))
        .collect::<Result<_>>()?;
    // Open/validate packs before creating any outputs. The experiment never repairs its input files.
    let mut ids = BTreeSet::new();
    let mut models = Vec::new();
    for path in &experiment.packs {
        let pack = Pack::open(path)?;
        checked_text(&pack.header.pgs_id)?;
        baseline.check(&pack)?;
        for source in scenarios
            .iter()
            .filter_map(|s| s.source.as_ref())
            .chain(references.iter().filter_map(|r| r.source.as_ref()))
        {
            source.check(&pack)?;
        }
        models.push(serde_json::json!({"pgs_id":pack.header.pgs_id,"terms":pack.columns.len(),
            "inventory":pack.header.inventory,"records_sha256":pack.header.records_sha256,"reference":pack.header.reference}));
        if !ids.insert(pack.header.pgs_id) {
            return invalid!("duplicate score IDs in experiment");
        }
    }
    Ok(Prepared {
        experiment,
        baseline,
        samples,
        scenarios,
        references,
        regions,
        metadata,
        input_sha256,
        models,
    })
}

fn verify_inputs(inputs: &BTreeMap<String, String>) -> Result<()> {
    for (path, hash) in inputs {
        if crate::digest::file_sha256(Path::new(path))? != *hash {
            return invalid!("input changed during experiment: {path}");
        }
    }
    Ok(())
}

/// All output is written under a newly created directory. report.json is the completion marker.
pub fn run(args: &RobustnessArgs) -> Result<()> {
    let mut progress = progress::Progress::new(args.progress, args.progress_interval);
    progress.stage("preparing");
    let prepared = prepare(&args.experiment)?;
    progress.stage("fingerprinting");
    verify_inputs(&prepared.input_sha256)?;
    if args.preflight {
        preflight::write(args, &prepared)?;
        progress.stage("preflight_complete");
        return Ok(());
    }
    let Prepared {
        experiment,
        baseline,
        samples,
        scenarios,
        references,
        regions,
        metadata,
        input_sha256,
        ..
    } = prepared;
    let n = samples.len();
    let options = Options::from(&experiment.scoring);
    std::fs::create_dir(&args.out).map_err(Error::io(&args.out))?;
    let mut scores = writer(&args.out, "samples.tsv")?;
    writeln!(scores,"pgs_id\tscenario\tsample\tbaseline_terms\tretained_terms\tmasked_terms\tunavailable_alternate_terms\tfilled_terms\tunfilled_terms\tbaseline_score\tscenario_score\tdelta\tbaseline_percentile\tfixed_baseline_percentile\treranked_percentile\tchanged_terms\ttop_terms_delta\tother_terms_delta").map_err(Error::io(&args.out))?;
    let mut changes = writer(&args.out, "top-changes.jsonl")?;
    let mut region_out = writer(&args.out, "regions.tsv")?;
    writeln!(region_out, "pgs_id\tscenario\tsample\tregion\tdelta").map_err(Error::io(&args.out))?;
    let mut reference_out = writer(&args.out, "references.tsv")?;
    writeln!(reference_out,"pgs_id\tscenario\tsample\treference\treference_samples\tshared_sample_ids\tbaseline_percentile\tscenario_percentile\tpercentile_delta").map_err(Error::io(&args.out))?;
    let mut terms = if args.terms {
        Some(zstd::stream::write::Encoder::new(writer(&args.out, "terms.tsv.zst")?, 3).map_err(Error::io(&args.out))?)
    } else {
        None
    };
    if let Some(out) = &mut terms {
        writeln!(
            out,
            "pgs_id\tscenario\tsample\tordinal\tcontig\tpos1\treason\tbaseline\tscenario_contribution\tdelta"
        )
        .map_err(Error::io(&args.out))?;
    }
    let executable = std::env::current_exe().map_err(Error::io("<executable>"))?;
    let mut report=Report {schema:"pgsum-robustness-report-v1",pgsum_version:env!("CARGO_PKG_VERSION"),
        executable_sha256: crate::digest::file_sha256(&executable)?, worker_threads: rayon::current_num_threads(), full_term_audit: args.terms, experiment,
        input_sha256,source_metadata:metadata,mask_algorithm:config::MASK_ALGORITHM,outputs:BTreeMap::new(),scores:Vec::new(),notes:vec![
            "Available-data baseline: one frozen common-callable autosomal term set across baseline targets and all selected reference samples. Excluded terms remain counted; a subset is never a complete published PGS.".into(),
            "Experimental omissions contribute zero. Optional population-frequency fills are hypothetical reconstructions, including below production coverage gates; they never change production scoring rules.".into(),
            "Percentiles are numerical sensitivity diagnostics against frozen baseline distributions. A perturbed incomplete score is not a calibrated population percentile or a risk estimate.".into(),
            "Within-cohort reranking is reported separately. A constant or fewer-than-two reference distribution has no reported percentile. Cutoffs are descriptive, not clinical.".into(),
            "Per-term and regional deltas explain raw-score movement, not additive percentile changes. Overlapping regions must not be summed.".into(),
            "Frequency-source independence is declared by the user, not established by the software. Evaluation genotypes must not contribute to frequency estimation.".into(),
        ]};
    for (score_index, path) in report.experiment.packs.iter().enumerate() {
        let pack = Pack::open(path)?;
        baseline.check(&pack)?;
        for s in &scenarios {
            if let Some(t) = &s.source {
                t.check(&pack)?
            }
        }
        for r in &references {
            if let Some(t) = &r.source {
                t.check(&pack)?
            }
        }
        progress.start_score(
            &pack.header.pgs_id,
            score_index + 1,
            report.experiment.packs.len(),
            n,
            scenarios.len(),
            pack.columns.len(),
        );
        let mut base_reader = baseline.reader();
        let mut scenario_readers: Vec<_> = scenarios
            .iter()
            .map(|s| s.source.as_ref().map(Source::reader))
            .collect();
        let mut reference_readers: Vec<_> = references
            .iter()
            .map(|r| r.source.as_ref().map(Source::reader))
            .collect();
        let mut base_sums = vec![ExactSum::default(); n];
        let mut ref_sums: Vec<_> = references
            .iter()
            .map(|r| vec![ExactSum::default(); r.indices.len()])
            .collect();
        let mut mask_draws = RandomDrawCache::new(&samples);
        let mut accumulators: Vec<_> = scenarios.iter().map(|_| Accumulator::new(n, regions.len())).collect();
        let mut eligible = 0usize;
        let mut excluded = BTreeMap::new();
        let mut universe = Sha256::new();
        let mut total_weight = 0.0;
        let mut baseline_weight = 0.0;
        for (ordinal, item) in pack.terms().enumerate() {
            if ordinal % 256 == 0 {
                progress.checkpoint(ordinal, eligible, false);
            }
            let term = item?;
            let weight = crate::score::term_effect(&term).unwrap_or(0.0);
            total_weight += weight;
            let reason = if term.contig == 0 || term.pos == 0 {
                Some("unresolved_position")
            } else if term.contig > 22 {
                Some("non_autosomal_term")
            } else if matches!(crate::score::plan(&term, &options).0, Plan::Unscorable(_)) {
                Some("model_or_orientation_unscorable")
            } else {
                None
            };
            if let Some(reason) = reason {
                *excluded.entry(reason.to_string()).or_insert(0) += 1;
                continue;
            }
            let base = base_reader.contributions(&term, &options)?;
            let anchor = match crate::score::plan(&term, &options).0 {
                Plan::Sequence { variant, .. } => u32::try_from(variant.pos)
                    .map_err(|_| Error::Invalid("normalized target position exceeds supported coordinates".into()))?,
                _ => term.pos,
            };
            if base.iter().any(Option::is_none) {
                *excluded
                    .entry("baseline_target_not_common_callable".into())
                    .or_insert(0) += 1;
                continue;
            }
            let mut reference_values = Vec::new();
            let mut reference_missing = false;
            for (r, reader) in references.iter().zip(reference_readers.iter_mut()) {
                let full = if let Some(reader) = reader {
                    reader.contributions(&term, &options)?
                } else {
                    base.clone()
                };
                let selected: Vec<_> = r.indices.iter().map(|i| full[*i].clone()).collect();
                reference_missing |= selected.iter().any(Option::is_none);
                reference_values.push(selected);
            }
            if reference_missing {
                *excluded.entry("reference_not_common_callable".into()).or_insert(0) += 1;
                continue;
            }
            eligible += 1;
            baseline_weight += weight;
            universe.update((ordinal as u64 + 1).to_le_bytes());
            universe.update([term.contig]);
            universe.update(anchor.to_le_bytes());
            for (sum, c) in base_sums.iter_mut().zip(&base) {
                sum.add(c.as_ref().expect("eligible"));
            }
            for (sums, values) in ref_sums.iter_mut().zip(&reference_values) {
                for (sum, c) in sums.iter_mut().zip(values) {
                    sum.add(c.as_ref().expect("eligible"));
                }
            }
            for (j, ((runtime, reader), acc)) in scenarios
                .iter()
                .zip(&mut scenario_readers)
                .zip(&mut accumulators)
                .enumerate()
            {
                let available = if let Some(reader) = reader {
                    reader.contributions(&term, &options)?
                } else {
                    base.clone()
                };
                let mut fill = None;
                for s in 0..n {
                    let original = base[s].as_ref().expect("eligible");
                    let hidden = runtime.masks.hidden_cached(&mut mask_draws, term.contig, anchor, s);
                    let alternate = &available[runtime.indices[s]];
                    let mut reason = "alternate_genotype";
                    let replacement = if hidden || alternate.is_none() {
                        if fill.is_none()
                            && let Some(table) = &runtime.fill
                        {
                            fill = Some(crate::fill::resolve(
                                &term,
                                pack.alleles(ordinal),
                                table,
                                crate::score::term_contribution(&term, 2),
                            )?);
                        }
                        if hidden {
                            acc.masked[s] += 1;
                            reason = "masked";
                        } else {
                            acc.unavailable[s] += 1;
                            reason = "alternate_unavailable";
                        }
                        if let Some(crate::fill::Outcome::Filled { contribution, .. }) = &fill {
                            acc.filled[s] += 1;
                            reason = if hidden {
                                "masked_frequency_fill"
                            } else {
                                "alternate_frequency_fill"
                            };
                            contribution.clone()
                        } else {
                            acc.unfilled[s] += 1;
                            if let Some(crate::fill::Outcome::Omitted(why)) = &fill {
                                *acc.fill_omissions.entry((*why).into()).or_default() += 1;
                            }
                            zero()
                        }
                    } else {
                        acc.retained[s] += 1;
                        alternate.as_ref().expect("available").clone()
                    };
                    acc.sums[s].add(&replacement);
                    if replacement == *original {
                        continue;
                    }
                    let delta = difference(&replacement.to_decimal(), &original.to_decimal());
                    if delta.coefficient == "0" {
                        continue;
                    }
                    acc.changed[s] += 1;
                    let delta_c = decimal_contribution(&delta);
                    acc.deltas[s].add(&delta_c);
                    for (r, bed) in regions.iter().enumerate() {
                        if bed.contains(term.contig, anchor) {
                            acc.regions[r][s].add(&delta_c)
                        }
                    }
                    let (effect, other) = pack
                        .alleles(ordinal)
                        .map_or((None, None), |(e, o)| (Some(e.to_owned()), Some(o.to_owned())));
                    let change = Change {
                        ordinal: ordinal + 1,
                        contig: crate::term::CONTIGS[term.contig as usize - 1].into(),
                        pos1: anchor,
                        source_pos1: term.pos,
                        effect_allele: effect,
                        other_allele: other,
                        reason: reason.into(),
                        baseline: original.to_decimal().to_python_string(),
                        scenario: replacement.to_decimal().to_python_string(),
                        delta: delta.to_python_string(),
                        delta_exact: delta.clone(),
                        magnitude: Some(Decimal {
                            negative: false,
                            ..delta
                        }),
                    };
                    if let Some(out) = &mut terms {
                        writeln!(
                            out,
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            pack.header.pgs_id,
                            report.experiment.scenarios[j].id,
                            samples[s],
                            change.ordinal,
                            change.contig,
                            change.pos1,
                            change.reason,
                            change.baseline,
                            change.scenario,
                            change.delta
                        )
                        .map_err(Error::io(&args.out))?;
                    }
                    let top = &mut acc.top[s];
                    top.push(change);
                    top.sort_by(|a, b| {
                        b.magnitude
                            .as_ref()
                            .unwrap()
                            .numeric_cmp(a.magnitude.as_ref().unwrap())
                            .then(a.ordinal.cmp(&b.ordinal))
                    });
                    top.truncate(report.experiment.top_terms);
                }
            }
        }
        progress.checkpoint(pack.columns.len(), eligible, true);
        if eligible == 0 {
            return invalid!(
                "{}: no common-callable autosomal baseline terms; no completed report",
                pack.header.pgs_id
            );
        }
        progress.stage("summarizing");
        let base_values: Vec<_> = base_sums.iter().map(ExactSum::finish).collect();
        let base_pct = stats::percentiles(&base_values, &base_values);
        let reference_values: Vec<Vec<_>> = ref_sums
            .iter()
            .map(|sums| sums.iter().map(ExactSum::finish).collect())
            .collect();
        let mut summaries = Vec::new();
        for (j, acc) in accumulators.iter().enumerate() {
            let scenario = &report.experiment.scenarios[j];
            let values: Vec<_> = acc.sums.iter().map(ExactSum::finish).collect();
            let fixed = stats::percentiles(&values, &base_values);
            let reranked = stats::percentiles(&values, &values);
            for s in 0..n {
                let delta = difference(&values[s], &base_values[s]);
                if delta.numeric_cmp(&acc.deltas[s].finish()) != std::cmp::Ordering::Equal {
                    return invalid!("internal attribution reconciliation failed");
                }
                let mut top_sum = ExactSum::default();
                for c in &acc.top[s] {
                    top_sum.add(&decimal_contribution(&c.delta_exact));
                }
                let other = difference(&delta, &top_sum.finish());
                writeln!(
                    scores,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    pack.header.pgs_id,
                    scenario.id,
                    samples[s],
                    eligible,
                    acc.retained[s],
                    acc.masked[s],
                    acc.unavailable[s],
                    acc.filled[s],
                    acc.unfilled[s],
                    base_values[s].to_python_string(),
                    values[s].to_python_string(),
                    delta.to_python_string(),
                    optional(base_pct[s]),
                    optional(fixed[s]),
                    optional(reranked[s]),
                    acc.changed[s],
                    top_sum.finish().to_python_string(),
                    other.to_python_string()
                )
                .map_err(Error::io(&args.out))?;
                serde_json::to_writer(&mut changes,&serde_json::json!({"pgs_id":pack.header.pgs_id,"scenario":scenario.id,"sample":samples[s],"changed_terms":acc.changed[s],"top":acc.top[s],"top_terms_delta":top_sum.finish().to_python_string(),"other_terms_delta":other.to_python_string(),"delta":delta.to_python_string()})).map_err(|e|Error::Invalid(e.to_string()))?;
                writeln!(changes).map_err(Error::io(&args.out))?;
                for (r, region) in report.experiment.regions.iter().enumerate() {
                    writeln!(
                        region_out,
                        "{}\t{}\t{}\t{}\t{}",
                        pack.header.pgs_id,
                        scenario.id,
                        samples[s],
                        region.id,
                        acc.regions[r][s].finish().to_python_string()
                    )
                    .map_err(Error::io(&args.out))?;
                }
            }
            let mut reference_comparisons = Vec::new();
            for (r, reference) in references.iter().enumerate() {
                let a = stats::percentiles(&base_values, &reference_values[r]);
                let b = stats::percentiles(&values, &reference_values[r]);
                let shared = reference.samples.iter().filter(|x| samples.contains(x)).count();
                let pairs: Vec<_> = a.iter().zip(&b).filter_map(|(a, b)| a.zip(*b)).collect();
                let mean =
                    |values: Vec<f64>| (!values.is_empty()).then(|| values.iter().sum::<f64>() / values.len() as f64);
                reference_comparisons.push(ReferenceSummary {
                    id: reference.id.clone(),
                    reference_samples: reference.indices.len(),
                    shared_sample_ids: shared,
                    evaluable: pairs.len(),
                    mean_baseline_percentile: mean(pairs.iter().map(|(a, _)| *a).collect()),
                    mean_scenario_percentile: mean(pairs.iter().map(|(_, b)| *b).collect()),
                    mean_absolute_shift: mean(pairs.iter().map(|(a, b)| (b - a).abs()).collect()),
                    max_absolute_shift: pairs.iter().map(|(a, b)| (b - a).abs()).reduce(f64::max),
                });
                for s in 0..n {
                    writeln!(
                        reference_out,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        pack.header.pgs_id,
                        scenario.id,
                        samples[s],
                        reference.id,
                        reference.indices.len(),
                        shared,
                        optional(a[s]),
                        optional(b[s]),
                        optional(a[s].zip(b[s]).map(|(a, b)| b - a))
                    )
                    .map_err(Error::io(&args.out))?;
                }
            }
            summaries.push(ScenarioReport {
                id: scenario.id.clone(),
                masked_observations: acc.masked.iter().sum(),
                unavailable_alternate_observations: acc.unavailable.iter().sum(),
                filled_observations: acc.filled.iter().sum(),
                unfilled_observations: acc.unfilled.iter().sum(),
                fill_omissions: acc.fill_omissions.clone(),
                metrics: stats::summarize(&base_values, &values, &report.experiment.cutoffs)?,
                reference_comparisons,
            });
        }
        let weight_fraction = (total_weight > 0.0 && total_weight.is_finite() && baseline_weight.is_finite())
            .then_some(baseline_weight / total_weight);
        report.scores.push(ScoreReport {
            pgs_id: pack.header.pgs_id.clone(),
            total_model_terms: pack.columns.len(),
            baseline_terms: eligible,
            baseline_term_fraction: eligible as f64 / pack.columns.len() as f64,
            baseline_weight_fraction: weight_fraction,
            excluded_terms: excluded,
            baseline_term_sha256: crate::digest::hex(&universe.finalize()),
            inventory_consistent: pack.header.inventory.consistent,
            matches_publication: pack.header.matches_publication,
            license: pack.header.license.clone(),
            catalog_metadata: pack.header.catalog_metadata.clone(),
            scenarios: summaries,
        });
    }
    progress.writing();
    for out in [&mut scores, &mut changes, &mut region_out, &mut reference_out] {
        out.flush().map_err(Error::io(&args.out))?;
    }
    if let Some(out) = terms {
        out.finish().and_then(|mut w| w.flush()).map_err(Error::io(&args.out))?;
    }
    report::write(&args.out.join("report.html"), &report)?;
    for name in [
        "samples.tsv",
        "top-changes.jsonl",
        "regions.tsv",
        "references.tsv",
        "report.html",
    ] {
        report
            .outputs
            .insert(name.into(), crate::digest::file_sha256(&args.out.join(name))?);
    }
    if args.terms {
        report.outputs.insert(
            "terms.tsv.zst".into(),
            crate::digest::file_sha256(&args.out.join("terms.tsv.zst"))?,
        );
    }
    progress.stage("verifying_inputs");
    // Re-hash inputs so a concurrent replacement cannot leave a successful provenance claim.
    for (path, hash) in &report.input_sha256 {
        if crate::digest::file_sha256(Path::new(path))? != *hash {
            return invalid!("input changed during experiment: {path}");
        }
    }
    let marker = args.out.join(".report.json.tmp");
    write_json(&marker, &report)?;
    std::fs::rename(&marker, args.out.join("report.json")).map_err(Error::io(&marker))?;
    progress.complete(&args.out.join("report.html").to_string_lossy());
    Ok(())
}
fn writer(dir: &Path, name: &str) -> Result<BufWriter<File>> {
    let path = dir.join(name);
    Ok(BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .map_err(Error::io(&path))?,
    ))
}

/// Authenticate the completed native output inventory before downstream analysis.
fn verify_report_outputs(
    results: &Path,
    report: &serde_json::Value,
    provenance: &mut BTreeMap<PathBuf, String>,
) -> Result<()> {
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
        let p = results.join(name);
        let actual = crate::digest::file_sha256(&p)?;
        if hash.as_str() != Some(actual.as_str()) {
            return invalid!("report output hash differs: {name}");
        }
        provenance.insert(p, actual);
    }
    Ok(())
}
