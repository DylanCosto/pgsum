//! Reproducible robustness manifests and masks. BED intervals select the one-base,
//! zero-based anchor of a target; they do not expand to an indel's allele span.

use std::collections::{BTreeMap, HashSet};
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::{Error, Result, invalid};

pub const EXPERIMENT_SCHEMA: &str = "pgsum-robustness-experiment-v1";
/// SHA-256 over this NUL-terminated domain, the little-endian seed, contig byte,
/// little-endian 1-based position, unit byte (0 site, 1 sample_site), and, only for
/// sample_site, the little-endian u64 UTF-8 sample length and sample bytes. The
/// first eight digest bytes, read little-endian and shifted right by 11, give a
/// 53-bit integer; dividing by 2^53 gives the draw compared strictly below rate.
pub const MASK_ALGORITHM: &str = "pgsum-robustness-mask-sha256-v1";

fn default_cutoffs() -> Vec<f64> {
    vec![90.0, 95.0]
}

fn default_top_terms() -> usize {
    10
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Experiment {
    pub schema: String,
    pub genotypes: PathBuf,
    pub packs: Vec<PathBuf>,
    #[serde(default)]
    pub scoring: ScoringOptions,
    #[serde(default = "default_cutoffs")]
    pub cutoffs: Vec<f64>,
    pub scenarios: Vec<Scenario>,
    #[serde(default)]
    pub reference_sets: Vec<ReferenceSet>,
    #[serde(default)]
    pub regions: Vec<Region>,
    #[serde(default = "default_top_terms")]
    pub top_terms: usize,
}

#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct ScoringOptions {
    pub allow_inferred_other_allele: bool,
    pub accept_informational_descriptions: bool,
    pub allow_inferred_palindromes: bool,
    pub allow_inferred_indels: bool,
}

impl From<&ScoringOptions> for crate::score::Options {
    fn from(options: &ScoringOptions) -> Self {
        Self {
            allow_inferred_other_allele: options.allow_inferred_other_allele,
            accept_informational_descriptions: options.accept_informational_descriptions,
            allow_inferred_palindromes: options.allow_inferred_palindromes,
            allow_inferred_indels: options.allow_inferred_indels,
        }
    }
}

impl From<ScoringOptions> for crate::score::Options {
    fn from(options: ScoringOptions) -> Self {
        Self::from(&options)
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Scenario {
    pub id: String,
    #[serde(default)]
    pub masks: Vec<Mask>,
    #[serde(default)]
    pub fill_frequencies: Option<PathBuf>,
    #[serde(default)]
    pub frequency_source_note: Option<String>,
    /// Optional independently extracted cohort, for example at a stricter GQ.
    #[serde(default)]
    pub genotypes: Option<PathBuf>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceSet {
    pub id: String,
    #[serde(default)]
    pub genotypes: Option<PathBuf>,
    /// None uses every sample; Some requires a nonempty, duplicate-free list.
    #[serde(default)]
    pub samples: Option<Vec<String>>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Region {
    pub id: String,
    pub path: PathBuf,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Mask {
    Random { rate: f64, seed: u64, unit: RandomUnit },
    Bed { path: PathBuf, mode: BedMode },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum RandomUnit {
    SampleSite,
    Site,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum BedMode {
    Keep,
    Drop,
}

fn validate_id(scope: &str, id: &str) -> Result<()> {
    if id.is_empty() || id.trim() != id || id.chars().any(char::is_control) {
        return invalid!("{scope}: IDs must be nonempty, trimmed, and contain no control characters");
    }
    Ok(())
}

fn validate_ids<'a>(scope: &str, ids: impl Iterator<Item = &'a str>) -> Result<()> {
    let mut seen = HashSet::new();
    for id in ids {
        validate_id(scope, id)?;
        if !seen.insert(id) {
            return invalid!("{scope}: duplicate ID {id:?}");
        }
    }
    Ok(())
}

fn validate_path(name: &str, path: &Path) -> Result<()> {
    if path.as_os_str().is_empty() {
        return invalid!("{name}: path must not be empty");
    }
    Ok(())
}

fn validate_mask(mask: &Mask) -> Result<()> {
    match mask {
        Mask::Random { rate, .. } if !rate.is_finite() || !(0.0..=1.0).contains(rate) => {
            invalid!("random mask rate must be finite and between 0 and 1 inclusive")
        }
        Mask::Bed { path, .. } => validate_path("BED mask", path),
        _ => Ok(()),
    }
}

impl Experiment {
    pub fn validate(&self) -> Result<()> {
        if self.schema != EXPERIMENT_SCHEMA {
            return invalid!(
                "unsupported robustness schema {:?}; expected {EXPERIMENT_SCHEMA}",
                self.schema
            );
        }
        validate_path("genotypes", &self.genotypes)?;
        if self.packs.is_empty() {
            return invalid!("robustness experiment needs at least one pack");
        }
        for pack in &self.packs {
            validate_path("pack", pack)?;
        }
        if self.scenarios.is_empty() {
            return invalid!("robustness experiment needs at least one scenario");
        }
        if !(1..=1000).contains(&self.top_terms) {
            return invalid!("top_terms must be between 1 and 1000 inclusive");
        }
        for (i, cutoff) in self.cutoffs.iter().enumerate() {
            if !cutoff.is_finite() || !(0.0..=100.0).contains(cutoff) {
                return invalid!("cutoffs must be finite percentiles between 0 and 100 inclusive");
            }
            if self.cutoffs[..i].contains(cutoff) {
                return invalid!("duplicate percentile cutoff {cutoff}");
            }
        }
        validate_ids("scenarios", self.scenarios.iter().map(|s| s.id.as_str()))?;
        validate_ids("reference_sets", self.reference_sets.iter().map(|r| r.id.as_str()))?;
        validate_ids("regions", self.regions.iter().map(|r| r.id.as_str()))?;
        for scenario in &self.scenarios {
            for mask in &scenario.masks {
                validate_mask(mask)?;
            }
            if let Some(path) = &scenario.genotypes {
                validate_path("scenario genotypes", path)?;
            }
            if let Some(path) = &scenario.fill_frequencies {
                validate_path("fill_frequencies", path)?;
                if scenario
                    .frequency_source_note
                    .as_deref()
                    .is_none_or(|note| note.trim().is_empty())
                {
                    return invalid!(
                        "scenario {:?}: fill_frequencies requires a nonempty frequency_source_note describing the frozen source and evaluation-sample overlap",
                        scenario.id
                    );
                }
            }
        }
        for reference in &self.reference_sets {
            if let Some(path) = &reference.genotypes {
                validate_path("reference genotypes", path)?;
            }
            if let Some(samples) = &reference.samples {
                if samples.is_empty() {
                    return invalid!("reference {:?}: samples must not be empty", reference.id);
                }
                validate_ids("reference samples", samples.iter().map(String::as_str))?;
            }
        }
        for region in &self.regions {
            validate_path("region", &region.path)?;
        }
        Ok(())
    }

    fn resolve_paths(&mut self, directory: &Path) {
        let resolve = |path: &mut PathBuf| {
            if path.is_relative() {
                *path = directory.join(&*path);
            }
        };
        resolve(&mut self.genotypes);
        self.packs.iter_mut().for_each(&resolve);
        for scenario in &mut self.scenarios {
            for mask in &mut scenario.masks {
                if let Mask::Bed { path, .. } = mask {
                    resolve(path);
                }
            }
            if let Some(path) = &mut scenario.fill_frequencies {
                resolve(path);
            }
            if let Some(path) = &mut scenario.genotypes {
                resolve(path);
            }
        }
        for reference in &mut self.reference_sets {
            if let Some(path) = &mut reference.genotypes {
                resolve(path);
            }
        }
        for region in &mut self.regions {
            resolve(&mut region.path);
        }
    }
}

/// Read and validate a JSON manifest, resolving relative paths against its parent
/// directory. Referenced inputs are opened by the experiment, not by this loader.
pub fn load(path: &Path) -> Result<Experiment> {
    let bytes = std::fs::read(path).map_err(Error::io(path))?;
    let mut experiment: Experiment =
        serde_json::from_slice(&bytes).map_err(|error| Error::Invalid(format!("{}: {error}", path.display())))?;
    experiment.validate()?;
    let absolute = if path.is_absolute() {
        path.to_owned()
    } else {
        std::env::current_dir().map_err(Error::io("."))?.join(path)
    };
    let directory = absolute
        .parent()
        .ok_or_else(|| Error::Invalid("manifest has no parent directory".into()))?;
    experiment.resolve_paths(directory);
    Ok(experiment)
}

/// Sorted, merged BED3 intervals, keyed by pgsum's canonical contig code.
/// Additional BED columns are accepted but do not affect anchor selection.
#[derive(Clone, Debug, Default)]
pub struct Bed {
    intervals: BTreeMap<u8, Vec<(u32, u32)>>,
}

impl Bed {
    pub fn load(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path).map_err(Error::io(path))?;
        Self::parse(&text).map_err(|error| Error::Invalid(format!("{}: {error}", path.display())))
    }

    fn parse(text: &str) -> Result<Self> {
        let mut bed = Self::default();
        for (index, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<_> = line.split_whitespace().collect();
            if matches!(fields[0], "track" | "browser") {
                continue;
            }
            let number = index + 1;
            if fields.len() < 3 {
                return invalid!("BED line {number}: expected at least three columns");
            }
            let Some(contig) = crate::term::contig_code_of_name(fields[0]) else {
                return invalid!("BED line {number}: unsupported contig {:?}", fields[0]);
            };
            let coordinate = |value: &str| -> Result<u32> {
                if !value.bytes().all(|byte| byte.is_ascii_digit()) {
                    return invalid!("BED line {number}: coordinates must be unsigned integers");
                }
                value.parse().map_err(|_| {
                    Error::Invalid(format!(
                        "BED line {number}: coordinate {value:?} is outside the supported 0..={} range",
                        u32::MAX
                    ))
                })
            };
            let (start, end) = (coordinate(fields[1])?, coordinate(fields[2])?);
            if start >= end {
                return invalid!("BED line {number}: interval must have start < end (0-based, half-open)");
            }
            bed.intervals.entry(contig).or_default().push((start, end));
        }
        for intervals in bed.intervals.values_mut() {
            intervals.sort_unstable();
            let mut merged: Vec<(u32, u32)> = Vec::with_capacity(intervals.len());
            for &(start, end) in intervals.iter() {
                if let Some(last) = merged.last_mut().filter(|last| start <= last.1) {
                    last.1 = last.1.max(end);
                } else {
                    merged.push((start, end));
                }
            }
            *intervals = merged;
        }
        Ok(bed)
    }

    /// Whether the 1-based target anchor lies in a BED interval. Position zero
    /// and unsupported contigs never match; callers pass normalized pack targets.
    pub fn contains(&self, contig: u8, pos1: u32) -> bool {
        let Some(anchor) = pos1.checked_sub(1) else {
            return false;
        };
        let Some(intervals) = self.intervals.get(&contig) else {
            return false;
        };
        let index = intervals.partition_point(|&(start, _)| start <= anchor);
        index > 0 && anchor < intervals[index - 1].1
    }
}

#[derive(Clone, Debug)]
enum PreparedMask {
    Random { rate: f64, seed: u64, unit: RandomUnit },
    Bed { intervals: Bed, mode: BedMode },
}

/// Masks combine by union: a target is hidden if any mask selects it for removal.
#[derive(Clone, Debug, Default)]
pub struct PreparedMasks {
    masks: Vec<PreparedMask>,
}

impl PreparedMasks {
    pub fn load(masks: &[Mask]) -> Result<Self> {
        let mut prepared = Vec::with_capacity(masks.len());
        for mask in masks {
            validate_mask(mask)?;
            prepared.push(match mask {
                Mask::Random { rate, seed, unit } => PreparedMask::Random {
                    rate: *rate,
                    seed: *seed,
                    unit: *unit,
                },
                Mask::Bed { path, mode } => PreparedMask::Bed {
                    intervals: Bed::load(path)?,
                    mode: *mode,
                },
            });
        }
        Ok(Self { masks: prepared })
    }

    /// A locus-level decision shared by all alleles/score terms at that position.
    /// It is independent of scenario names, pack order, sample order, and threads.
    pub fn hidden(&self, contig: u8, pos1: u32, sample: &str) -> bool {
        self.masks.iter().any(|mask| match mask {
            PreparedMask::Random { rate, seed, unit } => {
                *rate == 1.0 || (*rate > 0.0 && random_draw(*seed, contig, pos1, sample, *unit) < *rate)
            }
            PreparedMask::Bed { intervals, mode } => match mode {
                BedMode::Keep => !intervals.contains(contig, pos1),
                BedMode::Drop => intervals.contains(contig, pos1),
            },
        })
    }
}

/// Per-experiment, one-locus cache. Sample identities are frozen at construction.
/// Sharing a seed/unit reuses the draw across policies and rates, not their decisions.
pub(super) struct RandomDrawCache<'a> {
    samples: &'a [String],
    locus: Option<(u8, u32)>,
    draws: BTreeMap<(u64, bool), Vec<Option<f64>>>,
}
impl<'a> RandomDrawCache<'a> {
    pub(super) fn new(samples: &'a [String]) -> Self {
        Self {
            samples,
            locus: None,
            draws: BTreeMap::new(),
        }
    }
    fn draw(&mut self, seed: u64, unit: RandomUnit, contig: u8, pos1: u32, sample: usize) -> f64 {
        if self.locus != Some((contig, pos1)) {
            for draws in self.draws.values_mut() {
                draws.fill(None);
            }
            self.locus = Some((contig, pos1));
        }
        let per_sample = unit == RandomUnit::SampleSite;
        let index = if per_sample { sample } else { 0 };
        let name = &self.samples[sample];
        let size = if per_sample { self.samples.len() } else { 1 };
        let draws = self.draws.entry((seed, per_sample)).or_insert_with(|| vec![None; size]);
        *draws[index].get_or_insert_with(|| random_draw(seed, contig, pos1, name, unit))
    }
}

impl PreparedMasks {
    pub(super) fn hidden_cached(&self, cache: &mut RandomDrawCache<'_>, contig: u8, pos1: u32, sample: usize) -> bool {
        self.masks.iter().any(|mask| match mask {
            PreparedMask::Random { rate, seed, unit } => {
                *rate == 1.0 || (*rate > 0.0 && cache.draw(*seed, *unit, contig, pos1, sample) < *rate)
            }
            PreparedMask::Bed { intervals, mode } => match mode {
                BedMode::Keep => !intervals.contains(contig, pos1),
                BedMode::Drop => intervals.contains(contig, pos1),
            },
        })
    }
}

fn random_draw(seed: u64, contig: u8, pos1: u32, sample: &str, unit: RandomUnit) -> f64 {
    let mut hasher = Sha256::new();
    hasher.update(MASK_ALGORITHM.as_bytes());
    hasher.update([0]);
    hasher.update(seed.to_le_bytes());
    hasher.update([contig]);
    hasher.update(pos1.to_le_bytes());
    hasher.update([u8::from(unit == RandomUnit::SampleSite)]);
    if unit == RandomUnit::SampleSite {
        hasher.update((sample.len() as u64).to_le_bytes());
        hasher.update(sample.as_bytes());
    }
    let digest = hasher.finalize();
    let bits = u64::from_le_bytes(digest[..8].try_into().expect("eight digest bytes")) >> 11;
    bits as f64 / ((1u64 << 53) as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::{Value, json};
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_DIRECTORY: AtomicU64 = AtomicU64::new(0);

    struct Directory(PathBuf);

    impl Directory {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "pgsum-robustness-config-{}-{}",
                std::process::id(),
                NEXT_DIRECTORY.fetch_add(1, Ordering::Relaxed)
            ));
            std::fs::create_dir_all(&path).unwrap();
            Self(path)
        }

        fn file(&self, name: &str, contents: &str) -> PathBuf {
            let path = self.0.join(name);
            std::fs::write(&path, contents).unwrap();
            path
        }
    }

    impl Drop for Directory {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    fn manifest() -> Value {
        json!({
            "schema": EXPERIMENT_SCHEMA,
            "genotypes": "cohort.pgsc",
            "packs": ["score.pgsp"],
            "scenarios": [{"id": "unchanged"}]
        })
    }

    fn experiment(value: Value) -> Experiment {
        serde_json::from_value(value).unwrap()
    }

    #[test]
    fn defaults_and_scoring_options_match_standard_scoring() {
        let parsed = experiment(manifest());
        parsed.validate().unwrap();
        assert_eq!(parsed.cutoffs, vec![90.0, 95.0]);
        assert_eq!(parsed.top_terms, 10);
        assert_eq!(
            crate::score::Options::from(&parsed.scoring),
            crate::score::Options::default()
        );
        let enabled = ScoringOptions {
            allow_inferred_other_allele: true,
            accept_informational_descriptions: true,
            allow_inferred_palindromes: true,
            allow_inferred_indels: true,
        };
        let converted = crate::score::Options::from(enabled);
        assert!(converted.allow_inferred_other_allele);
        assert!(converted.accept_informational_descriptions);
        assert!(converted.allow_inferred_palindromes);
        assert!(converted.allow_inferred_indels);
    }

    #[test]
    fn paths_resolve_against_manifest_and_absolute_paths_stay_absolute() {
        let directory = Directory::new();
        let absolute_pack = directory.0.join("absolute.pgsp");
        let mut value = manifest();
        value["packs"] = json!(["relative.pgsp", absolute_pack]);
        value["scenarios"] = json!([{
            "id": "frequency",
            "genotypes": "strict.pgsc",
            "fill_frequencies": "freq.pgsf",
            "frequency_source_note": "Frozen independent panel; no evaluation samples.",
            "masks": [{"kind": "bed", "path": "mask.bed", "mode": "drop"}]
        }]);
        value["reference_sets"] = json!([{"id": "reference", "genotypes": "reference.pgsc"}]);
        value["regions"] = json!([{"id": "region", "path": "region.bed"}]);
        let path = directory.file("manifest.json", &value.to_string());
        let parsed = load(&path).unwrap();
        assert_eq!(parsed.genotypes, directory.0.join("cohort.pgsc"));
        assert_eq!(parsed.packs, vec![directory.0.join("relative.pgsp"), absolute_pack]);
        assert_eq!(parsed.scenarios[0].genotypes, Some(directory.0.join("strict.pgsc")));
        assert_eq!(
            parsed.scenarios[0].fill_frequencies,
            Some(directory.0.join("freq.pgsf"))
        );
        let Mask::Bed { path, .. } = &parsed.scenarios[0].masks[0] else {
            panic!("BED mask")
        };
        assert_eq!(path, &directory.0.join("mask.bed"));
        assert_eq!(
            parsed.reference_sets[0].genotypes,
            Some(directory.0.join("reference.pgsc"))
        );
        assert_eq!(parsed.regions[0].path, directory.0.join("region.bed"));
    }

    #[test]
    fn rejects_unknown_fields_at_each_manifest_level() {
        for value in [
            {
                let mut v = manifest();
                v["typo"] = json!(1);
                v
            },
            {
                let mut v = manifest();
                v["scoring"] = json!({"typo": true});
                v
            },
            {
                let mut v = manifest();
                v["scenarios"][0]["typo"] = json!(1);
                v
            },
            {
                let mut v = manifest();
                v["reference_sets"] = json!([{"id": "r", "typo": 1}]);
                v
            },
            {
                let mut v = manifest();
                v["regions"] = json!([{"id": "r", "path": "x", "typo": 1}]);
                v
            },
            {
                let mut v = manifest();
                v["scenarios"][0]["masks"] =
                    json!([{"kind": "random", "rate": 0.1, "seed": 1, "unit": "site", "typo": 1}]);
                v
            },
        ] {
            assert!(serde_json::from_value::<Experiment>(value).is_err());
        }
    }

    #[test]
    fn validates_schema_ids_paths_and_frequency_provenance() {
        let invalid_values = [
            {
                let mut v = manifest();
                v["schema"] = json!("future");
                v
            },
            {
                let mut v = manifest();
                v["genotypes"] = json!("");
                v
            },
            {
                let mut v = manifest();
                v["packs"] = json!([]);
                v
            },
            {
                let mut v = manifest();
                v["packs"] = json!([""]);
                v
            },
            {
                let mut v = manifest();
                v["scenarios"] = json!([]);
                v
            },
            {
                let mut v = manifest();
                v["scenarios"] = json!([{"id": "same"}, {"id": "same"}]);
                v
            },
            {
                let mut v = manifest();
                v["scenarios"][0]["id"] = json!(" bad ");
                v
            },
            {
                let mut v = manifest();
                v["scenarios"][0]["fill_frequencies"] = json!("freq");
                v
            },
            {
                let mut v = manifest();
                v["reference_sets"] = json!([{"id": "r", "samples": []}]);
                v
            },
            {
                let mut v = manifest();
                v["reference_sets"] = json!([{"id": "r", "samples": ["x", "x"]}]);
                v
            },
            {
                let mut v = manifest();
                v["regions"] = json!([{"id": "r", "path": "a"}, {"id": "r", "path": "b"}]);
                v
            },
        ];
        for value in invalid_values {
            assert!(experiment(value.clone()).validate().is_err(), "accepted {value}");
        }
        let mut value = manifest();
        value["reference_sets"] = json!([{"id": "unchanged", "samples": ["x"]}]);
        value["regions"] = json!([{"id": "unchanged", "path": "region.bed"}]);
        experiment(value).validate().unwrap(); // IDs are scoped, not global.
    }

    #[test]
    fn rejects_invalid_numeric_parameters_including_nonfinite_values() {
        let mut parsed = experiment(manifest());
        for top in [0, 1001, usize::MAX] {
            parsed.top_terms = top;
            assert!(parsed.validate().is_err());
        }
        parsed.top_terms = 10;
        for cutoff in [-1.0, 101.0, f64::NAN, f64::INFINITY] {
            parsed.cutoffs = vec![cutoff];
            assert!(parsed.validate().is_err());
        }
        parsed.cutoffs = vec![95.0, 95.0];
        assert!(parsed.validate().is_err());
        parsed.cutoffs = vec![0.0, 100.0];
        parsed.validate().unwrap();
        for rate in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
            assert!(
                PreparedMasks::load(&[Mask::Random {
                    rate,
                    seed: 1,
                    unit: RandomUnit::Site
                }])
                .is_err()
            );
        }
    }

    #[test]
    fn bed_normalizes_merges_and_uses_half_open_anchor_boundaries() {
        let bed =
            Bed::parse("# regions\ntrack name=regions\n1\t10\t20\tone\nchr1\t0\t5\n1\t5\t12\nchrX\t9\t10\nMT\t0\t1\n")
                .unwrap();
        assert_eq!(bed.intervals[&1], vec![(0, 20)]);
        assert!(!bed.contains(1, 0));
        assert!(bed.contains(1, 1));
        assert!(bed.contains(1, 20));
        assert!(!bed.contains(1, 21));
        assert!(!bed.contains(23, 9));
        assert!(bed.contains(23, 10));
        assert!(!bed.contains(23, 11));
        assert!(bed.contains(25, 1));
        assert!(!bed.contains(2, 1));
        assert!(!bed.contains(0, 1));
        let last = Bed::parse(&format!("1 {} {}\n", u32::MAX - 1, u32::MAX)).unwrap();
        assert!(last.contains(1, u32::MAX));
        assert!(!last.contains(1, u32::MAX - 1));
    }

    #[test]
    fn bed_rejects_malformed_lines_instead_of_silently_ignoring_them() {
        for text in [
            "1 0",
            "chrUn 0 10",
            "1 -1 10",
            "1 +1 10",
            "1 1.5 10",
            "1 10 10",
            "1 11 10",
            "1 0 4294967296",
        ] {
            assert!(Bed::parse(text).is_err(), "accepted {text}");
        }
        assert!(!Bed::parse("\n# empty\n").unwrap().contains(1, 1));
    }

    #[test]
    fn bed_keep_drop_and_multiple_masks_have_explicit_union_semantics() {
        let directory = Directory::new();
        let path = directory.file("region.bed", "1\t9\t20\n");
        let keep = PreparedMasks::load(&[Mask::Bed {
            path: path.clone(),
            mode: BedMode::Keep,
        }])
        .unwrap();
        let drop = PreparedMasks::load(&[Mask::Bed {
            path: path.clone(),
            mode: BedMode::Drop,
        }])
        .unwrap();
        let union = PreparedMasks::load(&[
            Mask::Bed {
                path: path.clone(),
                mode: BedMode::Keep,
            },
            Mask::Bed {
                path,
                mode: BedMode::Drop,
            },
        ])
        .unwrap();
        for pos in 1..=30 {
            assert_eq!(keep.hidden(1, pos, "sample"), !(10..=20).contains(&pos));
            assert_eq!(drop.hidden(1, pos, "sample"), (10..=20).contains(&pos));
            assert!(union.hidden(1, pos, "sample"));
        }
        assert!(!PreparedMasks::load(&[]).unwrap().hidden(1, 1, "sample"));
    }

    #[test]
    fn seeded_masks_are_order_independent_nested_and_respect_units() {
        let prepare = |rate, unit| {
            PreparedMasks::load(&[Mask::Random {
                rate,
                seed: 12345,
                unit,
            }])
            .unwrap()
        };
        let low = prepare(0.2, RandomUnit::SampleSite);
        let high = prepare(0.8, RandomUnit::SampleSite);
        let sites = prepare(0.5, RandomUnit::Site);
        let none = prepare(0.0, RandomUnit::SampleSite);
        let all = prepare(1.0, RandomUnit::SampleSite);
        let mut samples_differ = false;
        let forward: Vec<_> = (1..=100).map(|pos| low.hidden(1, pos, "a")).collect();
        let mut reverse: Vec<_> = (1..=100).rev().map(|pos| low.hidden(1, pos, "a")).collect();
        reverse.reverse();
        assert_eq!(forward, reverse);
        for pos in 1..=100 {
            assert!(!low.hidden(1, pos, "a") || high.hidden(1, pos, "a"));
            assert_eq!(sites.hidden(1, pos, "a"), sites.hidden(1, pos, "b"));
            samples_differ |= low.hidden(1, pos, "a") != low.hidden(1, pos, "b");
            assert!(!none.hidden(1, pos, "a"));
            assert!(all.hidden(1, pos, "a"));
        }
        assert!(samples_differ);
    }

    #[test]
    fn cached_masks_match_uncached_masks_across_rates_units_samples_and_loci() {
        let samples = vec!["a".into(), "NA12878".into(), "sample α".into()];
        let mut masks = Vec::new();
        for seed in [0, 12345, u64::MAX] {
            for unit in [RandomUnit::Site, RandomUnit::SampleSite] {
                for rate in [0.0, 0.2, 0.5, 0.8, 1.0] {
                    masks.push(PreparedMasks::load(&[Mask::Random { rate, seed, unit }]).unwrap());
                }
            }
        }
        for mode in [BedMode::Keep, BedMode::Drop] {
            masks.push(PreparedMasks {
                masks: vec![
                    PreparedMask::Bed {
                        intervals: Bed::parse("1 1 3\n22 0 10\n").unwrap(),
                        mode,
                    },
                    PreparedMask::Random {
                        rate: 0.4,
                        seed: 12345,
                        unit: RandomUnit::SampleSite,
                    },
                    PreparedMask::Random {
                        rate: 0.2,
                        seed: 0,
                        unit: RandomUnit::Site,
                    },
                ],
            });
        }
        let mut cache = RandomDrawCache::new(&samples);
        for (contig, pos) in [(1, 1000), (1, 2), (1, 2), (2, 2), (1, 1000), (22, u32::MAX), (1, 1)] {
            for mask in masks.iter().rev() {
                for sample in [2, 0, 1] {
                    assert_eq!(
                        mask.hidden_cached(&mut cache, contig, pos, sample),
                        mask.hidden(contig, pos, &samples[sample])
                    );
                }
            }
        }
        // Verify the cached values against independent SHA-256 vectors, including UTF-8 names.
        assert_eq!(
            cache.draw(12345, RandomUnit::Site, 1, 1000, 0),
            3_076_170_996_652_624_f64 / ((1u64 << 53) as f64)
        );
        assert_eq!(
            cache.draw(12345, RandomUnit::SampleSite, 1, 1000, 1),
            3_112_782_897_755_573_f64 / ((1u64 << 53) as f64)
        );
        assert_eq!(
            cache.draw(12345, RandomUnit::SampleSite, 1, 1000, 2),
            2_833_699_031_287_665_f64 / ((1u64 << 53) as f64)
        );
        // Only one slot per seed/unit exists, not one slot per rate/policy or genomic term.
        assert_eq!(cache.draws.len(), 6);
    }

    #[test]
    fn masking_hash_matches_independently_computed_versioned_vectors() {
        // Computed independently with Python hashlib and struct.pack("<QBI").
        // Pin the bit-level algorithm so refactoring cannot silently change runs.
        for (unit, sample, numerator) in [
            (RandomUnit::Site, "ignored sample", 3_076_170_996_652_624u64),
            (RandomUnit::SampleSite, "NA12878", 3_112_782_897_755_573),
            (RandomUnit::SampleSite, "sample α", 2_833_699_031_287_665),
        ] {
            assert_eq!(
                random_draw(12345, 1, 1000, sample, unit),
                numerator as f64 / ((1u64 << 53) as f64)
            );
        }
    }
}
