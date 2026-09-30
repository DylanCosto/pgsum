//! Resumable sample-sheet execution with content-addressed, independently verifiable sample runs.
use crate::{Error, Result, invalid};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashSet};
use std::fs::{File, OpenOptions};
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Sample {
    pub sample_id: String,
    pub input: PathBuf,
    pub sample: Option<String>,
}

pub fn read_sheet(path: &Path) -> Result<Vec<Sample>> {
    let text = std::fs::read_to_string(path).map_err(Error::io(path))?;
    let mut lines = text
        .lines()
        .enumerate()
        .filter(|(_, s)| !s.trim().is_empty() && !s.starts_with('#'));
    let (_, header) = lines
        .next()
        .ok_or_else(|| Error::Invalid("empty sample sheet".into()))?;
    let columns: Vec<_> = header.split('\t').collect();
    let allowed = ["sample_id", "input", "sample"];
    let unique: HashSet<_> = columns.iter().collect();
    if unique.len() != columns.len() || columns.iter().any(|c| !allowed.contains(c)) {
        return invalid!("sample sheet columns must be sample_id, input and optional sample");
    }
    let col = |name| {
        columns
            .iter()
            .position(|c| *c == name)
            .ok_or_else(|| Error::Invalid(format!("sample sheet needs {name}")))
    };
    let (id, input) = (col("sample_id")?, col("input")?);
    let selected = columns.iter().position(|c| *c == "sample");
    let base = path.parent().unwrap_or(Path::new("."));
    let mut seen = HashSet::new();
    let mut samples = Vec::new();
    for (line, text) in lines {
        let fields: Vec<_> = text.split('\t').collect();
        if fields.len() != columns.len() {
            return invalid!(
                "sample sheet line {} has {} columns, expected {}",
                line + 1,
                fields.len(),
                columns.len()
            );
        }
        let name = fields[id];
        if name.is_empty()
            || name.len() > 128
            || name.starts_with('.')
            || !name.bytes().all(|b| b.is_ascii_alphanumeric() || b"._-".contains(&b))
            || !seen.insert(name.to_owned())
        {
            return invalid!(
                "sample sheet line {} has an invalid or duplicate sample_id {name:?}",
                line + 1
            );
        }
        if fields[input].is_empty() {
            return invalid!("sample sheet line {} has no input", line + 1);
        }
        let input = base.join(fields[input]);
        let input = std::fs::canonicalize(&input).map_err(Error::io(&input))?;
        if !input.is_file() {
            return invalid!("{} is not an input file", input.display());
        }
        samples.push(Sample {
            sample_id: name.into(),
            input,
            sample: selected.map(|i| fields[i]).filter(|v| !v.is_empty()).map(str::to_owned),
        });
    }
    if samples.is_empty() {
        return invalid!("sample sheet contains no samples");
    }
    Ok(samples)
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Settings {
    pub dosage_field: crate::dosage::Field,
    pub accept_missing_quality: bool,
    pub haploid_xy_as_homozygous: bool,
    pub skip_structural_alleles: bool,
    pub merge_split_records: bool,
    pub allow_inferred_other_allele: bool,
    pub accept_informational_descriptions: bool,
    pub allow_inferred_palindromes: bool,
    pub allow_inferred_indels: bool,
    pub terms: bool,
    pub bundle: bool,
}

pub struct Options<'a> {
    pub sheet: &'a Path,
    pub reference: &'a Path,
    pub packs: &'a [PathBuf],
    pub out: &'a Path,
    pub settings: Settings,
    pub jobs: usize,
    pub threads: usize,
    pub resume: bool,
}

#[derive(Serialize, Deserialize)]
struct Manifest {
    schema: String,
    signature: String,
    request: serde_json::Value,
    outputs: BTreeMap<String, String>,
}

#[derive(Serialize)]
pub struct SampleResult {
    pub sample_id: String,
    pub status: String,
    pub directory: PathBuf,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

fn json(path: &Path, value: &impl Serialize) -> Result<()> {
    let data = serde_json::to_vec_pretty(value).map_err(|e| Error::Invalid(e.to_string()))?;
    let temporary = path.with_extension("json.tmp");
    std::fs::write(&temporary, data).map_err(Error::io(&temporary))?;
    File::open(&temporary)
        .and_then(|f| f.sync_all())
        .map_err(Error::io(&temporary))?;
    std::fs::rename(&temporary, path).map_err(Error::io(path))
}

fn signature(value: &impl Serialize) -> Result<String> {
    Ok(crate::pack::records_sha256(
        &serde_json::to_vec(value).map_err(|e| Error::Invalid(e.to_string()))?,
    ))
}

fn reusable(directory: &Path, expected: &str) -> bool {
    let Ok(bytes) = std::fs::read(directory.join("run-manifest.json")) else {
        return false;
    };
    let Ok(manifest) = serde_json::from_slice::<Manifest>(&bytes) else {
        return false;
    };
    manifest.schema == "pgsum-sample-run-v1"
        && manifest.signature == expected
        && signature(&manifest.request).ok().as_deref() == Some(expected)
        && !manifest.outputs.is_empty()
        && manifest.outputs.iter().all(|(name, digest)| {
            let p = Path::new(name);
            p.components().count() == 1
                && matches!(p.components().next(), Some(std::path::Component::Normal(_)))
                && crate::digest::file_sha256(&directory.join(p)).ok().as_ref() == Some(digest)
        })
}

/// Combine successful sample summaries without retaining the score matrix in memory.
fn aggregate(root: &Path, results: &[SampleResult]) -> Result<String> {
    let path = root.join("batch-scores.tsv");
    let temporary = root.join("batch-scores.tsv.tmp");
    let mut output = std::io::BufWriter::new(File::create(&temporary).map_err(Error::io(&temporary))?);
    let mut header = None;
    for sample in results.iter().filter(|s| s.error.is_none()) {
        let input = sample.directory.join("scores.tsv");
        let mut lines = std::io::BufReader::new(File::open(&input).map_err(Error::io(&input))?).lines();
        let first = lines
            .next()
            .transpose()
            .map_err(Error::io(&input))?
            .ok_or_else(|| Error::Invalid(format!("{} is empty", input.display())))?;
        if let Some(expected) = &header {
            if expected != &first {
                return invalid!("sample {} has an incompatible score summary", sample.sample_id);
            }
        } else {
            writeln!(output, "sample_id\t{first}").map_err(Error::io(&temporary))?;
            header = Some(first);
        }
        for line in lines {
            writeln!(output, "{}\t{}", sample.sample_id, line.map_err(Error::io(&input))?)
                .map_err(Error::io(&temporary))?;
        }
    }
    if header.is_none() {
        writeln!(output, "sample_id\tpgs_id\tstatus").map_err(Error::io(&temporary))?;
    }
    output.flush().map_err(Error::io(&temporary))?;
    output.get_ref().sync_all().map_err(Error::io(&temporary))?;
    drop(output);
    std::fs::rename(&temporary, &path).map_err(Error::io(&path))?;
    crate::digest::file_sha256(&path)
}

pub fn run(options: &Options<'_>) -> Result<Vec<SampleResult>> {
    if options.jobs == 0 || options.threads == 0 || options.jobs > options.threads {
        return invalid!("batch --jobs must be between 1 and the total --threads budget");
    }
    let samples = read_sheet(options.sheet)?;
    let reference_path = std::fs::canonicalize(options.reference).map_err(Error::io(options.reference))?;
    let reference = crate::reference::Reference::open(&reference_path)?;
    let mut identity = crate::compile::reference_identity(&reference)?;
    // Resume verification must inspect content, even if size and timestamps were preserved.
    identity.fasta_sha256 = crate::digest::file_sha256(&reference_path)?;
    crate::digest::remember_file_sha256(&reference_path, &identity.fasta_sha256);
    let packs = options
        .packs
        .iter()
        .map(|p| std::fs::canonicalize(p).map_err(Error::io(p)))
        .collect::<Result<Vec<_>>>()?;
    let pack_hashes = packs
        .iter()
        .map(|p| Ok((p.clone(), crate::digest::file_sha256(p)?)))
        .collect::<Result<BTreeMap<_, _>>>()?;
    let executable = std::env::current_exe().map_err(Error::io("<executable>"))?;
    let executable_hash = crate::digest::file_sha256(&executable)?;
    let source_hashes = samples
        .iter()
        .map(|s| &s.input)
        .collect::<HashSet<_>>()
        .into_iter()
        .map(|p| Ok((p.clone(), crate::digest::file_sha256(p)?)))
        .collect::<Result<BTreeMap<_, _>>>()?;
    std::fs::create_dir_all(options.out).map_err(Error::io(options.out))?;
    // OS-backed lock: automatically released after a crash, without interpreting stale PID files.
    let lock_path = options.out.join(".batch.lock");
    let lock = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(&lock_path)
        .map_err(Error::io(&lock_path))?;
    fs2::FileExt::try_lock_exclusive(&lock)
        .map_err(|e| Error::Invalid(format!("another batch owns {}: {e}", options.out.display())))?;
    let root = std::fs::canonicalize(options.out).map_err(Error::io(options.out))?;
    let cache = root.join("targets.pgst");
    // Build once under the lock before any worker opens the shared target index.
    let expected_packs = {
        let (targets, _) = crate::targets::targets(&packs, &identity, Some(&cache), false)?;
        targets.packs
    };
    drop(reference);
    for (path, digest) in &source_hashes {
        crate::digest::remember_file_sha256(path, digest);
    }
    let settings = &options.settings;
    let per_sample_threads = options.threads / options.jobs;
    let context = serde_json::json!({ "pgsum_version": env!("CARGO_PKG_VERSION"), "executable_sha256": executable_hash,
        "reference": identity, "pack_sha256": pack_hashes, "settings": settings, "threads": per_sample_threads });
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.jobs)
        .build()
        .map_err(|e| Error::Invalid(e.to_string()))?;
    let results = pool.install(|| samples.par_iter().map(|sample| {
        let request = serde_json::json!({ "context": context, "sample": sample, "input_sha256": source_hashes[&sample.input] });
        let key = signature(&request).expect("JSON value serializes");
        let parent = root.join("samples").join(&sample.sample_id);
        let destination = parent.join(key.trim_start_matches("sha256:"));
        let execute = || -> Result<bool> {
            if destination.exists() {
                if !options.resume { return invalid!("{} already exists; use --resume to verify and reuse it", destination.display()); }
                if reusable(&destination, &key) { return Ok(true); }
            }
            std::fs::create_dir_all(&parent).map_err(Error::io(&parent))?;
            let nonce = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_err(|e| Error::Invalid(e.to_string()))?.as_nanos();
            let stage = parent.join(format!(".pending-{}-{nonce}", std::process::id()));
            std::fs::create_dir(&stage).map_err(Error::io(&stage))?;
            let log_path = stage.join("run.log");
            let log = File::create(&log_path).map_err(Error::io(&log_path))?;
            let mut cmd = Command::new(&executable);
            cmd.arg("run").arg("--gvcf").arg(&sample.input).arg("--reference").arg(&reference_path)
                .arg("--out").arg(&stage).arg("--targets-cache").arg(&cache)
                .arg("--threads").arg(per_sample_threads.to_string())
                .arg("--dosage-field").arg(settings.dosage_field.label().to_lowercase());
            for pack in &packs { cmd.arg("--pack").arg(pack); }
            if let Some(name) = &sample.sample { cmd.arg("--sample").arg(name); }
            for (enabled, flag) in [
                (settings.accept_missing_quality, "--accept-missing-quality"),
                (settings.haploid_xy_as_homozygous, "--haploid-xy-as-homozygous"),
                (settings.skip_structural_alleles, "--skip-structural-alleles"),
                (settings.merge_split_records, "--merge-split-records"),
                (settings.allow_inferred_other_allele, "--allow-inferred-other-allele"),
                (settings.accept_informational_descriptions, "--accept-informational-descriptions"),
                (settings.allow_inferred_palindromes, "--allow-inferred-palindromes"),
                (settings.allow_inferred_indels, "--allow-inferred-indels"),
                (settings.terms, "--terms"), (settings.bundle, "--bundle"),
            ] { if enabled { cmd.arg(flag); } }
            let status = cmd.stdin(Stdio::null()).stdout(log.try_clone().map_err(Error::io(&log_path))?).stderr(log).status().map_err(Error::io(&executable))?;
            if !status.success() { return invalid!("sample {} failed ({status}); inspect {}", sample.sample_id, log_path.display()); }
            let table = crate::genotypes::GenotypeTable::open(&stage.join("genotypes.pgsg"))?;
            if table.header.gvcf.sha256 != source_hashes[&sample.input] || table.header.reference != identity || table.header.packs != expected_packs {
                return invalid!("sample {} inputs changed during execution; incomplete result retained in {}", sample.sample_id, stage.display());
            }
            let mut outputs = BTreeMap::new();
            for entry in std::fs::read_dir(&stage).map_err(Error::io(&stage))? {
                let path = entry.map_err(Error::io(&stage))?.path();
                if path.is_file() {
                    outputs.insert(path.file_name().unwrap().to_string_lossy().into_owned(), crate::digest::file_sha256(&path)?);
                }
            }
            let manifest = Manifest { schema: "pgsum-sample-run-v1".into(), signature: key.clone(), request, outputs };
            json(&stage.join("run-manifest.json"), &manifest)?;
            // Keep a damaged previous generation for diagnosis instead of destroying it.
            if destination.exists() { std::fs::rename(&destination, parent.join(format!(".invalid-{nonce}"))).map_err(Error::io(&destination))?; }
            std::fs::rename(&stage, &destination).map_err(Error::io(&destination))?;
            Ok(false)
        };
        match execute() {
            Ok(reused) => SampleResult { sample_id: sample.sample_id.clone(), status: if reused { "reused" } else { "completed" }.into(), directory: destination, error: None },
            Err(e) => SampleResult { sample_id: sample.sample_id.clone(), status: "failed".into(), directory: destination, error: Some(e.to_string()) },
        }
    }).collect::<Vec<_>>());
    let aggregate_sha256 = aggregate(&root, &results)?;
    let report = serde_json::json!({ "schema": "pgsum-batch-v1", "samplesheet_sha256": crate::digest::file_sha256(options.sheet)?, "jobs": options.jobs, "threads_per_sample": per_sample_threads,
        "scores": { "path": "batch-scores.tsv", "sha256": aggregate_sha256 },
        "complete": results.iter().all(|r| r.error.is_none()), "samples": results });
    json(&root.join("batch-results.json"), &report)?;
    let failures = results.iter().filter(|r| r.status == "failed").count();
    if failures > 0 {
        return invalid!(
            "{failures} samples failed; see {}",
            root.join("batch-results.json").display()
        );
    }
    Ok(results)
}
