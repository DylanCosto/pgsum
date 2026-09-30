//! Verified, immutable Catalog packs for `run --ids` and `batch --ids`.
use crate::pack::{Pack, ReferenceIdentity};
use crate::reference::Reference;
use crate::{Error, Result, invalid};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

pub struct Options {
    pub directory: Option<PathBuf>,
    pub offline: bool,
    pub refresh: bool,
    pub jobs: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            directory: None,
            offline: false,
            refresh: false,
            jobs: 4,
        }
    }
}

#[derive(Serialize, Deserialize)]
struct Receipt {
    schema: String,
    id: String,
    pack_sha256: String,
    fetched_unix_seconds: u64,
}

pub fn default_directory() -> Result<PathBuf> {
    let env = |name| std::env::var_os(name).filter(|v| !v.is_empty()).map(PathBuf::from);
    env("PGSUM_CACHE_DIR")
        .or_else(|| env("XDG_CACHE_HOME").map(|p| p.join("pgsum")))
        .or_else(|| env("HOME").map(|p| p.join(".cache/pgsum")))
        .ok_or_else(|| Error::Invalid("no cache directory available; supply --catalog-cache".into()))
}

pub fn resolve(ids: &[String], reference: &Path, options: &Options) -> Result<Vec<PathBuf>> {
    let agent = crate::fetch::agent();
    resolve_with(ids, reference, options, &agent, |ids| {
        crate::fetch::catalog_records(&agent, ids)
    })
}

/// Injectable record provider keeps the actual download/compile/cache path independently testable.
pub fn resolve_with(
    ids: &[String],
    reference: &Path,
    options: &Options,
    agent: &ureq::Agent,
    records: impl FnOnce(&[String]) -> Result<Vec<serde_json::Value>>,
) -> Result<Vec<PathBuf>> {
    if ids.is_empty() || ids.iter().any(|id| !crate::scoring_file::is_pgs_id(id)) {
        return invalid!("automatic downloads require PGS Catalog IDs; use --pack for custom scores");
    }
    if options.offline && options.refresh {
        return invalid!("--offline cannot be combined with --refresh-catalog");
    }
    if options.jobs == 0 {
        return invalid!("Catalog download jobs must be positive");
    }
    let mut ids = ids.to_vec();
    ids.sort();
    ids.dedup();
    let reference = Reference::open(reference)?;
    let fresh = crate::digest::file_sha256(&reference.path)?;
    crate::digest::remember_file_sha256(&reference.path, &fresh);
    let mut identity = crate::compile::reference_identity(&reference)?;
    identity.fasta_sha256 = fresh;
    let key = crate::pack::records_sha256(
        &serde_json::to_vec(&serde_json::json!({
            "reference": identity, "schema": crate::pack::SCHEMA, "rules": crate::pack::COMPILE_RULES,
            "build": "GRCh38"
        }))
        .map_err(|e| Error::Invalid(e.to_string()))?,
    );
    let base = options.directory.clone().map(Ok).unwrap_or_else(default_directory)?;
    let root = base.join("catalog/GRCh38").join(&key[7..]);
    std::fs::create_dir_all(&root).map_err(Error::io(&root))?;
    let lock_path = root.join(".lock");
    let lock = std::fs::OpenOptions::new()
        .create(true)
        .truncate(false)
        .read(true)
        .write(true)
        .open(&lock_path)
        .map_err(Error::io(&lock_path))?;
    fs2::FileExt::try_lock_exclusive(&lock).map_err(|e| {
        Error::Invalid(format!(
            "Catalog cache is in use ({}): {e}; retry after the other run finishes",
            root.display()
        ))
    })?;
    let mut paths = std::collections::BTreeMap::new();
    let mut needed = Vec::new();
    for id in &ids {
        if !options.refresh
            && let Some(path) = cached(&root.join(id), id, &identity)
        {
            paths.insert(id.clone(), path);
        } else {
            needed.push(id.clone());
        }
    }
    if options.offline && !needed.is_empty() {
        return invalid!(
            "offline cache has no valid packs for {}; run once without --offline",
            needed.join(", ")
        );
    }
    if !needed.is_empty() {
        let records = records(&needed)?;
        let returned: std::collections::BTreeSet<_> = records.iter().filter_map(|r| r["id"].as_str()).collect();
        let expected: std::collections::BTreeSet<_> = needed.iter().map(String::as_str).collect();
        if records.len() != needed.len() || returned != expected {
            return invalid!("Catalog response IDs do not match requested scores");
        }
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_err(|e| Error::Invalid(e.to_string()))?;
        let stage = root.join(format!(".pending-{}-{}", std::process::id(), nonce.as_nanos()));
        std::fs::create_dir(&stage).map_err(Error::io(&stage))?;
        let downloads = stage.join("downloads");
        let fetch_options = crate::fetch::Options {
            reference: &reference,
            identity: &identity,
            out: &stage,
            downloads: &downloads,
            keep_downloads: false,
            max_variants: None,
            public: None,
        };
        let log = crate::fetch::fetch_all(agent, &records, &fetch_options, options.jobs, &|l| {
            eprintln!("Catalog {}: {:?}", l.id, l.outcome);
        })?;
        let publish = (|| -> Result<()> {
            for item in log {
                match item.outcome {
                    crate::fetch::Outcome::Compiled { .. } => {}
                    other => return invalid!("Catalog {}: {other:?}", item.id),
                }
                let name = format!("{}.pgsp", item.id);
                let source = stage.join(&name);
                let pack = Pack::open(&source)?;
                if pack.header.pgs_id != item.id || !pack.header.inventory.consistent {
                    return invalid!("Catalog {} compiled with an inconsistent term inventory", item.id);
                }
                drop(pack);
                let hash = crate::digest::file_sha256(&source)?;
                let parent = root.join(&item.id);
                let generation = parent.join(&hash[7..]);
                std::fs::create_dir_all(&generation).map_err(Error::io(&generation))?;
                let destination = generation.join(&name);
                // Same-content refreshes can repair a damaged cached file atomically.
                std::fs::rename(&source, &destination).map_err(Error::io(&destination))?;
                let receipt = Receipt {
                    schema: "pgsum-catalog-cache-v1".into(),
                    id: item.id.clone(),
                    pack_sha256: hash,
                    fetched_unix_seconds: nonce.as_secs(),
                };
                let tmp = parent.join("current.json.tmp");
                std::fs::write(
                    &tmp,
                    serde_json::to_vec_pretty(&receipt).map_err(|e| Error::Invalid(e.to_string()))?,
                )
                .map_err(Error::io(&tmp))?;
                std::fs::rename(&tmp, parent.join("current.json")).map_err(Error::io(&tmp))?;
                paths.insert(item.id, destination);
            }
            Ok(())
        })();
        let _ = std::fs::remove_dir_all(&stage);
        publish?;
    }
    Ok(ids.iter().map(|id| paths[id].clone()).collect())
}

fn cached(root: &Path, id: &str, reference: &ReferenceIdentity) -> Option<PathBuf> {
    let receipt: Receipt = serde_json::from_slice(&std::fs::read(root.join("current.json")).ok()?).ok()?;
    let hash = receipt.pack_sha256.strip_prefix("sha256:")?;
    if receipt.schema != "pgsum-catalog-cache-v1"
        || receipt.id != id
        || hash.len() != 64
        || !hash.bytes().all(|b| b.is_ascii_hexdigit())
    {
        return None;
    }
    let path = root.join(hash).join(format!("{id}.pgsp"));
    if crate::digest::file_sha256(&path).ok()? != receipt.pack_sha256 {
        return None;
    }
    let header = Pack::open_header(&path).ok()?;
    (header.pgs_id == id
        && header.reference == *reference
        && header.schema == crate::pack::SCHEMA
        && header.compile_rules == crate::pack::COMPILE_RULES
        && header.public_variants.is_none()
        && header.inventory.consistent)
        .then_some(path)
}
