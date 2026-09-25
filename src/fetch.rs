//! Download PGS Catalog scores and compile each into a pack as it arrives.
//!
//! For each score: fetch its Catalog REST record, skip it if a pack compiled under the current format and
//! rules already exists with an identical record, otherwise download the harmonized GRCh38 scoring file,
//! compile it, and delete the download (unless kept).
//! Only compiled packs need to stay on disk, which is what makes the whole Catalog fit.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::compile::compile_file;
use crate::pack::{self, Pack, ReferenceIdentity};
use crate::reference::Reference;
use crate::scoring_file::is_pgs_id;
use crate::{Error, Result, invalid};

pub const REST: &str = "https://www.pgscatalog.org/rest";
const ATTEMPTS: u32 = 8;
/// Above this many IDs, records come from the paginated listing (a few dozen requests) instead of one request
/// per score, which the Catalog's API rate-limits.
const LISTING_ABOVE: usize = 25;

pub fn agent() -> ureq::Agent {
    ureq::Agent::config_builder()
        .user_agent(concat!("pgsum/", env!("CARGO_PKG_VERSION")))
        .timeout_connect(Some(Duration::from_secs(30)))
        .timeout_recv_body(Some(Duration::from_secs(600)))
        .build()
        .new_agent()
}

/// Retry transient failures (network errors, 429 and 5xx) with exponential backoff; fail at once otherwise.
fn with_retries<T>(what: &str, mut f: impl FnMut() -> std::result::Result<T, ureq::Error>) -> Result<T> {
    let mut last = String::new();
    for attempt in 0..ATTEMPTS {
        match f() {
            Ok(v) => return Ok(v),
            Err(ureq::Error::StatusCode(code)) if code != 429 && code < 500 => {
                return invalid!("{what}: HTTP {code}");
            }
            Err(e) => last = e.to_string(),
        }
        if attempt + 1 < ATTEMPTS {
            // 1 s, 2 s, 4 s … 64 s: about two minutes in all, enough for a rate limit to clear.
            std::thread::sleep(Duration::from_secs(2u64.pow(attempt)));
        }
    }
    invalid!("{what}: {last} (after {ATTEMPTS} attempts)")
}

fn get_json(agent: &ureq::Agent, url: &str) -> Result<serde_json::Value> {
    let text = with_retries(url, || {
        agent
            .get(url)
            .call()?
            .body_mut()
            .with_config()
            .limit(256 << 20)
            .read_to_string()
    })?;
    serde_json::from_str(&text).map_err(|e| Error::Invalid(format!("{url}: {e}")))
}

/// Catalog records for the given IDs, or for every score when `ids` is empty.
pub fn catalog_records(agent: &ureq::Agent, ids: &[String]) -> Result<Vec<serde_json::Value>> {
    if let Some(bad) = ids.iter().find(|id| !is_pgs_id(id)) {
        return invalid!("{bad:?} is not a PGS ID");
    }
    if ids.len() > LISTING_ABOVE {
        let wanted: std::collections::HashSet<&str> = ids.iter().map(String::as_str).collect();
        let records: Vec<_> = all_records(agent)?
            .into_iter()
            .filter(|r| r["id"].as_str().is_some_and(|id| wanted.contains(id)))
            .collect();
        let found: std::collections::HashSet<&str> = records.iter().filter_map(|r| r["id"].as_str()).collect();
        let missing: Vec<_> = ids.iter().filter(|id| !found.contains(id.as_str())).collect();
        if !missing.is_empty() {
            return invalid!("not in the Catalog listing: {missing:?}");
        }
        return Ok(records);
    }
    if !ids.is_empty() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .map_err(|e| Error::Invalid(e.to_string()))?;
        return pool.install(|| {
            ids.par_iter()
                .map(|id| get_json(agent, &format!("{REST}/score/{id}")))
                .collect()
        });
    }
    all_records(agent)
}

/// Every score's record, from the paginated listing.
fn all_records(agent: &ureq::Agent) -> Result<Vec<serde_json::Value>> {
    let mut records = Vec::new();
    let mut url = Some(format!("{REST}/score/all?limit=250"));
    while let Some(u) = url {
        let page = get_json(agent, &u)?;
        records.extend(page["results"].as_array().cloned().unwrap_or_default());
        url = page["next"].as_str().map(str::to_owned);
    }
    Ok(records)
}

pub fn harmonized_url(record: &serde_json::Value) -> Option<&str> {
    record["ftp_harmonized_scoring_files"]["GRCh38"]["positions"].as_str()
}

#[derive(Clone, Debug, PartialEq)]
pub enum Outcome {
    Compiled { terms: u64, source_bytes: u64 },
    Unchanged,
    Skipped(String),
    Failed(String),
}

pub struct Options<'a> {
    pub reference: &'a Reference,
    pub identity: &'a ReferenceIdentity,
    pub out: &'a Path,
    pub downloads: &'a Path,
    pub keep_downloads: bool,
    pub max_variants: Option<u64>,
}

/// Fetch and compile one score.
pub fn fetch_one(agent: &ureq::Agent, record: &serde_json::Value, o: &Options<'_>) -> Outcome {
    let Some(id) = record["id"].as_str().filter(|id| is_pgs_id(id)) else {
        return Outcome::Failed("record has no valid id".into());
    };
    let pack_path = o.out.join(format!("{id}.{}", pack::EXTENSION));
    if let Ok(existing) = Pack::open_header(&pack_path)
        && existing.schema == pack::SCHEMA
        && existing.compile_rules == pack::COMPILE_RULES
        && &existing.catalog_metadata == record
        && &existing.reference == o.identity
    {
        return Outcome::Unchanged;
    }
    if let (Some(max), Some(n)) = (o.max_variants, record["variants_number"].as_u64())
        && n > max
    {
        return Outcome::Skipped(format!("{n} variants exceeds --max-variants {max}"));
    }
    let Some(url) = harmonized_url(record) else {
        return Outcome::Skipped("no harmonized GRCh38 scoring file".into());
    };
    let scoring = o.downloads.join(format!("{id}_hmPOS_GRCh38.txt.gz"));
    let metadata = o.downloads.join(format!("{id}.metadata.json"));
    let result = (|| -> Result<Outcome> {
        let json = serde_json::to_vec_pretty(record).map_err(|e| Error::Invalid(e.to_string()))?;
        std::fs::write(&metadata, json).map_err(Error::io(&metadata))?;
        download(agent, url, &scoring)?;
        let header = compile_file(&scoring, Some(&metadata), o.reference, o.identity, o.out)?;
        Ok(Outcome::Compiled {
            terms: header.inventory.actual_terms,
            source_bytes: header.source.bytes,
        })
    })();
    if !o.keep_downloads {
        let _ = std::fs::remove_file(&scoring);
        let _ = std::fs::remove_file(&metadata);
    }
    result.unwrap_or_else(|e| Outcome::Failed(e.to_string()))
}

/// Download to `<path>.part`, then rename into place.
fn download(agent: &ureq::Agent, url: &str, path: &Path) -> Result<()> {
    let part = PathBuf::from(format!("{}.part", path.display()));
    with_retries(url, || {
        let response = agent.get(url).call()?;
        let mut reader = response.into_body().into_reader();
        let mut file = std::io::BufWriter::new(std::fs::File::create(&part)?);
        std::io::copy(&mut reader, &mut file)?;
        file.flush()?;
        Ok(())
    })
    .inspect_err(|_| {
        let _ = std::fs::remove_file(&part);
    })?;
    std::fs::rename(&part, path).map_err(Error::io(path))
}

/// One line of the fetch log.
pub struct Logged {
    pub id: String,
    pub outcome: Outcome,
    pub seconds: f64,
}

/// Fetch `records` with `jobs` scores in flight at once, reporting each as it finishes.
pub fn fetch_all(
    agent: &ureq::Agent,
    records: &[serde_json::Value],
    o: &Options<'_>,
    jobs: usize,
    report: &(dyn Fn(&Logged) + Sync),
) -> Result<Vec<Logged>> {
    std::fs::create_dir_all(o.out).map_err(Error::io(o.out))?;
    std::fs::create_dir_all(o.downloads).map_err(Error::io(o.downloads))?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(jobs.max(1))
        .build()
        .map_err(|e| Error::Invalid(e.to_string()))?;
    Ok(pool.install(|| {
        records
            .par_iter()
            .map(|record| {
                let started = Instant::now();
                let logged = Logged {
                    id: record["id"].as_str().unwrap_or("?").to_owned(),
                    outcome: fetch_one(agent, record, o),
                    seconds: started.elapsed().as_secs_f64(),
                };
                report(&logged);
                logged
            })
            .collect()
    }))
}
