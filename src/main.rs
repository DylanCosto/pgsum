use clap::{Args, Parser, Subcommand};
use pgsum::compile::{compile_file, reference_identity};
use pgsum::genotypes::GenotypeTable;
use pgsum::pack::Pack;
use pgsum::reference::Reference;
use pgsum::{Error, Result};
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

/// Polygenic score calculation from single-sample gVCFs.
#[derive(Parser)]
#[command(version, about)]
struct Cli {
    /// Worker threads for every step (default: all cores).
    #[arg(long, global = true)]
    threads: Option<usize>,
    #[command(subcommand)]
    command: Command,
}

/// Which packs to use. `--pack` takes files or directories (every `*.pgsp` inside) and can be repeated;
/// `--pack-list` adds paths from a file; `--ids` keeps only the listed scores.
#[derive(Args, Clone)]
struct PackSelection {
    /// A pack, or a directory of packs. Repeatable.
    #[arg(long = "pack")]
    packs: Vec<PathBuf>,
    /// A file listing pack paths or directories, one per line (`#` comments allowed).
    #[arg(long)]
    pack_list: Option<PathBuf>,
    /// Keep only these score IDs (comma-separated or repeated), e.g. `PGS000001,PGS000013,MY_SCORE`.
    #[arg(long, value_delimiter = ',')]
    ids: Vec<String>,
}

impl PackSelection {
    fn resolve(&self) -> Result<Vec<PathBuf>> {
        let mut roots = self.packs.clone();
        if let Some(list) = &self.pack_list {
            let text = std::fs::read_to_string(list).map_err(Error::io(list))?;
            let base = list.parent().unwrap_or(Path::new("."));
            for line in text
                .lines()
                .map(str::trim)
                .filter(|l| !l.is_empty() && !l.starts_with('#'))
            {
                roots.push(base.join(line));
            }
        }
        let mut paths = Vec::new();
        for root in roots {
            if root.is_dir() {
                for entry in std::fs::read_dir(&root).map_err(Error::io(&root))? {
                    let path = entry.map_err(Error::io(&root))?.path();
                    if path.extension().is_some_and(|e| e == pgsum::pack::EXTENSION) && !is_hidden(&path) {
                        paths.push(path);
                    }
                }
            } else {
                paths.push(root);
            }
        }
        if !self.ids.is_empty() {
            for id in &self.ids {
                if !pgsum::scoring_file::is_pgs_id(id) && !pgsum::scoring_file::is_custom_id(id) {
                    return Err(Error::Invalid(format!(
                        "{id:?} is not a PGS Catalog or custom score ID"
                    )));
                }
            }
            paths.retain(|p| {
                p.file_stem()
                    .and_then(|s| s.to_str())
                    .is_some_and(|s| self.ids.iter().any(|id| id == s))
            });
            let found: Vec<_> = paths.iter().filter_map(|p| p.file_stem()?.to_str()).collect();
            let missing: Vec<_> = self.ids.iter().filter(|id| !found.contains(&id.as_str())).collect();
            if !missing.is_empty() {
                return Err(Error::Invalid(format!("no pack for {missing:?}")));
            }
        }
        paths.sort();
        paths.dedup();
        if paths.is_empty() {
            return Err(Error::Invalid("no packs selected; use --pack or --pack-list".into()));
        }
        Ok(paths)
    }
}

#[derive(Subcommand)]
enum Command {
    /// Compile PGS Catalog harmonized scoring files into packs.
    Compile {
        /// Scoring files: PGS Catalog harmonized files (`*_hmPOS_GRCh38.txt.gz`, or directories of them, each with
        /// its `<PGS_ID>.metadata.json` alongside), or custom scoring files (see DESIGN.md, "Custom scores").
        #[arg(required = true)]
        scoring_files: Vec<PathBuf>,
        /// GRCh38 reference FASTA (with `.fai`), used to resolve SNV orientation.
        #[arg(long)]
        reference: PathBuf,
        /// Directory to write packs into.
        #[arg(long)]
        out: PathBuf,
    },
    /// Download PGS Catalog scores and compile each into a pack, deleting downloads as it goes.
    ///
    /// Scores whose pack already exists with an identical Catalog record are skipped, so an interrupted
    /// run can be restarted. A log of every score is written to `<out>/fetch.tsv`.
    Fetch {
        /// These PGS IDs (comma-separated or repeated).
        #[arg(long, value_delimiter = ',', required_unless_present = "all", conflicts_with = "all")]
        ids: Vec<String>,
        /// Every score in the Catalog.
        #[arg(long)]
        all: bool,
        /// GRCh38 reference FASTA (with `.fai`).
        #[arg(long)]
        reference: PathBuf,
        /// Directory to write packs into.
        #[arg(long)]
        out: PathBuf,
        /// Where downloads go while they are compiled (default: `<out>/downloads`).
        #[arg(long)]
        downloads: Option<PathBuf>,
        /// Keep scoring files and Catalog records after compiling.
        #[arg(long)]
        keep_downloads: bool,
        /// Skip scores with more variants than this.
        #[arg(long)]
        max_variants: Option<u64>,
        /// Scores downloaded and compiled at once.
        #[arg(long, default_value_t = 4)]
        jobs: usize,
    },
    /// Print a pack's terms as TSV, or its header as JSON.
    Inspect {
        /// A pack, or a genotype table with `--header`.
        path: PathBuf,
        /// Print the header instead of the terms.
        #[arg(long)]
        header: bool,
        /// Add each term's status, call state and effect dosage from this genotype table.
        #[arg(long)]
        genotypes: Option<PathBuf>,
        /// With `--genotypes`: allow inferred other alleles, as `score` does.
        #[arg(long)]
        allow_inferred_other_allele: bool,
    },
    /// Read genotypes from a gVCF at every site the packs need.
    Extract {
        /// Bgzipped single-sample gVCF with a tabix or CSI index.
        #[arg(long)]
        gvcf: PathBuf,
        /// GRCh38 reference FASTA (with `.fai`).
        #[arg(long)]
        reference: PathBuf,
        #[command(flatten)]
        packs: PackSelection,
        /// Output genotype table (`.pgsg`).
        #[arg(long)]
        out: PathBuf,
        /// Target index to reuse when it matches the selected packs, or to create (`.pgst`).
        #[arg(long)]
        targets_cache: Option<PathBuf>,
    },
    /// Score packs against an extracted genotype table.
    Score {
        /// Genotype table from `extract`.
        #[arg(long)]
        genotypes: PathBuf,
        #[command(flatten)]
        packs: PackSelection,
        /// Output directory for one JSON result per score.
        #[arg(long)]
        out: PathBuf,
        /// Also write a per-term TSV for each score.
        #[arg(long)]
        terms: bool,
        /// Also score terms without an author other allele, using the orientation their pack inferred.
        #[arg(long)]
        allow_inferred_other_allele: bool,
    },
    /// Extract then score in one step.
    Run {
        #[arg(long)]
        gvcf: PathBuf,
        #[arg(long)]
        reference: PathBuf,
        #[command(flatten)]
        packs: PackSelection,
        #[arg(long)]
        out: PathBuf,
        #[arg(long)]
        terms: bool,
        /// Target index to reuse when it matches the selected packs, or to create (`.pgst`).
        #[arg(long)]
        targets_cache: Option<PathBuf>,
        /// Also score terms without an author other allele, using the orientation their pack inferred.
        #[arg(long)]
        allow_inferred_other_allele: bool,
    },
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let threads = threads(cli.threads);
    if let Err(e) = rayon::ThreadPoolBuilder::new().num_threads(threads).build_global() {
        eprintln!("pgsum: {e}");
        return ExitCode::FAILURE;
    }
    let result: Result<()> = match cli.command {
        Command::Compile {
            scoring_files,
            reference,
            out,
        } => compile(&scoring_files, &reference, &out),
        Command::Inspect {
            path,
            header,
            genotypes,
            allow_inferred_other_allele,
        } => inspect(
            &path,
            header,
            genotypes.as_deref(),
            &pgsum::score::Options {
                allow_inferred_other_allele,
            },
        ),
        Command::Fetch {
            ids,
            all: _,
            reference,
            out,
            downloads,
            keep_downloads,
            max_variants,
            jobs,
        } => fetch(&ids, &reference, &out, downloads, keep_downloads, max_variants, jobs),
        Command::Extract {
            gvcf,
            reference,
            packs,
            out,
            targets_cache,
        } => packs
            .resolve()
            .and_then(|packs| extract(&gvcf, &reference, &packs, &out, targets_cache.as_deref(), threads)),
        Command::Score {
            genotypes,
            packs,
            out,
            terms,
            allow_inferred_other_allele,
        } => packs.resolve().and_then(|packs| {
            let options = pgsum::score::Options {
                allow_inferred_other_allele,
            };
            score(&GenotypeTable::open(&genotypes)?, &packs, &out, terms, &options)
        }),
        Command::Run {
            gvcf,
            reference,
            packs,
            out,
            terms,
            targets_cache,
            allow_inferred_other_allele,
        } => packs.resolve().and_then(|packs| {
            let options = pgsum::score::Options {
                allow_inferred_other_allele,
            };
            std::fs::create_dir_all(&out).map_err(Error::io(&out))?;
            let table_path = out.join("genotypes.pgsg");
            extract(
                &gvcf,
                &reference,
                &packs,
                &table_path,
                targets_cache.as_deref(),
                threads,
            )?;
            score(&GenotypeTable::open(&table_path)?, &packs, &out, terms, &options)
        }),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("pgsum: {e}");
            ExitCode::FAILURE
        }
    }
}

fn threads(requested: Option<usize>) -> usize {
    requested
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()))
        .max(1)
}

/// Dot files, including the `._name` metadata files macOS writes next to files on exFAT and FAT drives.
fn is_hidden(path: &Path) -> bool {
    path.file_name()
        .and_then(|n| n.to_str())
        .is_some_and(|n| n.starts_with('.'))
}

/// Scoring files given directly, plus every `*_hmPOS_GRCh38.txt.gz` in given directories.
fn scoring_file_paths(inputs: &[PathBuf]) -> Result<Vec<PathBuf>> {
    let mut paths = Vec::new();
    for input in inputs {
        if input.is_dir() {
            for entry in std::fs::read_dir(input).map_err(Error::io(input))? {
                let path = entry.map_err(Error::io(input))?.path();
                if path.to_str().is_some_and(|p| p.ends_with("_hmPOS_GRCh38.txt.gz")) && !is_hidden(&path) {
                    paths.push(path);
                }
            }
        } else {
            paths.push(input.clone());
        }
    }
    paths.sort();
    paths.dedup();
    Ok(paths)
}

fn compile(inputs: &[PathBuf], reference: &Path, out: &Path) -> Result<()> {
    std::fs::create_dir_all(out).map_err(Error::io(out))?;
    let scoring_files = scoring_file_paths(inputs)?;
    let reference = Reference::open(reference)?;
    let identity = reference_identity(&reference)?;
    let failures: Vec<String> = scoring_files
        .par_iter()
        .filter_map(|path| {
            let started = std::time::Instant::now();
            let result = compile_file(path, None, &reference, &identity, out);
            match result {
                Ok(h) => {
                    eprintln!(
                        "{}: {} terms, {} resolved, inventory {} ({:.1}s)",
                        h.pgs_id,
                        h.inventory.actual_terms,
                        h.counts.orientation.get("resolved").copied().unwrap_or(0),
                        if h.inventory.consistent {
                            "consistent"
                        } else {
                            "INCONSISTENT"
                        },
                        started.elapsed().as_secs_f64()
                    );
                    None
                }
                Err(e) => Some(e.to_string()),
            }
        })
        .collect();
    if failures.is_empty() {
        Ok(())
    } else {
        Err(Error::Invalid(failures.join("\n")))
    }
}

fn fetch(
    ids: &[String],
    reference: &Path,
    out: &Path,
    downloads: Option<PathBuf>,
    keep_downloads: bool,
    max_variants: Option<u64>,
    jobs: usize,
) -> Result<()> {
    use pgsum::fetch::Outcome;
    let started = std::time::Instant::now();
    let reference = Reference::open(reference)?;
    let identity = reference_identity(&reference)?;
    let agent = pgsum::fetch::agent();
    let records = pgsum::fetch::catalog_records(&agent, ids)?;
    let variants: u64 = records.iter().filter_map(|r| r["variants_number"].as_u64()).sum();
    eprintln!("{} scores, {variants} variants in the Catalog records", records.len());
    let downloads = downloads.unwrap_or_else(|| out.join("downloads"));
    let options = pgsum::fetch::Options {
        reference: &reference,
        identity: &identity,
        out,
        downloads: &downloads,
        keep_downloads,
        max_variants,
    };
    let done = std::sync::atomic::AtomicUsize::new(0);
    let total = records.len();
    let report = |l: &pgsum::fetch::Logged| {
        let n = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
        let what = match &l.outcome {
            Outcome::Compiled { terms, .. } => format!("compiled, {terms} terms"),
            Outcome::Unchanged => "unchanged".into(),
            Outcome::Skipped(why) => format!("skipped: {why}"),
            Outcome::Failed(why) => format!("FAILED: {why}"),
        };
        eprintln!("[{n}/{total}] {}: {what} ({:.1}s)", l.id, l.seconds);
    };
    let mut log = pgsum::fetch::fetch_all(&agent, &records, &options, jobs, &report)?;
    log.sort_by(|a, b| a.id.cmp(&b.id));
    let mut tsv = String::from("pgs_id\toutcome\tterms\tsource_bytes\tseconds\tdetail\n");
    let mut failed = 0;
    for l in &log {
        let (outcome, terms, bytes, detail) = match &l.outcome {
            Outcome::Compiled { terms, source_bytes } => ("compiled", terms.to_string(), source_bytes.to_string(), ""),
            Outcome::Unchanged => ("unchanged", String::new(), String::new(), ""),
            Outcome::Skipped(why) => ("skipped", String::new(), String::new(), why.as_str()),
            Outcome::Failed(why) => {
                failed += 1;
                ("failed", String::new(), String::new(), why.as_str())
            }
        };
        tsv.push_str(&format!(
            "{}\t{outcome}\t{terms}\t{bytes}\t{:.1}\t{}\n",
            l.id,
            l.seconds,
            detail.replace(['\t', '\n'], " ")
        ));
    }
    let log_path = out.join("fetch.tsv");
    std::fs::write(&log_path, tsv).map_err(Error::io(&log_path))?;
    eprintln!(
        "done in {:.0}s; log in {}",
        started.elapsed().as_secs_f64(),
        log_path.display()
    );
    if failed > 0 {
        return Err(Error::Invalid(format!(
            "{failed} of {total} scores failed; see {}",
            log_path.display()
        )));
    }
    Ok(())
}

fn extract(
    gvcf: &Path,
    reference: &Path,
    packs: &[PathBuf],
    out: &Path,
    targets_cache: Option<&Path>,
    threads: usize,
) -> Result<()> {
    let started = std::time::Instant::now();
    let reference = Reference::open(reference)?;
    let identity = reference_identity(&reference)?;
    let loaded = started.elapsed().as_secs_f64();
    let (table, timings) = pgsum::extract::extract(gvcf, &reference, &identity, packs, targets_cache, threads)?;
    let t = std::time::Instant::now();
    table.write(out)?;
    eprintln!(
        "  timings: reference {loaded:.1}s, targets {:.1}s{}, scan {:.1}s, assess {:.1}s, table {:.1}s, write {:.1}s",
        timings.targets_s,
        if timings.targets_from_cache {
            " (from index)"
        } else {
            ""
        },
        timings.scan_s,
        timings.assess_s,
        timings.table_s,
        t.elapsed().as_secs_f64()
    );
    let h = &table.header;
    eprintln!(
        "{}: {} records scanned, {} kept, {} targets in {:.1}s",
        h.sample.sample_id,
        h.records_scanned,
        h.records_kept,
        h.targets,
        started.elapsed().as_secs_f64()
    );
    for (state, n) in &h.states {
        eprintln!("  {state}: {n}");
    }
    Ok(())
}

/// Score every pack against the table: packs in parallel, and each pack's terms in parallel chunks.
fn score(
    table: &GenotypeTable,
    packs: &[PathBuf],
    out: &Path,
    terms: bool,
    options: &pgsum::score::Options,
) -> Result<()> {
    std::fs::create_dir_all(out).map_err(Error::io(out))?;
    let started = std::time::Instant::now();
    let outcomes: Vec<(&PathBuf, Result<pgsum::score::ScoreResult>)> = packs
        .par_iter()
        .map(|path| {
            let result = (|| -> Result<pgsum::score::ScoreResult> {
                let pack = Pack::open(path)?;
                let id = &pack.header.pgs_id;
                let result = if terms {
                    let tsv = out.join(format!("{id}.terms.tsv"));
                    let mut w = BufWriter::new(std::fs::File::create(&tsv).map_err(Error::io(&tsv))?);
                    let r = pgsum::score::score(&pack, table, options, Some(&mut w))?;
                    w.flush().map_err(Error::io(&tsv))?;
                    r
                } else {
                    pgsum::score::score(&pack, table, options, None)?
                };
                let json_path = out.join(format!("{id}.score.json"));
                let mut json = serde_json::to_vec_pretty(&result).map_err(|e| Error::Invalid(e.to_string()))?;
                json.push(b'\n');
                std::fs::write(&json_path, json).map_err(Error::io(&json_path))?;
                Ok(result)
            })();
            (path, result)
        })
        .collect();
    let mut results = Vec::new();
    let mut failures = Vec::new();
    for (path, r) in outcomes {
        match r {
            Ok(r) => results.push(r),
            Err(e) => failures.push(format!("{}: {e}", path.display())),
        }
    }
    results.sort_by(|a, b| a.pgs_id.cmp(&b.pgs_id));
    let mut summary = String::from(
        "pgs_id\tstatus\traw_score\tpartial_raw_score\tscorable_terms\ttotal_terms\tterm_coverage\tweight_coverage\n",
    );
    for r in &results {
        let c = &r.partial.coverage;
        summary.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\n",
            r.pgs_id,
            r.status,
            r.raw_score.as_deref().unwrap_or(""),
            r.partial.raw_score,
            c.scorable_terms,
            c.total_terms,
            c.term_fraction,
            c.weight_fraction
        ));
        eprintln!(
            "{}: {} (partial {} over {:.2}% of terms, {:.2}% of weight)",
            r.pgs_id,
            r.status,
            r.partial.raw_score,
            100.0 * c.term_fraction,
            100.0 * c.weight_fraction
        );
    }
    let summary_path = out.join("scores.tsv");
    std::fs::write(&summary_path, summary).map_err(Error::io(&summary_path))?;
    let terms_scored: u64 = results.iter().map(|r| r.required_terms).sum();
    eprintln!(
        "scored {} packs, {} terms, in {:.1}s",
        results.len(),
        terms_scored,
        started.elapsed().as_secs_f64()
    );
    if failures.is_empty() {
        Ok(())
    } else {
        Err(Error::Invalid(failures.join("\n")))
    }
}

fn inspect(path: &Path, header: bool, genotypes: Option<&Path>, options: &pgsum::score::Options) -> Result<()> {
    let mut out = BufWriter::new(std::io::stdout().lock());
    let json = |out: &mut BufWriter<_>, value: &dyn erased::Json| value.write(out);
    let is_table = std::fs::File::open(path)
        .and_then(|mut f| {
            let mut magic = [0u8; 8];
            std::io::Read::read_exact(&mut f, &mut magic).map(|_| &magic == pgsum::genotypes::MAGIC)
        })
        .map_err(Error::io(path))?;
    if is_table {
        let table = GenotypeTable::open(path)?;
        json(&mut out, &table.header)?;
    } else {
        let pack = Pack::open(path)?;
        if header {
            json(&mut out, &pack.header)?;
        } else {
            match genotypes {
                Some(g) => {
                    pgsum::score::score(&pack, &GenotypeTable::open(g)?, options, Some(&mut out))?;
                }
                None => pack.write_terms_tsv(&mut out)?,
            }
        }
    }
    out.flush().map_err(Error::io("<stdout>"))
}

mod erased {
    use pgsum::{Error, Result};

    pub trait Json {
        fn write(&self, out: &mut dyn std::io::Write) -> Result<()>;
    }

    impl<T: serde::Serialize> Json for T {
        fn write(&self, out: &mut dyn std::io::Write) -> Result<()> {
            serde_json::to_writer_pretty(&mut *out, self).map_err(|e| Error::Invalid(e.to_string()))?;
            writeln!(out).map_err(Error::io("<stdout>"))
        }
    }
}
