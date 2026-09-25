use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use clap::{Parser, Subcommand};
use pgsum::compile::{compile_file, metadata_path, reference_identity};
use pgsum::genotypes::GenotypeTable;
use pgsum::pack::Pack;
use pgsum::reference::Reference;
use pgsum::{Error, Result};

/// Polygenic score calculation from single-sample gVCFs.
#[derive(Parser)]
#[command(version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Compile PGS Catalog harmonized scoring files into packs.
    Compile {
        /// Harmonized scoring files (`*_hmPOS_GRCh38.txt.gz`).
        #[arg(required = true)]
        scoring_files: Vec<PathBuf>,
        /// GRCh38 reference FASTA (with `.fai`), used to resolve SNV orientation.
        #[arg(long)]
        reference: PathBuf,
        /// Directory to write packs into.
        #[arg(long)]
        out: PathBuf,
        /// Scoring files compiled at once (default: all cores).
        #[arg(long)]
        threads: Option<usize>,
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
    },
    /// Read genotypes from a gVCF at every site the packs need.
    Extract {
        /// Bgzipped single-sample gVCF with a tabix or CSI index.
        #[arg(long)]
        gvcf: PathBuf,
        /// GRCh38 reference FASTA (with `.fai`).
        #[arg(long)]
        reference: PathBuf,
        /// Packs whose sites to extract.
        #[arg(long = "pack", required = true)]
        packs: Vec<PathBuf>,
        /// Output genotype table (`.pgsg`).
        #[arg(long)]
        out: PathBuf,
        /// Worker threads (default: all cores).
        #[arg(long)]
        threads: Option<usize>,
    },
    /// Score packs against an extracted genotype table.
    Score {
        /// Genotype table from `extract`.
        #[arg(long)]
        genotypes: PathBuf,
        #[arg(long = "pack", required = true)]
        packs: Vec<PathBuf>,
        /// Output directory for one JSON result per score.
        #[arg(long)]
        out: PathBuf,
        /// Also write a per-term TSV for each score.
        #[arg(long)]
        terms: bool,
    },
    /// Extract then score in one step.
    Run {
        #[arg(long)]
        gvcf: PathBuf,
        #[arg(long)]
        reference: PathBuf,
        #[arg(long = "pack", required = true)]
        packs: Vec<PathBuf>,
        #[arg(long)]
        out: PathBuf,
        #[arg(long)]
        threads: Option<usize>,
        #[arg(long)]
        terms: bool,
    },
}

fn main() -> ExitCode {
    let result: Result<()> = match Cli::parse().command {
        Command::Compile {
            scoring_files,
            reference,
            out,
            threads,
        } => compile(&scoring_files, &reference, &out, threads),
        Command::Inspect {
            path,
            header,
            genotypes,
        } => inspect(&path, header, genotypes.as_deref()),
        Command::Extract {
            gvcf,
            reference,
            packs,
            out,
            threads: t,
        } => extract(&gvcf, &reference, &packs, &out, t),
        Command::Score { .. } => Err(Error::NotImplemented("score")),
        Command::Run { .. } => Err(Error::NotImplemented("run")),
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

fn compile(scoring_files: &[PathBuf], reference: &Path, out: &Path, requested: Option<usize>) -> Result<()> {
    std::fs::create_dir_all(out).map_err(Error::io(out))?;
    let reference = Reference::open(reference)?;
    let identity = reference_identity(&reference)?;
    let next = AtomicUsize::new(0);
    let failures = Mutex::new(Vec::new());
    std::thread::scope(|scope| {
        for _ in 0..threads(requested).min(scoring_files.len()) {
            scope.spawn(|| {
                while let Some(path) = scoring_files.get(next.fetch_add(1, Ordering::Relaxed)) {
                    let started = std::time::Instant::now();
                    let result = metadata_path(path)
                        .ok_or_else(|| Error::Invalid(format!("{}: cannot name its metadata file", path.display())))
                        .and_then(|metadata| compile_file(path, &metadata, &reference, &identity, out));
                    match result {
                        Ok(h) => eprintln!(
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
                        ),
                        Err(e) => failures.lock().unwrap().push(e.to_string()),
                    }
                }
            });
        }
    });
    let failures = failures.into_inner().unwrap();
    if failures.is_empty() {
        Ok(())
    } else {
        Err(Error::Invalid(failures.join("\n")))
    }
}

fn extract(gvcf: &Path, reference: &Path, packs: &[PathBuf], out: &Path, requested: Option<usize>) -> Result<()> {
    let started = std::time::Instant::now();
    let reference = Reference::open(reference)?;
    let identity = reference_identity(&reference)?;
    let loaded = started.elapsed().as_secs_f64();
    let (table, timings) = pgsum::extract::extract(gvcf, &reference, &identity, packs, threads(requested))?;
    let t = std::time::Instant::now();
    table.write(out)?;
    eprintln!(
        "  timings: reference {loaded:.1}s, packs+targets {:.1}s, scan {:.1}s, assess {:.1}s, table {:.1}s, write {:.1}s",
        timings.targets_s,
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

fn inspect(path: &Path, header: bool, genotypes: Option<&Path>) -> Result<()> {
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
            let table = genotypes.map(GenotypeTable::open).transpose()?;
            pack.write_terms_tsv(&mut out, table.as_ref())?;
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
