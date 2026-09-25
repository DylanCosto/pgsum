use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use clap::{Parser, Subcommand};
use pgsum::compile::{compile_file, metadata_path, reference_identity};
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
        pack: PathBuf,
        /// Print the header instead of the terms.
        #[arg(long)]
        header: bool,
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
        /// Output genotype table.
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
        Command::Inspect { pack, header } => inspect(&pack, header),
        Command::Extract { .. } => Err(Error::NotImplemented("extract")),
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

fn inspect(path: &Path, header: bool) -> Result<()> {
    let pack = Pack::open(path)?;
    let mut out = BufWriter::new(std::io::stdout().lock());
    if header {
        serde_json::to_writer_pretty(&mut out, &pack.header).map_err(|e| Error::Invalid(e.to_string()))?;
        writeln!(out).map_err(Error::io("<stdout>"))?;
    } else {
        pack.write_terms_tsv(&mut out)?;
    }
    out.flush().map_err(Error::io("<stdout>"))
}
