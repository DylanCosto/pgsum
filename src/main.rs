use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Parser, Subcommand};
use pgsum::Error;

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
    let result: Result<(), Error> = match Cli::parse().command {
        Command::Compile { .. } => Err(Error::NotImplemented("compile")),
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
