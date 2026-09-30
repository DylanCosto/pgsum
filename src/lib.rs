//! Polygenic score calculation from single-sample gVCFs.
//!
//! The pipeline is `compile` (scoring file → pack), `extract` (gVCF → genotype table at the pack sites)
//! and `score` (packs × genotypes → one result per score). See `DESIGN.md` for the rules.

pub mod alleles;
pub mod batch;
pub mod bcf;
pub mod catalog_cache;
pub mod cohort;
mod cohort_diagnostics;
pub mod compare;
pub mod compile;
pub mod decimal;
pub mod digest;
pub mod dosage;
pub mod evidence;
pub mod extract;
pub mod fetch;
pub mod fill;
pub mod frequencies;
pub mod genotype;
pub mod genotypes;
pub mod gvcf;
pub mod index;
pub mod orient;
pub mod pack;
pub mod panel;
pub mod placement;
pub mod public;
pub mod reference;
pub mod score;
pub mod scoring_file;
pub mod targets;
pub mod term;

/// Errors surfaced by the library and the CLI.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("{0} is not implemented yet")]
    NotImplemented(&'static str),
    #[error("{path}: {source}")]
    Io {
        path: std::path::PathBuf,
        source: std::io::Error,
    },
    #[error("{0}")]
    Invalid(String),
}

impl Error {
    pub fn io(path: impl Into<std::path::PathBuf>) -> impl FnOnce(std::io::Error) -> Error {
        let path = path.into();
        move |source| Error::Io { path, source }
    }
}

pub type Result<T> = std::result::Result<T, Error>;

/// `Err(Error::Invalid(format!(...)))`.
#[macro_export]
macro_rules! invalid {
    ($($arg:tt)*) => { Err($crate::Error::Invalid(format!($($arg)*))) };
}

pub mod workflow;

pub mod missingness;
