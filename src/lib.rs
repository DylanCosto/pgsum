//! Polygenic score calculation from single-sample gVCFs.
//!
//! The pipeline is `compile` (scoring file → pack), `extract` (gVCF → genotype table at the pack sites)
//! and `score` (packs × genotypes → one result per score). See `DESIGN.md` for the rules.

pub mod genotype;

/// Errors surfaced by the library and the CLI.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("{0} is not implemented yet")]
    NotImplemented(&'static str),
}
