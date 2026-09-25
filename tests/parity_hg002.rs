//! Parity with an external reference implementation of the same rules, on GIAB HG002 (GRCh38).
//!
//! Needs local data, so it is ignored by default. Set `PGSUM_PARITY_DIR` to a directory holding the HG002
//! gVCF and index, the GRCh38 FASTA and `.fai`, the scoring files, and the reference implementation's
//! per-term output, then run `cargo test --release -- --ignored parity`.
//!
//! Pass criteria: for every score, each term's state and effect dosage match, and the raw score is
//! numerically equal (string formatting may differ).

#[test]
#[ignore = "needs local HG002 data; set PGSUM_PARITY_DIR"]
fn parity_hg002() {
    let dir = std::env::var("PGSUM_PARITY_DIR").expect("PGSUM_PARITY_DIR is not set");
    assert!(std::path::Path::new(&dir).is_dir(), "{dir} is not a directory");
    todo!("run pgsum on the HG002 inputs and compare with the reference output");
}
