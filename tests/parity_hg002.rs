//! Parity with an external reference implementation of the same rules, on real data (GIAB HG002).
//!
//! Needs local data, so it is ignored by default. Set `PGSUM_PARITY_DIR` to a directory holding:
//!
//! - `genotypes.pgsg`: a genotype table extracted with the packs below;
//! - `packs/<pgs_id>.pgsp`: the packs;
//! - `<pgs_id>.expected.tsv`: the reference implementation's per-term output with the columns ordinal, status,
//!   call_state, effect_dosage, contribution;
//! - optionally `<pgs_id>.expected.partial`: its exact partial sum.
//!
//! Then run `cargo test --release -- --ignored parity`. Every term and every sum must match exactly.

use std::path::Path;

use pgsum::genotypes::GenotypeTable;
use pgsum::pack::Pack;

#[test]
#[ignore = "needs local HG002 data; set PGSUM_PARITY_DIR"]
fn parity_hg002() {
    let dir = std::env::var("PGSUM_PARITY_DIR").expect("PGSUM_PARITY_DIR is not set");
    let dir = Path::new(&dir);
    let table = GenotypeTable::open(&dir.join("genotypes.pgsg")).unwrap();
    let mut checked = 0;
    for entry in std::fs::read_dir(dir.join("packs")).unwrap() {
        let path = entry.unwrap().path();
        let pack = Pack::open(&path).unwrap();
        let id = pack.header.pgs_id.clone();
        let Ok(expected) = std::fs::read_to_string(dir.join(format!("{id}.expected.tsv"))) else {
            continue;
        };
        let mut tsv = Vec::new();
        let result = pgsum::score::score(&pack, &table, Some(&mut tsv)).unwrap();
        let tsv = String::from_utf8(tsv).unwrap();
        let mut terms = 0;
        for (i, (a, e)) in tsv.lines().zip(expected.lines()).enumerate() {
            let f: Vec<&str> = a.split('\t').collect();
            assert_eq!([f[0], f[13], f[14], f[15], f[16]].join("\t"), e, "{id} line {}", i + 1);
            terms += 1;
        }
        assert_eq!(tsv.lines().count(), expected.lines().count(), "{id} term count");
        if let Ok(partial) = std::fs::read_to_string(dir.join(format!("{id}.expected.partial"))) {
            assert_eq!(result.partial.raw_score, partial.trim(), "{id} partial sum");
        }
        eprintln!("{id}: {} terms identical", terms - 1);
        checked += 1;
    }
    assert!(checked > 0, "no <pgs_id>.expected.tsv files in {}", dir.display());
}
