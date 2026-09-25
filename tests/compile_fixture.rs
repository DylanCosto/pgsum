//! Compile the synthetic scoring file and compare every term with the expected output.
//!
//! `tests/fixtures/PGS999999_hmPOS_GRCh38.txt.gz` has one or more rows per description and orientation rule
//! (see `DESIGN.md`), against a synthetic 60 bp chr1 and 30 bp chrX. `PGS999999.expected.tsv` is the output
//! of an independent reference implementation of the same rules.

use std::path::Path;

use pgsum::compile::{compile_file, reference_identity};
use pgsum::pack::Pack;
use pgsum::reference::Reference;

#[test]
fn synthetic_scoring_file_matches_reference_implementation() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let out = std::env::temp_dir().join(format!("pgsum-test-{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let reference = Reference::open(&fixtures.join("synthetic.fa")).unwrap();
    let identity = reference_identity(&reference).unwrap();
    let header = compile_file(
        &fixtures.join("PGS999999_hmPOS_GRCh38.txt.gz"),
        &fixtures.join("PGS999999.metadata.json"),
        &reference,
        &identity,
        &out,
    )
    .unwrap();
    assert!(header.inventory.consistent);
    assert_eq!(header.counts.exact_duplicate_terms, 2);
    assert_eq!(header.counts.centred_terms, 1);

    let pack = Pack::open(&out.join("PGS999999.pgsp")).unwrap();
    assert_eq!(pack.header, header);
    let mut tsv = Vec::new();
    pack.write_terms_tsv(&mut tsv).unwrap();
    let expected = std::fs::read_to_string(fixtures.join("PGS999999.expected.tsv")).unwrap();
    let actual = String::from_utf8(tsv).unwrap();
    for (i, (a, e)) in actual.lines().zip(expected.lines()).enumerate() {
        assert_eq!(a, e, "line {}", i + 1);
    }
    assert_eq!(actual.lines().count(), expected.lines().count());
    std::fs::remove_dir_all(&out).unwrap();
}
