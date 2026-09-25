//! Extract the synthetic gVCF and compare every term's status, call state and effect dosage with the
//! expected output.
//!
//! `tests/fixtures/synthetic.g.vcf.gz` has one scenario per genotype rule (see `make_genotype_fixture.py`) and
//! `PGS999998` has one SNV term at each. `PGS999998.expected.tsv` is the output of an independent reference
//! implementation of the same rules.

use std::path::Path;

use pgsum::compile::{compile_file, reference_identity};
use pgsum::extract::extract;
use pgsum::genotypes::GenotypeTable;
use pgsum::pack::Pack;
use pgsum::reference::Reference;

#[test]
fn synthetic_gvcf_matches_reference_implementation() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let out = std::env::temp_dir().join(format!("pgsum-extract-test-{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let reference = Reference::open(&fixtures.join("synthetic.fa")).unwrap();
    let identity = reference_identity(&reference).unwrap();
    compile_file(
        &fixtures.join("PGS999998_hmPOS_GRCh38.txt.gz"),
        &fixtures.join("PGS999998.metadata.json"),
        &reference,
        &identity,
        &out,
    )
    .unwrap();
    let pack_path = out.join("PGS999998.pgsp");
    let (table, _) = extract(
        &fixtures.join("synthetic.g.vcf.gz"),
        &reference,
        &identity,
        std::slice::from_ref(&pack_path),
        2,
    )
    .unwrap();
    let pack = Pack::open(&pack_path).unwrap();
    assert_eq!(table.header.sample.sample_id, "SYNTH");
    assert!(table.header.sample.refcall_defined);
    assert_eq!(table.header.records_scanned, 23);

    // The table survives a write and read.
    let path = out.join("synthetic.pgsg");
    table.write(&path).unwrap();
    let table = GenotypeTable::open(&path).unwrap();

    let mut tsv = Vec::new();
    pack.write_terms_tsv(&mut tsv, Some(&table)).unwrap();
    let actual: Vec<String> = String::from_utf8(tsv)
        .unwrap()
        .lines()
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            [f[0], f[13], f[14], f[15]].join("\t")
        })
        .collect();
    let expected = std::fs::read_to_string(fixtures.join("PGS999998.expected.tsv")).unwrap();
    let expected: Vec<&str> = expected.lines().collect();
    assert_eq!(actual, expected);
    std::fs::remove_dir_all(&out).unwrap();
}
