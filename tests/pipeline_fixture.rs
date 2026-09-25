//! Extract the synthetic gVCF, score it, and compare every term's status, call state, effect dosage and
//! contribution, and the scores, with the expected output.
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
    pgsum::score::score(&pack, &table, Some(&mut tsv)).unwrap();
    let actual: Vec<String> = String::from_utf8(tsv)
        .unwrap()
        .lines()
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            [f[0], f[13], f[14], f[15], f[16]].join("\t")
        })
        .collect();
    let expected = std::fs::read_to_string(fixtures.join("PGS999998.expected.tsv")).unwrap();
    let expected: Vec<&str> = expected.lines().collect();
    assert_eq!(actual, expected);
    std::fs::remove_dir_all(&out).unwrap();
}

/// PGS999997 has every weight model on passing calls, so the complete score is reported. PGS999998 is
/// withheld; its partial score sums its 10 scorable terms. Expected sums come from the reference
/// implementation's exact decimal reduction.
#[test]
fn synthetic_scores_match_reference_implementation() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let out = std::env::temp_dir().join(format!("pgsum-score-test-{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let reference = Reference::open(&fixtures.join("synthetic.fa")).unwrap();
    let identity = reference_identity(&reference).unwrap();
    let mut paths = Vec::new();
    for id in ["PGS999997", "PGS999998"] {
        compile_file(
            &fixtures.join(format!("{id}_hmPOS_GRCh38.txt.gz")),
            &fixtures.join(format!("{id}.metadata.json")),
            &reference,
            &identity,
            &out,
        )
        .unwrap();
        paths.push(out.join(format!("{id}.pgsp")));
    }
    let (table, _) = extract(&fixtures.join("synthetic.g.vcf.gz"), &reference, &identity, &paths, 2).unwrap();

    let complete = Pack::open(&paths[0]).unwrap();
    let mut tsv = Vec::new();
    let result = pgsum::score::score(&complete, &table, Some(&mut tsv)).unwrap();
    assert_eq!(result.status, "complete_uncalibrated_score");
    assert_eq!(result.raw_score.as_deref(), Some("1332.64750200000000000000000"));
    assert_eq!(result.partial.raw_score, "1332.64750200000000000000000");
    assert!(result.withheld_because.is_empty());
    let actual: Vec<String> = String::from_utf8(tsv)
        .unwrap()
        .lines()
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            [f[0], f[13], f[14], f[15], f[16]].join("\t")
        })
        .collect();
    let expected = std::fs::read_to_string(fixtures.join("PGS999997.expected.tsv")).unwrap();
    assert_eq!(actual, expected.lines().collect::<Vec<_>>());

    let withheld = pgsum::score::score(&Pack::open(&paths[1]).unwrap(), &table, None).unwrap();
    assert_eq!(withheld.status, "score_withheld");
    assert_eq!(withheld.raw_score, None);
    assert_eq!(withheld.partial.raw_score, "7.3");
    assert_eq!((withheld.scorable_terms, withheld.required_terms), (10, 26));
    std::fs::remove_dir_all(&out).unwrap();
}
