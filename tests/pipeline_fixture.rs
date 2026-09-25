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
        None,
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
    pgsum::score::score(&pack, &table, &Default::default(), Some(&mut tsv)).unwrap();
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
    let (table, _) = extract(
        &fixtures.join("synthetic.g.vcf.gz"),
        &reference,
        &identity,
        &paths,
        None,
        2,
    )
    .unwrap();

    let complete = Pack::open(&paths[0]).unwrap();
    let mut tsv = Vec::new();
    let result = pgsum::score::score(&complete, &table, &Default::default(), Some(&mut tsv)).unwrap();
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

    let withheld = pgsum::score::score(&Pack::open(&paths[1]).unwrap(), &table, &Default::default(), None).unwrap();
    assert_eq!(withheld.status, "score_withheld");
    assert_eq!(withheld.raw_score, None);
    assert_eq!(withheld.partial.raw_score, "7.3");
    assert_eq!((withheld.scorable_terms, withheld.required_terms), (10, 26));
    std::fs::remove_dir_all(&out).unwrap();
}

/// A target index is written on the first run, reused when the packs are unchanged (giving an identical
/// table), and ignored when the pack set changes.
#[test]
fn target_index_is_reused_only_for_the_same_packs() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let out = std::env::temp_dir().join(format!("pgsum-targets-test-{}", std::process::id()));
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
    let gvcf = fixtures.join("synthetic.g.vcf.gz");
    let cache = out.join("targets.pgst");
    let (first, t1) = extract(&gvcf, &reference, &identity, &paths, Some(&cache), 2).unwrap();
    assert!(!t1.targets_from_cache && cache.exists());
    let (second, t2) = extract(&gvcf, &reference, &identity, &paths, Some(&cache), 2).unwrap();
    assert!(t2.targets_from_cache);
    assert_eq!(first.header, second.header);
    let (third, t3) = extract(&gvcf, &reference, &identity, &paths[..1], Some(&cache), 2).unwrap();
    assert!(!t3.targets_from_cache);
    assert_eq!(third.header.packs.len(), 1);
    std::fs::remove_dir_all(&out).unwrap();
}

/// Scores without an author other allele: unscorable by default; with `allow_inferred_other_allele`, oriented
/// by the score's reference convention (PGS999996 effect is ALT, PGS999995 effect is REF) or, with no
/// convention, by the Catalog's inferred allele when it names one base that orients unambiguously (PGS999994).
/// Expected outcomes follow from the synthetic gVCF scenarios at each position.
#[test]
fn inferred_other_alleles_are_opt_in() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let out = std::env::temp_dir().join(format!("pgsum-infer-test-{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let reference = Reference::open(&fixtures.join("synthetic.fa")).unwrap();
    let identity = reference_identity(&reference).unwrap();
    let ids = ["PGS999996", "PGS999995", "PGS999994"];
    let mut paths = Vec::new();
    for id in ids {
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
    let (table, _) = extract(
        &fixtures.join("synthetic.g.vcf.gz"),
        &reference,
        &identity,
        &paths,
        None,
        2,
    )
    .unwrap();
    let allow = pgsum::score::Options {
        allow_inferred_other_allele: true,
    };

    // Per score: reference convention, then per term "status, effect dosage, contribution, inference method",
    // then the strict status and the partial sum.
    let expected = [
        (
            "PGS999996",
            Some("effect_is_alt"),
            vec![
                "scorable_observation\t0\t0.0\treference_anchored_effect_is_alt",
                "scorable_observation\t1\t0.2\treference_anchored_effect_is_alt",
                "scorable_observation\t2\t0.6\treference_anchored_effect_is_alt",
                "other_called_allele\t\t\treference_anchored_effect_is_alt",
                "scorable_observation\t1\t0.5\treference_anchored_effect_is_alt",
                "scorable_observation\t0\t0.0\treference_anchored_effect_is_alt",
            ],
            "score_withheld",
            "1.3",
        ),
        (
            "PGS999995",
            Some("effect_is_ref"),
            vec![
                "scorable_observation\t2\t0.2\treference_anchored_effect_is_ref",
                "scorable_observation\t1\t0.2\treference_anchored_effect_is_ref",
                "scorable_observation\t0\t0.0\treference_anchored_effect_is_ref",
                "scorable_observation\t2\t0.8\treference_anchored_effect_is_ref",
                "scorable_observation\t1\t0.5\treference_anchored_effect_is_ref",
                "scorable_observation\t1\t0.6\treference_anchored_effect_is_ref",
            ],
            "complete_uncalibrated_score",
            "2.3",
        ),
        (
            "PGS999994",
            None,
            vec![
                "scorable_observation\t1\t0.2\tcatalog_inferred_other_allele",
                "scorable_observation\t0\t0.0\tcatalog_inferred_other_allele",
                "model_term_requires_review\t\t\t",
                "model_term_requires_review\t\t\t",
                "model_term_requires_review\t\t\t",
            ],
            "score_withheld",
            "0.2",
        ),
    ];
    for (path, (id, convention, rows, status, partial)) in paths.iter().zip(expected) {
        let pack = Pack::open(path).unwrap();
        assert_eq!(pack.header.pgs_id, id);
        let default = pgsum::score::score(&pack, &table, &Default::default(), None).unwrap();
        assert_eq!(default.scorable_terms, 0, "{id}: nothing is scorable by default");
        assert_eq!(default.partial.raw_score, "0");

        let mut tsv = Vec::new();
        let result = pgsum::score::score(&pack, &table, &allow, Some(&mut tsv)).unwrap();
        assert_eq!(
            result.inferred_other_allele.reference_convention.as_deref(),
            convention,
            "{id}"
        );
        assert_eq!(result.status, status, "{id}");
        assert_eq!(result.partial.raw_score, partial, "{id}");
        let actual: Vec<String> = String::from_utf8(tsv)
            .unwrap()
            .lines()
            .skip(1)
            .map(|line| {
                let f: Vec<&str> = line.split('\t').collect();
                [f[13], f[15], f[16], f[17]].join("\t")
            })
            .collect();
        assert_eq!(actual, rows, "{id}");
    }
    std::fs::remove_dir_all(&out).unwrap();
}
