//! Extract the synthetic gVCF, score it, and compare every term's status, call state, effect dosage and
//! contribution, and the scores, with the expected output.
//!
//! `tests/fixtures/synthetic.g.vcf.gz` has one scenario per genotype rule (see `make_genotype_fixture.py`) and
//! `PGS999998` has one SNV term at each. `PGS999998.expected.tsv` is the output of an independent reference
//! implementation of the same rules.

use std::path::Path;

use pgsum::compile::{compile_file, reference_identity};
use pgsum::extract::{Options, ScanMode, extract};

fn opts(targets_cache: Option<&Path>, haploid_xy_as_homozygous: bool, scan: ScanMode) -> Options<'_> {
    Options {
        targets_cache,
        haploid_xy_as_homozygous,
        sample: None,
        scan,
        threads: 2,
    }
}
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
        Some(&fixtures.join("PGS999998.metadata.json")),
        &reference,
        &identity,
        None,
        &out,
    )
    .unwrap();
    let pack_path = out.join("PGS999998.pgsp");
    let (table, _) = extract(
        &fixtures.join("synthetic.g.vcf.gz"),
        &reference,
        &identity,
        std::slice::from_ref(&pack_path),
        &opts(None, false, ScanMode::Full),
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
            Some(&fixtures.join(format!("{id}.metadata.json"))),
            &reference,
            &identity,
            None,
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
        &opts(None, false, ScanMode::Full),
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
            Some(&fixtures.join(format!("{id}.metadata.json"))),
            &reference,
            &identity,
            None,
            &out,
        )
        .unwrap();
        paths.push(out.join(format!("{id}.pgsp")));
    }
    let gvcf = fixtures.join("synthetic.g.vcf.gz");
    let cache = out.join("targets.pgst");
    let (first, t1) = extract(
        &gvcf,
        &reference,
        &identity,
        &paths,
        &opts(Some(&cache), false, ScanMode::Full),
    )
    .unwrap();
    assert!(!t1.targets_from_cache && cache.exists());
    let (second, t2) = extract(
        &gvcf,
        &reference,
        &identity,
        &paths,
        &opts(Some(&cache), false, ScanMode::Full),
    )
    .unwrap();
    assert!(t2.targets_from_cache);
    assert_eq!(first.header, second.header);
    let (third, t3) = extract(
        &gvcf,
        &reference,
        &identity,
        &paths[..1],
        &opts(Some(&cache), false, ScanMode::Full),
    )
    .unwrap();
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
            Some(&fixtures.join(format!("{id}.metadata.json"))),
            &reference,
            &identity,
            None,
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
        &opts(None, false, ScanMode::Full),
    )
    .unwrap();
    let allow = pgsum::score::Options {
        allow_inferred_other_allele: true,
        ..Default::default()
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

/// A custom scoring file (plain text, `chr1` and `1` both used, no Catalog metadata) compiles, extracts and
/// scores like a Catalog one; its strict score needs no Catalog publication record.
#[test]
fn custom_scoring_file() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let out = std::env::temp_dir().join(format!("pgsum-custom-test-{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let reference = Reference::open(&fixtures.join("synthetic.fa")).unwrap();
    let identity = reference_identity(&reference).unwrap();
    let header = compile_file(&fixtures.join("MY_SCORE.tsv"), None, &reference, &identity, None, &out).unwrap();
    assert_eq!(header.origin, pgsum::scoring_file::Origin::Custom);
    assert_eq!(header.license.as_deref(), Some("CC0 1.0 (synthetic test fixture)"));
    assert!(header.inventory.consistent);
    let path = out.join("MY_SCORE.pgsp");
    let (table, _) = extract(
        &fixtures.join("synthetic.g.vcf.gz"),
        &reference,
        &identity,
        std::slice::from_ref(&path),
        &opts(None, false, ScanMode::Full),
    )
    .unwrap();
    let result = pgsum::score::score(&Pack::open(&path).unwrap(), &table, &Default::default(), None).unwrap();
    // 0.5×1 − 0.25×2 + 1.5×2 + 2×0 + 0.1×1, with Decimal's exponent: the smallest is −2 (from −0.25×2).
    assert_eq!(
        result.status, "complete_uncalibrated_score",
        "{:?}",
        result.withheld_because
    );
    assert_eq!(result.raw_score.as_deref(), Some("3.10"));

    // Rejected with a clear reason: another build, and an ID that is not a valid custom ID.
    for (name, text, expected) in [
        ("grch37", "#pgs_id=OLD\n#genome_build=GRCh37\n", "genome_build=GRCh38"),
        (
            "bad-id",
            "#pgs_id=my score\n#genome_build=GRCh38\n",
            "custom score needs #pgs_id=",
        ),
    ] {
        let file = out.join(format!("{name}.tsv"));
        std::fs::write(
            &file,
            format!("{text}chr_name\tchr_position\teffect_allele\teffect_weight\n1\t4\tC\t0.5\n"),
        )
        .unwrap();
        let err = compile_file(&file, None, &reference, &identity, None, &out)
            .expect_err("rejected")
            .to_string();
        assert!(err.contains(expected), "{name}: {err}");
    }
    std::fs::remove_dir_all(&out).unwrap();
}

/// Informational `variant_description`s are scored only with `accept_informational_descriptions`; a
/// description that changes the model stays under review either way.
#[test]
fn informational_descriptions_are_opt_in() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let out = std::env::temp_dir().join(format!("pgsum-annotated-test-{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let reference = Reference::open(&fixtures.join("synthetic.fa")).unwrap();
    let identity = reference_identity(&reference).unwrap();
    compile_file(&fixtures.join("ANNOTATED.tsv"), None, &reference, &identity, None, &out).unwrap();
    let path = out.join("ANNOTATED.pgsp");
    let (table, _) = extract(
        &fixtures.join("synthetic.g.vcf.gz"),
        &reference,
        &identity,
        std::slice::from_ref(&path),
        &opts(None, false, ScanMode::Full),
    )
    .unwrap();
    let pack = Pack::open(&path).unwrap();
    let statuses = |options: &pgsum::score::Options| {
        let mut tsv = Vec::new();
        let result = pgsum::score::score(&pack, &table, options, Some(&mut tsv)).unwrap();
        let rows: Vec<String> = String::from_utf8(tsv)
            .unwrap()
            .lines()
            .skip(1)
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                [f[13], f[18]].join("\t")
            })
            .collect();
        (result, rows)
    };
    let (default, rows) = statuses(&Default::default());
    assert_eq!(
        rows,
        [
            "model_term_requires_review\t",
            "model_term_requires_review\t",
            "scorable_observation\t"
        ]
    );
    assert_eq!(default.partial.raw_score, "0.4");
    let accept = pgsum::score::Options {
        accept_informational_descriptions: true,
        ..Default::default()
    };
    let (accepted, rows) = statuses(&accept);
    assert_eq!(
        rows,
        [
            "scorable_observation\taccepted",
            "model_term_requires_review\t",
            "scorable_observation\t"
        ]
    );
    // 0.5 × 1 (chr1:4 heterozygous) + 0.2 × 2 (chr1:5 homozygous ALT)
    assert_eq!(accepted.partial.raw_score, "0.9");
    assert_eq!(accepted.informational_descriptions.scorable_terms, 1);
    std::fs::remove_dir_all(&out).unwrap();
}

/// Palindromic SNVs are read on the forward strand only with `allow_inferred_palindromes`, and only in a
/// score whose other SNVs are (at least 99.9%, with at least 100 of them) on the forward strand.
#[test]
fn inferred_palindromes_need_a_forward_strand_score() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let out = std::env::temp_dir().join(format!("pgsum-palindrome-test-{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let reference = Reference::open(&fixtures.join("synthetic.fa")).unwrap();
    let identity = reference_identity(&reference).unwrap();
    let mut paths = Vec::new();
    for id in ["PAL_OK", "PAL_MIXED"] {
        let h = compile_file(
            &fixtures.join(format!("{id}.tsv")),
            None,
            &reference,
            &identity,
            None,
            &out,
        )
        .unwrap();
        let p = h.palindromes.expect("palindromes summarised");
        assert_eq!(p.applied, id == "PAL_OK", "{id}: {p:?}");
        paths.push(out.join(format!("{id}.pgsp")));
    }
    let (table, _) = extract(
        &fixtures.join("synthetic.g.vcf.gz"),
        &reference,
        &identity,
        &paths,
        &opts(None, false, ScanMode::Full),
    )
    .unwrap();
    let allow = pgsum::score::Options {
        allow_inferred_palindromes: true,
        ..Default::default()
    };
    let ok = Pack::open(&paths[0]).unwrap();
    let default = pgsum::score::score(&ok, &table, &Default::default(), None).unwrap();
    assert_eq!(
        (default.scorable_terms, default.partial.raw_score.as_str()),
        (100, "5.050")
    );
    // Plus chr1:1 A/T (effect is the reference, homozygous reference: 1.0 × 2) and chr1:2 G/C (effect is the ALT,
    // homozygous reference: 0.5 × 0); chr1:5 A/T is not scorable (the sample carries G).
    let inferred = pgsum::score::score(&ok, &table, &allow, None).unwrap();
    assert_eq!(
        (inferred.scorable_terms, inferred.partial.raw_score.as_str()),
        (102, "7.050")
    );
    assert_eq!(inferred.inferred_palindromes.scorable_terms, 2);
    assert_eq!(inferred.states.get("other_called_allele"), Some(&1));

    let mixed = pgsum::score::score(&Pack::open(&paths[1]).unwrap(), &table, &allow, None).unwrap();
    assert_eq!(mixed.inferred_palindromes.strand_consistent, Some(false));
    assert_eq!(mixed.inferred_palindromes.scorable_terms, 0);
    assert_eq!(mixed.states.get("unresolved_orientation"), Some(&1));
    std::fs::remove_dir_all(&out).unwrap();
}

/// Indels are scored only with `allow_inferred_indels`: oriented by reference fit or a public pair, matched to
/// gVCF records after normalization (including a record written with extra context), and called homozygous
/// reference when passing blocks cover the whole span. A both-fit indel without a public record stays
/// unresolved.
#[test]
fn inferred_indels() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let out = std::env::temp_dir().join(format!("pgsum-indel-test-{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let reference = Reference::open(&fixtures.join("synthetic.fa")).unwrap();
    let identity = reference_identity(&reference).unwrap();
    let public = pgsum::public::PublicVariants::open(&fixtures.join("public_variants.tsv")).unwrap();
    let header = compile_file(
        &fixtures.join("INDELS.tsv"),
        None,
        &reference,
        &identity,
        Some(&public),
        &out,
    )
    .unwrap();
    let s = header.sequences.expect("sequences summarised");
    assert_eq!(
        (s.terms, s.reference_fit, s.public_pair, s.unresolved_both_fit),
        (5, 1, 3, 1)
    );
    let path = out.join("INDELS.pgsp");
    let (table, _) = extract(
        &fixtures.join("synthetic_indel.g.vcf.gz"),
        &reference,
        &identity,
        std::slice::from_ref(&path),
        &opts(None, false, ScanMode::Full),
    )
    .unwrap();
    let pack = Pack::open(&path).unwrap();
    let rows = |options: &pgsum::score::Options| {
        let mut tsv = Vec::new();
        let r = pgsum::score::score(&pack, &table, options, Some(&mut tsv)).unwrap();
        let rows: Vec<String> = String::from_utf8(tsv)
            .unwrap()
            .lines()
            .skip(1)
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                [f[13], f[14], f[15], f[16], f[17]].join(" ")
            })
            .collect();
        (r, rows)
    };
    let (default, _) = rows(&Default::default());
    assert_eq!(default.scorable_terms, 0);
    let allow = pgsum::score::Options {
        allow_inferred_indels: true,
        ..Default::default()
    };
    let (result, rows) = rows(&allow);
    assert_eq!(
        rows,
        [
            "scorable_observation observed_variant 1 0.5 public_sequence_pair",
            "scorable_observation observed_variant 1 0.25 public_sequence_pair",
            "scorable_observation observed_reference 0 0 public_sequence_pair",
            "scorable_observation observed_reference 0 0.0 reference_fit_sequence",
            "unresolved_orientation    ",
        ]
    );
    assert_eq!(result.partial.raw_score, "0.75");
    assert_eq!(
        result.inferred_indels.scorable_terms.get("public_sequence_pair"),
        Some(&3)
    );
    std::fs::remove_dir_all(&out).unwrap();
}

/// The synthetic gVCF's haploid chrX call (`1` at chrX:2) is unsupported by default and read as `1/1` with
/// `haploid_xy_as_homozygous`, which also changes the policy ID in the table.
#[test]
fn haploid_sex_chromosome_calls() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let out = std::env::temp_dir().join(format!("pgsum-haploid-test-{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let reference = Reference::open(&fixtures.join("synthetic.fa")).unwrap();
    let identity = reference_identity(&reference).unwrap();
    compile_file(
        &fixtures.join("PGS999998_hmPOS_GRCh38.txt.gz"),
        Some(&fixtures.join("PGS999998.metadata.json")),
        &reference,
        &identity,
        None,
        &out,
    )
    .unwrap();
    let path = out.join("PGS999998.pgsp");
    let pack = Pack::open(&path).unwrap();
    let last_row = |haploid: bool| {
        let (table, _) = extract(
            &fixtures.join("synthetic.g.vcf.gz"),
            &reference,
            &identity,
            std::slice::from_ref(&path),
            &opts(None, haploid, ScanMode::Full),
        )
        .unwrap();
        let mut tsv = Vec::new();
        pgsum::score::score(&pack, &table, &Default::default(), Some(&mut tsv)).unwrap();
        let text = String::from_utf8(tsv).unwrap();
        let f: Vec<String> = text.lines().last().unwrap().split('\t').map(str::to_owned).collect();
        (
            table.header.policy.id.clone(),
            [f[1].clone(), f[13].clone(), f[15].clone()].join(" "),
        )
    };
    assert_eq!(
        last_row(false),
        (
            "pgsum-diploid-dp10-gq20-pass-v1".into(),
            "chrX unsupported_ploidy ".into()
        )
    );
    assert_eq!(
        last_row(true),
        (
            "pgsum-dp10-gq20-pass-haploid-xy-homozygous-v1".into(),
            "chrX scorable_observation 2".into()
        )
    );
    std::fs::remove_dir_all(&out).unwrap();
}

/// Reading only the index chunks gives the same table body as reading the whole file, with a tabix or a CSI
/// index, for SNV and indel targets; `--scan indexed` without an index is an error.
#[test]
fn indexed_reads_match_full_scan() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let out = std::env::temp_dir().join(format!("pgsum-indexed-test-{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let reference = Reference::open(&fixtures.join("synthetic.fa")).unwrap();
    let identity = reference_identity(&reference).unwrap();
    let public = pgsum::public::PublicVariants::open(&fixtures.join("public_variants.tsv")).unwrap();
    let mut packs = Vec::new();
    for id in ["PGS999997", "PGS999998"] {
        compile_file(
            &fixtures.join(format!("{id}_hmPOS_GRCh38.txt.gz")),
            Some(&fixtures.join(format!("{id}.metadata.json"))),
            &reference,
            &identity,
            None,
            &out,
        )
        .unwrap();
        packs.push(out.join(format!("{id}.pgsp")));
    }
    compile_file(
        &fixtures.join("INDELS.tsv"),
        None,
        &reference,
        &identity,
        Some(&public),
        &out,
    )
    .unwrap();
    let indel_pack = vec![out.join("INDELS.pgsp")];
    for (gvcf, packs) in [
        ("synthetic.g.vcf.gz", &packs),
        ("synthetic_indel.g.vcf.gz", &indel_pack),
    ] {
        for ext in ["tbi", "csi"] {
            let dir = out.join(format!("{gvcf}-{ext}"));
            std::fs::create_dir_all(&dir).unwrap();
            let copy = dir.join(gvcf);
            std::fs::copy(fixtures.join(gvcf), &copy).unwrap();
            let run = |mode| {
                extract(&copy, &reference, &identity, packs, &opts(None, false, mode))
                    .unwrap()
                    .0
            };
            assert!(
                extract(
                    &copy,
                    &reference,
                    &identity,
                    packs,
                    &opts(None, false, ScanMode::Indexed)
                )
                .is_err(),
                "{gvcf}: indexed without an index"
            );
            std::fs::copy(
                fixtures.join(format!("{gvcf}.{ext}")),
                dir.join(format!("{gvcf}.{ext}")),
            )
            .unwrap();
            // An index older than the file is not used.
            let index = dir.join(format!("{gvcf}.{ext}"));
            let old = std::fs::metadata(&copy).unwrap().modified().unwrap() - std::time::Duration::from_secs(3600);
            std::fs::File::options()
                .write(true)
                .open(&index)
                .unwrap()
                .set_modified(old)
                .unwrap();
            assert!(
                extract(
                    &copy,
                    &reference,
                    &identity,
                    packs,
                    &opts(None, false, ScanMode::Indexed)
                )
                .is_err()
            );
            std::fs::File::options()
                .write(true)
                .open(&index)
                .unwrap()
                .set_modified(std::time::SystemTime::now())
                .unwrap();
            let (full, indexed, auto) = (run(ScanMode::Full), run(ScanMode::Indexed), run(ScanMode::Auto));
            assert!(full.header.targets > 0);
            for other in [&indexed, &auto] {
                assert_eq!(other.header.body_sha256, full.header.body_sha256, "{gvcf} with .{ext}");
                assert_eq!(other.header.gvcf, full.header.gvcf, "{gvcf} with .{ext}");
                assert_eq!(other.header.states, full.header.states, "{gvcf} with .{ext}");
            }
        }
    }
    std::fs::remove_dir_all(&out).unwrap();
}

/// The same calls come out of the synthetic gVCF when it is uncompressed, plain gzip, has a second sample
/// (read with `sample`), or names contigs `1` and `X` against a reference that does too.
#[test]
fn input_variants_read_alike() {
    use std::io::Write;
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let out = std::env::temp_dir().join(format!("pgsum-inputs-test-{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let text = {
        let mut s = String::new();
        std::io::Read::read_to_string(
            &mut flate2::read::MultiGzDecoder::new(std::fs::File::open(fixtures.join("synthetic.g.vcf.gz")).unwrap()),
            &mut s,
        )
        .unwrap();
        s
    };
    let compile = |reference: &Reference, dir: &Path| {
        let identity = reference_identity(reference).unwrap();
        std::fs::create_dir_all(dir).unwrap();
        compile_file(
            &fixtures.join("PGS999998_hmPOS_GRCh38.txt.gz"),
            Some(&fixtures.join("PGS999998.metadata.json")),
            reference,
            &identity,
            None,
            dir,
        )
        .unwrap();
        (identity, vec![dir.join("PGS999998.pgsp")])
    };
    let reference = Reference::open(&fixtures.join("synthetic.fa")).unwrap();
    let (identity, packs) = compile(&reference, &out.join("packs"));
    let calls = |gvcf: &Path, reference: &Reference, identity, packs: &[std::path::PathBuf], sample: Option<&str>| {
        let options = Options {
            sample,
            ..opts(None, false, ScanMode::Full)
        };
        let (table, _) = extract(gvcf, reference, identity, packs, &options)?;
        let pack = Pack::open(&packs[0]).unwrap();
        let mut tsv = Vec::new();
        pgsum::score::score(&pack, &table, &Default::default(), Some(&mut tsv)).unwrap();
        Ok::<_, pgsum::Error>((table.header.states.clone(), table.header.body_sha256.clone(), tsv))
    };
    let expected = calls(
        &fixtures.join("synthetic.g.vcf.gz"),
        &reference,
        &identity,
        &packs,
        None,
    )
    .unwrap();

    let plain = out.join("synthetic.vcf");
    std::fs::write(&plain, &text).unwrap();
    assert_eq!(
        calls(&plain, &reference, &identity, &packs, None).unwrap(),
        expected,
        "uncompressed"
    );

    let gzip = out.join("synthetic.vcf.gz");
    let mut encoder = flate2::write::GzEncoder::new(std::fs::File::create(&gzip).unwrap(), Default::default());
    encoder.write_all(text.as_bytes()).unwrap();
    encoder.finish().unwrap();
    assert_eq!(
        calls(&gzip, &reference, &identity, &packs, None).unwrap(),
        expected,
        "plain gzip"
    );

    // A second sample after SYNTH; records are kept with SYNTH's column only, so even the body is unchanged.
    let two: String = text
        .lines()
        .map(|l| match l {
            _ if l.starts_with("##") => format!("{l}\n"),
            _ if l.starts_with("#CHROM") => format!("{l}\tOTHER\n"),
            _ => format!("{l}\t./.\n"),
        })
        .collect();
    let multi = out.join("two_samples.vcf");
    std::fs::write(&multi, two).unwrap();
    assert_eq!(
        calls(&multi, &reference, &identity, &packs, Some("SYNTH")).unwrap(),
        expected,
        "--sample"
    );
    assert!(
        calls(&multi, &reference, &identity, &packs, None).is_err(),
        "two samples, none chosen"
    );
    assert!(calls(&multi, &reference, &identity, &packs, Some("NOBODY")).is_err());

    // Another assembly's contig lengths, and BCF, are refused.
    let other = out.join("other_assembly.vcf");
    std::fs::write(
        &other,
        text.replace("##contig=<ID=chr1,length=60>", "##contig=<ID=chr1,length=61>"),
    )
    .unwrap();
    let err = calls(&other, &reference, &identity, &packs, None)
        .unwrap_err()
        .to_string();
    assert!(err.contains("different assembly"), "{err}");
    let bcf = out.join("synthetic.bcf");
    std::fs::write(&bcf, b"BCF\x02\x02rest").unwrap();
    let err = calls(&bcf, &reference, &identity, &packs, None)
        .unwrap_err()
        .to_string();
    assert!(err.contains("BCF is not supported"), "{err}");

    // `1` and `X` in the VCF and the reference FASTA and index.
    let unprefixed = |s: &str| s.replace("\nchr", "\n").replace("ID=chr", "ID=").replace(">chr", ">");
    assert_eq!(
        fasta_index(&std::fs::read_to_string(fixtures.join("synthetic.fa")).unwrap()),
        std::fs::read_to_string(fixtures.join("synthetic.fa.fai")).unwrap()
    );
    let fasta = out.join("renamed.fa");
    let renamed_fasta = std::fs::read_to_string(fixtures.join("synthetic.fa"))
        .unwrap()
        .replace(">chr", ">");
    std::fs::write(&fasta, &renamed_fasta).unwrap();
    std::fs::write(out.join("renamed.fa.fai"), fasta_index(&renamed_fasta)).unwrap();
    let renamed_reference = Reference::open(&fasta).unwrap();
    let (renamed_identity, renamed_packs) = compile(&renamed_reference, &out.join("renamed-packs"));
    let renamed = out.join("renamed.vcf");
    std::fs::write(&renamed, unprefixed(&text)).unwrap();
    let (states, _, tsv) = calls(&renamed, &renamed_reference, &renamed_identity, &renamed_packs, None).unwrap();
    assert_eq!((states, tsv), (expected.0, expected.2), "contigs named 1 and X");
    std::fs::remove_dir_all(&out).unwrap();
}

/// A `.fai` for FASTA text whose sequence lines all have the same width.
fn fasta_index(fasta: &str) -> String {
    let mut out = String::new();
    let mut offset = 0usize;
    let mut current: Option<(String, usize, usize, usize, usize)> = None; // name, length, start, bases, width
    let finish = |c: Option<(String, usize, usize, usize, usize)>, out: &mut String| {
        if let Some((name, length, start, bases, width)) = c {
            out.push_str(&format!("{name}\t{length}\t{start}\t{bases}\t{width}\n"));
        }
    };
    for line in fasta.split_inclusive('\n') {
        if let Some(name) = line.strip_prefix('>') {
            finish(current.take(), &mut out);
            let name = name.split_whitespace().next().unwrap().to_owned();
            current = Some((name, 0, offset + line.len(), 0, 0));
        } else if let Some(c) = current.as_mut() {
            let bases = line.trim_end().len();
            if c.3 == 0 {
                (c.3, c.4) = (bases, line.len());
            }
            c.1 += bases;
        }
        offset += line.len();
    }
    finish(current, &mut out);
    out
}
