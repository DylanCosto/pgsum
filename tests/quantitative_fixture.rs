//! End-to-end expectations are hand calculated; no scorer is used to construct expected scores.
use std::path::{Path, PathBuf};
use std::process::Command;

fn run(args: &[&str]) -> std::process::Output {
    let output = Command::new(env!("CARGO_BIN_EXE_pgsum")).args(args).output().unwrap();
    assert!(
        output.status.success(),
        "{:?}: {}",
        args,
        String::from_utf8_lossy(&output.stderr)
    );
    output
}
fn text(p: &Path) -> &str {
    p.to_str().unwrap()
}
fn fixture(name: &str) -> (PathBuf, PathBuf, PathBuf) {
    let dir = std::env::temp_dir().join(format!("pgsum-quantitative-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/synthetic.fa");
    let scoring = dir.join("QUANT.tsv");
    std::fs::write(&scoring, "#pgs_id=QUANT\n#genome_build=GRCh38\n#variants_number=2\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n1\t2\tT\tC\t0.2\n1\t3\tG\tT\t-0.5\n").unwrap();
    run(&[
        "compile",
        text(&scoring),
        "--reference",
        text(&reference),
        "--out",
        text(&dir),
    ]);
    let vcf = dir.join("sample.vcf");
    std::fs::write(&vcf, "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DS:GP\t1/1:0.125:0.2,0.3,0.5\n1\t3\t.\tG\tT\t.\tPASS\t.\tGT:DS:GP\t0/0:1.25:0.5,0.3,0.2\n").unwrap();
    (dir, reference, vcf)
}
#[test]
fn explicit_fields_survive_tables_and_use_correct_alleles() {
    let (dir, reference, vcf) = fixture("fields");
    for (field, expected, n) in [("gt", "-0.6", 0), ("ds", "-0.3500", 2), ("gp", "-0.39", 2)] {
        let out = dir.join(field);
        run(&[
            "run",
            "--gvcf",
            text(&vcf),
            "--reference",
            text(&reference),
            "--pack",
            text(&dir),
            "--out",
            text(&out),
            "--dosage-field",
            field,
            "--accept-missing-quality",
            "--terms",
        ]);
        let result: serde_json::Value =
            serde_json::from_slice(&std::fs::read(out.join("QUANT.score.json")).unwrap()).unwrap();
        assert_eq!(result["raw_score"], expected, "{field}: {result}");
        assert_eq!(result["expected_genotype_terms"], n);
        assert_eq!(result["imputation"]["observed_scores"], n > 0);
        let table = pgsum::genotypes::GenotypeTable::open(&out.join("genotypes.pgsg")).unwrap();
        assert_eq!(table.header.schema, "pgsum-genotypes-v5");
        table.preload().unwrap();
        let rewrite = out.join("rewrite.pgsg");
        table.write(&rewrite).unwrap();
        let rewritten = pgsum::genotypes::GenotypeTable::open(&rewrite).unwrap();
        assert_eq!(
            table.entries().unwrap().collect::<Vec<_>>(),
            rewritten.entries().unwrap().collect::<Vec<_>>()
        );
        let inspect = run(&["inspect", text(&rewrite)]);
        if field == "ds" {
            assert!(String::from_utf8(inspect.stdout).unwrap().contains("0.125"));
        }
    }
}
#[test]
fn ds_only_inputs_do_not_need_gt_and_never_fall_back_to_gt() {
    let (dir, reference, vcf) = fixture("missing");
    std::fs::write(&vcf, "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n1\t2\t.\tC\tT\t.\tPASS\t.\tDS\t0.125\n1\t3\t.\tG\tT\t.\tPASS\t.\tGT:DS\t0/0:.\n").unwrap();
    let out = dir.join("results");
    run(&[
        "run",
        "--gvcf",
        text(&vcf),
        "--reference",
        text(&reference),
        "--pack",
        text(&dir),
        "--out",
        text(&out),
        "--dosage-field",
        "ds",
        "--accept-missing-quality",
    ]);
    let result: serde_json::Value =
        serde_json::from_slice(&std::fs::read(out.join("QUANT.score.json")).unwrap()).unwrap();
    assert!(result["raw_score"].is_null());
    assert_eq!(result["partial"]["raw_score"], "0.0250");
    assert_eq!(result["states"]["dosage_missing"], 1);
    assert_eq!(result["expected_genotype_terms"], 1);
}

#[test]
fn filters_invalid_values_and_reference_blocks_keep_distinct_states() {
    let target = pgsum::genotype::Target {
        pos: 2,
        ref_allele: "C",
        alt: Some("T"),
        sequence: false,
        sex_chromosome: false,
    };
    let policy = pgsum::genotype::Policy {
        dosage_field: pgsum::dosage::Field::Ds,
        accept_missing_quality: true,
        ..Default::default()
    };
    let cases = [
        ("C\tT\t.\tPASS\t.\tDS\t2.00000000000000000001", "invalid_dosage"),
        ("C\tT\t.\tLowQual\t.\tDS\t0.4", "filtered"),
        ("C\tT\t.\tPASS\t.\tDS:FT\t0.4:FAIL", "genotype_filtered"),
        ("C\tT\t.\tPASS\t.\tDS:DP:GQ\t0.4:5:30", "low_quality"),
        ("C\tT\t.\tPASS\t.\tDS:DP\t0.4:30", "quality_missing"),
        ("C\tT,G\t.\tPASS\t.\tDS\t0.4,0.6", "unsupported_allele_representation"),
        ("C\tT\t.\tPASS\t.\tGT:DS\t1:0.4", "unsupported_ploidy"),
        (
            "C\t<NON_REF>\t.\tPASS\tEND=4\tGT:MIN_DP:GQ\t0/0:30:30",
            "observed_reference",
        ),
        ("C\t<NON_REF>\t.\tPASS\tEND=4\tGT:MIN_DP:GQ\t0/0:5:30", "low_quality"),
    ];
    for (record, expected) in cases {
        let r = pgsum::gvcf::parse_record(&format!("1\t2\t.\t{record}")).unwrap();
        let call = pgsum::genotype::assess(&[r], &target, &policy, |_, _| Some("C".into()));
        assert_eq!(call.state.as_str(), expected, "{record}");
        if expected == "observed_reference" {
            assert_eq!(call.alt_dosage, Some(0));
            assert!(call.measurement.is_none());
        }
    }
}

#[test]
fn gp_scores_nonadditive_models_and_ds_withholds_them() {
    let (dir, reference, vcf) = fixture("models");
    let scoring = dir.join("QUANT.tsv");
    std::fs::write(&scoring, "#pgs_id=QUANT\n#genome_build=GRCh38\n#variants_number=2\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\tis_dominant\tis_recessive\n1\t2\tT\tC\t2\tTrue\tFalse\n1\t3\tG\tT\t-0.5\tFalse\tTrue\n").unwrap();
    run(&[
        "compile",
        text(&scoring),
        "--reference",
        text(&reference),
        "--out",
        text(&dir),
    ]);
    for (field, score, status) in [
        ("gp", "1.35", "complete_uncalibrated_score"),
        ("ds", "0", "score_withheld"),
    ] {
        let out = dir.join(field);
        run(&[
            "run",
            "--gvcf",
            text(&vcf),
            "--reference",
            text(&reference),
            "--pack",
            text(&dir),
            "--out",
            text(&out),
            "--dosage-field",
            field,
            "--accept-missing-quality",
        ]);
        let result: serde_json::Value =
            serde_json::from_slice(&std::fs::read(out.join("QUANT.score.json")).unwrap()).unwrap();
        assert_eq!(result["status"], status, "{result}");
        assert_eq!(result["partial"]["raw_score"], score, "{result}");
        if field == "ds" {
            assert_eq!(result["states"]["genotype_probabilities_required"], 2);
        }
    }
}

#[test]
fn legacy_v4_files_read_and_upgrade_without_changing_calls() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/legacy-v4.pgsg");
    let old = pgsum::genotypes::GenotypeTable::open(&fixture).unwrap();
    assert_eq!(old.header.schema, "pgsum-genotypes-v4");
    assert_eq!(old.header.policy.dosage_field, pgsum::dosage::Field::Gt);
    let upgraded = std::env::temp_dir().join(format!("pgsum-v4-upgrade-{}.pgsg", std::process::id()));
    old.write(&upgraded).unwrap();
    let new = pgsum::genotypes::GenotypeTable::open(&upgraded).unwrap();
    assert_eq!(new.header.schema, "pgsum-genotypes-v5");
    assert_eq!(
        old.entries().unwrap().collect::<Vec<_>>(),
        new.entries().unwrap().collect::<Vec<_>>()
    );
    for e in old.entries().unwrap() {
        assert_eq!(old.records(e.key).unwrap(), new.records(e.key).unwrap());
    }
    std::fs::remove_file(upgraded).unwrap();
}

#[test]
fn quantitative_cohorts_match_selected_samples_across_threads() {
    let (dir, reference_path, vcf) = fixture("cohort");
    let sample_names: Vec<_> = (0..9).map(|i| format!("S{i}")).collect();
    let fields = [
        "0/0:0.125:0.2,0.3,0.5",
        "0/1:0:1,0,0",
        "1/1:2:0,0,1",
        "./.:.:.",
        "0/0:0.5:0.6,0.3,0.1",
        "0/1:1.5:0.1,0.3,0.6",
        "1/1:0.00001:0.99999,0.00001,0",
        "0/0:2.1:0.5,0.5,0.1",
        "0/1:1:0,1,0",
    ];
    let mut text_vcf = format!(
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
        sample_names.join("\t")
    );
    for (pos, r, a) in [(2, "C", "T"), (3, "G", "T")] {
        text_vcf.push_str(&format!(
            "1\t{pos}\t.\t{r}\t{a}\t.\tPASS\t.\tGT:DS:GP\t{}\n",
            fields.join("\t")
        ));
    }
    std::fs::write(&vcf, text_vcf).unwrap();
    let reference = pgsum::reference::Reference::open(&reference_path).unwrap();
    let identity = pgsum::compile::reference_identity(&reference).unwrap();
    let pack_path = dir.join("QUANT.pgsp");
    let packs = [pack_path.clone()];
    let pack = pgsum::pack::Pack::open(&pack_path).unwrap();
    for field in [pgsum::dosage::Field::Ds, pgsum::dosage::Field::Gp] {
        for threads in [1, 4] {
            let cohort_path = dir.join(format!("{}-{threads}.pgsc", field.label()));
            let options = pgsum::cohort::CohortOptions {
                dosage_field: field,
                accept_missing_quality: true,
                threads,
                ..Default::default()
            };
            pgsum::cohort::extract_cohort(&vcf, &reference, &identity, &packs, &options, &cohort_path).unwrap();
            let cohort = pgsum::cohort::CohortTable::open(&cohort_path).unwrap();
            assert_eq!(cohort.header.schema, "pgsum-cohort-v3");
            let scored = pgsum::cohort::score_cohort(&pack, &cohort, &Default::default()).unwrap();
            for (i, name) in sample_names.iter().enumerate() {
                let options = pgsum::extract::Options {
                    dosage_field: field,
                    sample: Some(name),
                    accept_missing_quality: true,
                    threads,
                    ..Default::default()
                };
                let (single, _) = pgsum::extract::extract(&vcf, &reference, &identity, &packs, &options).unwrap();
                let expected = pgsum::score::score(&pack, &single, &Default::default(), None).unwrap();
                assert_eq!(
                    scored.text(i),
                    expected.partial.raw_score,
                    "{} / {threads} / {name}",
                    field.label()
                );
                assert_eq!(scored.scorable[i], expected.scorable_terms);
                assert_eq!(scored.expected_genotype_terms[i], expected.expected_genotype_terms);
            }
            assert_eq!(scored.scorable[3], 0);
            assert_eq!(scored.scorable[7], 0);
            assert_eq!(scored.text(1), "-1.0");
            assert_eq!(scored.text(2), "0.4");
            if threads == 4 {
                let cli_table = dir.join(format!("{}-cli.pgsc", field.label()));
                run(&[
                    "extract",
                    "--gvcf",
                    text(&vcf),
                    "--reference",
                    text(&reference_path),
                    "--pack",
                    text(&pack_path),
                    "--out",
                    text(&cli_table),
                    "--all-samples",
                    "--dosage-field",
                    &field.label().to_lowercase(),
                    "--accept-missing-quality",
                    "--threads",
                    "4",
                ]);
                let output = dir.join(format!("{}-scores", field.label()));
                run(&[
                    "score",
                    "--genotypes",
                    text(&cli_table),
                    "--pack",
                    text(&pack_path),
                    "--out",
                    text(&output),
                ]);
                let metadata: serde_json::Value =
                    serde_json::from_slice(&std::fs::read(output.join("cohort-scores.metadata.json")).unwrap())
                        .unwrap();
                assert_eq!(metadata["dosage_field"], field.label());
                assert_eq!(metadata["scored_pgs_ids"], serde_json::json!(["QUANT"]));
                let tsv = String::from_utf8(
                    zstd::decode_all(std::fs::File::open(output.join("cohort-scores.tsv.zst")).unwrap()).unwrap(),
                )
                .unwrap();
                assert!(
                    tsv.lines()
                        .next()
                        .unwrap()
                        .ends_with("dosage_field\texpected_genotype_terms")
                );
                let rows: Vec<_> = tsv.lines().skip(1).collect();
                assert_eq!(rows.len(), 9);
                assert!(rows[3].ends_with(&format!("{}\t0", field.label())));
                assert!(rows[0].ends_with(&format!("{}\t2", field.label())));
            }
        }
    }
}
