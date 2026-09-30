//! Cohort reasons, numeric reports and source references must match independently extracted samples.
use pgsum::cohort::{CohortOptions, CohortTable};
use pgsum::dosage::Field;
use pgsum::genotype::State;
use pgsum::genotypes::target_key;
use pgsum::pack::Pack;
use serde_json::Value;
use std::path::{Path, PathBuf};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures").join(name)
}
fn root(name: &str) -> PathBuf {
    let root = std::env::temp_dir().join(format!("pgsum-cohort-diagnostics-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    root
}
fn compile(root: &Path) -> (pgsum::reference::Reference, pgsum::pack::ReferenceIdentity, PathBuf) {
    let reference = pgsum::reference::Reference::open(&fixture("synthetic.fa")).unwrap();
    let identity = pgsum::compile::reference_identity(&reference).unwrap();
    let scoring = root.join("DIAGNOSTICS.tsv");
    std::fs::write(&scoring, "#pgs_id=DIAGNOSTICS\n#genome_build=GRCh38\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\tis_dominant\n1\t2\tT\tC\t0.2\tFalse\n1\t3\tA\tG\t-2\tTrue\n1\t4\tC\tT\t1\tFalse\n").unwrap();
    pgsum::compile::compile_file(&scoring, None, &reference, &identity, None, root).unwrap();
    (reference, identity, root.join("DIAGNOSTICS.pgsp"))
}

#[test]
fn all_fields_preserve_reasons_and_match_single_sample_reports() {
    let root = root("parity");
    let (reference, identity, path) = compile(&root);
    let pack = Pack::open(&path).unwrap();
    let input = root.join("cohort.vcf");
    std::fs::write(&input, concat!(
        "##fileformat=VCFv4.3\n##contig=<ID=1,length=60>\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\n",
        // An ignored contig is still a record in the source ordinal sequence.
        "unplaced\t1\tignored\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\n",
        "1\t2\tfirst\tC\tT\t.\tPASS\t.\tGT:DS:GP:DP:GQ\t0/1:0.4:0.25,0.5,0.25:30:60\t./.:.:.:30:60\t0/1:0.4:0.25,0.5,0.25:5:60\t1:1:0,1,0:30:60\n",
        "1\t3\tsecond\tG\tA\t.\tPASS\t.\tGT:DS:GP:DP:GQ:FT\t0/1:0.4:0.25,0.5,0.25:30:60:PASS\t1/1:2:0,0,1:30:60:FAIL\t0/1:0.4:0.25,0.5,0.25:30:60:PASS\t0/1:0.4:0.25,0.5,0.25:30:60:PASS\n"
    )).unwrap();
    for field in [Field::Gt, Field::Ds, Field::Gp] {
        let mut prior_report = None;
        for threads in [1, 4] {
            let cohort_path = root.join(format!("{}-{threads}.pgsc", field.label()));
            pgsum::cohort::extract_cohort(
                &input,
                &reference,
                &identity,
                std::slice::from_ref(&path),
                &CohortOptions {
                    dosage_field: field,
                    threads,
                    ..Default::default()
                },
                &cohort_path,
            )
            .unwrap();
            let cohort = CohortTable::open(&cohort_path).unwrap();
            assert_eq!(cohort.header.schema, pgsum::cohort::DIAGNOSTIC_SCHEMA);
            cohort.validate_diagnostics().unwrap();
            let key = target_key(1, 2, b'C', b'T');
            assert_eq!(cohort.source_records(key).unwrap(), Some(vec![2]));
            assert_eq!(cohort.call_state(key, 2).unwrap(), Some(State::LowQuality));
            assert_eq!(cohort.call_state(key, 3).unwrap(), Some(State::UnsupportedPloidy));
            assert_eq!(
                cohort.source_records(target_key(1, 4, b'T', b'C')).unwrap(),
                Some(vec![])
            );
            let scores = pgsum::cohort::score_cohort(&pack, &cohort, &Default::default()).unwrap();
            let mut output = Vec::new();
            pgsum::missingness::write_cohort(&pack, &cohort, &Default::default(), &scores, &mut output).unwrap();
            if let Some(prior) = &prior_report {
                assert_eq!(&output, prior);
            }
            prior_report = Some(output.clone());
            let reports: Vec<Value> = String::from_utf8(output)
                .unwrap()
                .lines()
                .map(|l| serde_json::from_str(l).unwrap())
                .collect();
            for (sample, name) in ["A", "B", "C", "D"].iter().enumerate() {
                let (table, _) = pgsum::extract::extract(
                    &input,
                    &reference,
                    &identity,
                    std::slice::from_ref(&path),
                    &pgsum::extract::Options {
                        dosage_field: field,
                        sample: Some(name),
                        threads,
                        ..Default::default()
                    },
                )
                .unwrap();
                let single = pgsum::score::score(&pack, &table, &Default::default(), None).unwrap();
                assert_eq!(reports[sample]["partial_raw_score"], single.partial.raw_score);
                let single_report = serde_json::to_value(single.missingness.unwrap()).unwrap();
                let actual = &reports[sample]["missingness"];
                for field in [
                    "missing_terms",
                    "bounded_terms",
                    "unbounded_terms",
                    "by_status",
                    "bounded_missing_contribution",
                    "completion_score_bounds",
                    "completion_bounds_withheld_because",
                ] {
                    assert_eq!(actual[field], single_report[field], "{name} / {field}");
                }
                for field in ["ranked_missing_terms", "unbounded_examples"] {
                    let mut a = actual[field].clone();
                    let mut b = single_report[field].clone();
                    for row in a.as_array_mut().unwrap().iter_mut().chain(b.as_array_mut().unwrap()) {
                        row.as_object_mut().unwrap().remove("evidence");
                    }
                    assert_eq!(a, b, "{name} / {field}");
                }
                for entry in table.entries().unwrap() {
                    assert_eq!(cohort.call_state(entry.key, sample).unwrap(), Some(entry.state));
                }
            }
            // Exercise the public CLI and compressed report publication once.
            if field == Field::Gt && threads == 1 {
                let out = root.join("reports");
                let result = std::process::Command::new(env!("CARGO_BIN_EXE_pgsum"))
                    .args([
                        "score",
                        "--genotypes",
                        cohort_path.to_str().unwrap(),
                        "--pack",
                        path.to_str().unwrap(),
                        "--out",
                        out.to_str().unwrap(),
                        "--missingness",
                        "--threads",
                        "2",
                    ])
                    .output()
                    .unwrap();
                assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
                let text =
                    zstd::decode_all(std::fs::File::open(out.join("cohort-missingness.jsonl.zst")).unwrap()).unwrap();
                assert_eq!(String::from_utf8(text).unwrap().lines().count(), 4);
                let metadata: Value =
                    serde_json::from_slice(&std::fs::read(out.join("cohort-scores.metadata.json")).unwrap()).unwrap();
                assert_eq!(
                    metadata["missingness_report"]["sha256"],
                    pgsum::digest::file_sha256(&out.join("cohort-missingness.jsonl.zst")).unwrap()
                );
                let rerun = std::process::Command::new(env!("CARGO_BIN_EXE_pgsum"))
                    .arg("score")
                    .arg("--genotypes")
                    .arg(&cohort_path)
                    .arg("--pack")
                    .arg(&path)
                    .arg("--out")
                    .arg(&out)
                    .output()
                    .unwrap();
                assert!(rerun.status.success());
                let metadata: Value =
                    serde_json::from_slice(&std::fs::read(out.join("cohort-scores.metadata.json")).unwrap()).unwrap();
                assert!(metadata["missingness_report"].is_null());
                let duplicate = root.join("DUPLICATE.pgsp");
                std::fs::copy(&path, &duplicate).unwrap();
                let rejected = std::process::Command::new(env!("CARGO_BIN_EXE_pgsum"))
                    .arg("score")
                    .arg("--genotypes")
                    .arg(&cohort_path)
                    .arg("--pack")
                    .arg(&path)
                    .arg("--pack")
                    .arg(&duplicate)
                    .arg("--out")
                    .arg(root.join("duplicate-report"))
                    .arg("--missingness")
                    .output()
                    .unwrap();
                assert!(!rejected.status.success());
                assert!(String::from_utf8_lossy(&rejected.stderr).contains("one pack per score ID"));
                std::fs::remove_file(duplicate).unwrap();
            }
        }
    }
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn legacy_cohorts_keep_scores_without_inventing_missing_reasons() {
    let root = root("legacy");
    let reference = pgsum::reference::Reference::open(&fixture("synthetic.fa")).unwrap();
    let identity = pgsum::compile::reference_identity(&reference).unwrap();
    pgsum::compile::compile_file(
        &fixture("PGS999997_hmPOS_GRCh38.txt.gz"),
        Some(&fixture("PGS999997.metadata.json")),
        &reference,
        &identity,
        None,
        &root,
    )
    .unwrap();
    let pack = Pack::open(&root.join("PGS999997.pgsp")).unwrap();
    for version in [1, 2] {
        let table = CohortTable::open(&fixture(&format!("legacy-cohort-v{version}.pgsc"))).unwrap();
        let key = target_key(1, 2, b'C', b'T');
        assert_eq!(table.call_state(key, 0).unwrap(), Some(State::ObservedVariant));
        assert_eq!(table.call_state(key, 1).unwrap(), None);
        assert_eq!(table.source_records(key).unwrap(), None);
        let scores = pgsum::cohort::score_cohort(&pack, &table, &Default::default()).unwrap();
        let mut out = Vec::new();
        pgsum::missingness::write_cohort(&pack, &table, &Default::default(), &scores, &mut out).unwrap();
        let rows: Vec<Value> = String::from_utf8(out)
            .unwrap()
            .lines()
            .map(|l| serde_json::from_str(l).unwrap())
            .collect();
        assert_eq!(rows[1]["missingness"]["by_status"]["missing_reason_not_retained"], 10);
        assert!(rows[1]["missingness"]["completion_score_bounds"].is_null());
    }
    std::fs::remove_dir_all(root).unwrap();
}

// Mutate the independently compressed diagnostic frame while keeping the score calls intact.
fn alter_diagnostics(source: &Path, dest: &Path, change: impl FnOnce(&mut Value), bind_digest: bool) {
    use sha2::{Digest, Sha256};
    let bytes = std::fs::read(source).unwrap();
    let end = bytes.len();
    let len = u64::from_le_bytes(bytes[end - 16..end - 8].try_into().unwrap()) as usize;
    let start = end - 16 - len;
    let mut header: pgsum::cohort::Header = serde_json::from_slice(&bytes[start..end - 16]).unwrap();
    let frame = &mut header.diagnostics[0];
    let first = 8 + frame.offset as usize;
    let decoded = zstd::decode_all(&bytes[first..first + frame.compressed_bytes as usize]).unwrap();
    let mut rows: Value = serde_json::from_slice(&decoded).unwrap();
    change(&mut rows);
    let decoded = serde_json::to_vec(&rows).unwrap();
    let compressed = zstd::bulk::compress(&decoded, 3).unwrap();
    let mut output = bytes[..start].to_vec();
    frame.offset = (output.len() - 8) as u64;
    frame.compressed_bytes = compressed.len() as u64;
    if bind_digest {
        frame.sha256 = format!(
            "sha256:{}",
            Sha256::digest(&decoded)
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect::<String>()
        );
        let mut digest = Sha256::new();
        for value in std::iter::once(&header.sequence_frame.sha256)
            .chain(header.blocks.iter().map(|b| &b.frame.sha256))
            .chain(header.diagnostics.iter().map(|d| &d.sha256))
        {
            digest.update(value.as_bytes());
            digest.update(b"\n");
        }
        header.body_sha256 = format!(
            "sha256:{}",
            digest.finalize().iter().map(|b| format!("{b:02x}")).collect::<String>()
        );
    }
    output.extend(compressed);
    let json = serde_json::to_vec(&header).unwrap();
    output.extend(&json);
    output.extend((json.len() as u64).to_le_bytes());
    output.extend(b"PGSUMCOE");
    std::fs::write(dest, output).unwrap();
}

#[test]
fn corrupt_or_inconsistent_diagnostics_are_rejected_when_read() {
    let root = root("corrupt");
    let (reference, identity, path) = compile(&root);
    let pack = Pack::open(&path).unwrap();
    let vcf = root.join("source.vcf");
    std::fs::write(&vcf, "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n1\t2\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\n").unwrap();
    let source = root.join("source.pgsc");
    let options = CohortOptions {
        accept_missing_quality: true,
        ..Default::default()
    };
    pgsum::cohort::extract_cohort(
        &vcf,
        &reference,
        &identity,
        std::slice::from_ref(&path),
        &options,
        &source,
    )
    .unwrap();
    let original = CohortTable::open(&source).unwrap();
    let expected = pgsum::cohort::score_cohort(&pack, &original, &Default::default()).unwrap();
    for case in 0..5 {
        let dest = root.join(format!("bad-{case}.pgsc"));
        alter_diagnostics(
            &source,
            &dest,
            |rows| {
                match case {
                    0 | 1 => rows[0]["source_records"] = serde_json::json!([2]), // beyond source record count
                    2 => rows[0]["missing"] = serde_json::json!({"Uniform": 200}), // unknown state
                    3 => rows[0]["missing"] = serde_json::json!({"Uniform": State::NoCall as u8}), // passing packed call
                    _ => {
                        rows.as_array_mut().unwrap().pop();
                    } // incorrect row count
                }
            },
            case != 0,
        );
        let table = CohortTable::open(&dest).unwrap();
        // The numeric scoring path reads no diagnostic bytes.
        let scored = pgsum::cohort::score_cohort(&pack, &table, &Default::default()).unwrap();
        assert_eq!(scored.text(0), expected.text(0));
        assert!(table.validate_diagnostics().is_err(), "case {case}");
        assert!(
            pgsum::missingness::write_cohort(&pack, &table, &Default::default(), &scored, &mut Vec::new()).is_err()
        );
    }
    // The cohort path must apply the same reference-length safeguard as selected samples.
    let wrong = std::fs::read_to_string(&vcf).unwrap().replace(
        "##fileformat=VCFv4.3",
        "##fileformat=VCFv4.3\n##contig=<ID=1,length=61>",
    );
    std::fs::write(&vcf, wrong).unwrap();
    assert!(
        pgsum::cohort::extract_cohort(&vcf, &reference, &identity, &[path], &options, &root.join("wrong.pgsc"))
            .is_err()
    );
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn source_ordinals_count_ignored_records_across_scan_chunks() {
    let root = root("chunks");
    let (reference, identity, path) = compile(&root);
    let vcf = root.join("source.vcf");
    let mut text = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n1\t2\t.\tC\tT\t.\tPASS\t.\tGT\t./.\n".to_string();
    let ignored = format!("unplaced\t1\t{}\tA\tG\t.\tPASS\t.\tGT\t0/0\n", "x".repeat(1500));
    for _ in 0..3000 {
        text.push_str(&ignored);
    }
    text.push_str("1\t3\t.\tG\tA\t.\tPASS\t.\tGT\t./.\n");
    std::fs::write(&vcf, text).unwrap();
    for threads in [1, 4] {
        let dest = root.join(format!("source-{threads}.pgsc"));
        pgsum::cohort::extract_cohort(
            &vcf,
            &reference,
            &identity,
            std::slice::from_ref(&path),
            &CohortOptions {
                threads,
                accept_missing_quality: true,
                ..Default::default()
            },
            &dest,
        )
        .unwrap();
        let table = CohortTable::open(&dest).unwrap();
        assert_eq!(table.header.records_scanned, 3002);
        assert_eq!(
            table.source_records(target_key(1, 2, b'C', b'T')).unwrap(),
            Some(vec![1])
        );
        assert_eq!(
            table.source_records(target_key(1, 3, b'G', b'A')).unwrap(),
            Some(vec![3002])
        );
    }
    std::fs::remove_dir_all(root).unwrap();
}
