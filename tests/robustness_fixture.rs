//! Hand-calculated reconstruction experiments exercised through the public CLI.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};

use pgsum::decimal::Decimal;
use serde_json::{Value, json};

const SITES: [(u32, &str, &str, &str, &str, &str); 4] = [
    // Position, REF, ALT, effect allele, other allele, weight.
    (2, "C", "T", "T", "C", "0.2"),
    (3, "G", "T", "G", "T", "-0.5"),
    (4, "T", "C", "C", "T", "1.25"),
    (5, "T", "C", "T", "C", "2"),
];
static SERIAL: AtomicUsize = AtomicUsize::new(0);

type Row = BTreeMap<String, String>;

struct Fixture {
    dir: PathBuf,
    reference: PathBuf,
    terms: usize,
}

impl Fixture {
    fn new(name: &str, terms: usize) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "pgsum-robustness-{}-{}-{name}",
            std::process::id(),
            SERIAL.fetch_add(1, AtomicOrdering::Relaxed)
        ));
        std::fs::create_dir(&dir).unwrap();
        let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/synthetic.fa");
        let mut scoring = format!(
            "#pgs_id=ROBUST\n#genome_build=GRCh38\n#variants_number={terms}\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n"
        );
        for &(pos, _, _, effect, other, weight) in &SITES[..terms] {
            scoring.push_str(&format!("1\t{pos}\t{effect}\t{other}\t{weight}\n"));
        }
        let source = dir.join("ROBUST.tsv");
        std::fs::write(&source, scoring).unwrap();
        run(&[
            "compile",
            text(&source),
            "--reference",
            text(&reference),
            "--out",
            text(&dir),
        ]);
        Self { dir, reference, terms }
    }

    fn path(&self, name: &str) -> PathBuf {
        self.dir.join(name)
    }

    fn table(&self, name: &str, samples: &[&str], format: &str, fields: &[Vec<&str>], dosage: &str) -> PathBuf {
        assert_eq!(fields.len(), self.terms);
        let mut vcf = format!(
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
            samples.join("\t")
        );
        for (&(pos, reference, alt, _, _, _), values) in SITES.iter().zip(fields) {
            assert_eq!(values.len(), samples.len());
            vcf.push_str(&format!(
                "1\t{pos}\t.\t{reference}\t{alt}\t.\tPASS\t.\t{format}\t{}\n",
                values.join("\t")
            ));
        }
        let source = self.path(&format!("{name}.vcf"));
        let out = self.path(&format!("{name}.pgsc"));
        std::fs::write(&source, vcf).unwrap();
        run(&[
            "extract",
            "--gvcf",
            text(&source),
            "--reference",
            text(&self.reference),
            "--pack",
            text(&self.path("ROBUST.pgsp")),
            "--out",
            text(&out),
            "--all-samples",
            "--dosage-field",
            dosage,
            "--accept-missing-quality",
            "--threads",
            "2",
        ]);
        out
    }

    fn standard_table(&self) -> PathBuf {
        self.table("baseline", &["A", "B", "C"], "GT", &gt_fields(), "gt")
    }

    fn manifest(&self, genotypes: &Path, scenarios: Value) -> Value {
        json!({
            "schema": "pgsum-robustness-experiment-v1",
            "genotypes": genotypes,
            "packs": ["ROBUST.pgsp"],
            "cutoffs": [50.0, 90.0],
            "top_terms": 1,
            "scenarios": scenarios
        })
    }

    fn experiment(&self, name: &str, manifest: Value) -> (PathBuf, Value) {
        let source = self.path(&format!("{name}.json"));
        let out = self.path(name);
        std::fs::write(&source, serde_json::to_vec_pretty(&manifest).unwrap()).unwrap();
        run(&[
            "robustness",
            "--experiment",
            text(&source),
            "--out",
            text(&out),
            "--terms",
        ]);
        let report = serde_json::from_slice(&std::fs::read(out.join("report.json")).unwrap()).unwrap();
        (out, report)
    }

    fn bed(&self, name: &str, contents: &str) -> PathBuf {
        let path = self.path(name);
        std::fs::write(&path, contents).unwrap();
        path
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        if !std::thread::panicking() {
            std::fs::remove_dir_all(&self.dir).unwrap();
        }
    }
}

fn text(path: &Path) -> &str {
    path.to_str().unwrap()
}

fn run(args: &[&str]) -> Output {
    let output = Command::new(env!("CARGO_BIN_EXE_pgsum")).args(args).output().unwrap();
    assert!(
        output.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn gt_fields() -> Vec<Vec<&'static str>> {
    vec![
        vec!["0/0", "0/1", "1/1"],
        vec!["0/0", "0/1", "1/1"],
        vec!["0/1", "0/0", "1/1"],
        vec!["1/1", "0/1", "0/0"],
    ]
}

#[test]
fn preflight_reports_dimensions_and_rejects_invalid_inputs_before_scoring() {
    let fixture = Fixture::new("preflight", 4);
    let table = fixture.standard_table();
    let bed = fixture.bed("coverage.bed", "chr1\t0\t8\n");
    let manifest_path = fixture.path("experiment.json");
    let mut manifest = fixture.manifest(
        &table,
        json!([
            {"id":"control"}, {"id":"covered", "masks":[{"kind":"bed","mode":"keep","path":bed}]}
        ]),
    );
    manifest["reference_sets"] = json!([{"id":"all"},{"id":"one","samples":["A"]}]);
    manifest["regions"] = json!([{"id":"region","path":bed}]);
    let invoke = |m: &Value, name: &str, preflight: bool| {
        std::fs::write(&manifest_path, serde_json::to_vec(m).unwrap()).unwrap();
        let mut command = Command::new(env!("CARGO_BIN_EXE_pgsum"));
        command.args([
            "robustness",
            "--experiment",
            text(&manifest_path),
            "--out",
            text(&fixture.path(name)),
            "--progress",
            "json",
            "--terms",
        ]);
        if preflight {
            command.arg("--preflight");
        }
        command.output().unwrap()
    };
    let result = invoke(&manifest, "checked", true);
    assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
    let events: Vec<Value> = String::from_utf8(result.stderr)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(events.last().unwrap()["phase"], "preflight_complete");
    assert!(
        events
            .iter()
            .all(|e| e["phase"] != "scoring" && e["phase"] != "complete" && e["score"].is_null())
    );
    let out = fixture.path("checked");
    assert!(!out.join("report.json").exists() && !out.join("samples.tsv").exists());
    let report: Value = serde_json::from_slice(&std::fs::read(out.join("preflight.json")).unwrap()).unwrap();
    assert_eq!(report["schema"], "pgsum-robustness-preflight-v1");
    assert_eq!(report["scoring_performed"], false);
    for (key, expected) in [
        ("model_terms", 4),
        ("candidate_term_sample_scenario_cells", 24),
        ("samples_tsv_rows_if_complete", 6),
        ("reference_tsv_rows_if_complete", 12),
        ("region_tsv_rows_if_complete", 6),
        ("full_term_audit_max_rows", 24),
    ] {
        assert_eq!(report["dimensions"][key], expected, "{key}");
    }
    assert_eq!(report["reference_sets"][0]["shared_target_sample_count"], 3);
    assert_eq!(report["reference_sets"][1]["sample_count"], 1);
    assert!(
        report["observations"]
            .as_array()
            .unwrap()
            .iter()
            .any(|s| s.as_str().unwrap().contains("fewer than two"))
    );
    for (path, hash) in report["input_sha256"].as_object().unwrap() {
        assert_eq!(
            hash.as_str().unwrap(),
            pgsum::digest::file_sha256(Path::new(path)).unwrap()
        );
    }
    assert_eq!(
        report["outputs"]["preflight.txt"],
        pgsum::digest::file_sha256(&out.join("preflight.txt")).unwrap()
    );
    let original = std::fs::read(out.join("preflight.json")).unwrap();
    assert!(!invoke(&manifest, "checked", true).status.success());
    assert_eq!(original, std::fs::read(out.join("preflight.json")).unwrap());
    let mut invalid = manifest.clone();
    invalid["reference_sets"][1]["samples"] = json!(["absent"]);
    for preflight in [true, false] {
        assert!(!invoke(&invalid, "unknown-sample", preflight).status.success());
        assert!(!fixture.path("unknown-sample").exists());
    }
    // A valid second model that was not part of extraction must fail before any score runs.
    let other_source = fixture.path("OTHER.tsv");
    std::fs::write(
        &other_source,
        std::fs::read_to_string(fixture.path("ROBUST.tsv"))
            .unwrap()
            .replace("#pgs_id=ROBUST", "#pgs_id=OTHER"),
    )
    .unwrap();
    run(&[
        "compile",
        text(&other_source),
        "--reference",
        text(&fixture.reference),
        "--out",
        text(&fixture.dir),
    ]);
    invalid = manifest.clone();
    invalid["packs"] = json!(["ROBUST.pgsp", "OTHER.pgsp"]);
    for preflight in [true, false] {
        let result = invoke(&invalid, "incompatible-pack", preflight);
        assert!(!result.status.success());
        assert!(String::from_utf8_lossy(&result.stderr).contains("not extracted with"));
        assert!(!fixture.path("incompatible-pack").exists());
    }
    let mut large_header = pgsum::frequencies::MAGIC.to_vec();
    large_header.extend(u64::MAX.to_le_bytes());
    for bytes in [pgsum::frequencies::MAGIC.to_vec(), large_header] {
        std::fs::write(fixture.path("truncated.pgsf"), bytes).unwrap();
        invalid = manifest.clone();
        invalid["scenarios"][1]["fill_frequencies"] = json!("truncated.pgsf");
        invalid["scenarios"][1]["frequency_source_note"] = json!("independent synthetic frequencies");
        let result = invoke(&invalid, "invalid-frequency", true);
        assert!(!result.status.success());
        assert!(!String::from_utf8_lossy(&result.stderr).contains("panicked"));
        assert!(!fixture.path("invalid-frequency").exists());
    }
    // Preflight does not decode all cohort blocks; this limit must not become a score-success claim.
    let opened = pgsum::cohort::CohortTable::open(&table).unwrap();
    let offset = 8 + opened.header.blocks[0].frame.offset as usize;
    drop(opened);
    let mut bytes = std::fs::read(&table).unwrap();
    bytes[offset] ^= 0xff;
    std::fs::write(&table, bytes).unwrap();
    assert!(invoke(&manifest, "lazy-block-check", true).status.success());
    assert!(!invoke(&manifest, "damaged-block-scoring", false).status.success());
    assert!(!fixture.path("damaged-block-scoring/report.json").exists());
}

#[test]
fn progress_logs_account_for_exclusions_and_preserve_all_scientific_outputs() {
    let fixture = Fixture::new("progress", 4);
    let mut fields = gt_fields();
    fields[1][0] = "./.";
    let table = fixture.table("partial", &["A", "B", "C"], "GT", &fields, "gt");
    let source = fixture.path("experiment.json");
    let manifest = fixture.manifest(
        &table,
        json!([
            {"id":"control"},
            {"id":"mask", "masks":[{"kind":"random", "rate":0.2, "seed":42, "unit":"sample_site"}]}
        ]),
    );
    std::fs::write(&source, serde_json::to_vec(&manifest).unwrap()).unwrap();
    for format in ["quiet", "json", "human"] {
        let out = fixture.path(format);
        let output = run(&[
            "robustness",
            "--experiment",
            text(&source),
            "--out",
            text(&out),
            "--terms",
            "--progress",
            format,
            "--progress-interval",
            "1",
        ]);
        assert!(output.stdout.is_empty());
        let stderr = String::from_utf8(output.stderr).unwrap();
        if format == "quiet" {
            assert!(stderr.is_empty());
        } else if format == "human" {
            assert!(stderr.contains("4/4 terms processed, 3 baseline, 1 excluded"));
            assert!(stderr.contains("robustness: complete;"));
        } else {
            let events: Vec<Value> = stderr.lines().map(|line| serde_json::from_str(line).unwrap()).collect();
            let mut phases: Vec<_> = events.iter().map(|event| event["phase"].as_str().unwrap()).collect();
            phases.dedup();
            assert_eq!(
                phases,
                [
                    "preparing",
                    "fingerprinting",
                    "scoring",
                    "summarizing",
                    "writing",
                    "verifying_inputs",
                    "complete"
                ]
            );
            assert!(events.iter().all(|e| e["schema"] == "pgsum-robustness-progress-v1"));
            assert!(
                events
                    .windows(2)
                    .all(|pair| pair[0]["elapsed_seconds"].as_f64() <= pair[1]["elapsed_seconds"].as_f64())
            );
            let score = &events.iter().rev().find(|event| event["phase"] == "scoring").unwrap()["score"];
            for (key, expected) in [
                ("terms_total", 4),
                ("terms_processed", 4),
                ("baseline_terms", 3),
                ("excluded_terms", 1),
                ("score_index", 1),
                ("scores_total", 1),
                ("samples", 3),
                ("scenarios", 2),
            ] {
                assert_eq!(score[key], expected, "{key}");
            }
            assert_eq!(events.last().unwrap()["report"], text(&out.join("report.html")));
        }
        for name in [
            "report.json",
            "report.html",
            "samples.tsv",
            "regions.tsv",
            "references.tsv",
            "top-changes.jsonl",
            "terms.tsv.zst",
        ] {
            assert_eq!(
                std::fs::read(out.join(name)).unwrap(),
                std::fs::read(fixture.path("quiet").join(name)).unwrap(),
                "{format}/{name}"
            );
        }
    }
}

#[test]
fn a_finished_term_scan_does_not_claim_a_failed_experiment_is_complete() {
    let fixture = Fixture::new("progress-error", 4);
    let mut fields = gt_fields();
    for row in &mut fields {
        row[0] = "./.";
    }
    let table = fixture.table("missing", &["A", "B", "C"], "GT", &fields, "gt");
    let source = fixture.path("experiment.json");
    std::fs::write(
        &source,
        serde_json::to_vec(&fixture.manifest(&table, json!([{"id":"control"}]))).unwrap(),
    )
    .unwrap();
    let out = fixture.path("failed");
    let output = Command::new(env!("CARGO_BIN_EXE_pgsum"))
        .args([
            "robustness",
            "--experiment",
            text(&source),
            "--out",
            text(&out),
            "--progress",
            "json",
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(!out.join("report.json").exists());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("no common-callable"));
    let events: Vec<Value> = stderr
        .lines()
        .filter(|line| line.starts_with('{'))
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert!(events.iter().all(|event| event["phase"] != "complete"));
    assert_eq!(events.last().unwrap()["score"]["terms_processed"], 4);
    assert_eq!(events.last().unwrap()["score"]["excluded_terms"], 4);
    let invalid = Command::new(env!("CARGO_BIN_EXE_pgsum"))
        .args([
            "robustness",
            "--experiment",
            text(&source),
            "--out",
            text(&fixture.path("invalid")),
            "--progress-interval",
            "0",
        ])
        .output()
        .unwrap();
    assert!(!invalid.status.success());
    assert!(!fixture.path("invalid").exists());
}

fn rows(path: &Path) -> Vec<Row> {
    let contents = std::fs::read_to_string(path).unwrap();
    parse_rows(&contents)
}

fn parse_rows(contents: &str) -> Vec<Row> {
    let mut lines = contents.lines();
    let headers: Vec<_> = lines.next().unwrap().split('\t').collect();
    lines
        .map(|line| {
            let values: Vec<_> = line.split('\t').collect();
            assert_eq!(values.len(), headers.len());
            headers
                .iter()
                .zip(values)
                .map(|(name, value)| ((*name).into(), value.into()))
                .collect()
        })
        .collect()
}

fn row<'a>(rows: &'a [Row], scenario: &str, sample: &str) -> &'a Row {
    rows.iter()
        .find(|row| row["scenario"] == scenario && row["sample"] == sample)
        .unwrap()
}

fn exact(actual: &str, expected: &str) {
    assert_eq!(
        Decimal::parse(actual)
            .unwrap()
            .numeric_cmp(&Decimal::parse(expected).unwrap()),
        Ordering::Equal,
        "{actual} != {expected}"
    );
}

fn approximately(actual: &str, expected: f64) {
    let value: f64 = actual.parse().unwrap();
    assert!((value - expected).abs() < 1e-10, "{actual} != {expected}");
}

#[test]
fn explicit_groups_preserve_global_cutoffs_reference_distributions_and_small_group_limits() {
    let fixture = Fixture::new("groups", 4);
    let table = fixture.standard_table();
    let mut manifest = fixture.manifest(
        &table,
        json!([
            {"id":"control"}, {"id":"all","masks":[{"kind":"random","rate":1.0,"seed":42,"unit":"sample_site"}]}
        ]),
    );
    manifest["reference_sets"] = json!([{"id":"same_cohort","samples":["A","B","C"]}]);
    let (native, _) = fixture.experiment("native", manifest);
    let mapping = fixture.path("groups.tsv");
    std::fs::write(&mapping, "sample\tgroup\nC\t<script>\nA\tsingleton\nB\t<script>\n").unwrap();
    let out = fixture.path("grouped");
    run(&[
        "robustness-groups",
        "--results",
        text(&native),
        "--groups",
        text(&mapping),
        "--out",
        text(&out),
    ]);
    let report: Value = serde_json::from_slice(&std::fs::read(out.join("groups.json")).unwrap()).unwrap();
    assert_eq!(report["percentile_basis"], "full_target_baseline");
    let group_rows = report["scores"][0]["scenarios"].as_array().unwrap();
    let upper = group_rows
        .iter()
        .find(|r| r["scenario"] == "all" && r["group"] == "<script>")
        .unwrap();
    let lower = group_rows
        .iter()
        .find(|r| r["scenario"] == "all" && r["group"] == "singleton")
        .unwrap();
    approximately(&upper["metrics"]["mae"].to_string(), 4.3);
    approximately(
        &upper["metrics"]["rmse"].to_string(),
        ((1.7_f64.powi(2) + 6.9_f64.powi(2)) / 2.0).sqrt(),
    );
    approximately(&upper["metrics"]["mean_abs_percentile_shift"].to_string(), 200.0 / 3.0);
    assert_eq!(upper["metrics"]["crossings"][0]["down"], 2);
    assert_eq!(upper["metrics"]["crossings"][0]["evaluable"], 2);
    assert!(upper["metrics"]["spearman"].is_null());
    assert!(lower["metrics"]["pearson"].is_null());
    assert!(lower["metrics"]["unavailable_reasons"]["pearson"].is_string());
    approximately(&lower["metrics"]["mean_abs_percentile_shift"].to_string(), 100.0 / 6.0);
    let reference = &upper["reference_comparisons"][0];
    assert_eq!(reference["reference_samples"], 3);
    assert_eq!(reference["source_cohort_shared_sample_ids"], 3);
    assert_eq!(reference["evaluable"], 2);
    approximately(&reference["mean_absolute_shift"].to_string(), 200.0 / 3.0);
    let html = std::fs::read_to_string(out.join("groups.html")).unwrap();
    assert!(html.contains("&lt;script&gt;"));
    assert!(!html.contains("<td><script>"));
    assert!(!html.contains("data-group=\"<script>\""));
    for (name, hash) in report["outputs"].as_object().unwrap() {
        assert_eq!(
            pgsum::digest::file_sha256(&out.join(name)).unwrap(),
            hash.as_str().unwrap()
        );
    }
    let repeat = Command::new(env!("CARGO_BIN_EXE_pgsum"))
        .args([
            "robustness-groups",
            "--results",
            text(&native),
            "--groups",
            text(&mapping),
            "--out",
            text(&out),
        ])
        .output()
        .unwrap();
    assert!(!repeat.status.success());
    for (i, contents) in [
        "sample\tgroup\nA\tx\nB\tx\n",
        "sample\tgroup\nA\tx\nB\tx\nD\ty\n",
        "sample\tgroup\nA\tx\nA\tx\nB\ty\nC\ty\n",
        "sample\tgroup\nA\t\nB\tx\nC\ty\n",
    ]
    .iter()
    .enumerate()
    {
        let bad = fixture.path(&format!("bad-{i}.tsv"));
        std::fs::write(&bad, contents).unwrap();
        let dest = fixture.path(&format!("bad-{i}"));
        let output = Command::new(env!("CARGO_BIN_EXE_pgsum"))
            .args([
                "robustness-groups",
                "--results",
                text(&native),
                "--groups",
                text(&bad),
                "--out",
                text(&dest),
            ])
            .output()
            .unwrap();
        assert!(!output.status.success());
        assert!(!dest.exists());
    }
    let constant = fixture.table("constant", &["A", "B", "C"], "GT", &vec![vec!["0/0"; 3]; 4], "gt");
    let (constant_native, _) = fixture.experiment(
        "constant-native",
        fixture.manifest(&constant, json!([{"id":"control"}])),
    );
    let constant_out = fixture.path("constant-groups");
    run(&[
        "robustness-groups",
        "--results",
        text(&constant_native),
        "--groups",
        text(&mapping),
        "--out",
        text(&constant_out),
    ]);
    for row in rows(&constant_out.join("cutoffs.tsv")) {
        assert_eq!(row["evaluable"], "0");
        for key in ["up", "down", "unchanged", "changed_fraction"] {
            assert_eq!(row[key], "");
        }
        assert!(row["unavailable_reason"].contains("constant"));
    }
    let original = std::fs::read(native.join("samples.tsv")).unwrap();
    std::fs::write(native.join("samples.tsv"), b"corrupted").unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_pgsum"))
        .args([
            "robustness-groups",
            "--results",
            text(&native),
            "--groups",
            text(&mapping),
            "--out",
            text(&fixture.path("corrupt")),
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("output hash differs"));
    assert!(!fixture.path("corrupt").exists());
    std::fs::write(native.join("samples.tsv"), original).unwrap();
}

#[test]
fn no_mask_matches_ordinary_scoring_and_all_mask_changes_reconcile_exactly() {
    let fixture = Fixture::new("identity", 4);
    let table = fixture.standard_table();
    let (out, report) = fixture.experiment(
        "experiment",
        fixture.manifest(
            &table,
            json!([
                {"id":"identity"},
                {"id":"all", "masks":[{"kind":"random","rate":1.0,"seed":9,"unit":"sample_site"}]}
            ]),
        ),
    );
    let scored = pgsum::cohort::score_cohort(
        &pgsum::pack::Pack::open(&fixture.path("ROBUST.pgsp")).unwrap(),
        &pgsum::cohort::CohortTable::open(&table).unwrap(),
        &Default::default(),
    )
    .unwrap();
    let samples = rows(&out.join("samples.tsv"));
    for (index, (sample, expected, delta, top, other)) in [
        ("A", "0.25", "-0.25", "-1.25", "1"),
        ("B", "1.70", "-1.70", "-2", "0.30"),
        ("C", "6.90", "-6.90", "-4", "-2.90"),
    ]
    .iter()
    .enumerate()
    {
        let identity = row(&samples, "identity", sample);
        exact(&identity["baseline_score"], expected);
        exact(&identity["scenario_score"], expected);
        exact(&identity["scenario_score"], &scored.text(index));
        exact(&identity["delta"], "0");
        assert_eq!(identity["baseline_terms"], "4");
        assert_eq!(identity["retained_terms"], "4");
        let masked = row(&samples, "all", sample);
        exact(&masked["scenario_score"], "0");
        exact(&masked["delta"], delta);
        exact(&masked["top_terms_delta"], top);
        exact(&masked["other_terms_delta"], other);
        assert_eq!(masked["masked_terms"], "4");
        assert_eq!(masked["unfilled_terms"], "4");
    }
    assert_eq!(report["scores"][0]["baseline_terms"], 4);
    assert_eq!(report["scores"][0]["scenarios"][0]["metrics"]["mae"], 0.0);
    assert_eq!(report["scores"][0]["scenarios"][1]["masked_observations"], 12);
    assert!(report["scores"][0]["scenarios"][1]["metrics"]["spearman"].is_null());
    let terms =
        String::from_utf8(zstd::decode_all(std::fs::File::open(out.join("terms.tsv.zst")).unwrap()).unwrap()).unwrap();
    let changed = parse_rows(&terms);
    assert_eq!(changed.len(), 8); // A has two nonzero terms; B and C have three each.
    assert!(changed.iter().all(|row| row["scenario"] == "all"));
    for (name, digest) in report["outputs"].as_object().unwrap() {
        assert_eq!(
            pgsum::digest::file_sha256(&out.join(name)).unwrap(),
            digest.as_str().unwrap()
        );
    }
}

#[test]
fn ds_and_gp_measurements_preserve_ref_effect_contributions() {
    let fixture = Fixture::new("dosage", 2);
    let fields = vec![
        vec!["1/1:0.125:0.2,0.3,0.5", "0/0:1.5:0.1,0.3,0.6", "0/1:0.5:0.6,0.3,0.1"],
        vec!["0/0:1.25:0.5,0.3,0.2", "1/1:0.5:0.6,0.3,0.1", "0/1:1:0,1,0"],
    ];
    let bed = fixture.bed("ref-effect.bed", "1\t2\t3\n");
    for (field, expected, retained) in [
        ("ds", ["-0.3500", "-0.45", "-0.4"], ["0.0250", "0.3", "0.1"]),
        ("gp", ["-0.39", "-0.45", "-0.4"], ["0.26", "0.3", "0.1"]),
    ] {
        let table = fixture.table(field, &["A", "B", "C"], "GT:DS:GP", &fields, field);
        let (out, report) = fixture.experiment(
            &format!("experiment-{field}"),
            fixture.manifest(
                &table,
                json!([
                    {"id":"identity"},
                    {"id":"drop-ref", "masks":[{"kind":"bed","path":bed,"mode":"drop"}]}
                ]),
            ),
        );
        let samples = rows(&out.join("samples.tsv"));
        for (index, sample) in ["A", "B", "C"].iter().enumerate() {
            exact(&row(&samples, "identity", sample)["scenario_score"], expected[index]);
            exact(&row(&samples, "drop-ref", sample)["scenario_score"], retained[index]);
            assert_eq!(row(&samples, "drop-ref", sample)["masked_terms"], "1");
        }
        assert_eq!(report["scores"][0]["baseline_terms"], 2);
    }
}

#[test]
fn bed_keep_drop_and_region_attribution_obey_half_open_boundaries() {
    let fixture = Fixture::new("bed-boundaries", 4);
    let table = fixture.standard_table();
    // Includes one-based positions 3 and 4; excludes positions 2 and 5.
    fixture.bed("middle.bed", "chr1\t2\t4\n");
    let mut manifest = fixture.manifest(
        &table,
        json!([
            {"id":"drop", "masks":[{"kind":"bed","path":"middle.bed","mode":"drop"}]},
            {"id":"keep", "masks":[{"kind":"bed","path":"middle.bed","mode":"keep"}]}
        ]),
    );
    manifest["regions"] = json!([{"id":"middle","path":"middle.bed"}]);
    let (out, _) = fixture.experiment("experiment", manifest);
    let samples = rows(&out.join("samples.tsv"));
    let regions = rows(&out.join("regions.tsv"));
    for (sample, dropped, kept, region_delta) in [
        ("A", "0", "0.25", "-0.25"),
        ("B", "2.2", "-0.5", "0.5"),
        ("C", "4.4", "2.5", "-2.5"),
    ] {
        exact(&row(&samples, "drop", sample)["scenario_score"], dropped);
        exact(&row(&samples, "keep", sample)["scenario_score"], kept);
        assert_eq!(row(&samples, "drop", sample)["masked_terms"], "2");
        assert_eq!(row(&samples, "keep", sample)["masked_terms"], "2");
        exact(&row(&regions, "drop", sample)["delta"], region_delta);
        exact(&row(&regions, "keep", sample)["delta"], "0");
    }
}

#[test]
fn seeded_masks_follow_sample_identity_across_cohort_reordering() {
    let fixture = Fixture::new("seed-reorder", 4);
    let table = fixture.standard_table();
    let fields: Vec<_> = gt_fields()
        .iter()
        .map(|values| vec![values[2], values[0], values[1]])
        .collect();
    let reordered = fixture.table("reordered", &["C", "A", "B"], "GT", &fields, "gt");
    let scenarios = json!([
        {"id":"sample-site", "masks":[{"kind":"random","rate":0.5,"seed":12345,"unit":"sample_site"}]},
        {"id":"site", "masks":[{"kind":"random","rate":0.5,"seed":12345,"unit":"site"}]}
    ]);
    let (first, first_report) = fixture.experiment("first", fixture.manifest(&table, scenarios.clone()));
    let (second, second_report) = fixture.experiment("second", fixture.manifest(&reordered, scenarios));
    let a = rows(&first.join("samples.tsv"));
    let b = rows(&second.join("samples.tsv"));
    for scenario in ["sample-site", "site"] {
        for sample in ["A", "B", "C"] {
            assert_eq!(row(&a, scenario, sample), row(&b, scenario, sample));
        }
    }
    assert_eq!(
        first_report["scores"][0]["baseline_term_sha256"],
        second_report["scores"][0]["baseline_term_sha256"]
    );
    let total = first_report["scores"][0]["scenarios"][0]["masked_observations"]
        .as_u64()
        .unwrap();
    assert!(
        total > 0 && total < 12,
        "mask should exercise both hidden and observed terms"
    );
    assert_eq!(
        row(&a, "site", "A")["masked_terms"],
        row(&a, "site", "B")["masked_terms"]
    );
    assert_eq!(
        row(&a, "site", "B")["masked_terms"],
        row(&a, "site", "C")["masked_terms"]
    );
}

#[test]
fn external_frequency_fill_uses_the_effect_allele_and_reconciles_changes() {
    let fixture = Fixture::new("frequency-fill", 4);
    let table = fixture.standard_table();
    let population = fixture.path("population.vcf.gz");
    let mut gzip =
        flate2::write::GzEncoder::new(std::fs::File::create(&population).unwrap(), flate2::Compression::fast());
    gzip.write_all(b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t2\t.\tC\tT\t.\tPASS\tAF=0.25\n1\t3\t.\tG\tT\t.\tPASS\tAF=0.25\n1\t4\t.\tT\tC\t.\tPASS\tAF=0.5\n1\t5\t.\tT\tC\t.\tPASS\tAF=0.75\n").unwrap();
    gzip.finish().unwrap();
    let frequencies = fixture.path("population.pgsf");
    run(&[
        "frequencies",
        "--vcf",
        text(&population),
        "--field",
        "AF",
        "--out",
        text(&frequencies),
    ]);
    let zero_population = fixture.path("zero.vcf.gz");
    let mut gzip = flate2::write::GzEncoder::new(
        std::fs::File::create(&zero_population).unwrap(),
        flate2::Compression::fast(),
    );
    gzip.write_all(b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        .unwrap();
    for &(pos, reference, alt, _, _, _) in &SITES {
        writeln!(gzip, "1\t{pos}\t.\t{reference}\t{alt}\t.\tPASS\tAF=0").unwrap();
    }
    gzip.finish().unwrap();
    run(&[
        "frequencies",
        "--vcf",
        text(&zero_population),
        "--field",
        "AF",
        "--out",
        text(&fixture.path("zero.pgsf")),
    ]);
    let (out, report) = fixture.experiment(
        "experiment",
        fixture.manifest(
            &table,
            json!([{
                "id":"fill", "masks":[{"kind":"random","rate":1.0,"seed":5,"unit":"site"}],
                "fill_frequencies":"population.pgsf", "frequency_source_note":"independent synthetic population; fixed before masking"
            }, {
                "id":"different_input", "masks":[{"kind":"random","rate":1.0,"seed":5,"unit":"site"}],
                "fill_frequencies":"zero.pgsf", "frequency_source_note":"second synthetic population"
            }, {
                "id":"same_input_again", "masks":[{"kind":"random","rate":1.0,"seed":5,"unit":"site"}],
                "fill_frequencies":"population.pgsf", "frequency_source_note":"same frozen population reused"
            }]),
        ),
    );
    let samples = rows(&out.join("samples.tsv"));
    // Expected contributions: 0.1 - 0.75 + 1.25 + 1 = 1.6 for every masked sample.
    for (sample, delta, top, other) in [
        ("A", "1.35", "1", "0.35"),
        ("B", "-0.10", "1.25", "-1.35"),
        ("C", "-5.30", "-3", "-2.30"),
    ] {
        let sample = row(&samples, "fill", sample);
        exact(&sample["scenario_score"], "1.6");
        exact(&sample["delta"], delta);
        exact(&sample["top_terms_delta"], top);
        exact(&sample["other_terms_delta"], other);
        assert_eq!(sample["filled_terms"], "4");
        assert_eq!(sample["unfilled_terms"], "0");
        assert_eq!(sample["retained_terms"], "0");
    }
    for name in ["A", "B", "C"] {
        // Repeated inputs share storage, while a different input retains its own frequencies.
        exact(&row(&samples, "same_input_again", name)["scenario_score"], "1.6");
        exact(&row(&samples, "different_input", name)["scenario_score"], "3");
        for key in [
            "delta",
            "filled_terms",
            "changed_terms",
            "top_terms_delta",
            "other_terms_delta",
        ] {
            assert_eq!(
                row(&samples, "fill", name)[key],
                row(&samples, "same_input_again", name)[key]
            );
        }
    }
    assert_eq!(report["scores"][0]["scenarios"][0]["filled_observations"], 12);
    assert!(
        report["input_sha256"]
            .as_object()
            .unwrap()
            .contains_key(text(&frequencies))
    );
}

#[test]
fn alternate_genotype_tables_align_samples_and_keep_missingness_separate() {
    let fixture = Fixture::new("alternate", 4);
    let table = fixture.standard_table();
    let mut fields = gt_fields();
    fields[0][0] = "1/1"; // A: +0.4 at ALT effect term.
    fields[1][1] = "./."; // B: remove -0.5 at REF effect term.
    let fields: Vec<_> = fields
        .iter()
        .map(|values| vec![values[2], values[0], values[1]])
        .collect();
    let alternate = fixture.table("alternate", &["C", "A", "B"], "GT", &fields, "gt");
    let (out, report) = fixture.experiment(
        "experiment",
        fixture.manifest(&table, json!([{"id":"alternate", "genotypes":alternate}])),
    );
    let samples = rows(&out.join("samples.tsv"));
    for (sample, expected, delta, missing) in [
        ("A", "0.65", "0.4", "0"),
        ("B", "2.2", "0.5", "1"),
        ("C", "6.9", "0", "0"),
    ] {
        let sample = row(&samples, "alternate", sample);
        exact(&sample["scenario_score"], expected);
        exact(&sample["delta"], delta);
        assert_eq!(sample["unavailable_alternate_terms"], missing);
        assert_eq!(sample["masked_terms"], "0");
    }
    assert_eq!(report["scores"][0]["baseline_terms"], 4);
    assert_eq!(
        report["scores"][0]["scenarios"][0]["unavailable_alternate_observations"],
        1
    );
}

#[test]
fn reference_distributions_are_frozen_and_never_reranked_with_the_scenario() {
    let fixture = Fixture::new("references", 4);
    let table = fixture.standard_table();
    fixture.bed("ref-term.bed", "1\t2\t3\n");
    let mut manifest = fixture.manifest(
        &table,
        json!([{"id":"drop", "masks":[{"kind":"bed","path":"ref-term.bed","mode":"drop"}]}]),
    );
    manifest["reference_sets"] = json!([
        {"id":"low", "samples":["A","B"]},
        {"id":"high", "samples":["B","C"]}
    ]);
    let (out, report) = fixture.experiment("experiment", manifest);
    let references = rows(&out.join("references.tsv"));
    let samples = rows(&out.join("samples.tsv"));
    assert_eq!(report["scores"][0]["baseline_terms"], 4);
    for (sample, expected) in [("A", "1.25"), ("B", "2.2"), ("C", "6.9")] {
        exact(&row(&samples, "drop", sample)["scenario_score"], expected);
    }
    for (reference, sample, before, after) in [
        ("low", "A", 25.0, 50.0),
        ("low", "B", 75.0, 100.0),
        ("low", "C", 100.0, 100.0),
        ("high", "A", 0.0, 0.0),
        ("high", "B", 25.0, 50.0),
        ("high", "C", 75.0, 75.0),
    ] {
        let value = references
            .iter()
            .find(|row| row["reference"] == reference && row["sample"] == sample)
            .unwrap();
        approximately(&value["baseline_percentile"], before);
        approximately(&value["scenario_percentile"], after);
        approximately(&value["percentile_delta"], after - before);
        assert_eq!(value["reference_samples"], "2");
        assert_eq!(value["shared_sample_ids"], "2");
    }
}

#[test]
fn incomplete_external_reference_freezes_one_shared_variant_universe() {
    let fixture = Fixture::new("reference-common-terms", 4);
    let table = fixture.standard_table();
    let fields = vec![
        vec!["0/0", "1/1"],
        vec!["0/0", "1/1"],
        vec!["0/1", "./."], // This term must disappear from baseline and all comparisons.
        vec!["1/1", "0/0"],
    ];
    let reference = fixture.table("external", &["R1", "R2"], "GT", &fields, "gt");
    let mut manifest = fixture.manifest(
        &table,
        json!([
            {"id":"identity"},
            {"id":"all", "masks":[{"kind":"random","rate":1.0,"seed":6,"unit":"site"}]}
        ]),
    );
    manifest["reference_sets"] = json!([{"id":"external","genotypes":reference}]);
    let (out, report) = fixture.experiment("experiment", manifest);
    assert_eq!(report["scores"][0]["baseline_terms"], 3);
    assert_eq!(report["scores"][0]["baseline_term_fraction"], 0.75);
    assert_eq!(
        report["scores"][0]["excluded_terms"]["reference_not_common_callable"],
        1
    );
    let samples = rows(&out.join("samples.tsv"));
    for (sample, expected, delta) in [("A", "-1", "1"), ("B", "1.7", "-1.7"), ("C", "4.4", "-4.4")] {
        let identity = row(&samples, "identity", sample);
        exact(&identity["baseline_score"], expected);
        exact(&identity["scenario_score"], expected);
        assert_eq!(identity["baseline_terms"], "3");
        let hidden = row(&samples, "all", sample);
        exact(&hidden["baseline_score"], expected);
        exact(&hidden["scenario_score"], "0");
        exact(&hidden["delta"], delta);
        assert_eq!(hidden["masked_terms"], "3");
    }
    for reference in rows(&out.join("references.tsv")) {
        assert_eq!(reference["reference_samples"], "2");
        assert_eq!(reference["shared_sample_ids"], "0");
    }
}

#[test]
fn single_pgsg_preserves_exact_scores_and_withholds_unsupported_rank_statistics() {
    let fixture = Fixture::new("single-pgsg", 4);
    fixture.standard_table();
    let single = fixture.path("single.pgsg");
    run(&[
        "extract",
        "--gvcf",
        text(&fixture.path("baseline.vcf")),
        "--reference",
        text(&fixture.reference),
        "--pack",
        text(&fixture.path("ROBUST.pgsp")),
        "--sample",
        "A",
        "--out",
        text(&single),
        "--accept-missing-quality",
    ]);
    // Confirm the real single-sample format is used, rather than a one-row cohort.
    let table = pgsum::genotypes::GenotypeTable::open(&single).unwrap();
    assert_eq!(table.header.sample.sample_id, "A");
    let ordinary = pgsum::score::score(
        &pgsum::pack::Pack::open(&fixture.path("ROBUST.pgsp")).unwrap(),
        &table,
        &Default::default(),
        None,
    )
    .unwrap();
    let (out, report) = fixture.experiment(
        "experiment",
        fixture.manifest(
            &single,
            json!([
                {"id":"identity"},
                {"id":"all", "masks":[{"kind":"random","rate":1.0,"seed":1,"unit":"site"}]}
            ]),
        ),
    );
    let samples = rows(&out.join("samples.tsv"));
    assert_eq!(samples.len(), 2);
    let identity = row(&samples, "identity", "A");
    exact(&identity["baseline_score"], "0.25");
    exact(&identity["scenario_score"], &ordinary.partial.raw_score);
    exact(&identity["delta"], "0");
    let masked = row(&samples, "all", "A");
    exact(&masked["baseline_score"], "0.25");
    exact(&masked["scenario_score"], "0");
    exact(&masked["delta"], "-0.25");
    exact(&masked["top_terms_delta"], "-1.25");
    exact(&masked["other_terms_delta"], "1");
    for sample in &samples {
        for column in [
            "baseline_percentile",
            "fixed_baseline_percentile",
            "reranked_percentile",
        ] {
            assert_eq!(sample[column], "");
        }
    }
    for (scenario, expected_error) in report["scores"][0]["scenarios"]
        .as_array()
        .unwrap()
        .iter()
        .zip([0.0, 0.25])
    {
        let metrics = &scenario["metrics"];
        assert_eq!(metrics["n"], 1);
        for metric in ["mae", "rmse", "max_abs_change"] {
            assert_eq!(metrics[metric], expected_error);
        }
        for metric in [
            "pearson",
            "spearman",
            "mean_abs_percentile_shift",
            "max_abs_percentile_shift",
        ] {
            assert!(metrics[metric].is_null());
            assert!(!metrics["unavailable_reasons"][metric].as_str().unwrap().is_empty());
        }
        for crossing in metrics["crossings"].as_array().unwrap() {
            for count in ["up", "down", "unchanged", "evaluable"] {
                assert_eq!(crossing[count], 0);
            }
        }
    }
}

#[test]
fn completed_experiment_outputs_are_not_overwritten_by_a_rerun() {
    let fixture = Fixture::new("no-overwrite", 4);
    let table = fixture.standard_table();
    let (out, _) = fixture.experiment("experiment", fixture.manifest(&table, json!([{"id":"identity"}])));
    let snapshot = || -> BTreeMap<_, _> {
        std::fs::read_dir(&out)
            .unwrap()
            .map(|entry| {
                let entry = entry.unwrap();
                (entry.file_name(), std::fs::read(entry.path()).unwrap())
            })
            .collect()
    };
    let before = snapshot();
    assert!(before.contains_key(std::ffi::OsStr::new("report.json")));
    let rejected = Command::new(env!("CARGO_BIN_EXE_pgsum"))
        .args([
            "robustness",
            "--experiment",
            text(&fixture.path("experiment.json")),
            "--out",
            text(&out),
            "--terms",
        ])
        .output()
        .unwrap();
    assert!(!rejected.status.success(), "existing outputs must not be replaced");
    assert!(String::from_utf8_lossy(&rejected.stderr).contains(text(&out)));
    assert_eq!(
        snapshot(),
        before,
        "a rejected rerun must leave every completed output intact"
    );
}

#[test]
fn planned_sweeps_execute_paired_masks_and_preserve_controls_across_directories() {
    let fixture = Fixture::new("planned-sweep", 4);
    let table = fixture.standard_table();
    let template = fixture.path("policies.json");
    let mut manifest = fixture.manifest(&table, json!([{"id":"omit"},{"id":"same-policy"}]));
    manifest["cutoffs"] = json!([50, 90]);
    manifest["reference_sets"] = json!([{"id":"cohort","samples":["A","B","C"]},{"id":"singleton","samples":["A"]}]);
    std::fs::write(&template, serde_json::to_vec(&manifest).unwrap()).unwrap();
    let planned = fixture.path("planned experiments");
    let command = || {
        Command::new(env!("CARGO_BIN_EXE_pgsum"))
            .current_dir(std::env::temp_dir())
            .args([
                "robustness-plan",
                "--template",
                text(&template),
                "--rates",
                "0,0.5,1",
                "--seeds",
                "42,43",
                "--out",
                text(&planned),
            ])
            .output()
            .unwrap()
    };
    let output = command();
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    let experiment = planned.join("experiment.json");
    let before = std::fs::read(&experiment).unwrap();
    assert!(!command().status.success());
    assert_eq!(before, std::fs::read(&experiment).unwrap());
    let plan: Value = serde_json::from_slice(&std::fs::read(planned.join("plan.json")).unwrap()).unwrap();
    assert_eq!(plan["scenarios"].as_array().unwrap().len(), 14);
    assert_eq!(
        plan["experiment_sha256"],
        pgsum::digest::file_sha256(&experiment).unwrap()
    );
    let loaded = pgsum::robustness::config::load(&experiment).unwrap();
    assert!(loaded.genotypes.is_absolute());
    let out = fixture.path("sweep results");
    run(&["robustness", "--experiment", text(&experiment), "--out", text(&out)]);
    let summarized = fixture.path("study summary");
    run(&[
        "robustness-summarize",
        "--plan",
        text(&planned),
        "--results",
        text(&out),
        "--reference-policy",
        "omit",
        "--out",
        text(&summarized),
    ]);
    let summary: Value = serde_json::from_slice(&std::fs::read(summarized.join("summary.json")).unwrap()).unwrap();
    assert_eq!(summary["scores"][0]["summaries"].as_array().unwrap().len(), 42);
    for comparison in summary["scores"][0]["comparisons"].as_array().unwrap() {
        if comparison["distribution"]["available"].as_u64().unwrap() > 0 {
            assert_eq!(comparison["distribution"]["mean"], 0.0);
        }
    }
    let percentile_rows = summary["scores"][0]["percentile_comparisons"].as_array().unwrap();
    assert!(!percentile_rows.is_empty());
    for row in percentile_rows {
        if row["reference_set"] == "singleton" {
            assert_eq!(row["distribution"]["available"], 0);
            assert!(row["distribution"]["mean"].is_null());
        } else {
            assert_eq!(row["distribution"]["available"], 2);
            assert_eq!(row["distribution"]["mean"], 0.0);
        }
    }
    for name in ["percentile_summaries.tsv", "percentile_comparisons.tsv"] {
        assert_eq!(
            summary["outputs"][name],
            pgsum::digest::file_sha256(&summarized.join(name)).unwrap()
        );
    }
    let mapping = fixture.path("groups.tsv");
    std::fs::write(&mapping, "sample\tgroup\nA\tsingle\nB\t<script>\nC\t<script>\n").unwrap();
    let group_out = fixture.path("group study");
    run(&[
        "robustness-summarize",
        "--plan",
        text(&planned),
        "--results",
        text(&out),
        "--groups",
        text(&mapping),
        "--reference-policy",
        "omit",
        "--out",
        text(&group_out),
    ]);
    let grouped: Value = serde_json::from_slice(&std::fs::read(group_out.join("summary.json")).unwrap()).unwrap();
    assert_eq!(grouped["schema"], "pgsum-robustness-group-study-v1");
    assert_eq!(grouped["percentile_basis"], "full_target_baseline");
    assert_eq!(grouped["scores"].as_array().unwrap().len(), 2);
    for score in grouped["scores"].as_array().unwrap() {
        assert_eq!(score["source_samples"], 3);
        assert_eq!(score["reference_sets"][0]["source_cohort_shared_sample_ids"], 3);
        assert!(score["reference_sets"][0].get("shared_sample_ids").is_none());
        for row in score["comparisons"].as_array().unwrap() {
            if score["group"] == "single" && ["pearson", "spearman"].contains(&row["metric"].as_str().unwrap()) {
                assert_eq!(row["distribution"]["available"], 0);
                assert_eq!(row["distribution"]["unavailable"], 2);
            } else if row["distribution"]["available"].as_u64().unwrap() > 0 {
                assert_eq!(row["distribution"]["mean"], 0.0);
            }
        }
        if score["group"] == "<script>" {
            let row = score["percentile_summaries"]
                .as_array()
                .unwrap()
                .iter()
                .find(|r| {
                    r["policy"] == "omit" && r["rate"] == 1.0 && r["cutoff"] == 50.0 && r["metric"] == "down_count"
                })
                .unwrap();
            assert_eq!(row["distribution"]["mean"], 2.0);
            assert_eq!(row["seeds"][0]["evaluable"], 2);
        }
    }
    let html = std::fs::read_to_string(group_out.join("summary.html")).unwrap();
    assert!(html.contains("&lt;script&gt;") && !html.contains("<td><script>"));
    let tabular = rows(&group_out.join("summaries.tsv"));
    assert_eq!(tabular[0]["source_samples"], "3");
    assert_eq!(tabular.len(), 84);
    for (name, hash) in grouped["outputs"].as_object().unwrap() {
        assert_eq!(
            hash.as_str().unwrap(),
            pgsum::digest::file_sha256(&group_out.join(name)).unwrap()
        );
    }
    assert!(
        grouped["source_sha256"]
            .as_object()
            .unwrap()
            .values()
            .any(|h| h.as_str() == Some(pgsum::digest::file_sha256(&mapping).unwrap().as_str()))
    );
    // Corrupted output cannot be summarized even if report.json still exists.
    std::fs::write(out.join("references.tsv"), "corrupted").unwrap();
    let rejected = fixture.path("rejected study");
    let result = Command::new(env!("CARGO_BIN_EXE_pgsum"))
        .args([
            "robustness-summarize",
            "--plan",
            text(&planned),
            "--results",
            text(&out),
            "--reference-policy",
            "omit",
            "--out",
            text(&rejected),
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(!rejected.exists());
    let samples = rows(&out.join("samples.tsv"));
    assert_eq!(samples.len(), 42);
    for sample in ["A", "B", "C"] {
        for seed in [42, 43] {
            for rate in ["0", "0.5", "1"] {
                let a = row(&samples, &format!("omit/rate={rate}/seed={seed}"), sample);
                let b = row(&samples, &format!("same-policy/rate={rate}/seed={seed}"), sample);
                for key in ["masked_terms", "scenario_score", "delta", "baseline_score"] {
                    assert_eq!(a[key], b[key]);
                }
                if rate == "0" {
                    assert_eq!(a["scenario_score"], row(&samples, "omit", sample)["scenario_score"]);
                    assert_eq!(a["masked_terms"], "0");
                } else if rate == "1" {
                    exact(&a["scenario_score"], "0");
                    assert_eq!(a["masked_terms"], "4");
                }
            }
        }
    }
}
