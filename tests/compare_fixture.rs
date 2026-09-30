//! Migration gates with hand-computed sums, external-format fixtures and real pgsum exports.
use serde_json::Value;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};
static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "pgsum-compare-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
    fn write(&self, name: &str, text: &str) {
        std::fs::write(self.path(name), text).unwrap();
    }
    fn compare(&self, extra: &[&str]) -> Output {
        Command::new(env!("CARGO_BIN_EXE_pgsum"))
            .current_dir(&self.0)
            .args(["compare", "--left", "left", "--right", "right", "--out", "out"])
            .args(extra)
            .output()
            .unwrap()
    }
    fn report(&self) -> Value {
        serde_json::from_slice(&std::fs::read(self.path("out/comparison.json")).unwrap()).unwrap()
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
fn ok(o: Output) {
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
}
fn err(o: Output, text: &str) {
    assert!(!o.status.success());
    assert!(
        String::from_utf8_lossy(&o.stderr).contains(text),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
}
const HEADER: &str = "sample\tpgs_id\tpartial_raw_score\n";

#[test]
fn exact_decimal_tolerance_full_outer_join_and_unavailable() {
    let f = Fixture::new();
    f.write("left", &format!("{HEADER}S\tEXACT\t9007199254740993.0100\nS\tABS\t-1.1\nS\tREL\t100\nS\tBAD\t1e-1000\nS\tLEFT\t3\nS\tNA\t\nS\tZERO\t-0.00\nS\tSCALE\t1\n"));
    f.write("right", &format!("{HEADER}S\tEXACT\t900719925474099301e-2\nS\tABS\t-1.09\nS\tREL\t101\nS\tBAD\t1\nS\tRIGHT\t3\nS\tNA\tNA\nS\tZERO\t0\nS\tSCALE\t1.010000000000000000000000000000000000000000000000001\n"));
    err(
        f.compare(&[
            "--absolute-tolerance",
            "0.01",
            "--relative-tolerance",
            "0.01",
            "--fail-on-difference",
        ]),
        "comparison differs",
    );
    let r = f.report();
    assert_eq!(
        r["counts"],
        serde_json::json!({"rows":9,"exact":2,"within_tolerance":3,"different":1,"left_only":1,"right_only":1,"unavailable":1,"identity_mismatch":0})
    );
    assert_eq!(r["agrees"], false);
    let tsv = std::fs::read_to_string(f.path("out/scores.tsv")).unwrap();
    assert!(tsv.contains("S\tABS\twithin_tolerance\t-1.1\t-1.09\t1e-2"));
    assert!(tsv.contains("S\tREL\twithin_tolerance\t100\t101\t1e0"));
    assert_eq!(
        r["outputs"]["scores.tsv"],
        pgsum::digest::file_sha256(&f.path("out/scores.tsv")).unwrap()
    );
    assert!(!f.compare(&[]).status.success()); // Never overwrite the completed report.
    assert_eq!(f.report(), r);
}

#[test]
fn tolerance_boundary_is_not_rounded_through_float() {
    let f = Fixture::new();
    f.write("left", &format!("{HEADER}S\tA\t0\nS\tB\t0\nS\tC\t1e-1000\n"));
    f.write(
        "right",
        &format!("{HEADER}S\tA\t0.01\nS\tB\t0.01000000000000000000000000000000000001\nS\tC\t2e-1000\n"),
    );
    ok(f.compare(&["--absolute-tolerance", "0.01"]));
    assert_eq!(f.report()["counts"]["within_tolerance"], 2);
    assert_eq!(f.report()["counts"]["different"], 1);
    let g = Fixture::new();
    g.write("left", &format!("{HEADER}S\tA\t1e-1000\n"));
    g.write("right", &format!("{HEADER}S\tA\t2e-1000\n"));
    ok(g.compare(&[]));
    assert!(
        std::fs::read_to_string(g.path("out/scores.tsv"))
            .unwrap()
            .contains("\t1e-1000\t")
    );
    assert_eq!(g.report()["counts"]["different"], 1);
}

#[test]
fn plink_and_pgsc_adapters_compression_and_identity() {
    let f = Fixture::new();
    let text = b"#FID IID ALLELE_CT NAMED_ALLELE_DOSAGE_SUM PGS1_SUM PGS1_AVG\nF S 4 2 0.300000 0.075\n";
    std::fs::write(f.path("left"), zstd::encode_all(&text[..], 1).unwrap()).unwrap();
    let text = b"sampleset FID IID PGS SUM DENOM AVG\ntarget F S PGS1 0.3 4 0.075\nreference F S PGS1 8 4 2\n";
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    std::io::Write::write_all(&mut gz, text).unwrap();
    std::fs::write(f.path("right"), gz.finish().unwrap()).unwrap();
    let options = [
        "--left-format",
        "plink2",
        "--left-score-column",
        "PGS1_SUM",
        "--left-score-id",
        "PGS1",
        "--right-format",
        "pgsc-calc",
        "--right-sampleset",
        "target",
        "--fail-on-difference",
    ];
    ok(f.compare(&options));
    let r = f.report();
    assert_eq!(r["counts"]["exact"], 1);
    assert_eq!(r["right"]["filtered_rows"], 1);
    assert_eq!(r["left"]["value_column"], "PGS1_SUM");
    assert_eq!(
        r["left"]["sha256"],
        pgsum::digest::file_sha256(&f.path("left")).unwrap()
    );

    let g = Fixture::new();
    g.write("left", "#FID IID X_SUM\nF S 3\n");
    g.write("right", "#FID IID X_SUM\nOTHER S 3\n");
    err(
        g.compare(&[
            "--left-format",
            "plink2",
            "--right-format",
            "plink2",
            "--left-score-column",
            "X_SUM",
            "--right-score-column",
            "X_SUM",
            "--left-score-id",
            "X",
            "--right-score-id",
            "X",
            "--fail-on-difference",
        ]),
        "comparison differs",
    );
    assert_eq!(g.report()["counts"]["identity_mismatch"], 1);
}

#[test]
fn malformed_ambiguous_and_non_sum_inputs_fail_without_a_completion_marker() {
    let cases: Vec<(&str, &str, Vec<&str>, &str)> = vec![
        (
            "sample\tpgs_id\tpartial_raw_score\nS\tP\t1\nS\tP\t2\n",
            "",
            vec![],
            "duplicate sample/score",
        ),
        ("pgs_id\tpartial_raw_score\nP\t1\n", "", vec![], "requires an explicit"),
        (
            "sample\tpgs_id\tpartial_raw_score\nS\tP\tinf\n",
            "",
            vec![],
            "invalid finite",
        ),
        (
            "sample\tpgs_id\tpartial_raw_score\nS\tP\t1e999999999\n",
            "",
            vec![],
            "invalid finite",
        ),
        (
            "sample\tpgs_id\tpartial_raw_score\nS\tP\t1\tX\n",
            "",
            vec![],
            "expected 3 columns",
        ),
        (
            "sample\tpgs_id\tpartial_raw_score\tpgs_id\n",
            "",
            vec![],
            "duplicate column",
        ),
        (
            "#IID X_AVG\nS 1\n",
            "",
            vec![
                "--left-format",
                "plink2",
                "--left-score-column",
                "X_AVG",
                "--left-score-id",
                "P",
            ],
            "averages are not sums",
        ),
        (
            "#FID IID X_SUM\nF S 1\nG S 2\n",
            "",
            vec![
                "--left-format",
                "plink2",
                "--left-score-column",
                "X_SUM",
                "--left-score-id",
                "P",
            ],
            "different FID/SID",
        ),
        (
            "sampleset IID PGS SUM\na A P 1\nb B P 2\n",
            "",
            vec!["--left-format", "pgsc-calc"],
            "multiple samplesets",
        ),
        ("sample\tpgs_id\tpartial_raw_score\n", "", vec![], "no score rows"),
    ];
    for (left, right, options, message) in cases {
        let f = Fixture::new();
        f.write("left", left);
        f.write("right", right);
        err(f.compare(&options), message);
        assert!(!f.path("out").exists());
    }
    let f = Fixture::new();
    f.write("left", "#IID NAMED_ALLELE_DOSAGE_SUM\nS 2\n");
    f.write("right", "");
    err(
        f.compare(&[
            "--left-format",
            "plink2",
            "--left-score-column",
            "NAMED_ALLELE_DOSAGE_SUM",
            "--left-score-id",
            "P",
        ]),
        "averages are not sums",
    );
    assert!(!f.path("out").exists());
    for bad in ["-0.1", "NaN", "inf", "1e999999", "1e", "--1", "."] {
        let f = Fixture::new();
        f.write("left", "");
        f.write("right", "");
        let o = f.compare(&[&format!("--absolute-tolerance={bad}")]);
        assert!(!o.status.success());
        assert!(!f.path("out").exists());
    }
}

const TERMS: &str =
    "ordinal\tcontig\tpos\tmodel\tref\talt\teffect_is_alt\tweights\tstatus\tcall_state\teffect_dosage\tcontribution\n";
const TERM_ARGS: &[&str] = &[
    "--left-terms",
    "lt",
    "--right-terms",
    "rt",
    "--term-sample",
    "S",
    "--term-score",
    "P",
];
#[test]
fn term_evidence_explains_differences_and_reconciles_both_sums() {
    let f = Fixture::new();
    f.write("left", &format!("{HEADER}S\tP\t0.4\n"));
    f.write("right", &format!("{HEADER}S\tP\t0.6\n"));
    f.write("lt", &format!("{TERMS}1\t1\t10\tadditive\tA\tG\t1\t0.2\tscorable_observation\thet\t1\t0.2\n2\t1\t11\tadditive\tC\tT\t1\t0.1\tscorable_observation\thom_alt\t2\t0.2\n3\t1\t12\tadditive\tC\tG\t1\t0.4\tmissing\tmissing\t\t\n"));
    f.write("rt", &format!("{TERMS}1\t1\t10\tadditive\tA\tG\t1\t0.2\tscorable_observation\thom_alt\t2\t0.4\n2\t1\t11\tadditive\tC\tT\t1\t0.1\tscorable_observation\thom_alt\t2\t0.2\n3\t1\t12\tadditive\tC\tG\t1\t0.4\tmissing\tmissing\t\t\n"));
    ok(f.compare(TERM_ARGS));
    let r = f.report();
    assert_eq!(r["terms"]["left"]["score_reconciliation"], "exact");
    assert_eq!(r["terms"]["right"]["score_reconciliation"], "exact");
    assert_eq!(r["terms"]["rows_with_changed_evidence"], 1);
    assert_eq!(r["terms"]["jointly_unscored"], 1);
    let rows: Vec<Value> = std::fs::read_to_string(f.path("out/terms.jsonl"))
        .unwrap()
        .lines()
        .map(|s| serde_json::from_str(s).unwrap())
        .collect();
    assert_eq!(rows[0]["delta_right_minus_left"], "2e-1");
    assert_eq!(
        rows[0]["changed_fields"],
        serde_json::json!(["call_state", "effect_dosage"])
    );
    assert!(rows[0]["explanations"][0].as_str().unwrap().contains("dosage differs"));
    assert_eq!(rows[2]["left_contribution"], "");
}

#[test]
fn same_totals_do_not_hide_term_changes_or_incomplete_exports() {
    for (right, expected) in [
        (
            format!(
                "{TERMS}1\t1\t10\tadditive\tA\tG\t1\t2\tok\thet\t1\t2\n2\t1\t11\tadditive\tC\tT\t1\t1\tok\thet\t1\t1\n"
            ),
            "exact",
        ),
        (
            format!("{TERMS}1\t1\t10\tadditive\tA\tG\t1\t1\tok\thet\t1\t1\n"),
            "different",
        ),
    ] {
        let f = Fixture::new();
        f.write("left", &format!("{HEADER}S\tP\t3\n"));
        f.write("right", &format!("{HEADER}S\tP\t3\n"));
        f.write(
            "lt",
            &format!(
                "{TERMS}1\t1\t10\tadditive\tA\tG\t1\t1\tok\thet\t1\t1\n2\t1\t11\tadditive\tC\tT\t1\t2\tok\thet\t1\t2\n"
            ),
        );
        f.write("rt", &right);
        let mut args = TERM_ARGS.to_vec();
        args.push("--fail-on-difference");
        err(f.compare(&args), "comparison differs");
        let r = f.report();
        assert_eq!(r["counts"]["exact"], 1);
        assert_eq!(r["terms"]["agrees"], false);
        assert_eq!(r["terms"]["right"]["score_reconciliation"], expected);
    }
}

#[test]
fn out_of_order_terms_cleanup_and_original_input_is_never_overwritten() {
    let f = Fixture::new();
    f.write("left", &format!("{HEADER}S\tP\t1\n"));
    f.write("right", &format!("{HEADER}S\tP\t1\n"));
    let row = "1\t1\t10\tadditive\tA\tG\t1\t1\tok\thet\t1\t1\n";
    f.write("lt", &format!("{TERMS}{row}{row}"));
    f.write("rt", &format!("{TERMS}{row}"));
    err(f.compare(TERM_ARGS), "strictly increasing");
    assert!(!f.path("out").exists());
    let original = std::fs::read(f.path("left")).unwrap();
    let o = Command::new(env!("CARGO_BIN_EXE_pgsum"))
        .current_dir(&f.0)
        .args(["compare", "--left", "left", "--right", "right", "--out", "left"])
        .output()
        .unwrap();
    assert!(!o.status.success());
    assert_eq!(std::fs::read(f.path("left")).unwrap(), original);
}

#[test]
fn native_run_exports_round_trip_through_the_comparison_gate() {
    let f = Fixture::new();
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let o = Command::new(env!("CARGO_BIN_EXE_pgsum"))
        .args(["compile"])
        .arg(fixtures.join("PGS999998_hmPOS_GRCh38.txt.gz"))
        .arg("--reference")
        .arg(fixtures.join("synthetic.fa"))
        .arg("--out")
        .arg(f.path("packs"))
        .output()
        .unwrap();
    ok(o);
    let o = Command::new(env!("CARGO_BIN_EXE_pgsum"))
        .arg("run")
        .arg("--gvcf")
        .arg(fixtures.join("synthetic.g.vcf.gz"))
        .arg("--reference")
        .arg(fixtures.join("synthetic.fa"))
        .arg("--pack")
        .arg(f.path("packs"))
        .arg("--out")
        .arg(f.path("run"))
        .arg("--terms")
        .output()
        .unwrap();
    ok(o);
    for name in ["left", "right"] {
        std::fs::copy(f.path("run/scores.tsv"), f.path(name)).unwrap();
    }
    for name in ["lt", "rt"] {
        std::fs::copy(f.path("run/PGS999998.terms.tsv"), f.path(name)).unwrap();
    }
    ok(f.compare(&[
        "--left-sample",
        "S",
        "--right-sample",
        "S",
        "--left-terms",
        "lt",
        "--right-terms",
        "rt",
        "--term-sample",
        "S",
        "--term-score",
        "PGS999998",
        "--fail-on-difference",
    ]));
    let r = f.report();
    assert_eq!(r["agrees"], true);
    assert_eq!(r["terms"]["left"]["score_reconciliation"], "exact");
}
