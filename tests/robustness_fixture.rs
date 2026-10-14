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
    let (out, report) = fixture.experiment(
        "experiment",
        fixture.manifest(
            &table,
            json!([{
                "id":"fill", "masks":[{"kind":"random","rate":1.0,"seed":5,"unit":"site"}],
                "fill_frequencies":"population.pgsf", "frequency_source_note":"independent synthetic population; fixed before masking"
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
    let manifest = fixture.manifest(&table, json!([{"id":"omit"},{"id":"same-policy"}]));
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
