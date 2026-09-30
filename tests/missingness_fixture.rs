//! Bounds checked against every possible completion of independently specified contribution tables.
use serde_json::Value;
use std::path::PathBuf;
use std::process::Command;
struct Fixture {
    root: PathBuf,
}
impl Fixture {
    fn new(name: &str, extra: &str) -> Self {
        let root = std::env::temp_dir().join(format!("pgsum-missingness-{}-{name}", std::process::id()));
        std::fs::create_dir_all(&root).unwrap();
        let f = Self { root };
        f.write("reference.fa", ">chr1\nACGTTGCA\n");
        f.write("reference.fa.fai", "chr1\t8\t6\t8\t9\n");
        f.write("score.tsv", &format!("#pgs_id=BOUNDS\n#genome_build=GRCh38\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\tis_dominant\tis_recessive\tdosage_0_weight\tdosage_1_weight\tdosage_2_weight\tis_interaction\n1\t1\tG\tA\t-2\tFalse\tFalse\t\t\t\tFalse\n1\t2\tT\tC\t0.5\tFalse\tFalse\t\t\t\tFalse\n1\t3\tA\tG\t3\tTrue\tFalse\t\t\t\tFalse\n1\t4\tC\tT\t-1\tFalse\tTrue\t\t\t\tFalse\n1\t5\tG\tT\t\tFalse\tFalse\t-3\t7\t1\tFalse\n{extra}"));
        f.write("sample.vcf", "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n1\t1\tfiltered-site\tA\tG\t.\tLowQual\t.\tGT:DP:GQ\t0/1:30:60\n1\t2\tcalled-site\tC\tT\t.\tPASS\t.\tGT:DP:GQ\t0/1:30:60\n1\t4\tpartial-call\tT\tC\t.\tPASS\t.\tGT:DP:GQ\t0/.:30:60\n1\t5\tlow-quality\tT\tG\t.\tPASS\t.\tGT:DP:GQ\t1/1:30:10\n");
        f.cli(&["compile", "score.tsv", "--reference", "reference.fa", "--out", "packs"]);
        f
    }
    fn write(&self, file: &str, text: &str) {
        std::fs::write(self.root.join(file), text).unwrap();
    }
    fn cli(&self, args: &[&str]) {
        let result = Command::new(env!("CARGO_BIN_EXE_pgsum"))
            .current_dir(&self.root)
            .args(args)
            .args(["--threads", "2"])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{args:?}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
    fn score(&self) -> Value {
        self.cli(&[
            "run",
            "--gvcf",
            "sample.vcf",
            "--reference",
            "reference.fa",
            "--pack",
            "packs",
            "--terms",
            "--out",
            "out",
        ]);
        serde_json::from_slice(&std::fs::read(self.root.join("out/BOUNDS.score.json")).unwrap()).unwrap()
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.root);
    }
}

#[test]
fn every_completion_fits_exact_bounds_and_sources_explain_missing_terms() {
    let f = Fixture::new("envelope", "");
    let score = f.score();
    let report = &score["missingness"];
    assert_eq!(score["partial"]["raw_score"], "0.5");
    assert_eq!(report["missing_terms"], 4);
    assert_eq!(report["bounded_terms"], 4);
    assert_eq!(report["unbounded_terms"], 0);
    assert_eq!(
        report["bounded_missing_contribution"],
        serde_json::json!({"lower":"-8", "upper":"10"})
    );
    assert_eq!(
        report["completion_score_bounds"],
        serde_json::json!({"lower":"-7.5", "upper":"10.5"})
    );
    // Values in half-units, written independently from the scorer: 0.5 observed, then four missing terms.
    let mut completions = Vec::new();
    for a in [0, -4, -8] {
        for d in [0, 6, 6] {
            for r in [0, 0, -2] {
                for w in [-6, 14, 2] {
                    completions.push(1 + a + d + r + w);
                }
            }
        }
    }
    assert_eq!(completions.iter().min(), Some(&-15));
    assert_eq!(completions.iter().max(), Some(&21));
    assert!(completions.iter().all(|v| (-15..=21).contains(v)));
    let rows = report["ranked_missing_terms"].as_array().unwrap();
    assert_eq!(
        rows.iter().map(|r| r["ordinal"].as_u64().unwrap()).collect::<Vec<_>>(),
        [5, 1, 3, 4]
    );
    assert_eq!(rows[0]["status"], "low_quality");
    assert_eq!(rows[1]["status"], "filtered");
    assert_eq!(rows[2]["status"], "unknown_no_record");
    assert_eq!(rows[3]["status"], "partial_no_call");
    assert_eq!(rows[2]["evidence"]["availability"], "no_overlapping_record");
    let retained = rows[1]["evidence"]["records"][0]["text"].as_str().unwrap();
    assert!(retained.contains("filtered-site"));
    assert_eq!(
        rows[1]["evidence"]["records"][0]["sha256"],
        pgsum::pack::records_sha256(retained.as_bytes())
    );
    let tsv = std::fs::read_to_string(f.root.join("out/BOUNDS.terms.tsv")).unwrap();
    let lines: Vec<_> = tsv.lines().map(|l| l.split('\t').collect::<Vec<_>>()).collect();
    let lo = lines[0]
        .iter()
        .position(|v| *v == "missing_contribution_lower")
        .unwrap();
    assert_eq!(&lines[5][lo..lo + 3], ["-3", "7", ""]);
    assert_eq!(&lines[2][lo..lo + 3], ["", "", ""]);
    // Thread count cannot change rankings, bounds, source selection or the raw sum.
    let output = Command::new(env!("CARGO_BIN_EXE_pgsum"))
        .current_dir(&f.root)
        .args([
            "score",
            "--genotypes",
            "out/genotypes.pgsg",
            "--pack",
            "packs",
            "--out",
            "one-thread",
            "--threads",
            "1",
        ])
        .output()
        .unwrap();
    assert!(output.status.success());
    let single: Value =
        serde_json::from_slice(&std::fs::read(f.root.join("one-thread/BOUNDS.score.json")).unwrap()).unwrap();
    assert_eq!(single, score);
}

#[test]
fn unknown_models_prevent_a_misleading_whole_score_interval() {
    let f = Fixture::new("interaction", "1\t6\tA\tG\t100\tFalse\tFalse\t\t\t\tTrue\n");
    let score = f.score();
    let report = &score["missingness"];
    assert_eq!(report["missing_terms"], 5);
    assert_eq!(report["unbounded_terms"], 1);
    assert!(report["completion_score_bounds"].is_null());
    let example = &report["unbounded_examples"][0];
    assert_eq!(example["ordinal"], 6);
    assert_eq!(
        example["bounds_unavailable_because"],
        "ambiguous_or_unsupported_contribution_model"
    );
    assert_eq!(example["evidence"]["availability"], "position_evidence_not_extracted");
    assert!(example["contribution_bounds"].is_null());
}

#[test]
fn metadata_uncertainty_withholds_bounds_and_old_json_remains_readable() {
    let f = Fixture::new("metadata", "");
    let mut result = f.score();
    result.as_object_mut().unwrap().remove("missingness");
    let old: pgsum::score::ScoreResult = serde_json::from_value(result).unwrap();
    assert!(old.missingness.is_none());
    let mut pack = pgsum::pack::Pack::open(&f.root.join("packs/BOUNDS.pgsp")).unwrap();
    let genotypes = pgsum::genotypes::GenotypeTable::open(&f.root.join("out/genotypes.pgsg")).unwrap();
    pack.header.inventory.consistent = false;
    let r = pgsum::score::score(&pack, &genotypes, &Default::default(), None).unwrap();
    let report = r.missingness.unwrap();
    assert!(report.completion_score_bounds.is_none());
    assert_eq!(
        report.completion_bounds_withheld_because,
        ["inconsistent_term_inventory"]
    );
    assert_eq!(report.bounded_terms, 4);
    pack.header.inventory.consistent = true;
    pack.header.origin = pgsum::scoring_file::Origin::PgsCatalog;
    pack.header.matches_publication = Some(false);
    let r = pgsum::score::score(&pack, &genotypes, &Default::default(), None).unwrap();
    assert_eq!(
        r.missingness.unwrap().completion_bounds_withheld_because,
        ["scoring_file_not_confirmed_to_match_publication"]
    );
}
