//! Public, self-contained validation with hand-calculated expectations. No external genomes or oracle.
//! Run with `cargo test --locked --test reference_validation` (see bench/README.md).

use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

use serde_json::Value;

struct Fixture {
    dir: PathBuf,
}

impl Fixture {
    fn new(name: &str, two_terms: bool) -> Self {
        let dir = std::env::temp_dir().join(format!("pgsum-reference-{name}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let f = Self { dir };
        f.write("ref.fa", ">chr1\nACGT\n");
        f.write("ref.fa.fai", "chr1\t4\t6\t4\t5\n");
        let first = if two_terms { "chr1\t1\tG\tA\t1\n" } else { "" };
        f.write(
            "score.tsv",
            &format!(
                "#pgs_id=VALIDATION\n#genome_build=GRCh38\n\
                 chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n\
                 {first}chr1\t4\tC\tT\t1\n"
            ),
        );
        f.run(&["compile", "score.tsv", "--reference", "ref.fa", "--out", "packs"]);
        f
    }

    fn write(&self, name: &str, text: &str) {
        std::fs::write(self.dir.join(name), text).unwrap();
    }

    fn run(&self, args: &[&str]) {
        let output = Command::new(env!("CARGO_BIN_EXE_pgsum"))
            .current_dir(&self.dir)
            .args(args)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{args:?}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    fn vcf(names: &[String], first: &[&str], second: &[&str]) -> String {
        let fields = |gts: &[&str]| {
            gts.iter()
                .map(|gt| format!("{gt}:30:60"))
                .collect::<Vec<_>>()
                .join("\t")
        };
        format!(
            "##fileformat=VCFv4.2\n\
             ##contig=<ID=chr1,length=4>\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
             ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n\
             chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t{}\n\
             chr1\t4\t.\tT\tC\t.\tPASS\t.\tGT:DP:GQ\t{}\n",
            names.join("\t"),
            fields(first),
            fields(second)
        )
    }

    fn panel(&self, called: usize, varying_first: bool) {
        let names: Vec<_> = (0..100).map(|i| format!("R{i}")).collect();
        let first: Vec<_> = (0..100)
            .map(|i| {
                if varying_first {
                    ["0/0", "0/1", "1/1"][i % 3]
                } else {
                    "0/1"
                }
            })
            .collect();
        let second: Vec<_> = (0..100).map(|i| if i < called { "0/1" } else { "./." }).collect();
        self.write("panel.vcf", &Self::vcf(&names, &first, &second));
        self.write(
            "groups.tsv",
            &format!(
                "sample\tsuper_pop\n{}",
                names.iter().map(|n| format!("{n}\tX\n")).collect::<String>()
            ),
        );
        self.run(&[
            "extract",
            "--gvcf",
            "panel.vcf",
            "--all-samples",
            "--reference",
            "ref.fa",
            "--pack",
            "packs",
            "--out",
            "panel.pgsc",
        ]);
    }

    fn score(&self, second_gt: &str, bundle: bool) -> Value {
        self.write("sample.vcf", &Self::vcf(&["S".into()], &["0/1"], &[second_gt]));
        self.run(&[
            "extract",
            "--gvcf",
            "sample.vcf",
            "--reference",
            "ref.fa",
            "--pack",
            "packs",
            "--out",
            "sample.pgsg",
        ]);
        let mut args = vec![
            "score",
            "--genotypes",
            "sample.pgsg",
            "--pack",
            "packs",
            "--reference-panel",
            "panel.pgsc",
            "--reference-groups",
            "groups.tsv",
            "--out",
            "results",
        ];
        if bundle {
            args.push("--bundle");
        }
        self.run(&args);
        let bytes = if bundle {
            zstd::decode_all(std::fs::File::open(self.dir.join("results/results.jsonl.zst")).unwrap()).unwrap()
        } else {
            std::fs::read(self.dir.join("results/VALIDATION.score.json")).unwrap()
        };
        serde_json::from_slice(&bytes).unwrap()
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

#[test]
fn low_call_rate_terms_cannot_change_percentiles_or_inflate_coverage() {
    let f = Fixture::new("threshold", true);
    for called in [100, 99, 98, 0] {
        f.panel(called, true);
        let r = f.score("1/1", false);
        let p = &r["reference"];
        assert_eq!(r["schema"], "pgsum-score-v3");
        assert_eq!(r["partial"]["raw_score"], "3");
        let matched = if called >= 99 { 2 } else { 1 };
        assert_eq!(p["matched_terms"], matched);
        assert_eq!(p["excluded_low_call_rate_terms"], 2 - matched);
        assert_eq!(p["matched_term_fraction"], matched as f64 / 2.0);
        assert_eq!(p["matched_weight_fraction"], matched as f64 / 2.0);
        assert_eq!(p["meets_coverage_guideline"], called >= 99);
        // Position 1 has 34 zero, 33 one and 33 two dosages. Retained position 4 adds one to each
        // panel score and two to ours: 100 * (67 + 33/2) / 100 = 83.5. Without it: (34 + 33/2) = 50.5.
        assert_eq!(p["score"], if called >= 99 { "3" } else { "1" });
        assert_eq!(p["percentile_all"], if called >= 99 { 83.5 } else { 50.5 });
        assert_eq!(p["fills"]["genotypes_filled"], if called == 99 { 1 } else { 0 });
        assert_eq!(r["imputation_performed"], called == 99);
        assert_eq!(r["imputation"]["reference_panel"], called == 99);
        assert_eq!(r["imputation"]["sample_fill"], false);
        assert_eq!(r["imputation"]["observed_scores"], false);
        assert_eq!(f.score("1/1", true), r, "bundle metadata must match loose JSON");
        if called < 99 {
            let changed = f.score("0/0", false);
            assert_eq!(changed["partial"]["raw_score"], "1");
            assert_eq!(changed["reference"], *p, "excluded genotype must not move placement");
        }
    }
}

#[test]
fn no_matched_terms_withholds_percentiles() {
    let f = Fixture::new("unmatched", false);
    f.panel(98, false);
    let r = f.score("1/1", false);
    let p = &r["reference"];
    assert_eq!(p["matched_terms"], 0);
    assert_eq!(p["excluded_low_call_rate_terms"], 1);
    assert_eq!(p["matched_weight_fraction"], 0.0);
    assert_eq!(p["meets_coverage_guideline"], false);
    assert!(p.get("percentile_all").is_none());
    assert!(p.get("nearest_group_percentile").is_none());
    assert_eq!(p["groups"], serde_json::json!([]));
    assert_eq!(r["imputation_performed"], false);
}

#[test]
fn zero_variance_has_no_z_score() {
    let f = Fixture::new("constant", true);
    f.panel(98, false);
    let r = f.score("1/1", false);
    let group = &r["reference"]["groups"][0];
    assert_eq!(r["reference"]["score"], "1");
    assert_eq!(group["percentile"], 50.0);
    assert_eq!(group["sd"], 0.0);
    assert!(group.get("z").is_some_and(Value::is_null));
}

#[test]
fn zero_frequency_fill_is_reported_even_when_the_filled_score_is_withheld() {
    let f = Fixture::new("sample-fill", true);
    f.panel(100, true);
    let observed = f.score("./.", false);
    assert_eq!(observed["imputation_performed"], false);
    let mut frequencies = flate2::write::GzEncoder::new(
        std::fs::File::create(f.dir.join("frequencies.vcf.gz")).unwrap(),
        flate2::Compression::fast(),
    );
    frequencies
        .write_all(
            b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
          chr1\t4\t.\tT\tC\t.\tPASS\tAF=0\n",
        )
        .unwrap();
    frequencies.finish().unwrap();
    f.run(&[
        "frequencies",
        "--vcf",
        "frequencies.vcf.gz",
        "--field",
        "AF",
        "--out",
        "frequencies.pgsf",
    ]);
    f.run(&[
        "score",
        "--genotypes",
        "sample.pgsg",
        "--pack",
        "packs",
        "--fill-frequencies",
        "frequencies.pgsf",
        "--out",
        "filled",
    ]);
    let r: Value = serde_json::from_slice(&std::fs::read(f.dir.join("filled/VALIDATION.score.json")).unwrap()).unwrap();
    assert_eq!(
        r["partial"], observed["partial"],
        "filling must not change observed scores"
    );
    assert_eq!(r["fill"]["filled"]["terms"], 1);
    assert_eq!(r["fill"]["filled"]["sum"], "0");
    assert!(
        r["fill"]["filled_score"].is_null(),
        "one of two terms is below the coverage gate"
    );
    assert_eq!(r["imputation_performed"], true);
    assert_eq!(r["imputation"]["sample_fill"], true);
    assert_eq!(r["imputation"]["observed_scores"], false);
    assert_eq!(r["imputation"]["reference_panel"], false);
}
