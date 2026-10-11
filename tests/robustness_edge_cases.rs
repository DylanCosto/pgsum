//! Regression checks for missing-only frequency lookup and canonical indel masks.
//! Expected scores and regional deltas are hand calculated from two diploid calls.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use serde_json::{Value, json};

static NEXT_DIRECTORY: AtomicU64 = AtomicU64::new(0);

struct Fixture {
    root: PathBuf,
    reference: PathBuf,
}

fn run(args: &[&str]) {
    let result = Command::new(env!("CARGO_BIN_EXE_pgsum")).args(args).output().unwrap();
    assert!(
        result.status.success(),
        "{args:?}\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
}

fn text(path: &Path) -> &str {
    path.to_str().unwrap()
}

impl Fixture {
    fn new() -> Self {
        let root = std::env::temp_dir().join(format!(
            "pgsum-robustness-edges-{}-{}",
            std::process::id(),
            NEXT_DIRECTORY.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&root).unwrap();
        let reference = root.join("reference.fa");
        std::fs::write(&reference, ">chr1\nGACGATGCACTAGCATCGTA\n").unwrap();
        std::fs::write(root.join("reference.fa.fai"), "chr1\t20\t6\t20\t21\n").unwrap();
        Self { root, reference }
    }

    fn file(&self, name: &str, contents: &str) -> PathBuf {
        let path = self.root.join(name);
        std::fs::write(&path, contents).unwrap();
        path
    }

    fn gzip(&self, name: &str, contents: &str) -> PathBuf {
        let path = self.root.join(name);
        let mut out =
            flate2::write::GzEncoder::new(std::fs::File::create(&path).unwrap(), flate2::Compression::default());
        out.write_all(contents.as_bytes()).unwrap();
        out.finish().unwrap();
        path
    }

    fn compile(&self, id: &str, pos: u32, effect: &str, other: &str) -> PathBuf {
        let scoring = self.file(
            &format!("{id}.tsv"),
            &format!(
                "#pgs_id={id}\n#genome_build=GRCh38\n#variants_number=1\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n1\t{pos}\t{effect}\t{other}\t1\n"
            ),
        );
        run(&[
            "compile",
            text(&scoring),
            "--reference",
            text(&self.reference),
            "--out",
            text(&self.root.join("packs")),
        ]);
        self.root.join("packs").join(format!("{id}.pgsp"))
    }

    fn extract(&self, alt: &str) {
        let vcf = self.file(
            "cohort.vcf",
            &format!(
                "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n1\t3\t.\tC\t{alt}\t.\tPASS\t.\tGT\t0/1\t1/1\n"
            ),
        );
        run(&[
            "extract",
            "--gvcf",
            text(&vcf),
            "--reference",
            text(&self.reference),
            "--pack",
            text(&self.root.join("packs")),
            "--all-samples",
            "--accept-missing-quality",
            "--out",
            text(&self.root.join("cohort.pgsc")),
        ]);
    }

    fn experiment(&self, manifest: Value) -> PathBuf {
        let path = self.file("experiment.json", &manifest.to_string());
        let out = self.root.join("results");
        run(&[
            "robustness",
            "--experiment",
            text(&path),
            "--out",
            text(&out),
            "--terms",
        ]);
        assert!(out.join("report.json").exists());
        out
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.root);
    }
}

fn rows(path: &Path) -> Vec<BTreeMap<String, String>> {
    let data = std::fs::read_to_string(path).unwrap();
    let mut lines = data.lines();
    let header: Vec<_> = lines.next().unwrap().split('\t').collect();
    lines
        .map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            assert_eq!(fields.len(), header.len());
            header
                .iter()
                .zip(fields)
                .map(|(key, value)| ((*key).to_owned(), value.to_owned()))
                .collect()
        })
        .collect()
}

#[test]
fn frequency_records_are_not_consulted_for_retained_observations() {
    let fixture = Fixture::new();
    fixture.compile("RETAINED", 3, "T", "C");
    fixture.extract("T");
    // Such an unsplit record deliberately cannot be used by the frequency-fill
    // policy. Every call is retained, so its unusable AF must not affect scoring.
    let frequency_vcf = fixture.gzip(
        "frequencies.vcf.gz",
        "##fileformat=VCFv4.3\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"ALT frequencies\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t3\t.\tC\tT,G\t.\tPASS\tAF=0.2,0.3\n",
    );
    run(&[
        "frequencies",
        "--vcf",
        text(&frequency_vcf),
        "--field",
        "AF",
        "--out",
        text(&fixture.root.join("frequencies.pgsf")),
    ]);
    let out = fixture.experiment(json!({
        "schema": "pgsum-robustness-experiment-v1",
        "genotypes": "cohort.pgsc",
        "packs": ["packs/RETAINED.pgsp"],
        "scenarios": [
            {"id": "control"},
            {
                "id": "unused_fill",
                "masks": [{"kind": "random", "rate": 0.0, "seed": 1, "unit": "site"}],
                "fill_frequencies": "frequencies.pgsf",
                "frequency_source_note": "Independent synthetic frequencies; no evaluation samples."
            }
        ]
    }));
    let samples = rows(&out.join("samples.tsv"));
    assert_eq!(samples.len(), 4);
    for row in samples {
        let expected = if row["sample"] == "A" { "1" } else { "2" };
        assert_eq!(row["baseline_score"], expected);
        assert_eq!(row["scenario_score"], expected);
        assert_eq!(row["delta"], "0");
        assert_eq!(row["retained_terms"], "1");
        assert_eq!(row["masked_terms"], "0");
        assert_eq!(row["filled_terms"], "0");
    }
}

#[test]
fn normalized_indel_anchor_controls_bed_random_masks_and_region_attribution() {
    let fixture = Fixture::new();
    // The padded change at 2, AC>ATT, and the minimal change at 3, C>TT,
    // normalize to exactly the same indel. Both effect alleles are the ALT.
    let padded = fixture.compile("PADDED", 2, "ATT", "AC");
    let normalized = fixture.compile("NORMALIZED", 3, "TT", "C");
    for path in [&padded, &normalized] {
        let pack = pgsum::pack::Pack::open(path).unwrap();
        let term = pack.terms().next().unwrap().unwrap();
        let variant = term.inferred_sequence.expect("inferred normalized indel");
        assert_eq!(
            (variant.pos, variant.ref_allele.as_str(), variant.alt.as_str()),
            (3, "C", "TT")
        );
    }
    fixture.extract("TT");
    fixture.file("anchor.bed", "chr1\t2\t3\n");
    let out = fixture.experiment(json!({
        "schema": "pgsum-robustness-experiment-v1",
        "genotypes": "cohort.pgsc",
        "packs": ["packs/PADDED.pgsp", "packs/NORMALIZED.pgsp"],
        "scoring": {"allow_inferred_indels": true},
        "regions": [{"id": "normalized_anchor", "path": "anchor.bed"}],
        "scenarios": [
            {"id": "bed_drop", "masks": [{"kind": "bed", "mode": "drop", "path": "anchor.bed"}]},
            {"id": "bed_keep", "masks": [{"kind": "bed", "mode": "keep", "path": "anchor.bed"}]},
            // Independently computed SHA-256 draws: seed 2 at rate 0.5 hides
            // canonical position 3 and retains padded source position 2.
            {"id": "random", "masks": [{"kind": "random", "rate": 0.5, "seed": 2, "unit": "site"}]}
        ]
    }));
    let samples = rows(&out.join("samples.tsv"));
    assert_eq!(samples.len(), 12);
    for row in samples {
        let baseline = if row["sample"] == "A" { "1" } else { "2" };
        let kept = row["scenario"] == "bed_keep";
        assert_eq!(row["baseline_score"], baseline);
        assert_eq!(row["scenario_score"], if kept { baseline } else { "0" }, "{row:?}");
        assert_eq!(row["masked_terms"], if kept { "0" } else { "1" }, "{row:?}");
        assert_eq!(row["delta"], if kept { "0".into() } else { format!("-{baseline}") });
    }
    let regions = rows(&out.join("regions.tsv"));
    assert_eq!(regions.len(), 12);
    for row in regions {
        let expected = if row["scenario"] == "bed_keep" {
            "0"
        } else if row["sample"] == "A" {
            "-1"
        } else {
            "-2"
        };
        assert_eq!(row["delta"], expected, "{row:?}");
    }
}
