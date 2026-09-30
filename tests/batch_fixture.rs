//! Resume must validate content and complete outputs, not merely the presence of a directory.
use serde_json::Value;
use std::path::{Path, PathBuf};
use std::process::Command;

fn cli(args: &[&str], success: bool) -> std::process::Output {
    let out = Command::new(env!("CARGO_BIN_EXE_pgsum")).args(args).output().unwrap();
    assert_eq!(
        out.status.success(),
        success,
        "{args:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out
}
fn s(path: &Path) -> &str {
    path.to_str().unwrap()
}
fn read(path: &Path) -> Value {
    serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap()
}
fn input(path: &Path, gt: &str) {
    std::fs::write(path, format!("##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tINPUT\n1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DS\t{gt}:0.125\n")).unwrap();
}
fn setup(name: &str) -> (PathBuf, PathBuf, PathBuf) {
    let root = std::env::temp_dir().join(format!("pgsum-batch-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/synthetic.fa");
    let reference = root.join("reference.fa");
    std::fs::copy(&source, &reference).unwrap();
    std::fs::copy(source.with_extension("fa.fai"), reference.with_extension("fa.fai")).unwrap();
    let score = root.join("BATCH.tsv");
    std::fs::write(&score, "#pgs_id=BATCH\n#genome_build=GRCh38\n#variants_number=1\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n1\t2\tT\tC\t0.2\n").unwrap();
    cli(
        &["compile", s(&score), "--reference", s(&reference), "--out", s(&root)],
        true,
    );
    input(&root.join("a.vcf"), "0/1");
    input(&root.join("b.vcf"), "1/1");
    let sheet = root.join("samples.tsv");
    std::fs::write(&sheet, "sample_id\tinput\nA\ta.vcf\nB\tb.vcf\n").unwrap();
    (root, reference, sheet)
}
fn batch(root: &Path, reference: &Path, sheet: &Path, extra: &[&str], success: bool) -> Value {
    let out = root.join("out");
    let mut args = vec![
        "batch",
        "--samplesheet",
        s(sheet),
        "--reference",
        s(reference),
        "--pack",
        s(root),
        "--out",
        s(&out),
        "--accept-missing-quality",
        "--threads",
        "2",
        "--jobs",
        "2",
    ];
    args.extend_from_slice(extra);
    cli(&args, success);
    read(&out.join("batch-results.json"))
}
fn directory(report: &Value, i: usize) -> PathBuf {
    PathBuf::from(report["samples"][i]["directory"].as_str().unwrap())
}

#[test]
fn resume_invalidation_integrity_and_incremental_samples() {
    let (root, reference, sheet) = setup("resume");
    let first = batch(&root, &reference, &sheet, &[], true);
    assert_eq!(first["samples"][0]["status"], "completed");
    let aggregate = root.join("out/batch-scores.tsv");
    assert_eq!(std::fs::read_to_string(&aggregate).unwrap().lines().count(), 3);
    assert_eq!(
        first["scores"]["sha256"],
        pgsum::digest::file_sha256(&aggregate).unwrap()
    );
    assert_eq!(first["complete"], true);
    assert_eq!(read(&directory(&first, 0).join("BATCH.score.json"))["raw_score"], "0.2");
    assert_eq!(read(&directory(&first, 1).join("BATCH.score.json"))["raw_score"], "0.4");
    let repeat = batch(&root, &reference, &sheet, &["--resume"], true);
    assert!(
        repeat["samples"]
            .as_array()
            .unwrap()
            .iter()
            .all(|r| r["status"] == "reused")
    );
    // Change content while preserving size and mtime: the actual hash must invalidate A only.
    let path = root.join("a.vcf");
    let modified = std::fs::metadata(&path).unwrap().modified().unwrap();
    input(&path, "0/0");
    std::fs::File::options()
        .write(true)
        .open(&path)
        .unwrap()
        .set_times(std::fs::FileTimes::new().set_modified(modified))
        .unwrap();
    let changed = batch(&root, &reference, &sheet, &["--resume"], true);
    assert_eq!(changed["samples"][0]["status"], "completed");
    assert_eq!(changed["samples"][1]["status"], "reused");
    assert_ne!(directory(&changed, 0), directory(&first, 0));
    assert_eq!(
        read(&directory(&changed, 0).join("BATCH.score.json"))["raw_score"],
        "0.0"
    );
    // A damaged output cannot be accepted just because its manifest exists.
    std::fs::write(directory(&changed, 1).join("BATCH.score.json"), "broken").unwrap();
    let repaired = batch(&root, &reference, &sheet, &["--resume"], true);
    assert_eq!(repaired["samples"][1]["status"], "completed");
    assert_eq!(
        read(&directory(&repaired, 1).join("BATCH.score.json"))["raw_score"],
        "0.4"
    );
    // Add a sample: existing valid samples retain their completed generations.
    std::fs::write(&sheet, "sample_id\tinput\nA\ta.vcf\nB\tb.vcf\nC\tb.vcf\n").unwrap();
    let added = batch(&root, &reference, &sheet, &["--resume"], true);
    assert_eq!(added["samples"][0]["status"], "reused");
    assert_eq!(added["samples"][1]["status"], "reused");
    assert_eq!(added["samples"][2]["status"], "completed");
    let settings = batch(&root, &reference, &sheet, &["--resume", "--dosage-field", "ds"], true);
    assert!(
        settings["samples"]
            .as_array()
            .unwrap()
            .iter()
            .all(|r| r["status"] == "completed")
    );
    assert_eq!(
        read(&directory(&settings, 0).join("BATCH.score.json"))["raw_score"],
        "0.0250"
    );
}

#[test]
fn failed_samples_do_not_prevent_other_samples_finishing() {
    let (root, reference, sheet) = setup("failure");
    std::fs::write(
        &sheet,
        "sample_id\tinput\tsample\nA\ta.vcf\tNOT_PRESENT\nB\tb.vcf\tINPUT\n",
    )
    .unwrap();
    let failed = batch(&root, &reference, &sheet, &[], false);
    assert_eq!(failed["samples"][0]["status"], "failed");
    assert_eq!(failed["complete"], false);
    assert_eq!(
        std::fs::read_to_string(root.join("out/batch-scores.tsv"))
            .unwrap()
            .lines()
            .count(),
        2
    );
    assert_eq!(failed["samples"][1]["status"], "completed");
    assert!(!directory(&failed, 0).join("run-manifest.json").exists());
    std::fs::write(&sheet, "sample_id\tinput\tsample\nA\ta.vcf\tINPUT\nB\tb.vcf\tINPUT\n").unwrap();
    let fixed = batch(&root, &reference, &sheet, &["--resume"], true);
    assert_eq!(fixed["samples"][0]["status"], "completed");
    assert_eq!(fixed["samples"][1]["status"], "reused");
}

#[test]
fn sheet_rejects_duplicate_ids_and_path_traversal() {
    let (root, _, sheet) = setup("sheet");
    for contents in [
        "sample_id\tinput\nA\ta.vcf\nA\tb.vcf\n",
        "sample_id\tinput\n../escape\ta.vcf\n",
        "sample_id\tinput\nA\t\n",
        "sample_id\tinput\tinput\nA\ta.vcf\tb.vcf\n",
    ] {
        std::fs::write(&sheet, contents).unwrap();
        assert!(pgsum::batch::read_sheet(&sheet).is_err());
    }
    assert!(!root.parent().unwrap().join("escape").exists());
}

#[test]
fn changed_packs_and_reference_require_new_runs() {
    let (root, reference, sheet) = setup("packs");
    let initial = batch(&root, &reference, &sheet, &[], true);
    let score = root.join("BATCH.tsv");
    let changed = std::fs::read_to_string(&score).unwrap().replace("\t0.2\n", "\t0.3\n");
    std::fs::write(&score, changed).unwrap();
    cli(
        &["compile", s(&score), "--reference", s(&reference), "--out", s(&root)],
        true,
    );
    let updated = batch(&root, &reference, &sheet, &["--resume"], true);
    assert_ne!(directory(&initial, 0), directory(&updated, 0));
    assert_eq!(
        read(&directory(&updated, 0).join("BATCH.score.json"))["raw_score"],
        "0.3"
    );
    let modified = std::fs::metadata(&reference).unwrap().modified().unwrap();
    let content = std::fs::read_to_string(&reference).unwrap().replacen("NN", "NA", 1);
    std::fs::write(&reference, content).unwrap();
    std::fs::File::options()
        .write(true)
        .open(&reference)
        .unwrap()
        .set_times(std::fs::FileTimes::new().set_modified(modified))
        .unwrap();
    let out = root.join("out");
    let rejected = cli(
        &[
            "batch",
            "--samplesheet",
            s(&sheet),
            "--reference",
            s(&reference),
            "--pack",
            s(&root),
            "--out",
            s(&out),
            "--resume",
        ],
        false,
    );
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("reference"));
    cli(
        &["compile", s(&score), "--reference", s(&reference), "--out", s(&root)],
        true,
    );
    let rebuilt = batch(&root, &reference, &sheet, &["--resume"], true);
    assert_ne!(directory(&updated, 0), directory(&rebuilt, 0));
    assert_eq!(rebuilt["samples"][0]["status"], "completed");
}

#[test]
fn output_lock_releases_without_manual_cleanup() {
    let (root, reference, sheet) = setup("lock");
    let out = root.join("out");
    std::fs::create_dir_all(&out).unwrap();
    let lock = std::fs::File::create(out.join(".batch.lock")).unwrap();
    fs2::FileExt::try_lock_exclusive(&lock).unwrap();
    let result = cli(
        &[
            "batch",
            "--samplesheet",
            s(&sheet),
            "--reference",
            s(&reference),
            "--pack",
            s(&root),
            "--out",
            s(&out),
        ],
        false,
    );
    assert!(String::from_utf8_lossy(&result.stderr).contains("another batch owns"));
    drop(lock);
    let completed = batch(&root, &reference, &sheet, &[], true);
    let manifest = directory(&completed, 0).join("run-manifest.json");
    // Simulate an incomplete publication: no completion manifest means no resumable result.
    std::fs::remove_file(&manifest).unwrap();
    let repaired = batch(&root, &reference, &sheet, &["--resume"], true);
    assert_eq!(repaired["samples"][0]["status"], "completed");
    assert_eq!(repaired["samples"][1]["status"], "reused");
}
