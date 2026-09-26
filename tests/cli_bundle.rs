//! `pgsum run --bundle` writes the same results as one file per score, in one compressed JSON-lines file.

use std::path::Path;
use std::process::Command;

fn pgsum(args: &[&str]) {
    let status = Command::new(env!("CARGO_BIN_EXE_pgsum")).args(args).status().unwrap();
    assert!(status.success(), "pgsum {args:?}");
}

#[test]
fn bundle_matches_per_score_files() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let dir = std::env::temp_dir().join(format!("pgsum-bundle-test-{}", std::process::id()));
    let (packs, loose, bundled) = (dir.join("packs"), dir.join("loose"), dir.join("bundled"));
    std::fs::create_dir_all(&packs).unwrap();
    let f = |name: &str| fixtures.join(name).to_str().unwrap().to_owned();
    pgsum(&[
        "compile",
        &f("PGS999997_hmPOS_GRCh38.txt.gz"),
        &f("PGS999998_hmPOS_GRCh38.txt.gz"),
        "--reference",
        &f("synthetic.fa"),
        "--out",
        packs.to_str().unwrap(),
    ]);
    let (gvcf, fasta) = (f("synthetic.g.vcf.gz"), f("synthetic.fa"));
    for (out, extra) in [(&loose, None), (&bundled, Some("--bundle"))] {
        let mut args = vec![
            "run",
            "--gvcf",
            &gvcf,
            "--reference",
            &fasta,
            "--pack",
            packs.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
        ];
        args.extend(extra);
        pgsum(&args);
    }

    let text =
        String::from_utf8(zstd::decode_all(std::fs::File::open(bundled.join("results.jsonl.zst")).unwrap()).unwrap())
            .unwrap();
    let lines: Vec<serde_json::Value> = text.lines().map(|l| serde_json::from_str(l).unwrap()).collect();
    let ids: Vec<&str> = lines.iter().map(|v| v["pgs_id"].as_str().unwrap()).collect();
    assert_eq!(ids, ["PGS999997", "PGS999998"]);
    for v in &lines {
        let id = v["pgs_id"].as_str().unwrap();
        let file: serde_json::Value =
            serde_json::from_slice(&std::fs::read(loose.join(format!("{id}.score.json"))).unwrap()).unwrap();
        assert_eq!(v, &file);
        assert!(!bundled.join(format!("{id}.score.json")).exists());
    }
    assert_eq!(
        std::fs::read(loose.join("scores.tsv")).unwrap(),
        std::fs::read(bundled.join("scores.tsv")).unwrap()
    );
    std::fs::remove_dir_all(&dir).unwrap();
}
