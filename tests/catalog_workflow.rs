//! Exercise real HTTP downloads, compilation, verified cache reuse, and CLI execution.
use pgsum::catalog_cache::{Options, resolve_with};
use serde_json::Value;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures").join(name)
}
fn setup(name: &str) -> (PathBuf, PathBuf, Options) {
    let root = std::env::temp_dir().join(format!("pgsum-catalog-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let reference = root.join("reference.fa");
    std::fs::copy(fixture("synthetic.fa"), &reference).unwrap();
    std::fs::copy(fixture("synthetic.fa.fai"), reference.with_extension("fa.fai")).unwrap();
    let options = Options {
        directory: Some(root.join("cache")),
        jobs: 2,
        ..Default::default()
    };
    (root, reference, options)
}
fn record() -> Value {
    serde_json::from_slice(&std::fs::read(fixture("PGS999999.metadata.json")).unwrap()).unwrap()
}
fn serve(bytes: Vec<u8>) -> (String, std::thread::JoinHandle<()>) {
    let socket = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
    let url = format!("http://{}/score.gz", socket.local_addr().unwrap());
    socket.set_nonblocking(true).unwrap();
    let task = std::thread::spawn(move || {
        let start = std::time::Instant::now();
        let mut stream = loop {
            match socket.accept() {
                Ok((stream, _)) => break stream,
                Err(e) if e.kind() == std::io::ErrorKind::WouldBlock && start.elapsed().as_secs() < 10 => {
                    std::thread::sleep(std::time::Duration::from_millis(10))
                }
                Err(e) => panic!("test HTTP server: {e}"),
            }
        };
        stream
            .set_read_timeout(Some(std::time::Duration::from_secs(5)))
            .unwrap();
        let mut request = Vec::new();
        let mut byte = [0];
        while !request.ends_with(b"\r\n\r\n") {
            stream.read_exact(&mut byte).unwrap();
            request.push(byte[0]);
        }
        write!(
            stream,
            "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
            bytes.len()
        )
        .unwrap();
        stream.write_all(&bytes).unwrap();
    });
    (url, task)
}
fn download(reference: &Path, options: &Options) -> Vec<PathBuf> {
    let (url, server) = serve(std::fs::read(fixture("PGS999999_hmPOS_GRCh38.txt.gz")).unwrap());
    let mut r = record();
    r["ftp_harmonized_scoring_files"] = serde_json::json!({"GRCh38": {"positions": url}});
    let paths = resolve_with(
        &["PGS999999".into(), "PGS999999".into()],
        reference,
        options,
        &pgsum::fetch::agent(),
        |ids| {
            assert_eq!(ids, ["PGS999999"]);
            Ok(vec![r])
        },
    )
    .unwrap();
    server.join().unwrap();
    paths
}
fn no_network(reference: &Path, options: &Options) -> pgsum::Result<Vec<PathBuf>> {
    resolve_with(
        &["PGS999999".into()],
        reference,
        options,
        &pgsum::fetch::agent(),
        |_| panic!("unexpected Catalog request"),
    )
}
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
fn json(path: &Path) -> Value {
    serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap()
}

#[test]
fn verified_cache_download_reuse_refresh_and_corruption() {
    let (root, reference, mut options) = setup("cache");
    options.offline = true;
    assert!(
        no_network(&reference, &options)
            .unwrap_err()
            .to_string()
            .contains("offline cache")
    );
    options.offline = false;
    let first = download(&reference, &options);
    assert_eq!(first.len(), 1);
    assert_eq!(no_network(&reference, &options).unwrap(), first);
    options.offline = true;
    assert_eq!(no_network(&reference, &options).unwrap(), first);
    let saved = std::fs::read(&first[0]).unwrap();
    std::fs::write(&first[0], b"damaged").unwrap();
    assert!(no_network(&reference, &options).is_err());
    std::fs::write(&first[0], &saved).unwrap();
    // Even preserved size/timestamp cannot make changed reference bytes reuse old packs.
    let modified = std::fs::metadata(&reference).unwrap().modified().unwrap();
    let original = std::fs::read_to_string(&reference).unwrap();
    std::fs::write(&reference, original.replace("NN", "NA")).unwrap();
    std::fs::File::options()
        .write(true)
        .open(&reference)
        .unwrap()
        .set_times(std::fs::FileTimes::new().set_modified(modified))
        .unwrap();
    assert!(no_network(&reference, &options).is_err());
    std::fs::write(&reference, original).unwrap();
    options.offline = false;
    options.refresh = true;
    let second = download(&reference, &options);
    assert!(first[0].exists());
    assert!(second[0].exists());
    // Metadata's download URL changed; both immutable generations must remain available.
    assert_ne!(first, second);
    options.refresh = false;
    options.offline = true;
    assert_eq!(no_network(&reference, &options).unwrap(), second);
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn catalog_ids_run_preflight_manifests_and_batch() {
    let (root, reference, options) = setup("cli");
    let paths = download(&reference, &options);
    let cache = options.directory.as_ref().unwrap();
    let input = fixture("synthetic.g.vcf.gz");
    let out = root.join("run");
    let args = [
        "run",
        "--ids",
        "PGS999999",
        "--gvcf",
        s(&input),
        "--reference",
        s(&reference),
        "--catalog-cache",
        s(cache),
        "--offline",
        "--out",
        s(&out),
        "--threads",
        "2",
    ];
    let mut check_args = args.to_vec();
    check_args.push("--preflight");
    cli(&check_args, true);
    assert!(out.join("preflight.json").exists());
    assert!(!out.join("genotypes.pgsg").exists());
    cli(&args, true);
    let manifest = json(&out.join("execution-manifest.json"));
    assert_eq!(manifest["schema"], "pgsum-execution-v1");
    assert_eq!(manifest["complete"], true);
    assert_eq!(
        manifest["packs"][0]["sha256"],
        pgsum::digest::file_sha256(&paths[0]).unwrap()
    );
    for (name, digest) in manifest["outputs"].as_object().unwrap() {
        assert_eq!(digest, &pgsum::digest::file_sha256(&out.join(name)).unwrap());
    }
    let explicit = root.join("explicit");
    cli(
        &[
            "run",
            "--pack",
            s(&paths[0]),
            "--gvcf",
            s(&input),
            "--reference",
            s(&reference),
            "--out",
            s(&explicit),
            "--threads",
            "2",
        ],
        true,
    );
    assert_eq!(
        std::fs::read(out.join("scores.tsv")).unwrap(),
        std::fs::read(explicit.join("scores.tsv")).unwrap()
    );
    let sheet = root.join("samples.tsv");
    std::fs::write(&sheet, format!("sample_id\tinput\nA\t{}\n", input.display())).unwrap();
    let batch = root.join("batch");
    cli(
        &[
            "batch",
            "--ids",
            "PGS999999",
            "--samplesheet",
            s(&sheet),
            "--reference",
            s(&reference),
            "--catalog-cache",
            s(cache),
            "--offline",
            "--out",
            s(&batch),
            "--threads",
            "2",
        ],
        true,
    );
    assert_eq!(json(&batch.join("batch-results.json"))["complete"], true);
    // A failed rerun must invalidate the previous completion marker.
    let bad = root.join("bad.vcf");
    std::fs::write(&bad, "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n1\tBAD_POSITION\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\n").unwrap();
    cli(
        &[
            "run",
            "--pack",
            s(&paths[0]),
            "--gvcf",
            s(&bad),
            "--reference",
            s(&reference),
            "--out",
            s(&out),
        ],
        false,
    );
    assert!(!out.join("execution-manifest.json").exists());
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn malformed_download_and_wrong_catalog_ids_never_publish() {
    let (root, reference, mut options) = setup("failure");
    let ids = vec!["PGS999999".into()];
    let agent = pgsum::fetch::agent();
    let error = resolve_with(&ids, &reference, &options, &agent, |_| {
        Ok(vec![serde_json::json!({"id":"PGS999998"})])
    })
    .unwrap_err();
    assert!(error.to_string().contains("response IDs"));
    let (url, server) = serve(b"not a scoring file".to_vec());
    let mut r = record();
    r["ftp_harmonized_scoring_files"] = serde_json::json!({"GRCh38":{"positions":url}});
    assert!(resolve_with(&ids, &reference, &options, &agent, |_| Ok(vec![r])).is_err());
    server.join().unwrap();
    options.offline = true;
    assert!(no_network(&reference, &options).is_err());
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn bad_input_preflight_stops_before_catalog_access() {
    let (root, reference, options) = setup("preflight");
    let cache = options.directory.as_ref().unwrap();
    let input = root.join("bad.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.3\n##contig=<ID=1,length=1234>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n",
    )
    .unwrap();
    let output = cli(
        &[
            "run",
            "--ids",
            "PGS999999",
            "--gvcf",
            s(&input),
            "--reference",
            s(&reference),
            "--catalog-cache",
            s(cache),
            "--out",
            s(&root.join("out")),
            "--preflight",
        ],
        false,
    );
    assert!(String::from_utf8_lossy(&output.stderr).contains("different assembly"));
    assert!(!cache.exists());
    std::fs::remove_dir_all(root).unwrap();
}
