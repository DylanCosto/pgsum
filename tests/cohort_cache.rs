//! Eviction must release retained blocks without changing exact scores or invalidating live row handles.
use pgsum::cohort::{CohortOptions, CohortTable};
use std::io::Write;

#[test]
fn eviction_preserves_live_rows_and_scores_across_readers() {
    let root = std::env::temp_dir().join(format!("pgsum-cohort-cache-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let n = pgsum::cohort::BLOCK_TARGETS * 2 + 3;
    let reference_path = root.join("ref.fa");
    std::fs::write(&reference_path, format!(">1\n{}\n", "A".repeat(n))).unwrap();
    std::fs::write(root.join("ref.fa.fai"), format!("1\t{n}\t3\t{n}\t{}\n", n + 1)).unwrap();
    let reference = pgsum::reference::Reference::open(&reference_path).unwrap();
    let identity = pgsum::compile::reference_identity(&reference).unwrap();
    let scoring = root.join("CACHE.tsv");
    let vcf = root.join("cohort.vcf");
    let mut score = std::io::BufWriter::new(std::fs::File::create(&scoring).unwrap());
    let mut source = std::io::BufWriter::new(std::fs::File::create(&vcf).unwrap());
    writeln!(
        score,
        "#pgs_id=CACHE\n#genome_build=GRCh38\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight"
    )
    .unwrap();
    writeln!(
        source,
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD"
    )
    .unwrap();
    let mut total = 0i64;
    for i in 0..n {
        let (weight, units) = [("0.1", 100), ("0.002", 2), ("-0.03", -30)][i % 3];
        total += units;
        writeln!(score, "1\t{}\tA\tG\t{weight}", i + 1).unwrap();
        writeln!(source, "1\t{}\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t./.", i + 1).unwrap();
    }
    drop(score);
    drop(source);
    pgsum::compile::compile_file(&scoring, None, &reference, &identity, None, &root).unwrap();
    let pack_path = root.join("CACHE.pgsp");
    let path = root.join("cohort.pgsc");
    pgsum::cohort::extract_cohort(
        &vcf,
        &reference,
        &identity,
        std::slice::from_ref(&pack_path),
        &CohortOptions {
            accept_missing_quality: true,
            threads: 2,
            ..Default::default()
        },
        &path,
    )
    .unwrap();
    assert!(CohortTable::open_with_cache_bytes(&path, 0).is_err());
    let small = CohortTable::open_with_cache_bytes(&path, 1).unwrap();
    assert_eq!(small.header.blocks.len(), 3);
    let mut reader = small.reader();
    let first = small
        .row(pgsum::genotypes::target_key(1, 1, b'A', b'G'))
        .unwrap()
        .unwrap();
    let last = reader
        .row(pgsum::genotypes::target_key(1, n as u32, b'A', b'G'))
        .unwrap()
        .unwrap();
    assert_eq!(&*first, &*last);
    assert_eq!(pgsum::cohort::code(&first, 1), 1);
    assert_eq!(small.cache_stats().unwrap().resident_blocks, 1);
    assert_eq!(small.cache_stats().unwrap().evictions, 1);
    drop((first, last));
    drop(reader);
    let large = CohortTable::open(&path).unwrap();
    let pack = pgsum::pack::Pack::open(&pack_path).unwrap();
    let full = pgsum::cohort::score_cohort(&pack, &large, &Default::default()).unwrap();
    for (sample, dosage) in [2, 1, 0, 0].iter().enumerate() {
        let expected = pgsum::decimal::Decimal::parse(&format!("{}e-3", total * dosage)).unwrap();
        let actual = pgsum::decimal::Decimal::parse(&full.text(sample)).unwrap();
        assert_eq!(actual.numeric_cmp(&expected), std::cmp::Ordering::Equal);
    }
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..3)
            .map(|_| scope.spawn(|| pgsum::cohort::score_cohort(&pack, &small, &Default::default()).unwrap()))
            .collect();
        for handle in handles {
            let scored = handle.join().unwrap();
            for s in 0..4 {
                assert_eq!(scored.text(s), full.text(s));
            }
            assert_eq!(scored.scorable, full.scorable);
            assert_eq!(scored.effect_scorable, full.effect_scorable);
        }
    });
    let stats = small.cache_stats().unwrap();
    assert!(stats.evictions >= 3);
    assert!(stats.resident_bytes <= stats.largest_block_bytes.max(stats.capacity_bytes));
    assert_eq!(stats.resident_blocks, 1);
    let mut bounded_report = Vec::new();
    let mut full_report = Vec::new();
    pgsum::missingness::write_cohort(&pack, &small, &Default::default(), &full, &mut bounded_report).unwrap();
    pgsum::missingness::write_cohort(&pack, &large, &Default::default(), &full, &mut full_report).unwrap();
    assert_eq!(bounded_report, full_report);
    let out = root.join("cli");
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_pgsum"))
        .arg("score")
        .arg("--genotypes")
        .arg(&path)
        .arg("--pack")
        .arg(&pack_path)
        .arg("--out")
        .arg(&out)
        .args(["--cohort-cache-mib", "1"])
        .output()
        .unwrap();
    assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
    let metadata: serde_json::Value =
        serde_json::from_slice(&std::fs::read(out.join("cohort-scores.metadata.json")).unwrap()).unwrap();
    assert_eq!(metadata["cohort_cache"]["stats"]["capacity_bytes"], 1024 * 1024);
    assert!(metadata["cohort_cache"]["stats"]["evictions"].as_u64().unwrap() >= 1);
    std::fs::remove_dir_all(root).unwrap();
}

// Split a genuine extraction into smaller, valid v3 frames so cache/report regressions span many
// blocks without millions of probabilities. Retain all original calls, metadata and source references.
fn small_frames(source: &std::path::Path, dest: &std::path::Path) {
    use sha2::{Digest, Sha256};
    fn hash(bytes: &[u8]) -> String {
        format!(
            "sha256:{}",
            Sha256::digest(bytes)
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect::<String>()
        )
    }
    fn append(out: &mut Vec<u8>, bytes: &[u8]) -> pgsum::genotypes::Frame {
        let compressed = zstd::bulk::compress(bytes, 3).unwrap();
        let frame = pgsum::genotypes::Frame {
            offset: (out.len() - 8) as u64,
            compressed_bytes: compressed.len() as u64,
            sha256: hash(bytes),
        };
        out.extend(compressed);
        frame
    }
    let table = CohortTable::open(source).unwrap();
    let mut header = table.header.clone();
    assert_eq!(header.blocks.len(), 1);
    let bytes = std::fs::read(source).unwrap();
    let end = bytes.len();
    let header_len = u64::from_le_bytes(bytes[end - 16..end - 8].try_into().unwrap()) as usize;
    let mut out = bytes[..end - 16 - header_len].to_vec();
    let frame = &header.blocks[0].frame;
    let raw = zstd::decode_all(&bytes[8 + frame.offset as usize..8 + (frame.offset + frame.compressed_bytes) as usize])
        .unwrap();
    let n = u64::from_le_bytes(raw[..8].try_into().unwrap()) as usize;
    let width = pgsum::cohort::row_bytes(header.samples.len());
    let codes = 8 + n * 8;
    let tail = codes + n * width;
    let measurements: Vec<(usize, serde_json::Value)> = serde_json::from_slice(&raw[tail..]).unwrap();
    let frame = &header.diagnostics[0];
    let data =
        zstd::decode_all(&bytes[8 + frame.offset as usize..8 + (frame.offset + frame.compressed_bytes) as usize])
            .unwrap();
    let diagnostics: Vec<serde_json::Value> = serde_json::from_slice(&data).unwrap();
    header.blocks.clear();
    header.diagnostics.clear();
    for (i, diagnostic) in diagnostics.iter().enumerate() {
        let key = u64::from_le_bytes(raw[8 + i * 8..16 + i * 8].try_into().unwrap());
        let mut packed = 1u64.to_le_bytes().to_vec();
        packed.extend(key.to_le_bytes());
        packed.extend(&raw[codes + i * width..codes + (i + 1) * width]);
        let measured: Vec<_> = measurements
            .iter()
            .filter(|(r, _)| *r == i)
            .map(|(_, v)| (0, v))
            .collect();
        packed.extend(serde_json::to_vec(&measured).unwrap());
        let frame = append(&mut out, &packed);
        header.blocks.push(pgsum::genotypes::BlockInfo {
            first_key: key,
            last_key: key,
            targets: 1,
            frame,
        });
        header
            .diagnostics
            .push(append(&mut out, &serde_json::to_vec(&[diagnostic]).unwrap()));
    }
    let digests = std::iter::once(&header.sequence_frame.sha256)
        .chain(header.blocks.iter().map(|b| &b.frame.sha256))
        .chain(header.diagnostics.iter().map(|d| &d.sha256))
        .map(|s| format!("{s}\n"))
        .collect::<String>();
    header.body_sha256 = hash(digests.as_bytes());
    let header = serde_json::to_vec(&header).unwrap();
    out.extend(&header);
    out.extend((header.len() as u64).to_le_bytes());
    out.extend(b"PGSUMCOE");
    std::fs::write(dest, out).unwrap();
}

#[test]
fn quantitative_shared_pass_and_unordered_fallback_preserve_scores_and_bound_report_reads() {
    use pgsum::dosage::Field;
    use pgsum::pack::Pack;
    let root = std::env::temp_dir().join(format!("pgsum-cohort-shared-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let reference_path = root.join("ref.fa");
    std::fs::write(&reference_path, format!(">1\n{}\n", "A".repeat(24))).unwrap();
    std::fs::write(root.join("ref.fa.fai"), "1\t24\t3\t24\t25\n").unwrap();
    let reference = pgsum::reference::Reference::open(&reference_path).unwrap();
    let identity = pgsum::compile::reference_identity(&reference).unwrap();
    let mut paths = Vec::new();
    for (name, weight, dominant, reverse) in [
        ("ORDERED", "0.1", false, false),
        ("REVERSE", "-0.03", false, true),
        ("DOMINANT", "0.002", true, false),
    ] {
        let score_path = root.join(format!("{name}.tsv"));
        let mut text = format!(
            "#pgs_id={name}\n#genome_build=GRCh38\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\tis_dominant\n"
        );
        let positions: Vec<_> = if reverse {
            (1..=24).rev().collect()
        } else {
            (1..=24).collect()
        };
        if dominant {
            text.push_str("1\t\tA\tG\t1.234567890123\tTrue\n");
        }
        for pos in positions {
            text.push_str(&format!(
                "1\t{pos}\tA\tG\t{weight}\t{}\n",
                if dominant { "True" } else { "False" }
            ));
        }
        if dominant {
            text.push_str("1\t\tA\tG\t-9.87654\tTrue\n");
        }
        std::fs::write(&score_path, text).unwrap();
        pgsum::compile::compile_file(&score_path, None, &reference, &identity, None, &root).unwrap();
        paths.push(root.join(format!("{name}.pgsp")));
    }
    let input = root.join("source.vcf");
    let mut text =
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\n".to_owned();
    for pos in 1..=24 {
        text.push_str(&format!(
            "1\t{pos}\t.\tA\tG\t.\tPASS\t.\tGT:DS:GP\t0/1:1:0.25,0.5,0.25\t1/1:2:0,0,1\t0/1:1:0,1,0\t./.:.:.\n"
        ));
    }
    std::fs::write(&input, text).unwrap();
    let packs = paths.iter().map(|p| Pack::open(p).unwrap()).collect::<Vec<_>>();
    for field in [Field::Ds, Field::Gp] {
        let original = root.join(format!("original-{}.pgsc", field.label()));
        pgsum::cohort::extract_cohort(
            &input,
            &reference,
            &identity,
            &paths,
            &CohortOptions {
                dosage_field: field,
                accept_missing_quality: true,
                threads: 2,
                ..Default::default()
            },
            &original,
        )
        .unwrap();
        let path = root.join(format!("split-{}.pgsc", field.label()));
        small_frames(&original, &path);
        let wrong = root.join(format!("wrong-{}.pgsc", field.label()));
        let mut bad_header = CohortTable::open(&path).unwrap().header;
        bad_header.policy.dosage_field = if field == Field::Gp { Field::Ds } else { Field::Gp };
        let mut bytes = std::fs::read(&path).unwrap();
        let n = bytes.len();
        let len = u64::from_le_bytes(bytes[n - 16..n - 8].try_into().unwrap()) as usize;
        bytes.truncate(n - 16 - len);
        let json = serde_json::to_vec(&bad_header).unwrap();
        bytes.extend(&json);
        bytes.extend((json.len() as u64).to_le_bytes());
        bytes.extend(b"PGSUMCOE");
        std::fs::write(&wrong, bytes).unwrap();
        let bad = CohortTable::open(&wrong).unwrap();
        for _ in 0..2 {
            assert!(bad.row(pgsum::genotypes::target_key(1, 1, b'A', b'G')).is_err());
        }
        assert_eq!(bad.cache_stats().unwrap().loads, 0);

        for (threads, capacity) in [(1, 1), (4, 1), (1, 65536), (4, 65536)] {
            let table = CohortTable::open_with_cache_bytes(&path, capacity).unwrap();
            let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
            let batch = pool.install(|| pgsum::cohort::score_cohorts(&packs, &table, &Default::default()).unwrap());
            assert_eq!((batch.shared_packs, batch.independent_packs), (2, 1));
            assert_eq!(table.cache_stats().unwrap().loads, 48); // 24 shared once, 24 for reverse-order fallback.
            assert_eq!(
                table.cache_stats().unwrap().streaming_peak_blocks,
                if capacity == 1 { 1 } else { threads }
            );
            let expected = [
                ["2.4", "0", "2.4", "0"],
                ["-0.72", "0", "-0.72", "0"],
                if field == Field::Gp {
                    ["0.036", "0", "0.048", "0"]
                } else {
                    ["0", "0", "0", "0"]
                },
            ];
            for (i, scored) in batch.scores.iter().enumerate() {
                let independent = pgsum::cohort::score_cohort(&packs[i], &table, &Default::default()).unwrap();
                for sample in 0..4 {
                    let actual = pgsum::decimal::Decimal::parse(&scored.text(sample)).unwrap();
                    let expected = pgsum::decimal::Decimal::parse(expected[i][sample]).unwrap();
                    assert_eq!(actual.numeric_cmp(&expected), std::cmp::Ordering::Equal);
                    assert_eq!(scored.text(sample), independent.text(sample));
                }
                assert_eq!(scored.total_terms, if i == 2 { 26 } else { 24 });
                assert_eq!(scored.effect_all, independent.effect_all);
                assert_eq!(scored.expected_genotype_terms, independent.expected_genotype_terms);
                assert_eq!(scored.scorable, independent.scorable);
                assert_eq!(scored.effect_scorable, independent.effect_scorable);
                let report_table = CohortTable::open_with_cache_bytes(&path, 1).unwrap();
                let mut report = Vec::new();
                pgsum::missingness::write_cohort(&packs[i], &report_table, &Default::default(), scored, &mut report)
                    .unwrap();
                assert_eq!(String::from_utf8(report).unwrap().lines().count(), 4);
                assert!(
                    report_table.cache_stats().unwrap().loads <= 24,
                    "previously validated source evidence must not re-decode probability blocks"
                );
            }
        }
    }
    std::fs::remove_dir_all(root).unwrap();
}
