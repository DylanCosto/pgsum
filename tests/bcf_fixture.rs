//! Native BCF and textual VCF must apply the same rules to representable input values.
use noodles_vcf::variant::io::Write as _;
use pgsum::dosage::Field;
use std::path::Path;

fn encode(vcf: &Path, bcf: &Path, compressed: bool) {
    let mut source = noodles_vcf::io::Reader::new(std::io::BufReader::new(std::fs::File::open(vcf).unwrap()));
    let header = source.read_header().unwrap();
    let output: Box<dyn std::io::Write> = if compressed {
        Box::new(noodles_bgzf::io::Writer::new(std::fs::File::create(bcf).unwrap()))
    } else {
        Box::new(std::fs::File::create(bcf).unwrap())
    };
    let mut writer = noodles_bcf::io::Writer::from(output);
    writer.write_header(&header).unwrap();
    for record in source.record_bufs(&header) {
        writer.write_variant_record(&header, &record.unwrap()).unwrap();
    }
    std::io::Write::flush(writer.get_mut()).unwrap();
}

#[test]
fn bcf_streams_hard_calls_dosages_probabilities_and_reference_blocks() {
    let root = std::env::temp_dir().join(format!("pgsum-bcf-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let reference = pgsum::reference::Reference::open(&fixtures.join("synthetic.fa")).unwrap();
    let identity = pgsum::compile::reference_identity(&reference).unwrap();
    let score_path = root.join("BCF_TEST.tsv");
    std::fs::write(&score_path, "#pgs_id=BCF_TEST\n#genome_build=GRCh38\n#variants_number=4\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n1\t1\tA\tG\t1\n1\t2\tT\tC\t0.2\n1\t3\tG\tT\t-0.5\n1\t4\tC\tT\t0.1\n").unwrap();
    pgsum::compile::compile_file(&score_path, None, &reference, &identity, None, &root).unwrap();
    let pack_path = root.join("BCF_TEST.pgsp");
    let pack = pgsum::pack::Pack::open(&pack_path).unwrap();
    let paths = [pack_path];
    let vcf = root.join("source.vcf");
    std::fs::write(
        &vcf,
        concat!(
            "##fileformat=VCFv4.3\n##contig=<ID=1,length=60>\n",
            "##FILTER=<ID=LowQual,Description=\"Low quality\">\n",
            "##FILTER=<ID=RefCall,Description=\"Genotyping model thinks this site is reference.\">\n",
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n",
            "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Block depth\">\n",
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">\n",
            "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Dosage\">\n",
            "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Probabilities\">\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n",
            "1\t1\t.\tA\t<NON_REF>\t.\tRefCall\tEND=1\tGT:MIN_DP:GQ\t0/0:30:60\t0/0:30:60\n",
            "1\t2\t.\tC\tT\t.\tPASS\t.\tGT:DS:GP\t1/1:0.125:0.25,0.25,0.5\t0/1:1.25:0.125,0.5,0.375\n",
            "1\t3\t.\tG\tT\t.\tPASS\t.\tGT:DS:GP\t0/0:1.25:0.5,0.25,0.25\t./.:.:.\n",
            "1\t4\t.\tT\tC\t.\tLowQual\t.\tGT:DS:GP\t0/1:1:0,1,0\t1/1:2:0,0,1\n"
        ),
    )
    .unwrap();
    assert!(!pgsum::bcf::is_bcf(&vcf).unwrap());
    for compressed in [false, true] {
        // A nonstandard suffix verifies content-based detection.
        let bcf = root.join(format!("data-{compressed}.input"));
        encode(&vcf, &bcf, compressed);
        assert!(pgsum::bcf::is_bcf(&bcf).unwrap());
        let checks = pgsum::workflow::preflight(&bcf, &reference.path, Some("A")).unwrap();
        assert_eq!(checks["format"], "BCF");
        assert_eq!(checks["sample"], "A");
        assert_eq!(checks["matching_contig_lengths"], 1);
        assert!(pgsum::workflow::preflight(&bcf, &reference.path, None).is_err());
        for field in [Field::Gt, Field::Ds, Field::Gp] {
            for sample in ["A", "B"] {
                let options = pgsum::extract::Options {
                    dosage_field: field,
                    sample: Some(sample),
                    accept_missing_quality: true,
                    threads: 2,
                    ..Default::default()
                };
                let (text, _) = pgsum::extract::extract(&vcf, &reference, &identity, &paths, &options).unwrap();
                let (binary, _) = pgsum::extract::extract(&bcf, &reference, &identity, &paths, &options).unwrap();
                assert_eq!(text.header.states, binary.header.states);
                assert_eq!(text.header.records_scanned, binary.header.records_scanned);
                assert_eq!(binary.header.gvcf.sha256, pgsum::digest::file_sha256(&bcf).unwrap());
                assert_eq!(binary.header.gvcf.bytes, std::fs::metadata(&bcf).unwrap().len());
                let a = pgsum::score::score(&pack, &text, &Default::default(), None).unwrap();
                let b = pgsum::score::score(&pack, &binary, &Default::default(), None).unwrap();
                assert_eq!(a.partial, b.partial);
                assert_eq!(a.states, b.states);
                assert_eq!(a.expected_genotype_terms, b.expected_genotype_terms);
            }
            let options = pgsum::cohort::CohortOptions {
                dosage_field: field,
                accept_missing_quality: true,
                threads: 2,
                ..Default::default()
            };
            let mut results = Vec::new();
            for (source, name) in [(&vcf, "vcf"), (&bcf, "bcf")] {
                let dest = root.join(format!("{name}-{}.pgsc", field.label()));
                pgsum::cohort::extract_cohort(source, &reference, &identity, &paths, &options, &dest).unwrap();
                let table = pgsum::cohort::CohortTable::open(&dest).unwrap();
                table.validate_diagnostics().unwrap();
                for (position, reference, alternate) in
                    [(1, b'A', b'G'), (2, b'C', b'T'), (3, b'G', b'T'), (4, b'T', b'C')]
                {
                    let key = pgsum::genotypes::target_key(1, position, reference, alternate);
                    assert_eq!(table.source_records(key).unwrap(), Some(vec![position as u64]));
                }
                assert_eq!(
                    table
                        .call_state(pgsum::genotypes::target_key(1, 4, b'T', b'C'), 1)
                        .unwrap(),
                    Some(pgsum::genotype::State::Filtered)
                );
                let result = pgsum::cohort::score_cohort(&pack, &table, &Default::default()).unwrap();
                results.push((
                    result.text(0),
                    result.text(1),
                    result.scorable,
                    result.expected_genotype_terms,
                ));
            }
            assert_eq!(results[0], results[1]);
        }
    }
}

#[test]
fn truncated_bcf_and_unsupported_version_are_errors() {
    let root = std::env::temp_dir().join(format!("pgsum-invalid-bcf-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    for (name, bytes) in [
        ("truncated", b"BCF\x02\x02".as_slice()),
        ("version", b"BCF\x03\x09".as_slice()),
    ] {
        let input = root.join(name);
        std::fs::write(&input, bytes).unwrap();
        let result = pgsum::extract::scan(&input, &[], &[], &Default::default());
        assert!(result.is_err());
    }
}
