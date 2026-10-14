//! Reference ambiguities stay unknown through compilation and optional allele inference.

use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use pgsum::compile::{compile_file, reference_identity};
use pgsum::extract::{Options as ExtractOptions, ScanMode, extract};
use pgsum::orient::Status;
use pgsum::pack::{Pack, ReferenceIdentity};
use pgsum::reference::Reference;
use pgsum::score::{Options as ScoreOptions, score};
use sha2::{Digest, Sha256};

static NEXT_DIRECTORY: AtomicU64 = AtomicU64::new(0);

struct Fixture {
    root: PathBuf,
    reference: Reference,
}

impl Fixture {
    fn new(sequence: &str, line_bases: usize) -> Self {
        let root = std::env::temp_dir().join(format!(
            "pgsum-reference-ambiguity-{}-{}",
            std::process::id(),
            NEXT_DIRECTORY.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&root).unwrap();
        let path = root.join("reference.fa");
        let line_bases = line_bases.min(sequence.len());
        let mut fasta = b">chr1\n".to_vec();
        for chunk in sequence.as_bytes().chunks(line_bases) {
            fasta.extend_from_slice(chunk);
            fasta.push(b'\n');
        }
        std::fs::write(&path, fasta).unwrap();
        std::fs::write(
            root.join("reference.fa.fai"),
            format!("chr1\t{}\t6\t{line_bases}\t{}\n", sequence.len(), line_bases + 1),
        )
        .unwrap();
        Self {
            root,
            reference: Reference::open(&path).unwrap(),
        }
    }

    fn compile(&self, rows: &str) -> (Pack, PathBuf, ReferenceIdentity) {
        let path = self.root.join("score.tsv");
        std::fs::write(
            &path,
            format!(
                "#pgs_id=AMBIGUITY\n#genome_build=GRCh38\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\thm_inferOtherAllele\n{rows}"
            ),
        )
        .unwrap();
        let identity = reference_identity(&self.reference).unwrap();
        compile_file(&path, None, &self.reference, &identity, None, &self.root).unwrap();
        let path = self.root.join("AMBIGUITY.pgsp");
        (Pack::open(&path).unwrap(), path, identity)
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.root);
    }
}

#[test]
fn all_dna_iupac_ambiguities_are_unknown_without_changing_reference_hashes() {
    let sequence = "ACGTNRYSWKMBDHVacgtnryswkmbdhv";
    let fixture = Fixture::new(sequence, 7);
    let original = std::fs::read(&fixture.reference.path).unwrap();
    let original_fai = std::fs::read(&fixture.reference.fai_path).unwrap();
    let expected = "ACGTNNNNNNNNNNNACGTNNNNNNNNNNN";
    assert_eq!(
        fixture.reference.fetch("chr1", 0, sequence.len() as u64).unwrap(),
        expected
    );
    // Slices spanning FASTA line boundaries use the same ambiguity handling.
    assert_eq!(fixture.reference.fetch("chr1", 5, 21).unwrap(), &expected[5..21]);
    let identity = reference_identity(&fixture.reference).unwrap();
    let sha = |bytes: &[u8]| format!("sha256:{}", pgsum::digest::hex(&Sha256::digest(bytes)));
    assert_eq!(identity.fasta_sha256, sha(&original));
    assert_eq!(identity.fai_sha256, sha(&original_fai));
    assert_eq!(std::fs::read(&fixture.reference.path).unwrap(), original);
    assert_eq!(std::fs::read(&fixture.reference.fai_path).unwrap(), original_fai);
}

#[test]
fn non_dna_iupac_characters_still_fail() {
    for invalid in ['U', 'u', 'X', 'x', 'Z', 'z', '-', '.', '?', '*', '0', ' '] {
        let fixture = Fixture::new(&format!("A{invalid}G"), 3);
        assert_eq!(fixture.reference.fetch("chr1", 0, 1).unwrap(), "A");
        let error = fixture.reference.fetch("chr1", 1, 2).unwrap_err().to_string();
        assert!(error.contains("unsupported base"), "{invalid:?}: {error}");
    }
    let fixture = Fixture::new("AXG", 3);
    let score_path = fixture.root.join("invalid.tsv");
    std::fs::write(
        &score_path,
        "#pgs_id=INVALID\n#genome_build=GRCh38\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n1\t2\tG\tA\t1\n",
    )
    .unwrap();
    let identity = reference_identity(&fixture.reference).unwrap();
    let error = compile_file(&score_path, None, &fixture.reference, &identity, None, &fixture.root)
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("term 1") && error.contains("unsupported base"),
        "{error}"
    );
}

#[test]
fn ambiguities_never_supply_snv_palindrome_or_sequence_inference() {
    let ambiguities = "RYSWKMBDHVNryswkmbdhvn";
    // 100 concrete forward-strand terms enable the palindrome convention. The next concrete A enables
    // missing-other-allele inference, so the ambiguous rows cannot hide behind inactive opt-in policies.
    let sequence = format!("{}{ambiguities}A", "A".repeat(101));
    let fixture = Fixture::new(&sequence, 17);
    let mut rows: String = (1..=100).map(|p| format!("1\t{p}\tG\tA\t1\t\n")).collect();
    rows.push_str("1\t101\tG\t\t2\tA\n");
    for pos in 102..102 + ambiguities.len() {
        for (effect, other, weight) in [("G", "A", 10), ("G", "", 20), ("A", "T", 30), ("GA", "AT", 40)] {
            rows.push_str(&format!("1\t{pos}\t{effect}\t{other}\t{weight}\tA\n"));
        }
    }
    rows.push_str(&format!("1\t{}\tT\tA\t3\t\n", sequence.len()));
    let (pack, path, identity) = fixture.compile(&rows);
    let terms: Vec<_> = pack.terms().map(Result::unwrap).collect();
    assert!(
        terms[100].inferred.is_some(),
        "concrete missing-other control must infer"
    );
    assert!(
        terms.last().unwrap().inferred.is_some(),
        "concrete palindrome control must infer"
    );
    for term in &terms[101..terms.len() - 1] {
        assert_ne!(term.orientation.status, Status::Resolved, "position {}", term.pos);
        assert!(term.inferred.is_none(), "ambiguity inferred at {}", term.pos);
        assert!(
            term.inferred_sequence.is_none(),
            "ambiguity supplied sequence at {}",
            term.pos
        );
    }
    let mut vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n".to_owned();
    for pos in 1..=sequence.len() {
        let alt = if pos == sequence.len() { "T" } else { "G" };
        vcf.push_str(&format!("1\t{pos}\t.\tA\t{alt}\t.\tPASS\t.\tGT\t0/1\n"));
    }
    let vcf_path = fixture.root.join("calls.vcf");
    std::fs::write(&vcf_path, vcf).unwrap();
    let (table, _) = extract(
        &vcf_path,
        &fixture.reference,
        &identity,
        &[path],
        &ExtractOptions {
            accept_missing_quality: true,
            scan: ScanMode::Full,
            ..ExtractOptions::default()
        },
    )
    .unwrap();
    let default = score(&pack, &table, &ScoreOptions::default(), None).unwrap();
    assert_eq!(
        (default.scorable_terms, default.partial.raw_score.as_str()),
        (100, "100")
    );
    let opt_in = ScoreOptions {
        allow_inferred_other_allele: true,
        allow_inferred_palindromes: true,
        allow_inferred_indels: true,
        accept_informational_descriptions: true,
    };
    let inferred = score(&pack, &table, &opt_in, None).unwrap();
    assert_eq!(
        (inferred.scorable_terms, inferred.partial.raw_score.as_str()),
        (102, "105")
    );
    for term in &terms[101..terms.len() - 1] {
        assert!(
            pgsum::score::outcome(term, &table, &opt_in)
                .unwrap()
                .contribution
                .is_none()
        );
    }
}

#[test]
fn sequence_inference_cannot_normalize_through_an_ambiguous_upstream_anchor() {
    for ambiguous in ['N', 'Y', 'r'] {
        let fixture = Fixture::new(&format!("A{ambiguous}TAGA"), 6);
        // TA fits at position 3; A does not. Deleting T needs the preceding reference base as its
        // canonical VCF anchor. That base is unknown and must not produce NT>N or an IUPAC allele.
        let (pack, path, identity) = fixture.compile("1\t3\tA\tTA\t10\t\n1\t6\tG\tA\t1\t\n");
        let term = pack.terms().next().unwrap().unwrap();
        assert_eq!(term.orientation.status, Status::SequenceReferenceEvidenceRequired);
        assert!(term.inferred.is_none());
        assert!(term.inferred_sequence.is_none());
        let targets = pgsum::targets::collect(&[path], &identity, false).unwrap();
        assert!(targets.sequences.is_empty());
        assert_eq!(targets.keys.len(), 1, "only the concrete SNV remains a target");
    }
}
