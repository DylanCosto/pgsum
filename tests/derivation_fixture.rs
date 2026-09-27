//! Filling missing terms from allele frequencies and placing a score within a group of a PLINK panel, on the
//! synthetic gVCF and PGS999998 (one additive SNV term per genotype rule; see `pipeline_fixture.rs`).
//!
//! Synthetic chr1 reference bases 1–10: A C G T T G C A A C. PGS999998's scorable terms include
//! 1 (pos 1, effect G on REF A: effect is ALT, dosage 0), 2 (pos 2, effect C = REF, dosage 2) and
//! 3 (pos 3, effect G = REF, dosage 2); term 7 (pos 7, effect T on REF C, weight 0.7) and term 8 (pos 8,
//! effect G on REF A, weight 0.8) are missing (a filtered genotype and a low-quality call).

use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use pgsum::compile::{compile_file, reference_identity};
use pgsum::decimal::Decimal;
use pgsum::extract::{Options, ScanMode, extract};
use pgsum::fill::Outcome;
use pgsum::frequencies::FrequencyTable;
use pgsum::pack::Pack;
use pgsum::reference::Reference;
use pgsum::score::{Extras, Options as ScoreOptions, score_with};

fn gz(path: &Path, text: &str) {
    let mut e = flate2::write::GzEncoder::new(std::fs::File::create(path).unwrap(), flate2::Compression::fast());
    e.write_all(text.as_bytes()).unwrap();
    e.finish().unwrap();
}

/// A population VCF: pos 7 C>T (EUR_AF 0.25), pos 8 A>G flagged multi-allelic (0.1), pos 9 two records,
/// pos 10 without a frequency, pos 11 C>A (0.4, biallelic) for an effect-is-REF lookup.
fn frequency_table(dir: &Path) -> PathBuf {
    let vcf = dir.join("pop.vcf.gz");
    gz(
        &vcf,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
         chr1\t7\t.\tC\tT\t.\tPASS\tAF=0.2;EUR_AF=0.25\tGT\t0|1\n\
         chr1\t8\t.\tA\tG\t.\tPASS\tEUR_AF=0.1;MULTI_ALLELIC\tGT\t0|0\n\
         chr1\t9\t.\tA\tC\t.\tPASS\tEUR_AF=0.3\tGT\t0|0\n\
         chr1\t9\t.\tA\tT\t.\tPASS\tEUR_AF=0.01\tGT\t0|0\n\
         chr1\t10\t.\tC\tT\t.\tPASS\tAF=0.5\tGT\t0|0\n\
         chr1\t11\t.\tC\tA\t.\tPASS\tEUR_AF=0.4\tGT\t0|0\n",
    );
    let table = dir.join("pop.pgsf");
    let header = pgsum::frequencies::build(&[vcf], "EUR_AF", "MULTI_ALLELIC", None, &table).unwrap();
    assert_eq!(header.records, 6);
    table
}

fn setup(name: &str) -> (PathBuf, Reference, Pack, pgsum::genotypes::GenotypeTable) {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let out = std::env::temp_dir().join(format!("pgsum-derivation-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let reference = Reference::open(&fixtures.join("synthetic.fa")).unwrap();
    let identity = reference_identity(&reference).unwrap();
    compile_file(
        &fixtures.join("PGS999998_hmPOS_GRCh38.txt.gz"),
        Some(&fixtures.join("PGS999998.metadata.json")),
        &reference,
        &identity,
        None,
        &out,
    )
    .unwrap();
    let pack_path = out.join("PGS999998.pgsp");
    let options = Options {
        targets_cache: None,
        haploid_xy_as_homozygous: false,
        accept_missing_quality: false,
        skip_structural_alleles: false,
        merge_split_records: false,
        term_positions: false,
        sample: None,
        scan: ScanMode::Full,
        threads: 2,
    };
    let (table, _) = extract(
        &fixtures.join("synthetic.g.vcf.gz"),
        &reference,
        &identity,
        std::slice::from_ref(&pack_path),
        &options,
    )
    .unwrap();
    (out, reference, Pack::open(&pack_path).unwrap(), table)
}

#[test]
fn missing_terms_are_filled_or_omitted_by_the_rules() {
    let (out, _, pack, table) = setup("fill");
    let frequencies = FrequencyTable::open(&frequency_table(&out)).unwrap();
    let extras = Extras {
        fill: Some(&frequencies),
        sex: None,
    };
    let mut tsv = Vec::new();
    let result = score_with(&pack, &table, &ScoreOptions::default(), &extras, Some(&mut tsv)).unwrap();
    let fill = result.fill.unwrap();
    // Term 7: 0.7 × 2 × 0.25. Term 8: effect is the ALT, so the multi-allelic flag does not matter: 0.8 × 2 × 0.1.
    // Term 9 (pos 9, effect C on A): one of the two records there carries C/A: 0.9 × 2 × 0.3.
    assert_eq!(fill.filled.sum, "1.050");
    assert_eq!(fill.filled.part.terms, 3);
    assert_eq!(
        fill.filled.by_method,
        BTreeMap::from([("orientation_direct:effect_is_alt".into(), 3)])
    );
    assert_eq!(fill.omitted.by_reason["frequency_missing_or_invalid"].terms, 1);
    assert_eq!(fill.omitted.by_reason["no_frequency_source_for_contig"].terms, 1);
    assert!(fill.omitted.by_reason["site_absent"].terms > 0);
    assert_eq!(
        fill.filled.part.terms + fill.omitted.part.terms,
        fill.terms - fill.scorable_terms
    );
    // 10 of 26 terms are scorable: far below 99%, so no filled score.
    assert!(!fill.terms_at_least_99_percent && fill.filled_score.is_none());
    assert_eq!(fill.observed_sum, result.partial.raw_score);
    let rows: Vec<Vec<String>> = String::from_utf8(tsv)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect();
    let row = |ordinal: usize| &rows[ordinal - 1];
    let last = rows[0].len();
    assert_eq!(
        row(7)[last - 3..],
        ["filled:orientation_direct:effect_is_alt", "0.25", "0.350"]
    );
    assert_eq!(
        row(8)[last - 3..],
        ["filled:orientation_direct:effect_is_alt", "0.1", "0.16"]
    );
    assert_eq!(
        row(9)[last - 3..],
        ["filled:orientation_direct:effect_is_alt", "0.3", "0.54"]
    );
    assert_eq!(row(10)[last - 3], "omitted:frequency_missing_or_invalid");
    assert_eq!(row(1)[last - 3], "");
    std::fs::remove_dir_all(&out).unwrap();
}

#[test]
fn allele_resolution() {
    let out = std::env::temp_dir().join(format!("pgsum-derivation-alleles-{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let t = FrequencyTable::open(&frequency_table(&out)).unwrap();
    let w2 = pgsum::score::Contribution {
        negative: false,
        coefficient: pgsum::score::Coefficient::Small(2),
        exponent: -1,
    };
    let fill =
        |pos, effect, other| pgsum::fill::resolve_alleles(1, pos, effect, other, "m", &t, Some(w2.clone())).unwrap();
    let omitted = |o: Outcome| match o {
        Outcome::Omitted(r) => r,
        Outcome::Filled { .. } => "filled",
    };
    // Effect on the REF of a biallelic record: 1 − f.
    match fill(11, "C", "A") {
        Outcome::Filled {
            method,
            frequency,
            contribution,
        } => {
            assert_eq!(method, "m:effect_is_ref");
            assert_eq!(frequency, "0.4");
            assert_eq!(contribution.to_decimal().to_python_string(), "0.12");
        }
        o => panic!("{o:?}"),
    }
    // Effect on the REF of a record flagged multi-allelic, or at a position with two records.
    assert_eq!(omitted(fill(8, "A", "G")), "multiallelic_reference_effect_allele");
    assert_eq!(omitted(fill(9, "A", "C")), "multiallelic_reference_effect_allele");
    assert_eq!(omitted(fill(9, "C", "T")), "alleles_differ");
    assert_eq!(omitted(fill(12, "C", "T")), "site_absent");
    assert_eq!(
        omitted(pgsum::fill::resolve_alleles(23, 7, "T", "C", "m", &t, Some(w2.clone())).unwrap()),
        "no_frequency_source_for_contig"
    );
    std::fs::remove_dir_all(&out).unwrap();
}

/// A PLINK 1 fileset: `lines` are (pos, REF, ALT, per-sample ALT dosage or None).
fn plink(dir: &Path, samples: &[&str], lines: &[(u32, &str, &str, Vec<Option<u8>>)]) -> PathBuf {
    let fam: String = samples.iter().map(|s| format!("0 {s} 0 0 0 -9\n")).collect();
    std::fs::write(dir.join("panel.fam"), fam).unwrap();
    let bim: String = lines
        .iter()
        .map(|(p, r, a, _)| format!("1\t1:{p}:{r}:{a}\t0\t{p}\t{a}\t{r}\n"))
        .collect();
    std::fs::write(dir.join("panel.bim"), bim).unwrap();
    let mut bed = vec![0x6c, 0x1b, 0x01];
    for (_, _, _, dosages) in lines {
        let mut bytes = vec![0u8; samples.len().div_ceil(4)];
        for (s, d) in dosages.iter().enumerate() {
            let code = match d {
                Some(2) => 0,
                None => 1,
                Some(1) => 2,
                _ => 3,
            };
            bytes[s / 4] |= code << (2 * (s % 4));
        }
        bed.extend(bytes);
    }
    std::fs::write(dir.join("panel.bed"), bed).unwrap();
    dir.join("panel.bed")
}

#[test]
fn placement_scores_every_panel_sample_exactly() {
    let (out, _, pack, table) = setup("place");
    let frequencies = FrequencyTable::open(&frequency_table(&out)).unwrap();
    // Scorable terms (ordinal: position, effect, weight, sample effect dosage): 1: pos 1, G on A, 0.1, 0;
    // 2: pos 2, C = REF, 0.2, 2; 3: pos 3, G = REF, -0.3, 2. Other scorable terms have no panel line and no
    // frequency record, so they are left out for everyone.
    let samples = ["R1", "R2", "R3", "O1"];
    let bed = plink(
        &out,
        &samples,
        &[
            (1, "A", "G", vec![Some(0), Some(1), Some(2), Some(1)]),
            // pos 2 is split multi-allelic: C>T and C>G.
            (2, "C", "T", vec![Some(0), Some(1), Some(0), Some(2)]),
            (2, "C", "G", vec![Some(1), Some(0), None, Some(0)]),
            (3, "G", "A", vec![Some(0), Some(2), Some(1), Some(1)]),
        ],
    );
    let groups = out.join("groups.tsv");
    std::fs::write(&groups, "sample\tsuper_pop\nR1\tREF\nR2\tREF\nR3\tREF\nO1\tOTHER\n").unwrap();
    let panel = pgsum::placement::PlinkPanel::open(&bed).unwrap();
    let groups = pgsum::placement::read_groups(&groups, "super_pop").unwrap();
    let raw = Decimal::parse("1.0").unwrap();
    let zero = Decimal::parse("0").unwrap();
    let p = pgsum::placement::place(
        &pack,
        &table,
        &ScoreOptions::default(),
        &panel,
        &groups,
        "super_pop",
        "REF",
        Some(&frequencies),
        &raw,
        &zero,
        &[],
        None,
    )
    .unwrap();
    assert_eq!(p.samples_scored, 4);
    assert_eq!(p.terms.matched_in_panel, 3);
    assert_eq!(p.terms.matched_by["effect_is_panel_alt"], 1);
    assert_eq!(p.terms.matched_by["effect_is_panel_ref_multiallelic"], 1);
    assert_eq!(p.terms.matched_by["effect_is_panel_ref"], 1);
    // Term 1: 0.1 × ALT(pos 1). Term 2: 2 × 0.2 − 0.2 × (ALT(C>T) + ALT(C>G)); a missing C>G genotype takes
    // 2 × its REF-group frequency, 1/4. Term 3: 2 × (−0.3) − (−0.3) × ALT(pos 3).
    let expected = |a1: f64, t: f64, g: f64, a3: f64| 0.1 * a1 + (0.4 - 0.2 * (t + g)) + (-0.6 + 0.3 * a3);
    let scores: BTreeMap<String, f64> = p
        .scores
        .iter()
        .map(|(s, _, v)| (s.clone(), v.to_python_string().parse().unwrap()))
        .collect();
    assert!((scores["R1"] - expected(0.0, 0.0, 1.0, 0.0)).abs() < 1e-12);
    assert!((scores["R2"] - expected(1.0, 1.0, 0.0, 2.0)).abs() < 1e-12);
    assert_eq!(p.missing_genotypes.total, 1);
    assert!((scores["O1"] - expected(1.0, 2.0, 0.0, 1.0)).abs() < 1e-12);
    // R3's missing C>G genotype is imputed in its placement value (not in its exact score).
    assert_eq!(p.reference.n, 3);
    let r3 = expected(2.0, 0.0, 0.5, 1.0);
    let reference = [scores["R1"], scores["R2"], r3];
    let d = pgsum::placement::distribution(p.reference.mean + 0.0, &reference);
    assert!((p.reference.mean - d.mean).abs() < 1e-12);
    std::fs::remove_dir_all(&out).unwrap();
}
