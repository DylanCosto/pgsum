//! A cohort file gives every sample exactly what extracting and scoring that sample alone gives.

use std::path::{Path, PathBuf};

use pgsum::cohort::{CohortOptions, CohortTable, MISSING, code, extract_cohort, score_cohort};
use pgsum::compile::{compile_file, reference_identity};
use pgsum::extract::{Options, ScanMode, extract};
use pgsum::genotypes::key_position;
use pgsum::pack::Pack;
use pgsum::reference::Reference;

fn read_gz(path: &Path) -> String {
    let mut s = String::new();
    std::io::Read::read_to_string(
        &mut flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap()),
        &mut s,
    )
    .unwrap();
    s
}

#[test]
fn cohort_matches_single_samples() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let out = std::env::temp_dir().join(format!("pgsum-cohort-test-{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let reference = Reference::open(&fixtures.join("synthetic.fa")).unwrap();
    let identity = reference_identity(&reference).unwrap();
    let mut packs: Vec<PathBuf> = Vec::new();
    for id in ["PGS999997", "PGS999998"] {
        compile_file(
            &fixtures.join(format!("{id}_hmPOS_GRCh38.txt.gz")),
            Some(&fixtures.join(format!("{id}.metadata.json"))),
            &reference,
            &identity,
            None,
            &out,
        )
        .unwrap();
        packs.push(out.join(format!("{id}.pgsp")));
    }
    // Three samples: the fixture's own, one with every heterozygous GT made homozygous ALT, one with no calls.
    let text = read_gz(&fixtures.join("synthetic.g.vcf.gz"));
    let columns = |line: &str| -> [String; 3] {
        let sample = line.rsplit('\t').next().unwrap().to_owned();
        let alt = match sample.split_once(':') {
            Some((gt, rest)) if gt == "0/1" || gt == "0|1" => format!("1/1:{rest}"),
            _ => sample.clone(),
        };
        [sample, alt, "./.".into()]
    };
    let mut cohort_vcf = String::new();
    let mut single: [String; 3] = Default::default();
    for line in text.lines() {
        if line.starts_with("##") {
            cohort_vcf.push_str(line);
            cohort_vcf.push('\n');
            for s in &mut single {
                s.push_str(line);
                s.push('\n');
            }
        } else if line.starts_with("#CHROM") {
            let base = line.rsplit_once('\t').unwrap().0;
            cohort_vcf.push_str(&format!("{base}\tA\tB\tC\n"));
            for (s, name) in single.iter_mut().zip(["A", "B", "C"]) {
                s.push_str(&format!("{base}\t{name}\n"));
            }
        } else {
            let base = line.rsplit_once('\t').unwrap().0;
            let c = columns(line);
            cohort_vcf.push_str(&format!("{base}\t{}\t{}\t{}\n", c[0], c[1], c[2]));
            for (s, field) in single.iter_mut().zip(&c) {
                s.push_str(&format!("{base}\t{field}\n"));
            }
        }
    }
    let cohort_path = out.join("cohort.vcf");
    std::fs::write(&cohort_path, &cohort_vcf).unwrap();
    let options = CohortOptions {
        threads: 2,
        ..CohortOptions::default()
    };
    let pgsc = out.join("cohort.pgsc");
    let summary = extract_cohort(&cohort_path, &reference, &identity, &packs, &options, &pgsc).unwrap();
    assert_eq!(summary.header.samples, ["A", "B", "C"]);
    let cohort = CohortTable::open(&pgsc).unwrap();
    let scores: Vec<_> = packs
        .iter()
        .map(|p| score_cohort(&Pack::open(p).unwrap(), &cohort, &Default::default()).unwrap())
        .collect();
    for (s, vcf) in single.iter().enumerate() {
        let path = out.join(format!("sample{s}.vcf"));
        std::fs::write(&path, vcf).unwrap();
        let opts = Options {
            threads: 2,
            scan: ScanMode::Full,
            ..Options::default()
        };
        let (table, _) = extract(&path, &reference, &identity, &packs, &opts).unwrap();
        // Every target's code equals the single-sample call.
        let mut targets = 0;
        for e in table.entries().unwrap() {
            let expected = match e.alt_dosage {
                Some(d) if e.state.is_passing() => d,
                _ => MISSING,
            };
            let row = cohort.row(e.key).unwrap().expect("target in the cohort file");
            assert_eq!(code(row, s), expected, "sample {s} at {:?}", key_position(e.key));
            targets += 1;
        }
        assert!(targets > 20);
        // Every score's exact partial sum and scorable count equal the single-sample score.
        for (p, c) in packs.iter().zip(&scores) {
            let single = pgsum::score::score(&Pack::open(p).unwrap(), &table, &Default::default(), None).unwrap();
            assert_eq!(c.text(s), single.partial.raw_score, "sample {s} {}", c.pgs_id);
            assert_eq!(c.scorable[s], single.scorable_terms, "sample {s} {}", c.pgs_id);
        }
    }
    // Placing sample A among a panel of A and B: over the terms both have, A's sum is its own panel value.
    let two = cohort_vcf
        .lines()
        .map(|l| {
            if l.starts_with("##") {
                l.to_owned()
            } else {
                l.rsplit_once('\t').unwrap().0.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    let two_path = out.join("two.vcf");
    std::fs::write(&two_path, two).unwrap();
    let panel_path = out.join("two.pgsc");
    extract_cohort(&two_path, &reference, &identity, &packs, &options, &panel_path).unwrap();
    let panel = CohortTable::open(&panel_path).unwrap();
    let labels = out.join("labels.tsv");
    std::fs::write(&labels, "sample\tsuper_pop\nA\tX\nB\tX\n").unwrap();
    let groups = pgsum::panel::read_groups(&labels, &panel, "super_pop").unwrap();
    let a = std::fs::write(out.join("sample0.vcf"), &single[0])
        .map(|_| out.join("sample0.vcf"))
        .unwrap();
    let opts = Options {
        threads: 2,
        scan: ScanMode::Full,
        ..Options::default()
    };
    let (table_a, _) = extract(&a, &reference, &identity, &packs, &opts).unwrap();
    let pack = Pack::open(&packs[0]).unwrap();
    let placed = pgsum::panel::place(&pack, &table_a, &panel, "two", &groups, None, &Default::default()).unwrap();
    assert!(placed.matched_terms > 0);
    assert_eq!(
        placed.meets_coverage_guideline,
        placed.matched_term_fraction >= 0.99 && placed.matched_weight_fraction >= 0.99
    );
    assert!(placed.meets_coverage_guideline || placed.note.contains("below the 99% guideline"));
    let panel_scores = score_cohort(&pack, &panel, &Default::default()).unwrap();
    let (va, vb): (f64, f64) = (placed.score.parse().unwrap(), panel_scores.text(1).parse().unwrap());
    assert_ne!(va, vb);
    assert_eq!(placed.percentile_all, Some(if va < vb { 25.0 } else { 75.0 }));
    assert!((placed.groups[0].mean - (va + vb) / 2.0).abs() < 1e-9);

    assert_eq!(
        (placed.fills.terms_filled, placed.fills.terms_constant),
        (0, 0),
        "A and B are called everywhere"
    );

    // Panel calls missing for some samples are filled for the comparison. 101 samples: 100 copies of A and
    // one sample (Z) missing at the first record overlapping a target; 100/101 are called (at least 99%), so
    // Z's call there is filled and every other sample keeps its own.
    let mut filled_vcf = String::new();
    let mut first_record = true;
    for line in text.lines() {
        if line.starts_with("##") {
            filled_vcf.push_str(line);
            filled_vcf.push('\n');
            continue;
        }
        let (base, field) = line.rsplit_once('\t').unwrap();
        if line.starts_with("#CHROM") {
            let names: Vec<String> = (0..100).map(|i| format!("A{i}")).chain(["Z".into()]).collect();
            filled_vcf.push_str(&format!("{base}\t{}\n", names.join("\t")));
            continue;
        }
        let z = if first_record && !field.starts_with("0/0") {
            "./."
        } else {
            field
        };
        if z == "./." {
            first_record = false;
        }
        filled_vcf.push_str(&format!("{base}\t{}\t{z}\n", vec![field; 100].join("\t")));
    }
    let filled_path = out.join("filled.vcf");
    std::fs::write(&filled_path, filled_vcf).unwrap();
    let filled_panel_path = out.join("filled.pgsc");
    extract_cohort(
        &filled_path,
        &reference,
        &identity,
        &packs,
        &options,
        &filled_panel_path,
    )
    .unwrap();
    let filled_panel = CohortTable::open(&filled_panel_path).unwrap();
    let labels = out.join("labels101.tsv");
    let rows: String = (0..100).map(|i| format!("A{i}\tX\n")).collect();
    std::fs::write(&labels, format!("sample\tsuper_pop\n{rows}Z\tX\n")).unwrap();
    let groups101 = pgsum::panel::read_groups(&labels, &filled_panel, "super_pop").unwrap();
    let with_fill = pgsum::panel::place(
        &pack,
        &table_a,
        &filled_panel,
        "filled",
        &groups101,
        None,
        &Default::default(),
    )
    .unwrap();
    let f = &with_fill.fills;
    assert!(f.terms_filled >= 1 && f.genotypes_filled == f.terms_filled, "{f:?}");
    assert_eq!(
        (f.per_sample[100] as u64, f.max_per_sample as u64, f.samples_filled),
        (f.genotypes_filled, f.genotypes_filled, 1)
    );
    assert!(with_fill.note.contains("filled for the comparison only"));

    // With only 2 samples one missing call is under 99%: the term gives both the group's expectation.
    let one_missing: String = filled_path_two(&text);
    let two_missing_path = out.join("two_missing.vcf");
    std::fs::write(&two_missing_path, one_missing).unwrap();
    let two_missing_panel = out.join("two_missing.pgsc");
    extract_cohort(
        &two_missing_path,
        &reference,
        &identity,
        &packs,
        &options,
        &two_missing_panel,
    )
    .unwrap();
    let tm = CohortTable::open(&two_missing_panel).unwrap();
    std::fs::write(out.join("labels_az.tsv"), "sample\tsuper_pop\nA\tX\nZ\tX\n").unwrap();
    let groups_az = pgsum::panel::read_groups(&out.join("labels_az.tsv"), &tm, "super_pop").unwrap();
    let constant = pgsum::panel::place(&pack, &table_a, &tm, "two", &groups_az, None, &Default::default()).unwrap();
    assert!(
        constant.fills.terms_constant >= 1 && constant.fills.terms_filled == 0,
        "{:?}",
        constant.fills
    );

    // A panel sample the groups file does not list is left out.
    let only_a = out.join("only_a.tsv");
    std::fs::write(&only_a, "sample\tsuper_pop\nA\tX\n").unwrap();
    let groups_a = pgsum::panel::read_groups(&only_a, &panel, "super_pop").unwrap();
    assert_eq!(groups_a.members, [0]);
    let alone = pgsum::panel::place(&pack, &table_a, &panel, "two", &groups_a, None, &Default::default()).unwrap();
    assert_eq!((alone.panel_samples, alone.percentile_all), (1, Some(50.0)));

    // The samples differ, so the comparison above is not vacuous.
    assert_ne!(scores[0].text(0), scores[0].text(1));
    assert_eq!(scores[0].scorable[2], 0);
    std::fs::remove_dir_all(&out).unwrap();
}

/// A two-sample VCF (A and Z) from the fixture's text: Z copies A except for a no-call at the first record that
/// is not a hom-ref block.
fn filled_path_two(text: &str) -> String {
    let mut out = String::new();
    let mut first = true;
    for line in text.lines() {
        if line.starts_with("##") {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        let (base, field) = line.rsplit_once('\t').unwrap();
        if line.starts_with("#CHROM") {
            out.push_str(&format!("{base}\tA\tZ\n"));
            continue;
        }
        let z = if first && !field.starts_with("0/0") {
            "./."
        } else {
            field
        };
        if z == "./." {
            first = false;
        }
        out.push_str(&format!("{base}\t{field}\t{z}\n"));
    }
    out
}
