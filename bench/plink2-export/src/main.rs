//! Write the inputs plink2 needs to score exactly what `pgsum score` scores, for a like-for-like benchmark.
//!
//! Usage: plink2-export <genotypes.pgsg> <pack-list.txt> <out-dir>
//!
//! Only additive packs are exported (plink2 --score applies one model to all columns). For each pack, every
//! term with no review reasons and a resolved orientation is written, as `pgsum score` uses by default:
//!
//! - `genotypes.vcf`: one record per target used, `ID=chr:pos:ref:alt`, with pgsum's call as GT (`./.` unless
//!   the call passes).
//! - `alt.scores.tsv`, `ref.scores.tsv`: dense weight tables (ID, effect allele, one column per score) for
//!   terms whose effect allele is the ALT or the REF. Weights of terms repeated within a score are added.
//! - `columns.txt`: the PGS IDs, in column order.
//!
//! With `no-mean-imputation`, plink2's SCORE_SUM for a score is then the sum over its alt and ref tables,
//! which equals pgsum's partial raw score.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use pgsum::genotypes::{GenotypeTable, target_key, unpack_key};
use pgsum::orient::Status;
use pgsum::pack::Pack;
use pgsum::term::{CONTIGS, Model};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let [_, table, list, out] = &args[..] else {
        eprintln!("usage: plink2-export <genotypes.pgsg> <pack-list.txt> <out-dir>");
        std::process::exit(2);
    };
    let out = Path::new(out);
    std::fs::create_dir_all(out).unwrap();
    let started = std::time::Instant::now();
    let table = GenotypeTable::open(Path::new(table)).unwrap();
    let packs: Vec<PathBuf> =
        std::fs::read_to_string(list).unwrap().lines().filter(|l| !l.trim().is_empty()).map(PathBuf::from).collect();

    // (key, effect is ALT) -> per-score weight, by column.
    let mut columns: Vec<String> = Vec::new();
    let mut weights: BTreeMap<(u64, bool), HashMap<usize, f64>> = BTreeMap::new();
    for path in &packs {
        let pack = Pack::open(path).unwrap();
        if pack.columns.model.iter().any(|&m| m != Model::Additive as u8) {
            eprintln!("skipping {} (not purely additive)", pack.header.pgs_id);
            continue;
        }
        let column = columns.len();
        columns.push(pack.header.pgs_id.clone());
        for term in pack.terms() {
            let t = term.unwrap();
            let o = t.orientation;
            if !t.reasons.is_empty() || o.status != Status::Resolved {
                continue;
            }
            let Some(w) = t.weights[0].to_decimal() else { continue };
            let w: f64 = w.to_python_string().parse().unwrap();
            let key = target_key(t.contig, t.pos, o.ref_base, o.alt_base);
            *weights.entry((key, o.effect_is_alt)).or_default().entry(column).or_default() += w;
        }
    }
    eprintln!("{} scores, {} weight rows, read in {:.1}s", columns.len(), weights.len(), started.elapsed().as_secs_f64());

    let id = |key: u64| {
        let (c, p, r, a) = unpack_key(key);
        format!("{}:{}:{}:{}", CONTIGS[c as usize - 1], p, r as char, a as char)
    };
    let mut vcf = BufWriter::with_capacity(1 << 22, File::create(out.join("genotypes.vcf")).unwrap());
    writeln!(vcf, "##fileformat=VCFv4.2").unwrap();
    for c in CONTIGS {
        writeln!(vcf, "##contig=<ID={c}>").unwrap();
    }
    writeln!(vcf, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">").unwrap();
    writeln!(vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", table.header.sample.sample_id).unwrap();
    let mut keys: Vec<u64> = weights.keys().map(|(k, _)| *k).collect();
    keys.dedup();
    for &key in &keys {
        let (c, p, r, a) = unpack_key(key);
        let entry = table.get(key).expect("target in the genotype table");
        let gt = match entry.alt_dosage {
            Some(0) if entry.state.is_passing() => "0/0",
            Some(1) if entry.state.is_passing() => "0/1",
            Some(2) if entry.state.is_passing() => "1/1",
            _ => "./.",
        };
        writeln!(vcf, "{}\t{p}\t{}\t{}\t{}\t.\t.\t.\tGT\t{gt}", CONTIGS[c as usize - 1], id(key), r as char, a as char)
            .unwrap();
    }
    vcf.flush().unwrap();

    for (name, effect_is_alt) in [("alt", true), ("ref", false)] {
        let mut f = BufWriter::with_capacity(1 << 22, File::create(out.join(format!("{name}.scores.tsv"))).unwrap());
        write!(f, "ID\tA1").unwrap();
        for c in &columns {
            write!(f, "\t{c}").unwrap();
        }
        writeln!(f).unwrap();
        let mut row = vec![0f64; columns.len()];
        for ((key, side), per_score) in &weights {
            if *side != effect_is_alt {
                continue;
            }
            let (_, _, r, a) = unpack_key(*key);
            row.iter_mut().for_each(|v| *v = 0.0);
            for (&c, &w) in per_score {
                row[c] = w;
            }
            write!(f, "{}\t{}", id(*key), if effect_is_alt { a as char } else { r as char }).unwrap();
            for v in &row {
                if *v == 0.0 { write!(f, "\t0").unwrap() } else { write!(f, "\t{v}").unwrap() }
            }
            writeln!(f).unwrap();
        }
        f.flush().unwrap();
    }
    std::fs::write(out.join("columns.txt"), columns.join("\n") + "\n").unwrap();
    eprintln!("wrote {} variants in {:.1}s", keys.len(), started.elapsed().as_secs_f64());
}
