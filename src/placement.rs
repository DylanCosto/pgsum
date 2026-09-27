//! Placing a score within one reference group of a PLINK panel (`score --reference-panel <panel>.bed`).
//!
//! The panel is a PLINK 1 fileset (`.bed`, `.bim`, `.fam`) with biallelic lines, `A1` the ALT and `A2` the REF
//! allele, multi-allelic sites split into one line per ALT (as `plink2 --make-bed` writes them). A groups file
//! names each sample's group; every listed panel sample is scored, and the sample is placed among one group.
//!
//! The term set is the sample's own:
//! * each term scorable for the sample is matched to the panel line at its position with the same allele pair.
//!   An effect allele equal to the line's ALT adds `w × ALT dosage`; one equal to its REF adds
//!   `2w − w × (ALT dosage of every single-base line there with that REF)`, exact on split multi-allelic sites.
//!   A missing panel genotype takes the dosage `2 × f`, `f` the ALT frequency among the group's called samples.
//! * a scorable term with no panel line gets, for everyone, the expected value `w × 2 × f` from the frequency
//!   table (`--fill-frequencies`), replacing the sample's own contribution; one that cannot be filled is left
//!   out for everyone and reported, with a sensitivity result that gives it a homozygous-reference dosage in
//!   every panel sample instead.
//! * the sample's filled terms (`fill`) add the same amount to everyone, and its omitted terms nothing.
//!
//! Panel scores are exact sums. The placement uses their nearest floating-point values: the mid-rank
//! percentile `100 × (#{ref < x} + ½ #{ref = x}) / n`, the mean, the sample SD (`n − 1`), and
//! `Z = (x − mean) / SD`, with sums as Python's `math.fsum` computes them.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use num_bigint::BigInt;
use num_traits::{Signed, ToPrimitive, Zero};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::decimal::Decimal;
use crate::fill::Outcome as FillOutcome;
use crate::frequencies::FrequencyTable;
use crate::genotypes::{GenotypeTable, unpack_key};
use crate::pack::Pack;
use crate::score::{Coefficient, Contribution, ExactSum, Options, Plan};
use crate::term::{CONTIGS, contig_code_of_name};
use crate::{Error, Result, invalid};

/// A PLINK 1 panel, memory-mapped.
pub struct PlinkPanel {
    pub path: PathBuf,
    pub samples: Vec<String>,
    /// Per line: contig code, position, REF (`A2`), ALT (`A1`).
    lines: Vec<(u8, u32, String, String)>,
    by_position: HashMap<(u8, u32), Vec<u32>>,
    bed: memmap2::Mmap,
    bytes_per_line: usize,
}

impl PlinkPanel {
    pub fn open(bed: &Path) -> Result<PlinkPanel> {
        let with = |ext: &str| bed.with_extension(ext);
        let fam = std::fs::read_to_string(with("fam")).map_err(Error::io(with("fam")))?;
        let samples: Vec<String> = fam
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| l.split_whitespace().nth(1).unwrap_or_default().to_owned())
            .collect();
        let bim = std::fs::read_to_string(with("bim")).map_err(Error::io(with("bim")))?;
        let mut lines = Vec::new();
        let mut by_position: HashMap<(u8, u32), Vec<u32>> = HashMap::new();
        for (i, line) in bim.lines().filter(|l| !l.trim().is_empty()).enumerate() {
            let f: Vec<&str> = line.split_whitespace().collect();
            let [chrom, _id, _cm, pos, a1, a2] = f[..] else {
                return invalid!("{}: line {} does not have six fields", with("bim").display(), i + 1);
            };
            let contig = match chrom {
                "PAR1" | "PAR2" => 23,
                _ => contig_code_of_name(chrom).unwrap_or(0),
            };
            let pos: u32 = pos
                .parse()
                .map_err(|_| Error::Invalid(format!("{}: invalid position {pos}", with("bim").display())))?;
            by_position.entry((contig, pos)).or_default().push(i as u32);
            lines.push((contig, pos, a2.to_owned(), a1.to_owned()));
        }
        let file = File::open(bed).map_err(Error::io(bed))?;
        // SAFETY: the panel is read-only while it is used.
        let map = unsafe { memmap2::Mmap::map(&file) }.map_err(Error::io(bed))?;
        let bytes_per_line = samples.len().div_ceil(4);
        if map.get(..3) != Some([0x6c, 0x1b, 0x01].as_slice()) || map.len() != 3 + bytes_per_line * lines.len() {
            return invalid!(
                "{}: not a SNP-major PLINK 1 .bed of {} samples × {} lines",
                bed.display(),
                samples.len(),
                lines.len()
            );
        }
        Ok(PlinkPanel {
            path: bed.to_owned(),
            samples,
            lines,
            by_position,
            bed: map,
            bytes_per_line,
        })
    }

    /// ALT dosage of every sample at `line` (`None` when missing).
    fn dosages(&self, line: u32) -> impl Iterator<Item = Option<u8>> + '_ {
        let start = 3 + line as usize * self.bytes_per_line;
        let bytes = &self.bed[start..start + self.bytes_per_line];
        (0..self.samples.len()).map(move |s| match (bytes[s / 4] >> (2 * (s % 4))) & 3 {
            0 => Some(2),
            2 => Some(1),
            3 => Some(0),
            _ => None,
        })
    }
}

/// Each panel sample's group, from a TSV with a `sample` column and the group column.
pub fn read_groups(path: &Path, column: &str) -> Result<HashMap<String, String>> {
    let text = std::fs::read_to_string(path).map_err(Error::io(path))?;
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap_or_default().split('\t').map(str::trim).collect();
    let (Some(s), Some(g)) = (
        header.iter().position(|h| *h == "sample"),
        header.iter().position(|h| *h == column),
    ) else {
        return invalid!(
            "{}: needs a header with `sample` and `{column}` columns",
            path.display()
        );
    };
    Ok(lines
        .filter(|l| !l.trim().is_empty())
        .filter_map(|l| {
            let f: Vec<&str> = l.split('\t').map(str::trim).collect();
            Some(((*f.get(s)?).to_owned(), (*f.get(g)?).to_owned()))
        })
        .collect())
}

/// Python's `math.fsum`: the correctly rounded sum (Shewchuk's partials).
pub fn fsum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut partials: Vec<f64> = Vec::new();
    for mut x in values {
        let mut i = 0;
        for j in 0..partials.len() {
            let mut y = partials[j];
            if x.abs() < y.abs() {
                std::mem::swap(&mut x, &mut y);
            }
            let hi = x + y;
            let lo = y - (hi - x);
            if lo != 0.0 {
                partials[i] = lo;
                i += 1;
            }
            x = hi;
        }
        partials.truncate(i);
        partials.push(x);
    }
    // Round the partials' exact sum, as CPython does.
    let mut n = partials.len();
    if n == 0 {
        return 0.0;
    }
    n -= 1;
    let mut hi = partials[n];
    let mut lo = 0.0;
    while n > 0 {
        n -= 1;
        let x = hi;
        let y = partials[n];
        hi = x + y;
        let yr = hi - x;
        lo = y - yr;
        if lo != 0.0 {
            break;
        }
    }
    if n > 0 && ((lo < 0.0 && partials[n - 1] < 0.0) || (lo > 0.0 && partials[n - 1] > 0.0)) {
        let y = lo * 2.0;
        let x = hi + y;
        let yr = x - hi;
        if y == yr {
            hi = x;
        }
    }
    hi
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Distribution {
    pub n: u64,
    pub mean: f64,
    pub sd: f64,
    pub min: f64,
    pub max: f64,
    pub n_below: u64,
    pub n_equal: u64,
    pub n_above: u64,
    pub percentile: f64,
    pub z: f64,
}

/// The sample's place among `reference`.
pub fn distribution(x: f64, reference: &[f64]) -> Distribution {
    let n = reference.len();
    let below = reference.iter().filter(|v| **v < x).count();
    let equal = reference.iter().filter(|v| **v == x).count();
    let mean = fsum(reference.iter().copied()) / n as f64;
    let sd = (fsum(reference.iter().map(|v| (v - mean) * (v - mean))) / (n as f64 - 1.0)).sqrt();
    Distribution {
        n: n as u64,
        mean,
        sd,
        min: reference.iter().copied().fold(f64::INFINITY, f64::min),
        max: reference.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        n_below: below as u64,
        n_equal: equal as u64,
        n_above: (n - below - equal) as u64,
        percentile: 100.0 * (below as f64 + 0.5 * equal as f64) / n as f64,
        z: (x - mean) / sd,
    }
}

/// An omitted term: its ordinal, position, alleles, weight, the sample's contribution and the reason.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct OmittedTerm {
    pub ordinal: u64,
    pub contig: String,
    pub pos: u32,
    pub ref_allele: String,
    pub alt: String,
    pub effect_allele: String,
    pub weight: String,
    pub sample_contribution: String,
    pub reason: String,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct PlacementTerms {
    /// Terms scorable for the sample.
    pub observed: u64,
    pub matched_in_panel: u64,
    pub absent_from_panel: u64,
    pub absent_filled: u64,
    pub absent_omitted: u64,
    /// Matched terms by how they were expressed on panel lines.
    pub matched_by: BTreeMap<String, u64>,
    /// Panel lines with a non-zero coefficient.
    pub panel_lines_scored: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct GroupPlacement {
    pub panel: String,
    pub group_column: String,
    pub group: String,
    /// Panel samples scored (listed in the groups file).
    pub samples_scored: u64,
    /// The sample's score over the placement's term set, and what it is made of.
    pub sample_score: String,
    pub raw_score: String,
    pub sample_fill_sum: String,
    pub absent_fill: String,
    pub absent_filled_sample_contribution: String,
    pub terms: PlacementTerms,
    /// `2w` summed over matched terms whose effect allele is the panel line's REF.
    pub effect_is_ref_constant: String,
    pub reference: Distribution,
    pub omitted: Vec<OmittedTerm>,
    pub omitted_abs_weight: String,
    pub omitted_sample_contribution: String,
    /// The omitted terms kept for the sample and given a homozygous-reference dosage in every panel sample.
    pub sensitivity_absent_as_homozygous_reference: SensitivityResult,
    pub missing_genotypes: MissingGenotypes,
    /// Where the sample's deviation from the group mean comes from.
    pub contributions: Contributions,
    /// Each scored panel sample's exact score (see `write_scores`).
    #[serde(skip)]
    pub scores: Vec<(String, String, Decimal)>,
}

/// A named interval (BED: 0-based start, end).
#[derive(Clone, Debug, PartialEq)]
pub struct Region {
    pub name: String,
    pub contig: u8,
    pub start0: u64,
    pub end: u64,
}

/// Regions from a BED file (plain or gzip); each named by its fourth column, else the file's name.
pub fn read_regions(path: &Path) -> Result<Vec<Region>> {
    let bytes = std::fs::read(path).map_err(Error::io(path))?;
    let text = if bytes.starts_with(&[0x1f, 0x8b]) {
        let mut s = String::new();
        std::io::Read::read_to_string(&mut flate2::read::MultiGzDecoder::new(&bytes[..]), &mut s)
            .map_err(Error::io(path))?;
        s
    } else {
        String::from_utf8(bytes).map_err(|_| Error::Invalid(format!("{}: not text", path.display())))?
    };
    let stem = path
        .file_name()
        .map(|n| {
            n.to_string_lossy()
                .trim_end_matches(".gz")
                .trim_end_matches(".bed")
                .to_owned()
        })
        .unwrap_or_default();
    text.lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#') && !l.starts_with("track"))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let bad = || Error::Invalid(format!("{}: invalid line {l:?}", path.display()));
            Ok(Region {
                name: f.get(3).map_or(stem.clone(), |n| (*n).to_owned()),
                contig: contig_code_of_name(f.first().ok_or_else(bad)?).ok_or_else(bad)?,
                start0: f.get(1).and_then(|v| v.parse().ok()).ok_or_else(bad)?,
                end: f.get(2).and_then(|v| v.parse().ok()).ok_or_else(bad)?,
            })
        })
        .collect()
}

/// Window size for `Contributions::top_window` (fixed windows aligned to the contig start).
pub const WINDOW: u32 = 1_000_000;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct RegionContribution {
    pub name: String,
    pub contig: String,
    pub start0: u64,
    pub end: u64,
    pub contribution: f64,
    /// `contribution / deviation`.
    pub share: Option<f64>,
    pub lines: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Contributions {
    /// The sample's deviation from the group mean, summed over panel lines.
    pub deviation: f64,
    pub lines: u64,
    /// The 1 Mb window contributing most in the deviation's direction.
    pub top_window: Option<RegionContribution>,
    /// Each region given with `--contribution-region`.
    pub regions: Vec<RegionContribution>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SensitivityResult {
    pub reference_shift: String,
    pub percentile: f64,
    pub z: f64,
    pub mean: f64,
    pub sd: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct MissingGenotypes {
    pub samples_with_missing: u64,
    pub max_per_sample: u64,
    pub total: u64,
}

fn to_big(c: &Contribution) -> BigInt {
    let m = match &c.coefficient {
        Coefficient::Small(v) => BigInt::from(*v),
        Coefficient::Big(v) => v.clone(),
    };
    if c.negative { -m } else { m }
}

fn scaled(c: &Contribution, exponent: i64) -> BigInt {
    to_big(c) * BigInt::from(10u8).pow((c.exponent - exponent) as u32)
}

fn decimal_of(value: &BigInt, exponent: i64) -> Decimal {
    Decimal {
        negative: value.is_negative(),
        coefficient: value.abs().to_string(),
        exponent,
    }
}

fn contribution_of(d: &Decimal) -> Contribution {
    Contribution {
        negative: d.negative,
        coefficient: Coefficient::Big(d.coefficient.parse().expect("digits")),
        exponent: d.exponent,
    }
}

fn to_f64(d: &Decimal) -> f64 {
    d.to_python_string().parse().unwrap_or(f64::NAN)
}

/// A term scorable for the sample.
struct Observed {
    ordinal: u64,
    contig: u8,
    pos: u32,
    ref_allele: u8,
    alt: u8,
    effect: u8,
    /// The sample's effect-allele dosage.
    dosage: u8,
    /// `w` (multiplier 1) and the sample's contribution.
    weight: Contribution,
    contribution: Contribution,
}

/// Positions of the terms scorable for the sample (the lines a panel needs for `place`).
pub fn scorable_positions(
    pack: &Pack,
    genotypes: &GenotypeTable,
    options: &Options,
    _sex: Option<crate::score::Sex>,
) -> Result<Vec<(u8, u32)>> {
    let mut out = Vec::new();
    for term in pack.terms_from(0) {
        let term = term?;
        if crate::score::outcome(&term, genotypes, options)?.contribution.is_some() && term.contig != 0 {
            out.push((term.contig, term.pos));
        }
    }
    Ok(out)
}

/// Place `pack`'s score for the sample among `group` of `panel`. `fill` is the sample's fill result, if any.
#[allow(clippy::too_many_arguments)]
pub fn place(
    pack: &Pack,
    genotypes: &GenotypeTable,
    options: &Options,
    panel: &PlinkPanel,
    groups: &HashMap<String, String>,
    column: &str,
    group: &str,
    frequencies: Option<&FrequencyTable>,
    raw_score: &Decimal,
    sample_fill_sum: &Decimal,
    regions: &[Region],
    sex: Option<crate::score::Sex>,
) -> Result<GroupPlacement> {
    // The sample's scorable SNV terms.
    let mut observed = Vec::new();
    for (i, term) in pack.terms_from(0).enumerate() {
        let term = term?;
        let o = crate::score::outcome(&term, genotypes, options)?;
        let Some(c) = o.contribution else { continue };
        if term.model != crate::term::Model::Additive {
            return invalid!(
                "{}: term {} is not additive; placement takes additive scores",
                pack.header.pgs_id,
                i + 1
            );
        }
        let Plan::Snv { key, effect_is_alt } = crate::score::plan(&term, options).0 else {
            return invalid!(
                "{}: term {} is a scorable sequence term; placement takes SNVs",
                pack.header.pgs_id,
                i + 1
            );
        };
        let (contig, pos, r, a) = unpack_key(key);
        if sex.is_some_and(|s| s.xy) && contig == 23 && !crate::score::in_par(pos) {
            return invalid!(
                "{}: placing a score with hemizygous chrX terms for an XY sample is not supported yet",
                pack.header.pgs_id
            );
        }
        let weight = crate::score::contribution(crate::term::Model::Additive, &term.weights, 1)
            .ok_or_else(|| Error::Invalid(format!("{}: term {} has no weight", pack.header.pgs_id, i + 1)))?;
        observed.push(Observed {
            ordinal: i as u64 + 1,
            contig,
            pos,
            ref_allele: r,
            alt: a,
            effect: if effect_is_alt { a } else { r },
            dosage: o.effect_dosage.unwrap_or(0),
            weight,
            contribution: c,
        });
    }
    // Coefficients on panel lines, each held at the smallest exponent of the weights it sums, and what matched
    // how.
    let mut coef: HashMap<u32, (BigInt, i64)> = HashMap::new();
    let add = |coef: &mut HashMap<u32, (BigInt, i64)>, line: u32, w: &Contribution, negate: bool| {
        let entry = coef.entry(line).or_insert_with(|| (BigInt::zero(), w.exponent));
        if w.exponent < entry.1 {
            entry.0 *= BigInt::from(10u8).pow((entry.1 - w.exponent) as u32);
            entry.1 = w.exponent;
        }
        let v = scaled(w, entry.1);
        if negate { entry.0 -= v } else { entry.0 += v }
    };
    let mut constant = ExactSum::default();
    let mut matched_by: BTreeMap<String, u64> = BTreeMap::new();
    let (mut absent_fill, mut absent_actual) = (ExactSum::default(), ExactSum::default());
    let (mut omitted_abs, mut omitted_sample) = (ExactSum::default(), ExactSum::default());
    let mut homref_shift = ExactSum::default();
    let (mut matched, mut absent, mut filled) = (0u64, 0u64, 0u64);
    let mut omitted = Vec::new();
    // The sample's ALT dosage on each line a term touches (for `contributions`).
    let mut sample_alt: HashMap<u32, u8> = HashMap::new();
    let base = |b: u8| (b as char).to_string();
    for t in &observed {
        let here: &[u32] = panel.by_position.get(&(t.contig, t.pos)).map_or(&[], |v| v.as_slice());
        let pair = |l: u32| {
            let (_, _, r, a) = &panel.lines[l as usize];
            let (mut x, mut y) = ([r.as_str(), a.as_str()], [base(t.ref_allele), base(t.alt)]);
            x.sort_unstable();
            y.sort_unstable();
            x[0] == y[0] && x[1] == y[1]
        };
        let matches: Vec<u32> = here.iter().copied().filter(|&l| pair(l)).collect();
        if matches.len() > 1 {
            return invalid!("{}: several panel lines match term {}", pack.header.pgs_id, t.ordinal);
        }
        let Some(&line) = matches.first() else {
            absent += 1;
            let effect = base(t.effect);
            let other = if t.effect == t.ref_allele {
                base(t.alt)
            } else {
                base(t.ref_allele)
            };
            let w2 = Contribution {
                coefficient: Coefficient::Big(to_big(&t.weight).abs() * BigInt::from(2u8)),
                ..t.weight.clone()
            };
            let outcome = match frequencies {
                Some(table) => {
                    crate::fill::resolve_alleles(t.contig, t.pos, &effect, &other, "observed", table, Some(w2))?
                }
                None => FillOutcome::Omitted("no_frequency_table"),
            };
            match outcome {
                FillOutcome::Filled { contribution, .. } => {
                    filled += 1;
                    absent_fill.add(&contribution);
                    absent_actual.add(&t.contribution);
                }
                FillOutcome::Omitted(reason) => {
                    omitted_sample.add(&t.contribution);
                    omitted_abs.add(&Contribution {
                        negative: false,
                        ..t.weight.clone()
                    });
                    if t.effect == t.ref_allele {
                        homref_shift.add(&Contribution {
                            coefficient: Coefficient::Big(to_big(&t.weight).abs() * BigInt::from(2u8)),
                            ..t.weight.clone()
                        });
                    }
                    omitted.push(OmittedTerm {
                        ordinal: t.ordinal,
                        contig: CONTIGS[t.contig as usize - 1].to_owned(),
                        pos: t.pos,
                        ref_allele: base(t.ref_allele),
                        alt: base(t.alt),
                        effect_allele: effect,
                        weight: t.weight.to_decimal().to_python_string(),
                        sample_contribution: t.contribution.to_decimal().to_python_string(),
                        reason: reason.to_owned(),
                    });
                }
            }
            continue;
        };
        matched += 1;
        let (_, _, pref, palt) = &panel.lines[line as usize];
        if base(t.effect) == *palt {
            add(&mut coef, line, &t.weight, false);
            sample_alt.entry(line).or_insert(t.dosage);
            *matched_by.entry("effect_is_panel_alt".into()).or_default() += 1;
        } else {
            constant.add(&Contribution {
                coefficient: Coefficient::Big(to_big(&t.weight).abs() * BigInt::from(2u8)),
                ..t.weight.clone()
            });
            let others: Vec<u32> = here
                .iter()
                .copied()
                .filter(|&l| {
                    let (_, _, r, a) = &panel.lines[l as usize];
                    r == pref && r.len() == 1 && a.len() == 1
                })
                .collect();
            for &l in &others {
                add(&mut coef, l, &t.weight, true);
                let other = if t.effect == t.ref_allele { t.alt } else { t.ref_allele };
                let dose = if panel.lines[l as usize].3 == base(other) {
                    2 - t.dosage
                } else {
                    0
                };
                sample_alt.entry(l).or_insert(dose);
            }
            let kind = if others.len() > 1 {
                "effect_is_panel_ref_multiallelic"
            } else {
                "effect_is_panel_ref"
            };
            *matched_by.entry(kind.into()).or_default() += 1;
        }
        if *pref != base(t.ref_allele) {
            *matched_by
                .entry("panel_ref_differs_from_grch38_ref".into())
                .or_default() += 1;
        }
    }
    coef.retain(|_, (v, _)| !v.is_zero());
    let mut lines: Vec<(u32, BigInt, i64)> = coef.into_iter().map(|(l, (v, e))| (l, v, e)).collect();
    lines.sort_unstable_by_key(|(l, _, _)| *l);

    // Scored samples and the reference group.
    let scored: Vec<usize> = (0..panel.samples.len())
        .filter(|&s| groups.contains_key(&panel.samples[s]))
        .collect();
    let in_group: Vec<bool> = panel
        .samples
        .iter()
        .map(|s| groups.get(s).is_some_and(|g| g == group))
        .collect();
    if !scored.iter().any(|&s| in_group[s]) {
        return invalid!("no panel sample is in group {group}");
    }
    let n = panel.samples.len();
    let group_frequency = |line: u32| -> f64 {
        let (mut alt, mut called) = (0u64, 0u64);
        for (s, d) in panel.dosages(line).enumerate() {
            if let (true, Some(d)) = (in_group[s], d) {
                alt += d as u64;
                called += 1;
            }
        }
        if called == 0 {
            f64::NAN
        } else {
            alt as f64 / (2 * called) as f64
        }
    };
    // Per-sample sums, exact: lines are grouped by exponent and summed in 128-bit integers per group (a
    // coefficient too wide for that is summed in arbitrary precision). A missing genotype adds `coef × 2f` in
    // floating point, reported with `missing_genotypes`.
    let mut exponents: Vec<i64> = lines.iter().map(|(_, _, e)| *e).collect();
    exponents.sort_unstable();
    exponents.dedup();
    let bucket = |e: i64| exponents.binary_search(&e).expect("listed");
    let (mut small, mut wide) = (Vec::new(), Vec::new());
    for (l, v, e) in &lines {
        match v.to_i128().filter(|x| x.unsigned_abs() < 1 << 100) {
            Some(x) => small.push((*l, x, bucket(*e))),
            None => wide.push((*l, v.clone(), *e)),
        }
    }
    let buckets = exponents.len();
    type Part = (Vec<i128>, Vec<f64>, Vec<u64>, Vec<(u32, u64, u64)>);
    let parts: Vec<Part> = small
        .par_chunks(4096)
        .map(|chunk| {
            let (mut acc, mut imp, mut miss) = (vec![0i128; buckets * n], vec![0f64; n], vec![0u64; n]);
            let mut counts = Vec::with_capacity(chunk.len());
            for &(line, c, b) in chunk {
                let (mut group_alt, mut group_called) = (0u64, 0u64);
                let mut frequency = None;
                // PLINK codes: 0 two A1 (ALT), 1 missing, 2 one, 3 none.
                let values = [2 * c, 0, c, 0];
                let acc = &mut acc[b * n..(b + 1) * n];
                let start = 3 + line as usize * panel.bytes_per_line;
                for (i, &byte) in panel.bed[start..start + panel.bytes_per_line].iter().enumerate() {
                    for k in 0..4 {
                        let s = 4 * i + k;
                        if s >= n {
                            break;
                        }
                        let code = ((byte >> (2 * k)) & 3) as usize;
                        if code == 1 {
                            let f = *frequency.get_or_insert_with(|| group_frequency(line));
                            imp[s] += c as f64 * 10f64.powi(exponents[b] as i32) * 2.0 * f;
                            miss[s] += 1;
                        } else {
                            acc[s] += values[code];
                            if in_group[s] {
                                group_alt += [2, 0, 1, 0][code];
                                group_called += 1;
                            }
                        }
                    }
                }
                counts.push((line, group_alt, group_called));
            }
            (acc, imp, miss, counts)
        })
        .collect();
    let mut acc = vec![0i128; buckets * n];
    let (mut imputed, mut missing) = (vec![0f64; n], vec![0u64; n]);
    let mut group_counts: HashMap<u32, (u64, u64)> = HashMap::new();
    for (a, i, m, counts) in parts {
        group_counts.extend(counts.into_iter().map(|(l, x, y)| (l, (x, y))));
        for (x, y) in acc.iter_mut().zip(a) {
            *x += y;
        }
        for s in 0..n {
            imputed[s] += i[s];
            missing[s] += m[s];
        }
    }
    let mut wide_sums = vec![ExactSum::default(); n];
    for (line, c, e) in &wide {
        let (mut group_alt, mut group_called) = (0u64, 0u64);
        for (s, d) in panel.dosages(*line).enumerate() {
            if let (true, Some(d)) = (in_group[s], d) {
                group_alt += d as u64;
                group_called += 1;
            }
        }
        group_counts.insert(*line, (group_alt, group_called));
        let mut frequency = None;
        for (s, d) in panel.dosages(*line).enumerate() {
            match d {
                Some(d) => wide_sums[s].add(&contribution_of(&decimal_of(&(c * d), *e))),
                None => {
                    let f = *frequency.get_or_insert_with(|| group_frequency(*line));
                    imputed[s] += c.to_f64().unwrap_or(f64::NAN) * 10f64.powi(*e as i32) * 2.0 * f;
                    missing[s] += 1;
                }
            }
        }
    }
    // Everyone gets the constants: the effect-is-REF 2w, the absent-term fill and the sample's fill.
    let mut shift = constant.clone();
    shift.merge(absent_fill.clone());
    shift.add(&contribution_of(sample_fill_sum));
    let shift = shift.finish();
    let mut scores = Vec::with_capacity(scored.len());
    let (mut reference, mut reference_sensitivity) = (Vec::new(), Vec::new());
    let homref = homref_shift.finish();
    for &s in &scored {
        let mut exact = wide_sums[s].clone();
        for (b, e) in exponents.iter().enumerate() {
            exact.add(&contribution_of(&decimal_of(&BigInt::from(acc[b * n + s]), *e)));
        }
        exact.add(&contribution_of(&shift));
        let value = exact.finish();
        let float = to_f64(&value) + imputed[s];
        if in_group[s] {
            reference.push(float);
            let mut with = exact.clone();
            with.add(&contribution_of(&homref));
            reference_sensitivity.push(to_f64(&with.finish()) + imputed[s]);
        }
        scores.push((panel.samples[s].clone(), groups[&panel.samples[s]].clone(), value));
    }
    // The sample: its filled score with its own contribution at absent terms replaced by their fill, less its
    // contribution at terms left out for everyone.
    let mut comparison = ExactSum::default();
    comparison.add(&contribution_of(raw_score));
    comparison.add(&contribution_of(&omitted_sample.finish()).negated());
    comparison.add(&contribution_of(&absent_actual.finish()).negated());
    comparison.merge(absent_fill.clone());
    let comparison = comparison.finish();
    let placed = distribution(to_f64(&comparison), &reference);
    let sensitivity = distribution(to_f64(raw_score), &reference_sensitivity);
    let missing_in_scored: Vec<u64> = scored.iter().map(|&s| missing[s]).collect();
    // Each line's share of the sample's deviation from the group mean: `coef × (sample ALT dosage − group mean)`.
    let mut total = Vec::new();
    let mut windows: BTreeMap<(u8, u32), (Vec<f64>, u64)> = BTreeMap::new();
    let mut in_regions: Vec<(Vec<f64>, u64)> = vec![(Vec::new(), 0); regions.len()];
    for (line, v, e) in &lines {
        let (Some(&x), Some(&(alt, called))) = (sample_alt.get(line), group_counts.get(line)) else {
            continue;
        };
        if called == 0 {
            continue;
        }
        let (contig, pos, _, _) = &panel.lines[*line as usize];
        let d = to_f64(&decimal_of(v, *e)) * (x as f64 - alt as f64 / called as f64);
        total.push(d);
        let w = windows.entry((*contig, pos / WINDOW)).or_default();
        w.0.push(d);
        w.1 += 1;
        for (r, slot) in regions.iter().zip(in_regions.iter_mut()) {
            if r.contig == *contig && r.start0 < *pos as u64 && *pos as u64 <= r.end {
                slot.0.push(d);
                slot.1 += 1;
            }
        }
    }
    let deviation = fsum(total.iter().copied());
    let sign = if deviation >= 0.0 { 1.0 } else { -1.0 };
    let share = |x: f64| (deviation != 0.0).then(|| x / deviation);
    let top_window = windows
        .iter()
        .map(|((c, k), (parts, lines))| (*c, *k, fsum(parts.iter().copied()), *lines))
        .max_by(|a, b| (sign * a.2).total_cmp(&(sign * b.2)))
        .map(|(c, k, x, lines)| RegionContribution {
            name: format!("{}:{}-{} Mb", CONTIGS[c as usize - 1], k, k + 1),
            contig: CONTIGS[c as usize - 1].to_owned(),
            start0: k as u64 * WINDOW as u64,
            end: (k + 1) as u64 * WINDOW as u64,
            contribution: x,
            share: share(x),
            lines,
        });
    let contributions = Contributions {
        deviation,
        lines: total.len() as u64,
        top_window,
        regions: regions
            .iter()
            .zip(in_regions)
            .map(|(r, (parts, lines))| {
                let x = fsum(parts);
                RegionContribution {
                    name: r.name.clone(),
                    contig: CONTIGS[r.contig as usize - 1].to_owned(),
                    start0: r.start0,
                    end: r.end,
                    contribution: x,
                    share: share(x),
                    lines,
                }
            })
            .collect(),
    };
    Ok(GroupPlacement {
        panel: panel
            .path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default(),
        group_column: column.to_owned(),
        group: group.to_owned(),
        samples_scored: scored.len() as u64,
        sample_score: comparison.to_plain_string(),
        raw_score: raw_score.to_python_string(),
        sample_fill_sum: sample_fill_sum.to_python_string(),
        absent_fill: absent_fill.finish().to_python_string(),
        absent_filled_sample_contribution: absent_actual.finish().to_python_string(),
        terms: PlacementTerms {
            observed: observed.len() as u64,
            matched_in_panel: matched,
            absent_from_panel: absent,
            absent_filled: filled,
            absent_omitted: absent - filled,
            matched_by,
            panel_lines_scored: lines.len() as u64,
        },
        effect_is_ref_constant: constant.finish().to_python_string(),
        reference: placed,
        omitted,
        omitted_abs_weight: omitted_abs.finish().to_python_string(),
        omitted_sample_contribution: omitted_sample.finish().to_python_string(),
        sensitivity_absent_as_homozygous_reference: SensitivityResult {
            reference_shift: homref.to_python_string(),
            percentile: sensitivity.percentile,
            z: sensitivity.z,
            mean: sensitivity.mean,
            sd: sensitivity.sd,
        },
        missing_genotypes: MissingGenotypes {
            samples_with_missing: missing_in_scored.iter().filter(|m| **m > 0).count() as u64,
            max_per_sample: missing_in_scored.iter().copied().max().unwrap_or(0),
            total: missing_in_scored.iter().sum(),
        },
        contributions,
        scores,
    })
}

/// Write each scored panel sample's exact score as `sample`, `group`, `score`.
pub fn write_scores(placement: &GroupPlacement, path: &Path) -> Result<()> {
    let mut out = std::io::BufWriter::new(File::create(path).map_err(Error::io(path))?);
    let io = Error::io(path);
    let mut text = String::from("sample\tgroup\tscore\n");
    for (sample, group, score) in &placement.scores {
        text.push_str(&format!("{sample}\t{group}\t{}\n", score.to_plain_string()));
    }
    out.write_all(text.as_bytes()).map_err(io)?;
    out.flush().map_err(Error::io(path))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fsum_is_correctly_rounded() {
        assert_eq!(fsum([0.1; 10]), 1.0);
        assert_eq!(fsum([1e100, 1.0, -1e100, 1e-100, 1e50, -1.0, -1e50]), 1e-100);
        assert_eq!(
            fsum([2.0_f64.powi(53), -0.5, -2.0_f64.powi(-54)]),
            2.0_f64.powi(53) - 1.0
        );
    }

    #[test]
    fn distribution_matches_the_definitions() {
        let d = distribution(2.0, &[1.0, 2.0, 3.0, 4.0]);
        assert_eq!((d.n_below, d.n_equal, d.n_above), (1, 1, 2));
        assert_eq!(d.percentile, 37.5);
        assert_eq!(d.mean, 2.5);
        assert!((d.sd - (5.0f64 / 3.0).sqrt()).abs() < 1e-15);
    }
}
