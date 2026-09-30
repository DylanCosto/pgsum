//! Auditable comparisons of reported scores; never silently convert averages into sums.
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use clap::{Args, ValueEnum};
use num_bigint::BigInt;
use num_traits::{Signed, Zero};
use serde::Serialize;

use crate::{Error, Result, invalid};

#[derive(Clone, Copy, Debug, ValueEnum, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Format {
    Pgsum,
    Plink2,
    PgscCalc,
}

#[derive(Args, Debug)]
pub struct CompareArgs {
    #[arg(long)]
    pub left: PathBuf,
    #[arg(long)]
    pub right: PathBuf,
    #[arg(long, value_enum, default_value_t = Format::Pgsum)]
    pub left_format: Format,
    #[arg(long, value_enum, default_value_t = Format::Pgsum)]
    pub right_format: Format,
    /// Sample identity for a single-sample pgsum scores.tsv (never inferred from filenames).
    #[arg(long)]
    pub left_sample: Option<String>,
    #[arg(long)]
    pub right_sample: Option<String>,
    /// Explicit PLINK sum column, e.g. PGS000001_SUM; required for PLINK input.
    #[arg(long)]
    pub left_score_column: Option<String>,
    #[arg(long)]
    pub right_score_column: Option<String>,
    /// Score identity for the selected PLINK column; required for PLINK input.
    #[arg(long)]
    pub left_score_id: Option<String>,
    #[arg(long)]
    pub right_score_id: Option<String>,
    /// Select one pgsc_calc sampleset. Required when more than one is present.
    #[arg(long)]
    pub left_sampleset: Option<String>,
    #[arg(long)]
    pub right_sampleset: Option<String>,
    /// Exact nonnegative tolerance: abs(delta) <= absolute + relative * max(abs(left), abs(right)).
    #[arg(long, default_value = "0")]
    pub absolute_tolerance: String,
    #[arg(long, default_value = "0")]
    pub relative_tolerance: String,
    /// Optional paired pgsum --terms exports (or the documented same-source-ordinal TSV schema).
    #[arg(long, requires_all = ["right_terms", "term_sample", "term_score"])]
    pub left_terms: Option<PathBuf>,
    #[arg(long, requires = "left_terms")]
    pub right_terms: Option<PathBuf>,
    #[arg(long, requires = "left_terms")]
    pub term_sample: Option<String>,
    #[arg(long, requires = "left_terms")]
    pub term_score: Option<String>,
    /// New output directory: scores.tsv, optional terms.jsonl, and comparison.json completion marker.
    #[arg(long)]
    pub out: PathBuf,
    /// Return failure after writing the report if any row is missing, unavailable or outside tolerance,
    /// or any supplied term evidence differs or fails to reconcile with its selected score.
    #[arg(long)]
    pub fail_on_difference: bool,
}

/// Broad enough for exact pgsum sums (weights/dosages each have their own narrower input bounds).
/// These bounds limit pathological exponent expansion while preserving up to 32,768 significant digits.
#[derive(Clone, Debug)]
struct Number {
    coefficient: BigInt,
    exponent: i32,
}
impl Number {
    fn parse(s: &str) -> Result<Self> {
        let bad = || Error::Invalid(format!("invalid finite comparison decimal {s:?}"));
        if s.is_empty() || s.len() > 32768 {
            return Err(bad());
        }
        let (negative, rest) = if let Some(t) = s.strip_prefix('-') {
            (true, t)
        } else {
            (false, s.strip_prefix('+').unwrap_or(s))
        };
        let mut pieces = rest.split(['e', 'E']);
        let mantissa = pieces.next().unwrap();
        let exponent: i32 = pieces
            .next()
            .map(str::parse)
            .transpose()
            .map_err(|_| bad())?
            .unwrap_or(0);
        if pieces.next().is_some() {
            return Err(bad());
        }
        let mut digits = String::new();
        let mut places = 0;
        let mut point = false;
        for b in mantissa.bytes() {
            match b {
                b'.' if !point => point = true,
                b'0'..=b'9' => {
                    digits.push(b as char);
                    if point {
                        places += 1;
                    }
                }
                _ => return Err(bad()),
            }
        }
        if digits.is_empty() {
            return Err(bad());
        }
        let exponent = exponent
            .checked_sub(places)
            .filter(|e| (-16384..=16384).contains(e))
            .ok_or_else(bad)?;
        let mut coefficient = BigInt::parse_bytes(digits.as_bytes(), 10).ok_or_else(bad)?;
        if negative {
            coefficient = -coefficient;
        }
        Ok(Self { coefficient, exponent })
    }
    fn zero() -> Self {
        Self {
            coefficient: BigInt::zero(),
            exponent: 0,
        }
    }
    fn at(&self, exponent: i32) -> BigInt {
        &self.coefficient * BigInt::from(10).pow((self.exponent - exponent) as u32)
    }
    fn add(&self, other: &Self) -> Self {
        let exponent = self.exponent.min(other.exponent);
        Self {
            coefficient: self.at(exponent) + other.at(exponent),
            exponent,
        }
    }
    fn sub(&self, other: &Self) -> Self {
        let exponent = self.exponent.min(other.exponent);
        Self {
            coefficient: self.at(exponent) - other.at(exponent),
            exponent,
        }
    }
    fn abs(&self) -> Self {
        Self {
            coefficient: self.coefficient.abs(),
            exponent: self.exponent,
        }
    }
    fn mul(&self, other: &Self) -> Self {
        Self {
            coefficient: &self.coefficient * &other.coefficient,
            exponent: self.exponent + other.exponent,
        }
    }
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let exponent = self.exponent.min(other.exponent);
        self.at(exponent).cmp(&other.at(exponent))
    }
    fn text(&self) -> String {
        // Scientific form avoids enormous strings for small differences.
        if self.coefficient.is_zero() {
            return "0".into();
        }
        let mut digits = self.coefficient.abs().to_string();
        let mut exponent = self.exponent;
        while digits.ends_with('0') {
            digits.pop();
            exponent += 1;
        }
        format!(
            "{}{digits}e{exponent}",
            if self.coefficient.is_negative() { "-" } else { "" }
        )
    }
}

struct Tolerance {
    absolute: Number,
    relative: Number,
}
impl Tolerance {
    fn new(absolute: &str, relative: &str) -> Result<Self> {
        let (absolute, relative) = (Number::parse(absolute)?, Number::parse(relative)?);
        if absolute.coefficient.is_negative() || relative.coefficient.is_negative() {
            return invalid!("comparison tolerances must be nonnegative");
        }
        Ok(Self { absolute, relative })
    }
    fn classify(&self, left: Option<&Number>, right: Option<&Number>) -> (&'static str, String) {
        let (Some(left), Some(right)) = (left, right) else {
            return ("unavailable", String::new());
        };
        let delta = right.sub(left);
        if delta.coefficient.is_zero() {
            return ("exact", "0".into());
        }
        let (a, b) = (left.abs(), right.abs());
        let scale = if a.cmp(&b).is_gt() { a } else { b };
        let limit = self.absolute.add(&self.relative.mul(&scale));
        (
            if delta.abs().cmp(&limit).is_le() {
                "within_tolerance"
            } else {
                "different"
            },
            delta.text(),
        )
    }
}

/// Tabular reader retains tab-separated empty cells; whitespace mode is for external score tables.
struct Table {
    path: PathBuf,
    reader: Box<dyn BufRead>,
    header: Vec<String>,
    tabs: bool,
    line: usize,
}
impl Table {
    fn open(path: &Path) -> Result<Self> {
        let mut raw = BufReader::new(File::open(path).map_err(Error::io(path))?);
        let magic = raw.fill_buf().map_err(Error::io(path))?;
        let decoded: Box<dyn Read> = if magic.starts_with(&[0x1f, 0x8b]) {
            Box::new(flate2::bufread::MultiGzDecoder::new(raw))
        } else if magic.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]) {
            Box::new(zstd::Decoder::with_buffer(raw).map_err(Error::io(path))?)
        } else {
            Box::new(raw)
        };
        let mut reader: Box<dyn BufRead> = Box::new(BufReader::new(decoded));
        let mut line = String::new();
        reader.read_line(&mut line).map_err(Error::io(path))?;
        let first = line
            .trim_end_matches(['\r', '\n'])
            .trim_start_matches('\u{feff}')
            .trim_start_matches('#');
        let tabs = first.contains('\t');
        let header: Vec<String> = if tabs {
            first.split('\t').map(str::to_owned).collect()
        } else {
            first.split_whitespace().map(str::to_owned).collect()
        };
        if header.is_empty()
            || header.iter().any(String::is_empty)
            || header.iter().collect::<BTreeSet<_>>().len() != header.len()
        {
            return invalid!("{}: empty or duplicate column header", path.display());
        }
        Ok(Self {
            path: path.into(),
            reader,
            header,
            tabs,
            line: 1,
        })
    }
    fn column(&self, name: &str) -> Option<usize> {
        self.header.iter().position(|x| x == name)
    }
    fn required(&self, name: &str) -> Result<usize> {
        self.column(name)
            .ok_or_else(|| Error::Invalid(format!("{}: required column {name:?} is absent", self.path.display())))
    }
    fn next(&mut self) -> Result<Option<Vec<String>>> {
        let mut text = String::new();
        loop {
            if self.reader.read_line(&mut text).map_err(Error::io(&self.path))? == 0 {
                return Ok(None);
            }
            self.line += 1;
            let content = text.trim_end_matches(['\r', '\n']);
            if content.is_empty() {
                text.clear();
                continue;
            }
            let row: Vec<String> = if self.tabs {
                content.split('\t').map(str::to_owned).collect()
            } else {
                content.split_whitespace().map(str::to_owned).collect()
            };
            if row.len() != self.header.len() {
                return invalid!(
                    "{}:{}: expected {} columns, got {}",
                    self.path.display(),
                    self.line,
                    self.header.len(),
                    row.len()
                );
            }
            if row.iter().any(|s| s.chars().any(|c| c.is_control())) {
                return invalid!("{}:{}: control character in field", self.path.display(), self.line);
            }
            return Ok(Some(row));
        }
    }
}

type Key = (String, String);
#[derive(Debug)]
struct Value {
    text: String,
    number: Option<Number>,
    evidence: BTreeMap<String, String>,
}
#[derive(Serialize)]
struct Source {
    path: PathBuf,
    sha256: String,
    format: Format,
    value_column: String,
    sample_override: Option<String>,
    sampleset: Option<String>,
    score_id: Option<String>,
    rows: usize,
    filtered_rows: usize,
}
struct Input<'a> {
    path: &'a Path,
    format: Format,
    sample: Option<&'a str>,
    column: Option<&'a str>,
    score: Option<&'a str>,
    sampleset: Option<&'a str>,
}
fn number(text: &str) -> Result<Option<Number>> {
    if matches!(text, "" | "NA" | "N/A" | "." | "nan" | "NaN") {
        Ok(None)
    } else {
        Number::parse(text).map(Some)
    }
}
fn identifier(s: &str) -> Result<()> {
    if s.is_empty() || s.chars().any(char::is_control) {
        return invalid!("empty identity or control character in identity");
    }
    Ok(())
}
fn read_scores(input: Input<'_>) -> Result<(BTreeMap<Key, Value>, Source)> {
    let before = crate::digest::file_sha256(input.path)?;
    let mut table = Table::open(input.path)?;
    let (sample_column, score_column, value_column) = match input.format {
        Format::Pgsum => {
            if input.column.is_some() || input.score.is_some() || input.sampleset.is_some() {
                return invalid!(
                    "pgsum comparisons always select partial_raw_score; column, score-id and sampleset options apply to external formats"
                );
            }
            let sample = table.column("sample_id").or_else(|| table.column("sample"));
            if table.column("sample_id").is_some() && table.column("sample").is_some() {
                return invalid!("ambiguous pgsum sample columns");
            }
            if sample.is_some() && input.sample.is_some() {
                return invalid!("sample override is only allowed for pgsum files without a sample column");
            }
            if sample.is_none() && input.sample.is_none() {
                return invalid!("single-sample pgsum scores.tsv requires an explicit --left-sample or --right-sample");
            }
            (sample, Some(table.required("pgs_id")?), "partial_raw_score".to_owned())
        }
        Format::Plink2 => {
            if input.sample.is_some() || input.sampleset.is_some() {
                return invalid!("PLINK samples use IID; sample/sampleset overrides are not supported");
            }
            let column = input.column.ok_or_else(|| {
                Error::Invalid("PLINK input requires an explicit score-column ending in _SUM and score-id".into())
            })?;
            if !column.ends_with("_SUM") || column == "NAMED_ALLELE_DOSAGE_SUM" || input.score.is_none() {
                return invalid!(
                    "PLINK input requires an explicit score-column ending in _SUM and score-id; averages are not sums"
                );
            }
            (Some(table.required("IID")?), None, column.into())
        }
        Format::PgscCalc => {
            if input.sample.is_some() || input.column.is_some() || input.score.is_some() {
                return invalid!(
                    "pgsc_calc comparisons use IID, PGS and SUM; sample/column/score-id overrides are not supported"
                );
            }
            (Some(table.required("IID")?), Some(table.required("PGS")?), "SUM".into())
        }
    };
    let vcol = table.required(&value_column)?;
    let sampleset_column = table.column("sampleset");
    if input.sampleset.is_some() && sampleset_column.is_none() {
        return invalid!("sampleset filter requested but sampleset column is absent");
    }
    let mut samplesets = BTreeSet::new();
    let mut identities = BTreeMap::new();
    let mut values = BTreeMap::new();
    let mut filtered_rows = 0;
    while let Some(row) = table.next()? {
        if let Some(col) = sampleset_column {
            if input.sampleset.is_some_and(|s| s != row[col]) {
                filtered_rows += 1;
                continue;
            }
            samplesets.insert(row[col].clone());
        }
        let sample = sample_column.map(|i| row[i].as_str()).or(input.sample).unwrap();
        let score = score_column.map(|i| row[i].as_str()).or(input.score).unwrap();
        identifier(sample)?;
        identifier(score)?;
        // Do not conflate families or supplementary PLINK sample identifiers, even across scores.
        let identity: Vec<_> = ["FID", "SID"]
            .iter()
            .filter_map(|name| table.column(name).map(|i| (name.to_string(), row[i].clone())))
            .collect();
        if identities
            .insert(sample.to_owned(), identity.clone())
            .is_some_and(|old| old != identity)
        {
            return invalid!(
                "{}: IID {sample:?} occurs with different FID/SID identities; split or explicitly remap before comparison",
                input.path.display()
            );
        }
        let evidence = table
            .header
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != vcol)
            .map(|(i, k)| (k.clone(), row[i].clone()))
            .collect();
        let value = Value {
            text: row[vcol].clone(),
            number: number(&row[vcol])
                .map_err(|e| Error::Invalid(format!("{}:{}: {e}", input.path.display(), table.line)))?,
            evidence,
        };
        if values.insert((sample.to_owned(), score.to_owned()), value).is_some() {
            return invalid!(
                "{}: duplicate sample/score pair {sample:?}/{score:?}",
                input.path.display()
            );
        }
    }
    if samplesets.len() > 1 {
        return invalid!(
            "{}: multiple samplesets; select one explicitly with the corresponding --sampleset option",
            input.path.display()
        );
    }
    if values.is_empty() {
        return invalid!("{}: no score rows selected", input.path.display());
    }
    if before != crate::digest::file_sha256(input.path)? {
        return invalid!("score input changed during comparison: {}", input.path.display());
    }
    let source = Source {
        path: input.path.into(),
        sha256: before,
        format: input.format,
        value_column,
        sample_override: input.sample.map(str::to_owned),
        sampleset: samplesets.into_iter().next(),
        score_id: input.score.map(str::to_owned),
        rows: values.len(),
        filtered_rows,
    };
    Ok((values, source))
}

#[derive(Default, Serialize)]
struct Counts {
    rows: usize,
    exact: usize,
    within_tolerance: usize,
    different: usize,
    left_only: usize,
    right_only: usize,
    unavailable: usize,
    identity_mismatch: usize,
}
impl Counts {
    fn add(&mut self, status: &str) {
        self.rows += 1;
        match status {
            "exact" => self.exact += 1,
            "within_tolerance" => self.within_tolerance += 1,
            "different" => self.different += 1,
            "left_only" => self.left_only += 1,
            "right_only" => self.right_only += 1,
            "identity_mismatch" => self.identity_mismatch += 1,
            _ => self.unavailable += 1,
        }
    }
    fn agrees(&self) -> bool {
        self.rows > 0 && self.exact + self.within_tolerance == self.rows
    }
}
fn classify(left: Option<&Value>, right: Option<&Value>, tolerance: &Tolerance) -> (&'static str, String) {
    match (left, right) {
        (None, _) => ("right_only", String::new()),
        (_, None) => ("left_only", String::new()),
        (Some(a), Some(b)) => {
            for field in ["FID", "SID"] {
                if let (Some(x), Some(y)) = (a.evidence.get(field), b.evidence.get(field))
                    && x != y
                {
                    return ("identity_mismatch", String::new());
                }
            }
            tolerance.classify(a.number.as_ref(), b.number.as_ref())
        }
    }
}

const TERM_FIELDS: &[&str] = &[
    "contig",
    "pos",
    "model",
    "ref",
    "alt",
    "effect_is_alt",
    "weights",
    "status",
    "call_state",
    "effect_dosage",
];
struct Terms {
    table: Table,
    last: u64,
    ordinal: usize,
    contribution: usize,
    sum: Number,
}
impl Terms {
    fn open(path: &Path) -> Result<Self> {
        let table = Table::open(path)?;
        if !table.tabs {
            return invalid!(
                "{}: terms require tab-separated fields (empty fields are significant)",
                path.display()
            );
        }
        for name in TERM_FIELDS {
            table.required(name)?;
        }
        Ok(Self {
            ordinal: table.required("ordinal")?,
            contribution: table.required("contribution")?,
            table,
            last: 0,
            sum: Number::zero(),
        })
    }
    fn next(&mut self) -> Result<Option<(u64, Value)>> {
        let Some(row) = self.table.next()? else {
            return Ok(None);
        };
        let ordinal: u64 = row[self.ordinal]
            .parse()
            .map_err(|_| Error::Invalid("invalid term ordinal".into()))?;
        if ordinal <= self.last {
            return invalid!(
                "{}: term ordinals must be positive, unique and strictly increasing",
                self.table.path.display()
            );
        }
        self.last = ordinal;
        let number = number(&row[self.contribution])?;
        if let Some(n) = &number {
            self.sum = self.sum.add(n);
        }
        let evidence = self
            .table
            .header
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != self.ordinal && *i != self.contribution)
            .map(|(i, k)| (k.clone(), row[i].clone()))
            .collect();
        Ok(Some((
            ordinal,
            Value {
                text: row[self.contribution].clone(),
                number,
                evidence,
            },
        )))
    }
}
fn term_comparison(
    args: &CompareArgs,
    tolerance: &Tolerance,
    left_score: &Value,
    right_score: &Value,
) -> Result<serde_json::Value> {
    let (lp, rp) = (args.left_terms.as_ref().unwrap(), args.right_terms.as_ref().unwrap());
    let (lh, rh) = (crate::digest::file_sha256(lp)?, crate::digest::file_sha256(rp)?);
    let (mut left, mut right) = (Terms::open(lp)?, Terms::open(rp)?);
    let path = args.out.join("terms.jsonl");
    let mut output = BufWriter::new(File::create(&path).map_err(Error::io(&path))?);
    let (mut l, mut r) = (left.next()?, right.next()?);
    let mut counts = Counts::default();
    let mut evidence_changes = 0;
    let mut jointly_unscored = 0;
    while l.is_some() || r.is_some() {
        let ordinal = match (&l, &r) {
            (Some(a), Some(b)) => a.0.min(b.0),
            (Some(a), _) => a.0,
            (_, Some(b)) => b.0,
            _ => unreachable!(),
        };
        let lv = l.as_ref().filter(|a| a.0 == ordinal).map(|a| &a.1);
        let rv = r.as_ref().filter(|a| a.0 == ordinal).map(|a| &a.1);
        let (status, delta) = classify(lv, rv, tolerance);
        counts.add(status);
        let mut changes = Vec::new();
        if let (Some(a), Some(b)) = (lv, rv) {
            let names: BTreeSet<_> = a.evidence.keys().chain(b.evidence.keys()).collect();
            for name in names {
                if a.evidence.get(name) != b.evidence.get(name) {
                    changes.push(name.clone());
                }
            }
            if a.number.is_none() && b.number.is_none() {
                jointly_unscored += 1;
            }
        }
        if !changes.is_empty() {
            evidence_changes += 1;
        }
        let mut explanations = Vec::new();
        if changes
            .iter()
            .any(|f| ["contig", "pos", "ref", "alt", "effect_is_alt", "orientation", "method"].contains(&f.as_str()))
        {
            explanations.push("term identity or allele matching differs; this ordinal may no longer refer to an equivalent scored allele");
        }
        if changes.iter().any(|f| ["weights", "model"].contains(&f.as_str())) {
            explanations.push("published weight or genetic model differs");
        }
        if changes.iter().any(|f| {
            [
                "status",
                "review_reasons",
                "inferred_other_allele",
                "informational_description",
            ]
            .contains(&f.as_str())
        }) {
            explanations.push("term eligibility or recorded policy acceptance differs");
        }
        if changes
            .iter()
            .any(|f| ["call_state", "effect_dosage"].contains(&f.as_str()))
        {
            explanations.push("retained call state or effect-allele dosage differs");
        }
        if matches!(status, "left_only" | "right_only") {
            explanations.push("source ordinal absent from one term export");
        }
        if status == "different" && explanations.is_empty() {
            explanations.push("contribution differs without an explanation in the supplied evidence");
        }
        let record = serde_json::json!({"ordinal": ordinal, "explanations": explanations, "status": status, "delta_right_minus_left": delta,
            "left_contribution": lv.map(|v| &v.text), "right_contribution": rv.map(|v| &v.text),
            "changed_fields": changes, "left_evidence": lv.map(|v| &v.evidence), "right_evidence": rv.map(|v| &v.evidence)});
        serde_json::to_writer(&mut output, &record).map_err(|e| Error::Invalid(e.to_string()))?;
        writeln!(output).map_err(Error::io(&path))?;
        if lv.is_some() {
            l = left.next()?;
        }
        if rv.is_some() {
            r = right.next()?;
        }
    }
    output.flush().map_err(Error::io(&path))?;
    if lh != crate::digest::file_sha256(lp)? || rh != crate::digest::file_sha256(rp)? {
        return invalid!("term input changed during comparison");
    }
    let (ls, ld) = tolerance.classify(Some(&left.sum), left_score.number.as_ref());
    let (rs, rd) = tolerance.classify(Some(&right.sum), right_score.number.as_ref());
    let agrees = counts.rows > 0
        && counts.different + counts.left_only + counts.right_only == 0
        && counts.unavailable == jointly_unscored
        && evidence_changes == 0
        && matches!(ls, "exact" | "within_tolerance")
        && matches!(rs, "exact" | "within_tolerance");
    Ok(
        serde_json::json!({"sample": args.term_sample, "pgs_id": args.term_score,
        "left": {"path": lp, "sha256": lh, "contribution_sum": left.sum.text(), "score_reconciliation": ls, "score_minus_terms": ld},
        "right": {"path": rp, "sha256": rh, "contribution_sum": right.sum.text(), "score_reconciliation": rs, "score_minus_terms": rd},
        "counts": counts, "jointly_unscored": jointly_unscored, "rows_with_changed_evidence": evidence_changes, "agrees": agrees,
        "alignment": "source ordinal; changed term identity fields are reported, never silently harmonized",
        "explanation": "Changed fields are observed evidence, not proof of causality. Blank contributions are omitted from each partial sum; identical unscored terms are not numeric score agreement."}),
    )
}

/// Returns numeric/evidence agreement; the CLI may use it as a validation gate.
pub fn run(args: &CompareArgs) -> Result<bool> {
    if args.left_terms.is_some() != args.right_terms.is_some()
        || args.left_terms.is_some() != args.term_sample.is_some()
        || args.left_terms.is_some() != args.term_score.is_some()
    {
        return invalid!("term comparison requires both term files, term-sample and term-score together");
    }
    let tolerance = Tolerance::new(&args.absolute_tolerance, &args.relative_tolerance)?;
    let (left, ls) = read_scores(Input {
        path: &args.left,
        format: args.left_format,
        sample: args.left_sample.as_deref(),
        column: args.left_score_column.as_deref(),
        score: args.left_score_id.as_deref(),
        sampleset: args.left_sampleset.as_deref(),
    })?;
    let (right, rs) = read_scores(Input {
        path: &args.right,
        format: args.right_format,
        sample: args.right_sample.as_deref(),
        column: args.right_score_column.as_deref(),
        score: args.right_score_id.as_deref(),
        sampleset: args.right_sampleset.as_deref(),
    })?;
    let term_scores = if args.left_terms.is_some() {
        let key = (
            args.term_sample
                .clone()
                .ok_or_else(|| Error::Invalid("term sample is required".into()))?,
            args.term_score
                .clone()
                .ok_or_else(|| Error::Invalid("term score is required".into()))?,
        );
        if args.right_terms.is_none() {
            return invalid!("paired term files are required");
        }
        Some((
            left.get(&key)
                .ok_or_else(|| Error::Invalid("term sample/score absent from left scores".into()))?,
            right
                .get(&key)
                .ok_or_else(|| Error::Invalid("term sample/score absent from right scores".into()))?,
        ))
    } else {
        None
    };
    std::fs::create_dir(&args.out).map_err(Error::io(&args.out))?;
    let result = (|| {
        let path = args.out.join("scores.tsv");
        let mut output = BufWriter::new(File::create(&path).map_err(Error::io(&path))?);
        writeln!(
            output,
            "sample\tpgs_id\tstatus\tleft\tright\tdelta_right_minus_left\tleft_evidence_json\tright_evidence_json"
        )
        .map_err(Error::io(&path))?;
        let keys: BTreeSet<_> = left.keys().chain(right.keys()).collect();
        let mut counts = Counts::default();
        for key in keys {
            let (a, b) = (left.get(key), right.get(key));
            let (status, delta) = classify(a, b, &tolerance);
            counts.add(status);
            writeln!(
                output,
                "{}\t{}\t{status}\t{}\t{}\t{delta}\t{}\t{}",
                key.0,
                key.1,
                a.map_or("", |v| &v.text),
                b.map_or("", |v| &v.text),
                serde_json::to_string(&a.map(|v| &v.evidence)).unwrap(),
                serde_json::to_string(&b.map(|v| &v.evidence)).unwrap()
            )
            .map_err(Error::io(&path))?;
        }
        output.flush().map_err(Error::io(&path))?;
        let terms = term_scores
            .map(|(a, b)| term_comparison(args, &tolerance, a, b))
            .transpose()?;
        let agrees = counts.agrees() && terms.as_ref().is_none_or(|t| t["agrees"] == true);
        let mut outputs = BTreeMap::new();
        outputs.insert("scores.tsv", crate::digest::file_sha256(&path)?);
        if terms.is_some() {
            outputs.insert(
                "terms.jsonl",
                crate::digest::file_sha256(&args.out.join("terms.jsonl"))?,
            );
        }
        let report = serde_json::json!({"schema": "pgsum-comparison-v1", "pgsum_version": env!("CARGO_PKG_VERSION"),
            "left": ls, "right": rs, "absolute_tolerance": args.absolute_tolerance, "relative_tolerance": args.relative_tolerance,
            "counts": counts, "agrees": agrees, "terms": terms, "outputs": outputs,
            "comparison_scope": "Reported numeric values only. Agreement does not establish identical input variants, policies, score publication, assembly or clinical validity.",
            "policies": ["pgsum uses partial_raw_score: observed/expected contributions without frequency filling; withheld strict raw scores remain withheld.",
                "PLINK requires an explicitly named _SUM column. pgsc_calc uses SUM. No averages, denominator rescaling or ancestry-normalized values are substituted.",
                "Sample joins use sample_id/sample or IID. Duplicate sample/score pairs and mixed families/samplesets are rejected; Conflicting FID/SID values on joined rows fail comparison; absent family identifiers cannot be verified.",
                "Missing rows and unavailable values never become zero or agreement. Explicit tolerances cover reported precision only, not biological or policy equivalence.",
                "A numeric discrepancy alone has no proven cause. Inspect term evidence and align missing-call imputation, quality filters, dosage fields, allele matching, ploidy and models."]});
        let marker = args.out.join("comparison.json");
        let temporary = args.out.join("comparison.json.tmp");
        std::fs::write(
            &temporary,
            serde_json::to_vec_pretty(&report).map_err(|e| Error::Invalid(e.to_string()))?,
        )
        .map_err(Error::io(&temporary))?;
        std::fs::rename(&temporary, &marker).map_err(Error::io(&marker))?;
        eprintln!(
            "comparison: {} exact, {} within tolerance, {} different, {} unmatched, {} unavailable, {} identity conflicts",
            counts.exact,
            counts.within_tolerance,
            counts.different,
            counts.left_only + counts.right_only,
            counts.unavailable,
            counts.identity_mismatch
        );
        Ok(agrees)
    })();
    if result.is_err() {
        for name in ["scores.tsv", "terms.jsonl", "comparison.json.tmp", "comparison.json"] {
            let _ = std::fs::remove_file(args.out.join(name));
        }
        let _ = std::fs::remove_dir(&args.out);
    }
    result
}
