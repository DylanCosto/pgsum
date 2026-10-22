//! GenomeOS PGS evidence batches (`pgsum rows`).
//!
//! GenomeOS keeps per-term PGS evidence as batch documents of 2,000 rows. Each row merges pgsum's evidence row
//! (`evidence::for_each_row`) with the term's source row from GenomeOS's catalog pack (`score.sqlite`), its
//! description and its allele orientation; each batch document is digested, and each batch's contribution to the
//! score checkpoint is summed exactly. GenomeOS's `genomeos_api.pgs_pgsum_rows._build` builds the same batches in
//! Python; this module builds them byte for byte the same, so the Python code stays the reference and can recheck
//! any batch. The rules below mirror, in order:
//!
//! - `pgs_catalog_pack.PgsPack.terms`, `describe`, `number`, `centred_weight_description`, `finngen_identifier`;
//! - `pgs_scoring.orientation`, `weight`, `target_positions` and `CompleteScore.add`;
//! - `pgs_frequency.support` (no frequency evidence installed);
//! - `pgs_pgsum_rows._row`, `_source_query`, `_build`;
//! - `wgs_pgs.evidence_value`, `batch_delta`, `_checkpoint`;
//! - `artifact_primitives.canonical_json_bytes` (Python `json.dumps` with sorted keys, `ensure_ascii=False`,
//!   compact separators, then a newline) and `wgs_constraint_pack.digest`.
//!
//! Python's text rules (`str.strip`, `str.upper`, `\d`, `Decimal`) differ from ASCII ones only for non-ASCII text, so
//! a source line with any non-ASCII byte stops the run as unsupported (exit status 3) and GenomeOS builds that run
//! of batches in Python instead.

use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::path::Path;

use memmap2::Mmap;
use num_bigint::BigInt;
use serde::Deserialize;
use serde_json::{Map, Value, json};
use sha2::{Digest, Sha256};

use crate::decimal::Decimal;
use crate::digest::hex;
use crate::evidence::for_each_row;
use crate::genotypes::GenotypeTable;
use crate::pack::Pack;
use crate::reference::Reference;
use crate::{Error, Result, invalid};

/// What GenomeOS passes for one run of batches (`pgs_pgsum_rows._build`'s task, plus its constants).
#[derive(Deserialize)]
pub struct Task {
    pub first: usize,
    pub batches: usize,
    pub batch_size: usize,
    /// `wgs_pgs.DATA_BYTES - 1024**2`: the bound on a batch's summed row bytes.
    pub row_bytes_bound: usize,
    /// `pgs_pgsum_rows.TERM_PAGE` and `PgsPack.terms`' 16 MiB page bound.
    pub term_page: usize,
    pub page_bytes: usize,
    pub identity: Map<String, Value>,
    pub inputs_sha256: String,
    pub pack_sha256: String,
    #[serde(rename = "static")]
    pub static_fields: Map<String, Value>,
    pub keep: bool,
    /// `pgs_scoring.POLICY`, `PASSING` and `pgs_frequency.POLICY`.
    pub call_policy: String,
    pub passing: Vec<String>,
    pub frequency_policy: Value,
}

impl Task {
    /// The task from GenomeOS's JSON. Its numbers are written back into the digested documents, so each must read back
    /// exactly: a number that does not print to the same text Python wrote (a float serde_json reads one unit in the
    /// last place off, an integer beyond 64 bits) leaves the run to Python, as non-ASCII text does.
    pub fn parse(text: &str) -> Result<Task> {
        for token in number_tokens(text) {
            let value: Value = serde_json::from_str(token).map_err(|e| Error::Invalid(format!("task: {e}")))?;
            let mut printed = Vec::new();
            write_value(&mut printed, &value)?;
            if printed != token.as_bytes() {
                return Err(Error::Invalid(format!(
                    "unsupported: task number {token} does not read back exactly"
                )));
            }
        }
        serde_json::from_str(text).map_err(|e| Error::Invalid(format!("task: {e}")))
    }
}

/// The number tokens of a JSON text (outside strings).
fn number_tokens(text: &str) -> Vec<&str> {
    let bytes = text.as_bytes();
    let (mut tokens, mut i, mut in_string) = (Vec::new(), 0, false);
    while i < bytes.len() {
        let b = bytes[i];
        if in_string {
            match b {
                b'\\' => i += 1,
                b'"' => in_string = false,
                _ => {}
            }
        } else if b == b'"' {
            in_string = true;
        } else if b == b'-' || b.is_ascii_digit() {
            let start = i;
            while i < bytes.len() && matches!(bytes[i], b'-' | b'+' | b'.' | b'e' | b'E' | b'0'..=b'9') {
                i += 1;
            }
            tokens.push(&text[start..i]);
            continue;
        }
        i += 1;
    }
    tokens
}

const FLAGS: [&str; 5] = [
    "is_haplotype",
    "is_diplotype",
    "is_interaction",
    "is_dominant",
    "is_recessive",
];
const MAX_TARGETS: usize = 20_000;

fn canonical_contigs() -> Vec<String> {
    let mut v: Vec<String> = (1..=22).map(|n| format!("chr{n}")).collect();
    v.extend(["chrX", "chrY", "chrM"].map(String::from));
    v
}

fn chrom_of(hm_chr: &str) -> Option<String> {
    match hm_chr {
        "X" => Some("chrX".into()),
        "Y" => Some("chrY".into()),
        "MT" | "M" => Some("chrM".into()),
        _ => {
            let n: u32 = hm_chr.parse().ok()?;
            // `CHROMS` keys are "1".."22" exactly (no signs, no leading zeros).
            ((1..=22).contains(&n) && n.to_string() == hm_chr).then(|| format!("chr{n}"))
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------
// Canonical JSON (`canonical_json_bytes`) and digests.

/// Python's `repr(float)`: the shortest round-trip digits, positional for exponents −4..15, else `d.ddde±XX`.
pub fn py_float(x: f64) -> Result<String> {
    if !x.is_finite() {
        return invalid!("value is not finite canonical JSON");
    }
    if x == 0.0 {
        return Ok(if x.is_sign_negative() { "-0.0" } else { "0.0" }.into());
    }
    let s = format!("{x:e}");
    let (mantissa, exp) = s.split_once('e').expect("LowerExp writes an exponent");
    let exp: i32 = exp.parse().expect("LowerExp exponent");
    let (sign, mantissa) = match mantissa.strip_prefix('-') {
        Some(m) => ("-", m),
        None => ("", mantissa),
    };
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    Ok(if !(-4..16).contains(&exp) {
        let m = if digits.len() > 1 {
            format!("{}.{}", &digits[..1], &digits[1..])
        } else {
            digits
        };
        format!("{sign}{m}e{}{:02}", if exp < 0 { '-' } else { '+' }, exp.abs())
    } else if exp >= 0 {
        let e = exp as usize;
        if digits.len() <= e + 1 {
            format!("{sign}{digits}{}.0", "0".repeat(e + 1 - digits.len()))
        } else {
            format!("{sign}{}.{}", &digits[..e + 1], &digits[e + 1..])
        }
    } else {
        format!("{sign}0.{}{digits}", "0".repeat((-exp - 1) as usize))
    })
}

/// A JSON string as Python writes it with `ensure_ascii=False`: `"` `\` and control characters escaped, the rest
/// as UTF-8. Runs without escapes are copied whole.
fn write_str(out: &mut Vec<u8>, s: &str) {
    out.push(b'"');
    let bytes = s.as_bytes();
    let mut start = 0;
    for (i, &b) in bytes.iter().enumerate() {
        let short: &[u8] = match b {
            b'"' => b"\\\"",
            b'\\' => b"\\\\",
            b'\n' => b"\\n",
            b'\r' => b"\\r",
            b'\t' => b"\\t",
            0x08 => b"\\b",
            0x0c => b"\\f",
            0x00..=0x1f => b"",
            _ => continue,
        };
        out.extend_from_slice(&bytes[start..i]);
        if short.is_empty() {
            out.extend_from_slice(format!("\\u{b:04x}").as_bytes());
        } else {
            out.extend_from_slice(short);
        }
        start = i + 1;
    }
    out.extend_from_slice(&bytes[start..]);
    out.push(b'"');
}

fn write_value(out: &mut Vec<u8>, v: &Value) -> Result<()> {
    match v {
        Value::Null => out.extend_from_slice(b"null"),
        Value::Bool(true) => out.extend_from_slice(b"true"),
        Value::Bool(false) => out.extend_from_slice(b"false"),
        Value::Number(n) => {
            if let Some(i) = n.as_i64() {
                out.extend_from_slice(i.to_string().as_bytes());
            } else if let Some(u) = n.as_u64() {
                out.extend_from_slice(u.to_string().as_bytes());
            } else {
                out.extend_from_slice(py_float(n.as_f64().expect("a JSON number"))?.as_bytes());
            }
        }
        Value::String(s) => write_str(out, s),
        Value::Array(items) => {
            out.push(b'[');
            for (i, item) in items.iter().enumerate() {
                if i > 0 {
                    out.push(b',');
                }
                write_value(out, item)?;
            }
            out.push(b']');
        }
        // serde_json's map is ordered by key (byte order of UTF-8, which is code point order, as Python sorts).
        Value::Object(map) => write_map(out, map, None)?,
    }
    Ok(())
}

/// An object; with `rows`, the value of key `"rows"` is those already-written rows instead of the map's.
fn write_map(out: &mut Vec<u8>, map: &Map<String, Value>, rows: Option<&[Vec<u8>]>) -> Result<()> {
    out.push(b'{');
    for (i, (k, item)) in map.iter().enumerate() {
        if i > 0 {
            out.push(b',');
        }
        write_str(out, k);
        out.push(b':');
        match rows {
            Some(rows) if k == "rows" => {
                out.push(b'[');
                for (j, row) in rows.iter().enumerate() {
                    if j > 0 {
                        out.push(b',');
                    }
                    out.extend_from_slice(row);
                }
                out.push(b']');
            }
            _ => write_value(out, item)?,
        }
    }
    out.push(b'}');
    Ok(())
}

/// `canonical_json_bytes`: compact, sorted, non-ASCII kept, then a newline.
pub fn canonical(v: &Value) -> Result<Vec<u8>> {
    let mut out = Vec::with_capacity(256);
    write_value(&mut out, v)?;
    out.push(b'\n');
    Ok(out)
}

fn sha256_text(bytes: &[u8]) -> String {
    format!("sha256:{}", hex(&Sha256::digest(bytes)))
}

/// `wgs_constraint_pack.digest`.
pub fn digest(v: &Value) -> Result<String> {
    Ok(sha256_text(&canonical(v)?))
}

// ---------------------------------------------------------------------------------------------------------------
// Exact decimals (Python `Decimal` at 2,048 digits: every value here is exact, nothing is rounded).

/// `pgs_catalog_pack.number`: the pattern `[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?`, at most 128 characters, a
/// finite value whose adjusted exponent is in −1000..=100. pgsum's `Decimal::parse` applies the same rules.
fn number(text: &str) -> Option<Decimal> {
    Decimal::parse(text)
}

fn big(d: &Decimal) -> BigInt {
    let magnitude: BigInt = d.coefficient.parse().expect("decimal digits");
    if d.negative { -magnitude } else { magnitude }
}

fn from_big(value: &BigInt, exponent: i64, negative_zero: bool) -> Decimal {
    let negative =
        value.sign() == num_bigint::Sign::Minus || (value.sign() == num_bigint::Sign::NoSign && negative_zero);
    Decimal {
        negative,
        coefficient: value.magnitude().to_string(),
        exponent,
    }
}

fn scaled(value: &BigInt, from: i64, to: i64) -> BigInt {
    if from == to {
        value.clone()
    } else {
        value * BigInt::from(10u8).pow((from - to) as u32)
    }
}

/// `a` compared with `b` by value.
fn compare(a: &Decimal, b: &Decimal) -> std::cmp::Ordering {
    let e = a.exponent.min(b.exponent);
    scaled(&big(a), a.exponent, e).cmp(&scaled(&big(b), b.exponent, e))
}

/// `Decimal × int` for a non-negative integer: the exponent is kept and the sign is the decimal's (−0.25 × 0 is
/// −0.00, as in Python).
fn times(d: &Decimal, k: u64) -> Decimal {
    from_big(&(big(d) * BigInt::from(k)), d.exponent, d.negative)
}

/// `Decimal + Decimal`: the finer exponent; an exact zero sum is positive unless both operands are negative.
#[cfg(test)]
fn plus(a: &Decimal, b: &Decimal) -> Decimal {
    let e = a.exponent.min(b.exponent);
    let sum = scaled(&big(a), a.exponent, e) + scaled(&big(b), b.exponent, e);
    from_big(&sum, e, a.negative && b.negative)
}

/// A running `Decimal` sum that starts at `Decimal(0)`. Its zero is never negative: Python's sum of two zeros is
/// negative only when both are, and the running sum starts positive.
struct Sum {
    coefficient: BigInt,
    exponent: i64,
}

impl Sum {
    fn new() -> Sum {
        Sum {
            coefficient: BigInt::from(0),
            exponent: 0,
        }
    }

    fn add(&mut self, d: &Decimal) {
        let e = self.exponent.min(d.exponent);
        let current = scaled(&self.coefficient, self.exponent, e);
        self.coefficient = current + scaled(&big(d), d.exponent, e);
        self.exponent = e;
    }

    fn to_python_string(&self) -> String {
        from_big(&self.coefficient, self.exponent, false).to_python_string()
    }
}

// ---------------------------------------------------------------------------------------------------------------
// Source rows: `describe`, `orientation`, `weight`.

fn get<'a>(row: &'a Map<String, Value>, key: &str) -> &'a str {
    row.get(key).and_then(Value::as_str).unwrap_or("")
}

/// `str.strip()` on ASCII text (Python's whitespace includes the four information separators).
fn strip(s: &str) -> &str {
    s.trim_matches(|c: char| {
        matches!(
            c,
            ' ' | '\t' | '\n' | '\r' | '\x0b' | '\x0c' | '\x1c' | '\x1d' | '\x1e' | '\x1f'
        )
    })
}

fn literal(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

fn digits(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

fn centred_weight_description(value: &str) -> Option<Value> {
    let rest = value.strip_prefix("weight_type=")?;
    let (weight_type, rest) = ["beta", "OR", "HR", "NR"]
        .iter()
        .find_map(|t| rest.strip_prefix(t).map(|r| (*t, r)))?;
    let centre_text = rest.strip_prefix(";centre=")?;
    // The second group has no sign.
    if centre_text.starts_with(['+', '-']) {
        return None;
    }
    let centre = number(centre_text)?;
    let zero = Decimal::parse("0").expect("zero");
    let two = Decimal::parse("2").expect("two");
    if compare(&centre, &zero).is_lt() || compare(&centre, &two).is_gt() {
        return None;
    }
    Some(json!({ "weight_type": weight_type, "centre": centre_text }))
}

fn finngen_identifier(value: &str, row: &Map<String, Value>) -> bool {
    let Some(rest) = value.strip_prefix("FinnGen_VariantID=chr") else {
        return false;
    };
    let parts: Vec<&str> = rest.split('_').collect();
    let [chrom, pos, a, b] = parts[..] else {
        return false;
    };
    let chrom_ok = matches!(chrom, "X" | "Y") || (digits(chrom) && chrom.len() <= 2);
    if !chrom_ok || !digits(pos) || !literal(a) || !literal(b) {
        return false;
    }
    let chr_name = strip(get(row, "chr_name"));
    let chr_name = chr_name.strip_prefix("chr").unwrap_or(chr_name);
    let (effect, other) = (strip(get(row, "effect_allele")), strip(get(row, "other_allele")));
    let ours: BTreeSet<&str> = [a, b].into();
    let theirs: BTreeSet<&str> = [effect, other].into();
    chrom == chr_name && pos == strip(get(row, "chr_position")) && ours == theirs
}

/// `pgs_catalog_pack.describe`.
pub fn describe(row: &Map<String, Value>) -> Map<String, Value> {
    let mut reasons: Vec<String> = Vec::new();
    let mut flags: HashMap<&str, bool> = HashMap::new();
    for name in FLAGS {
        let raw = get(row, name);
        let upper = raw.to_ascii_uppercase();
        if upper == "TRUE" || upper == "FALSE" {
            flags.insert(name, upper == "TRUE");
        } else if matches!(raw, "" | "." | "NA") {
            flags.insert(name, false);
        } else {
            reasons.push(format!("invalid_model_flag:{name}"));
            flags.insert(name, false);
        }
    }
    for k in ["is_haplotype", "is_diplotype", "is_interaction"] {
        if flags[k] {
            reasons.push(format!("requires_special_model:{k}"));
        }
    }
    if !strip(get(row, "inclusion_criteria")).is_empty() {
        reasons.push("conditional_inclusion_requires_review".into());
    }
    let description = strip(get(row, "variant_description"));
    let centre = if description.is_empty() {
        None
    } else {
        centred_weight_description(description)
    };
    if !description.is_empty() && centre.is_none() && !finngen_identifier(description, row) {
        reasons.push("variant_description_requires_review".into());
    }
    if !strip(get(row, "imputation_method")).is_empty() {
        reasons.push("specific_imputation_method_requires_review".into());
    }
    let (dominant, recessive) = (flags["is_dominant"], flags["is_recessive"]);
    let weights = [
        get(row, "dosage_0_weight"),
        get(row, "dosage_1_weight"),
        get(row, "dosage_2_weight"),
    ];
    let explicit = weights.iter().any(|w| !w.is_empty());
    let model = if explicit {
        "dosage_weights"
    } else if dominant {
        "dominant"
    } else if recessive {
        "recessive"
    } else {
        "additive"
    };
    if dominant && recessive || explicit && (dominant || recessive) {
        reasons.push("conflicting_weight_models".into());
    }
    if explicit && weights.iter().any(|w| number(w).is_none()) {
        reasons.push("invalid_dosage_weights".into());
    }
    if !explicit && number(get(row, "effect_weight")).is_none() {
        reasons.push("invalid_effect_weight".into());
    }
    let chrom = chrom_of(get(row, "hm_chr"));
    let pos_raw = get(row, "hm_pos");
    let mut pos: Option<u64> = (pos_raw.len() <= 10
        && pos_raw.as_bytes().first().is_some_and(|b| (b'1'..=b'9').contains(b))
        && digits(pos_raw))
    .then(|| pos_raw.parse().expect("ten digits fit"));
    if pos.is_some_and(|p| p >= 1 << 31) {
        pos = None;
    }
    if chrom.is_none() || pos.is_none() {
        reasons.push("unresolved_harmonized_position".into());
    }
    for key in ["hm_match_chr", "hm_match_pos"] {
        let v = get(row, key).to_ascii_lowercase();
        if v == "false" {
            reasons.push(format!("harmonization_mismatch:{key}"));
        } else if !matches!(v.as_str(), "" | "true" | "na" | ".") {
            reasons.push(format!("invalid_harmonization_flag:{key}"));
        }
    }
    let (effect, other) = (get(row, "effect_allele"), get(row, "other_allele"));
    if other.is_empty() {
        reasons.push("author_other_allele_missing".into());
    }
    if !literal(effect) || !other.is_empty() && !literal(other) {
        reasons.push("nonliteral_alleles".into());
    }
    if effect == other {
        reasons.push("identical_effect_other_alleles".into());
    }
    let allele_kind = if effect.len() == 1 && other.len() == 1 && literal(effect) && literal(other) && effect != other {
        "literal_snv"
    } else if literal(effect) && literal(other) && effect != other {
        "literal_sequence"
    } else {
        "unresolved"
    };
    let palindromic = matches!((effect, other), ("A", "T") | ("T", "A") | ("C", "G") | ("G", "C"));
    let review: BTreeSet<String> = reasons.into_iter().collect();
    let mut result = Map::new();
    result.insert("contig".into(), chrom.map_or(Value::Null, Value::from));
    result.insert("pos1".into(), pos.map_or(Value::Null, Value::from));
    result.insert("model".into(), model.into());
    result.insert("allele_kind".into(), allele_kind.into());
    result.insert("palindromic".into(), palindromic.into());
    result.insert("review_reasons".into(), review.into_iter().collect::<Vec<_>>().into());
    if let Some(centre) = centre {
        result.insert("centred_weight".into(), centre);
    }
    result
}

/// `pgs_scoring.weight`.
fn weight(row: &Map<String, Value>, description: &Map<String, Value>, dosage: u64) -> Result<Option<Decimal>> {
    if dosage > 2 {
        return invalid!("Invalid PGS effect dosage");
    }
    if description["review_reasons"].as_array().is_some_and(|r| !r.is_empty()) {
        return Ok(None);
    }
    let model = description["model"].as_str().unwrap_or("");
    if model == "dosage_weights" {
        let key = format!("dosage_{dosage}_weight");
        let Some(text) = row.get(&key).and_then(Value::as_str) else {
            return invalid!("source row has no {key}");
        };
        return Ok(number(text));
    }
    let Some(w) = number(get(row, "effect_weight")) else {
        return Ok(None);
    };
    let factor = match model {
        "additive" => dosage,
        "dominant" => u64::from(dosage > 0),
        _ => u64::from(dosage == 2),
    };
    Ok(Some(times(&w, factor)))
}

/// The GRCh38 FASTA as GenomeOS reads it (`normalization._IndexedFasta`): the index's names exactly, upper case,
/// only A, C, G, T and N accepted.
pub struct Fasta {
    entries: HashMap<String, (u64, u64, u64, u64)>,
    map: Mmap,
}

impl Fasta {
    pub fn open(path: &Path) -> Result<Fasta> {
        let fai = format!("{}.fai", path.display());
        let text = std::fs::read_to_string(&fai).map_err(Error::io(&fai))?;
        let mut entries = HashMap::new();
        for line in text.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            let parsed: Option<Vec<u64>> = f
                .get(1..5)
                .and_then(|v| v.iter().map(|x| x.parse().ok()).collect::<Option<_>>());
            let Some([length, offset, bases, width]) = parsed.as_deref().and_then(|v| <[u64; 4]>::try_from(v).ok())
            else {
                return invalid!("{fai}: malformed or duplicate FAI entry");
            };
            if length < 1
                || bases < 1
                || width < bases
                || entries
                    .insert(f[0].to_owned(), (length, offset, bases, width))
                    .is_some()
            {
                return invalid!("{fai}: malformed or duplicate FAI entry");
            }
        }
        if entries.is_empty() {
            return invalid!("{fai}: no contigs");
        }
        let file = File::open(path).map_err(Error::io(path))?;
        let map = unsafe { Mmap::map(&file) }.map_err(Error::io(path))?;
        Ok(Fasta { entries, map })
    }

    fn length(&self, contig: &str) -> Option<u64> {
        self.entries.get(contig).map(|e| e.0)
    }

    fn fetch(&self, contig: &str, start: u64, end: u64) -> Result<String> {
        let Some(&(length, offset, bases, width)) = self.entries.get(contig) else {
            return invalid!("Normalized allele contig is absent from the reference: {contig}");
        };
        if end < start || end > length {
            return invalid!("Reference interval is outside {contig}: {start}-{end}");
        }
        let mut seq = Vec::with_capacity((end - start) as usize);
        let mut cursor = start;
        while cursor < end {
            let (line, in_line) = (cursor / bases, cursor % bases);
            let take = (end - cursor).min(bases - in_line);
            let at = (offset + line * width + in_line) as usize;
            let Some(chunk) = self.map.get(at..at + take as usize) else {
                return invalid!("Reference sequence could not be read for {contig}:{}-{end}", start + 1);
            };
            seq.extend_from_slice(chunk);
            cursor += take;
        }
        if !seq.is_ascii() {
            return invalid!("Reference sequence could not be read for {contig}:{}-{end}", start + 1);
        }
        seq.make_ascii_uppercase();
        if !seq.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N')) {
            return invalid!(
                "Reference sequence contains unsupported bases at {contig}:{}-{end}.",
                start + 1
            );
        }
        Ok(String::from_utf8(seq).expect("ASCII"))
    }
}

fn complement(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// `pgs_scoring.orientation` with the term's description.
fn orientation(row: &Map<String, Value>, d: &Map<String, Value>, fasta: &Fasta) -> Result<Map<String, Value>> {
    let mut result = Map::new();
    for key in ["method", "ref", "alt", "effect_allele", "other_allele"] {
        result.insert(key.into(), Value::Null);
    }
    result.insert("candidates".into(), Value::Array(Vec::new()));
    let with = |mut r: Map<String, Value>, status: &str| {
        r.insert("status".into(), status.into());
        Ok(r)
    };
    let (Some(chrom), Some(pos)) = (d["contig"].as_str(), d["pos1"].as_u64()) else {
        return with(result, "unresolved_position");
    };
    let (effect, other) = (get(row, "effect_allele"), get(row, "other_allele"));
    let kind = d["allele_kind"].as_str().unwrap_or("");
    if kind == "unresolved" {
        return with(result, "unresolved_source_alleles");
    }
    match fasta.length(chrom) {
        Some(length) if pos + effect.len().max(other.len()) as u64 - 1 <= length => {}
        _ => return with(result, "outside_reference"),
    }
    if kind != "literal_snv" {
        return with(result, "sequence_reference_evidence_required");
    }
    let r = fasta.fetch(chrom, pos - 1, pos)?;
    let mut candidates: Vec<Map<String, Value>> = Vec::new();
    for (method, e, o) in [
        ("direct", effect.to_owned(), other.to_owned()),
        ("complement", complement(effect), complement(other)),
    ] {
        if r == e || r == o {
            let alt = if e == r { o.clone() } else { e.clone() };
            let mut c = Map::new();
            c.insert("method".into(), method.into());
            c.insert("ref".into(), r.clone().into());
            c.insert("alt".into(), alt.into());
            c.insert("effect_allele".into(), e.into());
            c.insert("other_allele".into(), o.into());
            candidates.push(c);
        }
    }
    if candidates.len() != 1 {
        let palindromic = d["palindromic"].as_bool().unwrap_or(false);
        result.insert("candidates".into(), candidates.into_iter().map(Value::Object).collect());
        return with(
            result,
            if palindromic {
                "palindromic_orientation_unresolved"
            } else {
                "allele_representation_unresolved"
            },
        );
    }
    for (k, v) in &candidates[0] {
        result.insert(k.clone(), v.clone());
    }
    result.insert("candidates".into(), candidates.into_iter().map(Value::Object).collect());
    with(result, "resolved")
}

/// `pgs_frequency.support` with no frequency evidence installed.
fn support(d: &Map<String, Value>, candidates: usize, policy: &Value) -> Value {
    let palindromic = d["palindromic"].as_bool().unwrap_or(false);
    let status = if candidates != 2 || !palindromic {
        "not_applicable"
    } else {
        "reference_frequency_not_installed"
    };
    json!({
        "policy": policy, "status": status, "population": null, "source_frequency": null,
        "reference_counts": null, "reference_interval": null, "compatible_methods": [], "evidence_sha256": null,
        "selected_method": null, "clinical_use": false, "source_record_sha256": null, "excluded_samples": [],
        "target_genotype_used": false, "experimental": true,
    })
}

// ---------------------------------------------------------------------------------------------------------------
// Rows, batches and checkpoints.

struct Term {
    ordinal: u64,
    raw: Map<String, Value>,
    description: Map<String, Value>,
    source_row_sha256: String,
}

fn truthy(v: Option<&Value>) -> bool {
    match v {
        None | Some(Value::Null) | Some(Value::Bool(false)) => false,
        Some(Value::Bool(true)) => true,
        Some(Value::Number(n)) => n.as_f64() != Some(0.0),
        Some(Value::String(s)) => !s.is_empty(),
        Some(Value::Array(a)) => !a.is_empty(),
        Some(Value::Object(o)) => !o.is_empty(),
    }
}

/// `pgs_pgsum_rows._row`: the evidence row for `term` from pgsum's row, after checking them against each other.
fn build_row(term: Term, ours: Value, fasta: &Fasta, task: &Task) -> Result<Map<String, Value>> {
    let Value::Object(mut ours) = ours else {
        return invalid!("PGS term {}: pgsum row is not an object", term.ordinal);
    };
    let d = &term.description;
    if ours.get("ordinal").and_then(Value::as_u64) != Some(term.ordinal) || truthy(ours.get("pgsum_call_not_rebuilt")) {
        return invalid!("PGS term {}: pgsum row differs from the term", term.ordinal);
    }
    let mut o = orientation(&term.raw, d, fasta)?;
    if o["status"] == "palindromic_orientation_unresolved" {
        let n = o["candidates"].as_array().map_or(0, Vec::len);
        o.insert("frequency_support".into(), support(d, n, &task.frequency_policy));
    }
    let take = |ours: &mut Map<String, Value>, key: &str| ours.remove(key).unwrap_or(Value::Null);
    let status = take(&mut ours, "status");
    let call = take(&mut ours, "call");
    let effect_dosage = take(&mut ours, "effect_dosage");
    let contribution = take(&mut ours, "contribution");
    let review = d["review_reasons"].as_array().is_some_and(|r| !r.is_empty());
    let ok = if review || o["status"] != "resolved" {
        let expected = if review {
            "model_term_requires_review"
        } else {
            "unresolved_orientation"
        };
        status == expected && call.is_null() && effect_dosage.is_null() && contribution.is_null()
    } else {
        let alt = &o["alt"];
        let mut ok = call.is_object()
            && call.get("contig") == d.get("contig")
            && call.get("pos1") == d.get("pos1")
            && call.get("ref") == o.get("ref")
            && call
                .get("alts")
                .and_then(Value::as_array)
                .is_some_and(|a| a.len() == 1 && &a[0] == alt)
            && call.get("policy").and_then(Value::as_str) == Some(task.call_policy.as_str());
        let state = call.get("state").and_then(Value::as_str).unwrap_or("");
        if ok && task.passing.iter().any(|p| p == state) {
            let Some(first) = call["dosages"].get(0).and_then(Value::as_u64) else {
                return invalid!("PGS term {}: pgsum call has no dosage", term.ordinal);
            };
            let dose = if o["effect_allele"] == o["alt"] {
                first
            } else {
                2u64.checked_sub(first).unwrap_or(u64::MAX)
            };
            match weight(&term.raw, d, dose)? {
                None => ok = status == "model_term_requires_review" && effect_dosage.is_null(),
                Some(c) => {
                    ok = status == "scorable_observation"
                        && effect_dosage.as_u64() == Some(dose)
                        && contribution.as_str() == Some(c.to_python_string().as_str())
                }
            }
        } else if ok {
            ok = status.as_str() == Some(state) && effect_dosage.is_null() && contribution.is_null();
        }
        ok
    };
    if !ok {
        return invalid!("PGS term {}: pgsum row does not agree with the term", term.ordinal);
    }
    let mut row = Map::new();
    row.insert("ordinal".into(), term.ordinal.into());
    row.insert("source_row_sha256".into(), term.source_row_sha256.into());
    row.insert("source_records".into(), take(&mut ours, "source_records"));
    row.insert("orientation".into(), Value::Object(o));
    row.insert("status".into(), status);
    row.insert("call".into(), call);
    row.insert("effect_dosage".into(), effect_dosage);
    row.insert("contribution".into(), contribution);
    row.insert("clinical_use".into(), false.into());
    row.insert("source_term".into(), Value::Object(term.raw));
    row.insert("description".into(), Value::Object(term.description));
    Ok(row)
}

/// `wgs_pgs.batch_delta` (the reducer `CompleteScore.add`, then `_checkpoint`), written as JSON with the states in
/// the order they were first met, as Python's `Counter` keeps them.
struct Reducer {
    count: u64,
    total: u64,
    states: Vec<(String, u64)>,
    sum: Sum,
}

impl Reducer {
    /// The reducer's checks, including its independent arithmetic check: the stored contribution against the
    /// weight of the source row. Python describes the row again for that weight; a scorable row's description has
    /// no review reasons, so it equals that fresh description and is used as it is.
    fn add(&mut self, row: &Map<String, Value>) -> Result<()> {
        if row["ordinal"].as_u64() != Some(self.count + 1) {
            return invalid!("PGS terms are missing, duplicated or out of order");
        }
        self.count += 1;
        if self.count > self.total {
            return invalid!("PGS results exceed the source inventory");
        }
        let status = row["status"].as_str().unwrap_or("");
        match self.states.iter_mut().find(|(s, _)| s == status) {
            Some((_, n)) => *n += 1,
            None => self.states.push((status.to_owned(), 1)),
        }
        if status == "scorable_observation" {
            let Some(dose) = row["effect_dosage"].as_u64() else {
                return invalid!("Invalid PGS effect dosage");
            };
            let value = row["contribution"].as_str().and_then(number);
            let (Some(raw), Some(description)) = (row["source_term"].as_object(), row["description"].as_object())
            else {
                return invalid!("Invalid PGS contribution evidence");
            };
            let expected = weight(raw, description, dose)?;
            let call_state = row["call"]["state"].as_str().unwrap_or("");
            let resolved = row["orientation"]["status"] == "resolved";
            let reviewed = description["review_reasons"].as_array().is_some_and(|r| !r.is_empty());
            let (Some(value), Some(expected)) = (value, expected) else {
                return invalid!("Invalid PGS contribution evidence");
            };
            if !resolved
                || !matches!(call_state, "observed_reference" | "observed_variant")
                || reviewed
                || compare(&value, &expected).is_ne()
            {
                return invalid!("Invalid PGS contribution evidence");
            }
            self.sum.add(&value);
        } else if !row["contribution"].is_null() || !row["effect_dosage"].is_null() {
            return invalid!("Unresolved PGS term carries a contribution");
        }
        Ok(())
    }

    fn delta_json(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(format!("{{\"count\":{},\"states\":{{", self.count).as_bytes());
        for (i, (state, n)) in self.states.iter().enumerate() {
            if i > 0 {
                out.push(b',');
            }
            write_str(&mut out, state);
            out.extend_from_slice(format!(":{n}").as_bytes());
        }
        out.extend_from_slice(b"},\"sum\":");
        write_str(&mut out, &self.sum.to_python_string());
        out.extend_from_slice(
            b",\"frequency_supported_terms\":0,\"author_supported_terms\":0,\"author_sequence_supported_terms\":0}",
        );
        out
    }
}

/// The catalog pack's terms (`PgsPack.terms`): pages of `term_page` rows after an ordinal, each page's source lines
/// bounded in bytes.
struct Catalog {
    db: rusqlite::Connection,
    columns: Vec<String>,
    /// `{"columns":[...],"raw_line":` — the start of every term's source-row digest document.
    digest_prefix: Vec<u8>,
}

fn sql_error(e: rusqlite::Error) -> Error {
    Error::Invalid(e.to_string())
}

impl Catalog {
    fn page(&self, after: u64, limit: usize, page_bytes: usize) -> Result<Vec<(u64, String, Option<i64>)>> {
        let sql = "SELECT t.ordinal,t.raw_line,d.copies FROM terms t LEFT JOIN duplicates d USING(raw_line) \
                   WHERE ordinal>? ORDER BY ordinal LIMIT ?";
        let mut stmt = self.db.prepare_cached(sql).map_err(sql_error)?;
        let mut rows = stmt.query((after as i64, limit as i64)).map_err(sql_error)?;
        let mut out = Vec::new();
        let mut size = 0usize;
        while let Some(r) = rows.next().map_err(sql_error)? {
            let ordinal: i64 = r.get(0).map_err(sql_error)?;
            let line: String = r.get(1).map_err(sql_error)?;
            let copies: Option<i64> = r.get(2).map_err(sql_error)?;
            size += line.len();
            if size > page_bytes {
                return invalid!("PGS term page exceeds byte bound; request fewer terms");
            }
            out.push((ordinal as u64, line, copies));
        }
        Ok(out)
    }

    fn term(&self, ordinal: u64, line: &str, copies: Option<i64>) -> Result<Term> {
        if !line.is_ascii() {
            return Err(Error::Invalid(format!(
                "unsupported: non-ASCII source line at term {ordinal}"
            )));
        }
        let mut raw = Map::new();
        for (column, field) in self.columns.iter().zip(line.split('\t')) {
            raw.insert(column.clone(), Value::from(field));
        }
        let mut description = describe(&raw);
        if copies.is_some_and(|c| c != 0)
            && let Some(Value::Array(reasons)) = description.get_mut("review_reasons")
        {
            reasons.push("exact_duplicate_source_term".into());
        }
        let mut doc = self.digest_prefix.clone();
        write_str(&mut doc, line);
        doc.extend_from_slice(b"}\n");
        Ok(Term {
            ordinal,
            raw,
            description,
            source_row_sha256: sha256_text(&doc),
        })
    }
}

/// `_targets(sorted(targets))` digested: canonical contig order, then position.
fn targets_sha256(targets: &BTreeSet<(String, u64)>, contigs: &[String]) -> Result<String> {
    let mut list: Vec<(usize, u64, &str)> = Vec::new();
    for (i, (chrom, pos)) in targets.iter().enumerate() {
        if i >= MAX_TARGETS {
            return invalid!("Module targets exceed the bounded query limit");
        }
        let Some(index) = contigs.iter().position(|c| c == chrom) else {
            return invalid!("Module targets require canonical chromosomes and positive positions");
        };
        if !(1..1u64 << 31).contains(pos) {
            return invalid!("Module targets require canonical chromosomes and positive positions");
        }
        list.push((index, *pos, chrom));
    }
    list.sort();
    let value = Value::Array(list.into_iter().map(|(_, p, c)| json!([c, p])).collect());
    digest(&value)
}

/// Build batches `first..first+batches` of one model and write one JSON line per batch:
/// `{"batch", "result_sha256", "delta", "value"}` (`value` only when the task keeps documents).
#[allow(clippy::too_many_arguments)]
pub fn build(
    pack: &Pack,
    table: &GenotypeTable,
    reference: &Reference,
    fasta: &Fasta,
    catalog_dir: &Path,
    task: &Task,
    out: &mut impl std::io::Write,
) -> Result<u64> {
    let manifest_path = catalog_dir.join("manifest.json");
    let manifest: Value = serde_json::from_slice(&std::fs::read(&manifest_path).map_err(Error::io(&manifest_path))?)
        .map_err(|e| Error::Invalid(e.to_string()))?;
    let pgs_id = manifest["pgs_id"].as_str().unwrap_or("").to_owned();
    let Some(total) = manifest["counts"]["terms"].as_u64() else {
        return invalid!("{pgs_id}: manifest has no term count");
    };
    let columns: Vec<String> = manifest["source_header"]["columns"]
        .as_array()
        .map(|c| c.iter().filter_map(|v| v.as_str().map(String::from)).collect())
        .unwrap_or_default();
    let db_path = catalog_dir.join("score.sqlite");
    let uri = format!("file:{}?mode=ro&immutable=1", db_path.display());
    let db = rusqlite::Connection::open_with_flags(
        uri,
        rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY | rusqlite::OpenFlags::SQLITE_OPEN_URI,
    )
    .map_err(|e| Error::Invalid(format!("{}: {e}", db_path.display())))?;
    let mut digest_prefix = b"{\"columns\":".to_vec();
    write_value(&mut digest_prefix, &json!(columns))?;
    digest_prefix.extend_from_slice(b",\"raw_line\":");
    let catalog = Catalog {
        db,
        columns,
        digest_prefix,
    };

    let size = task.batch_size as u64;
    let first = task.first as u64;
    let last = (first + task.batches as u64).min(total.div_ceil(size));
    let after = first * size;
    let limit = (total.min(last * size)).saturating_sub(after);
    let span = Span {
        total,
        first,
        last,
        after,
        limit,
    };

    // pgsum's evidence rows for the same terms, from a second thread through a bounded queue.
    let policy = task.call_policy.as_str();
    std::thread::scope(|scope| {
        let (tx, rx) = std::sync::mpsc::sync_channel::<Value>(4096);
        let producer = scope.spawn(move || {
            for_each_row(
                pack,
                table,
                reference,
                policy,
                after as usize,
                Some(limit as usize),
                |row| tx.send(row).map_err(|_| Error::Invalid(STOPPED.into())),
            )
        });
        let built = batches(&catalog, rx, fasta, task, &pgs_id, &span, out);
        match producer.join().expect("the evidence thread does not panic") {
            Err(Error::Invalid(m)) if m == STOPPED => built,
            Err(e) => Err(e),
            Ok(_) => built,
        }
    })
}

const STOPPED: &str = "the batch builder stopped reading evidence rows";

struct Span {
    total: u64,
    first: u64,
    last: u64,
    after: u64,
    limit: u64,
}

fn batches(
    catalog: &Catalog,
    evidence: std::sync::mpsc::Receiver<Value>,
    fasta: &Fasta,
    task: &Task,
    pgs_id: &str,
    span: &Span,
    out: &mut impl std::io::Write,
) -> Result<u64> {
    let Span {
        total,
        first,
        last,
        after,
        limit,
    } = *span;
    let size = task.batch_size as u64;
    let contigs = canonical_contigs();
    let mut evidence = evidence.into_iter();
    let mut written = 0u64;
    let mut batch = first;
    let mut reducer: Option<Reducer> = None;
    let mut rows: Vec<Vec<u8>> = Vec::new();
    let mut first_ordinal = 0u64;
    let mut bytes = 0usize;
    let mut targets: BTreeSet<(String, u64)> = BTreeSet::new();
    let mut records: BTreeSet<String> = BTreeSet::new();
    let mut position = after;
    while position < after + limit {
        let page = catalog.page(
            position,
            task.term_page.min((after + limit - position) as usize),
            task.page_bytes,
        )?;
        let Some(&(last_ordinal, _, _)) = page.last() else {
            return invalid!("{pgs_id}: the catalog pack has fewer terms than its manifest");
        };
        for (ordinal, line, copies) in page {
            let term = catalog.term(ordinal, &line, copies)?;
            let Some(ours) = evidence.next() else {
                return invalid!("{pgs_id}: pgsum wrote fewer rows than the model has terms");
            };
            if let (Some(c), Some(p)) = (term.description["contig"].as_str(), term.description["pos1"].as_u64()) {
                targets.insert((c.to_owned(), p));
            }
            let row = build_row(term, ours, fasta, task)?;
            if let Some(list) = row["source_records"].as_array() {
                for r in list {
                    records.insert(r["query_record_sha256"].as_str().unwrap_or("").to_owned());
                }
            }
            let mut encoded = Vec::with_capacity(2048);
            write_map(&mut encoded, &row, None)?;
            bytes += encoded.len() + 1;
            if bytes > task.row_bytes_bound {
                return invalid!("PGS batch evidence exceeds byte bound");
            }
            if rows.is_empty() {
                first_ordinal = ordinal;
                reducer = Some(Reducer {
                    count: batch * size,
                    total,
                    states: Vec::new(),
                    sum: Sum::new(),
                });
            }
            reducer.as_mut().expect("a reducer for the open batch").add(&row)?;
            rows.push(encoded);
            if rows.len() as u64 == size || ordinal == total {
                if first_ordinal != batch * size + 1 {
                    return invalid!("PGS batch terms are out of order");
                }
                let source = if targets.is_empty() {
                    Value::Null
                } else {
                    let mut h = Sha256::new();
                    for r in &records {
                        h.update(r.as_bytes());
                        h.update(b"\n");
                    }
                    let mut q = Map::new();
                    q.insert("schema_version".into(), "genomeos-wgs-module-source-pgsum-v1".into());
                    for (k, v) in &task.static_fields {
                        q.insert(k.clone(), v.clone());
                    }
                    q.insert("targets_sha256".into(), targets_sha256(&targets, &contigs)?.into());
                    q.insert(
                        "source_records_sha256".into(),
                        format!("sha256:{}", hex(&h.finalize())).into(),
                    );
                    q.insert("absence_is_reference".into(), false.into());
                    q.insert("clinical_use".into(), false.into());
                    q.insert("callability_authority".into(), false.into());
                    Value::Object(q)
                };
                let mut value = task.identity.clone();
                value.insert("inputs_sha256".into(), task.inputs_sha256.clone().into());
                value.insert("pgs_id".into(), pgs_id.into());
                value.insert("pack_sha256".into(), task.pack_sha256.clone().into());
                value.insert("batch".into(), batch.into());
                value.insert("rows".into(), Value::Null);
                value.insert("source_query".into(), source);
                value.insert("assembly_projection".into(), Value::Null);
                // The document is written once with the rows already encoded; its digest is `result_sha256`.
                let mut doc = Vec::with_capacity(bytes + 1024);
                write_map(&mut doc, &value, Some(&rows))?;
                doc.push(b'\n');
                let result = sha256_text(&doc);
                let mut line = Vec::new();
                line.extend_from_slice(b"{\"batch\":");
                line.extend_from_slice(batch.to_string().as_bytes());
                line.extend_from_slice(b",\"result_sha256\":");
                write_str(&mut line, &result);
                line.extend_from_slice(b",\"delta\":");
                line.extend_from_slice(&reducer.take().expect("a reducer for the open batch").delta_json());
                line.extend_from_slice(b",\"value\":");
                if task.keep {
                    value.insert("result_sha256".into(), result.into());
                    write_map(&mut line, &value, Some(&rows))?;
                } else {
                    line.extend_from_slice(b"null");
                }
                line.extend_from_slice(b"}\n");
                out.write_all(&line).map_err(Error::io("<rows>"))?;
                written += 1;
                batch += 1;
                rows.clear();
                bytes = 0;
                targets.clear();
                records.clear();
            }
        }
        position = last_ordinal;
    }
    if evidence.next().is_some() {
        return invalid!("{pgs_id}: pgsum wrote more rows than asked for");
    }
    if written != last - first {
        return invalid!("{pgs_id}: built {written} of {} batches", last - first);
    }
    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn task_numbers_must_read_back_exactly() {
        let task = |policy: &str| {
            format!(
                r#"{{"first":0,"batches":1,"batch_size":2000,"row_bytes_bound":1,"term_page":20000,"page_bytes":1,
                "identity":{{"case_id":"a \"1.5\" b"}},"inputs_sha256":"x","pack_sha256":"y","static":{{}},"keep":false,
                "call_policy":"p","passing":[],"frequency_policy":{policy}}}"#
            )
        };
        assert!(Task::parse(&task(r#"{"z":1.959963984540054,"n":400,"low":"0.4"}"#)).is_ok());
        for bad in [
            "6.0753501176115204e-151",
            "10000000000000000000000000",
            "1.0e5",
            "1E+16",
        ] {
            let e = Task::parse(&task(&format!(r#"{{"z":{bad}}}"#))).err().expect("refused");
            assert!(
                matches!(&e, Error::Invalid(m) if m.starts_with("unsupported:")),
                "{bad}: {e}"
            );
        }
        assert_eq!(
            number_tokens(r#"{"a":"-1 \" 2","b":[-3.5e-07,4]}"#),
            vec!["-3.5e-07", "4"]
        );
    }

    #[test]
    fn floats_print_as_python_repr() {
        let cases = [
            (30.0, "30.0"),
            (0.5, "0.5"),
            (1.959963984540054, "1.959963984540054"),
            (1e16, "1e+16"),
            (1e15, "1000000000000000.0"),
            (1e-4, "0.0001"),
            (1e-5, "1e-05"),
            (-2.5e-7, "-2.5e-07"),
            (123456.789, "123456.789"),
            (-0.0, "-0.0"),
            (1.5e300, "1.5e+300"),
        ];
        for (x, want) in cases {
            assert_eq!(py_float(x).unwrap(), want, "{x}");
        }
    }

    #[test]
    fn strings_escape_as_python_json() {
        let v = json!({"b": "x\"\\\n\t\u{1}\u{7f}é\u{8}\u{c}\r", "a": [1, -2, 3.0, null, true]});
        assert_eq!(
            String::from_utf8(canonical(&v).unwrap()).unwrap(),
            "{\"a\":[1,-2,3.0,null,true],\"b\":\"x\\\"\\\\\\n\\t\\u0001\u{7f}é\\b\\f\\r\"}\n"
        );
    }

    #[test]
    fn decimal_arithmetic_keeps_python_exponents_and_signs() {
        let d = |t: &str| Decimal::parse(t).unwrap();
        assert_eq!(times(&d("-0.25"), 0).to_python_string(), "-0.00");
        assert_eq!(times(&d("0.25"), 2).to_python_string(), "0.50");
        assert_eq!(times(&d("1.5E-7"), 2).to_python_string(), "3.0E-7");
        assert_eq!(plus(&d("0"), &d("-0.00")).to_python_string(), "0.00");
        assert_eq!(plus(&d("0.25"), &d("-0.25")).to_python_string(), "0.00");
        assert_eq!(plus(&d("1E+2"), &d("1")).to_python_string(), "101");
        assert!(compare(&d("2.0"), &d("2")).is_eq());
        let mut s = Sum::new();
        assert_eq!(s.to_python_string(), "0");
        for t in ["-0.00", "0.25", "-0.5", "1E+2"] {
            s.add(&d(t));
        }
        assert_eq!(s.to_python_string(), "99.75");
    }

    fn row(pairs: &[(&str, &str)]) -> Map<String, Value> {
        pairs.iter().map(|(k, v)| ((*k).to_owned(), Value::from(*v))).collect()
    }

    #[test]
    fn describe_follows_the_catalog_rules() {
        let base = [
            ("effect_allele", "C"),
            ("other_allele", "A"),
            ("effect_weight", "0.25"),
            ("hm_chr", "1"),
            ("hm_pos", "1"),
        ];
        let d = describe(&row(&base));
        assert_eq!(
            Value::Object(d),
            json!({"contig": "chr1", "pos1": 1, "model": "additive", "allele_kind": "literal_snv",
                   "palindromic": false, "review_reasons": []})
        );
        let mut r = row(&base);
        r.insert("is_dominant".into(), "maybe".into());
        r.insert("hm_match_pos".into(), "False".into());
        r.insert("variant_description".into(), " weight_type=beta;centre=1.5 ".into());
        let d = describe(&r);
        assert_eq!(
            d["review_reasons"],
            json!(["harmonization_mismatch:hm_match_pos", "invalid_model_flag:is_dominant"])
        );
        assert_eq!(d["centred_weight"], json!({"weight_type": "beta", "centre": "1.5"}));
        let mut r = row(&base);
        r.insert("variant_description".into(), "FinnGen_VariantID=chr1_1_A_C".into());
        r.insert("chr_name".into(), "chr1".into());
        r.insert("chr_position".into(), "1".into());
        assert_eq!(describe(&r)["review_reasons"], json!([]));
        r.insert("chr_position".into(), "2".into());
        assert_eq!(
            describe(&r)["review_reasons"],
            json!(["variant_description_requires_review"])
        );
        let pal = describe(&row(&[
            ("effect_allele", "T"),
            ("other_allele", "A"),
            ("effect_weight", "x"),
            ("hm_chr", "MT"),
            ("hm_pos", "0"),
        ]));
        assert_eq!(pal["palindromic"], json!(true));
        assert_eq!(pal["contig"], json!("chrM"));
        assert_eq!(
            pal["review_reasons"],
            json!(["invalid_effect_weight", "unresolved_harmonized_position"])
        );
    }
}
