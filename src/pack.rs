//! Pack files: one compiled score.
//!
//! Layout: the magic bytes `PGSUMPK2`, a little-endian `u64` header length, the JSON header, then one zstd
//! frame holding the terms in source order as columns. Each column is a little-endian `u64` byte length
//! followed by its bytes, in this order:
//!
//! | column | encoding, one entry per term unless noted |
//! |---|---|
//! | contig code (0 = unresolved) | u8 |
//! | position (0 = unresolved) | zigzag varint of the difference from the previous term's position |
//! | model, allele kind | u8, u8 |
//! | flags: bit 0 palindromic, bit 1 effect allele is ALT | u8 |
//! | orientation status, method | u8, u8 |
//! | REF base, ALT base (0 unless resolved) | u8, u8 |
//! | review reasons | varint bit set |
//! | weight count (1, or 3 for dosage weights) | u8 |
//! | weight tag: 0 invalid, 1 positive, 3 negative, 2 text | u8 per weight |
//! | weight coefficient (tags 1 and 3) | varint per such weight |
//! | weight exponent (tags 1 and 3) | zigzag varint per such weight |
//! | weight text (tag 2) | u8 length then the text, per such weight |
//!
//! Columns keep like values together, which compresses about a third better than rows: most of what is left
//! is the weights' significant digits.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::decimal::Decimal;
use crate::digest::hex;
use crate::orient::{Method, Orientation, Status};
use crate::term::{AlleleKind, CONTIGS, Model, Reasons};
use crate::{Error, Result, invalid};

pub const MAGIC: &[u8; 8] = b"PGSUMPK2";
pub const SCHEMA: &str = "pgsum-pack-v2";
pub const EXTENSION: &str = "pgsp";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SourceFile {
    pub name: String,
    pub sha256: String,
    pub bytes: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ReferenceIdentity {
    pub fasta_name: String,
    pub fasta_sha256: String,
    pub fai_sha256: String,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Inventory {
    /// `variants_number` in the scoring file header.
    pub declared_terms: u64,
    /// `variants_number` in the Catalog metadata.
    pub catalog_terms: Option<u64>,
    pub actual_terms: u64,
    /// All three counts agree.
    pub consistent: bool,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct Counts {
    pub allele_kinds: BTreeMap<String, u64>,
    pub models: BTreeMap<String, u64>,
    pub review_reasons: BTreeMap<String, u64>,
    pub orientation: BTreeMap<String, u64>,
    pub palindromic_terms: u64,
    pub exact_duplicate_terms: u64,
    pub centred_terms: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct WeightType {
    pub scoring_file: Option<String>,
    pub catalog: Option<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Header {
    pub schema: String,
    pub pgsum_version: String,
    pub pgs_id: String,
    pub source: SourceFile,
    pub scoring_file_header: Vec<String>,
    /// Descriptive scoring-file header keys that appear more than once; each keeps all its lines above.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub duplicate_header_keys: Vec<String>,
    pub columns: Vec<String>,
    pub catalog_metadata_sha256: String,
    /// The Catalog REST record for the score, as downloaded.
    pub catalog_metadata: serde_json::Value,
    pub license: Option<String>,
    pub matches_publication: Option<bool>,
    pub weight_type: WeightType,
    pub reference: ReferenceIdentity,
    pub inventory: Inventory,
    pub counts: Counts,
    /// SHA-256 of the uncompressed record bytes.
    pub records_sha256: String,
}

/// A weight as stored: exact, with the original text kept only when the coefficient exceeds 19 digits.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Weight {
    Invalid,
    Small {
        negative: bool,
        coefficient: u64,
        exponent: i32,
    },
    Text(String),
}

impl Weight {
    pub fn from_decimal(d: Option<&Decimal>, text: &str) -> Weight {
        let Some(d) = d else { return Weight::Invalid };
        match (d.coefficient.parse::<u64>(), i32::try_from(d.exponent)) {
            (Ok(coefficient), Ok(exponent)) => Weight::Small {
                negative: d.negative,
                coefficient,
                exponent,
            },
            _ => Weight::Text(text.to_owned()),
        }
    }

    pub fn to_decimal(&self) -> Option<Decimal> {
        match self {
            Weight::Invalid => None,
            Weight::Small {
                negative,
                coefficient,
                exponent,
            } => Some(Decimal {
                negative: *negative,
                coefficient: coefficient.to_string(),
                exponent: *exponent as i64,
            }),
            Weight::Text(text) => Decimal::parse(text),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct TermRecord {
    pub contig: u8,
    pub pos: u32,
    pub model: Model,
    pub allele_kind: AlleleKind,
    pub palindromic: bool,
    pub orientation: Orientation,
    pub reasons: Reasons,
    pub weights: Vec<Weight>,
}

fn put_varint(out: &mut Vec<u8>, mut v: u64) {
    loop {
        let byte = (v & 0x7f) as u8;
        v >>= 7;
        if v == 0 {
            out.push(byte);
            return;
        }
        out.push(byte | 0x80);
    }
}

fn zigzag(v: i64) -> u64 {
    ((v << 1) ^ (v >> 63)) as u64
}

fn unzigzag(v: u64) -> i64 {
    (v >> 1) as i64 ^ -((v & 1) as i64)
}

/// Terms as columns, decoded. Built while compiling and read back when a pack is opened.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Columns {
    pub contig: Vec<u8>,
    pub pos: Vec<u32>,
    pub model: Vec<u8>,
    pub allele_kind: Vec<u8>,
    pub flags: Vec<u8>,
    pub status: Vec<u8>,
    pub method: Vec<u8>,
    pub ref_base: Vec<u8>,
    pub alt_base: Vec<u8>,
    pub reasons: Vec<u32>,
    /// Index of each term's first weight in `weights`; one extra entry for the end.
    pub weight_start: Vec<u32>,
    pub weights: Vec<Weight>,
}

impl Columns {
    pub fn len(&self) -> usize {
        self.contig.len()
    }

    pub fn is_empty(&self) -> bool {
        self.contig.is_empty()
    }

    pub fn push(&mut self, t: &TermRecord) {
        if self.weight_start.is_empty() {
            self.weight_start.push(0);
        }
        self.contig.push(t.contig);
        self.pos.push(t.pos);
        self.model.push(t.model as u8);
        self.allele_kind.push(t.allele_kind as u8);
        self.flags
            .push(t.palindromic as u8 | (t.orientation.effect_is_alt as u8) << 1);
        self.status.push(t.orientation.status as u8);
        self.method.push(t.orientation.method as u8);
        self.ref_base.push(t.orientation.ref_base);
        self.alt_base.push(t.orientation.alt_base);
        self.reasons.push(t.reasons.0);
        self.weights.extend(t.weights.iter().cloned());
        self.weight_start.push(self.weights.len() as u32);
    }

    /// The term at `i`.
    pub fn get(&self, i: usize) -> Result<TermRecord> {
        let bad = || Error::Invalid(format!("invalid pack term {}", i + 1));
        let (a, b) = (self.weight_start[i] as usize, self.weight_start[i + 1] as usize);
        Ok(TermRecord {
            contig: self.contig[i],
            pos: self.pos[i],
            model: Model::from_code(self.model[i]).ok_or_else(bad)?,
            allele_kind: AlleleKind::from_code(self.allele_kind[i]).ok_or_else(bad)?,
            palindromic: self.flags[i] & 1 != 0,
            orientation: Orientation {
                status: Status::from_code(self.status[i]).ok_or_else(bad)?,
                method: Method::from_code(self.method[i]).ok_or_else(bad)?,
                ref_base: self.ref_base[i],
                alt_base: self.alt_base[i],
                effect_is_alt: self.flags[i] & 2 != 0,
            },
            reasons: Reasons(self.reasons[i]),
            weights: self.weights[a..b].to_vec(),
        })
    }

    /// The uncompressed pack body.
    pub fn encode(&self) -> Vec<u8> {
        let mut sections: Vec<Vec<u8>> = Vec::with_capacity(14);
        let mut pos = Vec::with_capacity(self.len() * 2);
        let mut previous = 0i64;
        for &p in &self.pos {
            put_varint(&mut pos, zigzag(p as i64 - previous));
            previous = p as i64;
        }
        let mut reasons = Vec::with_capacity(self.len());
        for &r in &self.reasons {
            put_varint(&mut reasons, r as u64);
        }
        let counts: Vec<u8> = self.weight_start.windows(2).map(|w| (w[1] - w[0]) as u8).collect();
        let (mut tags, mut coefficients, mut exponents, mut texts) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for w in &self.weights {
            match w {
                Weight::Invalid => tags.push(0),
                Weight::Small {
                    negative,
                    coefficient,
                    exponent,
                } => {
                    tags.push(if *negative { 3 } else { 1 });
                    put_varint(&mut coefficients, *coefficient);
                    put_varint(&mut exponents, zigzag(*exponent as i64));
                }
                Weight::Text(text) => {
                    tags.push(2);
                    texts.push(text.len() as u8);
                    texts.extend_from_slice(text.as_bytes());
                }
            }
        }
        sections.extend([
            self.contig.clone(),
            pos,
            self.model.clone(),
            self.allele_kind.clone(),
            self.flags.clone(),
            self.status.clone(),
            self.method.clone(),
            self.ref_base.clone(),
            self.alt_base.clone(),
            reasons,
            counts,
            tags,
            coefficients,
            exponents,
            texts,
        ]);
        let mut out = Vec::with_capacity(sections.iter().map(|s| s.len() + 8).sum());
        for section in sections {
            out.extend_from_slice(&(section.len() as u64).to_le_bytes());
            out.extend_from_slice(&section);
        }
        out
    }

    #[cfg(test)]
    fn decode_one(bytes: &[u8]) -> u64 {
        let mut v = 0u64;
        for (i, b) in bytes.iter().enumerate() {
            v |= ((b & 0x7f) as u64) << (7 * i);
        }
        v
    }

    /// Decode a pack body.
    pub fn decode(body: &[u8]) -> Result<Columns> {
        let bad = |what: &str| Error::Invalid(format!("invalid pack body: {what}"));
        let mut at = 0usize;
        let mut sections = Vec::with_capacity(15);
        for _ in 0..15 {
            let len = body.get(at..at + 8).ok_or_else(|| bad("truncated"))?;
            let len = u64::from_le_bytes(len.try_into().expect("8 bytes")) as usize;
            at += 8;
            sections.push(body.get(at..at + len).ok_or_else(|| bad("truncated"))?);
            at += len;
        }
        if at != body.len() {
            return Err(bad("trailing bytes"));
        }
        let n = sections[0].len();
        struct Varints<'a>(&'a [u8], usize);
        impl Varints<'_> {
            fn next(&mut self) -> Option<u64> {
                let mut v = 0u64;
                for shift in (0..64).step_by(7) {
                    let b = *self.0.get(self.1)?;
                    self.1 += 1;
                    v |= ((b & 0x7f) as u64) << shift;
                    if b & 0x80 == 0 {
                        return Some(v);
                    }
                }
                None
            }
        }
        let mut pos = Vec::with_capacity(n);
        let mut it = Varints(sections[1], 0);
        let mut previous = 0i64;
        for _ in 0..n {
            previous += unzigzag(it.next().ok_or_else(|| bad("positions"))?);
            pos.push(u32::try_from(previous).map_err(|_| bad("positions"))?);
        }
        let mut reasons = Vec::with_capacity(n);
        let mut it = Varints(sections[9], 0);
        for _ in 0..n {
            reasons.push(u32::try_from(it.next().ok_or_else(|| bad("reasons"))?).map_err(|_| bad("reasons"))?);
        }
        let counts = sections[10];
        let mut weight_start = Vec::with_capacity(n + 1);
        weight_start.push(0u32);
        for &c in counts {
            weight_start.push(weight_start.last().copied().unwrap_or(0) + c as u32);
        }
        let tags = sections[11];
        let (mut coefficients, mut exponents) = (Varints(sections[12], 0), Varints(sections[13], 0));
        let mut texts = 0usize;
        let mut weights = Vec::with_capacity(tags.len());
        for &tag in tags {
            weights.push(match tag {
                0 => Weight::Invalid,
                1 | 3 => Weight::Small {
                    negative: tag == 3,
                    coefficient: coefficients.next().ok_or_else(|| bad("coefficients"))?,
                    exponent: i32::try_from(unzigzag(exponents.next().ok_or_else(|| bad("exponents"))?))
                        .map_err(|_| bad("exponents"))?,
                },
                2 => {
                    let len = *sections[14].get(texts).ok_or_else(|| bad("texts"))? as usize;
                    let text = sections[14]
                        .get(texts + 1..texts + 1 + len)
                        .ok_or_else(|| bad("texts"))?;
                    texts += 1 + len;
                    Weight::Text(String::from_utf8(text.to_vec()).map_err(|_| bad("texts"))?)
                }
                _ => return Err(bad("weight tag")),
            });
        }
        let fixed = [2usize, 3, 4, 5, 6, 7, 8];
        if fixed.iter().any(|&k| sections[k].len() != n)
            || counts.len() != n
            || weight_start.last().copied() != Some(tags.len() as u32)
        {
            return Err(bad("column lengths differ"));
        }
        Ok(Columns {
            contig: sections[0].to_vec(),
            pos,
            model: sections[2].to_vec(),
            allele_kind: sections[3].to_vec(),
            flags: sections[4].to_vec(),
            status: sections[5].to_vec(),
            method: sections[6].to_vec(),
            ref_base: sections[7].to_vec(),
            alt_base: sections[8].to_vec(),
            reasons,
            weight_start,
            weights,
        })
    }
}

/// Write a pack atomically: to a temporary name in the same directory, then renamed into place.
pub fn write(path: &Path, header: &Header, records: &[u8]) -> Result<()> {
    let tmp = path.with_extension(format!("{EXTENSION}.tmp"));
    let header_json = serde_json::to_vec(header).map_err(|e| Error::Invalid(e.to_string()))?;
    let result = (|| -> std::io::Result<()> {
        let mut out = BufWriter::new(File::create(&tmp)?);
        out.write_all(MAGIC)?;
        out.write_all(&(header_json.len() as u64).to_le_bytes())?;
        out.write_all(&header_json)?;
        let mut encoder = zstd::Encoder::new(out, 3)?;
        encoder.write_all(records)?;
        encoder.finish()?.into_inner().map_err(|e| e.into_error())?.sync_all()?;
        std::fs::rename(&tmp, path)
    })();
    result.map_err(|e| {
        let _ = std::fs::remove_file(&tmp);
        Error::Io {
            path: path.to_owned(),
            source: e,
        }
    })
}

pub fn records_sha256(records: &[u8]) -> String {
    format!("sha256:{}", hex(&Sha256::digest(records)))
}

/// An opened pack: its header and its terms as columns.
pub struct Pack {
    pub header: Header,
    pub columns: Columns,
}

impl Pack {
    /// Read only a pack's header.
    pub fn open_header(path: &Path) -> Result<Header> {
        let mut input = BufReader::new(File::open(path).map_err(Error::io(path))?);
        read_header(&mut input, path)
    }

    pub fn open(path: &Path) -> Result<Pack> {
        let mut input = BufReader::new(File::open(path).map_err(Error::io(path))?);
        let header = read_header(&mut input, path)?;
        let mut body = Vec::new();
        zstd::Decoder::new(input)
            .and_then(|mut d| d.read_to_end(&mut body))
            .map_err(Error::io(path))?;
        if records_sha256(&body) != header.records_sha256 {
            return invalid!("{}: pack records differ from their digest", path.display());
        }
        let columns = Columns::decode(&body).map_err(|e| Error::Invalid(format!("{}: {e}", path.display())))?;
        if columns.len() as u64 != header.inventory.actual_terms {
            return invalid!(
                "{}: pack holds {} terms, header says {}",
                path.display(),
                columns.len(),
                header.inventory.actual_terms
            );
        }
        Ok(Pack { header, columns })
    }

    /// Terms from index `start` on.
    pub fn terms_from(&self, start: usize) -> TermIter<'_> {
        TermIter {
            columns: &self.columns,
            at: start,
        }
    }

    /// `(start index, start index)` of every `every`-th term, starting with the first.
    pub fn chunk_starts(&self, every: usize) -> Result<Vec<(usize, usize)>> {
        Ok((0..self.columns.len()).step_by(every.max(1)).map(|i| (i, i)).collect())
    }

    pub fn terms(&self) -> TermIter<'_> {
        self.terms_from(0)
    }

    /// One line per term, in source order: description, review reasons, orientation and weights.
    pub fn write_terms_tsv(&self, out: &mut impl Write) -> Result<()> {
        out.write_all(TSV_HEADER.as_bytes()).map_err(io_error)?;
        for (i, term) in self.terms().enumerate() {
            write_term_line(out, i + 1, &term?)?;
            writeln!(out).map_err(io_error)?;
        }
        Ok(())
    }
}

fn io_error(e: std::io::Error) -> Error {
    Error::Io {
        path: "<output>".into(),
        source: e,
    }
}

/// The `TSV_HEADER` columns for one term, without a line ending.
pub fn write_term_line(out: &mut (impl Write + ?Sized), ordinal: usize, t: &TermRecord) -> Result<()> {
    let base = |b: u8| if b == 0 { String::new() } else { (b as char).to_string() };
    let o = t.orientation;
    let weights: Vec<String> = t
        .weights
        .iter()
        .map(|w| w.to_decimal().map(|d| d.to_python_string()).unwrap_or_default())
        .collect();
    write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        ordinal,
        if t.contig == 0 {
            ""
        } else {
            CONTIGS[t.contig as usize - 1]
        },
        if t.pos == 0 { String::new() } else { t.pos.to_string() },
        t.model.as_str(),
        t.allele_kind.as_str(),
        t.palindromic as u8,
        t.reasons.names().join(","),
        o.status.as_str(),
        o.method.as_str(),
        base(o.ref_base),
        base(o.alt_base),
        o.effect_is_alt as u8,
        weights.join(","),
    )
    .map_err(io_error)
}

pub const TSV_HEADER: &str = "ordinal\tcontig\tpos\tmodel\tallele_kind\tpalindromic\treview_reasons\torientation\tmethod\tref\talt\teffect_is_alt\tweights\n";

fn read_header(input: &mut impl Read, path: &Path) -> Result<Header> {
    let mut magic = [0u8; 8];
    let mut len = [0u8; 8];
    input.read_exact(&mut magic).map_err(Error::io(path))?;
    if &magic != MAGIC {
        return invalid!("{}: not a pgsum pack", path.display());
    }
    input.read_exact(&mut len).map_err(Error::io(path))?;
    let len = u64::from_le_bytes(len);
    if len > 64 << 20 {
        return invalid!("{}: pack header is too large", path.display());
    }
    let mut header_json = vec![0; len as usize];
    input.read_exact(&mut header_json).map_err(Error::io(path))?;
    let header: Header =
        serde_json::from_slice(&header_json).map_err(|e| Error::Invalid(format!("{}: {e}", path.display())))?;
    if header.schema != SCHEMA {
        return invalid!("{}: pack schema {} is not {SCHEMA}", path.display(), header.schema);
    }
    Ok(header)
}

pub struct TermIter<'a> {
    columns: &'a Columns,
    at: usize,
}

impl Iterator for TermIter<'_> {
    type Item = Result<TermRecord>;
    fn next(&mut self) -> Option<Self::Item> {
        (self.at < self.columns.len()).then(|| {
            self.at += 1;
            self.columns.get(self.at - 1)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn varints_round_trip() {
        for v in [0i64, 1, -1, 63, -64, 64, 1 << 31, -(1 << 31), i64::MAX, i64::MIN] {
            let mut b = Vec::new();
            put_varint(&mut b, zigzag(v));
            assert_eq!(unzigzag(Columns::decode_one(&b)), v);
        }
    }
}
