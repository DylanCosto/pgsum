//! Pack files: one compiled score.
//!
//! Layout: the magic bytes `PGSUMPK1`, a little-endian `u64` header length, the JSON header, then one zstd
//! frame holding every term record in source order. Records are little-endian:
//!
//! | field | type |
//! |---|---|
//! | contig code (0 = unresolved) | u8 |
//! | position (0 = unresolved) | u32 |
//! | model, allele kind | u8, u8 |
//! | flags: bit 0 palindromic, bit 1 effect allele is ALT | u8 |
//! | orientation status, method | u8, u8 |
//! | REF base, ALT base (0 unless resolved) | u8, u8 |
//! | review reasons | u32 bit set |
//! | weight count (1, or 3 for dosage weights) | u8 |
//! | each weight: tag, then payload | see `Weight` |

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

pub const MAGIC: &[u8; 8] = b"PGSUMPK1";
pub const SCHEMA: &str = "pgsum-pack-v1";
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

    fn write(&self, out: &mut Vec<u8>) {
        match self {
            Weight::Invalid => out.push(0),
            Weight::Small {
                negative,
                coefficient,
                exponent,
            } => {
                out.push(if *negative { 3 } else { 1 });
                out.extend_from_slice(&coefficient.to_le_bytes());
                out.extend_from_slice(&exponent.to_le_bytes());
            }
            Weight::Text(text) => {
                out.push(2);
                out.push(text.len() as u8);
                out.extend_from_slice(text.as_bytes());
            }
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

impl TermRecord {
    pub fn write(&self, out: &mut Vec<u8>) {
        out.push(self.contig);
        out.extend_from_slice(&self.pos.to_le_bytes());
        out.push(self.model as u8);
        out.push(self.allele_kind as u8);
        out.push(self.palindromic as u8 | (self.orientation.effect_is_alt as u8) << 1);
        out.push(self.orientation.status as u8);
        out.push(self.orientation.method as u8);
        out.push(self.orientation.ref_base);
        out.push(self.orientation.alt_base);
        out.extend_from_slice(&self.reasons.0.to_le_bytes());
        out.push(self.weights.len() as u8);
        for w in &self.weights {
            w.write(out);
        }
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
        let mut encoder = zstd::Encoder::new(out, 9)?;
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

/// An opened pack: its header and decompressed records.
pub struct Pack {
    pub header: Header,
    records: Vec<u8>,
}

impl Pack {
    pub fn open(path: &Path) -> Result<Pack> {
        let mut input = BufReader::new(File::open(path).map_err(Error::io(path))?);
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
        let mut records = Vec::new();
        zstd::Decoder::new(input)
            .and_then(|mut d| d.read_to_end(&mut records))
            .map_err(Error::io(path))?;
        if records_sha256(&records) != header.records_sha256 {
            return invalid!("{}: pack records differ from their digest", path.display());
        }
        Ok(Pack { header, records })
    }

    pub fn terms(&self) -> TermIter<'_> {
        TermIter {
            bytes: &self.records,
            at: 0,
        }
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

pub struct TermIter<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl TermIter<'_> {
    fn take<const N: usize>(&mut self) -> Result<[u8; N]> {
        let slice = self
            .bytes
            .get(self.at..self.at + N)
            .ok_or_else(|| Error::Invalid("truncated pack record".into()))?;
        self.at += N;
        Ok(slice.try_into().expect("length checked"))
    }

    fn u8(&mut self) -> Result<u8> {
        Ok(self.take::<1>()?[0])
    }

    fn record(&mut self) -> Result<TermRecord> {
        let bad = || Error::Invalid("invalid pack record".into());
        let contig = self.u8()?;
        let pos = u32::from_le_bytes(self.take()?);
        let model = Model::from_code(self.u8()?).ok_or_else(bad)?;
        let allele_kind = AlleleKind::from_code(self.u8()?).ok_or_else(bad)?;
        let flags = self.u8()?;
        let status = Status::from_code(self.u8()?).ok_or_else(bad)?;
        let method = Method::from_code(self.u8()?).ok_or_else(bad)?;
        let (ref_base, alt_base) = (self.u8()?, self.u8()?);
        let reasons = Reasons(u32::from_le_bytes(self.take()?));
        let count = self.u8()?;
        let mut weights = Vec::with_capacity(count as usize);
        for _ in 0..count {
            weights.push(match self.u8()? {
                0 => Weight::Invalid,
                tag @ (1 | 3) => Weight::Small {
                    negative: tag == 3,
                    coefficient: u64::from_le_bytes(self.take()?),
                    exponent: i32::from_le_bytes(self.take()?),
                },
                2 => {
                    let n = self.u8()? as usize;
                    let text = self.bytes.get(self.at..self.at + n).ok_or_else(bad)?;
                    self.at += n;
                    Weight::Text(String::from_utf8(text.to_vec()).map_err(|_| bad())?)
                }
                _ => return Err(bad()),
            });
        }
        Ok(TermRecord {
            contig,
            pos,
            model,
            allele_kind,
            palindromic: flags & 1 != 0,
            orientation: Orientation {
                status,
                method,
                ref_base,
                alt_base,
                effect_is_alt: flags & 2 != 0,
            },
            reasons,
            weights,
        })
    }
}

impl Iterator for TermIter<'_> {
    type Item = Result<TermRecord>;
    fn next(&mut self) -> Option<Self::Item> {
        (self.at < self.bytes.len()).then(|| self.record())
    }
}
