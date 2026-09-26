//! Genotype tables: one sample's assessed call at every target of a set of packs, with the gVCF records
//! each call was made from.
//!
//! Layout: the magic bytes `PGSUMGT1`, a little-endian `u64` header length, the JSON header, then one zstd
//! frame with four little-endian sections, each prefixed by its `u64` element count:
//!
//! 1. targets, sorted by key: key `u64`, state `u8`, ALT dosage `u8` (255 = none), flags `u8` (bit 0: RefCall
//!    convention applied, bit 1: phased), first record reference `u32`, record count `u16`;
//! 2. record references `u32` (indices into the records);
//! 3. record end offsets `u64` into the text;
//! 4. record text: each kept gVCF line verbatim, without its line ending, concatenated.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::alleles::Variant;
use crate::digest::hex;
use crate::extract::Extracted;
use crate::genotype::{Call, Policy, State};
use crate::gvcf::HeaderFacts;
use crate::pack::{ReferenceIdentity, SourceFile, records_sha256};
use crate::targets::TargetSet;
use crate::{Error, Result, invalid};
use sha2::{Digest, Sha256};

pub const MAGIC: &[u8; 8] = b"PGSUMGT1";
pub const SCHEMA: &str = "pgsum-genotypes-v3";
/// Earlier schemas still read. v1 has no sequence targets; v1 and v2 pack keys as contig << 40 | pos << 8 |
/// 8 low bits (sequence flag 0x80, so at most 128 sequence targets per position) and are converted on reading.
pub const SCHEMA_V2: &str = "pgsum-genotypes-v2";
pub const SCHEMA_V1: &str = "pgsum-genotypes-v1";

/// Target keys are contig << 48 | pos << 16 | 16 low bits, so keys sort by (contig, pos) and then by the low
/// bits: REF and ALT codes for an SNV, or the sequence flag and an index for a sequence target.
const CONTIG_SHIFT: u32 = 48;
const POS_SHIFT: u32 = 16;

/// Bit set in the key of a sequence (indel or multi-base) target; SNV keys never set it.
pub const SEQUENCE_FLAG: u64 = 0x8000;

/// Distinct sequence targets a position can hold.
pub const MAX_SEQUENCES_PER_POSITION: usize = 0x8000;

/// Key of the `index`-th distinct sequence target at a position.
pub fn sequence_key(contig: u8, pos: u32, index: u16) -> u64 {
    debug_assert!((index as usize) < MAX_SEQUENCES_PER_POSITION);
    (contig as u64) << CONTIG_SHIFT | (pos as u64) << POS_SHIFT | SEQUENCE_FLAG | index as u64
}

pub fn is_sequence_key(key: u64) -> bool {
    key & SEQUENCE_FLAG != 0
}

/// Contig code and position of any target key.
pub fn key_position(key: u64) -> (u8, u32) {
    ((key >> CONTIG_SHIFT) as u8, (key >> POS_SHIFT) as u32)
}

/// A key from a v1 or v2 genotype table in the current layout.
fn upgrade_v2_key(key: u64) -> u64 {
    let (contig, pos, low) = ((key >> 40) as u8, (key >> 8) as u32, key & 0xff);
    let low = if low & 0x80 != 0 { SEQUENCE_FLAG | (low & 0x7f) } else { low };
    (contig as u64) << CONTIG_SHIFT | (pos as u64) << POS_SHIFT | low
}

const BASES: [u8; 4] = *b"ACGT";

fn base_code(b: u8) -> u64 {
    BASES.iter().position(|&x| x == b).expect("resolved bases are ACGT") as u64
}

/// ALT byte of a target that counts any non-reference allele.
pub const ANY_ALT: u8 = b'*';

/// Target key: contig code, position, REF and ALT packed so that keys sort by (contig, pos, REF, ALT). An ALT
/// of `ANY_ALT` sets bit 2 and sorts after the specific ALTs.
pub fn target_key(contig: u8, pos: u32, ref_base: u8, alt_base: u8) -> u64 {
    let alt = if alt_base == ANY_ALT { 4 } else { base_code(alt_base) };
    (contig as u64) << CONTIG_SHIFT | (pos as u64) << POS_SHIFT | base_code(ref_base) << 4 | alt
}

/// `(contig, pos, REF, ALT)` from a target key; ALT is `ANY_ALT` for an any-allele target.
pub fn unpack_key(key: u64) -> (u8, u32, u8, u8) {
    let alt = if key & 4 != 0 {
        ANY_ALT
    } else {
        BASES[(key & 3) as usize]
    };
    (
        (key >> CONTIG_SHIFT) as u8,
        (key >> POS_SHIFT) as u32,
        BASES[(key >> 4 & 3) as usize],
        alt,
    )
}

fn index_sequences(sequences: &[(u64, Variant)]) -> HashMap<(u8, Variant), u64> {
    sequences
        .iter()
        .map(|(k, v)| ((key_position(*k).0, v.clone()), *k))
        .collect()
}

/// An assessed call, as stored.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CompactCall {
    pub state: State,
    pub alt_dosage: Option<u8>,
    pub refcall_adapted: bool,
    pub phased: bool,
}

impl From<Call> for CompactCall {
    fn from(c: Call) -> Self {
        CompactCall {
            state: c.state,
            alt_dosage: c.alt_dosage,
            refcall_adapted: c.refcall_adapted,
            phased: c.phase_set.is_some(),
        }
    }
}

/// Hashes everything written to it.
struct HashingWriter(Sha256);

impl Write for HashingWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0.update(buf);
        Ok(buf.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct PolicyInfo {
    pub id: String,
    pub min_depth: f64,
    pub min_gq: f64,
    pub refcall_is_reference: bool,
    #[serde(default)]
    pub haploid_xy_as_homozygous: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct PackRef {
    pub pgs_id: String,
    pub records_sha256: String,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Header {
    pub schema: String,
    pub pgsum_version: String,
    pub policy: PolicyInfo,
    pub sample: HeaderFacts,
    pub gvcf: SourceFile,
    pub reference: ReferenceIdentity,
    pub packs: Vec<PackRef>,
    pub records_scanned: u64,
    pub records_kept: u64,
    pub targets: u64,
    pub states: BTreeMap<String, u64>,
    /// SHA-256 of the uncompressed body.
    pub body_sha256: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Entry {
    pub key: u64,
    pub state: State,
    pub alt_dosage: Option<u8>,
    pub refcall_adapted: bool,
    pub phased: bool,
    pub first_ref: u32,
    pub ref_count: u16,
}

pub struct GenotypeTable {
    pub header: Header,
    pub entries: Vec<Entry>,
    refs: Vec<u32>,
    offsets: Vec<u64>,
    text: Vec<u8>,
    /// Sequence targets: key and normalized variant, sorted by key.
    sequences: Vec<(u64, Variant)>,
    sequence_index: HashMap<(u8, Variant), u64>,
}

impl GenotypeTable {
    /// Build a table from the scan and one call per target key (`calls[i]` is for `targets.keys[i]`).
    pub fn new(
        gvcf: &Path,
        reference: &ReferenceIdentity,
        policy: &Policy,
        targets: TargetSet,
        scanned: Extracted,
        calls: Vec<CompactCall>,
    ) -> GenotypeTable {
        let TargetSet { keys, sequences, packs } = targets;
        let mut states = BTreeMap::new();
        let mut refs = Vec::with_capacity(keys.len());
        let mut entries = Vec::with_capacity(keys.len());
        for (i, (key, c)) in keys.into_iter().zip(calls).enumerate() {
            *states.entry(c.state.as_str().to_owned()).or_insert(0) += 1;
            let first_ref = refs.len() as u32;
            if scanned.first[i] != u32::MAX {
                refs.push(scanned.first[i]);
                refs.extend(scanned.extra.get(&(i as u32)).into_iter().flatten().copied());
            }
            entries.push(Entry {
                key,
                state: c.state,
                alt_dosage: c.alt_dosage,
                refcall_adapted: c.refcall_adapted,
                phased: c.phased,
                first_ref,
                ref_count: (refs.len() as u32 - first_ref) as u16,
            });
        }
        let mut table = GenotypeTable {
            header: Header {
                schema: SCHEMA.into(),
                pgsum_version: env!("CARGO_PKG_VERSION").into(),
                policy: PolicyInfo {
                    id: policy.id.into(),
                    min_depth: policy.min_depth,
                    min_gq: policy.min_gq,
                    refcall_is_reference: policy.refcall_is_reference,
                    haploid_xy_as_homozygous: policy.haploid_xy_as_homozygous,
                },
                sample: scanned.header,
                gvcf: SourceFile {
                    name: gvcf
                        .file_name()
                        .map(|n| n.to_string_lossy().into_owned())
                        .unwrap_or_default(),
                    sha256: scanned.gvcf_sha256,
                    bytes: scanned.gvcf_bytes,
                },
                reference: reference.clone(),
                packs,
                records_scanned: scanned.records_scanned,
                records_kept: scanned.offsets.len() as u64 - 1,
                targets: entries.len() as u64,
                states,
                body_sha256: String::new(),
            },
            entries,
            refs,
            offsets: scanned.offsets,
            text: scanned.arena,
            sequence_index: index_sequences(&sequences),
            sequences,
        };
        let mut hasher = std::io::BufWriter::with_capacity(1 << 20, HashingWriter(Sha256::new()));
        table.write_body(&mut hasher).expect("hashing cannot fail");
        let hasher = hasher
            .into_inner()
            .map_err(|e| e.into_error())
            .expect("hashing cannot fail");
        table.header.body_sha256 = format!("sha256:{}", hex(&hasher.0.finalize()));
        table
    }

    /// Stream the body sections.
    fn write_body(&self, out: &mut impl Write) -> std::io::Result<()> {
        out.write_all(&(self.entries.len() as u64).to_le_bytes())?;
        for e in &self.entries {
            out.write_all(&e.key.to_le_bytes())?;
            out.write_all(&[
                e.state as u8,
                e.alt_dosage.unwrap_or(255),
                e.refcall_adapted as u8 | (e.phased as u8) << 1,
            ])?;
            out.write_all(&e.first_ref.to_le_bytes())?;
            out.write_all(&e.ref_count.to_le_bytes())?;
        }
        out.write_all(&(self.refs.len() as u64).to_le_bytes())?;
        for r in &self.refs {
            out.write_all(&r.to_le_bytes())?;
        }
        out.write_all(&(self.offsets.len() as u64).to_le_bytes())?;
        for o in &self.offsets {
            out.write_all(&o.to_le_bytes())?;
        }
        out.write_all(&(self.text.len() as u64).to_le_bytes())?;
        out.write_all(&self.text)?;
        out.write_all(&(self.sequences.len() as u64).to_le_bytes())?;
        for (key, v) in &self.sequences {
            out.write_all(&key.to_le_bytes())?;
            out.write_all(&(v.pos as u32).to_le_bytes())?;
            for a in [&v.ref_allele, &v.alt] {
                out.write_all(&(a.len() as u16).to_le_bytes())?;
                out.write_all(a.as_bytes())?;
            }
        }
        Ok(())
    }

    pub fn write(&self, path: &Path) -> Result<()> {
        let tmp = path.with_extension("tmp");
        let json = serde_json::to_vec(&self.header).map_err(|e| Error::Invalid(e.to_string()))?;
        let result = (|| -> std::io::Result<()> {
            let mut out = BufWriter::new(File::create(&tmp)?);
            out.write_all(MAGIC)?;
            out.write_all(&(json.len() as u64).to_le_bytes())?;
            out.write_all(&json)?;
            let mut zstd = zstd::Encoder::new(out, 6)?;
            zstd.multithread(std::thread::available_parallelism().map_or(1, |n| n.get()) as u32)?;
            let mut encoder = std::io::BufWriter::with_capacity(1 << 20, zstd);
            self.write_body(&mut encoder)?;
            encoder
                .into_inner()
                .map_err(|e| e.into_error())?
                .finish()?
                .into_inner()
                .map_err(|e| e.into_error())?
                .sync_all()?;
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

    pub fn open(path: &Path) -> Result<GenotypeTable> {
        let bad = || Error::Invalid(format!("{}: invalid genotype table", path.display()));
        let mut input = BufReader::new(File::open(path).map_err(Error::io(path))?);
        let mut magic = [0u8; 8];
        let mut len = [0u8; 8];
        input.read_exact(&mut magic).map_err(Error::io(path))?;
        if &magic != MAGIC {
            return invalid!("{}: not a pgsum genotype table", path.display());
        }
        input.read_exact(&mut len).map_err(Error::io(path))?;
        let len = u64::from_le_bytes(len);
        if len > 64 << 20 {
            return Err(bad());
        }
        let mut json = vec![0; len as usize];
        input.read_exact(&mut json).map_err(Error::io(path))?;
        let header: Header =
            serde_json::from_slice(&json).map_err(|e| Error::Invalid(format!("{}: {e}", path.display())))?;
        if ![SCHEMA, SCHEMA_V2, SCHEMA_V1].contains(&header.schema.as_str()) {
            return invalid!("{}: schema {} is not {SCHEMA}", path.display(), header.schema);
        }
        let mut body = Vec::new();
        zstd::Decoder::new(input)
            .and_then(|mut d| d.read_to_end(&mut body))
            .map_err(Error::io(path))?;
        if records_sha256(&body) != header.body_sha256 {
            return invalid!("{}: body differs from its digest", path.display());
        }
        let mut at = 0usize;
        let mut take = |n: usize| -> Result<&[u8]> {
            let s = body.get(at..at + n).ok_or_else(bad)?;
            at += n;
            Ok(s)
        };
        let count = |b: &[u8]| u64::from_le_bytes(b.try_into().expect("8 bytes")) as usize;
        let legacy = header.schema != SCHEMA;
        let upgrade = |key: u64| if legacy { upgrade_v2_key(key) } else { key };
        let n = count(take(8)?);
        let mut entries = Vec::with_capacity(n);
        for _ in 0..n {
            let b = take(17)?;
            entries.push(Entry {
                key: upgrade(u64::from_le_bytes(b[0..8].try_into().expect("8"))),
                state: State::from_code(b[8]).ok_or_else(bad)?,
                alt_dosage: (b[9] != 255).then_some(b[9]),
                refcall_adapted: b[10] & 1 != 0,
                phased: b[10] & 2 != 0,
                first_ref: u32::from_le_bytes(b[11..15].try_into().expect("4")),
                ref_count: u16::from_le_bytes(b[15..17].try_into().expect("2")),
            });
        }
        let n = count(take(8)?);
        let refs = take(n * 4)?
            .chunks_exact(4)
            .map(|c| u32::from_le_bytes(c.try_into().expect("4")))
            .collect();
        let n = count(take(8)?);
        let offsets = take(n * 8)?
            .chunks_exact(8)
            .map(|c| u64::from_le_bytes(c.try_into().expect("8")))
            .collect();
        let n = count(take(8)?);
        let text = take(n)?.to_vec();
        let mut sequences = Vec::new();
        if header.schema != SCHEMA_V1 {
            let n = count(take(8)?);
            for _ in 0..n {
                let key = upgrade(u64::from_le_bytes(take(8)?.try_into().expect("8")));
                let pos = u32::from_le_bytes(take(4)?.try_into().expect("4")) as u64;
                let mut alleles = [String::new(), String::new()];
                for a in &mut alleles {
                    let len = u16::from_le_bytes(take(2)?.try_into().expect("2")) as usize;
                    *a = String::from_utf8(take(len)?.to_vec()).map_err(|_| bad())?;
                }
                let [ref_allele, alt] = alleles;
                sequences.push((key, Variant { pos, ref_allele, alt }));
            }
        }
        Ok(GenotypeTable {
            header,
            entries,
            refs,
            offsets,
            text,
            sequence_index: index_sequences(&sequences),
            sequences,
        })
    }

    /// The entry for a sequence target, if the table has it.
    pub fn get_sequence(&self, contig: u8, variant: &Variant) -> Option<&Entry> {
        self.sequence_index
            .get(&(contig, variant.clone()))
            .and_then(|&key| self.get(key))
    }

    /// The normalized variant of a sequence target key.
    pub fn sequence(&self, key: u64) -> Option<&Variant> {
        self.sequences
            .binary_search_by_key(&key, |(k, _)| *k)
            .ok()
            .map(|i| &self.sequences[i].1)
    }

    /// The entry for a target, if the table has it.
    pub fn get(&self, key: u64) -> Option<&Entry> {
        self.entries
            .binary_search_by_key(&key, |e| e.key)
            .ok()
            .map(|i| &self.entries[i])
    }

    /// The gVCF lines an entry was assessed from.
    pub fn records(&self, entry: &Entry) -> impl Iterator<Item = &str> {
        let refs = &self.refs[entry.first_ref as usize..entry.first_ref as usize + entry.ref_count as usize];
        refs.iter().map(|&id| {
            let (a, b) = (
                self.offsets[id as usize] as usize,
                self.offsets[id as usize + 1] as usize,
            );
            std::str::from_utf8(&self.text[a..b]).unwrap_or("")
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keys_round_trip_and_sort() {
        let a = target_key(1, 69_516_650, b'C', b'T');
        assert_eq!(unpack_key(a), (1, 69_516_650, b'C', b'T'));
        assert!(target_key(1, 5, b'T', b'A') < target_key(1, 6, b'A', b'C'));
        assert!(target_key(1, u32::MAX >> 1, b'T', b'G') < target_key(2, 1, b'A', b'C'));
        let any = target_key(3, 100, b'G', ANY_ALT);
        assert_eq!(unpack_key(any), (3, 100, b'G', ANY_ALT));
        assert!(target_key(3, 100, b'G', b'T') < any && any < target_key(3, 101, b'A', b'C'));
    }

    #[test]
    fn sequence_keys_sort_after_snvs_and_hold_many_per_position() {
        let last = sequence_key(1, 109_219_262, (MAX_SEQUENCES_PER_POSITION - 1) as u16);
        assert!(target_key(1, 109_219_262, b'T', ANY_ALT) < sequence_key(1, 109_219_262, 0));
        assert!(sequence_key(1, 109_219_262, 128) < last && last < target_key(1, 109_219_263, b'A', b'C'));
        assert!(is_sequence_key(last) && !is_sequence_key(target_key(1, 5, b'A', b'C')));
        assert_eq!(key_position(last), (1, 109_219_262));
    }

    #[test]
    fn v2_keys_upgrade() {
        let v2_snv = 1u64 << 40 | 69_516_650u64 << 8 | 1 << 4 | 3;
        assert_eq!(upgrade_v2_key(v2_snv), target_key(1, 69_516_650, b'C', b'T'));
        let v2_sequence = 7u64 << 40 | 1_000u64 << 8 | 0x80 | 5;
        assert_eq!(upgrade_v2_key(v2_sequence), sequence_key(7, 1_000, 5));
    }

    #[test]
    fn state_codes_match_order() {
        for (i, s) in State::ALL.iter().enumerate() {
            assert_eq!(*s as u8 as usize, i);
            assert_eq!(State::from_code(i as u8), Some(*s));
        }
    }
}
