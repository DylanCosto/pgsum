//! Genotype tables: one sample's assessed call at every target of a set of packs, with the gVCF records
//! each call was made from.
//!
//! Layout (v4): the magic bytes `PGSUMGT1`, a little-endian `u64` header length, the JSON header, then zstd
//! frames back to back: the sequence targets, then one frame per block of up to `BLOCK_ENTRIES` targets.
//! The header lists each frame's offset (from the end of the header), compressed length and the SHA-256 of
//! its uncompressed bytes, and each block's first and last key, so a reader decompresses only the blocks it
//! needs. A block holds four little-endian sections, each prefixed by its `u64` element count:
//!
//! 1. targets, sorted by key: key `u64`, state `u8`, ALT dosage `u8` (255 = none), flags `u8` (bit 0: RefCall
//!    convention applied, bit 1: phased), first record reference `u32`, record count `u16`;
//! 2. record references `u32` (indices into this block's records);
//! 3. record end offsets `u64` into the text;
//! 4. record text: each kept gVCF line verbatim, without its line ending, concatenated (a record that
//!    overlaps targets in two blocks is stored in both).
//!
//! The sequence frame holds a `u64` count, then per sequence target its key `u64`, position `u32` and REF
//! and ALT (`u16` length and bytes). `body_sha256` is the SHA-256 of the frame digests, one per line, in file
//! order. Earlier tables (v1–v3) are one frame holding the block sections for every target, then (v2, v3)
//! the sequence targets; they are still read, all at once.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use memmap2::Mmap;
use rayon::prelude::*;
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
pub const SCHEMA: &str = "pgsum-genotypes-v4";
/// Earlier schemas still read. v3 is one frame; v1 has no sequence targets; v1 and v2 pack keys as
/// contig << 40 | pos << 8 | 8 low bits (sequence flag 0x80, so at most 128 sequence targets per position)
/// and are converted on reading.
pub const SCHEMA_V3: &str = "pgsum-genotypes-v3";
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
    let low = if low & 0x80 != 0 {
        SEQUENCE_FLAG | (low & 0x7f)
    } else {
        low
    };
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

pub(crate) fn index_sequences(sequences: &[(u64, Variant)]) -> HashMap<(u8, Variant), u64> {
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

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct PolicyInfo {
    pub id: String,
    pub min_depth: f64,
    pub min_gq: f64,
    pub refcall_is_reference: bool,
    #[serde(default)]
    pub haploid_xy_as_homozygous: bool,
    #[serde(default)]
    pub accept_missing_quality: bool,
    /// Records whose ALTs were all structural symbols were left out (`--skip-structural-alleles`); the ID then
    /// ends in `-skip-structural`.
    #[serde(default)]
    pub skip_structural_alleles: bool,
}

impl PolicyInfo {
    /// Record `--skip-structural-alleles` in the policy.
    pub fn with_structural_skipped(mut self, skipped: bool) -> Self {
        if skipped && !self.skip_structural_alleles {
            self.skip_structural_alleles = true;
            self.id.push_str("-skip-structural");
        }
        self
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct PackRef {
    pub pgs_id: String,
    pub records_sha256: String,
}

/// Where a zstd frame is and what it decompresses to.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct Frame {
    /// Byte offset from the end of the header.
    pub offset: u64,
    pub compressed_bytes: u64,
    /// SHA-256 of the uncompressed bytes.
    pub sha256: String,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct BlockInfo {
    pub first_key: u64,
    pub last_key: u64,
    pub targets: u64,
    pub frame: Frame,
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
    /// v4: SHA-256 of the frame digests, one per line in file order; earlier: of the uncompressed body.
    pub body_sha256: String,
    /// v4: the sequence-target frame.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sequence_frame: Option<Frame>,
    /// v4: the target blocks, in key order.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub blocks: Vec<BlockInfo>,
}

/// Targets per block: small enough that scoring a few scores decompresses little, large enough that the
/// header stays small (about 660 blocks for the whole Catalog).
pub const BLOCK_ENTRIES: usize = 1 << 16;

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

/// Targets with the gVCF records they were assessed from.
#[derive(Default)]
struct Block {
    entries: Vec<Entry>,
    refs: Vec<u32>,
    offsets: Vec<u64>,
    text: Vec<u8>,
}

impl Block {
    fn write(&self, out: &mut impl Write) -> std::io::Result<()> {
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
        out.write_all(&self.text)
    }

    fn bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.entries.len() * 17 + self.refs.len() * 4 + self.text.len() + 64);
        self.write(&mut out).expect("writing to memory");
        out
    }

    /// Parse the block sections from the start of `bytes`; returns the block and the bytes consumed.
    fn read(bytes: &[u8]) -> Option<(Block, usize)> {
        let mut at = 0usize;
        let take = |at: &mut usize, n: usize| -> Option<&[u8]> {
            let s = bytes.get(*at..at.checked_add(n)?)?;
            *at += n;
            Some(s)
        };
        let count = |at: &mut usize| -> Option<usize> {
            usize::try_from(u64::from_le_bytes(take(at, 8)?.try_into().ok()?)).ok()
        };
        let n = count(&mut at)?;
        let raw = take(&mut at, n.checked_mul(17)?)?;
        let mut entries = Vec::with_capacity(n);
        for b in raw.as_chunks::<17>().0 {
            entries.push(Entry {
                key: u64::from_le_bytes(b[0..8].try_into().expect("8")),
                state: State::from_code(b[8])?,
                alt_dosage: (b[9] != 255).then_some(b[9]),
                refcall_adapted: b[10] & 1 != 0,
                phased: b[10] & 2 != 0,
                first_ref: u32::from_le_bytes(b[11..15].try_into().expect("4")),
                ref_count: u16::from_le_bytes(b[15..17].try_into().expect("2")),
            });
        }
        let n = count(&mut at)?;
        let refs = take(&mut at, n.checked_mul(4)?)?
            .as_chunks::<4>()
            .0
            .iter()
            .map(|&c| u32::from_le_bytes(c))
            .collect();
        let n = count(&mut at)?;
        let offsets = take(&mut at, n.checked_mul(8)?)?
            .as_chunks::<8>()
            .0
            .iter()
            .map(|&c| u64::from_le_bytes(c))
            .collect();
        let n = count(&mut at)?;
        let text = take(&mut at, n)?.to_vec();
        let block = Block {
            entries,
            refs,
            offsets,
            text,
        };
        block.consistent().then_some((block, at))
    }

    /// Every reference and offset points inside the block.
    fn consistent(&self) -> bool {
        self.entries.windows(2).all(|w| w[0].key < w[1].key)
            && self
                .entries
                .iter()
                .all(|e| e.first_ref as usize + e.ref_count as usize <= self.refs.len())
            && self.refs.iter().all(|&r| (r as usize) + 1 < self.offsets.len().max(1))
            && self.offsets.windows(2).all(|w| w[0] <= w[1])
            && self.offsets.last().is_none_or(|&o| o as usize <= self.text.len())
    }

    fn get(&self, key: u64) -> Option<&Entry> {
        self.entries
            .binary_search_by_key(&key, |e| e.key)
            .ok()
            .map(|i| &self.entries[i])
    }
}

pub(crate) fn write_sequences(sequences: &[(u64, Variant)], out: &mut impl Write) -> std::io::Result<()> {
    out.write_all(&(sequences.len() as u64).to_le_bytes())?;
    for (key, v) in sequences {
        out.write_all(&key.to_le_bytes())?;
        out.write_all(&(v.pos as u32).to_le_bytes())?;
        for a in [&v.ref_allele, &v.alt] {
            out.write_all(&(a.len() as u16).to_le_bytes())?;
            out.write_all(a.as_bytes())?;
        }
    }
    Ok(())
}

pub(crate) fn read_sequences(bytes: &[u8], upgrade: impl Fn(u64) -> u64) -> Option<Vec<(u64, Variant)>> {
    let mut at = 0usize;
    let mut take = |n: usize| -> Option<&[u8]> {
        let s = bytes.get(at..at.checked_add(n)?)?;
        at += n;
        Some(s)
    };
    let n = usize::try_from(u64::from_le_bytes(take(8)?.try_into().ok()?)).ok()?;
    let mut sequences = Vec::with_capacity(n.min(1 << 24));
    for _ in 0..n {
        let key = upgrade(u64::from_le_bytes(take(8)?.try_into().ok()?));
        let pos = u32::from_le_bytes(take(4)?.try_into().ok()?) as u64;
        let mut alleles = [String::new(), String::new()];
        for a in &mut alleles {
            let len = u16::from_le_bytes(take(2)?.try_into().ok()?) as usize;
            *a = String::from_utf8(take(len)?.to_vec()).ok()?;
        }
        let [ref_allele, alt] = alleles;
        sequences.push((key, Variant { pos, ref_allele, alt }));
    }
    Some(sequences)
}

pub(crate) fn sha256(bytes: &[u8]) -> String {
    format!("sha256:{}", hex(&Sha256::digest(bytes)))
}

/// `body_sha256` of a v4 table: the digest of its frame digests, one per line.
pub(crate) fn digest_of_frames<'a>(digests: impl Iterator<Item = &'a str>) -> String {
    let mut h = Sha256::new();
    for d in digests {
        h.update(d.as_bytes());
        h.update(b"\n");
    }
    format!("sha256:{}", hex(&h.finalize()))
}

/// A memory-mapped v4 table file that blocks are decompressed from on first use.
struct Source {
    path: PathBuf,
    map: Mmap,
    body_start: usize,
}

pub struct GenotypeTable {
    pub header: Header,
    /// First key of each block, for finding a key's block.
    first_keys: Vec<u64>,
    last_keys: Vec<u64>,
    blocks: Vec<OnceLock<std::result::Result<Block, String>>>,
    source: Option<Source>,
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
        for c in &calls {
            *states.entry(c.state.as_str().to_owned()).or_insert(0) += 1;
        }
        let targets = keys.len() as u64;
        let records_kept = scanned.offsets.len() as u64 - 1;
        // Split into blocks, each with its own copy of the records its targets were assessed from.
        let starts: Vec<usize> = (0..keys.len()).step_by(BLOCK_ENTRIES).collect();
        let blocks: Vec<Block> = starts
            .par_iter()
            .map(|&start| {
                let end = (start + BLOCK_ENTRIES).min(keys.len());
                let mut block = Block {
                    offsets: vec![0],
                    ..Block::default()
                };
                let mut local: HashMap<u32, u32> = HashMap::new();
                for i in start..end {
                    let first_ref = block.refs.len() as u32;
                    if scanned.first[i] != u32::MAX {
                        let ids = std::iter::once(scanned.first[i])
                            .chain(scanned.extra.get(&(i as u32)).into_iter().flatten().copied());
                        for id in ids {
                            let l = *local.entry(id).or_insert_with(|| {
                                let (a, b) = (
                                    scanned.offsets[id as usize] as usize,
                                    scanned.offsets[id as usize + 1] as usize,
                                );
                                block.text.extend_from_slice(&scanned.arena[a..b]);
                                block.offsets.push(block.text.len() as u64);
                                (block.offsets.len() - 2) as u32
                            });
                            block.refs.push(l);
                        }
                    }
                    let c = calls[i];
                    block.entries.push(Entry {
                        key: keys[i],
                        state: c.state,
                        alt_dosage: c.alt_dosage,
                        refcall_adapted: c.refcall_adapted,
                        phased: c.phased,
                        first_ref,
                        ref_count: (block.refs.len() as u32 - first_ref) as u16,
                    });
                }
                block
            })
            .collect();
        let infos: Vec<BlockInfo> = blocks
            .par_iter()
            .map(|b| BlockInfo {
                first_key: b.entries.first().map_or(0, |e| e.key),
                last_key: b.entries.last().map_or(0, |e| e.key),
                targets: b.entries.len() as u64,
                frame: Frame {
                    sha256: sha256(&b.bytes()),
                    ..Frame::default()
                },
            })
            .collect();
        let mut sequence_bytes = Vec::new();
        write_sequences(&sequences, &mut sequence_bytes).expect("writing to memory");
        let sequence_frame = Frame {
            sha256: sha256(&sequence_bytes),
            ..Frame::default()
        };
        let body_sha256 = digest_of_frames(
            std::iter::once(sequence_frame.sha256.as_str()).chain(infos.iter().map(|b| b.frame.sha256.as_str())),
        );
        GenotypeTable {
            header: Header {
                schema: SCHEMA.into(),
                pgsum_version: env!("CARGO_PKG_VERSION").into(),
                policy: PolicyInfo {
                    id: policy.id.into(),
                    min_depth: policy.min_depth,
                    min_gq: policy.min_gq,
                    refcall_is_reference: policy.refcall_is_reference,
                    haploid_xy_as_homozygous: policy.haploid_xy_as_homozygous,
                    accept_missing_quality: policy.accept_missing_quality,
                    skip_structural_alleles: false,
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
                records_kept,
                targets,
                states,
                body_sha256,
                sequence_frame: Some(sequence_frame),
                blocks: infos,
            },
            first_keys: blocks.iter().map(|b| b.entries.first().map_or(0, |e| e.key)).collect(),
            last_keys: blocks.iter().map(|b| b.entries.last().map_or(0, |e| e.key)).collect(),
            blocks: blocks.into_iter().map(|b| OnceLock::from(Ok(b))).collect(),
            source: None,
            sequence_index: index_sequences(&sequences),
            sequences,
        }
    }

    /// Every block, loading any not yet in memory.
    fn loaded_blocks(&self) -> Result<Vec<&Block>> {
        (0..self.blocks.len()).map(|i| self.block(i)).collect()
    }

    pub fn write(&self, path: &Path) -> Result<()> {
        let blocks = self.loaded_blocks()?;
        let mut sequence_bytes = Vec::new();
        write_sequences(&self.sequences, &mut sequence_bytes).expect("writing to memory");
        let compress = |bytes: &[u8]| zstd::bulk::compress(bytes, 6);
        let frames: Vec<Vec<u8>> = std::iter::once(compress(&sequence_bytes))
            .chain(blocks.par_iter().map(|b| compress(&b.bytes())).collect::<Vec<_>>())
            .collect::<std::io::Result<_>>()
            .map_err(Error::io(path))?;
        let mut header = self.header.clone();
        let mut offset = 0u64;
        for (i, f) in frames.iter().enumerate() {
            let frame = if i == 0 {
                header.sequence_frame.as_mut().expect("v4 tables have a sequence frame")
            } else {
                &mut header.blocks[i - 1].frame
            };
            frame.offset = offset;
            frame.compressed_bytes = f.len() as u64;
            offset += f.len() as u64;
        }
        let tmp = path.with_extension("tmp");
        let json = serde_json::to_vec(&header).map_err(|e| Error::Invalid(e.to_string()))?;
        let result = (|| -> std::io::Result<()> {
            let mut out = BufWriter::new(File::create(&tmp)?);
            out.write_all(MAGIC)?;
            out.write_all(&(json.len() as u64).to_le_bytes())?;
            out.write_all(&json)?;
            for f in &frames {
                out.write_all(f)?;
            }
            out.into_inner().map_err(|e| e.into_error())?.sync_all()?;
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

    /// Open a table. v4 blocks are decompressed when first needed (see `preload` to load them all at once).
    pub fn open(path: &Path) -> Result<GenotypeTable> {
        let bad = |what: &str| Error::Invalid(format!("{}: invalid genotype table ({what})", path.display()));
        let file = File::open(path).map_err(Error::io(path))?;
        // SAFETY: the table is read-only input; changing the file while pgsum runs is unsupported.
        let map = unsafe { Mmap::map(&file) }.map_err(Error::io(path))?;
        if map.get(..8) != Some(MAGIC.as_slice()) {
            return invalid!("{}: not a pgsum genotype table", path.display());
        }
        let len = map
            .get(8..16)
            .map(|b| u64::from_le_bytes(b.try_into().expect("8")))
            .ok_or_else(|| bad("header"))?;
        if len > 64 << 20 {
            return Err(bad("header length"));
        }
        let body_start = 16 + len as usize;
        let header: Header = serde_json::from_slice(map.get(16..body_start).ok_or_else(|| bad("header"))?)
            .map_err(|e| Error::Invalid(format!("{}: {e}", path.display())))?;
        if header.schema != SCHEMA {
            return GenotypeTable::open_single_frame(path, header, &map[body_start..]);
        }
        let sequence_frame = header.sequence_frame.clone().ok_or_else(|| bad("no sequence frame"))?;
        let digest = digest_of_frames(
            std::iter::once(sequence_frame.sha256.as_str())
                .chain(header.blocks.iter().map(|b| b.frame.sha256.as_str())),
        );
        if digest != header.body_sha256
            || header.blocks.iter().map(|b| b.targets).sum::<u64>() != header.targets
            || header.blocks.windows(2).any(|w| w[0].last_key >= w[1].first_key)
        {
            return Err(bad("block list"));
        }
        let source = Source {
            path: path.to_owned(),
            map,
            body_start,
        };
        let bytes = source.frame(&sequence_frame)?;
        let sequences = read_sequences(&bytes, |k| k).ok_or_else(|| bad("sequence targets"))?;
        Ok(GenotypeTable {
            first_keys: header.blocks.iter().map(|b| b.first_key).collect(),
            last_keys: header.blocks.iter().map(|b| b.last_key).collect(),
            blocks: header.blocks.iter().map(|_| OnceLock::new()).collect(),
            source: Some(source),
            sequence_index: index_sequences(&sequences),
            sequences,
            header,
        })
    }

    /// v1–v3: one frame holding every target, read at once; keys are converted to the current layout.
    fn open_single_frame(path: &Path, header: Header, compressed: &[u8]) -> Result<GenotypeTable> {
        let bad = || Error::Invalid(format!("{}: invalid genotype table", path.display()));
        if ![SCHEMA_V3, SCHEMA_V2, SCHEMA_V1].contains(&header.schema.as_str()) {
            return invalid!("{}: schema {} is not {SCHEMA}", path.display(), header.schema);
        }
        let mut body = Vec::new();
        zstd::Decoder::new(compressed)
            .and_then(|mut d| d.read_to_end(&mut body))
            .map_err(Error::io(path))?;
        if records_sha256(&body) != header.body_sha256 {
            return invalid!("{}: body differs from its digest", path.display());
        }
        let legacy = header.schema != SCHEMA_V3;
        let upgrade = |key: u64| if legacy { upgrade_v2_key(key) } else { key };
        let (mut block, used) = Block::read(&body).ok_or_else(bad)?;
        for e in &mut block.entries {
            e.key = upgrade(e.key);
        }
        let sequences = if header.schema == SCHEMA_V1 {
            Vec::new()
        } else {
            read_sequences(&body[used..], upgrade).ok_or_else(bad)?
        };
        Ok(GenotypeTable {
            first_keys: vec![block.entries.first().map_or(0, |e| e.key)],
            last_keys: vec![block.entries.last().map_or(0, |e| e.key)],
            blocks: vec![OnceLock::from(Ok(block))],
            source: None,
            sequence_index: index_sequences(&sequences),
            sequences,
            header,
        })
    }

    /// Block `i`, decompressing and checking it on first use.
    fn block(&self, i: usize) -> Result<&Block> {
        self.blocks[i]
            .get_or_init(|| {
                let source = self.source.as_ref().ok_or("block not in memory")?;
                let info = &self.header.blocks[i];
                let bytes = source.frame(&info.frame).map_err(|e| e.to_string())?;
                let (block, used) = Block::read(&bytes)
                    .filter(|(b, used)| *used == bytes.len() && b.entries.len() as u64 == info.targets)
                    .ok_or_else(|| format!("{}: invalid genotype table block {i}", source.path.display()))?;
                debug_assert_eq!(used, bytes.len());
                Ok(block)
            })
            .as_ref()
            .map_err(|e| Error::Invalid(e.clone()))
    }

    /// Decompress every block now, in parallel, rather than as scoring reaches them.
    pub fn preload(&self) -> Result<()> {
        (0..self.blocks.len())
            .into_par_iter()
            .try_for_each(|i| self.block(i).map(|_| ()))
    }

    /// The entry for a sequence target, if the table has it.
    pub fn get_sequence(&self, contig: u8, variant: &Variant) -> Result<Option<&Entry>> {
        match self.sequence_index.get(&(contig, variant.clone())) {
            Some(&key) => self.get(key),
            None => Ok(None),
        }
    }

    /// The normalized variant of a sequence target key.
    pub fn sequence(&self, key: u64) -> Option<&Variant> {
        self.sequences
            .binary_search_by_key(&key, |(k, _)| *k)
            .ok()
            .map(|i| &self.sequences[i].1)
    }

    /// The entry for a target, if the table has it.
    pub fn get(&self, key: u64) -> Result<Option<&Entry>> {
        let i = self.first_keys.partition_point(|&k| k <= key);
        if i == 0 || key > self.last_keys[i - 1] {
            return Ok(None);
        }
        Ok(self.block(i - 1)?.get(key))
    }

    /// The gVCF lines the entry for `key` was assessed from.
    pub fn records(&self, key: u64) -> Result<Vec<&str>> {
        let i = self.first_keys.partition_point(|&k| k <= key);
        if i == 0 {
            return Ok(Vec::new());
        }
        let block = self.block(i - 1)?;
        let Some(entry) = block.get(key) else {
            return Ok(Vec::new());
        };
        let refs = &block.refs[entry.first_ref as usize..entry.first_ref as usize + entry.ref_count as usize];
        Ok(refs
            .iter()
            .map(|&id| {
                let (a, b) = (
                    block.offsets[id as usize] as usize,
                    block.offsets[id as usize + 1] as usize,
                );
                std::str::from_utf8(&block.text[a..b]).unwrap_or("")
            })
            .collect())
    }

    /// Every entry in key order, decompressing any blocks not yet loaded.
    pub fn entries(&self) -> Result<impl Iterator<Item = &Entry>> {
        self.preload()?;
        Ok(self.loaded_blocks()?.into_iter().flat_map(|b| b.entries.iter()))
    }

    /// Number of blocks decompressed so far.
    pub fn blocks_loaded(&self) -> usize {
        self.blocks.iter().filter(|b| b.get().is_some()).count()
    }
}

impl Source {
    /// Decompress a frame and check its digest.
    fn frame(&self, frame: &Frame) -> Result<Vec<u8>> {
        let bad = || Error::Invalid(format!("{}: invalid genotype table frame", self.path.display()));
        let start = self.body_start.checked_add(frame.offset as usize).ok_or_else(bad)?;
        let compressed = self
            .map
            .get(start..start.checked_add(frame.compressed_bytes as usize).ok_or_else(bad)?)
            .ok_or_else(bad)?;
        let mut bytes = Vec::new();
        zstd::Decoder::new(compressed)
            .and_then(|mut d| d.read_to_end(&mut bytes))
            .map_err(Error::io(&self.path))?;
        if sha256(&bytes) != frame.sha256 {
            return invalid!("{}: a frame differs from its digest", self.path.display());
        }
        Ok(bytes)
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

    /// A table of three blocks: written and reopened, only the blocks a lookup touches are decompressed, and
    /// every entry and record matches the table it was written from.
    #[test]
    fn blocks_load_on_demand() {
        let n = 2 * BLOCK_ENTRIES + 10;
        let keys: Vec<u64> = (0..n as u32)
            .map(|i| target_key(1, 1_000 + 2 * i, b'A', b'G'))
            .collect();
        // One record per target, the last record shared by the last two targets of block 0 and block 1's first.
        let mut arena = Vec::new();
        let mut offsets = vec![0u64];
        for i in 0..n {
            arena.extend_from_slice(format!("chr1\t{}\t.\tA\tG\t50\tPASS\t.\tGT\t0/1", 1_000 + 2 * i).as_bytes());
            offsets.push(arena.len() as u64);
        }
        let first: Vec<u32> = (0..n as u32).collect();
        let mut extra = HashMap::new();
        extra.insert(BLOCK_ENTRIES as u32, vec![(BLOCK_ENTRIES - 1) as u32]);
        let scanned = Extracted {
            header: HeaderFacts::default(),
            gvcf_sha256: "sha256:00".into(),
            gvcf_bytes: 0,
            records_scanned: n as u64,
            arena,
            offsets,
            first,
            extra,
        };
        let calls = vec![
            CompactCall {
                state: State::ObservedVariant,
                alt_dosage: Some(1),
                refcall_adapted: false,
                phased: false,
            };
            n
        ];
        let targets = TargetSet {
            keys: keys.clone(),
            sequences: Vec::new(),
            packs: Vec::new(),
        };
        let reference = ReferenceIdentity {
            fasta_name: "x.fa".into(),
            fasta_sha256: "sha256:00".into(),
            fai_sha256: "sha256:00".into(),
        };
        let table = GenotypeTable::new(
            Path::new("x.g.vcf.gz"),
            &reference,
            &Policy::for_gvcf(false, false, false),
            targets,
            scanned,
            calls,
        );
        assert_eq!(table.header.blocks.len(), 3);
        let path = std::env::temp_dir().join(format!("pgsum-blocks-test-{}.pgsg", std::process::id()));
        table.write(&path).unwrap();
        let read = GenotypeTable::open(&path).unwrap();
        assert_eq!(read.header.body_sha256, table.header.body_sha256);
        assert_eq!(read.blocks_loaded(), 0);
        let last = keys[n - 1];
        assert_eq!(read.get(last).unwrap(), table.get(last).unwrap());
        assert_eq!(read.blocks_loaded(), 1);
        assert_eq!(read.get(last + 1).unwrap(), None);
        assert_eq!(read.get(target_key(2, 1, b'A', b'C')).unwrap(), None);
        assert_eq!(read.blocks_loaded(), 1);
        let shared = keys[BLOCK_ENTRIES];
        assert_eq!(read.records(shared).unwrap().len(), 2);
        assert_eq!(read.records(shared).unwrap(), table.records(shared).unwrap());
        read.preload().unwrap();
        assert_eq!(read.blocks_loaded(), 3);
        for &k in keys.iter().step_by(997) {
            assert_eq!(read.get(k).unwrap(), table.get(k).unwrap());
            assert_eq!(read.records(k).unwrap(), table.records(k).unwrap());
        }
        std::fs::remove_file(&path).unwrap();
    }

    #[test]
    fn state_codes_match_order() {
        for (i, s) in State::ALL.iter().enumerate() {
            assert_eq!(*s as u8 as usize, i);
            assert_eq!(State::from_code(i as u8), Some(*s));
        }
    }
}
