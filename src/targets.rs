//! The union of target sites across packs, and a cache of it.
//!
//! Collecting targets reads every pack. For the whole Catalog that is ~32 GB, so the union can be saved to a
//! target index (`.pgst`) together with the identity of every pack it came from. A later run checks only the
//! packs' headers (a few KB each) and reuses the index if the pack set and reference are unchanged.
//!
//! Index layout: the magic bytes `PGSUMTI1`, a little-endian `u64` header length, the JSON header, then one
//! zstd frame: the sorted SNV keys as varints of the difference from the previous key, then (v2) the sequence
//! targets, each as contig code `u8`, position `u32` and REF/ALT (`u16` length and bytes), little-endian.
//!
//! Sequence targets get their keys from their sorted order (`genotypes::sequence_key`), so the same set always
//! gets the same keys.

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use std::collections::BTreeSet;

use crate::alleles::Variant;
use crate::genotypes::{MAX_SEQUENCES_PER_POSITION, PackRef, is_sequence_key, key_position, sequence_key};
use crate::pack::{Pack, ReferenceIdentity, records_sha256};
use crate::{Error, Result, invalid};

pub const MAGIC: &[u8; 8] = b"PGSUMTI1";
pub const SCHEMA: &str = "pgsum-targets-v3";

/// The targets of a set of packs.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TargetSet {
    /// SNV and sequence target keys, sorted.
    pub keys: Vec<u64>,
    /// Sequence targets by key, sorted by key.
    pub sequences: Vec<(u64, Variant)>,
    /// The packs, sorted by PGS ID.
    pub packs: Vec<PackRef>,
}

impl TargetSet {
    /// Assign keys to distinct sequence targets and merge them with the SNV keys.
    fn new(mut snv_keys: Vec<u64>, sequences: Vec<(u8, Variant)>, packs: Vec<PackRef>) -> Result<TargetSet> {
        let mut keyed = Vec::with_capacity(sequences.len());
        let mut previous: Option<(u8, u64)> = None;
        let mut index = 0usize;
        for (contig, v) in sequences {
            index = if previous == Some((contig, v.pos)) {
                index + 1
            } else {
                0
            };
            if index >= MAX_SEQUENCES_PER_POSITION {
                return invalid!(
                    "more than {MAX_SEQUENCES_PER_POSITION} distinct indel targets at {}:{}",
                    crate::term::CONTIGS[contig as usize - 1],
                    v.pos
                );
            }
            previous = Some((contig, v.pos));
            keyed.push((sequence_key(contig, v.pos as u32, index as u16), v));
        }
        snv_keys.extend(keyed.iter().map(|(k, _)| *k));
        snv_keys.sort_unstable();
        keyed.sort_by_key(|(k, _)| *k);
        Ok(TargetSet {
            keys: snv_keys,
            sequences: keyed,
            packs,
        })
    }

    fn sequence_list(&self) -> Vec<(u8, Variant)> {
        self.sequences
            .iter()
            .map(|(k, v)| (key_position(*k).0, v.clone()))
            .collect()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Header {
    pub schema: String,
    pub pgsum_version: String,
    pub reference: ReferenceIdentity,
    /// The packs the targets came from, sorted by PGS ID.
    pub packs: Vec<PackRef>,
    /// SNV targets.
    pub targets: u64,
    /// Sequence targets.
    #[serde(default)]
    pub sequence_targets: u64,
    /// SHA-256 of the uncompressed body.
    pub body_sha256: String,
    /// The keys include a position target for every term with a position.
    #[serde(default)]
    pub term_positions: bool,
}

/// Where the targets came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Source {
    Packs,
    Cache,
}

fn sorted_refs(mut refs: Vec<PackRef>) -> Vec<PackRef> {
    refs.sort_by(|a, b| a.pgs_id.cmp(&b.pgs_id));
    refs
}

fn check_reference(pgs_id: &str, pack: &ReferenceIdentity, reference: &ReferenceIdentity) -> Result<()> {
    if pack != reference {
        return invalid!(
            "{pgs_id} was compiled against {} ({}), not this reference",
            pack.fasta_name,
            pack.fasta_sha256
        );
    }
    Ok(())
}

/// Read every pack (in parallel, target columns only) and return its targets and identity.
pub fn collect(packs: &[PathBuf], reference: &ReferenceIdentity, positions: bool) -> Result<TargetSet> {
    // SNV keys from all packs accumulate here and are deduplicated whenever they have doubled.
    let union: Mutex<(Vec<u64>, usize)> = Mutex::new((Vec::new(), 0));
    let sequences: Mutex<BTreeSet<(u8, Variant)>> = Mutex::new(BTreeSet::new());
    let refs: Vec<PackRef> = packs
        .par_iter()
        .map(|path| -> Result<PackRef> {
            let (header, keys, seqs) = Pack::open_target_keys(path, positions)?;
            check_reference(&header.pgs_id, &header.reference, reference)?;
            sequences.lock().expect("no panics while holding the lock").extend(seqs);
            let mut guard = union.lock().expect("no panics while holding the lock");
            let (all, deduplicated) = &mut *guard;
            all.extend_from_slice(&keys);
            if all.len() > 2 * *deduplicated + (1 << 22) {
                all.sort_unstable();
                all.dedup();
                *deduplicated = all.len();
            }
            Ok(PackRef {
                pgs_id: header.pgs_id,
                records_sha256: header.records_sha256,
            })
        })
        .collect::<Result<_>>()?;
    let mut keys = union.into_inner().expect("no panics while holding the lock").0;
    keys.sort_unstable();
    keys.dedup();
    let sequences = sequences
        .into_inner()
        .expect("no panics while holding the lock")
        .into_iter()
        .collect();
    TargetSet::new(keys, sequences, sorted_refs(refs))
}

/// The identities of the packs, read from their headers only.
pub fn pack_refs(packs: &[PathBuf], reference: &ReferenceIdentity) -> Result<Vec<PackRef>> {
    let refs = packs
        .par_iter()
        .map(|path| -> Result<PackRef> {
            let header = Pack::open_header(path)?;
            check_reference(&header.pgs_id, &header.reference, reference)?;
            Ok(PackRef {
                pgs_id: header.pgs_id,
                records_sha256: header.records_sha256,
            })
        })
        .collect::<Result<_>>()?;
    Ok(sorted_refs(refs))
}

fn encode(keys: &[u64]) -> Vec<u8> {
    let mut out = Vec::with_capacity(keys.len() * 2);
    let mut previous = 0u64;
    for &k in keys {
        let mut v = k - previous;
        previous = k;
        loop {
            let byte = (v & 0x7f) as u8;
            v >>= 7;
            if v == 0 {
                out.push(byte);
                break;
            }
            out.push(byte | 0x80);
        }
    }
    out
}

/// The first `n` delta-encoded keys and the offset after them.
fn decode(body: &[u8], n: usize) -> Result<(Vec<u64>, usize)> {
    let bad = || Error::Invalid("invalid target index body".into());
    let mut keys = Vec::with_capacity(n);
    let (mut at, mut previous) = (0usize, 0u64);
    while keys.len() < n {
        let mut v = 0u64;
        let mut shift = 0;
        loop {
            let b = *body.get(at).ok_or_else(bad)?;
            at += 1;
            v |= ((b & 0x7f) as u64) << shift;
            if b & 0x80 == 0 {
                break;
            }
            shift += 7;
            if shift > 63 {
                return Err(bad());
            }
        }
        previous = previous.checked_add(v).ok_or_else(bad)?;
        keys.push(previous);
    }
    Ok((keys, at))
}

pub fn write(path: &Path, reference: &ReferenceIdentity, set: &TargetSet, positions: bool) -> Result<()> {
    let snv: Vec<u64> = set.keys.iter().copied().filter(|&k| !is_sequence_key(k)).collect();
    let mut body = encode(&snv);
    let sequences = set.sequence_list();
    for (contig, v) in &sequences {
        body.push(*contig);
        body.extend_from_slice(&(v.pos as u32).to_le_bytes());
        for a in [&v.ref_allele, &v.alt] {
            body.extend_from_slice(&(a.len() as u16).to_le_bytes());
            body.extend_from_slice(a.as_bytes());
        }
    }
    let header = Header {
        schema: SCHEMA.into(),
        pgsum_version: env!("CARGO_PKG_VERSION").into(),
        reference: reference.clone(),
        packs: set.packs.clone(),
        targets: snv.len() as u64,
        sequence_targets: sequences.len() as u64,
        body_sha256: records_sha256(&body),
        term_positions: positions,
    };
    let json = serde_json::to_vec(&header).map_err(|e| Error::Invalid(e.to_string()))?;
    let tmp = path.with_extension("pgst.tmp");
    let result = (|| -> std::io::Result<()> {
        let mut out = BufWriter::new(File::create(&tmp)?);
        out.write_all(MAGIC)?;
        out.write_all(&(json.len() as u64).to_le_bytes())?;
        out.write_all(&json)?;
        let mut zstd = zstd::Encoder::new(out, 3)?;
        zstd.write_all(&body)?;
        zstd.finish()?.into_inner().map_err(|e| e.into_error())?.sync_all()?;
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

pub fn read(path: &Path) -> Result<(Header, TargetSet)> {
    let mut input = BufReader::new(File::open(path).map_err(Error::io(path))?);
    let mut magic = [0u8; 8];
    let mut len = [0u8; 8];
    input.read_exact(&mut magic).map_err(Error::io(path))?;
    if &magic != MAGIC {
        return invalid!("{}: not a pgsum target index", path.display());
    }
    input.read_exact(&mut len).map_err(Error::io(path))?;
    let len = u64::from_le_bytes(len);
    if len > 256 << 20 {
        return invalid!("{}: target index header is too large", path.display());
    }
    let mut json = vec![0; len as usize];
    input.read_exact(&mut json).map_err(Error::io(path))?;
    let header: Header =
        serde_json::from_slice(&json).map_err(|e| Error::Invalid(format!("{}: {e}", path.display())))?;
    if header.schema != SCHEMA {
        return invalid!("{}: schema {} is not {SCHEMA}", path.display(), header.schema);
    }
    let mut body = Vec::new();
    zstd::Decoder::new(input)
        .and_then(|mut d| d.read_to_end(&mut body))
        .map_err(Error::io(path))?;
    if records_sha256(&body) != header.body_sha256 {
        return invalid!("{}: target index differs from its digest", path.display());
    }
    let bad = || Error::Invalid(format!("{}: invalid target index body", path.display()));
    let (keys, mut at) =
        decode(&body, header.targets as usize).map_err(|e| Error::Invalid(format!("{}: {e}", path.display())))?;
    let mut sequences = Vec::with_capacity(header.sequence_targets as usize);
    for _ in 0..header.sequence_targets {
        let contig = *body.get(at).ok_or_else(bad)?;
        let pos = u32::from_le_bytes(body.get(at + 1..at + 5).ok_or_else(bad)?.try_into().expect("4")) as u64;
        at += 5;
        let mut alleles = [String::new(), String::new()];
        for a in &mut alleles {
            let len = u16::from_le_bytes(body.get(at..at + 2).ok_or_else(bad)?.try_into().expect("2")) as usize;
            *a = String::from_utf8(body.get(at + 2..at + 2 + len).ok_or_else(bad)?.to_vec()).map_err(|_| bad())?;
            at += 2 + len;
        }
        let [ref_allele, alt] = alleles;
        sequences.push((contig, Variant { pos, ref_allele, alt }));
    }
    if at != body.len() {
        return Err(bad());
    }
    let set = TargetSet::new(keys, sequences, header.packs.clone())?;
    Ok((header, set))
}

/// Targets for `packs`: from `cache` when it matches the packs and reference exactly, otherwise collected
/// from the packs (and written to `cache`, if given).
pub fn targets(
    packs: &[PathBuf],
    reference: &ReferenceIdentity,
    cache: Option<&Path>,
    positions: bool,
) -> Result<(TargetSet, Source)> {
    if let Some(cache) = cache
        && cache.exists()
        && let Ok((header, set)) = read(cache)
    {
        let refs = pack_refs(packs, reference)?;
        if &header.reference == reference && header.packs == refs && header.term_positions == positions {
            return Ok((set, Source::Cache));
        }
    }
    let set = collect(packs, reference, positions)?;
    if let Some(cache) = cache {
        write(cache, reference, &set, positions)?;
    }
    Ok((set, Source::Packs))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keys_round_trip() {
        let keys = vec![0u64, 1, 2, 1 << 40, (1 << 40) + 5, u64::MAX >> 1];
        let body = encode(&keys);
        assert_eq!(decode(&body, keys.len()).unwrap(), (keys.clone(), body.len()));
        assert!(decode(&body, keys.len() + 1).is_err());
    }
}
