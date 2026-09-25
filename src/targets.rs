//! The union of target sites across packs, and a cache of it.
//!
//! Collecting targets reads every pack. For the whole Catalog that is ~32 GB, so the union can be saved to a
//! target index (`.pgst`) together with the identity of every pack it came from. A later run checks only the
//! packs' headers (a few KB each) and reuses the index if the pack set and reference are unchanged.
//!
//! Index layout: the magic bytes `PGSUMTI1`, a little-endian `u64` header length, the JSON header, then one
//! zstd frame of the sorted keys as varints of the difference from the previous key.

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::genotypes::PackRef;
use crate::pack::{Pack, ReferenceIdentity, records_sha256};
use crate::{Error, Result, invalid};

pub const MAGIC: &[u8; 8] = b"PGSUMTI1";
pub const SCHEMA: &str = "pgsum-targets-v1";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Header {
    pub schema: String,
    pub pgsum_version: String,
    pub reference: ReferenceIdentity,
    /// The packs the targets came from, sorted by PGS ID.
    pub packs: Vec<PackRef>,
    pub targets: u64,
    /// SHA-256 of the uncompressed key bytes.
    pub body_sha256: String,
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

/// Read every pack (in parallel, target columns only) and return the sorted, distinct target keys and the
/// packs' identities (sorted by PGS ID).
pub fn collect(packs: &[PathBuf], reference: &ReferenceIdentity) -> Result<(Vec<u64>, Vec<PackRef>)> {
    // Keys from all packs accumulate here and are deduplicated whenever they have doubled.
    let union: Mutex<(Vec<u64>, usize)> = Mutex::new((Vec::new(), 0));
    let refs: Vec<PackRef> = packs
        .par_iter()
        .map(|path| -> Result<PackRef> {
            let (header, keys) = Pack::open_target_keys(path)?;
            check_reference(&header.pgs_id, &header.reference, reference)?;
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
    Ok((keys, sorted_refs(refs)))
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

fn decode(body: &[u8], n: usize) -> Result<Vec<u64>> {
    let bad = || Error::Invalid("invalid target index body".into());
    let mut keys = Vec::with_capacity(n);
    let (mut at, mut previous) = (0usize, 0u64);
    while at < body.len() {
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
    if keys.len() != n {
        return Err(bad());
    }
    Ok(keys)
}

pub fn write(path: &Path, reference: &ReferenceIdentity, packs: &[PackRef], keys: &[u64]) -> Result<()> {
    let body = encode(keys);
    let header = Header {
        schema: SCHEMA.into(),
        pgsum_version: env!("CARGO_PKG_VERSION").into(),
        reference: reference.clone(),
        packs: packs.to_vec(),
        targets: keys.len() as u64,
        body_sha256: records_sha256(&body),
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

pub fn read(path: &Path) -> Result<(Header, Vec<u64>)> {
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
    let keys =
        decode(&body, header.targets as usize).map_err(|e| Error::Invalid(format!("{}: {e}", path.display())))?;
    Ok((header, keys))
}

/// Targets for `packs`: from `cache` when it matches the packs and reference exactly, otherwise collected
/// from the packs (and written to `cache`, if given).
pub fn targets(
    packs: &[PathBuf],
    reference: &ReferenceIdentity,
    cache: Option<&Path>,
) -> Result<(Vec<u64>, Vec<PackRef>, Source)> {
    if let Some(cache) = cache
        && cache.exists()
    {
        let (header, keys) = read(cache)?;
        let refs = pack_refs(packs, reference)?;
        if &header.reference == reference && header.packs == refs {
            return Ok((keys, refs, Source::Cache));
        }
    }
    let (keys, refs) = collect(packs, reference)?;
    if let Some(cache) = cache {
        write(cache, reference, &refs, &keys)?;
    }
    Ok((keys, refs, Source::Packs))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keys_round_trip() {
        let keys = vec![0u64, 1, 2, 1 << 40, (1 << 40) + 5, u64::MAX >> 1];
        assert_eq!(decode(&encode(&keys), keys.len()).unwrap(), keys);
        assert!(decode(&encode(&keys), keys.len() + 1).is_err());
    }
}
