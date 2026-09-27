//! Allele-frequency tables (`pgsum frequencies`).
//!
//! A table keeps, for every record of one or more population VCFs, its position, REF, ALT, one INFO frequency
//! field exactly as written (for example 1000 Genomes' `EUR_AF`) and whether the record carries a multi-allelic
//! flag (1000 Genomes' `MULTI_ALLELIC`, set on sites split into biallelic records). It is built once from the
//! VCFs, which may be large, unindexed and unsorted, and then answers position lookups directly. The header
//! records each source's size, SHA-256 and MD5, and whether the MD5 matches a published manifest (`name`,
//! `bytes`, `md5` per line).

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::digest::{Md5, hex};
use crate::term::contig_code_of_name;
use crate::{Error, Result, invalid};

pub const SCHEMA: &str = "pgsum-frequencies-v1";
pub const MAGIC: &[u8; 8] = b"PGSUMFQ1";
pub const EXTENSION: &str = "pgsf";
const BLOCK_RECORDS: usize = 1 << 16;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Source {
    pub name: String,
    pub bytes: u64,
    pub sha256: String,
    pub md5: String,
    /// Whether the MD5 equals the manifest's, when a manifest was given and lists the file.
    pub manifest_md5_match: Option<bool>,
    /// Records kept, by CHROM.
    pub records_by_contig: BTreeMap<String, u64>,
    /// Records whose CHROM is not a GRCh38 primary contig.
    pub skipped_records: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Block {
    pub first_key: u64,
    pub last_key: u64,
    pub records: u64,
    /// Byte offset of the zstd frame from the end of the header.
    pub offset: u64,
    pub compressed_bytes: u64,
    /// SHA-256 of the uncompressed block.
    pub sha256: String,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Header {
    pub schema: String,
    pub pgsum_version: String,
    /// The INFO field kept, and the INFO flag read as "multi-allelic".
    pub field: String,
    pub multiallelic_flag: String,
    pub sources: Vec<Source>,
    pub records: u64,
    pub blocks: Vec<Block>,
}

/// One population record at a position.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Record {
    pub ref_allele: String,
    pub alt: String,
    /// The INFO field's value as written; `None` when the record lacks the field.
    pub frequency: Option<String>,
    pub multiallelic: bool,
}

fn key(contig: u8, pos: u32) -> u64 {
    (contig as u64) << 32 | pos as u64
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

fn varint(bytes: &[u8], at: &mut usize) -> Option<u64> {
    let mut v = 0u64;
    for shift in (0..64).step_by(7) {
        let b = *bytes.get(*at)?;
        *at += 1;
        v |= ((b & 0x7f) as u64) << shift;
        if b & 0x80 == 0 {
            return Some(v);
        }
    }
    None
}

fn put_text(out: &mut Vec<u8>, text: &[u8]) {
    put_varint(out, text.len() as u64);
    out.extend_from_slice(text);
}

/// A record's payload (without its key): REF, ALT, flags (bit 0 multi-allelic, bit 1 frequency present) and
/// the frequency text.
fn encode_payload(out: &mut Vec<u8>, r: &[u8], a: &[u8], frequency: Option<&[u8]>, multiallelic: bool) {
    put_text(out, r);
    put_text(out, a);
    out.push(multiallelic as u8 | (frequency.is_some() as u8) << 1);
    if let Some(f) = frequency {
        put_text(out, f);
    }
}

fn decode_block(bytes: &[u8]) -> Option<Vec<(u64, Record)>> {
    let mut out = Vec::new();
    let (mut at, mut previous) = (0usize, 0u64);
    let text = |at: &mut usize| -> Option<String> {
        let len = varint(bytes, at)? as usize;
        let s = std::str::from_utf8(bytes.get(*at..*at + len)?).ok()?.to_owned();
        *at += len;
        Some(s)
    };
    while at < bytes.len() {
        previous += varint(bytes, &mut at)?;
        let ref_allele = text(&mut at)?;
        let alt = text(&mut at)?;
        let flags = *bytes.get(at)?;
        at += 1;
        let frequency = if flags & 2 != 0 { Some(text(&mut at)?) } else { None };
        out.push((
            previous,
            Record {
                ref_allele,
                alt,
                frequency,
                multiallelic: flags & 1 != 0,
            },
        ));
    }
    Some(out)
}

/// Reads the file and hashes what it reads.
struct Hashing {
    file: File,
    sha: Sha256,
    md5: Md5,
    bytes: u64,
}

impl Read for Hashing {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.file.read(buf)?;
        self.sha.update(&buf[..n]);
        self.md5.update(&buf[..n]);
        self.bytes += n as u64;
        Ok(n)
    }
}

/// One source's records: keys with offsets into one payload buffer.
struct Scanned {
    source: Source,
    entries: Vec<(u64, u32, u32)>,
    payload: Vec<u8>,
}

fn info_value<'a>(info: &'a [u8], name: &[u8]) -> (Option<&'a [u8]>, bool) {
    let mut found = None;
    for item in info.split(|&b| b == b';') {
        if item == name {
            return (None, true);
        }
        if item.len() > name.len() && item.starts_with(name) && item[name.len()] == b'=' {
            found = Some(&item[name.len() + 1..]);
        }
    }
    (found, false)
}

fn scan(path: &Path, field: &str, flag: &str, manifest: &BTreeMap<String, String>) -> Result<Scanned> {
    let file = File::open(path).map_err(Error::io(path))?;
    let hashing = Hashing {
        file,
        sha: Sha256::new(),
        md5: Md5::default(),
        bytes: 0,
    };
    let mut reader = BufReader::with_capacity(1 << 20, flate2::read::MultiGzDecoder::new(BufReader::new(hashing)));
    let (mut entries, mut payload) = (Vec::new(), Vec::new());
    let (mut by_contig, mut skipped) = (BTreeMap::<String, u64>::new(), 0u64);
    let mut line = Vec::with_capacity(1 << 16);
    let (field, flag) = (field.as_bytes(), flag.as_bytes());
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line).map_err(Error::io(path))? == 0 {
            break;
        }
        if line.first() == Some(&b'#') {
            continue;
        }
        let mut fields = line.splitn(9, |&b| b == b'\t');
        let mut next = || {
            fields
                .next()
                .ok_or_else(|| Error::Invalid(format!("{}: short record", path.display())))
        };
        let (chrom, pos, _id, r, a, _qual, _filter, info) =
            (next()?, next()?, next()?, next()?, next()?, next()?, next()?, next()?);
        let info = info.strip_suffix(b"\n").unwrap_or(info);
        let chrom = std::str::from_utf8(chrom).map_err(|_| Error::Invalid(format!("{}: CHROM", path.display())))?;
        let Some(contig) = contig_code_of_name(chrom) else {
            skipped += 1;
            continue;
        };
        let pos: u32 = std::str::from_utf8(pos)
            .ok()
            .and_then(|p| p.parse().ok())
            .ok_or_else(|| Error::Invalid(format!("{}: POS", path.display())))?;
        let (frequency, _) = info_value(info, field);
        let (_, multiallelic) = info_value(info, flag);
        let start = payload.len() as u32;
        encode_payload(&mut payload, r, a, frequency, multiallelic);
        entries.push((key(contig, pos), start, payload.len() as u32 - start));
        *by_contig.entry(chrom.to_owned()).or_default() += 1;
    }
    let hashing = reader.into_inner().into_inner().into_inner();
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let md5 = hashing.md5.finish();
    Ok(Scanned {
        source: Source {
            manifest_md5_match: manifest.get(&name).map(|m| *m == md5),
            name,
            bytes: hashing.bytes,
            sha256: format!("sha256:{}", hex(&hashing.sha.finalize())),
            md5,
            records_by_contig: by_contig,
            skipped_records: skipped,
        },
        entries,
        payload,
    })
}

/// A manifest of `name`, `bytes`, `md5` lines (tab-separated).
pub fn read_manifest(path: &Path) -> Result<BTreeMap<String, String>> {
    let text = std::fs::read_to_string(path).map_err(Error::io(path))?;
    Ok(text
        .lines()
        .filter_map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f.len() == 3).then(|| (f[0].to_owned(), f[2].to_owned()))
        })
        .collect())
}

/// Build a table from `vcfs` (gzip or bgzip), keeping INFO `field` and flag `flag`.
pub fn build(vcfs: &[PathBuf], field: &str, flag: &str, manifest: Option<&Path>, out: &Path) -> Result<Header> {
    let manifest = manifest.map(read_manifest).transpose()?.unwrap_or_default();
    let mut scanned: Vec<Scanned> = vcfs
        .par_iter()
        .map(|p| scan(p, field, flag, &manifest))
        .collect::<Result<_>>()?;
    // Every record, in position order; records at one position keep their source and file order.
    let mut order: Vec<(u64, u32, u32)> = Vec::new();
    for (i, s) in scanned.iter().enumerate() {
        order.extend(s.entries.iter().enumerate().map(|(j, e)| (e.0, i as u32, j as u32)));
    }
    order.par_sort_unstable();
    let mut frames: Vec<Vec<u8>> = Vec::new();
    let mut blocks = Vec::new();
    let mut offset = 0u64;
    for chunk in order.chunks(BLOCK_RECORDS) {
        let mut body = Vec::with_capacity(chunk.len() * 16);
        let mut previous = 0u64;
        for &(k, i, j) in chunk {
            put_varint(&mut body, k - previous);
            previous = k;
            let (_, start, len) = scanned[i as usize].entries[j as usize];
            body.extend_from_slice(&scanned[i as usize].payload[start as usize..(start + len) as usize]);
        }
        let frame = zstd::bulk::compress(&body, 6).map_err(Error::io(out))?;
        blocks.push(Block {
            first_key: chunk[0].0,
            last_key: chunk[chunk.len() - 1].0,
            records: chunk.len() as u64,
            offset,
            compressed_bytes: frame.len() as u64,
            sha256: format!("sha256:{}", hex(&Sha256::digest(&body))),
        });
        offset += frame.len() as u64;
        frames.push(frame);
    }
    let header = Header {
        schema: SCHEMA.into(),
        pgsum_version: env!("CARGO_PKG_VERSION").into(),
        field: field.into(),
        multiallelic_flag: flag.into(),
        sources: scanned.iter_mut().map(|s| s.source.clone()).collect(),
        records: order.len() as u64,
        blocks,
    };
    drop(scanned);
    let json = serde_json::to_vec(&header).map_err(|e| Error::Invalid(e.to_string()))?;
    let tmp = out.with_extension(format!("{EXTENSION}.tmp"));
    let written = (|| -> std::io::Result<()> {
        let mut f = std::io::BufWriter::new(File::create(&tmp)?);
        f.write_all(MAGIC)?;
        f.write_all(&(json.len() as u64).to_le_bytes())?;
        f.write_all(&json)?;
        for frame in &frames {
            f.write_all(frame)?;
        }
        f.into_inner().map_err(|e| e.into_error())?.sync_all()?;
        std::fs::rename(&tmp, out)
    })();
    written.map_err(|e| {
        let _ = std::fs::remove_file(&tmp);
        Error::io(out)(e)
    })?;
    Ok(header)
}

/// An opened table.
pub struct FrequencyTable {
    pub header: Header,
    pub path: PathBuf,
    body: memmap2::Mmap,
    body_start: usize,
    decoded: Vec<OnceLock<Vec<(u64, Record)>>>,
}

impl FrequencyTable {
    pub fn open(path: &Path) -> Result<FrequencyTable> {
        let file = File::open(path).map_err(Error::io(path))?;
        // SAFETY: the table is read-only once written; blocks are checked against their SHA-256 when decoded.
        let body = unsafe { memmap2::Mmap::map(&file) }.map_err(Error::io(path))?;
        if body.get(..8) != Some(MAGIC.as_slice()) {
            return invalid!("{} is not a pgsum frequency table", path.display());
        }
        let len = u64::from_le_bytes(body[8..16].try_into().expect("8 bytes")) as usize;
        let header: Header = serde_json::from_slice(
            body.get(16..16 + len)
                .ok_or_else(|| Error::Invalid(format!("{}: truncated", path.display())))?,
        )
        .map_err(|e| Error::Invalid(format!("{}: {e}", path.display())))?;
        if header.schema != SCHEMA {
            return invalid!("{}: schema {} is not {SCHEMA}", path.display(), header.schema);
        }
        let decoded = header.blocks.iter().map(|_| OnceLock::new()).collect();
        Ok(FrequencyTable {
            header,
            path: path.to_owned(),
            body,
            body_start: 16 + len,
            decoded,
        })
    }

    fn block(&self, i: usize) -> Result<&[(u64, Record)]> {
        if let Some(b) = self.decoded[i].get() {
            return Ok(b);
        }
        let info = &self.header.blocks[i];
        let bad = || Error::Invalid(format!("{}: block {i} is damaged", self.path.display()));
        let start = self.body_start + info.offset as usize;
        let frame = self
            .body
            .get(start..start + info.compressed_bytes as usize)
            .ok_or_else(bad)?;
        let bytes = zstd::bulk::decompress(frame, 1 << 28).map_err(|_| bad())?;
        if format!("sha256:{}", hex(&Sha256::digest(&bytes))) != info.sha256 {
            return Err(bad());
        }
        let records = decode_block(&bytes).ok_or_else(bad)?;
        Ok(self.decoded[i].get_or_init(|| records))
    }

    /// Every record at `contig`:`pos`, in source order.
    pub fn at(&self, contig: u8, pos: u32) -> Result<Vec<&Record>> {
        let k = key(contig, pos);
        let first = self.header.blocks.partition_point(|b| b.last_key < k);
        let mut out = Vec::new();
        for i in first..self.header.blocks.len() {
            if self.header.blocks[i].first_key > k {
                break;
            }
            let block = self.block(i)?;
            let from = block.partition_point(|(key, _)| *key < k);
            out.extend(block[from..].iter().take_while(|(key, _)| *key == k).map(|(_, r)| r));
        }
        Ok(out)
    }
}
