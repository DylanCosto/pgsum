//! Cohorts: every sample of a multi-sample VCF at every target, from one pass over the file, and scores for
//! all of them.
//!
//! Each sample's call at a target is assessed with exactly the single-sample rules (`genotype::assess`) on the
//! record lines reduced to that sample, and stored as a 2-bit code: the ALT dosage of a passing call (0, 1 or
//! 2), or `MISSING`. In a genotype-only VCF a site has only a few distinct sample fields (`0|0`, `0|1`, …), so
//! each distinct field is assessed once per target and the result reused for every sample that has it.
//!
//! File layout (`pgsum-cohort-v1`): the magic bytes `PGSUMCO1`, then zstd frames (the sequence targets, then
//! blocks of up to `BLOCK_TARGETS` targets: a `u64` count, the keys as `u64`, and each target's codes packed
//! four samples to a byte, sample 0 in the low bits), then the JSON header, its `u64` length and the magic
//! bytes `PGSUMCOE`. The header lists every frame (offset from byte 8, compressed length, SHA-256 of its
//! bytes). Blocks are written as the scan finishes each contig, so memory holds only the records still
//! needed.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use memmap2::Mmap;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::alleles::Variant;
use crate::extract::{ChunkScan, Collector, stream_records, target_for};
use crate::genotype::{self, Policy, State};
use crate::genotypes::{
    BlockInfo, Frame, PackRef, PolicyInfo, digest_of_frames, index_sequences, key_position, read_sequences, sha256,
    write_sequences,
};
use crate::gvcf::{self, HeaderFacts};
use crate::pack::{Pack, ReferenceIdentity, SourceFile};
use crate::reference::Reference;
use crate::score::{ExactSum, Options, Plan, plan, term_contribution, term_effect};
use crate::term::CONTIGS;
use crate::{Error, Result, invalid};

pub const MAGIC: &[u8; 8] = b"PGSUMCO1";
const MAGIC_END: &[u8; 8] = b"PGSUMCOE";
pub const SCHEMA: &str = "pgsum-cohort-v1";
/// The code of a sample without a passing call at a target.
pub const MISSING: u8 = 3;
pub const BLOCK_TARGETS: usize = 1 << 16;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Header {
    pub schema: String,
    pub pgsum_version: String,
    pub policy: PolicyInfo,
    pub samples: Vec<String>,
    pub gvcf: SourceFile,
    pub reference: ReferenceIdentity,
    pub packs: Vec<PackRef>,
    pub records_scanned: u64,
    pub targets: u64,
    /// Call states over every target and sample.
    pub states: BTreeMap<String, u64>,
    /// SHA-256 of the frame digests, one per line: the sequence frame, then the blocks in key order.
    pub body_sha256: String,
    pub sequence_frame: Frame,
    pub blocks: Vec<BlockInfo>,
}

/// How `extract_cohort` reads the VCF and calls genotypes.
#[derive(Clone, Debug, Default)]
pub struct CohortOptions<'a> {
    pub targets_cache: Option<&'a Path>,
    pub haploid_xy_as_homozygous: bool,
    pub accept_missing_quality: bool,
    /// Leave records whose ALTs are all structural-variant symbols out (see `extract::Options`).
    pub skip_structural_alleles: bool,
    /// Read split multi-allelic records as one (see `genotype::merge_split`).
    pub merge_split_records: bool,
    pub threads: usize,
}

/// Bytes of one target's packed codes.
pub fn row_bytes(samples: usize) -> usize {
    samples.div_ceil(4)
}

/// The code of sample `s` in a packed row.
pub fn code(row: &[u8], s: usize) -> u8 {
    (row[s / 4] >> ((s % 4) * 2)) & 3
}

/// Every sample `MISSING` (unused bits of the last byte are set too; readers never look at them).
fn missing_row(samples: usize) -> Vec<u8> {
    vec![0xff; row_bytes(samples)]
}

/// Streams finished blocks to the file.
struct Writer {
    out: BufWriter<File>,
    offset: u64,
    blocks: Vec<BlockInfo>,
}

impl Writer {
    fn frame(&mut self, bytes: &[u8]) -> Result<Frame> {
        let compressed = zstd::bulk::compress(bytes, 6).map_err(|e| Error::Invalid(e.to_string()))?;
        self.out
            .write_all(&compressed)
            .map_err(|e| Error::Invalid(e.to_string()))?;
        let frame = Frame {
            offset: self.offset,
            compressed_bytes: compressed.len() as u64,
            sha256: sha256(bytes),
        };
        self.offset += compressed.len() as u64;
        Ok(frame)
    }

    fn block(&mut self, keys: &[u64], rows: &[Vec<u8>]) -> Result<()> {
        let mut bytes = Vec::with_capacity(8 + keys.len() * 8 + rows.iter().map(Vec::len).sum::<usize>());
        bytes.extend_from_slice(&(keys.len() as u64).to_le_bytes());
        for k in keys {
            bytes.extend_from_slice(&k.to_le_bytes());
        }
        for r in rows {
            bytes.extend_from_slice(r);
        }
        let frame = self.frame(&bytes)?;
        self.blocks.push(BlockInfo {
            first_key: keys[0],
            last_key: *keys.last().expect("non-empty block"),
            targets: keys.len() as u64,
            frame,
        });
        Ok(())
    }
}

/// Targets of one contig waiting for their rows, in index order.
struct ContigRows {
    /// Index of the first target not yet written.
    written: usize,
    /// Rows of targets from `written` on, finished in order.
    rows: Vec<Vec<u8>>,
}

/// The scan's state between chunks: record lines still needed and targets not yet finished.
struct Window<'a> {
    keys: &'a [u64],
    sequences: &'a [(u64, Variant)],
    samples: usize,
    policy: Policy,
    reference: &'a Reference,
    /// Kept record lines by global id, with the number of unfinished targets that use each.
    lines: HashMap<u32, (Vec<u8>, u32)>,
    next_line: u32,
    /// Record lines of targets that have been hit and are not finished.
    hits: HashMap<u32, Vec<u32>>,
    /// Per contig code: the next target to finish.
    done: [usize; 26],
    finished: [bool; 26],
    rows: [Option<ContigRows>; 26],
    states: [u64; State::ALL.len()],
    writer: Writer,
    contig_ranges: [(usize, usize); 26],
}

impl<'a> Window<'a> {
    fn add(&mut self, scan: ChunkScan) {
        let base = self.next_line;
        for i in 0..scan.offsets.len() - 1 {
            let line = scan.arena[scan.offsets[i] as usize..scan.offsets[i + 1] as usize].to_vec();
            self.lines.insert(base + i as u32, (line, 0));
        }
        self.next_line += (scan.offsets.len() - 1) as u32;
        for (id, t) in scan.hits {
            let id = base + id;
            self.hits.entry(t).or_default().push(id);
            if let Some(l) = self.lines.get_mut(&id) {
                l.1 += 1;
            }
        }
    }

    /// Finish the targets of contig `code` before index `upto` (all of them when `upto` is its end).
    fn finish(&mut self, code: u8, upto: usize) -> Result<()> {
        let c = code as usize;
        let from = self.done[c].max(self.contig_ranges[c].0);
        if upto <= from {
            return Ok(());
        }
        let contig = CONTIGS[c - 1];
        let reference = self.reference;
        let fetch = |p: u64, len: usize| reference.fetch(contig, p - 1, p - 1 + len as u64).ok();
        let work: Vec<(usize, Vec<u32>)> = (from..upto)
            .map(|t| (t, self.hits.remove(&(t as u32)).unwrap_or_default()))
            .collect();
        let (keys, sequences, samples, policy, lines) =
            (self.keys, self.sequences, self.samples, &self.policy, &self.lines);
        let results: Vec<Result<(Vec<u8>, [u64; State::ALL.len()])>> = work
            .par_iter()
            .map(|(t, ids)| {
                let recs: Vec<&[u8]> = ids.iter().map(|id| lines[id].0.as_slice()).collect();
                target_row(keys[*t], sequences, &recs, samples, policy, &fetch)
                    .map_err(|e| Error::Invalid(format!("{contig}:{}: {e}", key_position(keys[*t]).1)))
            })
            .collect();
        let rows = self.rows[c].get_or_insert_with(|| ContigRows {
            written: from,
            rows: Vec::new(),
        });
        for r in results {
            let (row, states) = r?;
            rows.rows.push(row);
            for (a, b) in self.states.iter_mut().zip(states) {
                *a += b;
            }
        }
        for (_, ids) in &work {
            for id in ids {
                if let Some(l) = self.lines.get_mut(id) {
                    l.1 -= 1;
                    if l.1 == 0 {
                        self.lines.remove(id);
                    }
                }
            }
        }
        self.done[c] = upto;
        self.write_blocks(code, upto == self.contig_ranges[c].1)
    }

    /// Write every full block of the contig's finished rows (and the remainder when the contig is done).
    fn write_blocks(&mut self, code: u8, all: bool) -> Result<()> {
        let c = code as usize;
        let Some(rows) = self.rows[c].as_mut() else {
            return Ok(());
        };
        while rows.rows.len() >= BLOCK_TARGETS || all && !rows.rows.is_empty() {
            let n = rows.rows.len().min(BLOCK_TARGETS);
            let block: Vec<Vec<u8>> = rows.rows.drain(..n).collect();
            let keys = &self.keys[rows.written..rows.written + n];
            self.writer.block(keys, &block)?;
            rows.written += n;
        }
        Ok(())
    }

    /// Finish every target of a contig whose records have all been read.
    fn finish_contig(&mut self, code: u8) -> Result<()> {
        if !self.finished[code as usize] {
            self.finish(code, self.contig_ranges[code as usize].1)?;
            self.finished[code as usize] = true;
        }
        Ok(())
    }
}

/// One target's codes for every sample, and how many sample calls ended in each state.
fn target_row(
    key: u64,
    sequences: &[(u64, Variant)],
    lines: &[&[u8]],
    samples: usize,
    policy: &Policy,
    fetch: &(impl Fn(u64, usize) -> Option<String> + Sync),
) -> Result<(Vec<u8>, [u64; State::ALL.len()])> {
    let mut states = [0u64; State::ALL.len()];
    if lines.is_empty() {
        states[State::UnknownNoRecord as usize] = samples as u64;
        return Ok((missing_row(samples), states));
    }
    let (mut r, mut a) = (String::new(), String::new());
    let target = target_for(key, sequences, &mut r, &mut a)?;
    if let [line] = lines {
        return single_record_row(line, &target, samples, policy, fetch);
    }
    // Each line: its first nine columns and the start of every sample field.
    let split: Vec<(&[u8], Vec<&[u8]>)> = lines
        .iter()
        .map(|line| {
            let mut fields = line.split(|&b| b == b'\t');
            let mut shared_end = 0;
            for _ in 0..9 {
                let f = fields.next().unwrap_or_default();
                shared_end += f.len() + 1;
            }
            let shared = &line[..shared_end.saturating_sub(1).min(line.len())];
            (shared, fields.collect::<Vec<_>>())
        })
        .collect();
    for (_, fields) in &split {
        if fields.len() != samples {
            return invalid!("record has {} sample fields, expected {samples}", fields.len());
        }
    }
    let mut row = vec![0u8; row_bytes(samples)];
    let mut memo: Vec<(Vec<&[u8]>, u8, u8)> = Vec::new();
    let mut text = String::new();
    for s in 0..samples {
        let fields: Vec<&[u8]> = split.iter().map(|(_, f)| f[s]).collect();
        let (c, state) = match memo.iter().find(|(k, _, _)| *k == fields) {
            Some(&(_, c, st)) => (c, st),
            None => {
                let mut records = Vec::with_capacity(split.len());
                for ((shared, _), field) in split.iter().zip(&fields) {
                    text.clear();
                    text.push_str(std::str::from_utf8(shared).map_err(|_| Error::Invalid("invalid UTF-8".into()))?);
                    text.push('\t');
                    text.push_str(std::str::from_utf8(field).map_err(|_| Error::Invalid("invalid UTF-8".into()))?);
                    records.push(gvcf::parse_record(&text)?);
                }
                let call = genotype::assess(&records, &target, policy, fetch);
                let c = match call.alt_dosage {
                    Some(d) if call.state.is_passing() => d.min(2),
                    _ => MISSING,
                };
                memo.push((fields.clone(), c, call.state as u8));
                (c, call.state as u8)
            }
        };
        row[s / 4] |= c << ((s % 4) * 2);
        states[state as usize] += 1;
    }
    Ok((row, states))
}

/// `target_row` for the common case of one record: the sample fields are walked in place and each distinct
/// field is assessed once.
fn single_record_row(
    line: &[u8],
    target: &genotype::Target<'_>,
    samples: usize,
    policy: &Policy,
    fetch: &(impl Fn(u64, usize) -> Option<String> + Sync),
) -> Result<(Vec<u8>, [u64; State::ALL.len()])> {
    let mut states = [0u64; State::ALL.len()];
    let mut tabs = memchr::memchr_iter(b'\t', line);
    let shared_end = tabs
        .nth(8)
        .ok_or_else(|| Error::Invalid("record has too few fields".into()))?;
    let shared = std::str::from_utf8(&line[..shared_end]).map_err(|_| Error::Invalid("invalid UTF-8".into()))?;
    let mut row = vec![0u8; row_bytes(samples)];
    let mut memo: Vec<(&[u8], u8, u8)> = Vec::new();
    let mut text = String::new();
    let mut count = 0usize;
    for (s, field) in line[shared_end + 1..].split(|&b| b == b'\t').enumerate() {
        if s >= samples {
            return invalid!("record has more than {samples} sample fields");
        }
        count += 1;
        let (c, state) = match memo.iter().find(|(k, _, _)| *k == field) {
            Some(&(_, c, st)) => (c, st),
            None => {
                text.clear();
                text.push_str(shared);
                text.push('\t');
                text.push_str(std::str::from_utf8(field).map_err(|_| Error::Invalid("invalid UTF-8".into()))?);
                let record = gvcf::parse_record(&text)?;
                let call = genotype::assess(std::slice::from_ref(&record), target, policy, fetch);
                let c = match call.alt_dosage {
                    Some(d) if call.state.is_passing() => d.min(2),
                    _ => MISSING,
                };
                memo.push((field, c, call.state as u8));
                (c, call.state as u8)
            }
        };
        row[s / 4] |= c << ((s % 4) * 2);
        states[state as usize] += 1;
    }
    if count != samples {
        return invalid!("record has {count} sample fields, expected {samples}");
    }
    Ok((row, states))
}

/// What `extract_cohort` read and wrote.
pub struct CohortSummary {
    pub header: Header,
}

/// Read every sample of `vcf` at every target of `packs` and write the codes to `out`.
pub fn extract_cohort(
    vcf: &Path,
    reference: &Reference,
    identity: &ReferenceIdentity,
    packs: &[PathBuf],
    options: &CohortOptions,
    out: &Path,
) -> Result<CohortSummary> {
    let (set, _) = crate::targets::targets(packs, identity, options.targets_cache)?;
    let ends: Vec<u64> = set
        .keys
        .iter()
        .map(|&k| {
            let pos = key_position(k).1 as u64;
            match set.sequences.binary_search_by_key(&k, |(key, _)| *key) {
                Ok(i) => pos + set.sequences[i].1.ref_allele.len() as u64 - 1,
                Err(_) => pos,
            }
        })
        .collect();
    let tmp = out.with_extension("tmp");
    let file = File::create(&tmp).map_err(Error::io(&tmp))?;
    let mut writer = Writer {
        out: BufWriter::with_capacity(1 << 22, file),
        offset: 0,
        blocks: Vec::new(),
    };
    writer.out.write_all(MAGIC).map_err(Error::io(&tmp))?;
    let mut sequence_bytes = Vec::new();
    write_sequences(&set.sequences, &mut sequence_bytes).expect("writing to memory");
    let sequence_frame = writer.frame(&sequence_bytes)?;
    let header = HeaderFacts {
        all_samples: true,
        ..HeaderFacts::default()
    };
    // The window is created once the header (and so the sample count) is known.
    let mut window: Option<Window> = None;
    let mut writer = Some(writer);
    let mut policy_refcall = false;
    let mut failure: Option<Error> = None;
    let (header, collector, digest, bytes) = stream_records(
        vcf,
        &set.keys,
        &ends,
        header,
        options.threads,
        options.skip_structural_alleles,
        |collector: &mut Collector, scan: ChunkScan, lines_before: u64| {
            collector.check_order(&scan, lines_before)?;
            let w = window.get_or_insert_with(|| {
                let samples = collector.header_samples();
                policy_refcall = collector.refcall_defined();
                Window {
                    keys: &set.keys,
                    sequences: &set.sequences,
                    samples,
                    policy: Policy::for_gvcf(
                        policy_refcall,
                        options.haploid_xy_as_homozygous,
                        options.accept_missing_quality,
                    )
                    .with_merge_split(options.merge_split_records),
                    reference,
                    lines: HashMap::new(),
                    next_line: 0,
                    hits: HashMap::new(),
                    done: [0; 26],
                    finished: [false; 26],
                    rows: Default::default(),
                    states: [0; State::ALL.len()],
                    writer: writer.take().expect("one window"),
                    contig_ranges: collector.contig_ranges(),
                }
            });
            // Contigs this chunk leaves behind are complete.
            let segments = scan.segments.clone();
            w.add(scan);
            let result = (|| -> Result<()> {
                for (i, &(code, _, last, _)) in segments.iter().enumerate() {
                    if i + 1 < segments.len() {
                        w.finish_contig(code)?;
                        continue;
                    }
                    // In the contig still being read, targets ending before the last POS can't gain records.
                    let (lo, hi) = w.contig_ranges[code as usize];
                    let max_span = collector.max_span();
                    let upto = lo + w.keys[lo..hi].partition_point(|&k| key_position(k).1 as u64 + max_span < last);
                    w.finish(code, upto)?;
                }
                Ok(())
            })();
            result.map_err(|e| {
                let msg = e.to_string();
                failure = Some(e);
                (lines_before, Error::Invalid(msg))
            })
        },
    )
    .map_err(|e| failure.take().unwrap_or(e))?;
    let records_scanned = collector.records_scanned();
    let mut w = match window {
        Some(w) => w,
        None => return invalid!("{}: no records", vcf.display()),
    };
    for code in 1..=25u8 {
        w.finish_contig(code)?;
    }
    let samples = header.sample_names.clone();
    let mut states = BTreeMap::new();
    for (i, &n) in w.states.iter().enumerate() {
        if n > 0 {
            states.insert(State::ALL[i].as_str().to_owned(), n);
        }
    }
    let mut writer = w.writer;
    writer.blocks.sort_by_key(|b| b.first_key);
    let body_sha256 = digest_of_frames(
        std::iter::once(sequence_frame.sha256.as_str()).chain(writer.blocks.iter().map(|b| b.frame.sha256.as_str())),
    );
    let policy = Policy::for_gvcf(
        policy_refcall,
        options.haploid_xy_as_homozygous,
        options.accept_missing_quality,
    )
    .with_merge_split(options.merge_split_records);
    let header = Header {
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
            merge_split_records: false,
        }
        .with_structural_skipped(options.skip_structural_alleles)
        .with_split_merged(options.merge_split_records),
        samples,
        gvcf: SourceFile {
            name: vcf
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default(),
            sha256: digest,
            bytes,
        },
        reference: identity.clone(),
        packs: set.packs.clone(),
        records_scanned,
        targets: set.keys.len() as u64,
        states,
        body_sha256,
        sequence_frame,
        blocks: writer.blocks.clone(),
    };
    let json = serde_json::to_vec(&header).map_err(|e| Error::Invalid(e.to_string()))?;
    let io = |e| Error::Io {
        path: tmp.clone(),
        source: e,
    };
    writer.out.write_all(&json).map_err(io)?;
    writer.out.write_all(&(json.len() as u64).to_le_bytes()).map_err(io)?;
    writer.out.write_all(MAGIC_END).map_err(io)?;
    writer
        .out
        .into_inner()
        .map_err(|e| io(e.into_error()))?
        .sync_all()
        .map_err(io)?;
    std::fs::rename(&tmp, out).map_err(io)?;
    Ok(CohortSummary { header })
}

/// A decompressed block: its keys and packed rows.
type Block = (Vec<u64>, Vec<u8>);

/// A cohort file, its blocks decompressed on first use.
pub struct CohortTable {
    pub header: Header,
    path: PathBuf,
    map: Mmap,
    blocks: Vec<OnceLock<std::result::Result<Block, String>>>,
    sequence_index: HashMap<(u8, Variant), u64>,
    row: usize,
}

impl CohortTable {
    pub fn open(path: &Path) -> Result<CohortTable> {
        let bad = |what: &str| Error::Invalid(format!("{}: invalid cohort file ({what})", path.display()));
        let file = File::open(path).map_err(Error::io(path))?;
        // SAFETY: the file is read-only input; changing it while pgsum runs is unsupported.
        let map = unsafe { Mmap::map(&file) }.map_err(Error::io(path))?;
        let n = map.len();
        if n < 24 || &map[..8] != MAGIC || &map[n - 8..] != MAGIC_END {
            return invalid!("{}: not a pgsum cohort file", path.display());
        }
        let len = u64::from_le_bytes(map[n - 16..n - 8].try_into().expect("8")) as usize;
        let start = (n - 16)
            .checked_sub(len)
            .filter(|&s| s >= 8)
            .ok_or_else(|| bad("header length"))?;
        let header: Header = serde_json::from_slice(&map[start..n - 16])
            .map_err(|e| Error::Invalid(format!("{}: {e}", path.display())))?;
        if header.schema != SCHEMA {
            return invalid!("{}: schema {} is not {SCHEMA}", path.display(), header.schema);
        }
        let digest = digest_of_frames(
            std::iter::once(header.sequence_frame.sha256.as_str())
                .chain(header.blocks.iter().map(|b| b.frame.sha256.as_str())),
        );
        if digest != header.body_sha256 || header.blocks.windows(2).any(|w| w[0].last_key >= w[1].first_key) {
            return Err(bad("block list"));
        }
        let row = row_bytes(header.samples.len());
        let mut table = CohortTable {
            path: path.to_owned(),
            blocks: header.blocks.iter().map(|_| OnceLock::new()).collect(),
            sequence_index: HashMap::new(),
            row,
            map,
            header,
        };
        let bytes = table.frame(&table.header.sequence_frame.clone())?;
        let sequences = read_sequences(&bytes, |k| k).ok_or_else(|| bad("sequence targets"))?;
        table.sequence_index = index_sequences(&sequences);
        Ok(table)
    }

    fn frame(&self, frame: &Frame) -> Result<Vec<u8>> {
        let bad = || Error::Invalid(format!("{}: invalid cohort frame", self.path.display()));
        let start = 8usize.checked_add(frame.offset as usize).ok_or_else(bad)?;
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

    fn block(&self, i: usize) -> Result<&Block> {
        self.blocks[i]
            .get_or_init(|| {
                let info = &self.header.blocks[i];
                let bytes = self.frame(&info.frame).map_err(|e| e.to_string())?;
                let n = u64::from_le_bytes(bytes.get(..8).ok_or("short block")?.try_into().expect("8")) as usize;
                let keys_end = 8 + n * 8;
                if n as u64 != info.targets || bytes.len() != keys_end + n * self.row {
                    return Err(format!("{}: invalid cohort block {i}", self.path.display()));
                }
                let keys = bytes[8..keys_end]
                    .as_chunks::<8>()
                    .0
                    .iter()
                    .map(|&c| u64::from_le_bytes(c))
                    .collect();
                Ok((keys, bytes[keys_end..].to_vec()))
            })
            .as_ref()
            .map_err(|e| Error::Invalid(e.clone()))
    }

    /// Block `b`'s keys and packed rows (one row of `row_bytes(samples)` per key).
    pub fn block_rows(&self, b: usize) -> Result<(&[u64], &[u8])> {
        let (keys, rows) = self.block(b)?;
        Ok((keys, rows))
    }

    /// Decompress every block now, in parallel.
    pub fn preload(&self) -> Result<()> {
        (0..self.blocks.len())
            .into_par_iter()
            .try_for_each(|i| self.block(i).map(|_| ()))
    }

    /// A target's packed codes, if the file has the target.
    pub fn row(&self, key: u64) -> Result<Option<&[u8]>> {
        let i = self.header.blocks.partition_point(|b| b.first_key <= key);
        if i == 0 || key > self.header.blocks[i - 1].last_key {
            return Ok(None);
        }
        let (keys, rows) = self.block(i - 1)?;
        Ok(keys
            .binary_search(&key)
            .ok()
            .map(|j| &rows[j * self.row..(j + 1) * self.row]))
    }

    pub fn sequence_key(&self, contig: u8, variant: &Variant) -> Option<u64> {
        self.sequence_index.get(&(contig, variant.clone())).copied()
    }
}

/// One score for every sample of a cohort.
pub struct CohortScore {
    pub pgs_id: String,
    pub total_terms: u64,
    /// Per sample: the exact partial sum, scorable terms and scorable absolute effect.
    pub sums: Vec<ExactSum>,
    pub scorable: Vec<u64>,
    pub effect_scorable: Vec<f64>,
    pub effect_all: f64,
    /// Per sample: the exponent its sum is written with, as `score::score` would (the smallest exponent among
    /// the contributions it counts, and at most 0).
    pub exponents: Vec<i64>,
}

impl CohortScore {
    /// Sample `s`'s exact partial sum, written as `score::score` writes it.
    pub fn text(&self, s: usize) -> String {
        self.sums[s].finish().rescaled(self.exponents[s]).to_python_string()
    }
}

/// Score a pack for every sample: the same terms, rules and exact arithmetic as `score::score`, per sample.
///
/// A term's contribution at code 0 is added once to a shared baseline and removed for the samples that lack a
/// passing call, so the per-sample work is proportional to the samples that differ from the most common code.
pub fn score_cohort(pack: &Pack, cohort: &CohortTable, options: &Options) -> Result<CohortScore> {
    let h = &pack.header;
    if !cohort
        .header
        .packs
        .iter()
        .any(|p| p.pgs_id == h.pgs_id && p.records_sha256 == h.records_sha256)
    {
        return invalid!("the cohort file was not extracted with this {} pack", h.pgs_id);
    }
    let n = cohort.header.samples.len();
    let mut baseline = ExactSum::default();
    let mut sums = vec![ExactSum::default(); n];
    let mut scorable = vec![0u64; n];
    let mut effect_scorable = vec![0f64; n];
    let (mut base_scorable, mut base_effect, mut effect_all, mut total) = (0u64, 0f64, 0f64, 0u64);
    // Exponents: of the baseline's contributions (with counts), and per sample the smallest of its own
    // contributions and how many baseline contributions of each exponent it left.
    let mut base_exponents: BTreeMap<i64, u64> = BTreeMap::new();
    let mut own_min = vec![0i64; n];
    let mut left: Vec<Vec<(i64, u64)>> = vec![Vec::new(); n];
    for term in pack.terms_from(0) {
        let term = term?;
        total += 1;
        let effect = term_effect(&term);
        if let Some(e) = effect {
            effect_all += e;
        }
        let (plan, _, _) = plan(&term, options);
        let (key, effect_is_alt) = match plan {
            Plan::Unscorable(_) => continue,
            Plan::Snv { key, effect_is_alt } => (Some(key), effect_is_alt),
            Plan::Sequence { variant, effect_is_alt } => (cohort.sequence_key(term.contig, &variant), effect_is_alt),
        };
        let Some(row) = key.map(|k| cohort.row(k)).transpose()?.flatten() else {
            return invalid!(
                "the cohort file has no call for {}:{}; extract with this pack",
                CONTIGS[term.contig as usize - 1],
                term.pos
            );
        };
        // Contribution at each code (ALT dosage); None where the model gives no weight.
        let at = |c: u8| term_contribution(&term, if effect_is_alt { c } else { 2 - c });
        let contributions = [at(0), at(1), at(2)];
        let e = effect.unwrap_or(0.0);
        // Code 0 goes into the baseline when it has a contribution.
        let base = contributions[0].clone();
        if let Some(c) = &base {
            baseline.add(c);
            base_scorable += 1;
            base_effect += e;
            *base_exponents.entry(c.exponent).or_default() += 1;
        }
        for s in 0..n {
            let byte = row[s / 4];
            if byte == 0 {
                // Four samples at code 0: all in the baseline already.
                continue;
            }
            let code = (byte >> ((s % 4) * 2)) & 3;
            if code == 0 {
                continue;
            }
            // Leave the baseline for this sample...
            if let Some(c) = &base {
                sums[s].add(&c.negated());
                scorable[s] = scorable[s].wrapping_sub(1);
                effect_scorable[s] -= e;
                match left[s].iter_mut().find(|(x, _)| *x == c.exponent) {
                    Some((_, k)) => *k += 1,
                    None => left[s].push((c.exponent, 1)),
                }
            }
            // ...and add what its own code contributes.
            if code != MISSING
                && let Some(c) = &contributions[code as usize]
            {
                sums[s].add(c);
                scorable[s] = scorable[s].wrapping_add(1);
                effect_scorable[s] += e;
                own_min[s] = own_min[s].min(c.exponent);
            }
        }
    }
    let mut exponents = Vec::with_capacity(n);
    for s in 0..n {
        sums[s].merge(baseline.clone());
        scorable[s] = scorable[s].wrapping_add(base_scorable);
        effect_scorable[s] += base_effect;
        let kept_base = base_exponents
            .iter()
            .filter(|(e, count)| {
                let gone = left[s].iter().find(|(x, _)| x == *e).map_or(0, |(_, k)| *k);
                **count > gone
            })
            .map(|(e, _)| *e)
            .min()
            .unwrap_or(0);
        exponents.push(own_min[s].min(kept_base).min(0));
    }
    Ok(CohortScore {
        pgs_id: h.pgs_id.clone(),
        total_terms: total,
        sums,
        scorable,
        effect_scorable,
        effect_all,
        exponents,
    })
}
