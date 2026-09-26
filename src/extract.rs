//! Read a gVCF once and assess every target the packs need.
//!
//! A target is one oriented SNV `(contig, pos, REF, ALT)` from a pack term that has no review reasons and a
//! resolved orientation; other terms never need a genotype (see `targets`). The gVCF is decompressed on several threads and
//! scanned once in file order. Each record overlapping at least one target is kept verbatim, and after the
//! scan every target is assessed from its overlapping records, in parallel.

use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::time::Instant;

use rayon::prelude::*;

use crate::alleles::Variant;
use crate::digest::{HashingReader, cached_file_sha256, remember_file_sha256};
use crate::genotype::{self, Policy, Target};
use crate::genotypes::{CompactCall, GenotypeTable, is_sequence_key, key_position, unpack_key};
use crate::gvcf::{self, HeaderFacts};
use crate::index::{Chunk, Index, merge};
use crate::pack::ReferenceIdentity;
use crate::reference::Reference;
use crate::term::{CONTIGS, contig_code_of_name};
use crate::{Error, Result, invalid};

/// Records overlapping one target, beyond the first (rare).
type Extra = HashMap<u32, Vec<u32>>;

pub struct Extracted {
    pub header: HeaderFacts,
    pub gvcf_sha256: String,
    pub gvcf_bytes: u64,
    pub records_scanned: u64,
    /// Kept records, verbatim and without line endings, concatenated.
    pub arena: Vec<u8>,
    /// Start of each kept record in `arena` (one extra entry for the end).
    pub offsets: Vec<u64>,
    /// Per target: first overlapping kept record, or `u32::MAX`.
    pub first: Vec<u32>,
    pub extra: Extra,
}

/// How `scan` reads the gVCF.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum ScanMode {
    /// Use the gVCF's tabix or CSI index when the targets need less than half of the file, else read it
    /// whole.
    #[default]
    Auto,
    /// Read the whole file.
    Full,
    /// Read only the index chunks overlapping targets (fails without a usable index).
    Indexed,
}

/// How `extract` reads the gVCF and calls genotypes.
#[derive(Clone, Debug, Default)]
pub struct Options<'a> {
    /// Target index to reuse when it matches the packs, or to create (`.pgst`).
    pub targets_cache: Option<&'a Path>,
    /// Read haploid chrX/chrY calls (`1`) as homozygous (`1/1`).
    pub haploid_xy_as_homozygous: bool,
    /// Accept calls that report neither depth nor GQ (genotype-only VCFs).
    pub accept_missing_quality: bool,
    /// The sample to read from a multi-sample VCF.
    pub sample: Option<&'a str>,
    pub scan: ScanMode,
    /// Worker threads (at least one).
    pub threads: usize,
}

/// Share of the compressed file above which reading it whole beats seeking chunk by chunk.
const INDEXED_SHARE_LIMIT: f64 = 0.5;

/// Decompressed bytes per chunk of the full scan.
const CHUNK_BYTES: usize = 4 << 20;

/// Targets closer than this are read as one region (one tabix linear-index window).
const REGION_GAP: u64 = 1 << 14;

/// The per-record matching shared by both ways of reading: keeps each record that overlaps a target.
pub(crate) struct Collector<'a> {
    /// The sample's column, when the VCF has several: records are then kept with that column only.
    sample_column: Option<usize>,
    keys: &'a [u64],
    ends: &'a [u64],
    max_span: u64,
    /// Targets by contig code: the index range into `keys`.
    contig_ranges: [(usize, usize); 26],
    records_scanned: u64,
    canonical_records: u64,
    arena: Vec<u8>,
    offsets: Vec<u64>,
    first: Vec<u32>,
    extra: Extra,
    seen: [bool; 26],
    current: (u8, u64), // contig code, last POS
    cursor: usize,      // first target of the current contig not yet passed
    last_chrom: Vec<u8>,
    last_code: Option<u8>,
    /// Sample columns in the file, and whether its header defines DeepVariant's `RefCall` filter.
    samples: usize,
    refcall_defined: bool,
}

impl<'a> Collector<'a> {
    pub(crate) fn header_samples(&self) -> usize {
        self.samples
    }

    pub(crate) fn refcall_defined(&self) -> bool {
        self.refcall_defined
    }

    pub(crate) fn contig_ranges(&self) -> [(usize, usize); 26] {
        self.contig_ranges
    }

    pub(crate) fn max_span(&self) -> u64 {
        self.max_span
    }

    pub(crate) fn records_scanned(&self) -> u64 {
        self.records_scanned
    }

    pub(crate) fn new(keys: &'a [u64], ends: &'a [u64], header: &HeaderFacts) -> Self {
        let max_span = keys
            .iter()
            .zip(ends)
            .map(|(&k, &e)| e - key_position(k).1 as u64)
            .max()
            .unwrap_or(0);
        Collector {
            sample_column: (header.samples > 1 && !header.all_samples).then_some(header.sample_column),
            keys,
            ends,
            max_span,
            contig_ranges: contig_ranges(keys),
            records_scanned: 0,
            canonical_records: 0,
            arena: Vec::new(),
            offsets: vec![0],
            first: vec![u32::MAX; keys.len()],
            extra: Extra::new(),
            seen: [false; 26],
            current: (0, 0),
            cursor: 0,
            last_chrom: Vec::new(),
            last_code: None,
            samples: header.samples,
            refcall_defined: header.refcall_defined,
        }
    }

    /// One record line, without its line ending.
    fn record(&mut self, line: &[u8], text: &str) -> Result<()> {
        self.records_scanned += 1;
        let span = gvcf::interval(text)?;
        // Records arrive grouped by contig, so the name lookup only runs when the contig changes.
        if span.chrom.as_bytes() != self.last_chrom.as_slice() {
            self.last_chrom.clear();
            self.last_chrom.extend_from_slice(span.chrom.as_bytes());
            self.last_code = contig_code_of_name(span.chrom);
        }
        let Some(code) = self.last_code else {
            return Ok(());
        };
        self.canonical_records += 1;
        if code != self.current.0 {
            if self.seen[code as usize] {
                return invalid!("{} records are not contiguous", span.chrom);
            }
            self.seen[code as usize] = true;
            self.current = (code, 0);
            self.cursor = self.contig_ranges[code as usize].0;
        }
        if span.pos < self.current.1 {
            return invalid!("records are not sorted by position");
        }
        self.current.1 = span.pos;
        let (keys, ends) = (self.keys, self.ends);
        let end_of_contig = self.contig_ranges[code as usize].1;
        // Targets ending before this record's POS can't overlap it or any later record; with spans up to
        // `max_span`, any target starting more than `max_span` before it has ended.
        while self.cursor < end_of_contig && key_position(keys[self.cursor]).1 as u64 + self.max_span < span.pos {
            self.cursor += 1;
        }
        let mut t = self.cursor;
        let mut kept: Option<u32> = None;
        while t < end_of_contig && (key_position(keys[t]).1 as u64) <= span.end {
            if ends[t] >= span.pos {
                let id = *kept.get_or_insert_with(|| {
                    keep_line(self.sample_column, line, &mut self.arena);
                    self.offsets.push(self.arena.len() as u64);
                    (self.offsets.len() - 2) as u32
                });
                if self.first[t] == u32::MAX {
                    self.first[t] = id;
                } else {
                    self.extra.entry(t as u32).or_default().push(id);
                }
            }
            t += 1;
        }
        Ok(())
    }

    /// Append a chunk scanned by `scan_chunk`, checking order across chunks as `record` does. `lines_before`
    /// is the file's line count before the chunk; errors carry file line numbers.
    fn merge(&mut self, scan: ChunkScan, lines_before: u64) -> std::result::Result<(), (u64, Error)> {
        self.check_order(&scan, lines_before)?;
        let base_id = (self.offsets.len() - 1) as u32;
        let base = self.arena.len() as u64;
        self.arena.extend_from_slice(&scan.arena);
        self.offsets.extend(scan.offsets[1..].iter().map(|o| o + base));
        for (id, t) in scan.hits {
            let id = base_id + id;
            if self.first[t as usize] == u32::MAX {
                self.first[t as usize] = id;
            } else {
                self.extra.entry(t).or_default().push(id);
            }
        }
        Ok(())
    }

    /// Check a chunk's records continue the file in order (sorted, each contig contiguous) and count them.
    pub(crate) fn check_order(&mut self, scan: &ChunkScan, lines_before: u64) -> std::result::Result<(), (u64, Error)> {
        for &(code, first, last, line) in &scan.segments {
            if code != self.current.0 {
                if self.seen[code as usize] {
                    let name = CONTIGS[code as usize - 1];
                    return Err((
                        lines_before + line,
                        Error::Invalid(format!("{name} records are not contiguous")),
                    ));
                }
                self.seen[code as usize] = true;
                self.current = (code, 0);
            }
            if first < self.current.1 {
                return Err((
                    lines_before + line,
                    Error::Invalid("records are not sorted by position".into()),
                ));
            }
            self.current.1 = last;
        }
        self.records_scanned += scan.records;
        self.canonical_records += scan.canonical_records;
        Ok(())
    }

    fn finish(self, header: HeaderFacts, gvcf_sha256: String, gvcf_bytes: u64) -> Extracted {
        Extracted {
            header,
            gvcf_sha256,
            gvcf_bytes,
            records_scanned: self.records_scanned,
            arena: self.arena,
            offsets: self.offsets,
            first: self.first,
            extra: self.extra,
        }
    }
}

/// Append a kept record: the whole line, or for a multi-sample VCF its first nine columns and the chosen
/// sample's.
fn keep_line(sample_column: Option<usize>, line: &[u8], arena: &mut Vec<u8>) {
    let Some(column) = sample_column else {
        arena.extend_from_slice(line);
        return;
    };
    for (i, field) in line.split(|&b| b == b'\t').enumerate() {
        if i < 9 || i == column {
            if i > 0 {
                arena.push(b'\t');
            }
            arena.extend_from_slice(field);
        }
        if i >= column && i >= 9 {
            break;
        }
    }
}

/// Records of one chunk of the file, found overlapping targets by `scan_chunk`.
#[derive(Default)]
pub(crate) struct ChunkScan {
    /// Kept records (as `keep_line` writes them) and their end offsets in `arena`, from 0.
    pub(crate) arena: Vec<u8>,
    pub(crate) offsets: Vec<u64>,
    /// (kept record in this chunk, target) in file order.
    pub(crate) hits: Vec<(u32, u32)>,
    pub(crate) lines: u64,
    records: u64,
    canonical_records: u64,
    /// Runs of records on one contig: code, first and last POS, and the line (in the chunk) of the first.
    pub(crate) segments: Vec<(u8, u64, u64, u64)>,
}

/// Find the records of a chunk of whole lines that overlap targets, as `Collector::record` would. Errors
/// carry the line number within the chunk.
fn scan_chunk(chunk: &[u8], c: &Collector) -> std::result::Result<ChunkScan, (u64, Error)> {
    let mut out = ChunkScan {
        offsets: vec![0],
        ..ChunkScan::default()
    };
    let (keys, ends) = (c.keys, c.ends);
    let mut last_chrom: &[u8] = b"";
    let mut last_code: Option<u8> = None;
    let mut cursor = 0usize;
    let mut end_of_contig = 0usize;
    let mut start = 0usize;
    while start < chunk.len() {
        let stop = memchr::memchr(b'\n', &chunk[start..]).map_or(chunk.len(), |i| start + i);
        let mut line = &chunk[start..stop];
        start = stop + 1;
        while let [rest @ .., b'\r'] = line {
            line = rest;
        }
        out.lines += 1;
        let at = out.lines;
        let text = std::str::from_utf8(line).map_err(|_| (at, Error::Invalid("invalid UTF-8".into())))?;
        out.records += 1;
        let span = gvcf::interval(text).map_err(|e| (at, e))?;
        if span.chrom.as_bytes() != last_chrom {
            last_chrom = &line[..span.chrom.len()];
            last_code = contig_code_of_name(span.chrom);
        }
        let Some(code) = last_code else { continue };
        out.canonical_records += 1;
        match out.segments.last_mut() {
            Some(seg) if seg.0 == code => {
                if span.pos < seg.2 {
                    return Err((at, Error::Invalid("records are not sorted by position".into())));
                }
                seg.2 = span.pos;
            }
            _ => {
                out.segments.push((code, span.pos, span.pos, at));
                let (lo, hi) = c.contig_ranges[code as usize];
                cursor = lo + keys[lo..hi].partition_point(|&k| key_position(k).1 as u64 + c.max_span < span.pos);
                end_of_contig = hi;
            }
        }
        while cursor < end_of_contig && key_position(keys[cursor]).1 as u64 + c.max_span < span.pos {
            cursor += 1;
        }
        let mut t = cursor;
        let mut kept: Option<u32> = None;
        while t < end_of_contig && (key_position(keys[t]).1 as u64) <= span.end {
            if ends[t] >= span.pos {
                let id = *kept.get_or_insert_with(|| {
                    keep_line(c.sample_column, line, &mut out.arena);
                    out.offsets.push(out.arena.len() as u64);
                    (out.offsets.len() - 2) as u32
                });
                out.hits.push((id, t as u32));
            }
            t += 1;
        }
    }
    Ok(out)
}

fn contig_ranges(keys: &[u64]) -> [(usize, usize); 26] {
    let mut ranges = [(0usize, 0usize); 26];
    for code in 1..=25u8 {
        let lo = keys.partition_point(|&k| key_position(k).0 < code);
        let hi = keys.partition_point(|&k| key_position(k).0 <= code);
        ranges[code as usize] = (lo, hi);
    }
    ranges
}

fn trim_line_end(line: &mut Vec<u8>) {
    while matches!(line.last(), Some(b'\n' | b'\r')) {
        line.pop();
    }
}

/// Scan the gVCF and collect the records overlapping each target. `ends[i]` is the last base of target `i`
/// (its position for SNVs; indel targets span their REF allele).
pub fn scan(gvcf: &Path, keys: &[u64], ends: &[u64], options: &Options) -> Result<Extracted> {
    let mode = options.scan;
    if mode != ScanMode::Full {
        match plan_indexed(gvcf, keys, ends) {
            Ok(Some((index, chunks, share))) if mode == ScanMode::Indexed || share < INDEXED_SHARE_LIMIT => {
                eprintln!(
                    "  reading about {:.1}% of the gVCF through {}",
                    100.0 * share.min(1.0),
                    index.path.display()
                );
                return scan_indexed(gvcf, keys, ends, options.sample, &index, &chunks);
            }
            Ok(Some((_, _, share))) => eprintln!(
                "  reading the whole gVCF: the targets need about {:.1}% of it",
                100.0 * share.min(1.0)
            ),
            Ok(None) if mode == ScanMode::Auto && compression(gvcf)? == Compression::Bgzf => eprintln!(
                "  no .tbi or .csi index for {}; `tabix -p vcf` on it makes runs with few scores much faster",
                gvcf.display()
            ),
            Ok(_) if mode == ScanMode::Indexed => {
                return invalid!("{}: no .tbi or .csi index for --scan indexed", gvcf.display());
            }
            Ok(_) => {}
            Err(e) if mode == ScanMode::Indexed => return Err(e),
            Err(e) => eprintln!("pgsum: not using the gVCF index ({e}); reading the whole file"),
        }
    }
    scan_full(gvcf, keys, ends, options.sample, options.threads)
}

/// The index, the merged chunks the targets need, and their share of the compressed file; `None` without an
/// index.
fn plan_indexed(gvcf: &Path, keys: &[u64], ends: &[u64]) -> Result<Option<(Index, Vec<Chunk>, f64)>> {
    let Some(path) = Index::find(gvcf) else {
        return Ok(None);
    };
    let modified = |p: &Path| std::fs::metadata(p).and_then(|m| m.modified()).map_err(Error::io(p));
    if modified(&path)? < modified(gvcf)? {
        return invalid!("{} is older than the gVCF", path.display());
    }
    let index = Index::open(&path)?;
    let mut ids = [None; 26];
    for (id, name) in index.names().iter().enumerate() {
        if let Some(code) = contig_code_of_name(name) {
            if ids[code as usize].is_some() {
                return invalid!("{}: two sequences named like {name}", path.display());
            }
            ids[code as usize] = Some(id);
        }
    }
    let ranges = contig_ranges(keys);
    let mut chunks = Vec::new();
    for code in 1..=25usize {
        let (lo, hi) = ranges[code];
        let Some(id) = ids[code] else { continue };
        // Merge nearby targets into regions, 0-based half-open.
        let mut region: Option<(u64, u64)> = None;
        for t in lo..hi {
            let (beg, end) = (key_position(keys[t]).1 as u64 - 1, ends[t]);
            region = match region {
                Some((b, e)) if beg <= e + REGION_GAP => Some((b, e.max(end))),
                Some((b, e)) => {
                    chunks.extend(index.chunks(id, b, e));
                    Some((beg, end))
                }
                None => Some((beg, end)),
            };
        }
        if let Some((b, e)) = region {
            chunks.extend(index.chunks(id, b, e));
        }
    }
    let chunks = merge(chunks);
    let file_bytes = std::fs::metadata(gvcf).map_err(Error::io(gvcf))?.len().max(1);
    // Compressed bytes the chunks span, plus about one BGZF block (typically 10–20 kB compressed) per chunk.
    let read: u64 = chunks
        .iter()
        .map(|c| (c.end >> 16).saturating_sub(c.start >> 16) + (1 << 14))
        .sum();
    Ok(Some((index, chunks, read as f64 / file_bytes as f64)))
}

/// Read the header, then only the index chunks the targets need.
fn scan_indexed(
    gvcf: &Path,
    keys: &[u64],
    ends: &[u64],
    sample: Option<&str>,
    index: &Index,
    chunks: &[Chunk],
) -> Result<Extracted> {
    let file = File::open(gvcf).map_err(Error::io(gvcf))?;
    let mut reader = noodles_bgzf::io::Reader::new(std::io::BufReader::new(file));
    let mut header = HeaderFacts {
        requested_sample: sample.map(str::to_owned),
        ..HeaderFacts::default()
    };
    let mut line = Vec::with_capacity(1 << 16);
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line).map_err(Error::io(gvcf))? == 0 {
            return invalid!("{}: no #CHROM header line", gvcf.display());
        }
        trim_line_end(&mut line);
        let text =
            std::str::from_utf8(&line).map_err(|_| Error::Invalid(format!("{}: invalid UTF-8", gvcf.display())))?;
        if !text.starts_with('#') {
            return invalid!("{}: record before the #CHROM header line", gvcf.display());
        }
        if header.read_line(text)? {
            break;
        }
    }
    if !index.names().iter().any(|n| contig_code_of_name(n).is_some()) {
        return invalid!("{}: no records on chr1–chr22, chrX, chrY or chrM", index.path.display());
    }
    let mut collector = Collector::new(keys, ends, &header);
    for chunk in chunks {
        reader
            .seek(noodles_bgzf::VirtualPosition::from(chunk.start))
            .map_err(Error::io(gvcf))?;
        while u64::from(reader.virtual_position()) < chunk.end {
            line.clear();
            if reader.read_until(b'\n', &mut line).map_err(Error::io(gvcf))? == 0 {
                break;
            }
            trim_line_end(&mut line);
            let text =
                std::str::from_utf8(&line).map_err(|_| Error::Invalid(format!("{}: invalid UTF-8", gvcf.display())))?;
            if text.starts_with('#') || text.is_empty() {
                continue;
            }
            collector.record(&line, text).map_err(|e| {
                let at: String = text.split('\t').take(2).collect::<Vec<_>>().join(":");
                Error::Invalid(format!("{}: record at {at}: {e}", gvcf.display()))
            })?;
        }
    }
    let bytes = std::fs::metadata(gvcf).map_err(Error::io(gvcf))?.len();
    Ok(collector.finish(header, cached_file_sha256(gvcf)?, bytes))
}

/// A reader shared with the code that finishes hashing the file.
struct SharedReader<R>(std::sync::Arc<std::sync::Mutex<R>>);

impl<R: std::io::Read> std::io::Read for SharedReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.0
            .lock()
            .map_err(|_| std::io::Error::other("reader poisoned"))?
            .read(buf)
    }
}

/// How the file is compressed, from its first bytes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Compression {
    /// Gzip members with the `BC` extra field: decompressed on several threads.
    Bgzf,
    /// Other gzip: one thread.
    Gzip,
    None,
}

fn compression(path: &Path) -> Result<Compression> {
    let mut start = [0u8; 16];
    let mut file = File::open(path).map_err(Error::io(path))?;
    let n = std::io::Read::read(&mut file, &mut start).map_err(Error::io(path))?;
    let start = &start[..n];
    Ok(if !start.starts_with(&[0x1f, 0x8b]) {
        Compression::None
    } else if n >= 16 && start[3] & 4 != 0 && start[12..14] == *b"BC" {
        Compression::Bgzf
    } else {
        Compression::Gzip
    })
}

/// Read the whole file once: bgzipped VCFs are decompressed on `threads` threads, plain gzip and
/// uncompressed VCFs are read on one.
fn scan_full(gvcf: &Path, keys: &[u64], ends: &[u64], sample: Option<&str>, threads: usize) -> Result<Extracted> {
    let header = HeaderFacts {
        requested_sample: sample.map(str::to_owned),
        ..HeaderFacts::default()
    };
    let (header, collector, digest, bytes) =
        stream_records(gvcf, keys, ends, header, threads, |c, scan, lines_before| {
            c.merge(scan, lines_before)
        })?;
    Ok(collector.finish(header, digest, bytes))
}

/// Read a whole VCF once, handing each chunk of records, scanned for the targets they overlap, to `on_scan`
/// in file order together with the collector (which holds the targets and checks record order) and the number
/// of lines before the chunk. Returns the header facts, the collector, and the file's SHA-256 and size.
pub(crate) fn stream_records<'a>(
    gvcf: &Path,
    keys: &'a [u64],
    ends: &'a [u64],
    mut header: HeaderFacts,
    threads: usize,
    mut on_scan: impl FnMut(&mut Collector<'a>, ChunkScan, u64) -> std::result::Result<(), (u64, Error)>,
) -> Result<(HeaderFacts, Collector<'a>, String, u64)> {
    let file = File::open(gvcf).map_err(Error::io(gvcf))?;
    let hashing = HashingReader::new(file);
    let workers = NonZero::new(threads.max(1)).expect("at least one");
    let kind = compression(gvcf)?;
    let (mut reader, hashed): (Box<dyn BufRead>, _) = {
        // The hashing reader is shared so the digest of the whole file can be finished after reading.
        let shared = std::sync::Arc::new(std::sync::Mutex::new(hashing));
        let source = SharedReader(shared.clone());
        let reader: Box<dyn BufRead> = match kind {
            Compression::Bgzf => Box::new(std::io::BufReader::with_capacity(
                1 << 20,
                noodles_bgzf::io::MultithreadedReader::with_worker_count(workers, source),
            )),
            Compression::Gzip => Box::new(std::io::BufReader::with_capacity(
                1 << 20,
                flate2::read::MultiGzDecoder::new(source),
            )),
            Compression::None => Box::new(std::io::BufReader::with_capacity(1 << 20, source)),
        };
        (reader, shared)
    };
    let at = |line_no: u64, e: Error| Error::Invalid(format!("{}: line {line_no}: {e}", gvcf.display()));
    let mut line = Vec::with_capacity(1 << 16);
    let mut line_no = 0u64;
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line).map_err(Error::io(gvcf))? == 0 {
            return invalid!("{}: no #CHROM header line", gvcf.display());
        }
        line_no += 1;
        trim_line_end(&mut line);
        if line_no == 1 && line.starts_with(b"BCF") {
            return invalid!(
                "{}: BCF is not supported; convert it with `bcftools view -Oz -o out.vcf.gz` and index it",
                gvcf.display()
            );
        }
        let text = std::str::from_utf8(&line).map_err(|_| at(line_no, Error::Invalid("invalid UTF-8".into())))?;
        if !text.starts_with('#') {
            return Err(at(
                line_no,
                Error::Invalid("record before the #CHROM header line".into()),
            ));
        }
        if header.read_line(text).map_err(|e| at(line_no, e))? {
            break;
        }
    }
    // Records: chunks of whole lines, parsed in parallel a batch at a time and handed on in file order, so the
    // result is the same as reading line by line.
    let mut collector = Collector::new(keys, ends, &header);
    let batch_len = 2 * threads.max(1);
    let mut carry: Vec<u8> = Vec::new();
    let mut eof = false;
    while !eof {
        let mut batch: Vec<Vec<u8>> = Vec::with_capacity(batch_len);
        while batch.len() < batch_len && !eof {
            let mut chunk = std::mem::take(&mut carry);
            while chunk.len() < CHUNK_BYTES {
                let buf = reader.fill_buf().map_err(Error::io(gvcf))?;
                if buf.is_empty() {
                    eof = true;
                    break;
                }
                chunk.extend_from_slice(buf);
                let n = buf.len();
                reader.consume(n);
            }
            if !eof {
                match memchr::memrchr(b'\n', &chunk) {
                    Some(i) => carry = chunk.split_off(i + 1),
                    None => {
                        // A line longer than a chunk: keep reading.
                        carry = chunk;
                        continue;
                    }
                }
            }
            if !chunk.is_empty() {
                batch.push(chunk);
            }
        }
        let scans: Vec<_> = batch.par_iter().map(|chunk| scan_chunk(chunk, &collector)).collect();
        for scan in scans {
            let scan = scan.map_err(|(line, e)| at(line_no + line, e))?;
            let lines = scan.lines;
            on_scan(&mut collector, scan, line_no).map_err(|(line, e)| at(line, e))?;
            line_no += lines;
        }
    }
    if collector.canonical_records == 0 {
        return invalid!(
            "{}: no records on chr1–chr22, chrX, chrY or chrM (named chr1 or 1, …)",
            gvcf.display()
        );
    }
    // Finish hashing the file (a decompressor may stop before the final empty BGZF block).
    drop(reader);
    let mut hashing = std::sync::Arc::try_unwrap(hashed)
        .map_err(|_| Error::Invalid("reader still in use".into()))?
        .into_inner()
        .map_err(|_| Error::Invalid("reader poisoned".into()))?;
    std::io::copy(&mut hashing, &mut std::io::sink()).map_err(Error::io(gvcf))?;
    let gvcf_bytes = hashing.bytes;
    let digest = hashing.finish();
    remember_file_sha256(gvcf, &digest);
    Ok((header, collector, digest, gvcf_bytes))
}

/// The target of a key: a sequence target's normalized alleles, or an SNV's REF and ALT (`r` and `a` hold
/// their text).
pub(crate) fn target_for<'t>(
    key: u64,
    sequences: &'t [(u64, Variant)],
    r: &'t mut String,
    a: &'t mut String,
) -> Result<Target<'t>> {
    let (code, pos) = key_position(key);
    let contig = CONTIGS[code as usize - 1];
    let sex_chromosome = matches!(contig, "chrX" | "chrY");
    if is_sequence_key(key) {
        let i = sequences
            .binary_search_by_key(&key, |(k, _)| *k)
            .map_err(|_| Error::Invalid(format!("{contig}:{pos}: sequence target without alleles")))?;
        let v = &sequences[i].1;
        return Ok(Target {
            pos: v.pos,
            ref_allele: &v.ref_allele,
            alt: Some(&v.alt),
            sequence: true,
            sex_chromosome,
        });
    }
    let (_, _, ref_base, alt_base) = unpack_key(key);
    *r = (ref_base as char).to_string();
    *a = (alt_base as char).to_string();
    Ok(Target {
        pos: pos as u64,
        ref_allele: r.as_str(),
        alt: (alt_base != crate::genotypes::ANY_ALT).then_some(a.as_str()),
        sequence: false,
        sex_chromosome,
    })
}

/// Assess every target from its records, on `threads` threads.
pub fn assess(
    keys: &[u64],
    sequences: &[(u64, Variant)],
    scanned: &Extracted,
    reference: &Reference,
    policy: &Policy,
    threads: usize,
) -> Result<Vec<CompactCall>> {
    let record = |id: u32| -> &str {
        let (a, b) = (
            scanned.offsets[id as usize] as usize,
            scanned.offsets[id as usize + 1] as usize,
        );
        std::str::from_utf8(&scanned.arena[a..b]).expect("validated when scanned")
    };
    let chunk = keys.len().div_ceil(threads.max(1)).max(1);
    let results: Vec<Result<Vec<CompactCall>>> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..keys.len())
            .step_by(chunk)
            .map(|start| {
                let end = (start + chunk).min(keys.len());
                let record = &record;
                scope.spawn(move || {
                    let mut out = Vec::with_capacity(end - start);
                    // Neighbouring targets are often covered by the same records (one reference block covers
                    // many), so the last target's parsed records are reused when its record list repeats.
                    let mut ids = Vec::new();
                    let mut last_ids: Vec<u32> = Vec::new();
                    let mut records = Vec::new();
                    for (i, &key) in keys.iter().enumerate().take(end).skip(start) {
                        let (code, pos) = key_position(key);
                        let contig = CONTIGS[code as usize - 1];
                        ids.clear();
                        if scanned.first[i] != u32::MAX {
                            ids.push(scanned.first[i]);
                            ids.extend(scanned.extra.get(&(i as u32)).into_iter().flatten().copied());
                        }
                        if ids != last_ids {
                            records = ids
                                .iter()
                                .map(|&id| {
                                    gvcf::parse_record(record(id))
                                        .map_err(|e| Error::Invalid(format!("{contig}:{pos}: {e}")))
                                })
                                .collect::<Result<Vec<_>>>()?;
                            last_ids.clone_from(&ids);
                        }
                        let (mut r, mut a) = (String::new(), String::new());
                        let target = target_for(key, sequences, &mut r, &mut a)?;
                        let call = genotype::assess(&records, &target, policy, |p, len| {
                            reference.fetch(contig, p - 1, p - 1 + len as u64).ok()
                        });
                        out.push(CompactCall::from(call));
                    }
                    Ok(out)
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("assessment thread panicked"))
            .collect()
    });
    let mut calls = Vec::with_capacity(keys.len());
    for r in results {
        calls.extend(r?);
    }
    Ok(calls)
}

/// Extract genotypes for `packs` from `gvcf` into a table, with the time each phase took.
pub fn extract(
    gvcf: &Path,
    reference: &Reference,
    reference_identity: &ReferenceIdentity,
    packs: &[PathBuf],
    options: &Options,
) -> Result<(GenotypeTable, Timings)> {
    let mut timings = Timings::default();
    let t = Instant::now();
    let (set, source) = crate::targets::targets(packs, reference_identity, options.targets_cache)?;
    timings.targets_s = t.elapsed().as_secs_f64();
    timings.targets_from_cache = source == crate::targets::Source::Cache;
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
    let t = Instant::now();
    let scanned = scan(gvcf, &set.keys, &ends, options)?;
    timings.scan_s = t.elapsed().as_secs_f64();
    check_assembly(gvcf, &scanned.header, reference)?;
    let t = Instant::now();
    let policy = Policy::for_gvcf(
        scanned.header.refcall_defined,
        options.haploid_xy_as_homozygous,
        options.accept_missing_quality,
    );
    let calls = assess(&set.keys, &set.sequences, &scanned, reference, &policy, options.threads)?;
    timings.assess_s = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let table = GenotypeTable::new(gvcf, reference_identity, &policy, set, scanned, calls);
    timings.table_s = t.elapsed().as_secs_f64();
    Ok((table, timings))
}

/// A VCF whose `##contig` lengths differ from the reference's is on another assembly (e.g. GRCh37), and every
/// call would be compared with the wrong bases.
fn check_assembly(gvcf: &Path, header: &HeaderFacts, reference: &Reference) -> Result<()> {
    for &(code, length) in &header.contig_lengths {
        let name = CONTIGS[code as usize - 1];
        if let Some(expected) = reference.contig_length(name)
            && expected != length
        {
            return invalid!(
                "{}: {name} is {length} bp in the VCF header but {expected} bp in {} (a different assembly? pgsum \
                 needs GRCh38)",
                gvcf.display(),
                reference.path.display()
            );
        }
    }
    Ok(())
}

/// Wall-clock seconds per extract phase.
#[derive(Clone, Copy, Debug, Default)]
pub struct Timings {
    pub targets_s: f64,
    /// The targets came from a matching target index rather than from reading the packs.
    pub targets_from_cache: bool,
    pub scan_s: f64,
    pub assess_s: f64,
    pub table_s: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotypes::target_key;

    const RECORDS: &str = "chr1\t1\t.\tA\t<*>\t0\t.\tEND=5\tGT\t0/0
chr1\t6\t.\tC\tT,<*>\t9\tPASS\t.\tGT\t0/1
chrUn_x\t1\t.\tA\t<*>\t0\t.\tEND=9\tGT\t0/0
chr1\t7\t.\tG\t<*>\t0\t.\tEND=20\tGT\t0/0
chr2\t3\t.\tAT\tA,<*>\t9\tPASS\t.\tGT\t1/1
chr2\t4\t.\tT\t<*>\t0\t.\tEND=8\tGT\t0/0
chrX\t2\t.\tG\tA,<*>\t9\tPASS\t.\tGT\t1\r
chrX\t5\t.\tC\t<*>\t0\t.\tEND=30\tGT\t0
";

    fn targets() -> (Vec<u64>, Vec<u64>) {
        let keys = vec![
            target_key(1, 2, b'A', b'C'),
            target_key(1, 6, b'C', b'T'),
            target_key(1, 15, b'A', b'G'),
            target_key(2, 3, b'A', b'C'),
            target_key(2, 8, b'A', b'C'),
            target_key(23, 2, b'G', b'A'),
            target_key(23, 29, b'A', b'C'),
        ];
        let ends = keys.iter().map(|&k| key_position(k).1 as u64).collect();
        (keys, ends)
    }

    fn line_by_line(text: &str, keys: &[u64], ends: &[u64]) -> Result<Collector<'static>> {
        let keys: &'static [u64] = Box::leak(keys.to_vec().into_boxed_slice());
        let ends: &'static [u64] = Box::leak(ends.to_vec().into_boxed_slice());
        let mut c = Collector::new(keys, ends, &HeaderFacts::default());
        for line in text.lines() {
            let line = line.trim_end_matches('\r');
            c.record(line.as_bytes(), line)?;
        }
        Ok(c)
    }

    fn chunked(text: &str, cuts: &[usize], keys: &[u64], ends: &[u64]) -> std::result::Result<Collector<'static>, u64> {
        let keys: &'static [u64] = Box::leak(keys.to_vec().into_boxed_slice());
        let ends: &'static [u64] = Box::leak(ends.to_vec().into_boxed_slice());
        let mut c = Collector::new(keys, ends, &HeaderFacts::default());
        let mut lines_before = 0;
        let mut from = 0;
        for &to in cuts.iter().chain([text.len()].iter()) {
            let scan = scan_chunk(&text.as_bytes()[from..to], &c).map_err(|(l, _)| lines_before + l)?;
            let lines = scan.lines;
            c.merge(scan, lines_before).map_err(|(l, _)| l)?;
            lines_before += lines;
            from = to;
        }
        Ok(c)
    }

    type Summary = (u64, u64, Vec<u8>, Vec<u64>, Vec<u32>, Vec<(u32, Vec<u32>)>);

    fn summary(c: &Collector) -> Summary {
        let mut extra: Vec<_> = c.extra.iter().map(|(k, v)| (*k, v.clone())).collect();
        extra.sort();
        (
            c.records_scanned,
            c.canonical_records,
            c.arena.clone(),
            c.offsets.clone(),
            c.first.clone(),
            extra,
        )
    }

    /// Any split into chunks at line ends gives what reading line by line gives.
    #[test]
    fn chunks_merge_like_lines() {
        let (keys, ends) = targets();
        let expected = summary(&line_by_line(RECORDS, &keys, &ends).unwrap());
        assert!(expected.4.iter().all(|&f| f != u32::MAX), "every target has a record");
        let line_ends: Vec<usize> = RECORDS.match_indices('\n').map(|(i, _)| i + 1).collect();
        for mask in 0..1u32 << (line_ends.len() - 1) {
            let cuts: Vec<usize> = (0..line_ends.len() - 1)
                .filter(|i| mask & 1 << i != 0)
                .map(|i| line_ends[i])
                .collect();
            let c = chunked(RECORDS, &cuts, &keys, &ends).unwrap();
            assert_eq!(summary(&c), expected, "cuts {cuts:?}");
        }
    }

    /// Order errors are found at the same line whether or not a chunk boundary falls before it.
    #[test]
    fn order_errors_across_chunks() {
        let (keys, ends) = targets();
        let unsorted = "chr1\t6\t.\tC\tT\t9\tPASS\t.\tGT\t0/1\nchr1\t2\t.\tA\t<*>\t0\t.\t.\tGT\t0/0\n";
        let split = "chr1\t6\t.\tC\tT\t9\tPASS\t.\tGT\t0/1\nchr2\t2\t.\tA\t<*>\t0\t.\t.\tGT\t0/0\nchr1\t9\t.\tA\t<*>\t0\t.\t.\tGT\t0/0\n";
        for (text, line) in [(unsorted, 2), (split, 3)] {
            assert!(line_by_line(text, &keys, &ends).is_err());
            let ends_at: Vec<usize> = text.match_indices('\n').map(|(i, _)| i + 1).collect();
            for cut in [vec![], vec![ends_at[0]], vec![ends_at[line - 2]]] {
                assert_eq!(chunked(text, &cut, &keys, &ends).err(), Some(line as u64), "{cut:?}");
            }
        }
    }
}
