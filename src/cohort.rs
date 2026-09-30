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
//! needed. Quantitative v2 blocks append a JSON list of [row, [[sample, measurement], ...]] pairs.
//! Their packed slots are MISSING; callers must consult `measurements` as well as `row`.
//! New extraction writes v3, with separately compressed missing-call states and source-record ordinals.
//! Ordinary scoring does not load diagnostic frames. Quantitative blocks cap at 262,144 sample-target
//! cells and GT blocks at 8,388,608 cells (or one row for wider cohorts), with at most 65,536 targets.
//! Readers retain v1/v2 compatibility; legacy missing calls have no recoverable exclusion reason.

use crate::cohort_diagnostics::{MissingStates, PASSING, Row as Diagnostics};
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};

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
pub const QUANTITATIVE_SCHEMA: &str = "pgsum-cohort-v2";
pub const DIAGNOSTIC_SCHEMA: &str = "pgsum-cohort-v3";
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
    /// SHA-256 of frame digests, one per line: sequence, blocks in key order, then v3 diagnostics.
    pub body_sha256: String,
    pub sequence_frame: Frame,
    pub blocks: Vec<BlockInfo>,
    /// v3: one diagnostic frame per genotype block, in the same target order.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub diagnostics: Vec<Frame>,
}

/// How `extract_cohort` reads the VCF and calls genotypes.
#[derive(Clone, Debug, Default)]
pub struct CohortOptions<'a> {
    pub dosage_field: crate::dosage::Field,
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
#[derive(Default)]
struct AssessedRow {
    codes: Vec<u8>,
    measurements: Vec<(usize, crate::dosage::Measurement)>,
    diagnostics: Diagnostics,
}

struct Writer {
    quantitative: bool,
    out: BufWriter<File>,
    offset: u64,
    blocks: Vec<BlockInfo>,
    diagnostics: BTreeMap<u64, Frame>,
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

    fn block(&mut self, keys: &[u64], rows: &[AssessedRow]) -> Result<()> {
        let mut bytes = Vec::with_capacity(8 + keys.len() * 8 + rows.iter().map(|r| r.codes.len()).sum::<usize>());
        bytes.extend_from_slice(&(keys.len() as u64).to_le_bytes());
        for k in keys {
            bytes.extend_from_slice(&k.to_le_bytes());
        }
        for r in rows {
            bytes.extend_from_slice(&r.codes);
        }
        if self.quantitative {
            let extra: Vec<_> = rows
                .iter()
                .enumerate()
                .filter(|(_, r)| !r.measurements.is_empty())
                .map(|(i, r)| (i, &r.measurements))
                .collect();
            let json = serde_json::to_vec(&extra).map_err(|e| Error::Invalid(e.to_string()))?;
            bytes.extend_from_slice(&json);
        }
        let frame = self.frame(&bytes)?;
        let diagnostics = serde_json::to_vec(&rows.iter().map(|r| &r.diagnostics).collect::<Vec<_>>())
            .map_err(|e| Error::Invalid(e.to_string()))?;
        let diagnostic_frame = self.frame(&diagnostics)?;
        self.diagnostics.insert(keys[0], diagnostic_frame);
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
    rows: Vec<AssessedRow>,
}

/// The scan's state between chunks: record lines still needed and targets not yet finished.
struct Window<'a> {
    keys: &'a [u64],
    sequences: &'a [(u64, Variant)],
    samples: usize,
    policy: Policy,
    reference: &'a Reference,
    /// Kept record lines by global id, with the number of unfinished targets that use each.
    lines: HashMap<u32, (Vec<u8>, u32, u64)>,
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
    fn add(&mut self, scan: ChunkScan, records_before: u64) {
        let base = self.next_line;
        for i in 0..scan.offsets.len() - 1 {
            let line = scan.arena[scan.offsets[i] as usize..scan.offsets[i + 1] as usize].to_vec();
            self.lines
                .insert(base + i as u32, (line, 0, records_before + scan.source_records[i]));
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
        // Bound temporary quantitative genotype allocations independently of cohort width.
        let limit = self.block_limit();
        if upto - from > limit {
            for end in (from..upto).step_by(limit).map(|start| (start + limit).min(upto)) {
                self.finish(code, end)?;
            }
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
        let results: Vec<Result<(AssessedRow, [u64; State::ALL.len()])>> = work
            .par_iter()
            .map(|(t, ids)| {
                let recs: Vec<&[u8]> = ids.iter().map(|id| lines[id].0.as_slice()).collect();
                let result = if policy.dosage_field == crate::dosage::Field::Gt {
                    target_row(keys[*t], sequences, &recs, samples, policy, &fetch)
                } else {
                    quantitative_target_row(keys[*t], sequences, &recs, samples, policy, &fetch)
                };
                result
                    .map(|(mut row, states)| {
                        row.diagnostics.source_records = ids.iter().map(|id| lines[id].2).collect();
                        (row, states)
                    })
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

    fn block_limit(&self) -> usize {
        if self.policy.dosage_field == crate::dosage::Field::Gt {
            (8_388_608 / self.samples.max(1)).clamp(1, BLOCK_TARGETS)
        } else {
            (262_144 / self.samples.max(1)).clamp(1, BLOCK_TARGETS)
        }
    }

    /// Write every full block of the contig's finished rows (and the remainder when the contig is done).
    fn write_blocks(&mut self, code: u8, all: bool) -> Result<()> {
        let c = code as usize;
        let limit = self.block_limit();
        let Some(rows) = self.rows[c].as_mut() else {
            return Ok(());
        };
        while rows.rows.len() >= limit || all && !rows.rows.is_empty() {
            let n = rows.rows.len().min(limit);
            let block: Vec<AssessedRow> = rows.rows.drain(..n).collect();
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
) -> Result<(AssessedRow, [u64; State::ALL.len()])> {
    let mut states = [0u64; State::ALL.len()];
    if lines.is_empty() {
        states[State::UnknownNoRecord as usize] = samples as u64;
        return Ok((
            AssessedRow {
                codes: missing_row(samples),
                diagnostics: Diagnostics {
                    missing: MissingStates::Uniform(State::UnknownNoRecord as u8),
                    source_records: Vec::new(),
                },
                ..Default::default()
            },
            states,
        ));
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
    let mut missing = vec![PASSING; samples];
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
        if c == MISSING {
            missing[s] = state;
        }
    }
    Ok((
        AssessedRow {
            codes: row,
            diagnostics: Diagnostics {
                missing: MissingStates::from_codes(missing),
                source_records: Vec::new(),
            },
            ..Default::default()
        },
        states,
    ))
}

/// `target_row` for the common case of one record: the sample fields are walked in place and each distinct
/// field is assessed once.
fn single_record_row(
    line: &[u8],
    target: &genotype::Target<'_>,
    samples: usize,
    policy: &Policy,
    fetch: &(impl Fn(u64, usize) -> Option<String> + Sync),
) -> Result<(AssessedRow, [u64; State::ALL.len()])> {
    let mut states = [0u64; State::ALL.len()];
    let mut tabs = memchr::memchr_iter(b'\t', line);
    let shared_end = tabs
        .nth(8)
        .ok_or_else(|| Error::Invalid("record has too few fields".into()))?;
    let shared = std::str::from_utf8(&line[..shared_end]).map_err(|_| Error::Invalid("invalid UTF-8".into()))?;
    let mut row = vec![0u8; row_bytes(samples)];
    let mut missing = vec![PASSING; samples];
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
        if c == MISSING {
            missing[s] = state;
        }
    }
    if count != samples {
        return invalid!("record has {count} sample fields, expected {samples}");
    }
    Ok((
        AssessedRow {
            codes: row,
            diagnostics: Diagnostics {
                missing: MissingStates::from_codes(missing),
                source_records: Vec::new(),
            },
            ..Default::default()
        },
        states,
    ))
}

/// Quantitative rows retain exact measurements alongside packed hard calls. Reference blocks use GT.
fn quantitative_target_row(
    key: u64,
    sequences: &[(u64, Variant)],
    lines: &[&[u8]],
    samples: usize,
    policy: &Policy,
    fetch: &(impl Fn(u64, usize) -> Option<String> + Sync),
) -> Result<(AssessedRow, [u64; State::ALL.len()])> {
    let mut states = [0; State::ALL.len()];
    let mut row = AssessedRow {
        codes: missing_row(samples),
        ..Default::default()
    };
    if lines.is_empty() {
        states[State::UnknownNoRecord as usize] = samples as u64;
        row.diagnostics.missing = MissingStates::Uniform(State::UnknownNoRecord as u8);
        return Ok((row, states));
    }
    let mut missing = vec![PASSING; samples];
    let (mut r, mut a) = (String::new(), String::new());
    let target = target_for(key, sequences, &mut r, &mut a)?;
    let split = lines
        .iter()
        .map(|line| {
            let text = std::str::from_utf8(line).map_err(|_| Error::Invalid("invalid UTF-8".into()))?;
            let fields: Vec<_> = text.split('\t').collect();
            if fields.len() != 9 + samples {
                return invalid!("record has {} fields, expected {}", fields.len(), 9 + samples);
            }
            Ok((fields[..9].join("\t"), fields[9..].to_vec()))
        })
        .collect::<Result<Vec<_>>>()?;
    // Many samples share FORMAT strings; memoize assessed calls within a target.
    let mut memo: HashMap<Vec<&str>, genotype::Call> = HashMap::new();
    for sample in 0..samples {
        let fields: Vec<_> = split.iter().map(|(_, fields)| fields[sample]).collect();
        let call = if let Some(call) = memo.get(&fields) {
            call.clone()
        } else {
            let records = split
                .iter()
                .zip(&fields)
                .map(|((shared, _), field)| gvcf::parse_record(&format!("{shared}\t{field}")))
                .collect::<Result<Vec<_>>>()?;
            let call = genotype::assess(&records, &target, policy, fetch);
            memo.insert(fields, call.clone());
            call
        };
        states[call.state as usize] += 1;
        if !call.state.is_passing() {
            missing[sample] = call.state as u8;
        }
        if call.state.is_passing() {
            if let Some(m) = call.measurement {
                row.measurements.push((sample, *m));
            } else if let Some(d) = call.alt_dosage {
                let shift = (sample % 4) * 2;
                row.codes[sample / 4] = (row.codes[sample / 4] & !(3 << shift)) | (d << shift);
            }
        }
    }
    row.diagnostics.missing = MissingStates::from_codes(missing);
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
    let (set, _) = crate::targets::targets(packs, identity, options.targets_cache, false)?;
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
        quantitative: options.dosage_field != crate::dosage::Field::Gt,
        out: BufWriter::with_capacity(1 << 22, file),
        offset: 0,
        blocks: Vec::new(),
        diagnostics: BTreeMap::new(),
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
            let records_before = collector.records_scanned();
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
                    .with_merge_split(options.merge_split_records)
                    .with_dosage_field(options.dosage_field),
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
            w.add(scan, records_before);
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
    crate::extract::check_assembly(vcf, &header, reference)?;
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
    let diagnostics: Vec<_> = writer
        .blocks
        .iter()
        .map(|b| writer.diagnostics[&b.first_key].clone())
        .collect();
    let body_sha256 = digest_of_frames(
        std::iter::once(sequence_frame.sha256.as_str())
            .chain(writer.blocks.iter().map(|b| b.frame.sha256.as_str()))
            .chain(diagnostics.iter().map(|d| d.sha256.as_str())),
    );
    let policy = Policy::for_gvcf(
        policy_refcall,
        options.haploid_xy_as_homozygous,
        options.accept_missing_quality,
    )
    .with_merge_split(options.merge_split_records)
    .with_dosage_field(options.dosage_field);
    let header = Header {
        schema: DIAGNOSTIC_SCHEMA.into(),
        pgsum_version: env!("CARGO_PKG_VERSION").into(),
        policy: PolicyInfo {
            id: if policy.dosage_field == crate::dosage::Field::Gt {
                policy.id.into()
            } else {
                format!("{}-{}-strict-v1", policy.id, policy.dosage_field.label().to_lowercase())
            },
            dosage_field: policy.dosage_field,
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
        diagnostics,
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
type Measurements = BTreeMap<usize, Vec<(usize, crate::dosage::Measurement)>>;
type Block = (Vec<u64>, Vec<u8>, Measurements);

/// Default retained decoded-block payload budget. In-use handles and decoder workspace are additional.
pub const DEFAULT_CACHE_BYTES: usize = 128 * 1024 * 1024;

#[derive(Clone, Debug, Default, Serialize)]
pub struct CacheStats {
    pub capacity_bytes: usize,
    pub resident_bytes: usize,
    pub resident_blocks: usize,
    pub largest_block_bytes: usize,
    pub loads: u64,
    pub hits: u64,
    pub evictions: u64,
    pub streaming_peak_bytes: usize,
    pub streaming_peak_blocks: usize,
    pub key_prefix_reads: u64,
}

struct BlockCache {
    entries: std::collections::VecDeque<(usize, Arc<Block>, usize)>,
    stats: CacheStats,
}

// Count vector/string capacities and estimate map storage with two entries' space per row. This is a
// decoded payload estimate, not a promise about allocator metadata or total process RSS.
fn block_bytes(block: &Block) -> usize {
    let mut bytes = std::mem::size_of::<Block>() + block.0.capacity() * 8 + block.1.capacity();
    bytes += block.2.len() * 2 * std::mem::size_of::<(usize, Vec<(usize, crate::dosage::Measurement)>)>();
    for entries in block.2.values() {
        bytes += entries.capacity() * std::mem::size_of::<(usize, crate::dosage::Measurement)>();
        for (_, measurement) in entries {
            bytes += match measurement {
                crate::dosage::Measurement::DS(d) => d.coefficient.capacity(),
                crate::dosage::Measurement::GP(p) => p.iter().map(|d| d.coefficient.capacity()).sum(),
            };
        }
    }
    bytes
}

/// Shared ownership lets cache eviction release blocks as soon as active readers finish with them.
#[derive(Clone)]
pub struct CohortBlock(Arc<Block>);
impl CohortBlock {
    pub fn keys(&self) -> &[u64] {
        &self.0.0
    }
    pub fn rows(&self) -> &[u8] {
        &self.0.1
    }
}

#[derive(Clone)]
pub struct CohortRow {
    block: Arc<Block>,
    index: usize,
    width: usize,
}
impl std::ops::Deref for CohortRow {
    type Target = [u8];
    fn deref(&self) -> &[u8] {
        &self.block.1[self.index * self.width..(self.index + 1) * self.width]
    }
}
impl CohortRow {
    pub fn measurements(&self) -> &[(usize, crate::dosage::Measurement)] {
        self.block.2.get(&self.index).map_or(&[], Vec::as_slice)
    }
}

/// A borrowed view into a worker's current block; no atomic reference-count update per term.
pub struct CohortRowRef<'a> {
    block: &'a Block,
    index: usize,
    width: usize,
}
impl std::ops::Deref for CohortRowRef<'_> {
    type Target = [u8];
    fn deref(&self) -> &[u8] {
        &self.block.1[self.index * self.width..(self.index + 1) * self.width]
    }
}
impl CohortRowRef<'_> {
    pub fn measurements(&self) -> &[(usize, crate::dosage::Measurement)] {
        self.block.2.get(&self.index).map_or(&[], Vec::as_slice)
    }
}

/// One reader per scoring worker keeps its current block without locking on each term.
pub struct CohortReader<'a> {
    table: &'a CohortTable,
    current: Option<(usize, Arc<Block>)>,
}
impl CohortReader<'_> {
    pub fn row(&mut self, key: u64) -> Result<Option<CohortRowRef<'_>>> {
        let Some(index) = self.table.block_index(key) else {
            return Ok(None);
        };
        if self.current.as_ref().is_none_or(|(i, _)| *i != index) {
            // Release the old worker handle before decoding a replacement.
            self.current = None;
            self.current = Some((index, self.table.block(index)?));
        }
        let block = &self.current.as_ref().expect("loaded block").1;
        Ok(block.0.binary_search(&key).ok().map(|row| CohortRowRef {
            block: block.as_ref(),
            index: row,
            width: self.table.row,
        }))
    }
}

struct DiagnosticBlock {
    keys: Vec<u64>,
    rows: Vec<Diagnostics>,
}

/// A cohort file, its blocks decompressed on first use.
pub struct CohortTable {
    pub header: Header,
    path: PathBuf,
    map: Mmap,
    cache: Mutex<BlockCache>,
    validated_blocks: Vec<AtomicBool>,
    validated_diagnostics: Vec<AtomicBool>,
    sequence_index: HashMap<(u8, Variant), u64>,
    row: usize,
    diagnostics_cache: Mutex<std::collections::VecDeque<(usize, Arc<DiagnosticBlock>)>>,
}

impl CohortTable {
    pub fn open(path: &Path) -> Result<CohortTable> {
        Self::open_with_cache_bytes(path, DEFAULT_CACHE_BYTES)
    }

    /// Bound retained decoded payload; one oversized block is retained to avoid repeated decoding.
    /// Active worker handles, diagnostic frames, temporary decoding buffers and score sums are extra.
    pub fn open_with_cache_bytes(path: &Path, capacity_bytes: usize) -> Result<CohortTable> {
        if capacity_bytes == 0 {
            return invalid!("cohort cache capacity must be positive");
        }
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
        if header.schema != SCHEMA && header.schema != QUANTITATIVE_SCHEMA && header.schema != DIAGNOSTIC_SCHEMA {
            return invalid!("{}: schema {} is not {SCHEMA}", path.display(), header.schema);
        }
        if (header.schema == DIAGNOSTIC_SCHEMA && header.diagnostics.len() != header.blocks.len())
            || (header.schema != DIAGNOSTIC_SCHEMA && !header.diagnostics.is_empty())
        {
            return Err(bad("diagnostic frame list"));
        }
        let digest = digest_of_frames(
            std::iter::once(header.sequence_frame.sha256.as_str())
                .chain(header.blocks.iter().map(|b| b.frame.sha256.as_str()))
                .chain(header.diagnostics.iter().map(|d| d.sha256.as_str())),
        );
        if digest != header.body_sha256 || header.blocks.windows(2).any(|w| w[0].last_key >= w[1].first_key) {
            return Err(bad("block list"));
        }
        if header.samples.is_empty() || header.samples.len() > u32::MAX as usize {
            return Err(bad("sample count"));
        }
        let row = row_bytes(header.samples.len());
        let mut table = CohortTable {
            path: path.to_owned(),
            cache: Mutex::new(BlockCache {
                entries: Default::default(),
                stats: CacheStats {
                    capacity_bytes,
                    ..Default::default()
                },
            }),
            validated_blocks: (0..header.blocks.len()).map(|_| AtomicBool::new(false)).collect(),
            validated_diagnostics: (0..header.diagnostics.len()).map(|_| AtomicBool::new(false)).collect(),
            sequence_index: HashMap::new(),
            row,
            diagnostics_cache: Mutex::new(std::collections::VecDeque::new()),
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

    fn decode_block(&self, i: usize) -> Result<Block> {
        let validated = self.validated_blocks[i].load(Ordering::Acquire);
        let result = (|| -> std::result::Result<Block, String> {
            let info = &self.header.blocks[i];
            let bytes = self.frame(&info.frame).map_err(|e| e.to_string())?;
            let n = u64::from_le_bytes(bytes.get(..8).ok_or("short block")?.try_into().expect("8")) as usize;
            let keys_end = n
                .checked_mul(8)
                .and_then(|v| v.checked_add(8))
                .ok_or("cohort block dimensions overflow")?;
            let rows_end = n
                .checked_mul(self.row)
                .and_then(|v| v.checked_add(keys_end))
                .ok_or("cohort block dimensions overflow")?;
            if n == 0
                || n as u64 != info.targets
                || bytes.len() < rows_end
                || (self.header.policy.dosage_field == crate::dosage::Field::Gt && bytes.len() != rows_end)
            {
                return Err(format!("{}: invalid cohort block {i}", self.path.display()));
            }
            let keys: Vec<u64> = bytes[8..keys_end]
                .as_chunks::<8>()
                .0
                .iter()
                .map(|&c| u64::from_le_bytes(c))
                .collect();
            if keys.first() != Some(&info.first_key)
                || keys.last() != Some(&info.last_key)
                || keys.windows(2).any(|w| w[0] >= w[1])
            {
                return Err("invalid cohort block keys".into());
            }
            let measurements: Measurements = if self.header.policy.dosage_field != crate::dosage::Field::Gt {
                let values: Vec<(usize, Vec<(usize, crate::dosage::Measurement)>)> =
                    serde_json::from_slice(&bytes[rows_end..]).map_err(|e| e.to_string())?;
                let mut result = Measurements::new();
                for (r, entries) in values {
                    if r >= n || result.contains_key(&r) {
                        return Err("invalid quantitative row".into());
                    }
                    let mut seen = std::collections::HashSet::new();
                    for (s, m) in &entries {
                        if *s >= self.header.samples.len()
                            || !seen.insert(*s)
                            || (!validated && !m.valid())
                            || !matches!(
                                (self.header.policy.dosage_field, m),
                                (crate::dosage::Field::Ds, crate::dosage::Measurement::DS(_))
                                    | (crate::dosage::Field::Gp, crate::dosage::Measurement::GP(_))
                            )
                            || code(&bytes[keys_end + r * self.row..keys_end + (r + 1) * self.row], *s) != MISSING
                        {
                            return Err("invalid quantitative sample".into());
                        }
                    }
                    result.insert(r, entries);
                }
                result
            } else {
                Measurements::new()
            };
            Ok((keys, bytes[keys_end..rows_end].to_vec(), measurements))
        })()
        .map_err(Error::Invalid)?;
        // Every reload still verifies the full frame digest before trusting prior numeric validation.
        self.validated_blocks[i].store(true, Ordering::Release);
        Ok(result)
    }

    fn block(&self, i: usize) -> Result<Arc<Block>> {
        let mut cache = self
            .cache
            .lock()
            .map_err(|_| Error::Invalid("cohort cache poisoned".into()))?;
        if let Some(at) = cache.entries.iter().position(|(index, _, _)| *index == i) {
            let entry = cache.entries.remove(at).expect("found block");
            let block = entry.1.clone();
            cache.entries.push_back(entry);
            cache.stats.hits += 1;
            return Ok(block);
        }
        // Serialize decoding so concurrent packs do not allocate duplicate temporary copies.
        let block = Arc::new(self.decode_block(i)?);
        let bytes = block_bytes(&block);
        while !cache.entries.is_empty() && cache.stats.resident_bytes.saturating_add(bytes) > cache.stats.capacity_bytes
        {
            let (_, _, removed) = cache.entries.pop_front().expect("nonempty cache");
            cache.stats.resident_bytes -= removed;
            cache.stats.evictions += 1;
        }
        cache.stats.resident_bytes += bytes;
        cache.stats.largest_block_bytes = cache.stats.largest_block_bytes.max(bytes);
        cache.stats.loads += 1;
        cache.entries.push_back((i, block.clone(), bytes));
        cache.stats.resident_blocks = cache.entries.len();
        Ok(block)
    }

    pub fn cache_stats(&self) -> Result<CacheStats> {
        Ok(self
            .cache
            .lock()
            .map_err(|_| Error::Invalid("cohort cache poisoned".into()))?
            .stats
            .clone())
    }

    pub fn reader(&self) -> CohortReader<'_> {
        CohortReader {
            table: self,
            current: None,
        }
    }

    /// Shared block handle; drop it when finished so evicted storage can be reclaimed.
    pub fn block_rows(&self, b: usize) -> Result<CohortBlock> {
        if b >= self.header.blocks.len() {
            return invalid!("cohort block index out of range");
        }
        Ok(CohortBlock(self.block(b)?))
    }

    /// Validate all genotype frames. Only the cache budget's worth of blocks stays resident.
    pub fn preload(&self) -> Result<()> {
        for i in 0..self.header.blocks.len() {
            self.block(i)?;
        }
        Ok(())
    }

    fn block_index(&self, key: u64) -> Option<usize> {
        let i = self.header.blocks.partition_point(|b| b.first_key <= key);
        (i != 0 && key <= self.header.blocks[i - 1].last_key).then(|| i - 1)
    }

    pub fn row(&self, key: u64) -> Result<Option<CohortRow>> {
        let Some(i) = self.block_index(key) else {
            return Ok(None);
        };
        let block = self.block(i)?;
        Ok(block.0.binary_search(&key).ok().map(|index| CohortRow {
            block,
            index,
            width: self.row,
        }))
    }

    /// Only called after successful full validation; input mutation while mapped is unsupported.
    fn validated_key_prefix(&self, block: usize) -> Result<Vec<u64>> {
        let info = &self.header.blocks[block];
        let start = 8usize
            .checked_add(info.frame.offset as usize)
            .ok_or_else(|| Error::Invalid("invalid cohort offset".into()))?;
        let end = start
            .checked_add(info.frame.compressed_bytes as usize)
            .ok_or_else(|| Error::Invalid("invalid cohort size".into()))?;
        let compressed = self
            .map
            .get(start..end)
            .ok_or_else(|| Error::Invalid("invalid cohort frame".into()))?;
        let mut decoder = zstd::Decoder::new(compressed).map_err(Error::io(&self.path))?;
        let mut count = [0u8; 8];
        decoder.read_exact(&mut count).map_err(Error::io(&self.path))?;
        if u64::from_le_bytes(count) != info.targets {
            return invalid!("cohort block target count changed");
        }
        let len = usize::try_from(info.targets)
            .ok()
            .and_then(|n| n.checked_mul(8))
            .ok_or_else(|| Error::Invalid("cohort key dimensions overflow".into()))?;
        let mut bytes = vec![0; len];
        decoder.read_exact(&mut bytes).map_err(Error::io(&self.path))?;
        self.cache
            .lock()
            .map_err(|_| Error::Invalid("cohort cache poisoned".into()))?
            .stats
            .key_prefix_reads += 1;
        Ok(bytes
            .as_chunks::<8>()
            .0
            .iter()
            .map(|b| u64::from_le_bytes(*b))
            .collect())
    }

    fn diagnostic_block(&self, block: usize) -> Result<Arc<DiagnosticBlock>> {
        let mut cache = self
            .diagnostics_cache
            .lock()
            .map_err(|_| Error::Invalid("diagnostic cache poisoned".into()))?;
        if let Some(i) = cache.iter().position(|(b, _)| *b == block) {
            let entry = cache.remove(i).expect("found entry");
            let rows = entry.1.clone();
            cache.push_back(entry);
            return Ok(rows);
        }
        let bytes = self.frame(&self.header.diagnostics[block])?;
        let rows: Vec<Diagnostics> =
            serde_json::from_slice(&bytes).map_err(|e| Error::Invalid(format!("invalid cohort diagnostics: {e}")))?;
        let keys = if self.validated_diagnostics[block].load(Ordering::Acquire) {
            // The immutable frame has already been fully authenticated and cross-checked. Only its
            // key prefix is needed when rereading evicted diagnostic metadata; do not parse GP again.
            self.validated_key_prefix(block)?
        } else {
            let calls = self.block(block)?;
            let (keys, codes, measurements) = calls.as_ref();
            let samples = self.header.samples.len();
            if rows.len() != keys.len() {
                return invalid!("diagnostic row count does not match cohort block");
            }
            for (r, diagnostic) in rows.iter().enumerate() {
                if !diagnostic.missing.validate(samples)
                    || diagnostic
                        .source_records
                        .iter()
                        .any(|n| *n == 0 || *n > self.header.records_scanned)
                    || diagnostic.source_records.windows(2).any(|w| w[0] >= w[1])
                {
                    return invalid!("invalid cohort diagnostic state or source ordinal");
                }
                let measured: std::collections::HashSet<_> =
                    measurements.get(&r).into_iter().flatten().map(|(s, _)| *s).collect();
                let row = &codes[r * self.row..(r + 1) * self.row];
                for sample in 0..samples {
                    let missing = code(row, sample) == MISSING && !measured.contains(&sample);
                    if missing != diagnostic.missing.at(sample).is_some() {
                        return invalid!("cohort diagnostic state disagrees with packed call");
                    }
                }
            }
            self.validated_diagnostics[block].store(true, Ordering::Release);
            keys.clone()
        };
        let rows = Arc::new(DiagnosticBlock { keys, rows });
        if cache.len() == 8 {
            cache.pop_front();
        }
        cache.push_back((block, rows.clone()));
        Ok(rows)
    }

    /// Exact assessed call state. Legacy files return None for missing cells whose reasons were not retained.
    pub fn call_state(&self, key: u64, sample: usize) -> Result<Option<State>> {
        if sample >= self.header.samples.len() {
            return invalid!("cohort sample index out of range");
        }
        let Some(row) = self.row(key)? else {
            return Ok(None);
        };
        let c = code(&row, sample);
        if c != MISSING {
            return Ok(Some(if c == 0 {
                State::ObservedReference
            } else {
                State::ObservedVariant
            }));
        }
        if let Some((_, m)) = row.measurements().iter().find(|(s, _)| *s == sample) {
            return Ok(Some(if m.alt_dosage().coefficient == "0" {
                State::ObservedReference
            } else {
                State::ObservedVariant
            }));
        }
        self.missing_state(key, sample)
    }

    /// The caller has already checked that this is an existing missing cell, with no measurement.
    pub(crate) fn missing_state(&self, key: u64, sample: usize) -> Result<Option<State>> {
        if self.header.schema != DIAGNOSTIC_SCHEMA {
            return Ok(None);
        }
        let block = self.block_index(key).expect("existing row");
        let diagnostics = self.diagnostic_block(block)?;
        let row = diagnostics.keys.binary_search(&key).expect("existing row");
        Ok(diagnostics.rows[row].missing.at(sample))
    }

    /// Source references are one-based record ordinals, excluding headers, in the input identified by its hash.
    /// None means diagnostics were not retained; Some(empty) means there was no overlapping record.
    pub fn source_records(&self, key: u64) -> Result<Option<Vec<u64>>> {
        if self.header.schema != DIAGNOSTIC_SCHEMA {
            return Ok(None);
        }
        let Some(block) = self.block_index(key) else {
            return Ok(None);
        };
        let diagnostics = self.diagnostic_block(block)?;
        Ok(diagnostics
            .keys
            .binary_search(&key)
            .ok()
            .map(|row| diagnostics.rows[row].source_records.clone()))
    }

    /// Check every diagnostic frame without retaining more than eight at a time.
    pub fn validate_diagnostics(&self) -> Result<()> {
        for block in 0..self.header.diagnostics.len() {
            self.diagnostic_block(block)?;
        }
        Ok(())
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
    pub expected_genotype_terms: Vec<u64>,
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

fn check_pack(pack: &Pack, cohort: &CohortTable) -> Result<()> {
    if !cohort
        .header
        .packs
        .iter()
        .any(|p| p.pgs_id == pack.header.pgs_id && p.records_sha256 == pack.header.records_sha256)
    {
        return invalid!(
            "the cohort file was not extracted with this {} pack",
            pack.header.pgs_id
        );
    }
    Ok(())
}

fn term_key(term: &crate::pack::TermRecord, cohort: &CohortTable, options: &Options) -> Result<Option<(u64, bool)>> {
    let (key, alt) = match plan(term, options).0 {
        Plan::Unscorable(_) => return Ok(None),
        Plan::Snv { key, effect_is_alt } => (Some(key), effect_is_alt),
        Plan::Sequence { variant, effect_is_alt } => (cohort.sequence_key(term.contig, &variant), effect_is_alt),
    };
    key.map(|k| Some((k, alt))).ok_or_else(|| {
        Error::Invalid(format!(
            "the cohort file has no call for {}:{}; extract with this pack",
            CONTIGS[term.contig as usize - 1],
            term.pos
        ))
    })
}

struct Accumulator {
    baseline: ExactSum,
    sums: Vec<ExactSum>,
    scorable: Vec<u64>,
    expected_genotype_terms: Vec<u64>,
    effect_scorable: Vec<f64>,
    base_scorable: u64,
    base_effect: f64,
    effect_all: f64,
    total: u64,
    base_exponents: BTreeMap<i64, u64>,
    own_min: Vec<i64>,
    left: Vec<Vec<(i64, u64)>>,
}
impl Accumulator {
    fn new(n: usize) -> Self {
        Self {
            baseline: Default::default(),
            sums: vec![ExactSum::default(); n],
            scorable: vec![0; n],
            expected_genotype_terms: vec![0; n],
            effect_scorable: vec![0.0; n],
            base_scorable: 0,
            base_effect: 0.0,
            effect_all: 0.0,
            total: 0,
            base_exponents: BTreeMap::new(),
            own_min: vec![0; n],
            left: vec![Vec::new(); n],
        }
    }
    fn count(&mut self, term: &crate::pack::TermRecord) -> Option<f64> {
        self.total += 1;
        let effect = term_effect(term);
        if let Some(e) = effect {
            self.effect_all += e;
        }
        effect
    }
    fn add(
        &mut self,
        term: &crate::pack::TermRecord,
        row: &[u8],
        measurements: &[(usize, crate::dosage::Measurement)],
        effect_is_alt: bool,
    ) {
        let n = self.sums.len();
        let effect = self.count(term);
        // Contribution at each code (ALT dosage); None where the model gives no weight.
        let at = |c: u8| term_contribution(term, if effect_is_alt { c } else { 2 - c });
        let contributions = [at(0), at(1), at(2)];
        let e = effect.unwrap_or(0.0);
        // Code 0 goes into the baseline when it has a contribution.
        let base = contributions[0].clone();
        if let Some(c) = &base {
            self.baseline.add(c);
            self.base_scorable += 1;
            self.base_effect += e;
            *self.base_exponents.entry(c.exponent).or_default() += 1;
        }
        let negative_base = base.as_ref().map(|c| c.negated());
        let differences = contributions
            .each_ref()
            .map(|c| c.as_ref().zip(base.as_ref()).and_then(|(c, b)| c.small_difference(b)));
        for (byte_index, &byte) in row.iter().enumerate() {
            if byte == 0 {
                // Four homozygous-reference samples, already included in the baseline.
                continue;
            }
            let start = byte_index * 4;
            let mut codes = byte;
            for s in start..(start + 4).min(n) {
                let code = codes & 3;
                codes >>= 2;
                if code == 0 {
                    continue;
                }
                if code != MISSING
                    && let Some(delta) = &differences[code as usize]
                {
                    self.sums[s].add(delta);
                    // Both contributions have the same exponent: retain its baseline count.
                    // Keep the previous floating-point operation order for byte-identical coverage.
                    self.effect_scorable[s] -= e;
                    self.effect_scorable[s] += e;
                    continue;
                }
                // Wider coefficients, different dosage-weight exponents, or missing calls use
                // the general path, including precision tracking when a baseline term is removed.
                if let Some(c) = &negative_base {
                    self.sums[s].add(c);
                    self.scorable[s] = self.scorable[s].wrapping_sub(1);
                    self.effect_scorable[s] -= e;
                    match self.left[s].iter_mut().find(|(x, _)| *x == c.exponent) {
                        Some((_, k)) => *k += 1,
                        None => self.left[s].push((c.exponent, 1)),
                    }
                }
                if code != MISSING
                    && let Some(c) = &contributions[code as usize]
                {
                    self.sums[s].add(c);
                    self.scorable[s] = self.scorable[s].wrapping_add(1);
                    self.effect_scorable[s] += e;
                    self.own_min[s] = self.own_min[s].min(c.exponent);
                }
            }
        }
        {
            for (sample, measurement) in measurements {
                if let Some(c) = measurement.contribution(term.model, &term.weights, effect_is_alt) {
                    self.sums[*sample].add(&c);
                    self.expected_genotype_terms[*sample] += 1;
                    self.scorable[*sample] = self.scorable[*sample].wrapping_add(1);
                    self.effect_scorable[*sample] += e;
                    self.own_min[*sample] = self.own_min[*sample].min(c.exponent);
                }
            }
        }
    }
    fn finish(mut self, pgs_id: &str) -> CohortScore {
        let n = self.sums.len();
        let mut exponents = Vec::with_capacity(n);
        for s in 0..n {
            self.sums[s].merge(self.baseline.clone());
            self.scorable[s] = self.scorable[s].wrapping_add(self.base_scorable);
            self.effect_scorable[s] += self.base_effect;
            let kept_base = self
                .base_exponents
                .iter()
                .filter(|(e, count)| {
                    let gone = self.left[s].iter().find(|(x, _)| x == *e).map_or(0, |(_, k)| *k);
                    **count > gone
                })
                .map(|(e, _)| *e)
                .min()
                .unwrap_or(0);
            exponents.push(self.own_min[s].min(kept_base).min(0));
        }
        CohortScore {
            pgs_id: pgs_id.into(),
            total_terms: self.total,
            sums: self.sums,
            scorable: self.scorable,
            expected_genotype_terms: self.expected_genotype_terms,
            effect_scorable: self.effect_scorable,
            effect_all: self.effect_all,
            exponents,
        }
    }
}

/// Score a pack for every sample: the same terms, rules and exact arithmetic as `score::score`, per sample.
///
/// A term's contribution at code 0 is added once to a shared baseline. Packed groups of four code-0
/// samples are skipped together; other samples apply an exact precomputed difference where possible.
/// Code 0 is the homozygous-reference dosage, not necessarily the most common code.
pub fn score_cohort(pack: &Pack, cohort: &CohortTable, options: &Options) -> Result<CohortScore> {
    check_pack(pack, cohort)?;
    let mut sums = Accumulator::new(cohort.header.samples.len());
    let mut reader = cohort.reader();
    for term in pack.terms_from(0) {
        let term = term?;
        if let Some((key, effect_is_alt)) = term_key(&term, cohort, options)? {
            let row = reader
                .row(key)?
                .ok_or_else(|| Error::Invalid("cohort scoring target is absent; extract with this pack".into()))?;
            sums.add(&term, &row, row.measurements(), effect_is_alt);
        } else {
            sums.count(&term);
        }
    }
    Ok(sums.finish(&pack.header.pgs_id))
}

/// Result counts distinguish shared ordered passes from the source-order fallback.
pub struct CohortBatch {
    pub scores: Vec<CohortScore>,
    pub shared_packs: usize,
    pub independent_packs: usize,
}

struct PreparedTerm {
    term: crate::pack::TermRecord,
    key: u64,
    effect_is_alt: bool,
    block: usize,
}
struct PendingScore<'a> {
    pack: &'a Pack,
    terms: crate::pack::TermIter<'a>,
    next: Option<PreparedTerm>,
    sums: Accumulator,
}
impl PendingScore<'_> {
    fn advance(&mut self, cohort: &CohortTable, options: &Options) -> Result<()> {
        self.next = None;
        for term in self.terms.by_ref() {
            let term = term?;
            if let Some((key, effect_is_alt)) = term_key(&term, cohort, options)? {
                let block = cohort
                    .block_index(key)
                    .ok_or_else(|| Error::Invalid("cohort scoring target is absent; extract with this pack".into()))?;
                self.next = Some(PreparedTerm {
                    term,
                    key,
                    effect_is_alt,
                    block,
                });
                break;
            }
            self.sums.count(&term);
        }
        Ok(())
    }
}

fn ordered_blocks(pack: &Pack, cohort: &CohortTable, options: &Options) -> Result<Option<Vec<usize>>> {
    check_pack(pack, cohort)?;
    let mut blocks = Vec::new();
    for term in pack.terms_from(0) {
        if let Some((key, _)) = term_key(&term?, cohort, options)? {
            let block = cohort
                .block_index(key)
                .ok_or_else(|| Error::Invalid("cohort scoring target is absent; extract with this pack".into()))?;
            if blocks.last().is_some_and(|b| *b > block) {
                return Ok(None);
            }
            if blocks.last() != Some(&block) {
                blocks.push(block);
            }
        }
    }
    Ok(Some(blocks))
}

/// At most eight packs' accumulators are active together. Every pack retains its original term order.
pub fn score_cohorts(packs: &[Pack], cohort: &CohortTable, options: &Options) -> Result<CohortBatch> {
    // Packed GT blocks are cheap and scoring dominates; avoid an extra term-order pass there.
    if cohort.header.policy.dosage_field == crate::dosage::Field::Gt {
        let scores = packs
            .par_iter()
            .map(|p| score_cohort(p, cohort, options))
            .collect::<Result<Vec<_>>>()?;
        return Ok(CohortBatch {
            scores,
            shared_packs: 0,
            independent_packs: packs.len(),
        });
    }

    let mut result = CohortBatch {
        scores: Vec::new(),
        shared_packs: 0,
        independent_packs: 0,
    };
    for group in packs.chunks(8) {
        let mut slots: Vec<Option<CohortScore>> = (0..group.len()).map(|_| None).collect();
        let mut ordered = Vec::new();
        let mut independent = Vec::new();
        let mut needed = std::collections::BTreeSet::new();
        for (i, pack) in group.iter().enumerate() {
            if let Some(blocks) = ordered_blocks(pack, cohort, options)? {
                needed.extend(blocks);
                ordered.push(i);
            } else {
                independent.push(i);
            }
        }
        result.shared_packs += ordered.len();
        result.independent_packs += independent.len();
        let mut pending = ordered
            .iter()
            .map(|i| {
                let pack = &group[*i];
                let mut pending = PendingScore {
                    pack,
                    terms: pack.terms_from(0),
                    next: None,
                    sums: Accumulator::new(cohort.header.samples.len()),
                };
                pending.advance(cohort, options)?;
                Ok(pending)
            })
            .collect::<Result<Vec<_>>>()?;
        let capacity = cohort.cache_stats()?.capacity_bytes;
        let mut blocks = needed.into_iter().peekable();
        let mut largest = 0usize;
        while blocks.peek().is_some() {
            // Start with one block to measure decoded size, then bound each parallel wave by the
            // cache target and available threads. Decoder workspace remains additional.
            let width = capacity
                .checked_div(largest)
                .unwrap_or(1)
                .max(1)
                .min(rayon::current_num_threads());
            let indexes: Vec<_> = blocks.by_ref().take(width).collect();
            let decoded = indexes
                .par_iter()
                .map(|i| {
                    let block = Arc::new(cohort.decode_block(*i)?);
                    let bytes = block_bytes(&block);
                    Ok((*i, block, bytes))
                })
                .collect::<Result<Vec<_>>>()?;
            let bytes: usize = decoded.iter().map(|(_, _, bytes)| bytes).sum();
            largest = largest.max(decoded.iter().map(|(_, _, bytes)| *bytes).max().unwrap_or(0));
            {
                let mut cache = cohort
                    .cache
                    .lock()
                    .map_err(|_| Error::Invalid("cohort cache poisoned".into()))?;
                cache.stats.loads += decoded.len() as u64;
                cache.stats.largest_block_bytes = cache.stats.largest_block_bytes.max(largest);
                cache.stats.streaming_peak_bytes = cache.stats.streaming_peak_bytes.max(bytes);
                cache.stats.streaming_peak_blocks = cache.stats.streaming_peak_blocks.max(decoded.len());
            }
            pending.par_iter_mut().try_for_each(|pending| -> Result<()> {
                while let Some(term) = &pending.next {
                    let Ok(i) = decoded.binary_search_by_key(&term.block, |(b, _, _)| *b) else {
                        break;
                    };
                    let term = pending.next.take().expect("prepared term");
                    let block = &decoded[i].1;
                    let index = block.0.binary_search(&term.key).map_err(|_| {
                        Error::Invalid("cohort scoring target is absent; extract with this pack".into())
                    })?;
                    let row = CohortRowRef {
                        block,
                        index,
                        width: cohort.row,
                    };
                    pending
                        .sums
                        .add(&term.term, &row, row.measurements(), term.effect_is_alt);
                    pending.advance(cohort, options)?;
                }
                Ok(())
            })?;
        }
        for (i, pending) in ordered.into_iter().zip(pending) {
            if pending.next.is_some() {
                return invalid!("cohort shared pass did not consume every term");
            }
            slots[i] = Some(pending.sums.finish(&pending.pack.header.pgs_id));
        }
        let scored = independent
            .par_iter()
            .map(|i| Ok((*i, score_cohort(&group[*i], cohort, options)?)))
            .collect::<Result<Vec<_>>>()?;
        for (i, score) in scored {
            slots[i] = Some(score);
        }
        result
            .scores
            .extend(slots.into_iter().map(|s| s.expect("every pack was scored")));
    }
    Ok(result)
}
