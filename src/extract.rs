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

use crate::alleles::Variant;
use crate::digest::HashingReader;
use crate::genotype::{self, Policy, Target};
use crate::genotypes::{CompactCall, GenotypeTable, is_sequence_key, key_position, unpack_key};
use crate::gvcf::{self, HeaderFacts};
use crate::pack::ReferenceIdentity;
use crate::reference::Reference;
use crate::term::CONTIGS;
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

/// Scan the gVCF and collect the records overlapping each target. `ends[i]` is the last base of target `i`
/// (its position for SNVs; indel targets span their REF allele).
pub fn scan(gvcf: &Path, keys: &[u64], ends: &[u64], threads: usize) -> Result<Extracted> {
    let max_span = keys
        .iter()
        .zip(ends)
        .map(|(&k, &e)| e - key_position(k).1 as u64)
        .max()
        .unwrap_or(0);
    let file = File::open(gvcf).map_err(Error::io(gvcf))?;
    let hashing = HashingReader::new(file);
    let workers = NonZero::new(threads.max(1)).expect("at least one");
    let mut reader = std::io::BufReader::with_capacity(
        1 << 20,
        noodles_bgzf::io::MultithreadedReader::with_worker_count(workers, hashing),
    );
    let at = |line_no: u64, e: Error| Error::Invalid(format!("{}: line {line_no}: {e}", gvcf.display()));

    // Targets grouped by contig code: the index range into `keys`.
    let mut contig_ranges = [(0usize, 0usize); 26];
    for code in 1..=25u8 {
        let lo = keys.partition_point(|&k| unpack_key(k).0 < code);
        let hi = keys.partition_point(|&k| unpack_key(k).0 <= code);
        contig_ranges[code as usize] = (lo, hi);
    }
    let contig_codes: HashMap<&str, u8> = CONTIGS.iter().enumerate().map(|(i, c)| (*c, i as u8 + 1)).collect();

    let mut header = HeaderFacts::default();
    let mut in_header = true;
    let mut line = Vec::with_capacity(1 << 16);
    let mut line_no = 0u64;
    let mut records_scanned = 0u64;
    let mut arena = Vec::new();
    let mut offsets = vec![0u64];
    let mut first = vec![u32::MAX; keys.len()];
    let mut extra = Extra::new();
    let mut seen = [false; 26];
    let mut current: (u8, u64) = (0, 0); // contig code, last POS
    let mut cursor = 0usize; // first target of the current contig not yet passed
    let mut canonical_records = 0u64;
    let mut last_chrom: Vec<u8> = Vec::new();
    let mut last_code: Option<u8> = None;
    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line).map_err(Error::io(gvcf))?;
        if n == 0 {
            break;
        }
        line_no += 1;
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        let text = std::str::from_utf8(&line).map_err(|_| at(line_no, Error::Invalid("invalid UTF-8".into())))?;
        if in_header {
            if !text.starts_with('#') {
                return Err(at(
                    line_no,
                    Error::Invalid("record before the #CHROM header line".into()),
                ));
            }
            in_header = !header.read_line(text).map_err(|e| at(line_no, e))?;
            continue;
        }
        records_scanned += 1;
        let span = gvcf::interval(text).map_err(|e| at(line_no, e))?;
        // Records arrive grouped by contig, so the name lookup only runs when the contig changes.
        if span.chrom.as_bytes() != last_chrom.as_slice() {
            last_chrom.clear();
            last_chrom.extend_from_slice(span.chrom.as_bytes());
            last_code = contig_codes.get(span.chrom).copied();
        }
        let Some(code) = last_code else {
            continue;
        };
        canonical_records += 1;
        if code != current.0 {
            if seen[code as usize] {
                return Err(at(
                    line_no,
                    Error::Invalid(format!("{} records are not contiguous", span.chrom)),
                ));
            }
            seen[code as usize] = true;
            current = (code, 0);
            cursor = contig_ranges[code as usize].0;
        }
        if span.pos < current.1 {
            return Err(at(line_no, Error::Invalid("records are not sorted by position".into())));
        }
        current.1 = span.pos;
        let end_of_contig = contig_ranges[code as usize].1;
        // Targets ending before this record's POS can't overlap it or any later record; with spans up to
        // `max_span`, any target starting more than `max_span` before it has ended.
        while cursor < end_of_contig && key_position(keys[cursor]).1 as u64 + max_span < span.pos {
            cursor += 1;
        }
        let mut t = cursor;
        let mut kept: Option<u32> = None;
        while t < end_of_contig && (key_position(keys[t]).1 as u64) <= span.end {
            if ends[t] >= span.pos {
                let id = *kept.get_or_insert_with(|| {
                    arena.extend_from_slice(&line);
                    offsets.push(arena.len() as u64);
                    (offsets.len() - 2) as u32
                });
                if first[t] == u32::MAX {
                    first[t] = id;
                } else {
                    extra.entry(t as u32).or_default().push(id);
                }
            }
            t += 1;
        }
    }
    if in_header {
        return invalid!("{}: no #CHROM header line", gvcf.display());
    }
    if canonical_records == 0 {
        return invalid!(
            "{}: no records on chr1–chr22, chrX, chrY or chrM (are contigs named chr1, …?)",
            gvcf.display()
        );
    }
    // Finish hashing the compressed file (the reader may stop before the final empty BGZF block).
    let mut hashing = reader.into_inner().finish().map_err(Error::io(gvcf))?;
    std::io::copy(&mut hashing, &mut std::io::sink()).map_err(Error::io(gvcf))?;
    let gvcf_bytes = hashing.bytes;
    Ok(Extracted {
        header,
        gvcf_sha256: hashing.finish(),
        gvcf_bytes,
        records_scanned,
        arena,
        offsets,
        first,
        extra,
    })
}

/// Assess every target from its records, on `threads` threads.
pub fn assess(
    keys: &[u64],
    sequences: &[(u64, Variant)],
    scanned: &Extracted,
    reference: &Reference,
    threads: usize,
) -> Result<Vec<CompactCall>> {
    let policy = Policy {
        refcall_is_reference: scanned.header.refcall_defined,
        ..Policy::default()
    };
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
                let policy = &policy;
                scope.spawn(move || {
                    let mut out = Vec::with_capacity(end - start);
                    for (i, &key) in keys.iter().enumerate().take(end).skip(start) {
                        let (code, pos) = key_position(key);
                        let contig = CONTIGS[code as usize - 1];
                        let mut ids = Vec::new();
                        if scanned.first[i] != u32::MAX {
                            ids.push(scanned.first[i]);
                            ids.extend(scanned.extra.get(&(i as u32)).into_iter().flatten().copied());
                        }
                        let records = ids
                            .iter()
                            .map(|&id| {
                                gvcf::parse_record(record(id))
                                    .map_err(|e| Error::Invalid(format!("{contig}:{pos}: {e}")))
                            })
                            .collect::<Result<Vec<_>>>()?;
                        let (r, a, sequence);
                        let target = if is_sequence_key(key) {
                            let v = sequences
                                .binary_search_by_key(&key, |(k, _)| *k)
                                .map(|i| &sequences[i].1)
                                .map_err(|_| {
                                    Error::Invalid(format!("{contig}:{pos}: sequence target without alleles"))
                                })?;
                            sequence = v;
                            Target {
                                pos: sequence.pos,
                                ref_allele: &sequence.ref_allele,
                                alt: Some(&sequence.alt),
                                sequence: true,
                            }
                        } else {
                            let (_, _, ref_base, alt_base) = unpack_key(key);
                            (r, a) = ((ref_base as char).to_string(), (alt_base as char).to_string());
                            Target {
                                pos: pos as u64,
                                ref_allele: &r,
                                alt: (alt_base != crate::genotypes::ANY_ALT).then_some(a.as_str()),
                                sequence: false,
                            }
                        };
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
    targets_cache: Option<&Path>,
    threads: usize,
) -> Result<(GenotypeTable, Timings)> {
    let mut timings = Timings::default();
    let t = Instant::now();
    let (set, source) = crate::targets::targets(packs, reference_identity, targets_cache)?;
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
    let scanned = scan(gvcf, &set.keys, &ends, threads)?;
    timings.scan_s = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let calls = assess(&set.keys, &set.sequences, &scanned, reference, threads)?;
    timings.assess_s = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let table = GenotypeTable::new(
        gvcf,
        reference_identity,
        set.packs,
        set.keys,
        set.sequences,
        scanned,
        calls,
    );
    timings.table_s = t.elapsed().as_secs_f64();
    Ok((table, timings))
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
