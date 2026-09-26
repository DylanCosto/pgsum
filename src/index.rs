//! Tabix (`.tbi`) and CSI (`.csi`) indexes of a bgzipped VCF: which parts of the file can hold records
//! overlapping a region.
//!
//! Only what `extract` needs: for a region, the chunks (ranges of BGZF virtual positions) of every bin that
//! overlaps it, less those ending before the region's minimum offset (the tabix linear index, or for CSI the
//! `loffset` of the smallest bin holding the region start). Layouts follow the SAM/BAM specification
//! (section 5, "Indexing BAM"), the tabix format and the CSI v1 specification. Both are built by htslib
//! (`tabix -p vcf`, `bcftools index`), which places a record in the bins of its whole span, `INFO/END`
//! included, so a gVCF reference block is found from any position it covers.

use std::collections::HashMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;

use crate::{Error, Result, invalid};

/// A range of BGZF virtual positions, `[start, end)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Chunk {
    pub start: u64,
    pub end: u64,
}

#[derive(Default)]
struct Reference {
    /// Bin number → (CSI `loffset`, chunks).
    bins: HashMap<u32, (u64, Vec<Chunk>)>,
    /// Tabix linear index: the smallest offset of a record overlapping each 16 kb window.
    intervals: Vec<u64>,
}

pub struct Index {
    pub path: PathBuf,
    csi: bool,
    min_shift: u32,
    depth: u32,
    names: Vec<String>,
    references: Vec<Reference>,
}

/// Little-endian reader over the decompressed index.
struct Cursor<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl Cursor<'_> {
    fn take(&mut self, n: usize) -> Result<&[u8]> {
        let s = self
            .bytes
            .get(self.at..self.at + n)
            .ok_or_else(|| Error::Invalid("index ends early".into()))?;
        self.at += n;
        Ok(s)
    }
    fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().expect("4")))
    }
    fn count(&mut self) -> Result<usize> {
        usize::try_from(self.i32()?).map_err(|_| Error::Invalid("negative count in index".into()))
    }
    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().expect("4")))
    }
    fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().expect("8")))
    }
}

/// The tabix header fields shared by `.tbi` and the CSI auxiliary data: format, the sequence, begin and end
/// columns, meta character and skipped lines, then the NUL-separated sequence names.
fn read_names(c: &mut Cursor<'_>) -> Result<Vec<String>> {
    for _ in 0..6 {
        c.i32()?;
    }
    let len = c.count()?;
    let raw = c.take(len)?;
    raw.split(|&b| b == 0)
        .filter(|n| !n.is_empty())
        .map(|n| String::from_utf8(n.to_vec()).map_err(|_| Error::Invalid("sequence name is not UTF-8".into())))
        .collect()
}

/// Bins at every level overlapping the 0-based half-open `[beg, end)` (`reg2bins` generalised as in htslib).
fn region_bins(beg: u64, end: u64, min_shift: u32, depth: u32) -> Vec<u32> {
    let mut bins = Vec::new();
    if beg >= end {
        return bins;
    }
    let last = end - 1;
    for level in 0..=depth {
        let shift = min_shift + 3 * (depth - level);
        let first_bin = ((1u64 << (3 * level)) - 1) / 7;
        for b in first_bin + (beg >> shift)..=first_bin + (last >> shift) {
            bins.push(b as u32);
        }
    }
    bins
}

impl Index {
    /// The index next to `gvcf`: `<gvcf>.tbi`, else `<gvcf>.csi`.
    pub fn find(gvcf: &Path) -> Option<PathBuf> {
        ["tbi", "csi"]
            .iter()
            .map(|ext| PathBuf::from(format!("{}.{ext}", gvcf.display())))
            .find(|p| p.is_file())
    }

    pub fn open(path: &Path) -> Result<Index> {
        let mut bytes = Vec::new();
        MultiGzDecoder::new(std::fs::File::open(path).map_err(Error::io(path))?)
            .read_to_end(&mut bytes)
            .map_err(Error::io(path))?;
        Index::parse(&bytes)
            .map_err(|e| Error::Invalid(format!("{}: {e}", path.display())))
            .map(|mut i| {
                i.path = path.to_owned();
                i
            })
    }

    fn parse(bytes: &[u8]) -> Result<Index> {
        let mut c = Cursor { bytes, at: 0 };
        let magic = c.take(4)?;
        let (csi, min_shift, depth, names) = match magic {
            b"TBI\x01" => {
                let n_ref = c.count()?;
                let names = read_names(&mut c)?;
                if names.len() != n_ref {
                    return invalid!("{} sequence names for {n_ref} references", names.len());
                }
                (false, 14, 5, names)
            }
            b"CSI\x01" => {
                let min_shift = c.count()? as u32;
                let depth = c.count()? as u32;
                if !(1..=30).contains(&min_shift) || !(1..=10).contains(&depth) || min_shift + 3 * depth > 62 {
                    return invalid!("unsupported CSI dimensions (min_shift {min_shift}, depth {depth})");
                }
                let aux_len = c.count()?;
                let aux = c.take(aux_len)?;
                if aux_len < 28 {
                    return invalid!("CSI index without sequence names (not built for a VCF)");
                }
                let names = read_names(&mut Cursor { bytes: aux, at: 0 })?;
                (true, min_shift, depth, names)
            }
            _ => return invalid!("not a tabix or CSI index"),
        };
        let n_ref = if csi { c.count()? } else { names.len() };
        if n_ref != names.len() {
            return invalid!("{} sequence names for {n_ref} references", names.len());
        }
        let max_bin = (((1u64 << (3 * (depth + 1))) - 1) / 7) as u32;
        let mut references = Vec::with_capacity(n_ref);
        for _ in 0..n_ref {
            let mut r = Reference::default();
            for _ in 0..c.count()? {
                let bin = c.u32()?;
                let loffset = if csi { c.u64()? } else { 0 };
                let n_chunk = c.count()?;
                let mut chunks = Vec::with_capacity(n_chunk);
                for _ in 0..n_chunk {
                    let (start, end) = (c.u64()?, c.u64()?);
                    chunks.push(Chunk { start, end });
                }
                // The pseudo-bin past the last real bin holds statistics, not records.
                if bin < max_bin {
                    r.bins.insert(bin, (loffset, chunks));
                }
            }
            if !csi {
                for _ in 0..c.count()? {
                    r.intervals.push(c.u64()?);
                }
            }
            references.push(r);
        }
        Ok(Index {
            path: PathBuf::new(),
            csi,
            min_shift,
            depth,
            names,
            references,
        })
    }

    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Chunks that can hold records of reference `id` overlapping the 0-based half-open `[beg, end)`.
    pub fn chunks(&self, id: usize, beg: u64, end: u64) -> Vec<Chunk> {
        let Some(r) = self.references.get(id) else {
            return Vec::new();
        };
        let min_offset = if self.csi {
            // The loffset of the smallest existing bin holding `beg`.
            let mut bin = ((1u64 << (3 * self.depth)) - 1) / 7 + (beg >> self.min_shift);
            loop {
                if let Some((loffset, _)) = r.bins.get(&(bin as u32)) {
                    break *loffset;
                }
                if bin == 0 {
                    break 0;
                }
                bin = (bin - 1) >> 3;
            }
        } else {
            let window = (beg >> self.min_shift) as usize;
            r.intervals
                .get(window)
                .or(r.intervals.last().filter(|_| window >= r.intervals.len()))
                .copied()
                .unwrap_or(0)
        };
        region_bins(beg, end, self.min_shift, self.depth)
            .into_iter()
            .filter_map(|b| r.bins.get(&b))
            .flat_map(|(_, chunks)| chunks.iter().copied())
            .filter(|c| c.end > min_offset)
            .collect()
    }
}

/// Sort chunks and merge those that overlap or share a BGZF block, so each block is decompressed once.
pub fn merge(mut chunks: Vec<Chunk>) -> Vec<Chunk> {
    chunks.sort_unstable();
    let mut out: Vec<Chunk> = Vec::with_capacity(chunks.len());
    for c in chunks {
        match out.last_mut() {
            Some(last) if c.start >> 16 <= last.end >> 16 || c.start <= last.end => last.end = last.end.max(c.end),
            _ => out.push(c),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins_match_htslib() {
        // reg2bins(0, 1) for the BAI/tabix layout: bin 0 and the first bin of every level.
        assert_eq!(region_bins(0, 1, 14, 5), vec![0, 1, 9, 73, 585, 4681]);
        // One 16 kb window further along: the same parents, and leaf 4682.
        assert_eq!(region_bins(16_384, 16_385, 14, 5), vec![0, 1, 9, 73, 585, 4682]);
        // A region spanning two leaves.
        assert_eq!(region_bins(16_000, 17_000, 14, 5), vec![0, 1, 9, 73, 585, 4681, 4682]);
        assert!(region_bins(5, 5, 14, 5).is_empty());
    }

    #[test]
    fn merging_joins_chunks_in_one_block() {
        let c = |s: u64, e: u64| Chunk { start: s, end: e };
        let merged = merge(vec![
            c(5 << 16, (5 << 16) + 90),
            c(1 << 16, (1 << 16) + 10),
            c((5 << 16) + 200, 6 << 16),
        ]);
        assert_eq!(merged, vec![c(1 << 16, (1 << 16) + 10), c(5 << 16, 6 << 16)]);
    }
}
