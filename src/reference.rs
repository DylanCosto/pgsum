//! Memory-mapped reference FASTA with a `.fai` index.

use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};

use memmap2::Mmap;

use crate::{Error, Result, invalid};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Entry {
    length: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

pub struct Reference {
    pub path: PathBuf,
    pub fai_path: PathBuf,
    entries: HashMap<String, Entry>,
    map: Mmap,
}

impl Reference {
    pub fn open(path: &Path) -> Result<Reference> {
        let fai_path = PathBuf::from(format!("{}.fai", path.display()));
        if !fai_path.is_file() {
            return invalid!(
                "{}: no FASTA index; create it with `samtools faidx {}`",
                fai_path.display(),
                path.display()
            );
        }
        let text = std::fs::read_to_string(&fai_path).map_err(Error::io(&fai_path))?;
        let mut entries = HashMap::new();
        for line in text.lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            let numbers: Option<Vec<u64>> = fields
                .get(1..5)
                .map(|f| f.iter().filter_map(|v| v.parse().ok()).collect());
            let Some([length, offset, line_bases, line_width]) =
                numbers.as_deref().and_then(|n| <[u64; 4]>::try_from(n).ok())
            else {
                return invalid!("{}: malformed index line", fai_path.display());
            };
            if length < 1 || line_bases < 1 || line_width < line_bases {
                return invalid!("{}: invalid index dimensions for {}", fai_path.display(), fields[0]);
            }
            if entries
                .insert(
                    fields[0].to_owned(),
                    Entry {
                        length,
                        offset,
                        line_bases,
                        line_width,
                    },
                )
                .is_some()
            {
                return invalid!("{}: duplicate contig {}", fai_path.display(), fields[0]);
            }
        }
        if entries.is_empty() {
            return invalid!("{}: no contigs", fai_path.display());
        }
        // Reach `1`, `X` or `MT` (Ensembl, NCBI) under the names pgsum uses, unless the FASTA also has those.
        let aliases: Vec<(String, Entry)> = entries
            .iter()
            .filter_map(|(name, e)| {
                let canonical = crate::term::CONTIGS[crate::term::contig_code_of_name(name)? as usize - 1];
                (canonical != name && !entries.contains_key(canonical)).then(|| (canonical.to_owned(), *e))
            })
            .collect();
        entries.extend(aliases);
        let file = File::open(path).map_err(Error::io(path))?;
        // SAFETY: the reference is read-only input; changing the file while pgsum runs is unsupported.
        let map = unsafe { Mmap::map(&file) }.map_err(Error::io(path))?;
        Ok(Reference {
            path: path.to_owned(),
            fai_path,
            entries,
            map,
        })
    }

    pub fn contig_length(&self, contig: &str) -> Option<u64> {
        self.entries.get(contig).map(|e| e.length)
    }

    /// Upper-cased bases of the 0-based half-open interval `[start, end)`. Errors if the interval is outside
    /// the contig or a base is not `ACGTN`.
    pub fn fetch(&self, contig: &str, start: u64, end: u64) -> Result<String> {
        let Some(e) = self.entries.get(contig) else {
            return invalid!("{contig} is not in the reference");
        };
        if end < start || end > e.length {
            return invalid!("{contig}:{}-{end} is outside the reference", start + 1);
        }
        let mut out = String::with_capacity((end - start) as usize);
        let mut cursor = start;
        while cursor < end {
            let (line, in_line) = (cursor / e.line_bases, cursor % e.line_bases);
            let take = (end - cursor).min(e.line_bases - in_line);
            let from = (e.offset + line * e.line_width + in_line) as usize;
            let Some(chunk) = self.map.get(from..from + take as usize) else {
                return invalid!("{}: unexpected end of file", self.path.display());
            };
            for &b in chunk {
                let b = b.to_ascii_uppercase();
                if !matches!(b, b'A' | b'C' | b'G' | b'T' | b'N') {
                    return invalid!("{contig}:{}-{end} contains an unsupported base", start + 1);
                }
                out.push(b as char);
            }
            cursor += take;
        }
        Ok(out)
    }
}
