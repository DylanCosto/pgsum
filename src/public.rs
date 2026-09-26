//! A public variant set used to orient indels and multi-base terms (e.g. 1000 Genomes non-SNV PASS records).
//!
//! Format: tab-separated `CHROM POS REF ALT`, plain or gzipped, `#` lines ignored. Only biallelic records
//! are kept (an ALT with a comma is skipped). A term is oriented by it when exactly one record at the term's
//! harmonized position has exactly the term's two alleles; its REF is then the reference allele.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::digest::file_sha256;
use crate::pack::SourceFile;
use crate::{Error, Result, invalid};

pub struct PublicVariants {
    pub identity: SourceFile,
    pairs: HashMap<(u8, u32), Vec<(String, String)>>,
}

impl PublicVariants {
    pub fn open(path: &Path) -> Result<PublicVariants> {
        let identity = SourceFile {
            name: path
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default(),
            sha256: file_sha256(path)?,
            bytes: std::fs::metadata(path).map_err(Error::io(path))?.len(),
        };
        let file = std::fs::File::open(path).map_err(Error::io(path))?;
        let mut raw = BufReader::new(file);
        let gzip = raw.fill_buf().map_err(Error::io(path))?.starts_with(&[0x1f, 0x8b]);
        let reader: Box<dyn Read> = if gzip {
            Box::new(MultiGzDecoder::new(raw))
        } else {
            Box::new(raw)
        };
        let mut pairs: HashMap<(u8, u32), Vec<(String, String)>> = HashMap::new();
        for (n, line) in BufReader::with_capacity(1 << 20, reader).lines().enumerate() {
            let line = line.map_err(Error::io(path))?;
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            let (Some(&chrom), Some(pos), Some(&r), Some(&a)) = (
                f.first(),
                f.get(1).and_then(|p| p.parse::<u32>().ok()),
                f.get(2),
                f.get(3),
            ) else {
                return invalid!("{}: line {} is not CHROM POS REF ALT", path.display(), n + 1);
            };
            let Some(code) = crate::term::contig_code_of_name(chrom) else {
                continue;
            };
            if a.contains(',') {
                continue;
            }
            pairs.entry((code, pos)).or_default().push((r.to_owned(), a.to_owned()));
        }
        Ok(PublicVariants { identity, pairs })
    }

    /// The `(REF, ALT)` of the one record at `(contig, pos)` whose alleles are exactly `{a, b}`, if there is
    /// exactly one.
    pub fn exact_pair(&self, contig: u8, pos: u32, a: &str, b: &str) -> Option<(&str, &str)> {
        let mut matches = self
            .pairs
            .get(&(contig, pos))?
            .iter()
            .filter(|(r, alt)| (r == a && alt == b) || (r == b && alt == a));
        let (r, alt) = matches.next()?;
        matches.next().is_none().then_some((r.as_str(), alt.as_str()))
    }

    pub fn len(&self) -> usize {
        self.pairs.values().map(Vec::len).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }
}
