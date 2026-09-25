//! Streaming reader for PGS Catalog harmonized scoring files (format 2.0, `HmPOS_build=GRCh38`).

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::bufread::MultiGzDecoder;

use crate::digest::HashingReader;
use crate::{Error, Result, invalid};

pub const MAX_LINE: usize = 65536;
pub const MAX_HEADER_LINES: usize = 2000;
pub const MAX_COLUMNS: usize = 128;
/// Larger than any current Catalog score (13.1M variants in September 2026).
pub const MAX_TERMS: u64 = 50_000_000;
/// Header keys that must appear at most once: the ones pgsum reads, and the genome builds.
pub const STRICT_HEADER_KEYS: [&str; 6] = [
    "format_version",
    "pgs_id",
    "variants_number",
    "weight_type",
    "genome_build",
    "HmPOS_build",
];

/// The columns the rules read, by index into a row.
#[derive(Clone, Debug, Default)]
pub struct Columns {
    pub names: Vec<String>,
    pub effect_allele: usize,
    pub other_allele: Option<usize>,
    pub effect_weight: Option<usize>,
    /// `dosage_{0,1,2}_weight`; a column may be absent even when others are present.
    pub dosage_weights: [Option<usize>; 3],
    pub hm_chr: usize,
    pub hm_pos: usize,
    pub hm_match_chr: Option<usize>,
    pub hm_match_pos: Option<usize>,
    pub flags: [Option<usize>; 5],
    pub inclusion_criteria: Option<usize>,
    pub variant_description: Option<usize>,
    pub imputation_method: Option<usize>,
}

/// Model flag columns, in the order of `Columns::flags`.
pub const FLAGS: [&str; 5] = [
    "is_haplotype",
    "is_diplotype",
    "is_interaction",
    "is_dominant",
    "is_recessive",
];

impl Columns {
    fn new(names: Vec<String>) -> Result<Columns> {
        let find = |name: &str| names.iter().position(|c| c == name);
        let (Some(effect_allele), Some(hm_chr), Some(hm_pos), Some(_)) =
            (find("effect_allele"), find("hm_chr"), find("hm_pos"), find("hm_source"))
        else {
            return invalid!("required harmonized columns are missing");
        };
        let dosage_weights = [
            find("dosage_0_weight"),
            find("dosage_1_weight"),
            find("dosage_2_weight"),
        ];
        let effect_weight = find("effect_weight");
        if effect_weight.is_none() && dosage_weights.iter().any(Option::is_none) {
            return invalid!("neither effect_weight nor all dosage weight columns are present");
        }
        Ok(Columns {
            effect_allele,
            other_allele: find("other_allele"),
            effect_weight,
            dosage_weights,
            hm_chr,
            hm_pos,
            hm_match_chr: find("hm_match_chr"),
            hm_match_pos: find("hm_match_pos"),
            flags: FLAGS.map(find),
            inclusion_criteria: find("inclusion_criteria"),
            variant_description: find("variant_description"),
            imputation_method: find("imputation_method"),
            names,
        })
    }
}

/// An open scoring file positioned at its first term.
pub struct ScoringFile {
    pub path: PathBuf,
    pub header_lines: Vec<String>,
    /// `#key=value` header lines; a repeated key keeps its first value.
    pub metadata: BTreeMap<String, String>,
    /// Descriptive header keys that appear more than once (all their lines are in `header_lines`).
    pub duplicate_keys: Vec<String>,
    pub columns: Columns,
    pub pgs_id: String,
    pub declared_terms: u64,
    reader: BufReader<MultiGzDecoder<BufReader<HashingReader<File>>>>,
    line: Vec<u8>,
    terms: u64,
}

/// One term row: the line without its line ending, and the byte ranges of its fields.
pub struct Row<'a> {
    pub line: &'a str,
    fields: &'a [(usize, usize)],
}

impl<'a> Row<'a> {
    pub fn get(&self, index: usize) -> &'a str {
        let (a, b) = self.fields[index];
        &self.line[a..b]
    }

    pub fn opt(&self, index: Option<usize>) -> &'a str {
        index.map_or("", |i| self.get(i))
    }
}

impl ScoringFile {
    pub fn open(path: &Path) -> Result<ScoringFile> {
        let file = File::open(path).map_err(Error::io(path))?;
        let decoder = MultiGzDecoder::new(BufReader::with_capacity(1 << 20, HashingReader::new(file)));
        let mut reader = BufReader::with_capacity(1 << 20, decoder);
        let mut header_lines = Vec::new();
        let mut metadata = BTreeMap::new();
        let mut duplicate_keys = Vec::new();
        let mut line = Vec::new();
        let columns = loop {
            let text = read_line(&mut reader, &mut line, path)?
                .ok_or_else(|| Error::Invalid(format!("{}: no terms", path.display())))?;
            if text.starts_with('#') {
                header_lines.push(text.to_owned());
                if header_lines.len() > MAX_HEADER_LINES {
                    return invalid!("{}: header exceeds {MAX_HEADER_LINES} lines", path.display());
                }
                if !text.starts_with("##")
                    && let Some((key, value)) = text.strip_prefix('#').and_then(|t| t.split_once('='))
                {
                    if metadata.contains_key(key) {
                        // A repeated descriptive key (e.g. two citations) keeps its first value here; every
                        // line stays in `header_lines`. Keys the rules read must be unambiguous.
                        if STRICT_HEADER_KEYS.contains(&key) {
                            return invalid!("{}: duplicate header key {key}", path.display());
                        }
                        if !duplicate_keys.iter().any(|k| k == key) {
                            duplicate_keys.push(key.to_owned());
                        }
                    } else {
                        metadata.insert(key.to_owned(), value.to_owned());
                    }
                }
                continue;
            }
            if text.is_empty() {
                return invalid!("{}: unexpected empty header line", path.display());
            }
            let names: Vec<String> = text.split('\t').map(str::to_owned).collect();
            let unique: std::collections::HashSet<_> = names.iter().collect();
            if names.len() > MAX_COLUMNS
                || unique.len() != names.len()
                || !names
                    .iter()
                    .all(|c| !c.is_empty() && c.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_'))
            {
                return invalid!("{}: invalid columns", path.display());
            }
            break Columns::new(names).map_err(|e| Error::Invalid(format!("{}: {e}", path.display())))?;
        };
        if metadata.get("format_version").map(String::as_str) != Some("2.0")
            || metadata.get("HmPOS_build").map(String::as_str) != Some("GRCh38")
        {
            return invalid!("{}: requires harmonized GRCh38 format 2.0", path.display());
        }
        let pgs_id = metadata.get("pgs_id").cloned().unwrap_or_default();
        if !is_pgs_id(&pgs_id) {
            return invalid!("{}: invalid pgs_id {pgs_id:?}", path.display());
        }
        let declared_terms = metadata
            .get("variants_number")
            .filter(|n| !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()))
            .and_then(|n| n.parse::<u64>().ok())
            .filter(|n| (1..=MAX_TERMS).contains(n))
            .ok_or_else(|| Error::Invalid(format!("{}: invalid variants_number", path.display())))?;
        Ok(ScoringFile {
            path: path.to_owned(),
            header_lines,
            metadata,
            duplicate_keys,
            columns,
            pgs_id,
            declared_terms,
            reader,
            line,
            terms: 0,
        })
    }

    /// The next term, or `None` at the end of the file.
    pub fn next_row<'a>(&'a mut self, fields: &'a mut Vec<(usize, usize)>) -> Result<Option<Row<'a>>> {
        let Some(text) = read_line(&mut self.reader, &mut self.line, &self.path)? else {
            if self.terms == 0 {
                return invalid!("{}: no terms", self.path.display());
            }
            return Ok(None);
        };
        if text.is_empty() || text.starts_with('#') {
            return invalid!("{}: unexpected line inside terms", self.path.display());
        }
        fields.clear();
        let mut start = 0;
        for (i, b) in text.bytes().enumerate() {
            if b == b'\t' {
                fields.push((start, i));
                start = i + 1;
            }
        }
        fields.push((start, text.len()));
        if fields.len() != self.columns.names.len() {
            return invalid!(
                "{}: term {} has {} fields, expected {}",
                self.path.display(),
                self.terms + 1,
                fields.len(),
                self.columns.names.len()
            );
        }
        self.terms += 1;
        if self.terms > MAX_TERMS {
            return invalid!("{}: more than {MAX_TERMS} terms", self.path.display());
        }
        Ok(Some(Row { line: text, fields }))
    }

    /// Read to the end of the compressed file and return its SHA-256 and size.
    pub fn finish(self) -> Result<(String, u64)> {
        let path = self.path;
        let mut hashing = self.reader.into_inner().into_inner().into_inner();
        std::io::copy(&mut hashing, &mut std::io::sink()).map_err(Error::io(&path))?;
        let bytes = hashing.bytes;
        Ok((hashing.finish(), bytes))
    }
}

pub fn is_pgs_id(value: &str) -> bool {
    value.len() == 9 && value.starts_with("PGS") && value[3..].bytes().all(|b| b.is_ascii_digit())
}

/// One line without its trailing `\r`/`\n` characters, or `None` at EOF.
fn read_line<'a>(reader: &mut impl BufRead, line: &'a mut Vec<u8>, path: &Path) -> Result<Option<&'a str>> {
    line.clear();
    let n = reader
        .take(MAX_LINE as u64 + 1)
        .read_until(b'\n', line)
        .map_err(Error::io(path))?;
    if n == 0 {
        return Ok(None);
    }
    if n > MAX_LINE {
        return invalid!("{}: line exceeds {MAX_LINE} bytes", path.display());
    }
    while matches!(line.last(), Some(b'\n' | b'\r')) {
        line.pop();
    }
    std::str::from_utf8(line)
        .map(Some)
        .map_err(|_| Error::Invalid(format!("{}: invalid UTF-8", path.display())))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn file(name: &str, header: &[&str]) -> PathBuf {
        let path = std::env::temp_dir().join(format!("pgsum-{}-{name}.txt.gz", std::process::id()));
        let mut gz = flate2::write::GzEncoder::new(File::create(&path).unwrap(), flate2::Compression::fast());
        for line in header {
            writeln!(gz, "{line}").unwrap();
        }
        writeln!(
            gz,
            "effect_allele\teffect_weight\thm_source\thm_chr\thm_pos\nA\t0.1\tENSEMBL\t1\t100"
        )
        .unwrap();
        gz.finish().unwrap();
        path
    }

    const BASE: [&str; 5] = [
        "#format_version=2.0",
        "#pgs_id=PGS000001",
        "#variants_number=1",
        "#weight_type=beta",
        "#HmPOS_build=GRCh38",
    ];

    #[test]
    fn repeated_descriptive_key_is_kept() {
        let mut header = BASE.to_vec();
        header.extend(["#citation=First et al.", "#citation=Second et al."]);
        let path = file("citation", &header);
        let f = ScoringFile::open(&path).unwrap();
        assert_eq!(f.metadata["citation"], "First et al.");
        assert_eq!(f.duplicate_keys, ["citation"]);
        assert_eq!(f.header_lines.iter().filter(|l| l.starts_with("#citation=")).count(), 2);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn repeated_strict_key_is_rejected() {
        for key in STRICT_HEADER_KEYS {
            let mut header = BASE.to_vec();
            let line = format!("#{key}=other");
            header.push(&line);
            if !BASE.iter().any(|b| b.starts_with(&format!("#{key}="))) {
                header.push(&line);
            }
            let path = file(key, &header);
            let err = ScoringFile::open(&path).err().expect("rejected").to_string();
            assert!(err.contains("duplicate header key"), "{key}: {err}");
            std::fs::remove_file(path).unwrap();
        }
    }
}
