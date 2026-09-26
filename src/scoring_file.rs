//! Streaming reader for scoring files: PGS Catalog harmonized files (format 2.0, `HmPOS_build=GRCh38`), and
//! custom files in the same layout with the author's own GRCh38 positions (see `DESIGN.md`, "Custom scores").

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

/// Where a scoring file comes from.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Origin {
    /// A PGS Catalog harmonized scoring file.
    #[default]
    PgsCatalog,
    /// Any other score: positions from `chr_name`/`chr_position` in GRCh38.
    Custom,
}

/// A custom score ID: letters, digits, `_`, `.`, `-` (at most 64, starting with a letter or digit), and not a
/// Catalog ID, so a custom score can never be mistaken for a Catalog one.
pub fn is_custom_id(value: &str) -> bool {
    (1..=64).contains(&value.len())
        && value.as_bytes()[0].is_ascii_alphanumeric()
        && value
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.' | b'-'))
        && !is_pgs_id(value)
}

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
    pub hm_infer_other_allele: Option<usize>,
    pub flags: [Option<usize>; 5],
    pub inclusion_criteria: Option<usize>,
    pub variant_description: Option<usize>,
    pub imputation_method: Option<usize>,
    /// The author's `chr_name` and `chr_position` (read by `term::finngen_identifier`).
    pub chr_name: Option<usize>,
    pub chr_position: Option<usize>,
    /// Accept `chr1` as well as `1` in the position columns (custom files).
    pub strip_chr_prefix: bool,
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
    fn new(names: Vec<String>, origin: Origin) -> Result<Columns> {
        let find = |name: &str| names.iter().position(|c| c == name);
        let (effect_allele, hm_chr, hm_pos) = match origin {
            Origin::PgsCatalog => {
                let (Some(e), Some(c), Some(p), Some(_)) =
                    (find("effect_allele"), find("hm_chr"), find("hm_pos"), find("hm_source"))
                else {
                    return invalid!("required harmonized columns are missing");
                };
                (e, c, p)
            }
            Origin::Custom => {
                let (Some(e), Some(c), Some(p)) = (find("effect_allele"), find("chr_name"), find("chr_position"))
                else {
                    return invalid!("custom scores need chr_name, chr_position and effect_allele columns");
                };
                (e, c, p)
            }
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
            strip_chr_prefix: origin == Origin::Custom,
            hm_match_chr: find("hm_match_chr"),
            hm_match_pos: find("hm_match_pos"),
            hm_infer_other_allele: find("hm_inferOtherAllele"),
            flags: FLAGS.map(find),
            inclusion_criteria: find("inclusion_criteria"),
            variant_description: find("variant_description"),
            imputation_method: find("imputation_method"),
            chr_name: find("chr_name"),
            chr_position: find("chr_position"),
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
    pub origin: Origin,
    /// `variants_number`; required for Catalog files, optional for custom ones.
    pub declared_terms: Option<u64>,
    reader: Source,
    line: Vec<u8>,
    terms: u64,
}

/// The file's bytes, gunzipped when they start with the gzip magic number, hashed as they are read.
enum Source {
    Gzip(Box<BufReader<MultiGzDecoder<BufReader<HashingReader<File>>>>>),
    Plain(BufReader<HashingReader<File>>),
}

impl Read for Source {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            Source::Gzip(r) => r.read(buf),
            Source::Plain(r) => r.read(buf),
        }
    }
}

impl BufRead for Source {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        match self {
            Source::Gzip(r) => r.fill_buf(),
            Source::Plain(r) => r.fill_buf(),
        }
    }
    fn consume(&mut self, n: usize) {
        match self {
            Source::Gzip(r) => r.consume(n),
            Source::Plain(r) => r.consume(n),
        }
    }
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
        let mut raw = BufReader::with_capacity(1 << 20, HashingReader::new(file));
        let gzip = raw.fill_buf().map_err(Error::io(path))?.starts_with(&[0x1f, 0x8b]);
        let mut reader = if gzip {
            Source::Gzip(Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(raw))))
        } else {
            Source::Plain(raw)
        };
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
            break Columns::new(names, origin_of(&metadata))
                .map_err(|e| Error::Invalid(format!("{}: {e}", path.display())))?;
        };
        let origin = origin_of(&metadata);
        let get = |key: &str| metadata.get(key).map(String::as_str);
        let pgs_id = get("pgs_id").unwrap_or_default().to_owned();
        match origin {
            Origin::PgsCatalog => {
                if get("format_version") != Some("2.0") || get("HmPOS_build") != Some("GRCh38") {
                    return invalid!("{}: requires harmonized GRCh38 format 2.0", path.display());
                }
                if !is_pgs_id(&pgs_id) {
                    return invalid!("{}: invalid pgs_id {pgs_id:?}", path.display());
                }
            }
            Origin::Custom => {
                if !is_custom_id(&pgs_id) {
                    return invalid!(
                        "{}: a custom score needs #pgs_id= with letters, digits, _ . - (not a PGSnnnnnn Catalog ID); got {pgs_id:?}",
                        path.display()
                    );
                }
                if !get("genome_build")
                    .is_some_and(|b| b.eq_ignore_ascii_case("GRCh38") || b.eq_ignore_ascii_case("hg38"))
                {
                    return invalid!(
                        "{}: a custom score needs #genome_build=GRCh38 (positions on another build must be lifted over first)",
                        path.display()
                    );
                }
            }
        }
        let declared_terms = match get("variants_number") {
            Some(n) => Some(
                Some(n)
                    .filter(|n| !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()))
                    .and_then(|n| n.parse::<u64>().ok())
                    .filter(|n| (1..=MAX_TERMS).contains(n))
                    .ok_or_else(|| Error::Invalid(format!("{}: invalid variants_number", path.display())))?,
            ),
            None if origin == Origin::Custom => None,
            None => return invalid!("{}: invalid variants_number", path.display()),
        };
        Ok(ScoringFile {
            path: path.to_owned(),
            header_lines,
            metadata,
            duplicate_keys,
            columns,
            pgs_id,
            origin,
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
        let mut hashing = match self.reader {
            Source::Gzip(r) => r.into_inner().into_inner().into_inner(),
            Source::Plain(r) => r.into_inner(),
        };
        std::io::copy(&mut hashing, &mut std::io::sink()).map_err(Error::io(&path))?;
        let bytes = hashing.bytes;
        Ok((hashing.finish(), bytes))
    }
}

/// Catalog when the header names a harmonization build or a Catalog ID; custom otherwise.
fn origin_of(metadata: &BTreeMap<String, String>) -> Origin {
    if metadata.contains_key("HmPOS_build") || metadata.get("pgs_id").is_some_and(|id| is_pgs_id(id)) {
        Origin::PgsCatalog
    } else {
        Origin::Custom
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
