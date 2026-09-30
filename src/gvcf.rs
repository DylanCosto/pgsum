//! Single-sample gVCF text: header facts and record parsing.
//!
//! Every record is scanned for its interval (`CHROM`, `POS`, `REF`, `INFO/END`); only records overlapping a
//! target are fully parsed and validated.

use crate::genotype::Record;
use crate::{Error, Result, invalid};

/// The DeepVariant header line that defines the `RefCall` filter.
pub const REFCALL_DEFINITION: &str =
    "##FILTER=<ID=RefCall,Description=\"Genotyping model thinks this site is reference.\">";

#[derive(Clone, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct HeaderFacts {
    pub sample_id: String,
    /// The chosen sample's column (0-based) in a multi-sample VCF; 9 for a single-sample file.
    #[serde(skip)]
    pub sample_column: usize,
    /// Sample columns in the file.
    #[serde(skip)]
    pub samples: usize,
    /// `##contig` lengths by contig code, where the header gives them.
    #[serde(skip)]
    pub contig_lengths: Vec<(u8, u64)>,
    /// The sample to read from a multi-sample VCF (`--sample`), set before reading the header.
    #[serde(skip)]
    pub requested_sample: Option<String>,
    /// Read every sample (cohort extraction): the header's sample names are kept in `sample_names`.
    #[serde(skip)]
    pub all_samples: bool,
    #[serde(skip)]
    pub sample_names: Vec<String>,
    /// `##DeepVariant_version=`, when it is a plain `X.Y.Z` version.
    pub deepvariant_version: Option<String>,
    /// The header defines DeepVariant's `RefCall` filter.
    pub refcall_defined: bool,
}

impl HeaderFacts {
    /// Update from one header line; returns true once the `#CHROM` line has been read.
    pub fn read_line(&mut self, line: &str) -> Result<bool> {
        if let Some(version) = line.strip_prefix("##DeepVariant_version=") {
            let parts: Vec<&str> = version.split('.').collect();
            let plain = parts.len() == 3
                && parts
                    .iter()
                    .all(|p| !p.is_empty() && p.bytes().all(|b| b.is_ascii_digit()));
            self.deepvariant_version = plain.then(|| version.to_owned());
        }
        if let Some(fields) = line.strip_prefix("##contig=<").and_then(|l| l.strip_suffix('>')) {
            let value = |key: &str| {
                fields
                    .split(',')
                    .find_map(|f| f.strip_prefix(key).and_then(|v| v.strip_prefix('=')))
            };
            if let (Some(code), Some(Ok(length))) = (
                value("ID").and_then(crate::term::contig_code_of_name),
                value("length").map(str::parse::<u64>),
            ) {
                self.contig_lengths.push((code, length));
            }
        }
        // BCF headers commonly attach IDX and can reorder attributes. Match the known
        // DeepVariant definition semantically instead of comparing the whole header line.
        if let Some(fields) = line.strip_prefix("##FILTER=<").and_then(|s| s.strip_suffix('>')) {
            let fields: Vec<_> = fields.split(',').collect();
            if fields.contains(&"ID=RefCall") {
                self.refcall_defined = fields.iter().filter(|f| f.starts_with("Description=")).count() == 1
                    && fields.contains(&"Description=\"Genotyping model thinks this site is reference.\"");
            }
        }
        if line.starts_with("#CHROM\t") {
            let fields: Vec<&str> = line.split('\t').collect();
            let samples = fields.get(9..).unwrap_or_default();
            if self.all_samples {
                if samples.is_empty() {
                    return invalid!("the VCF has no sample columns");
                }
                let unique: std::collections::HashSet<_> = samples.iter().collect();
                if unique.len() != samples.len() {
                    return invalid!("the VCF names a sample more than once");
                }
                self.sample_names = samples.iter().map(|s| (*s).to_owned()).collect();
                self.sample_id = format!("{} samples", samples.len());
                self.sample_column = 9;
                self.samples = samples.len();
                return Ok(true);
            }
            let column = match (&self.requested_sample, samples) {
                (_, []) => return invalid!("the VCF has no sample columns"),
                (Some(name), _) => match samples.iter().position(|s| s == name) {
                    Some(i) if samples.iter().filter(|s| *s == name).count() == 1 => 9 + i,
                    Some(_) => return invalid!("sample {name} appears more than once"),
                    None => return invalid!("sample {name} is not in the VCF ({} samples)", samples.len()),
                },
                (None, [_]) => 9,
                (None, _) => {
                    return invalid!(
                        "the VCF has {} samples; choose one with --sample (e.g. --sample {})",
                        samples.len(),
                        samples[0]
                    );
                }
            };
            if fields[column].is_empty() {
                return invalid!("empty sample name");
            }
            self.sample_id = fields[column].to_owned();
            self.sample_column = column;
            self.samples = samples.len();
            return Ok(true);
        }
        Ok(false)
    }
}

/// The interval of a record: `CHROM`, `POS` and the inclusive end (`INFO/END`, else `POS + len(REF) - 1`).
pub struct Interval<'a> {
    pub chrom: &'a str,
    pub pos: u64,
    pub end: u64,
    pub alt: &'a str,
}

/// Whether every ALT allele of a record is a structural-variant symbol (`<DEL>`, `<INV>`, …) or a breakend;
/// the gVCF placeholders `<*>` and `<NON_REF>` are not.
pub fn is_structural(alt: &str) -> bool {
    !alt.is_empty()
        && alt
            .split(',')
            .all(|a| (a.starts_with('<') && a != "<*>" && a != "<NON_REF>") || a.contains('[') || a.contains(']'))
}

fn parse_u64(text: &str) -> Option<u64> {
    (!text.is_empty() && text.bytes().all(|b| b.is_ascii_digit()))
        .then(|| text.parse().ok())
        .flatten()
}

/// `INFO/END` when present (the last one, if repeated).
fn info_end(info: &str) -> Result<Option<u64>> {
    let mut end = None;
    if info == "." {
        return Ok(None);
    }
    for item in info.split(';') {
        if item == "END" {
            return invalid!("INFO/END has no value");
        }
        if let Some(value) = item.strip_prefix("END=") {
            end = Some(parse_u64(value).ok_or_else(|| Error::Invalid(format!("invalid INFO/END {value:?}")))?);
        }
    }
    Ok(end)
}

pub fn interval(line: &str) -> Result<Interval<'_>> {
    let mut fields = line.splitn(9, '\t');
    let mut next = || {
        fields
            .next()
            .ok_or_else(|| Error::Invalid("record has too few fields".into()))
    };
    let chrom = next()?;
    let pos_text = next()?;
    let _id = next()?;
    let ref_allele = next()?;
    let alt = next()?;
    let _qual = next()?;
    let _filter = next()?;
    let info = next()?;
    let pos = parse_u64(pos_text)
        .filter(|&p| p >= 1)
        .ok_or_else(|| Error::Invalid(format!("invalid POS {pos_text:?}")))?;
    let end = match info_end(info)? {
        Some(end) => end,
        None => pos + ref_allele.len() as u64 - 1,
    };
    if end < pos {
        return invalid!("record end {end} is before POS {pos}");
    }
    Ok(Interval { chrom, pos, end, alt })
}

fn valid_gt(gt: &str) -> bool {
    !gt.is_empty()
        && gt
            .split(['/', '|'])
            .all(|a| a == "." || !a.is_empty() && a.bytes().all(|b| b.is_ascii_digit()))
}

/// Fully parse and validate a record for the genotype rules.
pub fn parse_record(line: &str) -> Result<Record> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() != 10 {
        return invalid!("record has {} fields, expected 10", fields.len());
    }
    let span = interval(line)?;
    let keys: Vec<&str> = fields[8].split(':').collect();
    let values: Vec<&str> = fields[9].split(':').collect();
    let unique: std::collections::HashSet<_> = keys.iter().collect();
    if unique.len() != keys.len() || values.len() > keys.len() {
        return invalid!("invalid FORMAT or sample fields");
    }
    let sample = |key: &str| {
        keys.iter()
            .position(|k| *k == key)
            .and_then(|i| values.get(i))
            .map(|v| (*v).to_owned())
    };
    let raw_gt = sample("GT").unwrap_or_else(|| ".".to_owned());
    if !valid_gt(&raw_gt) {
        return invalid!("invalid GT {raw_gt:?}");
    }
    let alts: Vec<String> = fields[4].split(',').map(str::to_owned).collect();
    let mut gt = Vec::new();
    for allele in raw_gt.split(['/', '|']) {
        if allele == "." {
            gt.push(None);
        } else {
            let index: u32 = allele
                .parse()
                .map_err(|_| Error::Invalid(format!("invalid GT {raw_gt:?}")))?;
            if index as usize > alts.len() {
                return invalid!("GT {raw_gt:?} exceeds the allele list");
            }
            gt.push(Some(index));
        }
    }
    Ok(Record {
        pos: span.pos,
        end: span.end,
        ref_allele: fields[3].to_owned(),
        alts,
        filters: fields[6].split(';').map(str::to_owned).collect(),
        gt_phased: raw_gt.contains('|'),
        gt,
        phase_set: sample("PS"),
        ft: sample("FT"),
        dp: sample("DP"),
        min_dp: sample("MIN_DP"),
        gq: sample("GQ"),
        ds: sample("DS"),
        gp: sample("GP"),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn refcall_definition_allows_bcf_indices_without_trusting_other_descriptions() {
        let mut header = HeaderFacts::default();
        header
            .read_line("##FILTER=<IDX=7,Description=\"Genotyping model thinks this site is reference.\",ID=RefCall>")
            .unwrap();
        assert!(header.refcall_defined);
        header
            .read_line("##FILTER=<ID=RefCall,Description=\"An unrelated filter\",IDX=7>")
            .unwrap();
        assert!(!header.refcall_defined);
    }

    #[test]
    fn structural_alleles() {
        for alt in ["<DEL>", "<INS:ME:ALU>", "<DUP>,<INV>", "G]17:198982]", "]13:123456]T"] {
            assert!(is_structural(alt), "{alt}");
        }
        for alt in ["T", "<*>", "<NON_REF>", "C,<*>", "<DEL>,T", "."] {
            assert!(!is_structural(alt), "{alt}");
        }
    }

    #[test]
    fn block_and_variant_records() {
        let block = "chr1\t10002\t.\tA\t<*>\t0\t.\tEND=10010\tGT:GQ:MIN_DP:PL\t0/0:24:28:0,24,719";
        let r = parse_record(block).unwrap();
        assert!(r.is_reference_block());
        assert_eq!((r.pos, r.end, r.gt.clone()), (10002, 10010, vec![Some(0), Some(0)]));
        assert_eq!(
            (r.min_dp.as_deref(), r.dp.as_deref(), r.gq.as_deref()),
            (Some("28"), None, Some("24"))
        );

        let variant = "chr1\t10011\t.\tC\tCCT,<*>\t0.3\tNoCall\t.\tGT:GQ:DP:PS\t1|0:11:33:9876";
        let r = parse_record(variant).unwrap();
        assert!(!r.is_reference_block());
        assert_eq!(
            (r.end, r.gt_phased, r.phase_set.as_deref()),
            (10011, true, Some("9876"))
        );
        assert_eq!(r.filters, ["NoCall"]);
    }

    #[test]
    fn deletion_end_and_missing_gt() {
        let r = parse_record("chr1\t100\t.\tACGT\tA\t10\tPASS\t.\tDP\t20").unwrap();
        assert_eq!((r.end, r.gt.clone()), (103, vec![None]));
    }

    #[test]
    fn short_sample_field_is_allowed() {
        let r = parse_record("chr1\t100\t.\tA\tG\t10\tPASS\t.\tGT:GQ:DP\t0/1:30").unwrap();
        assert_eq!(r.dp, None);
    }

    #[test]
    fn rejects_invalid_records() {
        for bad in [
            "chr1\t100\t.\tA\tG\t10\tPASS\t.\tGT\t0/2",
            "chr1\t100\t.\tA\tG\t10\tPASS\t.\tGT\t0-1",
            "chr1\t100\t.\tA\tG\t10\tPASS\t.\tGT:GT\t0/1:0/1",
            "chr1\t100\t.\tA\tG\t10\tPASS\t.\tGT\t0/1:30",
            "chr1\t100\t.\tA\t<*>\t0\t.\tEND=99\tGT\t0/0",
            "chr1\t0\t.\tA\tG\t10\tPASS\t.\tGT\t0/1",
            "chr1\t100\t.\tA\tG\t10\tPASS\t.\tGT",
        ] {
            assert!(parse_record(bad).is_err(), "{bad}");
        }
    }
}
