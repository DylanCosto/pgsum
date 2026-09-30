//! Input preflight and provenance for the single-command workflow.
use crate::{Error, Result};
use std::io::BufRead;
use std::path::{Path, PathBuf};

/// Read only the VCF/BCF header and check sample selection and stated contig lengths.
/// This cannot establish assembly identity when the source omits assembly evidence.
pub fn preflight(input: &Path, reference: &Path, sample: Option<&str>) -> Result<serde_json::Value> {
    let reference = crate::reference::Reference::open(reference)?;
    // Ensure a previous metadata-based digest cache cannot hide changed reference bytes.
    let digest = crate::digest::file_sha256(&reference.path)?;
    crate::digest::remember_file_sha256(&reference.path, &digest);
    let mut raw = std::io::BufReader::new(std::fs::File::open(input).map_err(Error::io(input))?);
    let compressed = raw.fill_buf().map_err(Error::io(input))?.starts_with(&[31, 139]);
    let mut reader: Box<dyn BufRead> = if compressed {
        Box::new(std::io::BufReader::new(flate2::read::MultiGzDecoder::new(raw)))
    } else {
        Box::new(raw)
    };
    let bcf = reader.fill_buf().map_err(Error::io(input))?.starts_with(b"BCF");
    if bcf {
        reader = Box::new(crate::bcf::TextReader::new(reader).map_err(Error::io(input))?);
    }
    let mut header = crate::gvcf::HeaderFacts {
        requested_sample: sample.map(str::to_owned),
        ..Default::default()
    };
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line).map_err(Error::io(input))? == 0 {
            return invalid!("{}: no VCF sample header", input.display());
        }
        let line = line.trim_end_matches(['\r', '\n']);
        if !line.starts_with('#') {
            return invalid!("{}: record before VCF sample header", input.display());
        }
        if header.read_line(line)? {
            break;
        }
    }
    crate::extract::check_assembly(input, &header, &reference)?;
    let checked_contigs = header
        .contig_lengths
        .iter()
        .filter(|(code, _)| {
            reference
                .contig_length(crate::term::CONTIGS[*code as usize - 1])
                .is_some()
        })
        .count();
    Ok(
        serde_json::json!({ "input": input, "format": if bcf { "BCF" } else { "VCF" },
        "compressed": compressed, "sample": header.sample_id, "samples_in_file": header.samples,
        "matching_contig_lengths": checked_contigs,
        "assembly_evidence": if checked_contigs == 0 { "not_provided" } else { "matching_contig_lengths" },
        "note": "Header checks do not validate all records or prove assembly identity." }),
    )
}

pub fn pack_provenance(paths: &[PathBuf]) -> Result<Vec<serde_json::Value>> {
    paths
        .iter()
        .map(|path| {
            let h = crate::pack::Pack::open_header(path)?;
            Ok(serde_json::json!({ "path": path, "pgs_id": h.pgs_id,
            "sha256": crate::digest::file_sha256(path)?, "source": h.source,
            "catalog_metadata_sha256": h.catalog_metadata_sha256,
            "records_sha256": h.records_sha256, "reference": h.reference,
            "compile_rules": h.compile_rules }))
        })
        .collect()
}

pub fn write_json(path: &Path, value: &serde_json::Value) -> Result<()> {
    let temp = path.with_extension("json.tmp");
    std::fs::write(
        &temp,
        serde_json::to_vec_pretty(value).map_err(|e| Error::Invalid(e.to_string()))?,
    )
    .map_err(Error::io(&temp))?;
    std::fs::File::open(&temp)
        .and_then(|f| f.sync_all())
        .map_err(Error::io(&temp))?;
    std::fs::rename(&temp, path).map_err(Error::io(path))
}

/// Include only outputs produced by this invocation, never stale files from an earlier selection.
pub fn finish(
    out: &Path,
    packs: &[serde_json::Value],
    settings: serde_json::Value,
    preflight: serde_json::Value,
    terms: bool,
    bundle: bool,
) -> Result<()> {
    let table = crate::genotypes::GenotypeTable::open(&out.join("genotypes.pgsg"))?;
    for pack in packs {
        let path = pack["path"]
            .as_str()
            .ok_or_else(|| Error::Invalid("missing pack path".into()))?;
        if crate::digest::file_sha256(Path::new(path))? != pack["sha256"]
            || !table
                .header
                .packs
                .iter()
                .any(|p| p.pgs_id == pack["pgs_id"] && p.records_sha256 == pack["records_sha256"])
        {
            return invalid!("pack changed during execution: {path}; no completion manifest was published");
        }
    }
    let mut names = vec!["genotypes.pgsg".to_owned(), "scores.tsv".to_owned()];
    if bundle {
        names.push("results.jsonl.zst".into());
    }
    for pack in packs {
        let id = pack["pgs_id"]
            .as_str()
            .ok_or_else(|| Error::Invalid("missing pack ID".into()))?;
        if !bundle {
            names.push(format!("{id}.score.json"));
        }
        if terms {
            names.push(format!("{id}.terms.tsv"));
        }
    }
    let mut outputs = std::collections::BTreeMap::new();
    for name in names {
        outputs.insert(name.clone(), crate::digest::file_sha256(&out.join(name))?);
    }
    let executable = std::env::current_exe().map_err(|e| Error::Invalid(e.to_string()))?;
    let manifest = serde_json::json!({ "schema": "pgsum-execution-v1", "complete": true,
        "pgsum_version": env!("CARGO_PKG_VERSION"), "executable_sha256": crate::digest::file_sha256(&executable)?,
        "settings": settings, "preflight": preflight, "input": table.header.gvcf,
        "reference": table.header.reference, "policy": table.header.policy, "sample": table.header.sample,
        "packs": packs, "outputs": outputs });
    write_json(&out.join("execution-manifest.json"), &manifest)
}
