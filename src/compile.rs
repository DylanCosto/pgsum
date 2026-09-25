//! Compile a harmonized scoring file and its Catalog metadata into a pack.

use std::path::{Path, PathBuf};

use sha2::{Digest, Sha256};

use crate::digest::file_sha256;
use crate::orient::orient;
use crate::pack::{
    self, Columns, Counts, Header, Inventory, ReferenceIdentity, SourceFile, TermRecord, Weight, WeightType,
};
use crate::reference::Reference;
use crate::scoring_file::ScoringFile;
use crate::term::{Reason, describe};
use crate::{Error, Result, invalid};

/// The Catalog metadata file expected next to a scoring file: `PGS000001_hmPOS_GRCh38.txt.gz` →
/// `PGS000001.metadata.json`.
pub fn metadata_path(scoring_file: &Path) -> Option<PathBuf> {
    let name = scoring_file.file_name()?.to_str()?;
    let id = name.split('_').next()?;
    Some(scoring_file.with_file_name(format!("{id}.metadata.json")))
}

pub fn reference_identity(reference: &Reference) -> Result<ReferenceIdentity> {
    Ok(ReferenceIdentity {
        fasta_name: reference
            .path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default(),
        fasta_sha256: file_sha256(&reference.path)?,
        fai_sha256: file_sha256(&reference.fai_path)?,
    })
}

/// Compile one scoring file into `<out_dir>/<pgs_id>.pgsp` and return the pack header.
pub fn compile_file(
    scoring_path: &Path,
    metadata_path: &Path,
    reference: &Reference,
    reference_identity: &ReferenceIdentity,
    out_dir: &Path,
) -> Result<Header> {
    let metadata_bytes = std::fs::read(metadata_path).map_err(Error::io(metadata_path))?;
    let metadata: serde_json::Value = serde_json::from_slice(&metadata_bytes)
        .map_err(|e| Error::Invalid(format!("{}: {e}", metadata_path.display())))?;
    let mut file = ScoringFile::open(scoring_path)?;
    if metadata.get("id").and_then(|v| v.as_str()) != Some(file.pgs_id.as_str()) {
        return invalid!("{}: metadata is not for {}", metadata_path.display(), file.pgs_id);
    }
    let header_weight_type = file.metadata.get("weight_type").cloned();
    let columns = file.columns.clone();
    let mut terms = Columns::default();
    let mut line_hashes: Vec<([u8; 16], u32)> = Vec::with_capacity(file.declared_terms as usize);
    let mut fields = Vec::new();
    let mut centred_terms = 0;
    let mut ordinal: u32 = 0;
    while let Some(row) = file.next_row(&mut fields)? {
        let d = describe(&row, &columns);
        if let Some(c) = &d.centred {
            centred_terms += 1;
            if Some(&c.weight_type) != header_weight_type.as_ref() {
                return invalid!(
                    "{}: term {} weight type differs from the header",
                    scoring_path.display(),
                    ordinal + 1
                );
            }
        }
        let orientation = orient(&d, reference)
            .map_err(|e| Error::Invalid(format!("{}: term {}: {e}", scoring_path.display(), ordinal + 1)))?;
        let texts: Vec<&str> = if d.weights.len() == 3 {
            columns.dosage_weights.iter().map(|i| row.opt(*i)).collect()
        } else {
            vec![row.opt(columns.effect_weight)]
        };
        let weights = d
            .weights
            .iter()
            .zip(texts)
            .map(|(w, t)| Weight::from_decimal(w.as_ref(), t))
            .collect();
        let hash: [u8; 32] = Sha256::digest(row.line.as_bytes()).into();
        line_hashes.push((hash[..16].try_into().expect("16 bytes"), ordinal));
        terms.push(&TermRecord {
            contig: d.contig,
            pos: d.pos,
            model: d.model,
            allele_kind: d.allele_kind,
            palindromic: d.palindromic,
            orientation,
            reasons: d.reasons,
            weights,
        });
        ordinal += 1;
    }
    let declared_terms = file.declared_terms;
    let scoring_file_header = file.header_lines.clone();
    let duplicate_header_keys = file.duplicate_keys.clone();
    let pgs_id = file.pgs_id.clone();
    let (source_sha256, source_bytes) = file.finish()?;

    // Exact duplicate lines (by 128-bit SHA-256 prefix) are all marked; none can contribute twice.
    line_hashes.sort_unstable();
    let mut duplicates = 0;
    for group in line_hashes.chunk_by(|a, b| a.0 == b.0).filter(|g| g.len() > 1) {
        for &(_, i) in group {
            terms.reasons[i as usize] |= Reason::ExactDuplicate as u32;
            duplicates += 1;
        }
    }

    let mut counts = Counts {
        exact_duplicate_terms: duplicates,
        centred_terms,
        ..Counts::default()
    };
    for i in 0..terms.len() {
        let t = terms.get(i)?;
        *counts.allele_kinds.entry(t.allele_kind.as_str().into()).or_default() += 1;
        *counts.models.entry(t.model.as_str().into()).or_default() += 1;
        *counts
            .orientation
            .entry(t.orientation.status.as_str().into())
            .or_default() += 1;
        counts.palindromic_terms += t.palindromic as u64;
        for r in t.reasons.iter() {
            *counts.review_reasons.entry(r.as_str().into()).or_default() += 1;
        }
    }
    let actual_terms = terms.len() as u64;
    let records = terms.encode();
    drop(terms);
    let catalog_terms = metadata.get("variants_number").and_then(|v| v.as_u64());
    let header = Header {
        schema: pack::SCHEMA.into(),
        pgsum_version: env!("CARGO_PKG_VERSION").into(),
        pgs_id: pgs_id.clone(),
        source: SourceFile {
            name: scoring_path
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default(),
            sha256: source_sha256,
            bytes: source_bytes,
        },
        scoring_file_header,
        duplicate_header_keys,
        columns: columns.names.clone(),
        catalog_metadata_sha256: crate::digest::file_sha256(metadata_path)?,
        license: metadata.get("license").and_then(|v| v.as_str()).map(str::to_owned),
        matches_publication: metadata.get("matches_publication").and_then(|v| v.as_bool()),
        weight_type: WeightType {
            scoring_file: header_weight_type,
            catalog: metadata.get("weight_type").and_then(|v| v.as_str()).map(str::to_owned),
        },
        catalog_metadata: metadata,
        reference: reference_identity.clone(),
        inventory: Inventory {
            declared_terms,
            catalog_terms,
            actual_terms,
            consistent: actual_terms == declared_terms && catalog_terms == Some(actual_terms),
        },
        counts,
        records_sha256: pack::records_sha256(&records),
    };
    pack::write(
        &out_dir.join(format!("{pgs_id}.{}", pack::EXTENSION)),
        &header,
        &records,
    )?;
    Ok(header)
}
