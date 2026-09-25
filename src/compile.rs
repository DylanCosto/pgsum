//! Compile a harmonized scoring file and its Catalog metadata into a pack.

use std::collections::BTreeMap;
use std::path::Path;

use sha2::{Digest, Sha256};

use crate::digest::file_sha256;
use crate::genotypes::ANY_ALT;
use crate::orient::{Method, Orientation, Status, orient};
use crate::pack::{
    self, Columns, Counts, Header, Inference, InferenceSummary, Inventory, ReferenceIdentity, SourceFile, TermRecord,
    Weight, WeightType,
};
use crate::reference::Reference;
use crate::scoring_file::{Origin, ScoringFile};
use crate::term::{AlleleKind, CONTIGS, Description, Reason, describe};
use crate::{Error, Result, invalid};

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

/// Share of eligible terms that must follow one reference convention for it to be used.
pub const REFERENCE_CONVENTION_THRESHOLD: f64 = 0.99;

/// What is known about a term without an author `other_allele`.
#[derive(Clone, Copy, Debug)]
struct Candidate {
    effect: u8,
    /// The reference base at the term's position, when it is A, C, G or T.
    reference: Option<u8>,
    /// The orientation from the Catalog's inferred other allele, when it names one base and resolves.
    catalog: Option<Orientation>,
}

/// A term is a candidate when it has no author other allele and its effect allele is one base at a position
/// inside the reference.
fn inference_candidate(d: &Description<'_>, catalog_other: &str, reference: &Reference) -> Result<Option<Candidate>> {
    let effect = d.effect_allele.as_bytes();
    if !d.other_allele.is_empty() || effect.len() != 1 || !b"ACGT".contains(&effect[0]) || d.contig == 0 || d.pos == 0 {
        return Ok(None);
    }
    let contig = CONTIGS[d.contig as usize - 1];
    if reference.contig_length(contig).is_none_or(|len| d.pos as u64 > len) {
        return Ok(None);
    }
    let base = reference.fetch(contig, d.pos as u64 - 1, d.pos as u64)?.as_bytes()[0];
    let catalog = if catalog_other.len() == 1
        && b"ACGT".contains(&catalog_other.as_bytes()[0])
        && catalog_other != d.effect_allele
    {
        let mut with_other = d.clone();
        with_other.other_allele = catalog_other;
        with_other.allele_kind = AlleleKind::LiteralSnv;
        with_other.palindromic = matches!(
            (d.effect_allele, catalog_other),
            ("A", "T") | ("T", "A") | ("C", "G") | ("G", "C")
        );
        Some(orient(&with_other, reference)?).filter(|o| o.status == Status::Resolved)
    } else {
        None
    };
    Ok(Some(Candidate {
        effect: effect[0],
        reference: (base != b'N').then_some(base),
        catalog,
    }))
}

/// Decide the score's reference convention and fill in each candidate's inferred orientation: the
/// reference-anchored one when the term follows the convention, otherwise the Catalog's.
fn infer(terms: &mut Columns, candidates: &[(u32, Candidate)]) -> Option<InferenceSummary> {
    if candidates.is_empty() {
        return None;
    }
    let with_reference: Vec<_> = candidates
        .iter()
        .filter_map(|(_, c)| c.reference.map(|r| (c.effect, r)))
        .collect();
    let eligible = with_reference.len() as u64;
    let not_reference = with_reference.iter().filter(|(e, r)| e != r).count() as u64;
    let is_reference = eligible - not_reference;
    let reaches = |n: u64| eligible > 0 && n as f64 >= REFERENCE_CONVENTION_THRESHOLD * eligible as f64;
    let convention = if reaches(not_reference) {
        Some(Inference::ReferenceAnchoredEffectIsAlt)
    } else if reaches(is_reference) {
        Some(Inference::ReferenceAnchoredEffectIsRef)
    } else {
        None
    };
    let mut summary = InferenceSummary {
        eligible_terms: eligible,
        effect_not_reference: not_reference,
        effect_is_reference: is_reference,
        reference_convention: convention.map(|c| match c {
            Inference::ReferenceAnchoredEffectIsAlt => "effect_is_alt".to_owned(),
            _ => "effect_is_ref".to_owned(),
        }),
        threshold: REFERENCE_CONVENTION_THRESHOLD,
        terms: BTreeMap::new(),
    };
    let anchored = |c: &Candidate| -> Option<(Inference, Orientation)> {
        let r = c.reference?;
        let o = |alt_base, effect_is_alt| Orientation {
            status: Status::Resolved,
            method: Method::Direct,
            ref_base: r,
            alt_base,
            effect_is_alt,
        };
        match convention? {
            Inference::ReferenceAnchoredEffectIsAlt if c.effect != r => {
                Some((Inference::ReferenceAnchoredEffectIsAlt, o(c.effect, true)))
            }
            Inference::ReferenceAnchoredEffectIsRef if c.effect == r => {
                Some((Inference::ReferenceAnchoredEffectIsRef, o(ANY_ALT, false)))
            }
            _ => None,
        }
    };
    for (i, c) in candidates {
        let chosen = anchored(c).or(c.catalog.map(|o| (Inference::CatalogInferredOtherAllele, o)));
        if let Some((kind, o)) = chosen {
            let i = *i as usize;
            terms.inferred_kind[i] = kind as u8;
            terms.inferred_method[i] = o.method as u8;
            terms.inferred_ref[i] = o.ref_base;
            terms.inferred_alt[i] = o.alt_base;
            terms.flags[i] |= (o.effect_is_alt as u8) << 2;
            *summary.terms.entry(kind.as_str().to_owned()).or_default() += 1;
        }
    }
    Some(summary)
}

/// Compile one scoring file into `<out_dir>/<pgs_id>.pgsp` and return the pack header.
///
/// The metadata record is `metadata_path`, else `<pgs_id>.metadata.json` next to the scoring file. A Catalog
/// score needs one; a custom score without one gets a record built from its header.
pub fn compile_file(
    scoring_path: &Path,
    metadata_path: Option<&Path>,
    reference: &Reference,
    reference_identity: &ReferenceIdentity,
    out_dir: &Path,
) -> Result<Header> {
    let mut file = ScoringFile::open(scoring_path)?;
    let sibling = scoring_path.with_file_name(format!("{}.metadata.json", file.pgs_id));
    let metadata_path = metadata_path
        .map(Path::to_path_buf)
        .or_else(|| sibling.exists().then_some(sibling));
    let metadata: serde_json::Value = match (&metadata_path, file.origin) {
        (Some(path), _) => {
            let bytes = std::fs::read(path).map_err(Error::io(path))?;
            let value: serde_json::Value =
                serde_json::from_slice(&bytes).map_err(|e| Error::Invalid(format!("{}: {e}", path.display())))?;
            if value.get("id").and_then(|v| v.as_str()) != Some(file.pgs_id.as_str()) {
                return invalid!("{}: metadata is not for {}", path.display(), file.pgs_id);
            }
            value
        }
        (None, Origin::PgsCatalog) => {
            return invalid!(
                "{}: a Catalog score needs its metadata record ({}.metadata.json from the Catalog REST API)",
                scoring_path.display(),
                file.pgs_id
            );
        }
        (None, Origin::Custom) => {
            let get = |k: &str| file.metadata.get(k).cloned();
            serde_json::json!({
                "id": file.pgs_id,
                "name": get("pgs_name"),
                "weight_type": get("weight_type"),
                "license": get("license"),
                "genome_build": get("genome_build"),
                "source": "custom scoring file header",
            })
        }
    };
    let header_weight_type = file.metadata.get("weight_type").cloned();
    let columns = file.columns.clone();
    let mut terms = Columns::default();
    let mut line_hashes: Vec<([u8; 16], u32)> = Vec::with_capacity(file.declared_terms.unwrap_or(0) as usize);
    let mut fields = Vec::new();
    let mut centred_terms = 0;
    let mut candidates: Vec<(u32, Candidate)> = Vec::new();
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
            inferred: None,
        });
        if let Some(c) = inference_candidate(&d, row.opt(columns.hm_infer_other_allele), reference)
            .map_err(|e| Error::Invalid(format!("{}: term {}: {e}", scoring_path.display(), ordinal + 1)))?
        {
            candidates.push((ordinal, c));
        }
        ordinal += 1;
    }
    let inference = infer(&mut terms, &candidates);
    drop(candidates);
    let declared_terms = file.declared_terms;
    let origin = file.origin;
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
        origin,
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
        catalog_metadata_sha256: metadata_path.as_deref().map(crate::digest::file_sha256).transpose()?,
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
            consistent: match origin {
                Origin::PgsCatalog => declared_terms == Some(actual_terms) && catalog_terms == Some(actual_terms),
                Origin::Custom => {
                    declared_terms.is_none_or(|n| n == actual_terms) && catalog_terms.is_none_or(|n| n == actual_terms)
                }
            },
        },
        counts,
        inference,
        records_sha256: pack::records_sha256(&records),
    };
    pack::write(
        &out_dir.join(format!("{pgs_id}.{}", pack::EXTENSION)),
        &header,
        &records,
    )?;
    Ok(header)
}
