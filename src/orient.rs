//! Orient a term's effect and other alleles against the reference genome.
//!
//! Only SNVs are resolved: exactly one of {direct, complemented} allele pairs must contain the reference base.
//! Palindromic pairs (A/T, C/G) match both ways and stay unresolved. Longer alleles need sequence evidence
//! that the scoring file does not carry.

use crate::Result;
use crate::reference::Reference;
use crate::term::{AlleleKind, CONTIGS, Description};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum Status {
    Resolved = 0,
    UnresolvedPosition = 1,
    UnresolvedSourceAlleles = 2,
    OutsideReference = 3,
    SequenceReferenceEvidenceRequired = 4,
    PalindromicOrientationUnresolved = 5,
    AlleleRepresentationUnresolved = 6,
}

impl Status {
    const ALL: [Status; 7] = [
        Status::Resolved,
        Status::UnresolvedPosition,
        Status::UnresolvedSourceAlleles,
        Status::OutsideReference,
        Status::SequenceReferenceEvidenceRequired,
        Status::PalindromicOrientationUnresolved,
        Status::AlleleRepresentationUnresolved,
    ];
    pub fn as_str(self) -> &'static str {
        [
            "resolved",
            "unresolved_position",
            "unresolved_source_alleles",
            "outside_reference",
            "sequence_reference_evidence_required",
            "palindromic_orientation_unresolved",
            "allele_representation_unresolved",
        ][self as usize]
    }
    pub fn from_code(code: u8) -> Option<Status> {
        Status::ALL.get(code as usize).copied()
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum Method {
    None = 0,
    Direct = 1,
    Complement = 2,
}

impl Method {
    pub fn as_str(self) -> &'static str {
        ["", "direct", "complement"][self as usize]
    }
    pub fn from_code(code: u8) -> Option<Method> {
        [Method::None, Method::Direct, Method::Complement]
            .get(code as usize)
            .copied()
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Orientation {
    pub status: Status,
    pub method: Method,
    /// Reference and alternate base on the genome's forward strand; 0 unless resolved.
    pub ref_base: u8,
    pub alt_base: u8,
    /// The model's effect allele is the ALT (so effect dosage = ALT dosage); otherwise it is the REF.
    pub effect_is_alt: bool,
}

impl Orientation {
    fn status(status: Status) -> Orientation {
        Orientation {
            status,
            method: Method::None,
            ref_base: 0,
            alt_base: 0,
            effect_is_alt: false,
        }
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

pub fn orient(d: &Description<'_>, reference: &Reference) -> Result<Orientation> {
    if d.contig == 0 || d.pos == 0 {
        return Ok(Orientation::status(Status::UnresolvedPosition));
    }
    if d.allele_kind == AlleleKind::Unresolved {
        return Ok(Orientation::status(Status::UnresolvedSourceAlleles));
    }
    let contig = CONTIGS[d.contig as usize - 1];
    let span = d.effect_allele.len().max(d.other_allele.len()) as u64;
    match reference.contig_length(contig) {
        Some(length) if d.pos as u64 + span - 1 <= length => {}
        _ => return Ok(Orientation::status(Status::OutsideReference)),
    }
    if d.allele_kind != AlleleKind::LiteralSnv {
        return Ok(Orientation::status(Status::SequenceReferenceEvidenceRequired));
    }
    let ref_base = reference.fetch(contig, d.pos as u64 - 1, d.pos as u64)?.as_bytes()[0];
    let (effect, other) = (d.effect_allele.as_bytes()[0], d.other_allele.as_bytes()[0]);
    let mut found = None;
    let mut matches = 0;
    for (method, e, o) in [
        (Method::Direct, effect, other),
        (Method::Complement, complement(effect), complement(other)),
    ] {
        if ref_base == e || ref_base == o {
            matches += 1;
            found = Some(Orientation {
                status: Status::Resolved,
                method,
                ref_base,
                alt_base: if e == ref_base { o } else { e },
                effect_is_alt: e != ref_base,
            });
        }
    }
    Ok(match (matches, found) {
        (1, Some(o)) => o,
        _ if d.palindromic => Orientation::status(Status::PalindromicOrientationUnresolved),
        _ => Orientation::status(Status::AlleleRepresentationUnresolved),
    })
}
