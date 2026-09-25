//! Describe one scoring-file term: its position, weight model, allele kind and any reasons it needs review.
//!
//! The rules are listed in `DESIGN.md` under "Term description". A term with any review reason never
//! contributes to a score.

use crate::decimal::Decimal;
use crate::scoring_file::{Columns, Row};

/// Canonical contigs, in the order used by pack contig codes (code = index + 1; 0 = unresolved).
pub const CONTIGS: [&str; 25] = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
    "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM",
];

/// `hm_chr` value → contig code.
fn contig_code(hm_chr: &str) -> u8 {
    match hm_chr {
        "X" => 23,
        "Y" => 24,
        "MT" | "M" => 25,
        _ => match hm_chr.parse::<u8>() {
            Ok(n @ 1..=22) if hm_chr.as_bytes()[0] != b'0' && !hm_chr.starts_with('+') => n,
            _ => 0,
        },
    }
}

macro_rules! reasons {
    ($($variant:ident = $bit:literal => $name:literal,)*) => {
        /// Why a term needs review. Stored as a bit set in packs.
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        #[repr(u32)]
        pub enum Reason { $($variant = 1 << $bit,)* }

        impl Reason {
            pub const ALL: &[Reason] = &[$(Reason::$variant,)*];
            pub fn as_str(self) -> &'static str {
                match self { $(Reason::$variant => $name,)* }
            }
        }
    };
}

reasons! {
    InvalidFlagHaplotype = 0 => "invalid_model_flag:is_haplotype",
    InvalidFlagDiplotype = 1 => "invalid_model_flag:is_diplotype",
    InvalidFlagInteraction = 2 => "invalid_model_flag:is_interaction",
    InvalidFlagDominant = 3 => "invalid_model_flag:is_dominant",
    InvalidFlagRecessive = 4 => "invalid_model_flag:is_recessive",
    SpecialHaplotype = 5 => "requires_special_model:is_haplotype",
    SpecialDiplotype = 6 => "requires_special_model:is_diplotype",
    SpecialInteraction = 7 => "requires_special_model:is_interaction",
    ConditionalInclusion = 8 => "conditional_inclusion_requires_review",
    VariantDescription = 9 => "variant_description_requires_review",
    ImputationMethod = 10 => "specific_imputation_method_requires_review",
    ConflictingWeightModels = 11 => "conflicting_weight_models",
    InvalidDosageWeights = 12 => "invalid_dosage_weights",
    InvalidEffectWeight = 13 => "invalid_effect_weight",
    UnresolvedPosition = 14 => "unresolved_harmonized_position",
    MismatchChr = 15 => "harmonization_mismatch:hm_match_chr",
    MismatchPos = 16 => "harmonization_mismatch:hm_match_pos",
    InvalidMatchChr = 17 => "invalid_harmonization_flag:hm_match_chr",
    InvalidMatchPos = 18 => "invalid_harmonization_flag:hm_match_pos",
    OtherAlleleMissing = 19 => "author_other_allele_missing",
    NonliteralAlleles = 20 => "nonliteral_alleles",
    IdenticalAlleles = 21 => "identical_effect_other_alleles",
    ExactDuplicate = 22 => "exact_duplicate_source_term",
}

/// A set of review reasons.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Reasons(pub u32);

impl Reasons {
    pub fn insert(&mut self, reason: Reason) {
        self.0 |= reason as u32;
    }
    pub fn is_empty(self) -> bool {
        self.0 == 0
    }
    pub fn iter(self) -> impl Iterator<Item = Reason> {
        Reason::ALL.iter().copied().filter(move |r| self.0 & *r as u32 != 0)
    }
    /// Reason names, sorted.
    pub fn names(self) -> Vec<&'static str> {
        let mut names: Vec<_> = self.iter().map(Reason::as_str).collect();
        names.sort_unstable();
        names
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum Model {
    Additive = 0,
    Dominant = 1,
    Recessive = 2,
    DosageWeights = 3,
}

impl Model {
    pub fn as_str(self) -> &'static str {
        ["additive", "dominant", "recessive", "dosage_weights"][self as usize]
    }
    pub fn from_code(code: u8) -> Option<Model> {
        [Model::Additive, Model::Dominant, Model::Recessive, Model::DosageWeights]
            .get(code as usize)
            .copied()
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum AlleleKind {
    LiteralSnv = 0,
    LiteralSequence = 1,
    Unresolved = 2,
}

impl AlleleKind {
    pub fn as_str(self) -> &'static str {
        ["literal_snv", "literal_sequence", "unresolved"][self as usize]
    }
    pub fn from_code(code: u8) -> Option<AlleleKind> {
        [
            AlleleKind::LiteralSnv,
            AlleleKind::LiteralSequence,
            AlleleKind::Unresolved,
        ]
        .get(code as usize)
        .copied()
    }
}

/// A recognised `variant_description` of the form `weight_type=<type>;centre=<0..2>` (LDAK/MegaPRS).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Centred {
    pub weight_type: String,
    pub centre: String,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Description<'a> {
    /// Contig code (index into `CONTIGS` + 1), 0 when unresolved.
    pub contig: u8,
    /// 1-based harmonized position, 0 when unresolved.
    pub pos: u32,
    pub model: Model,
    pub allele_kind: AlleleKind,
    pub palindromic: bool,
    pub reasons: Reasons,
    pub effect_allele: &'a str,
    pub other_allele: &'a str,
    /// One weight (additive/dominant/recessive) or three (dosage 0, 1, 2); `None` where invalid.
    pub weights: Vec<Option<Decimal>>,
    pub centred: Option<Centred>,
}

fn is_acgt(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

/// `[+-]`-free decimal text as accepted inside a centred description.
fn unsigned_number(s: &str) -> bool {
    !s.starts_with(['+', '-']) && Decimal::parse(s).is_some()
}

fn centred(description: &str) -> Option<Centred> {
    let rest = description.strip_prefix("weight_type=")?;
    let (weight_type, centre) = rest.split_once(";centre=")?;
    if !matches!(weight_type, "beta" | "OR" | "HR" | "NR") || !unsigned_number(centre) {
        return None;
    }
    let value = Decimal::parse(centre)?;
    // 0 <= centre <= 2: non-negative by the grammar; compare the magnitude with 2.
    let at_most_two = value.coefficient == "0"
        || value.adjusted() < 0
        || value.adjusted() == 0 && {
            let digits = value.coefficient.as_bytes();
            digits[0] < b'2' || digits[0] == b'2' && digits[1..].iter().all(|&d| d == b'0')
        };
    at_most_two.then(|| Centred {
        weight_type: weight_type.to_owned(),
        centre: centre.to_owned(),
    })
}

pub fn describe<'a>(row: &Row<'a>, columns: &Columns) -> Description<'a> {
    let mut reasons = Reasons::default();
    let mut flags = [false; 5];
    const INVALID_FLAG: [Reason; 5] = [
        Reason::InvalidFlagHaplotype,
        Reason::InvalidFlagDiplotype,
        Reason::InvalidFlagInteraction,
        Reason::InvalidFlagDominant,
        Reason::InvalidFlagRecessive,
    ];
    for (i, column) in columns.flags.iter().enumerate() {
        let raw = row.opt(*column);
        if raw.eq_ignore_ascii_case("TRUE") {
            flags[i] = true;
        } else if raw.eq_ignore_ascii_case("FALSE") || matches!(raw, "" | "." | "NA") {
        } else {
            reasons.insert(INVALID_FLAG[i]);
        }
    }
    let [haplotype, diplotype, interaction, dominant, recessive] = flags;
    for (set, reason) in [
        (haplotype, Reason::SpecialHaplotype),
        (diplotype, Reason::SpecialDiplotype),
        (interaction, Reason::SpecialInteraction),
    ] {
        if set {
            reasons.insert(reason);
        }
    }
    if !row.opt(columns.inclusion_criteria).trim().is_empty() {
        reasons.insert(Reason::ConditionalInclusion);
    }
    let description = row.opt(columns.variant_description).trim();
    let centred = if description.is_empty() {
        None
    } else {
        centred(description)
    };
    if !description.is_empty() && centred.is_none() {
        reasons.insert(Reason::VariantDescription);
    }
    if !row.opt(columns.imputation_method).trim().is_empty() {
        reasons.insert(Reason::ImputationMethod);
    }
    let dosage_text: [&str; 3] = columns.dosage_weights.map(|i| row.opt(i));
    let explicit = dosage_text.iter().any(|w| !w.is_empty());
    let model = if explicit {
        Model::DosageWeights
    } else if dominant {
        Model::Dominant
    } else if recessive {
        Model::Recessive
    } else {
        Model::Additive
    };
    if dominant && recessive || explicit && (dominant || recessive) {
        reasons.insert(Reason::ConflictingWeightModels);
    }
    let weights: Vec<Option<Decimal>> = if explicit {
        let w: Vec<_> = dosage_text.iter().map(|t| Decimal::parse(t)).collect();
        if w.iter().any(Option::is_none) {
            reasons.insert(Reason::InvalidDosageWeights);
        }
        w
    } else {
        let w = Decimal::parse(row.opt(columns.effect_weight));
        if w.is_none() {
            reasons.insert(Reason::InvalidEffectWeight);
        }
        vec![w]
    };
    let contig = contig_code(row.get(columns.hm_chr));
    let pos_text = row.get(columns.hm_pos);
    let pos = if (1..=10).contains(&pos_text.len())
        && pos_text.as_bytes()[0] != b'0'
        && pos_text.bytes().all(|b| b.is_ascii_digit())
    {
        pos_text.parse::<u64>().ok().filter(|&p| p < 1 << 31).unwrap_or(0) as u32
    } else {
        0
    };
    if contig == 0 || pos == 0 {
        reasons.insert(Reason::UnresolvedPosition);
    }
    for (column, mismatch, invalid) in [
        (columns.hm_match_chr, Reason::MismatchChr, Reason::InvalidMatchChr),
        (columns.hm_match_pos, Reason::MismatchPos, Reason::InvalidMatchPos),
    ] {
        let value = row.opt(column).to_lowercase();
        if value == "false" {
            reasons.insert(mismatch);
        } else if !matches!(value.as_str(), "" | "true" | "na" | ".") {
            reasons.insert(invalid);
        }
    }
    let effect = row.get(columns.effect_allele);
    let other = row.opt(columns.other_allele);
    if other.is_empty() {
        reasons.insert(Reason::OtherAlleleMissing);
    }
    if !is_acgt(effect) || !other.is_empty() && !is_acgt(other) {
        reasons.insert(Reason::NonliteralAlleles);
    }
    if effect == other {
        reasons.insert(Reason::IdenticalAlleles);
    }
    let allele_kind = if effect.len() == 1 && other.len() == 1 && is_acgt(effect) && is_acgt(other) && effect != other {
        AlleleKind::LiteralSnv
    } else if is_acgt(effect) && is_acgt(other) && effect != other {
        AlleleKind::LiteralSequence
    } else {
        AlleleKind::Unresolved
    };
    let palindromic = matches!((effect, other), ("A", "T") | ("T", "A") | ("C", "G") | ("G", "C"));
    Description {
        contig,
        pos,
        model,
        allele_kind,
        palindromic,
        reasons,
        effect_allele: effect,
        other_allele: other,
        weights,
        centred,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contig_codes() {
        assert_eq!(contig_code("1"), 1);
        assert_eq!(contig_code("22"), 22);
        assert_eq!(contig_code("X"), 23);
        assert_eq!(contig_code("MT"), 25);
        for bad in ["", "0", "23", "01", "+1", "chr1", "x"] {
            assert_eq!(contig_code(bad), 0, "{bad:?}");
        }
    }

    #[test]
    fn centred_descriptions() {
        assert!(centred("weight_type=beta;centre=1.5").is_some());
        assert!(centred("weight_type=NR;centre=2").is_some());
        assert!(centred("weight_type=NR;centre=2.000").is_some());
        assert!(centred("weight_type=NR;centre=0").is_some());
        assert!(centred("weight_type=NR;centre=2.01").is_none());
        assert!(centred("weight_type=NR;centre=-0.1").is_none());
        assert!(centred("weight_type=NR;centre=+1").is_none());
        assert!(centred("weight_type=logOR;centre=1").is_none());
        assert!(centred("weight_type=beta;centre=1;x=2").is_none());
    }

    #[test]
    fn reason_names_sort() {
        let mut r = Reasons::default();
        r.insert(Reason::UnresolvedPosition);
        r.insert(Reason::ConflictingWeightModels);
        assert_eq!(
            r.names(),
            ["conflicting_weight_models", "unresolved_harmonized_position"]
        );
    }
}
