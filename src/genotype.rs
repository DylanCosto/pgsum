//! Assess one target site from the gVCF records that overlap it.
//!
//! The rules and their order are specified in `DESIGN.md` under "Genotype rules". The first failing check
//! decides the state; only `ObservedReference` and `ObservedVariant` carry a dosage.

/// A parsed gVCF record, reduced to what the genotype rules read.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    /// 1-based position.
    pub pos: u64,
    /// 1-based inclusive end: `INFO/END` when present, else `pos + len(ref) - 1`.
    pub end: u64,
    pub ref_allele: String,
    pub alts: Vec<String>,
    pub filters: Vec<String>,
    /// Allele indices from GT; `None` is a no-call (`.`).
    pub gt: Vec<Option<u32>>,
    /// Whether GT used `|`.
    pub gt_phased: bool,
    pub phase_set: Option<String>,
    pub ft: Option<String>,
    pub dp: Option<String>,
    pub min_dp: Option<String>,
    pub gq: Option<String>,
}

impl Record {
    /// Every ALT is `<NON_REF>`, `<*>` or `.`.
    pub fn is_reference_block(&self) -> bool {
        self.alts
            .iter()
            .all(|a| matches!(a.as_str(), "<NON_REF>" | "<*>" | "."))
    }
}

/// Thresholds for a passing call.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Policy {
    pub id: &'static str,
    pub min_depth: f64,
    pub min_gq: f64,
    /// Treat DeepVariant `FILTER=RefCall` on a `0/0` record as passing. Set only when the gVCF header
    /// defines the `RefCall` filter.
    pub refcall_is_reference: bool,
}

impl Default for Policy {
    fn default() -> Self {
        Policy {
            id: "pgsum-diploid-dp10-gq20-pass-v1",
            min_depth: 10.0,
            min_gq: 20.0,
            refcall_is_reference: false,
        }
    }
}

/// The outcome for one target site.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum State {
    GenotypeFiltered,
    ReferenceAnchorMismatch,
    ReferenceMismatch,
    UnknownNoRecord,
    AmbiguousOverlappingRecords,
    UnsupportedPloidy,
    NoCall,
    PartialNoCall,
    Filtered,
    QualityMissing,
    LowQuality,
    UnsupportedSymbolicGenotype,
    IncompleteReferenceSpan,
    UnsupportedAlleleRepresentation,
    OtherCalledAllele,
    ObservedReference,
    ObservedVariant,
}

impl State {
    pub const ALL: [State; 17] = [
        State::GenotypeFiltered,
        State::ReferenceAnchorMismatch,
        State::ReferenceMismatch,
        State::UnknownNoRecord,
        State::AmbiguousOverlappingRecords,
        State::UnsupportedPloidy,
        State::NoCall,
        State::PartialNoCall,
        State::Filtered,
        State::QualityMissing,
        State::LowQuality,
        State::UnsupportedSymbolicGenotype,
        State::IncompleteReferenceSpan,
        State::UnsupportedAlleleRepresentation,
        State::OtherCalledAllele,
        State::ObservedReference,
        State::ObservedVariant,
    ];

    pub fn from_code(code: u8) -> Option<State> {
        State::ALL.get(code as usize).copied()
    }

    pub fn as_str(self) -> &'static str {
        use State::*;
        match self {
            GenotypeFiltered => "genotype_filtered",
            ReferenceAnchorMismatch => "reference_anchor_mismatch",
            ReferenceMismatch => "reference_mismatch",
            UnknownNoRecord => "unknown_no_record",
            AmbiguousOverlappingRecords => "ambiguous_overlapping_records",
            UnsupportedPloidy => "unsupported_ploidy",
            NoCall => "no_call",
            PartialNoCall => "partial_no_call",
            Filtered => "filtered",
            QualityMissing => "quality_missing",
            LowQuality => "low_quality",
            UnsupportedSymbolicGenotype => "unsupported_symbolic_genotype",
            IncompleteReferenceSpan => "incomplete_reference_span",
            UnsupportedAlleleRepresentation => "unsupported_allele_representation",
            OtherCalledAllele => "other_called_allele",
            ObservedReference => "observed_reference",
            ObservedVariant => "observed_variant",
        }
    }

    pub fn is_passing(self) -> bool {
        matches!(self, State::ObservedReference | State::ObservedVariant)
    }
}

/// The oriented target: the term's REF and ALT at `pos` (1-based).
#[derive(Clone, Debug, PartialEq)]
pub struct Target<'a> {
    pub pos: u64,
    pub ref_allele: &'a str,
    pub alt: &'a str,
}

/// An assessed call.
#[derive(Clone, Debug, PartialEq)]
pub struct Call {
    pub state: State,
    /// Copies of the target ALT (0, 1 or 2); set only for passing states.
    pub alt_dosage: Option<u8>,
    /// Phase set, kept only when GT is phased and `PS` is numeric.
    pub phase_set: Option<String>,
    /// Whether the DeepVariant `RefCall` convention was applied.
    pub refcall_adapted: bool,
}

impl Call {
    fn state(state: State) -> Self {
        Call {
            state,
            alt_dosage: None,
            phase_set: None,
            refcall_adapted: false,
        }
    }
}

/// Assess `target` from the records overlapping it.
///
/// `reference` returns the FASTA bases for a 1-based start and a length; it is called for the target span
/// and for each reference block's REF base.
pub fn assess(
    records: &[Record],
    target: &Target<'_>,
    policy: &Policy,
    reference: impl Fn(u64, usize) -> Option<String>,
) -> Call {
    let mut refcall_adapted = false;
    let mut filters: Vec<&[String]> = Vec::with_capacity(records.len());
    for record in records {
        if record.ft.as_deref().is_some_and(|ft| ft != "PASS" && ft != ".") {
            return Call::state(State::GenotypeFiltered);
        }
        if record.is_reference_block()
            && reference(record.pos, record.ref_allele.len()).as_deref() != Some(record.ref_allele.as_str())
        {
            return Call::state(State::ReferenceAnchorMismatch);
        }
        if policy.refcall_is_reference && record.filters == ["RefCall"] && record.gt == [Some(0), Some(0)] {
            refcall_adapted = true;
            filters.push(&[]);
        } else {
            filters.push(&record.filters);
        }
    }
    let mut call = assess_site(records, &filters, target, policy, &reference);
    call.refcall_adapted = refcall_adapted;
    call
}

fn assess_site(
    records: &[Record],
    filters: &[&[String]],
    target: &Target<'_>,
    policy: &Policy,
    reference: &impl Fn(u64, usize) -> Option<String>,
) -> Call {
    if reference(target.pos, target.ref_allele.len()).as_deref() != Some(target.ref_allele) {
        return Call::state(State::ReferenceMismatch);
    }
    let (record, filters) = match records {
        [] => return Call::state(State::UnknownNoRecord),
        [record] => (record, filters[0]),
        _ => return Call::state(State::AmbiguousOverlappingRecords),
    };
    if record.gt.len() != 2 {
        return Call::state(State::UnsupportedPloidy);
    }
    let called = record.gt.iter().filter(|a| a.is_some()).count();
    if called == 0 {
        return Call::state(State::NoCall);
    }
    if called < record.gt.len() {
        return Call::state(State::PartialNoCall);
    }
    if filters.iter().any(|f| f != "PASS" && f != ".") {
        return Call::state(State::Filtered);
    }
    if record.ft.as_deref().is_some_and(|ft| ft != "PASS" && ft != ".") {
        return Call::state(State::Filtered);
    }
    let block = record.is_reference_block();
    let depth = number(if block { &record.min_dp } else { &record.dp });
    let gq = number(&record.gq);
    let (Some(depth), Some(gq)) = (depth, gq) else {
        return Call::state(State::QualityMissing);
    };
    if depth < policy.min_depth || gq < policy.min_gq {
        return Call::state(State::LowQuality);
    }
    let indices: Vec<u32> = record.gt.iter().map(|a| a.expect("no-calls handled above")).collect();
    let dosage = if block {
        if indices.iter().any(|&i| i != 0) {
            return Call::state(State::UnsupportedSymbolicGenotype);
        }
        if record.end < target.pos + target.ref_allele.len() as u64 - 1 {
            return Call::state(State::IncompleteReferenceSpan);
        }
        0
    } else {
        if record.pos != target.pos || record.ref_allele != target.ref_allele {
            return Call::state(State::UnsupportedAlleleRepresentation);
        }
        let mut dosage = 0u8;
        for &i in &indices {
            if i == 0 {
                continue;
            }
            match record.alts.get(i as usize - 1) {
                Some(alt) if alt == target.alt => dosage += 1,
                _ => return Call::state(State::OtherCalledAllele),
            }
        }
        dosage
    };
    let phased = record.gt_phased
        && record
            .phase_set
            .as_deref()
            .is_some_and(|ps| !ps.is_empty() && ps.chars().all(|c| c.is_ascii_digit()));
    Call {
        state: if dosage == 0 {
            State::ObservedReference
        } else {
            State::ObservedVariant
        },
        alt_dosage: Some(dosage),
        phase_set: if phased { record.phase_set.clone() } else { None },
        refcall_adapted: false,
    }
}

/// A finite, non-negative number, else `None`.
fn number(value: &Option<String>) -> Option<f64> {
    value
        .as_deref()?
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite() && *v >= 0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SEQ: &str = "ACGTACGTAC"; // positions 1..=10

    fn reference(pos: u64, len: usize) -> Option<String> {
        let start = pos as usize - 1;
        SEQ.get(start..start + len).map(str::to_owned)
    }

    fn variant(gt: [Option<u32>; 2]) -> Record {
        Record {
            pos: 3,
            end: 3,
            ref_allele: "G".into(),
            alts: vec!["T".into(), "<*>".into()],
            filters: vec!["PASS".into()],
            gt: gt.to_vec(),
            gt_phased: false,
            phase_set: None,
            ft: None,
            dp: Some("30".into()),
            min_dp: None,
            gq: Some("40".into()),
        }
    }

    fn block(pos: u64, end: u64) -> Record {
        Record {
            pos,
            end,
            ref_allele: reference(pos, 1).unwrap(),
            alts: vec!["<*>".into()],
            filters: vec![".".into()],
            gt: vec![Some(0), Some(0)],
            gt_phased: false,
            phase_set: None,
            ft: None,
            dp: None,
            min_dp: Some("25".into()),
            gq: Some("50".into()),
        }
    }

    const TARGET: Target<'static> = Target {
        pos: 3,
        ref_allele: "G",
        alt: "T",
    };

    fn state(records: &[Record]) -> State {
        assess(records, &TARGET, &Policy::default(), reference).state
    }

    #[test]
    fn het_and_hom_variant_dosages() {
        let het = assess(&[variant([Some(0), Some(1)])], &TARGET, &Policy::default(), reference);
        assert_eq!((het.state, het.alt_dosage), (State::ObservedVariant, Some(1)));
        let hom = assess(&[variant([Some(1), Some(1)])], &TARGET, &Policy::default(), reference);
        assert_eq!(hom.alt_dosage, Some(2));
    }

    #[test]
    fn covering_reference_block_is_homozygous_reference() {
        let call = assess(&[block(1, 8)], &TARGET, &Policy::default(), reference);
        assert_eq!((call.state, call.alt_dosage), (State::ObservedReference, Some(0)));
    }

    #[test]
    fn reference_block_uses_min_dp() {
        let mut r = block(1, 8);
        r.min_dp = Some("9".into());
        assert_eq!(state(&[r.clone()]), State::LowQuality);
        r.min_dp = None;
        r.dp = Some("50".into());
        assert_eq!(state(&[r]), State::QualityMissing);
    }

    #[test]
    fn no_record_is_unknown_not_reference() {
        assert_eq!(state(&[]), State::UnknownNoRecord);
    }

    #[test]
    fn overlapping_records_are_ambiguous() {
        assert_eq!(
            state(&[block(1, 8), variant([Some(0), Some(1)])]),
            State::AmbiguousOverlappingRecords
        );
    }

    #[test]
    fn target_reference_mismatch() {
        let target = Target {
            pos: 3,
            ref_allele: "A",
            alt: "T",
        };
        assert_eq!(
            assess(&[block(1, 8)], &target, &Policy::default(), reference).state,
            State::ReferenceMismatch
        );
    }

    #[test]
    fn block_anchor_mismatch() {
        let mut r = block(1, 8);
        r.ref_allele = "T".into();
        assert_eq!(state(&[r]), State::ReferenceAnchorMismatch);
    }

    #[test]
    fn ploidy_and_no_calls() {
        let mut haploid = variant([Some(1), None]);
        haploid.gt = vec![Some(1)];
        assert_eq!(state(&[haploid]), State::UnsupportedPloidy);
        assert_eq!(state(&[variant([None, None])]), State::NoCall);
        assert_eq!(state(&[variant([Some(0), None])]), State::PartialNoCall);
    }

    #[test]
    fn filters() {
        let mut r = variant([Some(0), Some(1)]);
        r.filters = vec!["LowQual".into()];
        assert_eq!(state(&[r]), State::Filtered);
        let mut r = variant([Some(0), Some(1)]);
        r.ft = Some("DP".into());
        assert_eq!(state(&[r]), State::GenotypeFiltered);
    }

    #[test]
    fn refcall_convention_only_when_enabled() {
        let mut r = variant([Some(0), Some(0)]);
        r.filters = vec!["RefCall".into()];
        assert_eq!(state(&[r.clone()]), State::Filtered);
        let policy = Policy {
            refcall_is_reference: true,
            ..Policy::default()
        };
        let call = assess(&[r], &TARGET, &policy, reference);
        assert_eq!((call.state, call.refcall_adapted), (State::ObservedReference, true));
    }

    #[test]
    fn low_quality() {
        let mut r = variant([Some(0), Some(1)]);
        r.gq = Some("19".into());
        assert_eq!(state(&[r]), State::LowQuality);
    }

    #[test]
    fn block_with_nonreference_gt() {
        let mut r = block(1, 8);
        r.gt = vec![Some(0), Some(1)];
        assert_eq!(state(&[r]), State::UnsupportedSymbolicGenotype);
    }

    #[test]
    fn block_must_cover_whole_target() {
        let target = Target {
            pos: 3,
            ref_allele: "GTA",
            alt: "G",
        };
        let call = assess(&[block(1, 4)], &target, &Policy::default(), reference);
        assert_eq!(call.state, State::IncompleteReferenceSpan);
    }

    #[test]
    fn variant_representation_and_other_allele() {
        let mut r = variant([Some(0), Some(1)]);
        r.pos = 2;
        r.ref_allele = "CG".into();
        assert_eq!(state(&[r]), State::UnsupportedAlleleRepresentation);
        let mut r = variant([Some(0), Some(1)]);
        r.alts = vec!["A".into()];
        assert_eq!(state(&[r]), State::OtherCalledAllele);
    }

    #[test]
    fn phase_kept_only_with_numeric_phase_set() {
        let mut r = variant([Some(0), Some(1)]);
        r.gt_phased = true;
        r.phase_set = Some("1234".into());
        let call = assess(&[r.clone()], &TARGET, &Policy::default(), reference);
        assert_eq!(call.phase_set.as_deref(), Some("1234"));
        r.phase_set = Some(".".into());
        assert_eq!(assess(&[r], &TARGET, &Policy::default(), reference).phase_set, None);
    }
}
