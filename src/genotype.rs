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
    /// Read a haploid call on chrX or chrY (`1`) as homozygous (`1/1`), as DeepVariant writes male chrX.
    pub haploid_xy_as_homozygous: bool,
    /// Accept a variant-record call that reports neither depth nor GQ (a genotype-only VCF: imputed,
    /// array or joint-called data) on its GT and FILTER alone. A record reporting either is still checked.
    pub accept_missing_quality: bool,
    /// Read several records starting at the target's position as one multi-allelic record (see
    /// `merge_split`), as when a panel splits multi-allelic sites into one record per ALT.
    pub merge_split_records: bool,
}

impl Default for Policy {
    fn default() -> Self {
        Policy {
            id: "pgsum-diploid-dp10-gq20-pass-v1",
            min_depth: 10.0,
            min_gq: 20.0,
            refcall_is_reference: false,
            haploid_xy_as_homozygous: false,
            accept_missing_quality: false,
            merge_split_records: false,
        }
    }
}

impl Policy {
    /// The default policy, with DeepVariant's `RefCall` convention when the gVCF defines it and, optionally,
    /// haploid chrX/chrY calls read as homozygous and genotype-only calls accepted (each changes the policy ID).
    pub fn for_gvcf(
        refcall_is_reference: bool,
        haploid_xy_as_homozygous: bool,
        accept_missing_quality: bool,
    ) -> Policy {
        Policy {
            id: match (haploid_xy_as_homozygous, accept_missing_quality) {
                (false, false) => Policy::default().id,
                (true, false) => "pgsum-dp10-gq20-pass-haploid-xy-homozygous-v1",
                (false, true) => "pgsum-diploid-dp10-gq20-pass-or-genotype-only-v1",
                (true, true) => "pgsum-dp10-gq20-pass-or-genotype-only-haploid-xy-homozygous-v1",
            },
            refcall_is_reference,
            haploid_xy_as_homozygous,
            accept_missing_quality,
            ..Policy::default()
        }
    }
}

impl Policy {
    /// The same policy, reading split multi-allelic records as one (`merge_split_records`).
    pub fn with_merge_split(mut self, on: bool) -> Policy {
        self.merge_split_records = on;
        self
    }
}

/// Whether every ALT of a record is an indel that keeps its first base, so the record says nothing about the
/// base at its position (an insertion or deletion anchored there).
fn after_its_first_base(r: &Record) -> bool {
    let first = r.ref_allele.as_bytes().first();
    !r.alts.is_empty()
        && r.alts
            .iter()
            .all(|a| a.len() != r.ref_allele.len() && !a.starts_with('<') && a.as_bytes().first() == first)
}

/// The smallest of a quality field across records, when every record reports it.
fn least(values: impl Iterator<Item = Option<String>>) -> Option<String> {
    let mut best: Option<(f64, String)> = None;
    for v in values {
        let v = v?;
        let x = v.parse::<f64>().ok()?;
        if best.as_ref().is_none_or(|(b, _)| x < *b) {
            best = Some((x, v));
        }
    }
    best.map(|(_, v)| v)
}

/// Records split from one multi-allelic site, read as that site: the records starting at the target's
/// position are merged into one record (REF the longest of their REFs, each ALT extended to it, as
/// `bcftools norm -m+` does), and for an SNV target the records there that are only indels after its base are
/// set aside. A haplotype with ALT alleles in two records becomes a no-call; FILTER and FT keep any failure,
/// and DP, MIN_DP and GQ the smallest value. Records starting elsewhere are kept as they are, so a deletion
/// spanning the target still makes it ambiguous. Nothing is merged when the records don't fit together
/// (different ploidy, a reference block or a symbolic ALT, REFs that are not prefixes of one another).
pub fn merge_split(records: &[Record], target: &Target<'_>) -> Vec<Record> {
    let (here, elsewhere): (Vec<&Record>, Vec<&Record>) = records.iter().partition(|r| r.pos == target.pos);
    let here: Vec<&Record> = if target.sequence {
        here
    } else {
        here.into_iter().filter(|r| !after_its_first_base(r)).collect()
    };
    let keep = |here: &[&Record]| elsewhere.iter().chain(here).map(|r| (*r).clone()).collect::<Vec<_>>();
    if here.len() < 2 {
        return keep(&here);
    }
    let longest = here
        .iter()
        .map(|r| &r.ref_allele)
        .max_by_key(|r| r.len())
        .expect("two records")
        .clone();
    let ploidy = here[0].gt.len();
    let fits = here.iter().all(|r| {
        longest.starts_with(r.ref_allele.as_str())
            && r.gt.len() == ploidy
            && !r.is_reference_block()
            && r.alts.iter().all(|a| !a.starts_with('<'))
    });
    if !fits {
        return keep(&here);
    }
    let mut alts: Vec<String> = Vec::new();
    let index: Vec<Vec<u32>> = here
        .iter()
        .map(|r| {
            let suffix = &longest[r.ref_allele.len()..];
            r.alts
                .iter()
                .map(|a| {
                    let extended = format!("{a}{suffix}");
                    let i = match alts.iter().position(|x| *x == extended) {
                        Some(i) => i,
                        None => {
                            alts.push(extended);
                            alts.len() - 1
                        }
                    };
                    i as u32 + 1
                })
                .collect()
        })
        .collect();
    let gt = (0..ploidy)
        .map(|h| {
            let mut allele = Some(0u32);
            for (r, map) in here.iter().zip(&index) {
                match r.gt[h] {
                    None => return None,
                    Some(0) => {}
                    Some(k) => {
                        let combined = *map.get(k as usize - 1)?;
                        allele = match allele {
                            Some(0) => Some(combined),
                            Some(prev) if prev == combined => Some(prev),
                            _ => return None,
                        };
                    }
                }
            }
            allele
        })
        .collect();
    let failing: Vec<String> = here
        .iter()
        .flat_map(|r| r.filters.iter())
        .filter(|f| *f != "PASS" && *f != ".")
        .cloned()
        .collect();
    let merged = Record {
        pos: target.pos,
        end: target.pos + longest.len() as u64 - 1,
        ref_allele: longest,
        alts,
        filters: if failing.is_empty() {
            here[0].filters.clone()
        } else {
            failing
        },
        gt,
        gt_phased: here.iter().all(|r| r.gt_phased),
        phase_set: here
            .iter()
            .all(|r| r.phase_set == here[0].phase_set)
            .then(|| here[0].phase_set.clone())
            .flatten(),
        ft: here
            .iter()
            .find_map(|r| r.ft.clone().filter(|f| f != "PASS" && f != "."))
            .or_else(|| here[0].ft.clone()),
        dp: least(here.iter().map(|r| r.dp.clone())),
        min_dp: least(here.iter().map(|r| r.min_dp.clone())),
        gq: least(here.iter().map(|r| r.gq.clone())),
    };
    let mut out: Vec<Record> = elsewhere.into_iter().cloned().collect();
    out.push(merged);
    out
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

/// The oriented target: the term's REF and ALT at `pos` (1-based). `alt: None` counts any non-reference
/// allele, for terms whose effect allele is the reference and whose other allele is not given.
#[derive(Clone, Debug, PartialEq)]
pub struct Target<'a> {
    pub pos: u64,
    pub ref_allele: &'a str,
    pub alt: Option<&'a str>,
    /// A normalized indel or multi-base target: records' alleles are compared after normalization, and
    /// several passing reference blocks may together cover it.
    pub sequence: bool,
    /// The target is on chrX or chrY, where `Policy::haploid_xy_as_homozygous` applies.
    pub sex_chromosome: bool,
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
    let doubled: Vec<Record>;
    let records = if policy.haploid_xy_as_homozygous && target.sex_chromosome && records.iter().any(|r| r.gt.len() == 1)
    {
        doubled = records
            .iter()
            .map(|r| {
                let mut r = r.clone();
                if r.gt.len() == 1 {
                    r.gt.push(r.gt[0]);
                }
                r
            })
            .collect();
        &doubled[..]
    } else {
        records
    };
    let merged: Vec<Record>;
    let records = if policy.merge_split_records && records.len() > 1 {
        merged = merge_split(records, target);
        &merged[..]
    } else {
        records
    };
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
        _ if target.sequence && blocks_cover(records, filters, target, policy) => {
            return Call {
                state: State::ObservedReference,
                alt_dosage: Some(0),
                phase_set: None,
                refcall_adapted: false,
            };
        }
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
    match (depth, gq) {
        (Some(depth), Some(gq)) if depth < policy.min_depth || gq < policy.min_gq => {
            return Call::state(State::LowQuality);
        }
        (Some(_), Some(_)) => {}
        (None, None) if policy.accept_missing_quality && !block => {}
        _ => return Call::state(State::QualityMissing),
    }
    let indices: Vec<u32> = record.gt.iter().map(|a| a.expect("no-calls handled above")).collect();
    let dosage = if block {
        if indices.iter().any(|&i| i != 0) {
            return Call::state(State::UnsupportedSymbolicGenotype);
        }
        // A sequence target spans several bases, so the block must also start at or before it.
        if record.end < target.pos + target.ref_allele.len() as u64 - 1 || record.pos > target.pos {
            return Call::state(State::IncompleteReferenceSpan);
        }
        0
    } else if target.sequence {
        match sequence_dosage(record, &indices, target, reference) {
            Ok(d) => d,
            Err(state) => return Call::state(state),
        }
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
                Some(alt) if target.alt.is_none_or(|t| alt == t) => dosage += 1,
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

/// Several reference blocks overlapping a sequence target count as homozygous reference only if every one is
/// a passing diploid `0/0` block and together they cover the target's REF span without a gap.
fn blocks_cover(records: &[Record], filters: &[&[String]], target: &Target<'_>, policy: &Policy) -> bool {
    let end = target.pos + target.ref_allele.len() as u64 - 1;
    let mut covered = target.pos;
    for (record, filters) in records.iter().zip(filters) {
        let passing = record.is_reference_block()
            && record.gt == [Some(0), Some(0)]
            && filters.iter().all(|f| f == "PASS" || f == ".")
            && record.ft.as_deref().is_none_or(|ft| ft == "PASS" || ft == ".")
            && number(&record.min_dp).is_some_and(|d| d >= policy.min_depth)
            && number(&record.gq).is_some_and(|q| q >= policy.min_gq);
        if !passing || record.pos > covered {
            return false;
        }
        covered = covered.max(record.end + 1);
    }
    covered > end
}

/// The target ALT dosage from one variant record, comparing each called allele with the target after
/// normalization. A REF call must span the whole target; any other called allele is not the target's.
fn sequence_dosage(
    record: &Record,
    indices: &[u32],
    target: &Target<'_>,
    reference: &impl Fn(u64, usize) -> Option<String>,
) -> Result<u8, State> {
    let base = |p: u64| reference(p, 1).and_then(|b| b.bytes().next());
    let wanted = crate::alleles::Variant {
        pos: target.pos,
        ref_allele: target.ref_allele.to_owned(),
        alt: target.alt.unwrap_or_default().to_owned(),
    };
    let end = target.pos + target.ref_allele.len() as u64 - 1;
    let mut dosage = 0u8;
    for &i in indices {
        if i == 0 {
            if record.pos > target.pos || record.end < end {
                return Err(State::UnsupportedAlleleRepresentation);
            }
            continue;
        }
        let alt = record.alts.get(i as usize - 1).ok_or(State::OtherCalledAllele)?;
        match crate::alleles::normalize(record.pos, &record.ref_allele, alt, base) {
            Some(v) if v == wanted => dosage += 1,
            _ => return Err(State::OtherCalledAllele),
        }
    }
    Ok(dosage)
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
        alt: Some("T"),
        sequence: false,
        sex_chromosome: false,
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

    /// A genotype-only record passes only with `accept_missing_quality`; a record reporting one of depth and
    /// GQ is still checked, and a reference block without quality never passes.
    #[test]
    fn genotype_only_calls_are_opt_in() {
        let accepting = Policy::for_gvcf(false, false, true);
        assert_ne!(accepting.id, Policy::default().id);
        let mut r = variant([Some(0), Some(1)]);
        (r.dp, r.gq) = (None, None);
        assert_eq!(state(&[r.clone()]), State::QualityMissing);
        let call = assess(std::slice::from_ref(&r), &TARGET, &accepting, reference);
        assert_eq!((call.state, call.alt_dosage), (State::ObservedVariant, Some(1)));
        r.gq = Some("5".into());
        assert_eq!(
            assess(std::slice::from_ref(&r), &TARGET, &accepting, reference).state,
            State::QualityMissing
        );
        r.dp = Some("30".into());
        assert_eq!(
            assess(std::slice::from_ref(&r), &TARGET, &accepting, reference).state,
            State::LowQuality
        );
        let mut b = block(1, 8);
        (b.min_dp, b.gq) = (None, None);
        assert_eq!(
            assess(&[b], &TARGET, &accepting, reference).state,
            State::QualityMissing
        );
    }

    fn split(pos: u64, r: &str, alt: &str, gt: [u32; 2]) -> Record {
        Record {
            pos,
            end: pos + r.len() as u64 - 1,
            ref_allele: r.into(),
            alts: vec![alt.into()],
            filters: vec!["PASS".into()],
            gt: gt.iter().map(|&a| Some(a)).collect(),
            gt_phased: true,
            phase_set: None,
            ft: None,
            dp: Some("30".into()),
            min_dp: None,
            gq: Some("40".into()),
        }
    }

    /// A multi-allelic site split into one record per ALT reads, with `merge_split_records`, as the site
    /// itself: the same call as the native multi-allelic record.
    #[test]
    fn split_records_read_as_one_site() {
        let merging = Policy::default().with_merge_split(true);
        let call = |records: &[Record], target: &Target<'_>, policy: &Policy| {
            let c = assess(records, target, policy, reference);
            (c.state, c.alt_dosage)
        };
        // G>T carried on one haplotype, G>A on neither.
        let site = [split(3, "G", "T", [0, 1]), split(3, "G", "A", [0, 0])];
        assert_eq!(
            call(&site, &TARGET, &Policy::default()).0,
            State::AmbiguousOverlappingRecords
        );
        assert_eq!(call(&site, &TARGET, &merging), (State::ObservedVariant, Some(1)));
        let native = Record {
            alts: vec!["T".into(), "A".into()],
            ..split(3, "G", "T", [0, 1])
        };
        assert_eq!(call(&site, &TARGET, &merging), call(&[native], &TARGET, &merging));
        // T on one haplotype, A on the other: another allele is called, as for the native record.
        let other = [split(3, "G", "T", [0, 1]), split(3, "G", "A", [1, 0])];
        assert_eq!(call(&other, &TARGET, &merging).0, State::OtherCalledAllele);
        // Two ALTs on one haplotype cannot both be true: no call.
        let conflict = [split(3, "G", "T", [1, 0]), split(3, "G", "A", [1, 0])];
        assert_eq!(call(&conflict, &TARGET, &merging).0, State::PartialNoCall);
        // A deletion anchored at the target's base leaves that base alone.
        let anchored = [split(3, "G", "T", [1, 1]), split(3, "GT", "G", [0, 1])];
        assert_eq!(call(&anchored, &TARGET, &merging), (State::ObservedVariant, Some(2)));
        // For the deletion itself (a sequence target) the records merge to REF GT, ALTs TT and G.
        let deletion = Target {
            pos: 3,
            ref_allele: "GT",
            alt: Some("G"),
            sequence: true,
            sex_chromosome: false,
        };
        // One haplotype carrying both the SNV and the deletion is a conflict in the merged record: no call.
        assert_eq!(call(&anchored, &deletion, &merging).0, State::PartialNoCall);
        // The SNV on one haplotype and the deletion on the other: another allele is called.
        let apart = [split(3, "G", "T", [1, 0]), split(3, "GT", "G", [0, 1])];
        assert_eq!(call(&apart, &deletion, &merging).0, State::OtherCalledAllele);
        let only_deletion = [split(3, "G", "T", [0, 0]), split(3, "GT", "G", [0, 1])];
        assert_eq!(
            call(&only_deletion, &deletion, &merging),
            (State::ObservedVariant, Some(1))
        );
        // A deletion starting before the target still spans it: ambiguous either way.
        let spanning = [split(2, "CG", "C", [0, 1]), split(3, "G", "T", [0, 1])];
        assert_eq!(call(&spanning, &TARGET, &merging).0, State::AmbiguousOverlappingRecords);
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
            alt: Some("T"),
            sequence: false,
            sex_chromosome: false,
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
    fn haploid_sex_chromosome_calls_only_with_the_policy() {
        let mut haploid = variant([Some(1), None]);
        haploid.gt = vec![Some(1)];
        let x = Target {
            sex_chromosome: true,
            ..TARGET
        };
        assert_eq!(
            assess(&[haploid.clone()], &x, &Policy::default(), reference).state,
            State::UnsupportedPloidy
        );
        let policy = Policy::for_gvcf(false, true, false);
        let call = assess(&[haploid.clone()], &x, &policy, reference);
        assert_eq!((call.state, call.alt_dosage), (State::ObservedVariant, Some(2)));
        // Autosomes stay diploid-only.
        assert_eq!(
            assess(&[haploid], &TARGET, &policy, reference).state,
            State::UnsupportedPloidy
        );
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
            alt: Some("G"),
            sequence: false,
            sex_chromosome: false,
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
    fn any_alt_target_counts_every_non_reference_allele() {
        let any = Target {
            pos: 3,
            ref_allele: "G",
            alt: None,
            sequence: false,
            sex_chromosome: false,
        };
        let mut r = variant([Some(1), Some(2)]);
        r.alts = vec!["T".into(), "C".into()];
        let call = assess(&[r], &any, &Policy::default(), reference);
        assert_eq!((call.state, call.alt_dosage), (State::ObservedVariant, Some(2)));
        let call = assess(&[block(1, 8)], &any, &Policy::default(), reference);
        assert_eq!((call.state, call.alt_dosage), (State::ObservedReference, Some(0)));
    }

    #[test]
    fn sequence_targets_need_gapless_passing_blocks() {
        let target = Target {
            pos: 1,
            ref_allele: "ACG",
            alt: Some("A"),
            sequence: true,
            sex_chromosome: false,
        };
        let covered = assess(&[block(1, 1), block(2, 5)], &target, &Policy::default(), reference);
        assert_eq!((covered.state, covered.alt_dosage), (State::ObservedReference, Some(0)));
        let gap = assess(&[block(1, 1), block(3, 5)], &target, &Policy::default(), reference);
        assert_eq!(gap.state, State::AmbiguousOverlappingRecords);
        let mut low = block(2, 5);
        low.gq = Some("5".into());
        let failing = assess(&[block(1, 1), low], &target, &Policy::default(), reference);
        assert_eq!(failing.state, State::AmbiguousOverlappingRecords);
        // A single block must also start at or before the target.
        assert_eq!(
            assess(&[block(2, 5)], &target, &Policy::default(), reference).state,
            State::IncompleteReferenceSpan
        );
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
