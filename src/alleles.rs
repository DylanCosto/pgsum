//! Allele normalization: the one representation of a variant that both scoring files and gVCFs are compared in.
//!
//! A variant `(pos, REF, ALT)` is normalized as in Tan, Abecasis & Kang (2015): common trailing bases are
//! removed and the variant is shifted left through the reference while the alleles allow it, then common
//! leading bases are removed while both alleles keep at least one base. The result is left-aligned and
//! parsimonious; two descriptions of the same change always normalize to the same triple.

/// A normalized variant.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Variant {
    /// 1-based position of the first REF base.
    pub pos: u64,
    pub ref_allele: String,
    pub alt: String,
}

/// Normalize `(pos, ref, alt)`. `base(p)` is the upper-case reference base at 1-based position `p`, used when
/// the variant shifts left. `None` if the alleles are equal, not ACGT, or the shift runs off the contig.
pub fn normalize(pos: u64, ref_allele: &str, alt: &str, base: impl Fn(u64) -> Option<u8>) -> Option<Variant> {
    let acgt = |s: &str| !s.is_empty() && s.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'));
    if !acgt(ref_allele) || !acgt(alt) || ref_allele == alt {
        return None;
    }
    let (mut pos, mut r, mut a) = (pos, ref_allele.as_bytes().to_vec(), alt.as_bytes().to_vec());
    loop {
        let mut changed = false;
        if !r.is_empty() && !a.is_empty() && r.last() == a.last() {
            r.pop();
            a.pop();
            changed = true;
        }
        if r.is_empty() || a.is_empty() {
            pos = pos.checked_sub(1).filter(|&p| p >= 1)?;
            let b = base(pos)?;
            r.insert(0, b);
            a.insert(0, b);
            changed = true;
        }
        if !changed {
            break;
        }
    }
    while r.len() >= 2 && a.len() >= 2 && r[0] == a[0] {
        r.remove(0);
        a.remove(0);
        pos += 1;
    }
    Some(Variant {
        pos,
        ref_allele: String::from_utf8(r).expect("ACGT"),
        alt: String::from_utf8(a).expect("ACGT"),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    //                     1234567890
    const SEQ: &[u8] = b"GCATTTTGAC";

    fn n(pos: u64, r: &str, a: &str) -> Option<(u64, String, String)> {
        normalize(pos, r, a, |p| SEQ.get(p as usize - 1).copied()).map(|v| (v.pos, v.ref_allele, v.alt))
    }

    fn v(pos: u64, r: &str, a: &str) -> Option<(u64, String, String)> {
        Some((pos, r.to_owned(), a.to_owned()))
    }

    #[test]
    fn deletions_in_a_run_left_align() {
        // Deleting one T from TTTT, written at the end of the run, anchored at the preceding A.
        assert_eq!(n(6, "TT", "T"), v(3, "AT", "A"));
        assert_eq!(n(3, "AT", "A"), v(3, "AT", "A"));
        assert_eq!(n(3, "ATTTTG", "ATTTG"), v(3, "AT", "A"));
    }

    #[test]
    fn insertions_in_a_run_left_align() {
        assert_eq!(n(7, "T", "TT"), v(3, "A", "AT"));
        assert_eq!(n(5, "TT", "TTT"), v(3, "A", "AT"));
    }

    #[test]
    fn multi_allelic_record_alleles_normalize() {
        // A record REF=ATTTT with ALT=ATTT (one T deleted) is the same as REF=AT ALT=A.
        assert_eq!(n(3, "ATTTT", "ATTT"), v(3, "AT", "A"));
    }

    #[test]
    fn substitutions_are_trimmed() {
        assert_eq!(n(2, "CA", "CG"), v(3, "A", "G"));
        assert_eq!(n(3, "AT", "GC"), v(3, "AT", "GC"));
        assert_eq!(n(3, "A", "G"), v(3, "A", "G"));
    }

    #[test]
    fn invalid_inputs() {
        assert_eq!(n(3, "A", "A"), None);
        assert_eq!(n(3, "A", "N"), None);
        assert_eq!(n(1, "GC", "C"), None, "shifting left past position 1");
    }
}
