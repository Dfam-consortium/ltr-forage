//! Per-candidate structural features, reported rather than used to filter.
//!
//! Each function reads the indexed text directly.  Coordinates are absolute
//! offsets into the concatenated text, inclusive; `seq_start` and `seq_end`
//! bound the sequence the candidate lies in so that no window crosses into a
//! neighbouring sequence.  The text is upper-case ACGTN (sufr upper-cases DNA
//! at index time), so byte comparisons are enough.

/// The longest exact direct repeat that flanks the element: it must end just
/// before the left LTR and begin just after the right LTR, allowing up to
/// `slack` bases of play at each end because the boundaries are seldom exact
/// to the base.  Returns the length, the repeat, and the play used on each
/// side; length 0 when nothing of at least `min_len` is found.
pub struct Tsd {
    pub len: usize,
    /// The left copy (the right copy differs from it when `mismatches` is 1).
    pub seq: String,
    pub shift_left: usize,
    pub shift_right: usize,
    /// 0 for an exact repeat; 1 for the fallback in `tsd_or_mm1`.
    pub mismatches: usize,
}

#[allow(clippy::too_many_arguments)]
pub fn tsd_at(
    text: &[u8],
    left_ltr_start: usize,
    right_ltr_end: usize,
    seq_start: usize,
    seq_end: usize,
    min_len: usize,
    max_len: usize,
    slack: usize,
) -> Tsd {
    let mut best = Tsd { len: 0, seq: String::new(), shift_left: 0, shift_right: 0, mismatches: 0 };
    for dl in 0..=slack {
        for dr in 0..=slack {
            // Left repeat occupies [L - dl - k, L - dl); right repeat [R + 1 + dr, R + 1 + dr + k).
            for k in (min_len..=max_len).rev() {
                if k <= best.len {
                    break;
                }
                if left_ltr_start < dl + k || left_ltr_start - dl - k < seq_start {
                    continue;
                }
                let ls = left_ltr_start - dl - k;
                let rs = right_ltr_end + 1 + dr;
                if rs + k - 1 > seq_end {
                    continue;
                }
                let a = &text[ls..ls + k];
                let b = &text[rs..rs + k];
                if a == b && a.iter().all(|&c| c != b'N') {
                    best = Tsd {
                        len: k,
                        seq: String::from_utf8_lossy(a).into_owned(),
                        shift_left: dl,
                        shift_right: dr,
                        mismatches: 0,
                    };
                    break;
                }
            }
        }
    }
    best
}

/// `tsd_at`, and when it finds nothing, a 6 or 5 bp repeat with one
/// mismatch lying directly against both LTR ends (no play).  A TSD picks up
/// point mutations after insertion like any other sequence: in lotus, 20% of
/// TG..CA elements with no exact TSD have one at the LTR ends, against 1% in
/// flanks shifted 37 bp out.  Allowing play as well would raise that
/// background several-fold, so the fallback is fixed in place.
#[allow(clippy::too_many_arguments)]
pub fn tsd_or_mm1(
    text: &[u8],
    left_ltr_start: usize,
    right_ltr_end: usize,
    seq_start: usize,
    seq_end: usize,
    min_len: usize,
    max_len: usize,
    slack: usize,
) -> Tsd {
    let t = tsd_at(text, left_ltr_start, right_ltr_end, seq_start, seq_end, min_len, max_len, slack);
    if t.len > 0 {
        return t;
    }
    let acgt = |c: &u8| matches!(c, b'A' | b'C' | b'G' | b'T');
    for k in [6usize, 5] {
        if left_ltr_start < seq_start + k || right_ltr_end + k > seq_end {
            continue;
        }
        let a = &text[left_ltr_start - k..left_ltr_start];
        let b = &text[right_ltr_end + 1..right_ltr_end + 1 + k];
        if !a.iter().chain(b).all(acgt) {
            continue;
        }
        if a.iter().zip(b).filter(|(x, y)| x != y).count() == 1 {
            return Tsd { len: k, seq: String::from_utf8_lossy(a).into_owned(), shift_left: 0, shift_right: 0, mismatches: 1 };
        }
    }
    t
}

/// The first two bases of the left LTR and the last two of the right LTR,
/// as a 4-character string (TGCA for the canonical termini).
pub fn motif_at(text: &[u8], left_ltr_start: usize, right_ltr_end: usize) -> String {
    let mut s = String::with_capacity(4);
    for i in [left_ltr_start, left_ltr_start + 1, right_ltr_end - 1, right_ltr_end] {
        s.push(*text.get(i).unwrap_or(&b'N') as char);
    }
    s
}

/// Both LTRs end in a poly-A tail (at least 12 A in their last 15 bp), or
/// both start with poly-T (the same on the other strand).  LTRs end in CA;
/// a pair of poly-A tails is two copies of a non-LTR retroposon's 3' end.
pub fn polya_ltrs(text: &[u8], ltrs: (usize, usize, usize, usize)) -> bool {
    let (l1s, l1e, l2s, l2e) = ltrs;
    let count = |a: usize, b: usize, base: u8| -> usize {
        (a..=b).filter(|&i| text.get(i).map_or(false, |c| c.to_ascii_uppercase() == base)).count()
    };
    if l1e < l1s + 14 || l2e < l2s + 14 {
        return false;
    }
    let tail = count(l1e - 14, l1e, b'A') >= 12 && count(l2e - 14, l2e, b'A') >= 12;
    let head = count(l1s, l1s + 14, b'T') >= 12 && count(l2s, l2s + 14, b'T') >= 12;
    tail || head
}

/// LTR_retriever's reading of the terminal motif: `TGCA` for the canonical
/// termini, `known` for the other dinucleotide pairs its rules accept
/// (TGCT, TACA, TACT, TGGA, TGGT, TATA, TGTA, TGCC), `other` for the rest.
pub fn motif_class(motif: &str) -> &'static str {
    match motif {
        "TGCA" => "TGCA",
        "TGCT" | "TACA" | "TACT" | "TGGA" | "TGGT" | "TATA" | "TGTA" | "TGCC" => "known",
        _ => "other",
    }
}

/// Polypurine tract candidate: the most purine-rich 15 bp window in the 30
/// bp before the right LTR.  Returns the purine count in that window (out of
/// 15) and the window's start relative to the right LTR start (negative).
pub fn ppt(text: &[u8], right_ltr_start: usize, seq_start: usize) -> (usize, i64) {
    let lo = right_ltr_start.saturating_sub(30).max(seq_start);
    if right_ltr_start < lo + 15 {
        return (0, 0);
    }
    let mut best = (0usize, 0i64);
    for w in lo..=right_ltr_start - 15 {
        let n = text[w..w + 15].iter().filter(|&&c| c == b'A' || c == b'G').count();
        if n > best.0 {
            best = (n, w as i64 - right_ltr_start as i64);
        }
    }
    best
}

/// Count of N in the element.
pub fn n_count(text: &[u8], start: usize, end: usize) -> usize {
    text[start..=end].iter().filter(|&&c| c == b'N').count()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn motif_classes() {
        assert_eq!(motif_class("TGCA"), "TGCA");
        assert_eq!(motif_class("TACA"), "known");
        assert_eq!(motif_class("TGAA"), "other");
        assert_eq!(motif_class("NNNN"), "other");
    }

    #[test]
    fn tsd_found_with_slack() {
        //            0123456789012345678901234567890123456789
        let text = b"GGGGGACGTATGXXXXXXXXXXXXXXXXCAGACGTAGGGGG";
        // Element 10..29 (TG ... CA); the flanking repeat is GACGTA at 4..9
        // and 30..35 (the G of the poly-G run extends ACGTA to six bases).
        let t = tsd_at(text, 10, 29, 0, text.len() - 1, 4, 20, 1);
        assert_eq!(t.len, 6);
        assert_eq!(t.seq, "GACGTA");
        assert_eq!((t.shift_left, t.shift_right), (0, 0));
        // One base off on the left: still found with slack 1.
        let t = tsd_at(text, 11, 29, 0, text.len() - 1, 4, 20, 1);
        assert_eq!(t.len, 6);
        assert_eq!(t.shift_left, 1);
        let t = tsd_at(text, 11, 29, 0, text.len() - 1, 4, 20, 0);
        assert_eq!(t.len, 0);
    }

    #[test]
    fn tsd_one_mismatch_fallback() {
        //            0123456789012345678901234567890
        let text = b"CCCAGGATTGXXXXXXXXXXCAAAGATCCC";
        // Element 8..21 (TG ... CA); flanks AGGAT / AAGAT differ at one base.
        let t = tsd_at(text, 8, 21, 0, text.len() - 1, 4, 20, 1);
        assert_eq!(t.len, 0);
        let t = tsd_or_mm1(text, 8, 21, 0, text.len() - 1, 4, 20, 1);
        assert_eq!((t.len, t.mismatches, t.seq.as_str()), (5, 1, "AGGAT"));
        // Two mismatches: nothing.
        let text = b"CCCCTGATTGXXXXXXXXXXCAAAGATCCC";
        let t = tsd_or_mm1(text, 8, 21, 0, text.len() - 1, 4, 20, 1);
        assert_eq!(t.len, 0);
        // An exact TSD wins and reports no mismatch.
        let text = b"CCCAAGATTGXXXXXXXXXXCAAAGATCCC";
        let t = tsd_or_mm1(text, 8, 21, 0, text.len() - 1, 4, 20, 1);
        assert_eq!((t.len, t.mismatches), (5, 0));
    }

    #[test]
    fn motif_and_ppt() {
        //            0123456789012345678901234567890123
        let text = b"AAAATGCCCTGGAAAAAAAAAAAAAAAGGGCAAAA";
        assert_eq!(motif_at(text, 4, 31), "TGCA");
        let (n, off) = ppt(text, 32, 0);
        assert!(n >= 13);
        assert!(off < 0);
        assert_eq!(n_count(text, 0, text.len() - 1), 0);
    }
}
