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
    pub seq: String,
    pub shift_left: usize,
    pub shift_right: usize,
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
    let mut best = Tsd { len: 0, seq: String::new(), shift_left: 0, shift_right: 0 };
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
                    };
                    break;
                }
            }
        }
    }
    best
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

/// Primer binding site candidate: offset of the first TGG (the complement of
/// a tRNA's CCA end) within 20 bp after the left LTR, or -1.
pub fn pbs(text: &[u8], left_ltr_end: usize, seq_end: usize) -> i64 {
    let lo = left_ltr_end + 1;
    let hi = (left_ltr_end + 20).min(seq_end);
    if hi < lo + 2 {
        return -1;
    }
    (lo..=hi - 2)
        .find(|&i| &text[i..i + 3] == b"TGG")
        .map_or(-1, |i| (i - lo) as i64)
}

/// Count of N in the element.
pub fn n_count(text: &[u8], start: usize, end: usize) -> usize {
    text[start..=end].iter().filter(|&&c| c == b'N').count()
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn motif_and_ppt_and_pbs() {
        //            0123456789012345678901234567890123
        let text = b"AAAATGCCCTGGAAAAAAAAAAAAAAAGGGCAAAA";
        assert_eq!(motif_at(text, 4, 31), "TGCA");
        let (n, off) = ppt(text, 32, 0);
        assert!(n >= 13);
        assert!(off < 0);
        assert_eq!(pbs(text, 8, text.len() - 1), 0);
        assert_eq!(n_count(text, 0, text.len() - 1), 0);
    }
}
