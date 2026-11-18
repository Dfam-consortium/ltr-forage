/// Target Site Duplication (TSD) detection and palindromic motif checking.
///
/// After X-drop boundary extension, LTRHarvest optionally searches for:
///   1. Short direct repeats (TSDs) flanking the element at the outer LTR
///      borders (5' of left LTR and 3' of right LTR).
///   2. A TG...CA palindromic motif immediately inside the TSDs.
///   3. The complementary motif (CA/TG) at the inner LTR borders (3' of
///      left LTR and 5' of right LTR).
///
/// The search window ("vicinity") is a configurable number of bp around each
/// boundary; the default in GenomeTools is 60 bp.

use crate::ltr::LtrCandidate;
use crate::xdrop::same_base;

// --------------------------------------------------
/// Motif configuration (default: TG...CA).
#[derive(Clone, Debug)]
pub struct MotifConfig {
    /// First char at the 5' end of each LTR (default 'T')
    pub firstleft: u8,
    /// Second char at the 5' end of each LTR (default 'G')
    pub secondleft: u8,
    /// First char at the 3' end of each LTR (default 'C')
    pub firstright: u8,
    /// Second char at the 3' end of each LTR (default 'A')
    pub secondright: u8,
    /// Maximum allowed mismatches per side (default 0)
    pub allowed_mismatches: u32,
}

impl Default for MotifConfig {
    fn default() -> Self {
        MotifConfig {
            firstleft: b'T',
            secondleft: b'G',
            firstright: b'C',
            secondright: b'A',
            allowed_mismatches: 0,
        }
    }
}

// --------------------------------------------------
/// Result of TSD + motif post-processing for one LTR candidate.
#[derive(Clone, Debug)]
pub struct TsdMotifInfo {
    /// A TSD was found at the outer borders.
    pub has_tsd: bool,
    /// TSD length (0 when has_tsd is false).
    pub tsd_len: usize,
    /// Outer motif (TG/CA at 5'/3' of the LTRs, near the TSD) was found.
    pub motif_near_tsd: bool,
    /// Inner motif (CA/TG at 3'/5' inner LTR borders) was found.
    pub motif_far_tsd: bool,
    /// Adjusted left LTR start (may differ from original when TSD/motif found).
    pub left_ltr_start: usize,
    /// Adjusted left LTR end (may differ when inner motif found).
    pub left_ltr_end: usize,
    /// Adjusted right LTR start (may differ when inner motif found).
    pub right_ltr_start: usize,
    /// Adjusted right LTR end (may differ from original when TSD/motif found).
    pub right_ltr_end: usize,
}

// --------------------------------------------------
/// Find all left-maximal exact matches of length ≥ `min_len` between the
/// `left` and `right` slices, mirroring GenomeTools `gt_sarrquerysubstringmatch`.
///
/// GT builds a suffix array of the left window (database) and searches each
/// position of the right window (query) against it, iterating SA hits in rank
/// order.  That outer-loop order (right position ascending, then left SA rank
/// ascending) is what determines tie-breaking in
/// `searchforbestTSDandormotifatborders` when two TSDs have equal total
/// deviation.  Matching that order exactly is required for GT-compatible output.
///
/// Returns absolute `(left_pos, right_pos, match_len)` tuples in GT's emission
/// order: right position ascending (primary), left SA rank ascending (secondary).
fn find_exact_matches_sa(
    left: &[u8],
    left_offset: usize,
    right: &[u8],
    right_offset: usize,
    min_len: usize,
) -> Vec<(usize, usize, usize)> {
    let llen = left.len();
    let rlen = right.len();

    if llen < min_len || rlen < min_len {
        return vec![];
    }

    // Build SA of the left window (the "database" in GT terminology).
    let mut sa: Vec<usize> = (0..llen).collect();
    sa.sort_unstable_by(|&a, &b| left[a..].cmp(&left[b..]));

    let mut results: Vec<(usize, usize, usize)> = Vec::new();

    // Outer loop: right window position ascending — GT's querysubstring.currentoffset.
    // This is the primary key for tie-breaking.
    for rpos in 0..=(rlen - min_len) {
        let pattern = &right[rpos..rpos + min_len];
        // A TSD never contains a wildcard (GT: special characters match nothing).
        if !pattern.iter().all(|&c| same_base(c, c)) {
            continue;
        }

        // Binary search for SA interval [lo, hi) where left[sa[k]..] shares
        // at least min_len characters with right[rpos..].  Mirrors GT's
        // gt_mmsearch binary search for the LCP-interval.
        let lo = sa.partition_point(|&i| {
            left[i..].iter().take(min_len).lt(pattern.iter())
        });
        let hi = sa.partition_point(|&i| {
            !left[i..].iter().take(min_len).gt(pattern.iter())
        });

        // Iterate SA[lo..hi] in rank order — GT's mmsearchiterator_next.
        // SA rank within the interval is the secondary tie-breaking key.
        for k in lo..hi {
            let lpos = sa[k];

            // Left-maximality check — GT's gt_mmsearch_isleftmaximal.
            if lpos > 0 && rpos > 0 && same_base(left[lpos - 1], right[rpos - 1]) {
                continue;
            }

            // Extend match past min_len — GT's gt_mmsearch_extendright.
            let mut ext = 0;
            while lpos + min_len + ext < llen
                && rpos + min_len + ext < rlen
                && same_base(left[lpos + min_len + ext], right[rpos + min_len + ext])
            {
                ext += 1;
            }

            results.push((left_offset + lpos, right_offset + rpos, min_len + ext));
        }
    }

    results
}

// --------------------------------------------------
/// Known LTR terminal 4-mers beyond TGCA, from LTR_retriever's curated list.
/// Each entry is [first2 of left LTR, last2 of right LTR].
const KNOWN_MOTIFS: &[[u8; 4]] = &[
    [b'T', b'G', b'C', b'T'],
    [b'T', b'A', b'C', b'A'],
    [b'T', b'A', b'C', b'T'],
    [b'T', b'G', b'G', b'A'],
    [b'T', b'A', b'T', b'A'],
    [b'T', b'G', b'T', b'A'],
];

/// Tiered preference for a TSD candidate, mirroring LTR_retriever's six-tier
/// ranking.  Lower tier = higher preference.
///
/// Tier 1: 5 bp TSD + TGCA motif (canonical)
/// Tier 2: 5 bp TSD + known non-TGCA motif
/// Tier 3: <5 bp TSD + TGCA
/// Tier 4: 5 bp TSD + unknown motif
/// Tier 5: <5 bp TSD + non-TGCA / unknown
/// Tier 6: >5 bp TSD (unusual; likely over-extended)
fn tsd_tier(left2: [u8; 2], right2: [u8; 2], tsd_len: usize) -> u8 {
    let is_tgca = left2 == [b'T', b'G'] && right2 == [b'C', b'A'];
    let is_known = is_tgca || KNOWN_MOTIFS.iter().any(|m| {
        left2[0] == m[0] && left2[1] == m[1] && right2[0] == m[2] && right2[1] == m[3]
    });
    match (tsd_len.cmp(&5), is_tgca, is_known) {
        (std::cmp::Ordering::Equal,   true,  _)     => 1,
        (std::cmp::Ordering::Equal,   false, true)  => 2,
        (std::cmp::Ordering::Less,    true,  _)     => 3,
        (std::cmp::Ordering::Equal,   false, false) => 4,
        (std::cmp::Ordering::Less,    false, _)     => 5,
        _                                           => 6,
    }
}

// --------------------------------------------------
/// Search for TSDs (and optionally the outer TG/CA motif) at the 5' border
/// of the left LTR and the 3' border of the right LTR.
///
/// Mirrors GenomeTools `searchforbestTSDandormotifatborders`.
///
/// On success returns adjusted `(left_ltr_start, right_ltr_end, tsd_len,
/// motif_near_tsd)`.  The boundaries may be unchanged if no suitable TSD is
/// found.
fn search_tsd_at_outer_borders(
    text: &[u8],
    cand: &LtrCandidate,
    seq_start: usize,
    seq_end: usize,   // last base of this sequence
    min_tsd: usize,
    max_tsd: usize,
    vicinity: usize,
    motif: &MotifConfig,
    check_motif: bool,
    tiered: bool,
) -> (usize, usize, usize, bool) // (new_left_start, new_right_end, tsd_len, motif_ok)
{
    let old_left = cand.left_ltr_start;
    let old_right = cand.right_ltr_end;

    // ---- Define left search window ----
    let start_left = if old_left >= seq_start + vicinity {
        old_left - vicinity
    } else {
        seq_start
    };
    // Do not extend into the LTR body (keep 2 bp gap for possible motif).
    let end_left = {
        let e = old_left + vicinity;
        if e + 2 > cand.left_ltr_end { cand.left_ltr_end.saturating_sub(2) } else { e }
    };

    // ---- Define right search window ----
    // Do not extend into the LTR body.
    let start_right = {
        let s = old_right.saturating_sub(vicinity);
        if s < cand.right_ltr_start + 2 { cand.right_ltr_start + 2 } else { s }
    };
    let end_right = if old_right + vicinity <= seq_end {
        old_right + vicinity
    } else {
        seq_end
    };

    if end_left < start_left || end_right < start_right {
        return (old_left, old_right, 0, false);
    }

    let left_slice  = &text[start_left..=end_left];
    let right_slice = &text[start_right..=end_right];

    // ---- Find all exact matches ≥ min_tsd between the two windows ----
    let repeats = find_exact_matches_sa(left_slice, start_left, right_slice, start_right, min_tsd);

    if repeats.is_empty() {
        return (old_left, old_right, 0, false);
    }

    // ---- Process each repeat to find the best TSD (± motif) ----
    let mut best_left     = old_left;
    let mut best_right    = old_right;
    let mut best_tsd_len  = 0usize;
    let mut best_motif_ok = false;
    let mut best_tier     = u8::MAX;
    let mut best_dev      = usize::MAX;
    let mut has_hit       = false;

    for (pos1, pos2, rep_len) in &repeats {
        let (pos1, pos2, rep_len) = (*pos1, *pos2, *rep_len);

        // motifpos1 = first position after left TSD = tentative new leftLTR_5.
        let motifpos1 = pos1 + rep_len;
        // motifpos2 = two positions before right TSD = tentative rightLTR_3 - 1.
        // (motifpos2 and motifpos2+1 are the two CA chars; rightLTR_3 = motifpos2+1)
        if pos2 < 2 { continue; }
        let motifpos2 = pos2 - 2;

        // Trim `back` off the right end of left TSD and `forward` off the left
        // end of right TSD.  tsd_len = rep_len - back - forward ≥ min_tsd.
        let max_trim = rep_len.saturating_sub(min_tsd);

        for back in 0..=max_trim {
            let fwd_limit = max_trim - back;
            for forward in 0..=fwd_limit {
                let tsd_len = rep_len - back - forward;
                if tsd_len > max_tsd {
                    continue;
                }

                // Candidate LTR boundaries after this (back, forward) pair.
                let new_left  = motifpos1.saturating_sub(back);
                // rightLTR_3 = motifpos2 + 1 + forward
                if motifpos2 + 1 + forward > seq_end {
                    continue;
                }
                let new_right = motifpos2 + 1 + forward;

                // Sanity: new boundaries must not flip or overflow.
                if new_left > cand.left_ltr_end || new_right < cand.right_ltr_start {
                    continue;
                }

                // Bounds check needed for both motif reading and tier computation.
                if new_left + 1 >= text.len() { continue; }
                if new_right < 1 || new_right >= text.len() { continue; }

                // Read the 4-mer at the TSD-LTR junction regardless of --motif.
                // [first 2 bp of left LTR] + [last 2 bp of right LTR].
                let left2  = [text[new_left], text[new_left + 1]];
                let right2 = [text[new_right - 1], text[new_right]];

                let motif_ok = if check_motif {
                    let mm_left = {
                        let mut mm = 0u32;
                        if text[new_left]     != motif.firstleft  { mm += 1; }
                        if text[new_left + 1] != motif.secondleft { mm += 1; }
                        mm
                    };
                    let mm_right = {
                        let mut mm = 0u32;
                        if text[new_right - 1] != motif.firstright  { mm += 1; }
                        if text[new_right]      != motif.secondright { mm += 1; }
                        mm
                    };
                    mm_left  <= motif.allowed_mismatches &&
                    mm_right <= motif.allowed_mismatches
                } else {
                    false
                };

                // When --motif is requested, only consider TGCA-passing candidates.
                if check_motif && !motif_ok {
                    continue;
                }

                // Rank by (tier, deviation): biologically preferred motif+length
                // combinations win over candidates that merely have lower deviation.
                let tier = tsd_tier(left2, right2, tsd_len);
                let dev = {
                    let d1 = if new_left  > old_left  { new_left  - old_left  }
                             else                      { old_left  - new_left  };
                    let d2 = if new_right > old_right  { new_right - old_right }
                             else                      { old_right - new_right };
                    d1 + d2
                };

                // GT ranks by deviation only; tiered mode (--tiered-tsd) prefers
                // biologically preferred motif+length combinations first.
                let is_better = !has_hit || if tiered {
                    (tier, dev) < (best_tier, best_dev)
                } else {
                    dev < best_dev
                };
                if is_better {
                    has_hit       = true;
                    best_left     = new_left;
                    best_right    = new_right;
                    best_tsd_len  = tsd_len;
                    best_motif_ok = motif_ok;
                    best_tier     = tier;
                    best_dev      = dev;
                }
            }
        }
    }

    (best_left, best_right, best_tsd_len, best_motif_ok)
}

// --------------------------------------------------
/// Search for the outer TG/CA motif only (no TSD), scanning the vicinity
/// of the outer LTR borders.
///
/// Mirrors GenomeTools `searchformotifonlyborders`.
///
/// Returns `(new_left_start, new_right_end, motif_found)`.
fn search_motif_only_borders(
    text: &[u8],
    cand: &LtrCandidate,
    seq_start: usize,
    seq_end: usize,
    vicinity: usize,
    motif: &MotifConfig,
) -> (usize, usize, bool) {
    let old_left  = cand.left_ltr_start;
    let old_right = cand.right_ltr_end;

    // Left search window for TG (5' of left LTR).
    let start_left = if old_left >= seq_start + vicinity {
        old_left - vicinity
    } else {
        seq_start
    };
    let end_left = {
        let e = old_left + vicinity;
        // endleftLTR must be ≥ 2 bp before leftLTR_3 (for possible motif space).
        if e + 2 > cand.left_ltr_end { cand.left_ltr_end.saturating_sub(2) } else { e }
    };

    // Right search window for CA (3' of right LTR).
    let start_right = {
        let s = old_right.saturating_sub(vicinity);
        if s < cand.right_ltr_start + 2 { cand.right_ltr_start + 2 } else { s }
    };
    let end_right = if old_right + vicinity <= seq_end {
        old_right + vicinity
    } else {
        seq_end
    };

    let mut new_left = old_left;
    let mut new_right = old_right;
    let mut motif1 = false;
    let mut motif2 = false;
    let mut best_dev_left  = usize::MAX;
    let mut best_dev_right = usize::MAX;

    // Search for TG around leftLTR_5.
    if end_left + 1 < text.len() {
        for idx in start_left..end_left {
            let mm = {
                let mut m = 0u32;
                if text[idx]     != motif.firstleft  { m += 1; }
                if text[idx + 1] != motif.secondleft { m += 1; }
                m
            };
            if mm <= motif.allowed_mismatches {
                let dev = if idx > old_left { idx - old_left } else { old_left - idx };
                if !motif1 || dev < best_dev_left {
                    motif1 = true;
                    new_left = idx;
                    best_dev_left = dev;
                }
            }
        }
    }

    // Search for CA around rightLTR_3.
    if end_right >= 1 && start_right + 1 < text.len() {
        for idx in (start_right + 1)..=end_right {
            if idx >= text.len() { break; }
            let mm = {
                let mut m = 0u32;
                if text[idx - 1] != motif.firstright  { m += 1; }
                if text[idx]     != motif.secondright { m += 1; }
                m
            };
            if mm <= motif.allowed_mismatches {
                let dev = if idx > old_right { idx - old_right } else { old_right - idx };
                if !motif2 || dev < best_dev_right {
                    motif2 = true;
                    new_right = idx;
                    best_dev_right = dev;
                }
            }
        }
    }

    (new_left, new_right, motif1 && motif2)
}

// --------------------------------------------------
/// Search for the inner CA/TG motif at the 3' border of the left LTR
/// and the 5' border of the right LTR.
///
/// Mirrors GenomeTools `searchformotifonlyinside`.
///
/// Returns `(new_left_end, new_right_start, motif_found)`.
fn search_motif_inside(
    text: &[u8],
    cand: &LtrCandidate,
    vicinity: usize,
    motif: &MotifConfig,
) -> (usize, usize, bool) {
    let old_left_end   = cand.left_ltr_end;
    let old_right_start = cand.right_ltr_start;

    // Search vicinity of 3' border of left LTR.
    // Do not go below leftLTR_5 + 2 or above rightLTR_5 - 1.
    let start_left = {
        let s = old_left_end.saturating_sub(vicinity - 1);
        if s < cand.left_ltr_start + 2 { cand.left_ltr_start + 2 } else { s }
    };
    let end_left = {
        let e = old_left_end + vicinity - 1;
        if e >= old_right_start { old_right_start - 1 } else { e }
    };

    // Search vicinity of 5' border of right LTR.
    // Do not go below leftLTR_3 + 1 or above rightLTR_3 - 2.
    let start_right = {
        let s = old_right_start.saturating_sub(vicinity - 1);
        if s <= old_left_end { old_left_end + 1 } else { s }
    };
    let end_right = {
        let e = old_right_start + vicinity - 1;
        if e + 2 > cand.right_ltr_end { cand.right_ltr_end.saturating_sub(2) } else { e }
    };

    let mut new_left_end    = old_left_end;
    let mut new_right_start = old_right_start;
    let mut motif1 = false;
    let mut motif2 = false;
    let mut best_dev_left  = usize::MAX;
    let mut best_dev_right = usize::MAX;

    // Search for CA (firstright/secondright) at 3' of left LTR.
    // idx is the position of secondright (the 'A').
    if end_left < text.len() && start_left + 1 <= end_left {
        for idx in (start_left + 1)..=end_left {
            if idx >= text.len() { break; }
            let mm = {
                let mut m = 0u32;
                if text[idx - 1] != motif.firstright  { m += 1; }
                if text[idx]     != motif.secondright { m += 1; }
                m
            };
            if mm <= motif.allowed_mismatches {
                let dev = if idx > old_left_end { idx - old_left_end } else { old_left_end - idx };
                if !motif1 || dev < best_dev_left {
                    motif1 = true;
                    new_left_end = idx;
                    best_dev_left = dev;
                }
            }
        }
    }

    // Search for TG (firstleft/secondleft) at 5' of right LTR.
    // idx is the position of firstleft (the 'T').
    if end_right + 1 < text.len() && start_right < end_right {
        for idx in start_right..end_right {
            if idx + 1 >= text.len() { break; }
            let mm = {
                let mut m = 0u32;
                if text[idx]     != motif.firstleft  { m += 1; }
                if text[idx + 1] != motif.secondleft { m += 1; }
                m
            };
            if mm <= motif.allowed_mismatches {
                let dev = if idx > old_right_start { idx - old_right_start } else { old_right_start - idx };
                if !motif2 || dev < best_dev_right {
                    motif2 = true;
                    new_right_start = idx;
                    best_dev_right = dev;
                }
            }
        }
    }

    (new_left_end, new_right_start, motif1 && motif2)
}

// --------------------------------------------------
/// Main entry point: run TSD and/or motif detection on an LTR candidate.
///
/// `min_tsd == 0` disables TSD search (GenomeTools: min_tsd == 1 means no
/// TSD, because min_tsd > 1U is the condition for searching).
///
/// `check_motif` enables the motif search (both outer and inner borders).
pub fn annotate_tsd_and_motif(
    cand: &LtrCandidate,
    text: &[u8],
    sequence_starts: &[usize],
    min_tsd: usize,
    max_tsd: usize,
    vicinity: usize,
    motif: &MotifConfig,
    check_motif: bool,
    tiered_tsd: bool,
) -> TsdMotifInfo {
    let seq_start = sequence_starts[cand.seq_num];
    // the sequence's last base (GT's seqendpos): two before the next sequence's start, since
    // one separator byte lies between, or two before the end of the text (the terminator)
    let seq_end = sequence_starts
        .get(cand.seq_num + 1)
        .copied()
        .unwrap_or(text.len())
        .saturating_sub(2);

    let mut info = TsdMotifInfo {
        has_tsd:          false,
        tsd_len:          0,
        motif_near_tsd:   false,
        motif_far_tsd:    false,
        left_ltr_start:   cand.left_ltr_start,
        left_ltr_end:     cand.left_ltr_end,
        right_ltr_start:  cand.right_ltr_start,
        right_ltr_end:    cand.right_ltr_end,
    };

    // ---- Outer borders: TSD and/or motif ----
    if min_tsd >= 2 {
        let (new_left, new_right, tsd_len, motif_ok) = search_tsd_at_outer_borders(
            text, cand, seq_start, seq_end,
            min_tsd, max_tsd, vicinity, motif, check_motif, tiered_tsd,
        );
        if tsd_len > 0 {
            info.has_tsd        = true;
            info.tsd_len        = tsd_len;
            info.motif_near_tsd = motif_ok;
            info.left_ltr_start = new_left;
            info.right_ltr_end  = new_right;
        }
    } else if check_motif {
        // No TSD search, motif only.
        // Create a temporary candidate with current boundaries for the search.
        let (new_left, new_right, motif_ok) = search_motif_only_borders(
            text, cand, seq_start, seq_end, vicinity, motif,
        );
        info.motif_near_tsd = motif_ok;
        if motif_ok {
            info.left_ltr_start = new_left;
            info.right_ltr_end  = new_right;
        }
    }

    // ---- Inner borders: motif only (if enabled) ----
    if check_motif && motif.allowed_mismatches < 4 {
        // Build a temporary candidate with potentially adjusted outer borders.
        let adjusted = LtrCandidate {
            seq_num:         cand.seq_num,
            left_ltr_start:  info.left_ltr_start,
            left_ltr_end:    info.left_ltr_end,
            right_ltr_start: info.right_ltr_start,
            right_ltr_end:   info.right_ltr_end,
            similarity:      cand.similarity,
        };
        let (new_le, new_rs, motif_ok) = search_motif_inside(
            text, &adjusted, vicinity, motif,
        );
        info.motif_far_tsd   = motif_ok;
        if motif_ok {
            info.left_ltr_end    = new_le;
            info.right_ltr_start = new_rs;
        }
    }

    info
}
