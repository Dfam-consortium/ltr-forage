/// LTR candidate seed filtering and X-drop extension into full LTR predictions.
///
/// Phase 1 (filter_seed): select maximal repeats that could be LTR seeds.
/// Phase 2 (extend_seed): X-drop extend each seed to find full LTR boundaries.
///   Length, distance and similarity filters run in main.rs after the TSD
///   search, matching GenomeTools' order.

use crate::xdrop::{XdropScores, same_base, xdrop_extend};

// --------------------------------------------------
/// A maximal repeat that passed all LTR distance/length filters.
#[derive(Debug, Clone)]
pub struct Seed {
    /// Absolute start of the left repeat occurrence in the concatenated text
    pub pos1: usize,
    /// Absolute start of the right repeat occurrence in the concatenated text
    pub pos2: usize,
    /// Length of the shared exact match
    pub len: usize,
    /// Index into the sequence_starts / sequence_names arrays
    pub seq_num: usize,
}

// --------------------------------------------------
/// A fully extended LTR retrotransposon candidate.
#[derive(Debug, Clone)]
pub struct LtrCandidate {
    /// Sequence index (contig)
    pub seq_num: usize,
    /// Left LTR: 5′ boundary (absolute, inclusive)
    pub left_ltr_start: usize,
    /// Left LTR: 3′ boundary (absolute, inclusive)
    pub left_ltr_end: usize,
    /// Right LTR: 5′ boundary (absolute, inclusive)
    pub right_ltr_start: usize,
    /// Right LTR: 3′ boundary (absolute, inclusive)
    pub right_ltr_end: usize,
    /// Percent identity between the two LTR copies
    pub similarity: f64,
}

// --------------------------------------------------
/// Return the sequence index that contains `pos`, using the sorted
/// `sequence_starts` slice (binary search).
pub fn seq_num_for(pos: usize, sequence_starts: &[usize]) -> usize {
    match sequence_starts.binary_search(&pos) {
        Ok(i) => i,
        Err(i) => i.saturating_sub(1),
    }
}

// --------------------------------------------------
/// Test whether a raw maximal repeat satisfies the LTR seed constraints.
pub fn filter_seed(
    pos1: usize,
    pos2: usize,
    len: usize,
    lmax: usize,
    dmin: usize,
    dmax: usize,
    sequence_starts: &[usize],
) -> Option<Seed> {
    if len > lmax {
        return None;
    }
    let dist = pos2 - pos1;
    if dist < dmin || dist > dmax {
        return None;
    }
    let sn1 = seq_num_for(pos1, sequence_starts);
    let sn2 = seq_num_for(pos2, sequence_starts);
    if sn1 != sn2 {
        return None;
    }
    if let Some(&next_start) = sequence_starts.get(sn1 + 1) {
        if pos2 + len > next_start {
            return None;
        }
    }
    Some(Seed { pos1, pos2, len, seq_num: sn1 })
}

// --------------------------------------------------
/// Compute unit-cost edit distance between two sequences using Ukkonen's
/// greedy front algorithm (O(n·d) time, O(d) space).
///
/// Returns the edit distance.  Used for similarity scoring after extension.
pub fn unit_edit_distance(u: &[u8], v: &[u8]) -> usize {
    let ulen = u.len() as i64;
    let vlen = v.len() as i64;

    // Two alternating front arrays: front[k] = furthest row on diagonal k.
    // Diagonal k: i - j = k, so j = i - k.
    let max_d = (ulen + vlen) as usize;
    let offset = max_d as i64; // shift so index = k + offset
    let width = 2 * max_d + 1;
    let mut prev = vec![-1i64; width];
    let mut curr = vec![-1i64; width];

    // Helper: extend matches from (i, j) on diagonal k.
    let extend = |i: i64, k: i64| -> i64 {
        let mut i = i;
        let j = i - k;
        while i < ulen && j >= 0 && (i - k) < vlen && same_base(u[i as usize], v[(i - k) as usize]) {
            i += 1;
        }
        i
    };

    // d = 0
    let i0 = extend(0, 0);
    prev[(0 + offset) as usize] = i0;
    if i0 == ulen && i0 - 0 == vlen {
        return 0;
    }

    for d in 1..=(max_d as i64) {
        for k in (-d..=d).step_by(1) {
            // Choose the best predecessor.
            let i_del = if k + 1 <= d {
                *prev.get((k + 1 + offset) as usize).unwrap_or(&-1)
            } else { -1 };
            let i_ins = if k - 1 >= -d {
                let v = *prev.get((k - 1 + offset) as usize).unwrap_or(&-1);
                if v >= 0 { v + 1 } else { -1 }
            } else { -1 };
            let i_sub = if k >= -(d-1) && k <= d-1 {
                let v = *prev.get((k + offset) as usize).unwrap_or(&-1);
                if v >= 0 { v + 1 } else { -1 }
            } else { -1 };

            let i_best = i_del.max(i_ins).max(i_sub);
            let i_ext = if i_best >= 0 { extend(i_best, k) } else { -1 };
            curr[(k + offset) as usize] = i_ext;

            if i_ext == ulen && i_ext - k == vlen {
                return d as usize;
            }
        }
        std::mem::swap(&mut prev, &mut curr);
        for v in curr.iter_mut() { *v = -1; }
    }

    max_d
}

// --------------------------------------------------
/// GenomeTools' `checklengthanddistanceconstraints`, applied to a set of
/// LTR boundaries.  Called after TSD/motif adjustment, as GT does, because
/// the TSD search can move the outer boundaries by up to `vicinity` bp.
pub fn constraints_ok(
    left_ltr_start: usize,
    left_ltr_end: usize,
    right_ltr_start: usize,
    right_ltr_end: usize,
    lmin: usize,
    lmax: usize,
    dmin: usize,
    dmax: usize,
) -> bool {
    if left_ltr_end >= right_ltr_start {
        return false;
    }
    let ulen = left_ltr_end - left_ltr_start + 1;
    let vlen = right_ltr_end - right_ltr_start + 1;
    let dist = right_ltr_start - left_ltr_start;
    ulen >= lmin && ulen <= lmax && vlen >= lmin && vlen <= lmax && dist >= dmin && dist <= dmax
}

// --------------------------------------------------
/// Extend a seed with X-drop alignment in both directions and return the
/// resulting LTR boundaries.
///
/// No length, distance or similarity filtering happens here.  GT applies
/// those after the TSD search has adjusted the outer boundaries, so the
/// caller runs `constraints_ok` and the similarity computation on the
/// post-TSD boundaries.  `similarity` is left at 0.0.
pub fn extend_seed(
    seed: &Seed,
    text: &[u8],
    sequence_starts: &[usize],
    lmax: usize,
    scores: &XdropScores,
) -> LtrCandidate {
    // Maximum extension length in either direction (seed already covers `len` bp)
    let max_ext = lmax.saturating_sub(seed.len);

    let seq_start = sequence_starts[seed.seq_num];
    // One past the sequence's last base: the separator before the next sequence, or the
    // terminator that ends the text.  Used as an exclusive end below.
    let seq_end = sequence_starts
        .get(seed.seq_num + 1)
        .map(|&s| s - 1)
        .unwrap_or(text.len() - 1);

    // ---- Left (backward) extension ----
    // Extend from pos1 (left LTR seed start) and pos2 (right LTR seed start)
    // leftward.  The extension window is bounded by seqstart and by max_ext.
    let left_avail_u = seed.pos1.saturating_sub(seq_start).min(max_ext);
    let left_avail_v = seed.pos2.saturating_sub(seq_start).min(max_ext);
    let left_len = left_avail_u.min(left_avail_v);

    let (left_ext_i, left_ext_j) = if left_len > 0 {
        let u_slice = &text[seed.pos1 - left_len..seed.pos1];
        let v_slice = &text[seed.pos2 - left_len..seed.pos2];
        xdrop_extend(u_slice, v_slice, false, scores)
    } else {
        (0, 0)
    };

    // ---- Right (forward) extension ----
    // Extend from the end of the seed rightward.
    // GT uses `seqend - (pos + len)` (exclusive end is seqend, not seqend+1).
    let seed_end1 = seed.pos1 + seed.len; // exclusive
    let seed_end2 = seed.pos2 + seed.len;

    let right_avail_u = seq_end.saturating_sub(seed_end1).min(max_ext);
    let right_avail_v = seq_end.saturating_sub(seed_end2).min(max_ext);
    let right_len = right_avail_u.min(right_avail_v);

    let (right_ext_i, right_ext_j) = if right_len > 0 {
        let u_slice = &text[seed_end1..seed_end1 + right_len];
        let v_slice = &text[seed_end2..seed_end2 + right_len];
        xdrop_extend(u_slice, v_slice, true, scores)
    } else {
        (0, 0)
    };

    // left_ext_i/j is how far we extended leftward into u/v (from seed start).
    // right_ext_i/j is how far we extended rightward from seed end.
    LtrCandidate {
        seq_num: seed.seq_num,
        left_ltr_start:  seed.pos1 - left_ext_i,
        right_ltr_start: seed.pos2 - left_ext_j,
        left_ltr_end:    seed.pos1 + seed.len - 1 + right_ext_i,
        right_ltr_end:   seed.pos2 + seed.len - 1 + right_ext_j,
        similarity: 0.0,
    }
}
