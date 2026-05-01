/// X-drop gapped extension, ported from GenomeTools `gt_evalxdroparbitscoresextend`.
///
/// Given two sequences `u` and `v` and a starting point at offset (0,0),
/// extends an alignment in the requested direction while the score stays
/// within `xdrop` of the running best.
///
/// Returns `(i_ext, j_ext)` — the number of characters consumed in `u` and
/// `v` respectively at the point of maximum score.
///
/// # Direction
/// * `forward = true`  — extend right: reads `u[0..]` and `v[0..]` forward.
/// * `forward = false` — extend left:  reads `u` and `v` from their ends
///   backward (equivalent to passing reversed slices and extending right).

// --------------------------------------------------
/// Alignment scores. Defaults match LTRHarvest's defaults.
#[derive(Clone, Copy, Debug)]
pub struct XdropScores {
    pub mat: i32,   // match award (positive)
    pub mis: i32,   // mismatch penalty (negative)
    pub ins: i32,   // insertion penalty (negative)
    pub del: i32,   // deletion penalty (negative)
    pub xdrop: i32, // drop threshold (positive)
}

impl Default for XdropScores {
    fn default() -> Self {
        XdropScores { mat: 2, mis: -2, ins: -3, del: -3, xdrop: 5 }
    }
}

// --------------------------------------------------
// Distance metric derived from scores (GenomeTools gt_calculatedistancesfromscores)
struct Dist {
    mis: i64,
    ins: i64,
    del: i64,
    gcd: i64,
    half_mat: i64,
}

fn gcd(a: i64, b: i64) -> i64 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

fn compute_dist(s: &XdropScores) -> Dist {
    // Double all scores if mat is odd (GenomeTools requirement)
    let (mat, mis, ins, del) = if s.mat % 2 != 0 {
        (s.mat as i64 * 2, s.mis as i64 * 2, s.ins as i64 * 2, s.del as i64 * 2)
    } else {
        (s.mat as i64, s.mis as i64, s.ins as i64, s.del as i64)
    };
    let g = gcd(gcd((mat - mis).abs(), (mat / 2 - ins).abs()), (mat / 2 - del).abs());
    Dist {
        mis: (mat - mis) / g,
        ins: (mat / 2 - ins) / g,
        del: (mat / 2 - del) / g,
        gcd: g,
        half_mat: mat / 2,
    }
}

// --------------------------------------------------
// Front storage: front[d*d + d + k] = furthest row i on diagonal k at distance d.
// j = i - k, so accessing (i, j) on diagonal k means i - j = k.

const MINUS_INF: i64 = i64::MIN / 4;

#[inline]
fn front_idx(d: i64, k: i64) -> usize {
    (d * d + d + k) as usize
}

fn get_front(fronts: &[i64], d: i64, k: i64) -> i64 {
    let idx = front_idx(d, k);
    *fronts.get(idx).unwrap_or(&MINUS_INF)
}

fn set_front(fronts: &mut Vec<i64>, d: i64, k: i64, val: i64) {
    let idx = front_idx(d, k);
    if idx >= fronts.len() {
        fronts.resize(idx + 128, MINUS_INF);
    }
    fronts[idx] = val;
}

// --------------------------------------------------
/// Compute LCP from positions (i, j) in u and v.
/// `forward=true`: reads u[i..] vs v[j..] left-to-right.
/// `forward=false`: reads u backward from u[ulen-1-i] vs v[vlen-1-j].
fn seq_lcp(u: &[u8], v: &[u8], i: usize, j: usize, forward: bool) -> usize {
    let mut n = 0;
    if forward {
        while i + n < u.len() && j + n < v.len() && u[i + n] == v[j + n] {
            n += 1;
        }
    } else {
        let ubase = u.len().saturating_sub(1 + i);
        let vbase = v.len().saturating_sub(1 + j);
        while n <= ubase && n <= vbase && u[ubase - n] == v[vbase - n] {
            n += 1;
        }
    }
    n
}

// --------------------------------------------------
/// X-drop extension entry point.
///
/// `u` and `v` are the two sequence slices to align, starting from the
/// boundary of the seed.  The direction (`forward`) determines whether we
/// extend right or left.
///
/// Returns `(i_ext, j_ext)` characters consumed in u and v at the best-
/// scoring alignment end-point.
pub fn xdrop_extend(
    u: &[u8],
    v: &[u8],
    forward: bool,
    scores: &XdropScores,
) -> (usize, usize) {
    if u.is_empty() || v.is_empty() {
        return (0, 0);
    }

    let dist = compute_dist(scores);
    let ulen = u.len() as i64;
    let vlen = v.len() as i64;
    let end_k = ulen - vlen; // diagonal of the endpoint (ulen, vlen)
    let xdrop = scores.xdrop as i64;

    // How many distance steps to look back for the X-drop comparison.
    let dback = (xdrop + dist.half_mat) / dist.gcd + 1;

    // Allowed consecutive all-MINUS_INF distance levels before termination.
    let allowed_inf_gen = dist.mis.max(dist.ins).max(dist.del) - 1;
    let mut inf_gen = 0i64;

    let mut fronts: Vec<i64> = Vec::new();
    // big_t[d] = running best score at distance ≤ d (monotone non-decreasing)
    let mut big_t: Vec<i64> = Vec::new();

    // ---- Phase 0: initial match run ----
    let init_lcp = seq_lcp(u, v, 0, 0, forward) as i64;
    set_front(&mut fronts, 0, 0, init_lcp);

    let init_score = init_lcp * 2 * dist.half_mat; // (i+j)*half_mat = 2*i*half_mat at d=0
    big_t.push(init_score);

    let mut best_score = init_score;
    let mut best_i = init_lcp as usize;
    let mut best_j = init_lcp as usize;

    // Determine initial diagonal bounds for the while loop.
    let (mut lbound, mut ubound) = if init_lcp >= ulen || init_lcp >= vlen {
        (1i64, -1i64) // alignment already finished; loop won't execute
    } else {
        (0i64, 0i64)
    };

    let mut currd = 0i64;

    while lbound <= ubound {
        currd += 1;
        let mut always_inf = true;

        for k in (lbound - 1)..=(ubound + 1) {
            let mut i = MINUS_INF;
            let mut any_case_applied = false;

            // Case 1: DELETION (advance in u only; come from diagonal k-1)
            //   condition: k > lbound  (so k-1 >= lbound)
            if k > lbound {
                let pd = currd - dist.del;
                if pd >= 0 && (k - 1).abs() <= pd {
                    let prev = get_front(&fronts, pd, k - 1);
                    if prev > MINUS_INF {
                        i = prev + 1;
                    }
                    any_case_applied = true;
                }
            }

            // Case 2: MISMATCH/REPLACEMENT (advance both; come from diagonal k)
            //   condition: k within [lbound, ubound]
            if k >= lbound && k <= ubound {
                let pd = currd - dist.mis;
                if pd >= 0 && k.abs() <= pd {
                    let row = get_front(&fronts, pd, k) + 1;
                    if row > i {
                        i = row;
                    }
                    any_case_applied = true;
                }
            }

            // Case 3: INSERTION (advance in v only; come from diagonal k+1)
            //   condition: k < ubound  (so k+1 <= ubound)
            if k < ubound {
                let pd = currd - dist.ins;
                if pd >= 0 && (k + 1).abs() <= pd {
                    let prev = get_front(&fronts, pd, k + 1);
                    if prev > i {
                        i = prev;
                    }
                    any_case_applied = true;
                }
            }

            if i < 0 {
                // No valid i: this diagonal is inactive.
                // If no case applied at all, GenomeTools sets always_inf=false
                // to avoid premature termination at boundary diagonals.
                if !any_case_applied {
                    always_inf = false;
                }
                set_front(&mut fronts, currd, k, MINUS_INF);
                continue;
            }

            let j = i - k;

            // X-drop check: if the score here is more than X below the
            // best score seen dback steps ago, prune this diagonal.
            let prev_d = currd - dback;
            if prev_d > 0 {
                if let Some(&prev_best) = big_t.get(prev_d as usize) {
                    let cur_score = (i + j) * dist.half_mat - currd * dist.gcd;
                    if cur_score < prev_best - xdrop {
                        set_front(&mut fronts, currd, k, MINUS_INF);
                        continue;
                    }
                }
            }

            // Extend along matching characters only if we made genuine progress
            // over the previous distance level on this diagonal.
            let prev_i = get_front(&fronts, currd - 1, k);
            if k <= -currd || k >= currd
                || (prev_i < i && i <= ulen.min(vlen + k))
            {
                // Genuine progress: extend along matches.
                if i < ulen && j < vlen {
                    let lcp_ext = seq_lcp(u, v, i as usize, j as usize, forward) as i64;
                    i += lcp_ext;
                }
                set_front(&mut fronts, currd, k, i);

                let j = i - k;
                let score = (i + j) * dist.half_mat - currd * dist.gcd;
                if score > best_score {
                    best_score = score;
                    best_i = i as usize;
                    best_j = j as usize;
                }
            } else {
                // No new progress; carry forward the old value.
                set_front(&mut fronts, currd, k, prev_i);
            }

            always_inf = false;
        }

        // Track consecutive all-inactive-front generations.
        if always_inf {
            inf_gen += 1;
            if inf_gen > allowed_inf_gen {
                break;
            }
        } else {
            inf_gen = 0;
        }

        // Append running best score to big_t for X-drop look-back.
        while big_t.len() <= currd as usize {
            big_t.push(i64::MIN / 2);
        }
        big_t[currd as usize] = best_score;

        // Early exit if we've reached the far corner of the DP matrix.
        if (-currd..=currd).contains(&end_k)
            && get_front(&fronts, currd, end_k) == ulen
        {
            break;
        }

        // ---- Update lbound / ubound for the next round ----

        // Prune lower bound: find smallest k with a valid front.
        for k in (lbound - 1)..=(ubound + 1) {
            if get_front(&fronts, currd, k) > MINUS_INF {
                lbound = k;
                break;
            }
        }

        // Prune upper bound: find largest k with a valid front.
        for k in ((lbound - 1)..=(ubound + 1)).rev() {
            if get_front(&fronts, currd, k) > MINUS_INF {
                ubound = k;
                break;
            }
        }

        // Boundary: if diagonal k has consumed all of v (front == vlen + k),
        // we cannot go below that diagonal.
        for k in (lbound..=0).rev() {
            if get_front(&fronts, currd, k) == vlen + k {
                lbound = k;
                break;
            }
        }

        // Boundary: if diagonal k has consumed all of u (front == ulen),
        // we cannot go above that diagonal.
        for k in 0..=ubound {
            if get_front(&fronts, currd, k) == ulen {
                ubound = k;
                break;
            }
        }
    }

    (best_i, best_j)
}
