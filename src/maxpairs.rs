/// Bottom-up suffix array traversal to enumerate all maximal repeats.
///
/// A repeat (pos1, pos2, len) is *maximal* when:
///   - text[pos1..pos1+len] == text[pos2..pos2+len]  (len-character match)
///   - The characters immediately to the left of each occurrence differ
///     (left-maximality), OR one occurrence starts at position 0
///   - len >= min_len
///
/// This is a direct port of the Kurtz bottom-up ESA traversal used in
/// GenomeTools' `gt_enumeratemaxpairs` (esa-maxpairs.c /
/// esa-bottomup-maxpairs.inc).  The key difference from GenomeTools is
/// that we use per-node position lists instead of a shared global array;
/// the algorithm is otherwise identical.

// --------------------------------------------------
// Alphabet encoding for left-extension characters
// A=0, C=1, G=2, T=3, everything else (separators, etc.) = SPECIAL
const SPECIAL: usize = 4;
const ALPHA: usize = 5; // 0..=4

// Sentinels stored in NodeInfo::common
const COMMON_INITIAL: u8 = u8::MAX;     // node created, no children yet
const COMMON_DIVERSE: u8 = u8::MAX - 1; // left-diverse (multiple left chars seen)

// Map a raw byte to a left-char index (0..=SPECIAL)
#[inline]
fn base_idx(c: u8) -> usize {
    match c {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => SPECIAL,
    }
}

// --------------------------------------------------
/// Per-node state accumulated during the bottom-up traversal.
#[derive(Clone, Debug)]
struct NodeInfo {
    /// COMMON_INITIAL  — no leaves processed yet
    /// COMMON_DIVERSE  — at least two different left chars seen
    /// 0..=SPECIAL     — all leaves so far share this left char
    common: u8,
    /// Positions of leaves grouped by their left-extension character.
    pos: [Vec<usize>; ALPHA],
}

impl NodeInfo {
    fn new() -> Self {
        Self {
            common: COMMON_INITIAL,
            pos: Default::default(),
        }
    }
}

// --------------------------------------------------
/// Update `common` when a new left-char index `new_c` (0..=SPECIAL as u8)
/// is observed.  Mirrors GenomeTools' CHECKCHAR macro:
///   if (father->commonchar != CC || CC >= ISLEFTDIVERSE) → ISLEFTDIVERSE
/// Here ISLEFTDIVERSE == SPECIAL == 4; anything >= 4 always goes diverse.
#[inline]
fn update_common(common: &mut u8, new_c: u8) {
    if *common != new_c || new_c >= SPECIAL as u8 {
        *common = COMMON_DIVERSE;
    }
}

// --------------------------------------------------
/// Emit all pairs (leaf_pos, existing_pos, depth) where the existing
/// position has a left-char different from `leaf_char`.
#[inline]
fn emit_leaf_pairs(
    father: &NodeInfo,
    leaf_char: usize,
    leaf_pos: usize,
    depth: usize,
    cb: &mut dyn FnMut(usize, usize, usize),
) {
    for ci in 0..ALPHA {
        if ci != leaf_char {
            for &p in &father.pos[ci] {
                let (a, b) = if leaf_pos < p { (leaf_pos, p) } else { (p, leaf_pos) };
                cb(a, b, depth);
            }
        }
    }
}

// --------------------------------------------------
/// Emit all cross-character pairs between father's positions and son's
/// positions (different left chars only).
#[inline]
fn emit_branch_pairs(
    father: &NodeInfo,
    son: &NodeInfo,
    depth: usize,
    cb: &mut dyn FnMut(usize, usize, usize),
) {
    for fc in 0..ALPHA {
        for sc in 0..ALPHA {
            if fc != sc {
                for &fp in &father.pos[fc] {
                    for &sp in &son.pos[sc] {
                        let (a, b) = if fp < sp { (fp, sp) } else { (sp, fp) };
                        cb(a, b, depth);
                    }
                }
            }
        }
    }
}

/// Sentinel value used in cleanup: a position that is NOT in the SA.
/// Must be >= text.len() to guarantee it's invalid.
pub const SENTINEL_POS: usize = usize::MAX;

// --------------------------------------------------
/// Process a leaf edge: suffix `leaf_pos` is being attached as a child of
/// `father` (whose LCP depth is `depth`).
///
/// `first_succ` is true only for the very first child of a newly-created
/// internal node.
///
/// Use `leaf_pos = SENTINEL_POS` for the end-of-input sentinel (matches GT's
/// wildcard at position 0): forces left-diversity but is NOT added to any
/// position list and therefore generates no spurious repeat pairs.
fn process_leaf(
    first_succ: bool,
    depth: usize,
    father: &mut NodeInfo,
    leaf_pos: usize,
    text: &[u8],
    min_len: usize,
    cb: &mut dyn FnMut(usize, usize, usize),
) {
    if depth < min_len {
        *father = NodeInfo::new();
        return;
    }

    // Determine left-extension character of this suffix.
    // SENTINEL_POS (= usize::MAX) is the wildcard sentinel (like GT's pos 0):
    // it forces left-diversity but is NOT a real suffix, so it is NOT added
    // to any position list and emits no pairs.
    // Position 0 is a real suffix; its left-char is text[-1] which doesn't
    // exist, so we treat it as SPECIAL (same as GT's INITIALCHAR).
    let is_sentinel = leaf_pos == SENTINEL_POS;
    let lc = if is_sentinel || leaf_pos == 0 {
        SPECIAL
    } else {
        base_idx(text[leaf_pos - 1])
    };
    // SENTINEL_POS and pos==0 force diversity (matches GT's INITIALCHAR logic)
    let force_diverse = is_sentinel || leaf_pos == 0;

    if first_succ {
        *father = NodeInfo::new();
        father.common = if force_diverse { COMMON_DIVERSE } else { lc as u8 };
        if !is_sentinel {
            father.pos[lc].push(leaf_pos);
        }
        return;
    }

    // Not first child: update diversity tracking
    if father.common != COMMON_DIVERSE {
        if force_diverse {
            father.common = COMMON_DIVERSE;
        } else {
            update_common(&mut father.common, lc as u8);
        }
    }

    if father.common == COMMON_DIVERSE && !is_sentinel {
        emit_leaf_pairs(father, lc, leaf_pos, depth, cb);
    }

    if !is_sentinel {
        father.pos[lc].push(leaf_pos);
    }
}

// --------------------------------------------------
/// Process a branching edge: a completed sub-interval `son` is being
/// attached as a child of `father` (whose LCP depth is `depth`).
///
/// When `first_succ` is true the son is the first child; GenomeTools
/// returns immediately in this case (father stays at COMMON_INITIAL).
/// `son` is `None` exactly when `first_succ` is true.
fn process_branch(
    first_succ: bool,
    depth: usize,
    father: &mut NodeInfo,
    son: Option<&NodeInfo>,
    min_len: usize,
    cb: &mut dyn FnMut(usize, usize, usize),
) {
    if depth < min_len {
        *father = NodeInfo::new();
        return;
    }

    if first_succ {
        // Father stays INITIAL; son's info is picked up on the next merge.
        return;
    }

    let son = son.expect("son must be Some when first_succ=false");

    // Update father's common char based on son's.
    // COMMON_INITIAL in son → CHECKCHAR(INITIALCHAR) → always diverse.
    if father.common != COMMON_DIVERSE {
        match son.common {
            COMMON_DIVERSE | COMMON_INITIAL => father.common = COMMON_DIVERSE,
            sc => update_common(&mut father.common, sc),
        }
    }

    if father.common == COMMON_DIVERSE {
        emit_branch_pairs(father, son, depth, cb);
    }

    // Merge son's position lists into father
    for i in 0..ALPHA {
        for &p in &son.pos[i] {
            father.pos[i].push(p);
        }
    }
}

// --------------------------------------------------
/// Enumerate all maximal repeats of length ≥ `min_len` in the text
/// represented by the given suffix array (`sa`) and LCP array (`lcp`).
///
/// For each repeat found, `callback(pos1, pos2, len)` is called with
/// `pos1 < pos2`.
///
/// # Arguments
/// * `sa`           — suffix array (length `num_suffixes`)
/// * `lcp`          — LCP array (length `num_suffixes`, `lcp[0]` = 0)
/// * `text`         — original text bytes
/// * `num_suffixes` — number of valid suffixes (may be < text.len() when
///                    non-ACGT positions are excluded)
/// * `min_len`      — minimum repeat length
/// * `callback`     — called for each maximal repeat `(pos1, pos2, len)`
pub fn enumerate_maximal_repeats(
    sa: &[usize],
    lcp: &[usize],
    text: &[u8],
    num_suffixes: usize,
    min_len: usize,
    cb: &mut dyn FnMut(usize, usize, usize),
) {
    if num_suffixes == 0 {
        return;
    }

    // Stack entries: (lcp_depth, left_boundary, node_info)
    let mut s_lcp: Vec<usize> = vec![0];
    let mut s_lb: Vec<usize> = vec![0];
    let mut s_info: Vec<NodeInfo> = vec![NodeInfo::new()];

    // A completed interval waiting to be attached to a new parent
    let mut last_iv: Option<(usize, usize, usize, NodeInfo)> = None; // (lcp, lb, rb, info)

    let mut first_from_root = true;

    for idx in 0..num_suffixes {
        let cur_sa = sa[idx];
        // The C ESA bottom-up traversal (SAIN path) uses a look-ahead LCP:
        // at step idx, lcpvalue = LCP(sa[idx], sa[idx+1]).
        // Sufr stores the standard look-back LCP: lcp[i] = LCP(sa[i-1], sa[i]).
        // Shift by one so each suffix is placed at the correct node depth.
        let cur_lcp = if idx + 1 < num_suffixes { lcp[idx + 1] } else { 0 };

        // ---- Part 1: add current suffix as leaf to current top, if depth matches ----
        let top_lcp = *s_lcp.last().unwrap();
        if cur_lcp <= top_lcp {
            let first_succ = top_lcp == 0 && first_from_root;
            if first_succ {
                first_from_root = false;
            }
            process_leaf(
                first_succ,
                top_lcp,
                s_info.last_mut().unwrap(),
                cur_sa,
                text,
                min_len,
                cb,
            );
        }

        debug_assert!(last_iv.is_none(), "last_iv must be None at loop start");

        // ---- Part 2: pop intervals whose depth exceeds cur_lcp ----
        while cur_lcp < *s_lcp.last().unwrap() {
            let p_lcp = s_lcp.pop().unwrap();
            let p_lb = s_lb.pop().unwrap();
            let p_info = s_info.pop().unwrap();
            let p_rb = idx;

            let top_lcp = *s_lcp.last().unwrap();
            if cur_lcp <= top_lcp {
                let first_succ = top_lcp == 0 && first_from_root;
                if first_succ {
                    first_from_root = false;
                }
                process_branch(
                    first_succ,
                    top_lcp,
                    s_info.last_mut().unwrap(),
                    Some(&p_info),
                    min_len,
                    cb,
                );
                last_iv = None;
            } else {
                // This popped interval becomes a child of the next push
                last_iv = Some((p_lcp, p_lb, p_rb, p_info));
            }
        }

        // ---- Part 3: push a new node when depth increases ----
        if cur_lcp > *s_lcp.last().unwrap() {
            if let Some((_li_lcp, li_lb, _li_rb, li_info)) = last_iv.take() {
                // In GenomeTools, the new node is pushed onto the same stack slot
                // just vacated by lastinterval, so it *inherits* lastinterval's
                // NodeInfo through memory reuse.  We replicate that by seeding the
                // new entry directly with li_info rather than a blank NodeInfo.
                s_lcp.push(cur_lcp);
                s_lb.push(li_lb);
                // If depth is below min_len the positions won't be useful anyway.
                if cur_lcp < min_len {
                    s_info.push(NodeInfo::new());
                } else {
                    s_info.push(li_info);
                }
                // process_branch(firstsucc=true) does nothing beyond the depth
                // check; we handled that above, so no call needed here.
            } else {
                // New node with current leaf as first child
                s_lcp.push(cur_lcp);
                s_lb.push(idx);
                s_info.push(NodeInfo::new());
                process_leaf(
                    true,
                    cur_lcp,
                    s_info.last_mut().unwrap(),
                    cur_sa,
                    text,
                    min_len,
                    cb,
                );
            }
        }
    }

    // ---- Cleanup: close any still-open interval at the top of the stack ----
    // In GenomeTools' SAIN path, `lastsuftabvalue` stays 0, which is the
    // wildcard/sentinel position (not a real suffix).  We use SENTINEL_POS
    // (= usize::MAX) as our equivalent: it forces left-diversity to flush
    // any remaining pairs, but is NOT added to the position list and
    // generates no spurious pairs.
    let top_lcp = *s_lcp.last().unwrap();
    if top_lcp > 0 {
        process_leaf(false, top_lcp, s_info.last_mut().unwrap(), SENTINEL_POS, text, min_len, cb);
    }
}
