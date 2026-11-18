//! Re-extension of LTR boundaries with a local alignment.
//!
//! X-drop extension with the LTRharvest defaults stops at the first short
//! indel or cluster of mismatches between the two LTR copies, so many real
//! elements come out with LTRs short by tens to hundreds of bases.  This
//! stage takes each surviving candidate, cuts a window around each LTR wide
//! enough for the longest LTR we are prepared to believe, aligns the left
//! window against the right window with the RepeatMasker scoring scheme
//! (rmblast-lib, in process), and takes the best plus-strand HSP as the LTR
//! pair.  The HSP ends are then trimmed to a clean terminal window so that
//! a few chance matches past the true end do not become part of the LTR.
//!
//! The window sizes, scoring parameters and trim settings are the ones that
//! measured best against RepeatMasker-annotated elements in human,
//! Arabidopsis and rice; see `ltrstruct_without_ltr_retriever.md` in the
//! project directory.

use std::io::BufReader;

use anyhow::{Context, Result};
use rmblast_lib::encoding::{encode_iupac, NUCL_SENTINEL};
use rmblast_lib::hits::Strand;
use rmblast_lib::matrix::ScoreMatrix;
use rmblast_lib::options::{MtMode, SearchParams, SeedMode};
use rmblast_lib::search::engine::search_query_vs_subject;

/// RepeatMasker's nucleotide comparison matrix, embedded so the stage needs
/// no file at run time.
const COMPARISON_MATRIX: &str = include_str!("data/comparison.matrix");

/// Window and trimming settings for one run.
#[derive(Debug, Clone)]
pub struct ReextendConfig {
    /// Longest LTR the windows allow for.  Each window reaches this far past
    /// the current LTR end, plus `vicinity`.  About 2 kb measured best; wider
    /// windows run into neighbouring elements in nested genomes.
    pub ltr_max: usize,
    /// Extra bases beyond `ltr_max` on the outer side, for the TSD search
    /// that follows.
    pub vicinity: usize,
    /// Terminal window, in alignment columns, that must be clean for an end
    /// to stand.
    pub trim_window: usize,
    /// Matches required inside that window.
    pub trim_matches: usize,
    /// An HSP ending within this many bases of a window edge (that is not a
    /// sequence boundary) counts as reaching it.
    pub edge_margin: usize,
}

/// A scoring scheme and search parameters shared by every candidate.
pub struct Aligner {
    params: SearchParams,
    matrix: ScoreMatrix,
}

impl Aligner {
    /// The settings RepeatMasker uses for nucleotide searches, minus DUST,
    /// which would mask any low-complexity stretch inside an LTR.
    pub fn new() -> Result<Self> {
        let matrix = ScoreMatrix::from_reader(
            "comparison.matrix",
            BufReader::new(COMPARISON_MATRIX.as_bytes()),
        )
        .context("parsing the embedded comparison matrix")?;
        let params = SearchParams {
            gap_open: 20,
            gap_extend: 5,
            matrix_name: "comparison.matrix".to_string(),
            word_size: 14,
            // RepeatMasker's translation of a minimum score of 200: twice
            // it for the ungapped X-drop, half for the preliminary gapped
            // one, four times for the final gapped alignment.  A final
            // X-drop of 200 cut divergent 3 kb LTR pairs to a few hundred
            // bases on rice.
            xdrop_ungap: 400,
            xdrop_gap: 100,
            xdrop_gap_final: 800,
            min_raw_gapped_score: 200,
            ungapped_cutoff: None,
            complexity_adjust: true,
            dust: false,
            mask_level: 101,
            num_threads: 1,
            mt_mode: MtMode::SplitByDb,
            seed_mode: SeedMode::Combined,
        };
        Ok(Aligner { params, matrix })
    }
}

/// The re-extended boundaries of one candidate.  Coordinates are absolute
/// offsets into the concatenated text, inclusive, like the input.
#[derive(Debug, Clone)]
#[allow(dead_code)] // score, columns and trimmed feed the per-candidate table (milestone 2)
pub struct Reextended {
    pub left_ltr_start: usize,
    pub left_ltr_end: usize,
    pub right_ltr_start: usize,
    pub right_ltr_end: usize,
    /// Raw alignment score of the HSP.
    pub score: i32,
    /// Matches over aligned columns of the trimmed alignment, in percent.
    /// Gaps count as differences.
    pub identity: f64,
    /// Aligned columns after trimming.
    pub columns: usize,
    /// Bases removed by the end trim: left LTR 5', right LTR 5', left LTR 3',
    /// right LTR 3'.
    pub trimmed: (usize, usize, usize, usize),
    /// The HSP reached within a few bases of a window edge that was not a
    /// sequence boundary.  The LTR may be longer than the window allowed, or
    /// the pair sits inside a larger duplication.
    pub at_window_edge: bool,
    /// Self-repeat coverage of the two LTRs (`self_repeat_cov`), when the
    /// two-pass re-extension already computed it.
    pub self_cov: Option<(f64, f64)>,
}

/// Bases dropped from each end of an alignment so that its terminal
/// `window` columns are gap-free with at least `matches` matches.
/// Returns (query left, subject left, query right, subject right), or None
/// when no such window exists.
fn trim_ends(q: &[u8], s: &[u8], window: usize, matches: usize) -> Option<(usize, usize, usize, usize)> {
    let n = q.len();
    if window == 0 || n < window {
        return None;
    }
    let clean = |a: usize| -> bool {
        let mut m = 0;
        for i in a..a + window {
            if q[i] == NUCL_SENTINEL || s[i] == NUCL_SENTINEL {
                return false;
            }
            if q[i] == s[i] && q[i] < 4 {
                m += 1;
            }
        }
        m >= matches
    };
    let mut left = 0;
    while left + window <= n && !clean(left) {
        left += 1;
    }
    if left + window > n {
        return None;
    }
    let mut right = n;
    while right >= left + window && !clean(right - window) {
        right -= 1;
    }
    if right < left + window {
        return None;
    }
    let bases = |v: &[u8]| v.iter().filter(|&&b| b != NUCL_SENTINEL).count();
    Some((bases(&q[..left]), bases(&s[..left]), bases(&q[right..]), bases(&s[right..])))
}

/// Fraction of `seq` covered by plus-strand alignments of the sequence
/// against itself off the main diagonal.  A satellite array or a
/// microsatellite aligns to itself shifted by one period, so its coverage is
/// close to 1; an LTR made of unique sequence has none.
pub fn self_repeat_cov(aligner: &Aligner, seq: &[u8]) -> f64 {
    if seq.len() < 2 * aligner.params.word_size.max(1) {
        return 0.0;
    }
    let enc = encode_iupac(seq);
    let results = search_query_vs_subject(&enc, "q", &enc, "s", &aligner.params, &aligner.matrix);
    let mut covered = vec![false; seq.len()];
    for r in results.iter().filter(|r| r.hsp.strand == Strand::Plus) {
        let h = &r.hsp;
        if h.q_start == h.s_start {
            continue;
        }
        for c in covered.iter_mut().take((h.q_end as usize).min(seq.len())).skip(h.q_start as usize) {
            *c = true;
        }
    }
    covered.iter().filter(|&&c| c).count() as f64 / seq.len() as f64
}

/// Re-extend one candidate.  `seq_start` and `seq_end` are the first and last
/// base of the sequence the candidate lies in (absolute, inclusive).  Returns
/// None when no plus-strand HSP scores above the cutoff, in which case the
/// caller keeps the X-drop boundaries.
#[allow(clippy::too_many_arguments)]
pub fn reextend(
    aligner: &Aligner,
    cfg: &ReextendConfig,
    text: &[u8],
    seq_start: usize,
    seq_end: usize,
    left_ltr_start: usize,
    left_ltr_end: usize,
    right_ltr_start: usize,
    right_ltr_end: usize,
) -> Option<Reextended> {
    if right_ltr_start <= left_ltr_end + 1 {
        return None;
    }
    let ltr_len = (left_ltr_end - left_ltr_start + 1).max(right_ltr_end - right_ltr_start + 1);
    let outward = if ltr_len < cfg.ltr_max { cfg.ltr_max - ltr_len + cfg.vicinity } else { cfg.vicinity };
    let internal = right_ltr_start - left_ltr_end - 1;
    let inward = outward.min(internal / 2);

    // Window bounds, absolute and inclusive.  The two windows never overlap
    // because the inward reach is at most half the internal region.
    let qs = left_ltr_start.saturating_sub(outward).max(seq_start);
    let qe = (left_ltr_end + inward).min(seq_end);
    let ss = right_ltr_start.saturating_sub(inward).max(seq_start);
    let se = (right_ltr_end + outward).min(seq_end);
    if qe <= qs || se <= ss || ss <= qe {
        return None;
    }

    let query = encode_iupac(&text[qs..=qe]);
    let subject = encode_iupac(&text[ss..=se]);
    let results = search_query_vs_subject(&query, "left", &subject, "right", &aligner.params, &aligner.matrix);

    // The X-drop LTRs, in window coordinates (0-based, end exclusive).  The
    // windows can reach a neighbouring copy of the same LTR family, whose
    // alignment may score higher than the pair we started from; only HSPs
    // that overlap both X-drop LTRs are candidates for the re-extension.
    let (q_lo, q_hi) = (left_ltr_start - qs, left_ltr_end - qs + 1);
    let (s_lo, s_hi) = (right_ltr_start - ss, right_ltr_end - ss + 1);
    let best = results
        .iter()
        .filter(|r| r.hsp.strand == Strand::Plus)
        .filter(|r| {
            let h = &r.hsp;
            (h.q_start as usize) < q_hi && q_lo < h.q_end as usize
                && (h.s_start as usize) < s_hi && s_lo < h.s_end as usize
        })
        .max_by_key(|r| r.hsp.score)?;
    let hsp = &best.hsp;

    // Offsets in the HSP are 0-based, end exclusive, relative to the window
    // without its sentinels.
    let (q_start, q_end) = (hsp.q_start as usize, hsp.q_end as usize);
    let (s_start, s_end) = (hsp.s_start as usize, hsp.s_end as usize);
    let (ql, sl, qr, sr) = trim_ends(&hsp.q_seq, &hsp.s_seq, cfg.trim_window, cfg.trim_matches)?;

    // Identity over the trimmed columns.
    let n = hsp.q_seq.len();
    let (mut lead, mut trail) = (0usize, n);
    {
        // Recover the column range the trim kept, by walking the same bases.
        let (mut qb, mut sb) = (0usize, 0usize);
        while lead < n && (qb < ql || sb < sl) {
            if hsp.q_seq[lead] != NUCL_SENTINEL { qb += 1; }
            if hsp.s_seq[lead] != NUCL_SENTINEL { sb += 1; }
            lead += 1;
        }
        let (mut qb, mut sb) = (0usize, 0usize);
        while trail > lead && (qb < qr || sb < sr) {
            trail -= 1;
            if hsp.q_seq[trail] != NUCL_SENTINEL { qb += 1; }
            if hsp.s_seq[trail] != NUCL_SENTINEL { sb += 1; }
        }
    }
    let columns = trail.saturating_sub(lead);
    let matches = (lead..trail)
        .filter(|&i| hsp.q_seq[i] == hsp.s_seq[i] && hsp.q_seq[i] < 4)
        .count();
    let identity = if columns > 0 { 100.0 * matches as f64 / columns as f64 } else { 0.0 };

    let new_l1s = qs + q_start + ql;
    let new_l1e = qs + q_end - 1 - qr;
    let new_l2s = ss + s_start + sl;
    let new_l2e = ss + s_end - 1 - sr;
    if new_l1e <= new_l1s || new_l2e <= new_l2s || new_l2s <= new_l1e {
        return None;
    }

    let edge = cfg.edge_margin;
    let q_len = qe - qs + 1;
    let s_len = se - ss + 1;
    let at_window_edge = (q_start <= edge && qs > seq_start)
        || (q_end + edge >= q_len && qe < seq_end)
        || (s_start <= edge && ss > seq_start)
        || (s_end + edge >= s_len && se < seq_end);

    Some(Reextended {
        left_ltr_start: new_l1s,
        left_ltr_end: new_l1e,
        right_ltr_start: new_l2s,
        right_ltr_end: new_l2e,
        score: hsp.score,
        identity,
        columns,
        trimmed: (ql, sl, qr, sr),
        at_window_edge,
        self_cov: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn self_repeat_flags_tandem_arrays_only() {
        let aligner = Aligner::new().unwrap();
        let unit = b"ATCGAATGGAATGGAATCGAATGGAATCGAATGGACTCGAATGGAATCGAATG";
        let tandem: Vec<u8> = unit.iter().cycle().take(600).copied().collect();
        assert!(self_repeat_cov(&aligner, &tandem) > 0.9);
        // a pseudo-random sequence from a linear congruential generator
        let mut x: u64 = 12345;
        let random: Vec<u8> = (0..600).map(|_| { x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); b"ACGT"[(x >> 33) as usize % 4] }).collect();
        assert!(self_repeat_cov(&aligner, &random) < 0.1);
    }

    #[test]
    fn trim_removes_ragged_ends() {
        // BLASTNA: A=0 C=1 G=2 T=3, gap=15.  Two chance matches, a gap, then
        // a clean run, then a mismatching tail.
        let q = [0, 1, 15, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let s = [0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 2, 3];
        let (ql, sl, qr, sr) = trim_ends(&q, &s, 4, 4).unwrap();
        assert_eq!((ql, sl), (2, 3));
        assert_eq!((qr, sr), (2, 2));
    }

    #[test]
    fn trim_fails_when_nothing_is_clean() {
        let q = [0, 1, 2, 3];
        let s = [3, 2, 1, 0];
        assert!(trim_ends(&q, &s, 4, 4).is_none());
    }
}
