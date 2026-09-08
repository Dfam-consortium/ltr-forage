//! Protein-domain evidence for each candidate, in process.
//!
//! Of LTR_retriever's tests, the search of each candidate's six-frame
//! translation against transposable-element protein profiles (Pfam domains
//! grouped as LTR, not-LTR, Gypsy and Copia) separates real elements from
//! chance pairs best in our measurements.  LTR_retriever runs HMMER for it.
//! This module does the same search with libnail's profile HMM Forward
//! algorithm, so it needs no external program.
//!
//! Seeding: nail's own pipeline seeds with MMseqs2 and then bounds the
//! Forward pass with its cloud search.  The profile set here is small
//! (about a hundred profiles) and each target is one candidate, so this
//! module scores every diagonal sharing a 3-mer with the profile consensus
//! as an ungapped alignment, and the best one, when it reaches a minimum
//! score, bounds a Forward pass on a rectangle of rows around it.  Scoring
//! is then nail's: Forward, composition-bias correction, E-value.  I tried
//! the cloud search first; it does not prune on chance seeds and computed
//! most of the matrix for every pair.
//!
//! Classification follows LTR_retriever's `annotate_TE.pl`: a candidate is
//! LTR-order when all its hits are to LTR-class profiles, notLTR when all are
//! to the not-LTR list, and "mixture" when the hits disagree (a hit outside
//! both lists counts as disagreement, as in the original).  Superfamily and
//! strand come from hit counts with the same ratios.
//!
//! The HMMER3 profile reader is adapted from nail's `io/p7hmm.rs`
//! (Jack Roddy, BSD-3-Clause).

use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::time::{Duration, Instant};

use anyhow::{bail, Context, Result};
use libnail::align::structs::{DpMatrixSparse, RowBounds};
use libnail::align::{backward, e_value, forward, null_one_score, null_two_score, p_value, posterior, Bits};
use libnail::structs::profile::{ProfileBuilder, Transition};
use libnail::structs::{Profile, Sequence};

// ---------------------------------------------------------------------------
// Profile classes, from LTR_retriever's annotate_TE.pl (accessions without
// their version suffix).

const NOT_LTR: &[&str] = &[
    "PF10551", "PF13966", "PF04434", "PF03108", "PF13359", "PF00872", "PF05699", "PF13963", "PF04827",
    "PF05970", "PF03101", "PF13960", "PF02992", "PF14214", "PF12776", "PF13952", "PF03372", "PF14372",
    "PF14303", "PF13604", "PF01498", "PF01609", "PF13837", "PF13538", "PF02892", "PF13384", "PF06839",
    "PF03004", "PF04937", "PF13191", "PF02902", "PF12728", "PF14291", "PF14529", "PF10446", "PF14111",
    "PF13873", "PF13245", "PF04147", "PF13542", "PF05285", "PF00376", "PF00249", "PF00646", "PF12116",
    "PF09322", "PF00437", "PF13401", "PF06869", "PF13857", "PF04967", "PF00564", "PF13412", "PF12796",
    "PF06072",
];

const LTR: &[&str] = &[
    "PF00665", "PF03732", "PF14223", "PF07727", "PF13650", "PF13975", "PF13976", "PF08284", "PF00385",
    "PF13683", "PF09337", "PF13961", "PF09668", "PF14244", "PF13917", "PF04195", "PF03578", "PF14787",
    "PF02160", "PF05754", "PF04094", "PF16025", "PF11835", "PF12435", "PF05585", "PF05377", "PF04582",
    "PF12353", "PF06882", "PF08847", "PF10046", "PF17150", "PF12384", "PF01519", "PF00077", "PF15450",
    "PF11180", "PF08149", "PF07889", "PF12425", "PF07197", "PF05837", "PF04513", "PF03154", "PF04111",
    "PF04102", "PF06827",
];

const GYPSY: &[&str] = &[
    "PF13650", "PF13975", "PF08284", "PF00385", "PF09337", "PF13456", "PF04195", "PF03578", "PF00075",
    "PF05754", "PF04094", "PF05585", "PF06882", "PF00077", "PF10536", "PF00078", "PF09668", "PF02160",
    "PF16025", "PF11835", "PF04582", "PF05377", "PF17150", "PF10046", "PF12384", "PF12353", "PF07197",
    "PF15450", "PF04513", "PF08149", "PF03732", "PF00665",
];

const COPIA: &[&str] = &[
    "PF07727", "PF13961", "PF14244", "PF03564", "PF13976", "PF14223", "PF08847", "PF14787", "PF13917",
];

// ---------------------------------------------------------------------------
// HMMER3 profile file reader.

/// Parse every profile in a HMMER3 (`HMMER3/f`) file.
pub fn read_hmmer3<P: AsRef<Path>>(path: P) -> Result<Vec<Profile>> {
    let bytes = std::fs::read(&path)
        .with_context(|| format!("reading {}", path.as_ref().display()))?;
    let text = String::from_utf8_lossy(&bytes);
    let mut profiles = Vec::new();
    for record in text.split("\n//") {
        let record = record.trim();
        if !record.contains("\nHMM ") && !record.starts_with("HMM ") {
            continue;
        }
        profiles.push(parse_record(record).with_context(|| {
            let name = record.lines().find(|l| l.starts_with("NAME")).unwrap_or("?");
            format!("parsing profile at '{}'", name)
        })?);
    }
    if profiles.is_empty() {
        bail!("no HMMER3 profiles found in {}", path.as_ref().display());
    }
    Ok(profiles)
}

fn parse_floats<const N: usize>(fields: &[&str]) -> Result<[f32; N]> {
    let mut out = [0.0f32; N];
    if fields.len() < N {
        bail!("expected {} values, found {}", N, fields.len());
    }
    for (o, f) in out.iter_mut().zip(fields.iter()) {
        *o = if *f == "*" { 0.0 } else { (-f.parse::<f32>()?).exp() };
    }
    Ok(out)
}

fn parse_record(record: &str) -> Result<Profile> {
    let mut prf = ProfileBuilder::default();
    let mut in_model = false;
    let mut seen_compo = false;
    let mut pos = 0usize;
    // 0: expecting a match line, 1: insert line, 2: transition line
    let mut phase = 0;
    for line in record.lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.is_empty() {
            continue;
        }
        if !in_model {
            match fields[0] {
                "NAME" => { prf.name(fields.get(1).unwrap_or(&"").to_string()); }
                "ACC" => { prf.accession(fields.get(1).unwrap_or(&"").to_string()); }
                "LENG" => { prf.length(fields.get(1).unwrap_or(&"0").parse()?); }
                "ALPH" => {
                    if fields.get(1).copied() != Some("amino") {
                        bail!("non-amino alphabet {:?}", fields.get(1));
                    }
                }
                "STATS" => {
                    if fields.get(2).copied() == Some("FORWARD") {
                        prf.fwd_tau(fields[3].parse()?);
                        prf.fwd_lambda(fields[4].parse()?);
                    }
                }
                "HMM" => { in_model = true; }
                _ => {}
            }
            continue;
        }
        if !seen_compo {
            // The "m->m m->i ..." header line, then COMPO.  Node 0 has no
            // match line: its insert line is skipped and its transition line
            // is read like any other node's.
            if fields[0] == "COMPO" {
                seen_compo = true;
                phase = 1;
            }
            continue;
        }
        match phase {
            0 => {
                let new_pos: usize = fields[0].parse().context("node number")?;
                if new_pos != pos + 1 {
                    bail!("non-linear node numbers {} -> {}", pos, new_pos);
                }
                pos = new_pos;
                let em = parse_floats::<20>(&fields[1..])?;
                prf.mat_emission(pos, em);
                phase = 1;
            }
            1 => { phase = 2; }
            _ => {
                let t = parse_floats::<7>(&fields)?;
                // file order: MM MI MD IM II DM DD; Profile order per Transition
                // (MM=0 IM=1 MI=2 DM=3 II=4 MD=5 DD=6).
                let mut swapped = [0.0f32; 7];
                let order = [Transition::MM as usize, Transition::MI as usize, Transition::MD as usize, Transition::IM as usize, Transition::II as usize, Transition::DM as usize, Transition::DD as usize];
                for (i, idx) in order.iter().enumerate() {
                    swapped[*idx] = t[i];
                }
                prf.transition(pos, swapped);
                phase = 0;
            }
        }
    }
    prf.build()
}

// ---------------------------------------------------------------------------
// Six-frame translation.

fn codon(c: &[u8]) -> u8 {
    let idx = |b: u8| -> Option<usize> {
        match b { b'T' => Some(0), b'C' => Some(1), b'A' => Some(2), b'G' => Some(3), _ => None }
    };
    match (idx(c[0]), idx(c[1]), idx(c[2])) {
        (Some(a), Some(b), Some(d)) => {
            // Standard code, TCAG order; stops as X so the profile search
            // treats them as unknown residues.
            const AA: &[u8; 64] = b"FFLLSSSSYYXXCCXWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
            AA[a * 16 + b * 4 + d]
        }
        _ => b'X',
    }
}

fn translate(seq: &[u8]) -> Vec<u8> {
    seq.chunks_exact(3).map(codon).collect()
}

/// The six frames of `seq`, in LTR_retriever's order: aa1, aa2, aa3 forward,
/// then rev_aa1, rev_aa2, rev_aa3 of the reverse complement.
pub fn six_frames(seq: &[u8]) -> [Vec<u8>; 6] {
    let rc: Vec<u8> = seq
        .iter()
        .rev()
        .map(|&c| match c { b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', o => o })
        .collect();
    [
        translate(seq),
        translate(&seq[1.min(seq.len())..]),
        translate(&seq[2.min(seq.len())..]),
        translate(&rc),
        translate(&rc[1.min(rc.len())..]),
        translate(&rc[2.min(rc.len())..]),
    ]
}

// ---------------------------------------------------------------------------
// Profile set with a 3-mer index of each consensus for seeding.

const SEED_K: usize = 3;
/// Minimum ungapped diagonal score, in bits, for a (frame, profile) pair to
/// get a Forward pass.
const SEED_MIN_BITS: f32 = 18.0;

pub struct ProfileEntry {
    pub profile: Profile,
    pub accession: String,
    pub name: String,
    kmers: HashMap<[u8; SEED_K], Vec<usize>>,
}

pub struct DomainDb {
    pub entries: Vec<ProfileEntry>,
}

/// The best ungapped segment on one diagonal: profile and target positions
/// (1-based, inclusive) and its score in bits.
#[derive(Debug, Clone, Copy)]
struct DiagSeed {
    seq_start: usize,
    seq_end: usize,
    prf_start: usize,
    prf_end: usize,
    bits: f32,
}

impl DomainDb {
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let entries = read_hmmer3(path)?
            .into_iter()
            .map(|profile| {
                let consensus: Vec<u8> = profile.consensus_seq_bytes_utf8.iter().map(|b| b.to_ascii_uppercase()).collect();
                let mut kmers: HashMap<[u8; SEED_K], Vec<usize>> = HashMap::new();
                // consensus_seq_bytes_utf8 is 1-based (index 0 is padding),
                // so the window index is the 1-based profile position.
                for (i, w) in consensus.windows(SEED_K).enumerate() {
                    if i == 0 {
                        continue;
                    }
                    kmers.entry([w[0], w[1], w[2]]).or_default().push(i);
                }
                let accession = profile.accession.split('.').next().unwrap_or("").to_string();
                let name = profile.name.clone();
                ProfileEntry { profile, accession, name, kmers }
            })
            .collect();
        Ok(DomainDb { entries })
    }

    /// Score every diagonal that shares a 3-mer between `frame` and the
    /// profile consensus as an ungapped local alignment against the
    /// profile's match emissions, and return the best segment.
    fn seed(entry: &ProfileEntry, frame: &[u8], target: &Sequence) -> Option<DiagSeed> {
        let prf = &entry.profile;
        let m = prf.length;
        let l = target.length;
        let mut diagonals: HashSet<i64> = HashSet::new();
        for (t, w) in frame.windows(SEED_K).enumerate() {
            if let Some(ps) = entry.kmers.get(&[w[0], w[1], w[2]]) {
                for &p in ps {
                    diagonals.insert((t + 1) as i64 - p as i64);
                }
            }
        }
        let mut best: Option<DiagSeed> = None;
        for d in diagonals {
            // walk (t, p) with t - p = d over the overlap of both sequences
            let p0 = if d >= 0 { 1 } else { (1 - d) as usize };
            let t0 = (p0 as i64 + d) as usize;
            let mut run = 0.0f32;
            let mut run_start = (t0, p0);
            let mut seg_best = 0.0f32;
            let mut seg = ((t0, p0), (t0, p0));
            let (mut t, mut p) = (t0, p0);
            while t <= l && p <= m {
                let sc = prf.match_score(target.digital_bytes[t] as usize, p);
                if run <= 0.0 {
                    run = 0.0;
                    run_start = (t, p);
                }
                run += sc;
                if run > seg_best {
                    seg_best = run;
                    seg = (run_start, (t, p));
                }
                t += 1;
                p += 1;
            }
            let bits = seg_best / std::f32::consts::LN_2;
            if bits >= SEED_MIN_BITS && best.map_or(true, |b| bits > b.bits) {
                best = Some(DiagSeed { seq_start: seg.0 .0, seq_end: seg.1 .0, prf_start: seg.0 .1, prf_end: seg.1 .1, bits });
            }
        }
        best
    }
}

// ---------------------------------------------------------------------------
// Search and classification.

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct DomainHit {
    pub accession: String,
    pub name: String,
    /// Frame index 0..6 in `six_frames` order.
    pub frame: usize,
    pub bits: f32,
    pub e_value: f64,
}

/// Counts of pairs entering each stage, for the run summary.
#[derive(Debug, Clone, Copy, Default)]
pub struct SearchStats {
    pub pairs: usize,
    pub seeded: usize,
    pub forward_passed: usize,
    pub seed_time: Duration,
    pub align_time: Duration,
}

impl std::ops::AddAssign for SearchStats {
    fn add_assign(&mut self, o: Self) {
        self.pairs += o.pairs;
        self.seeded += o.seeded;
        self.forward_passed += o.forward_passed;
        self.seed_time += o.seed_time;
        self.align_time += o.align_time;
    }
}

#[derive(Debug, Clone, Default)]
pub struct DomainCall {
    pub stats: SearchStats,
    pub hits: Vec<DomainHit>,
    /// LTR, notLTR, mixture, NA (hits only outside both lists), or none.
    pub order: String,
    /// Gypsy, Copia or unknown.
    pub superfamily: String,
    /// +, - or ?
    pub strand: String,
}

/// nail's Forward filter threshold (its CLI default).
const FORWARD_P_THRESHOLD: f64 = 1e-4;

/// Per-thread working state: profile copies (their special transitions are
/// set per target length) and the dynamic-programming matrices.
pub struct Searcher {
    profiles: Vec<Profile>,
    fwd: DpMatrixSparse,
    bwd: DpMatrixSparse,
    post: DpMatrixSparse,
    bounds: RowBounds,
}

impl Searcher {
    pub fn new(db: &DomainDb) -> Self {
        Searcher {
            profiles: db.entries.iter().map(|e| e.profile.clone()).collect(),
            fwd: DpMatrixSparse::default(),
            bwd: DpMatrixSparse::default(),
            post: DpMatrixSparse::default(),
            bounds: RowBounds::default(),
        }
    }

    /// Search the six frames of `seq` against every profile.  `num_targets`
    /// is the database size used for the E-value, as hmmsearch would count
    /// it (the number of frame sequences searched in the run).
    pub fn search(&mut self, db: &DomainDb, seq: &[u8], num_targets: usize, max_e: f64) -> Result<DomainCall> {
        let frames = six_frames(seq);
        let mut hits = Vec::new();
        let mut stats = SearchStats::default();
        for (fi, frame) in frames.iter().enumerate() {
            if frame.len() < 10 {
                continue;
            }
            let target = Sequence::from_utf8(frame).context("building the translated frame")?;
            for (entry, prf) in db.entries.iter().zip(self.profiles.iter_mut()) {
                stats.pairs += 1;
                let t0 = Instant::now();
                let seed = DomainDb::seed(entry, frame, &target);
                stats.seed_time += t0.elapsed();
                let Some(seed) = seed else { continue };
                stats.seeded += 1;
                let t0 = Instant::now();
                prf.configure_for_target_length(target.length);
                let scored = score_pair(prf, &target, &seed, &mut self.bounds, &mut self.fwd, &mut self.bwd, &mut self.post, num_targets);
                stats.align_time += t0.elapsed();
                let Some((bits, e)) = scored else { continue };
                stats.forward_passed += 1;
                if e <= max_e {
                    hits.push(DomainHit { accession: entry.accession.clone(), name: entry.name.clone(), frame: fi, bits, e_value: e });
                }
            }
        }
        let mut call = classify(hits);
        call.stats = stats;
        Ok(call)
    }
}

/// nail's align stage on a rectangle of rows around the seed: Forward, then
/// Backward, posterior and the composition-bias correction.  Returns the
/// corrected bit score and E-value, or None when the Forward filter rejected
/// the pair.
#[allow(clippy::too_many_arguments)]
fn score_pair(
    prf: &Profile,
    target: &Sequence,
    seed: &DiagSeed,
    bounds: &mut RowBounds,
    fwd: &mut DpMatrixSparse,
    bwd: &mut DpMatrixSparse,
    post: &mut DpMatrixSparse,
    num_targets: usize,
) -> Option<(f32, f64)> {
    // Rows: where the rest of the profile would land if the seed diagonal
    // continued, with slack for insertions in the target.
    let slack = prf.length / 2 + 20;
    let row_start = seed.seq_start.saturating_sub(seed.prf_start + slack).max(1);
    let row_end = (seed.seq_end + (prf.length - seed.prf_end) + slack).min(target.length);
    bounds.fill_rectangle(row_start, 1, row_end, prf.length);
    let null_one = null_one_score(target.length);
    fwd.reuse(target.length, prf.length, bounds);
    let forward_score = (forward(prf, target, fwd, bounds) - null_one).to_bits();
    if p_value(forward_score, prf.fwd_lambda, prf.fwd_tau) >= FORWARD_P_THRESHOLD {
        return None;
    }
    bwd.reuse(target.length, prf.length, bounds);
    post.reuse(target.length, prf.length, bounds);
    backward(prf, target, bwd, bounds);
    posterior(prf, fwd, bwd, post, bounds);
    let null_two = Bits(null_two_score(post, prf, target, bounds).to_bits().0.max(0.0));
    let bits = forward_score - null_two;
    let p = p_value(bits, prf.fwd_lambda, prf.fwd_tau);
    Some((bits.0, e_value(p, num_targets)))
}

fn classify(hits: Vec<DomainHit>) -> DomainCall {
    if hits.is_empty() {
        return DomainCall { stats: SearchStats::default(), hits, order: "none".into(), superfamily: "unknown".into(), strand: "?".into() };
    }
    let ltr: HashSet<&str> = LTR.iter().copied().collect();
    let not_ltr: HashSet<&str> = NOT_LTR.iter().copied().collect();
    let gypsy: HashSet<&str> = GYPSY.iter().copied().collect();
    let copia: HashSet<&str> = COPIA.iter().copied().collect();
    let mut order: Option<&str> = None;
    let (mut n_gypsy, mut n_copia, mut n_plus, mut n_minus) = (0usize, 0usize, 0usize, 0usize);
    for h in &hits {
        let acc = h.accession.as_str();
        let d = if ltr.contains(acc) { "LTR" } else if not_ltr.contains(acc) { "notLTR" } else { "NA" };
        order = Some(match order {
            None => d,
            Some(prev) if prev == d => prev,
            Some(_) => "mixture",
        });
        if gypsy.contains(acc) { n_gypsy += 1; }
        if copia.contains(acc) { n_copia += 1; }
        if h.frame < 3 { n_plus += 1; } else { n_minus += 1; }
    }
    let superfamily = if n_gypsy > n_copia && (n_copia as f64) / (n_gypsy as f64) < 0.3 {
        "Gypsy"
    } else if n_copia > n_gypsy && (n_gypsy as f64) / (n_copia as f64) < 0.9 {
        "Copia"
    } else {
        "unknown"
    };
    let strand = if n_plus > n_minus && (n_minus as f64) / (n_plus as f64) < 0.4 {
        "+"
    } else if n_minus > n_plus && (n_plus as f64) / (n_minus as f64) < 0.4 {
        "-"
    } else {
        "?"
    };
    DomainCall { stats: SearchStats::default(), hits, order: order.unwrap_or("none").to_string(), superfamily: superfamily.into(), strand: strand.into() }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translation_matches_the_standard_code() {
        assert_eq!(translate(b"ATGGCCTGA"), b"MAX".to_vec());
        let f = six_frames(b"ATGGCCTGAT");
        assert_eq!(f[0], b"MAX".to_vec());
        assert_eq!(f[1], b"WPD".to_vec());
        assert_eq!(f[3], translate(b"ATCAGGCCAT"));
    }

    #[test]
    fn classification_rules() {
        let mk = |acc: &str, frame: usize| DomainHit { accession: acc.into(), name: String::new(), frame, bits: 10.0, e_value: 1e-5 };
        let c = classify(vec![mk("PF00665", 0), mk("PF03732", 1)]);
        assert_eq!((c.order.as_str(), c.superfamily.as_str(), c.strand.as_str()), ("LTR", "Gypsy", "+"));
        let c = classify(vec![mk("PF00665", 0), mk("PF00078", 4)]);
        assert_eq!(c.order, "mixture");
        let c = classify(vec![mk("PF07727", 5)]);
        assert_eq!((c.order.as_str(), c.superfamily.as_str(), c.strand.as_str()), ("LTR", "Copia", "-"));
        assert_eq!(classify(vec![]).order, "none");
    }
}
