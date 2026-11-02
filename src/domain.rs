//! Protein-domain evidence for each candidate.
//!
//! Of LTR_retriever's tests, the search of each candidate's six-frame
//! translation against transposable-element protein profiles (Pfam domains
//! grouped as LTR, not-LTR, Gypsy and Copia) separates real elements from
//! chance pairs best in our measurements.  LTR_retriever runs HMMER for it;
//! here the frames are written out as proteins and searched with nail
//! (MMseqs2 seeds, cloud search, Forward), whose full pipeline is the fastest
//! of what we measured: rice's 39,924 frames against TEfam.hmm's 102 profiles
//! take nail 206 CPU seconds, hmmsearch 250, and the in-process libnail
//! search this replaced about 18,000, with the same candidates hit (4,043 /
//! 4,044 at E <= 1e-5).
//!
//! Classification follows LTR_retriever's `annotate_TE.pl`: a candidate is
//! LTR-order when all its hits are to LTR-class profiles, notLTR when all are
//! to the not-LTR list, and "mixture" when the hits disagree (a hit outside
//! both lists counts as disagreement, as in the original).  Superfamily and
//! strand come from hit counts with the same ratios.

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{bail, Context, Result};
use rayon::prelude::*;

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
// The nail search.

/// nail's seed filter (`-S`, P-value; nail's default 1e-4).  At 0.01 the
/// rice search takes twice as long (70 s against 36 on 8 threads) and misses
/// 47 of the candidates hmmsearch finds at E <= 0.05 instead of 91; the
/// difference is weak hits, which old elements such as gorilla's depend on.
const SEED_P: f64 = 0.01;

/// Where nail and MMseqs2 are, and where nail works.
/// Oldest nail the search runs with.
pub const NAIL_MIN: (u32, u32, u32) = (0, 7, 1);

pub struct Nail {
    pub nail: PathBuf,
    pub mmseqs: PathBuf,
    pub threads: usize,
    /// Scratch directory for the frames, the profiles and nail's own files;
    /// removed after a successful search.
    pub work: PathBuf,
}

/// An executable named on the command line (a path, or a name looked up on
/// PATH).
pub fn find_program(name: &str) -> Result<PathBuf> {
    let p = Path::new(name);
    if p.components().count() > 1 {
        if p.is_file() {
            return Ok(p.to_path_buf());
        }
        bail!("{} not found", name);
    }
    std::env::var_os("PATH")
        .and_then(|path| std::env::split_paths(&path).map(|d| d.join(name)).find(|c| c.is_file()))
        .with_context(|| format!("{} not found on PATH", name))
}

/// Profile names by accession (without its version), in file order.
fn profile_names(hmm: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut name = String::new();
    for line in hmm.lines() {
        let mut f = line.split_whitespace();
        match f.next() {
            Some("NAME") => name = f.next().unwrap_or("").to_string(),
            Some("ACC") => out.push((f.next().unwrap_or("").split('.').next().unwrap_or("").to_string(), name.clone())),
            _ => {}
        }
    }
    out
}

/// Search the six frames of every sequence in `seqs` (None: not searched)
/// against the profiles in `hmm` with one nail run, and classify each
/// sequence's hits.  `num_targets` is the database size used for the
/// E-value, as hmmsearch would count it.
pub fn search(nail: &Nail, hmm: &Path, seqs: &[Option<&[u8]>], num_targets: usize, max_e: f64) -> Result<Vec<Option<DomainCall>>> {
    std::fs::create_dir_all(&nail.work).with_context(|| format!("creating {}", nail.work.display()))?;
    // nail 0.7.1 panics on a blank line in the profile file (LTR_retriever's
    // TEfam.hmm starts with one), so it gets a copy without them
    let hmm_text = std::fs::read_to_string(hmm).with_context(|| format!("reading {}", hmm.display()))?;
    let profiles = profile_names(&hmm_text);
    if profiles.is_empty() {
        bail!("no profiles with an ACC line in {}", hmm.display());
    }
    let hmm_copy = nail.work.join("profiles.hmm");
    std::fs::write(&hmm_copy, hmm_text.lines().filter(|l| !l.trim().is_empty()).map(|l| format!("{l}\n")).collect::<String>())?;

    // Frames as proteins: stops are X, as nail's alphabet has no '*'; frames
    // under 10 residues are left out.
    let frames_path = nail.work.join("frames.fa");
    let mut n_frames = 0usize;
    {
        let mut out = std::io::BufWriter::new(std::fs::File::create(&frames_path)?);
        for (ci, chunk) in seqs.chunks(4096).enumerate() {
            let text: Vec<(usize, Vec<u8>)> = chunk
                .par_iter()
                .enumerate()
                .map(|(j, s)| {
                    let mut buf = Vec::new();
                    let mut n = 0;
                    if let Some(s) = s {
                        for (fi, frame) in six_frames(s).iter().enumerate() {
                            if frame.len() >= 10 {
                                buf.extend(format!(">{}_{}\n", ci * 4096 + j, fi).bytes());
                                buf.extend(frame);
                                buf.push(b'\n');
                                n += 1;
                            }
                        }
                    }
                    (n, buf)
                })
                .collect();
            for (n, buf) in text {
                n_frames += n;
                out.write_all(&buf)?;
            }
        }
        out.flush()?;
    }
    let mut calls: Vec<Option<DomainCall>> = seqs.iter().map(|s| s.map(|_| DomainCall::default())).collect();
    if n_frames == 0 {
        return Ok(calls.into_iter().map(|c| c.map(|_| classify(Vec::new()))).collect());
    }

    let tbl = nail.work.join("hits.tbl");
    let out = Command::new(&nail.nail)
        .arg("search")
        .args(["-t", &nail.threads.max(1).to_string()])
        .args(["-E", &max_e.to_string(), "-Z", &num_targets.max(1).to_string(), "-S", &SEED_P.to_string()])
        .args(["-X", "--use-accession"])
        .arg("--mmseqs-path").arg(&nail.mmseqs)
        .arg("--tmp-dir").arg(nail.work.join("tmp"))
        .arg("--tbl-out").arg(&tbl)
        .arg(&hmm_copy)
        .arg(&frames_path)
        .output()
        .with_context(|| format!("running {}", nail.nail.display()))?;
    if !out.status.success() {
        let err = String::from_utf8_lossy(&out.stderr);
        let tail: Vec<&str> = err.lines().rev().take(15).collect();
        bail!("nail search failed ({}); work files kept in {}:\n{}", out.status, nail.work.display(), tail.into_iter().rev().collect::<Vec<_>>().join("\n"));
    }

    // target query tstart tend qstart qend score bias evalue cellfrac; one
    // row per (frame, profile) pair
    let order: HashMap<&str, usize> = profiles.iter().enumerate().map(|(i, (a, _))| (a.as_str(), i)).collect();
    let mut hits: Vec<Vec<(usize, DomainHit)>> = vec![Vec::new(); seqs.len()];
    for line in std::fs::read_to_string(&tbl).with_context(|| format!("reading {}", tbl.display()))?.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 9 {
            bail!("unexpected nail table line: {}", line);
        }
        let (idx, fi) = f[0].rsplit_once('_').context("nail target name")?;
        let (idx, fi): (usize, usize) = (idx.parse()?, fi.parse()?);
        let e: f64 = f[8].parse()?;
        if e > max_e {
            continue;
        }
        let accession = f[1].split('.').next().unwrap_or("").to_string();
        let rank = order.get(accession.as_str()).copied().unwrap_or(usize::MAX);
        let name = profiles.get(rank).map(|p| p.1.clone()).unwrap_or_default();
        // amino acid t (1-based) of frame fi covers strand bases
        // fi%3 + 3(t-1) .. fi%3 + 3t; the reverse frames are read on the
        // reverse complement
        let len = seqs[idx].map_or(0, |s| s.len());
        let (ts, te): (usize, usize) = (f[2].parse()?, f[3].parse()?);
        let o = fi % 3;
        let (a, b) = (o + 3 * (ts - 1), (o + 3 * te).min(len));
        let (start, end) = if fi < 3 { (a, b) } else { (len - b, len - a) };
        hits[idx].push((rank, DomainHit { accession, name, frame: fi, bits: f[6].parse()?, e_value: e, start, end }));
    }
    for (c, mut h) in calls.iter_mut().zip(hits) {
        if let Some(c) = c {
            h.sort_by_key(|(rank, d)| (d.frame, *rank));
            *c = classify(h.into_iter().map(|(_, d)| d).collect());
        }
    }
    std::fs::remove_dir_all(&nail.work).ok();
    Ok(calls)
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
    /// The alignment's target span, as 0-based offsets into the searched
    /// sequence on its forward strand (half-open).
    pub start: usize,
    pub end: usize,
}

#[derive(Debug, Clone, Default)]
pub struct DomainCall {
    pub hits: Vec<DomainHit>,
    /// LTR, notLTR, mixture, NA (hits only outside both lists), or none.
    pub order: String,
    /// Gypsy, Copia or unknown.
    pub superfamily: String,
    /// +, - or ?
    pub strand: String,
    /// Absolute text offset of the searched sequence's first base (set by
    /// the caller), so hit positions can be reported in genome coordinates.
    pub origin: usize,
}

impl DomainCall {
    /// Hits to profiles in the LTR-order list, as `accession:E-value:frame`
    /// strings.  Inside an LTR record a hit to one of the polyprotein
    /// domains means the record is coding sequence, not an LTR; some
    /// profiles on the list (DUF1263 and its relatives) do occur inside
    /// LTRs, so the caller decides which accessions count.
    pub fn ltr_order_hits(&self) -> Vec<String> {
        let ltr: HashSet<&str> = LTR.iter().copied().collect();
        self.hits
            .iter()
            .filter(|h| ltr.contains(h.accession.as_str()))
            .map(|h| format!("{}:{:.1e}:{}", h.accession, h.e_value, FRAME_NAMES[h.frame]))
            .collect()
    }
}

pub const FRAME_NAMES: [&str; 6] = ["aa1", "aa2", "aa3", "rev_aa1", "rev_aa2", "rev_aa3"];

/// Polyprotein domains that decide the veto: a hit to one of these inside an
/// LTR means the "LTR" is internal sequence.  The rest of the LTR list
/// (DUF834 and DUF1263 above all) turns up inside real LTRs: 748 annotated
/// LTR pairs on rice.
pub const CORE: &[&str] = &[
    "PF00665", "PF13683", // integrase core
    "PF03732", "PF14223", "PF14244", "PF13976", // gag, gag-pre-integrase
    "PF07727", // reverse transcriptase (copia)
    "PF13650", "PF13975", "PF09668", "PF00077", "PF08284", "PF02160", "PF12384", // aspartyl proteases
    "PF00385", // chromo
];

pub fn is_core(hit: &str) -> bool {
    CORE.iter().any(|a| hit.starts_with(a))
}

/// E-value a core hit needs before it counts as evidence of nesting.
pub const NEST_MAX_E: f64 = 1e-5;
/// Distance between two hits of one core profile on one strand that marks
/// the domain as duplicated rather than split across frames by a frameshift.
pub const NEST_DUP_GAP: usize = 1000;

/// Signs that a candidate is two elements, one inside the other, read from
/// its core polyprotein hits (`CORE`, E <= `NEST_MAX_E`): `opp` when there
/// are core hits on both strands, `dup` when one core profile hits twice on
/// one strand at least `NEST_DUP_GAP` bp apart.  A single element encodes
/// each core domain once, on one strand.  None when neither holds.
pub fn nesting(call: &DomainCall) -> Option<&'static str> {
    let core: Vec<&DomainHit> = call.hits.iter().filter(|h| is_core(&h.accession) && h.e_value <= NEST_MAX_E).collect();
    let plus = core.iter().any(|h| h.frame < 3);
    let minus = core.iter().any(|h| h.frame >= 3);
    if plus && minus {
        return Some("opp");
    }
    for (i, a) in core.iter().enumerate() {
        for b in &core[i + 1..] {
            if a.accession == b.accession && (a.frame < 3) == (b.frame < 3) {
                let gap = if a.start >= b.end { a.start - b.end } else if b.start >= a.end { b.start - a.end } else { 0 };
                if gap >= NEST_DUP_GAP {
                    return Some("dup");
                }
            }
        }
    }
    None
}

fn classify(hits: Vec<DomainHit>) -> DomainCall {
    if hits.is_empty() {
        return DomainCall { hits, order: "none".into(), superfamily: "unknown".into(), strand: "?".into(), ..Default::default() };
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
    DomainCall { hits, order: order.unwrap_or("none").to_string(), superfamily: superfamily.into(), strand: strand.into(), ..Default::default() }
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
    fn nesting_rules() {
        let mk = |acc: &str, frame: usize, start: usize, e: f64| DomainHit { accession: acc.into(), name: String::new(), frame, bits: 30.0, e_value: e, start, end: start + 300 };
        let call = |h: Vec<DomainHit>| DomainCall { hits: h, ..Default::default() };
        // one element: integrase and RT forward
        assert_eq!(nesting(&call(vec![mk("PF00665", 0, 1000, 1e-20), mk("PF07727", 1, 3000, 1e-30)])), None);
        // core domains on both strands
        assert_eq!(nesting(&call(vec![mk("PF00665", 0, 1000, 1e-20), mk("PF07727", 4, 6000, 1e-30)])), Some("opp"));
        // a weak opposite-strand hit does not count
        assert_eq!(nesting(&call(vec![mk("PF00665", 0, 1000, 1e-20), mk("PF07727", 4, 6000, 1e-3)])), None);
        // non-core hit on the other strand (accessory ORF) does not count
        assert_eq!(nesting(&call(vec![mk("PF00665", 0, 1000, 1e-20), mk("PF03078", 4, 6000, 1e-30)])), None);
        // integrase twice on one strand, 4 kb apart
        assert_eq!(nesting(&call(vec![mk("PF00665", 0, 1000, 1e-20), mk("PF00665", 2, 5300, 1e-20)])), Some("dup"));
        // two frames of one integrase (frameshift): adjacent, not duplicated
        assert_eq!(nesting(&call(vec![mk("PF00665", 0, 1000, 1e-20), mk("PF00665", 2, 1200, 1e-20)])), None);
    }

    #[test]
    fn classification_rules() {
        let mk = |acc: &str, frame: usize| DomainHit { accession: acc.into(), name: String::new(), frame, bits: 10.0, e_value: 1e-5, start: 0, end: 300 };
        let c = classify(vec![mk("PF00665", 0), mk("PF03732", 1)]);
        assert_eq!((c.order.as_str(), c.superfamily.as_str(), c.strand.as_str()), ("LTR", "Gypsy", "+"));
        let c = classify(vec![mk("PF00665", 0), mk("PF00078", 4)]);
        assert_eq!(c.order, "mixture");
        let c = classify(vec![mk("PF07727", 5)]);
        assert_eq!((c.order.as_str(), c.superfamily.as_str(), c.strand.as_str()), ("LTR", "Copia", "-"));
        assert_eq!(classify(vec![]).order, "none");
    }
}
