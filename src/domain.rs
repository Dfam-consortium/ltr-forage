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
//! to not-LTR profiles, and "mixture" when the hits disagree (a hit to a
//! profile of no class counts as disagreement, as in the original).
//! Superfamily and strand come from hit counts with the same ratios.
//!
//! Each profile's class and role come from the profile file itself when it
//! is annotated (see `ProfileTable`): the NAME carries a RepeatMasker
//! class, `rve#LTR` or `RVT_2#LTR/Copia`, and the DESC of the LTR
//! polyprotein domains ends in `[role=core]`.  A file without any class
//! (LTR_retriever's TEfam.hmm) gets annotate_TE.pl's lists below.

use std::collections::HashMap;
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

/// Polyprotein domains of the legacy lists (`CORE_ROLE` in an annotated
/// file): a hit to one of these inside an LTR means the "LTR" is internal
/// sequence.  The rest of the LTR list (DUF834 and DUF1263 above all) turns
/// up inside real LTRs: 748 annotated LTR pairs on rice.
const CORE: &[&str] = &[
    "PF00665", "PF13683", // integrase core
    "PF03732", "PF14223", "PF14244", "PF13976", // gag, gag-pre-integrase
    "PF07727", // reverse transcriptase (copia)
    "PF13650", "PF13975", "PF09668", "PF00077", "PF08284", "PF02160", "PF12384", // aspartyl proteases
    "PF00385", // chromo
];

// ---------------------------------------------------------------------------
// Profile classes from the profile file.

/// What a hit to a profile says about the order of the candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Order {
    Ltr,
    NotLtr,
    /// A profile of no class: counts as disagreement.
    Unclassified,
}

#[derive(Debug, Clone)]
pub struct ProfileClass {
    pub order: Order,
    /// Gypsy or Copia when the class says so ( the only superfamilies the
    /// call names ), otherwise None.
    pub superfamily: Option<&'static str>,
    /// An LTR polyprotein domain: evidence of coding sequence inside an LTR
    /// and of nesting.
    pub core: bool,
}

/// The role tag in a profile's DESC line.
const CORE_ROLE: &str = "[role=core]";

/// Profile classes by key (the accession without its version, or the NAME
/// before its '#' for a profile without an ACC line).
///
/// An annotated profile file names each profile `<name>#<class>`, the
/// class as RepeatMasker writes it ( `LTR`, `LTR/Gypsy`, `LTR/Copia`,
/// `LINE`, `DNA/hAT`, `RC/Helitron` ... ): a type of `LTR` is evidence for
/// the candidate, any other type evidence against it, and a subtype of
/// Gypsy or Copia counts toward the superfamily call.  An LTR profile whose
/// DESC ends in `[role=core]` is a polyprotein domain ( see `nesting` and
/// the LTR veto in decide.rs ); any other LTR profile is accessory
/// evidence.  HMMER and nail read both lines unchanged.
#[derive(Debug, Clone)]
pub struct ProfileTable {
    by_key: HashMap<String, ProfileClass>,
    /// From the file's own classes rather than the legacy lists.
    pub annotated: bool,
}

impl ProfileTable {
    /// annotate_TE.pl's lists, for an unannotated profile file.
    pub fn legacy() -> Self {
        let mut by_key = HashMap::new();
        for (list, order) in [(LTR, Order::Ltr), (NOT_LTR, Order::NotLtr)] {
            for &a in list {
                by_key.insert(a.to_string(), ProfileClass { order, superfamily: None, core: false });
            }
        }
        for (list, sf) in [(GYPSY, "Gypsy"), (COPIA, "Copia")] {
            for &a in list {
                by_key.entry(a.to_string()).or_insert(ProfileClass { order: Order::Unclassified, superfamily: None, core: false }).superfamily = Some(sf);
            }
        }
        for &a in CORE {
            if let Some(c) = by_key.get_mut(a) {
                c.core = true;
            }
        }
        ProfileTable { by_key, annotated: false }
    }

    /// The table for one profile file and its profiles' (key, name) in file
    /// order.  Annotated when any NAME carries a class; a profile without
    /// one in an annotated file is unclassified.
    pub fn from_hmm(hmm: &str) -> Result<(Self, Vec<(String, String)>)> {
        struct Rec { name: String, acc: Option<String>, desc: String }
        let mut recs: Vec<Rec> = Vec::new();
        let mut cur: Option<Rec> = None;
        for line in hmm.lines() {
            let tag = line.split_whitespace().next().unwrap_or("");
            let rest = line.get(tag.len()..).unwrap_or("").trim();
            match tag {
                "NAME" => cur = Some(Rec { name: rest.to_string(), acc: None, desc: String::new() }),
                "ACC" => if let Some(r) = cur.as_mut() { r.acc = Some(strip_version(rest).to_string()) },
                "DESC" => if let Some(r) = cur.as_mut() { r.desc = rest.to_string() },
                _ if tag.starts_with("//") => if let Some(r) = cur.take() { recs.push(r) },
                _ => {}
            }
        }
        if recs.is_empty() {
            bail!("no profiles in the profile file");
        }
        let annotated = recs.iter().any(|r| r.name.contains('#'));
        let mut table = if annotated { ProfileTable { by_key: HashMap::new(), annotated } } else { Self::legacy() };
        let mut keys = Vec::new();
        let mut unclassified = Vec::new();
        for r in &recs {
            let (stem, class) = match r.name.split_once('#') { Some((n, c)) => (n, Some(c)), None => (r.name.as_str(), None) };
            let key = r.acc.clone().unwrap_or_else(|| stem.to_string());
            if keys.iter().any(|(k, _): &(String, String)| *k == key) {
                bail!("profile {} appears twice in the profile file", key);
            }
            if annotated {
                let core = r.desc.ends_with(CORE_ROLE);
                if r.desc.contains("[role=") && !core {
                    bail!("profile {}: the only role tag is {} at the end of DESC (got '{}')", r.name, CORE_ROLE, r.desc);
                }
                let pc = match class {
                    None => {
                        unclassified.push(r.name.clone());
                        ProfileClass { order: Order::Unclassified, superfamily: None, core: false }
                    }
                    Some(c) => {
                        let (ty, sub) = c.split_once('/').unwrap_or((c, ""));
                        if ty.is_empty() || c.contains(char::is_whitespace) {
                            bail!("profile {}: bad class '{}'", r.name, c);
                        }
                        let order = if ty == "LTR" { Order::Ltr } else { Order::NotLtr };
                        if core && order != Order::Ltr {
                            bail!("profile {}: {} on a profile of class {} ( only LTR profiles are core )", r.name, CORE_ROLE, c);
                        }
                        let superfamily = match (order, sub) { (Order::Ltr, "Gypsy") => Some("Gypsy"), (Order::Ltr, "Copia") => Some("Copia"), _ => None };
                        ProfileClass { order, superfamily, core }
                    }
                };
                table.by_key.insert(key.clone(), pc);
            }
            keys.push((key, stem.to_string()));
        }
        if !unclassified.is_empty() {
            log_warn!("[8b protein domains] {} profile(s) without a class count as disagreement: {}",
                      unclassified.len(), unclassified.join(", "));
        }
        Ok((table, keys))
    }

    pub fn get(&self, key: &str) -> Option<&ProfileClass> {
        self.by_key.get(key)
    }
}

/// `PF00665.27` -> `PF00665`; a name without a numeric version unchanged.
fn strip_version(acc: &str) -> &str {
    match acc.rsplit_once('.') {
        Some((a, v)) if !v.is_empty() && v.bytes().all(|b| b.is_ascii_digit()) => a,
        _ => acc,
    }
}

/// The table of the run's profile file, set by `search`; the legacy lists
/// before that (and in unit tests).
static TABLE: std::sync::OnceLock<ProfileTable> = std::sync::OnceLock::new();
static LEGACY: std::sync::LazyLock<ProfileTable> = std::sync::LazyLock::new(ProfileTable::legacy);

pub fn table() -> &'static ProfileTable {
    TABLE.get().unwrap_or(&LEGACY)
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
// The nail search.

/// nail's seed filter (`-S`, P-value; nail's default 1e-4).  At 0.01 the
/// rice search takes twice as long (70 s against 36 on 8 threads) and misses
/// 47 of the candidates hmmsearch finds at E <= 0.05 instead of 91; the
/// difference is weak hits, which old elements such as gorilla's depend on.
const SEED_P: f64 = 0.01;

/// nail's `--prog-f`: its progressive seeding aligns a profile's prefilter
/// results in rounds and stops when under this fraction of a round are
/// hits (nail's default 0.01).  On a large target set that stops before
/// strong hits further down the list: on gorilla (2.7 M frames) 0 brings
/// back 97 hits, 30 of them under E = 1e-5, loses none, and costs no time
/// or memory.
const PROG_F: &str = "0";

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

/// Search the six frames of every sequence in `seqs` (None: not searched)
/// against the profiles in `hmm` with one nail run, and classify each
/// sequence's hits.  `num_targets` is the database size used for the
/// E-value (nail's -Z).
pub fn search(nail: &Nail, hmm: &Path, seqs: &[Option<&[u8]>], num_targets: usize, max_e: f64) -> Result<Vec<Option<DomainCall>>> {
    std::fs::create_dir_all(&nail.work).with_context(|| format!("creating {}", nail.work.display()))?;
    // nail 0.7.1 panics on a blank line in the profile file (LTR_retriever's
    // TEfam.hmm starts with one), so it gets a copy without them
    let hmm_text = std::fs::read_to_string(hmm).with_context(|| format!("reading {}", hmm.display()))?;
    let (tbl_classes, profiles) = ProfileTable::from_hmm(&hmm_text).with_context(|| format!("reading {}", hmm.display()))?;
    log_detail!("[8b protein domains] {} profiles, classes {}", profiles.len(),
              if tbl_classes.annotated { "from the profile file" } else { "from LTR_retriever's lists ( unannotated profile file )" });
    if TABLE.set(tbl_classes).is_err() {
        bail!("domain::search ran twice in one process; the profile classes are set once");
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
        std::fs::remove_dir_all(&nail.work).ok();
        return Ok(calls.into_iter().map(|c| c.map(|_| classify(Vec::new()))).collect());
    }

    let tbl = nail.work.join("hits.tbl");
    let out = Command::new(&nail.nail)
        .arg("search")
        .args(["-t", &nail.threads.max(1).to_string()])
        .args(["-E", &max_e.to_string(), "-Z", &num_targets.max(1).to_string(), "-S", &SEED_P.to_string()])
        .args(["-X", "--use-accession", "--prog-f", PROG_F])
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
        if idx >= seqs.len() || fi >= 6 {
            bail!("nail table names a target that was not searched: {}", line);
        }
        let e: f64 = f[8].parse()?;
        if e > max_e {
            continue;
        }
        // the accession, or for a profile without one its NAME ( class dropped )
        let accession = strip_version(f[1].split('#').next().unwrap_or("")).to_string();
        let rank = order.get(accession.as_str()).copied().unwrap_or(usize::MAX);
        let name = profiles.get(rank).map(|p| p.1.clone()).unwrap_or_default();
        // amino acid t (1-based) of frame fi covers strand bases
        // fi%3 + 3(t-1) .. fi%3 + 3t; the reverse frames are read on the
        // reverse complement
        let len = seqs[idx].map_or(0, |s| s.len());
        let (ts, te): (usize, usize) = (f[2].parse()?, f[3].parse()?);
        if ts == 0 || te < ts {
            bail!("nail table line with target coordinates {}..{}: {}", ts, te, line);
        }
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
        let t = table();
        self.hits
            .iter()
            .filter(|h| t.get(&h.accession).map_or(false, |c| c.order == Order::Ltr))
            .map(|h| format!("{}:{:.1e}:{}", h.accession, h.e_value, FRAME_NAMES[h.frame]))
            .collect()
    }
}

pub const FRAME_NAMES: [&str; 6] = ["aa1", "aa2", "aa3", "rev_aa1", "rev_aa2", "rev_aa3"];

/// A hit (`accession` or `accession:E-value:frame`) to a core polyprotein
/// profile: such a hit inside an LTR means the "LTR" is internal sequence.
pub fn is_core(hit: &str) -> bool {
    let key = hit.split(':').next().unwrap_or("");
    table().get(key).map_or(false, |c| c.core)
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
/// nail reports one row per (frame, profile), its best, so `dup` sees two
/// copies only when they lie in different reading frames of the strand
/// (two of three nested same-strand copies, by chance of frame).
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
    let t = table();
    let mut order: Option<&str> = None;
    let (mut n_gypsy, mut n_copia, mut n_plus, mut n_minus) = (0usize, 0usize, 0usize, 0usize);
    for h in &hits {
        let pc = t.get(&h.accession);
        let d = match pc.map(|c| c.order) { Some(Order::Ltr) => "LTR", Some(Order::NotLtr) => "notLTR", _ => "NA" };
        order = Some(match order {
            None => d,
            Some(prev) if prev == d => prev,
            Some(_) => "mixture",
        });
        match pc.and_then(|c| c.superfamily) {
            Some("Gypsy") => n_gypsy += 1,
            Some("Copia") => n_copia += 1,
            _ => {}
        }
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

    fn profile(name: &str, acc: Option<&str>, desc: &str) -> String {
        let acc = acc.map(|a| format!("ACC   {a}\n")).unwrap_or_default();
        format!("HMMER3/f [3.4 | Aug 2023]\nNAME  {name}\n{acc}DESC  {desc}\nLENG  10\nALPH  amino\n//\n")
    }

    #[test]
    fn profile_table_from_annotated_file() {
        let hmm = [
            profile("rve#LTR", Some("PF00665.27"), "Integrase core domain [role=core]"),
            profile("RVT_2#LTR/Copia", Some("PF07727.19"), "Reverse transcriptase (RNA-dependent DNA polymerase) [role=core]"),
            profile("DUF834#LTR/Gypsy", Some("PF05754.15"), "Domain of unknown function (DUF834)"),
            profile("MULE#DNA/MULE-MuDR", Some("PF10551.13"), "MULE transposase domain"),
            profile("myRT#LTR/Pao", None, "custom profile [role=core]"),
            profile("Other", Some("PF99999.1"), "no class"),
        ].concat();
        let (t, keys) = ProfileTable::from_hmm(&hmm).unwrap();
        assert!(t.annotated);
        assert_eq!(keys.iter().map(|k| k.0.as_str()).collect::<Vec<_>>(), ["PF00665", "PF07727", "PF05754", "PF10551", "myRT", "PF99999"]);
        assert_eq!(keys[0].1, "rve");
        let c = |k: &str| t.get(k).cloned().unwrap();
        assert!(c("PF00665").core && c("PF00665").order == Order::Ltr && c("PF00665").superfamily.is_none());
        assert_eq!(c("PF07727").superfamily, Some("Copia"));
        assert!(!c("PF05754").core && c("PF05754").superfamily == Some("Gypsy"));
        assert!(c("PF10551").order == Order::NotLtr && !c("PF10551").core);
        assert!(c("myRT").core && c("myRT").superfamily.is_none());
        assert_eq!(c("PF99999").order, Order::Unclassified);
        // not in the file: no class at all ( unlike the legacy lists )
        assert!(t.get("PF13650").is_none());
    }

    #[test]
    fn profile_table_errors_and_legacy() {
        assert!(ProfileTable::from_hmm(&profile("MULE#DNA", Some("PF10551.1"), "x [role=core]")).is_err());
        assert!(ProfileTable::from_hmm(&profile("rve#LTR", Some("PF00665.1"), "x [role=accessory]")).is_err());
        assert!(ProfileTable::from_hmm(&[profile("a#LTR", Some("PF1.1"), ""), profile("b#LTR", Some("PF1.2"), "")].concat()).is_err());
        // no class anywhere: annotate_TE.pl's lists
        let (t, _) = ProfileTable::from_hmm(&profile("rve", Some("PF00665.27"), "Integrase core domain")).unwrap();
        assert!(!t.annotated);
        assert!(t.get("PF00665").unwrap().core && t.get("PF00665").unwrap().superfamily == Some("Gypsy"));
        assert_eq!(t.get("PF00078").unwrap().order, Order::Unclassified);
        assert_eq!(t.get("PF10551").unwrap().order, Order::NotLtr);
        assert_eq!(strip_version("PF00665.27"), "PF00665");
        assert_eq!(strip_version("Gag_p24.v2"), "Gag_p24.v2");
    }
}
