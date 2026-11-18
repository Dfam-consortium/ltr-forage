//! tRNA gene detection: a Rust port of the cytosolic tRNA search in
//! ARAGORN 1.2.41 (`aragorn -t -gc1 -l`, with `-i` intron search).
//!
//! ARAGORN is Copyright (C) 2003-Present Dean Laslett and is distributed under
//! the GNU General Public License version 3.  This file is a translation of
//! its `tmioptimise`, `find_tstems`, `find_astem5`, `vloop_stability`,
//! `aatail`, `ti_genedetected`, `find_slot`/`nearest_trna_gene`,
//! `remove_overlapping_trna`, `gene_sort`, `name` and the `iopt_fastafile`
//! window driver, and is likewise licensed under the GNU GPL v3 (or later).
//! Scores, thresholds and tie-breaks are kept identical so the calls match
//! the C program; tmRNA, mitochondrial tRNA and all output formatting are
//! left out.
//!
//! Reference: Laslett, D. and Canback, B. (2004) ARAGORN, a program for the
//! detection of transfer RNA and transfer-messenger RNA genes in nucleotide
//! sequences.  Nucleic Acids Research, 32:11-16.
//!
//! The C program scans a sequence in 20 kb windows that overlap by twice the
//! longest gene; detections only reach the gene store through `find_slot`.
//! Here the windows are searched in parallel and the detections replayed
//! through `find_slot` in the C program's order, which gives the same result
//! as a single-threaded run.

use rayon::prelude::*;
use std::collections::BTreeSet;

const ADENINE: u8 = 0;
const CYTOSINE: u8 = 1;
const GUANINE: u8 = 2;
const THYMINE: u8 = 3;
const AMBIG: u8 = 4;
const NOBASE: u8 = 5;

const ATBOND: f64 = 2.5;
const INACTIVE: f64 = 2.0e35;
const ASTEM2_EXT: i64 = 9;
const MINTSTEM_DIST: i64 = 17 + ASTEM2_EXT;
const MAXTSTEM_DIST: i64 = 26 + ASTEM2_EXT;
const MINCTRNALEN: i64 = 62;
const MAXCTRNALEN: i64 = 110;
const MINTRNALEN: i64 = MINCTRNALEN + 1;
const MAXTRNALEN: i64 = MAXCTRNALEN + ASTEM2_EXT;
/// ARAGORN's `MAXINTRONLEN`, the `-i` default.
pub const MAXINTRONLEN: i64 = 3000;
const MAXTAGDIST: i64 = 102;
const VARMIN: i64 = 3;
const VARDIFF: i64 = 23;
const NA: usize = MAXINTRONLEN as usize;
const ND: usize = 100;
const NTH: usize = 3000;
const NC: usize = 5000;
const LSEQ: i64 = 20000;
const PROXIMITY: i64 = 7 * MINCTRNALEN / 10;

// Search switches (ARAGORN defaults with -t -gc1 -l).
const TRNATHRESH: f64 = 132.0;
const TTSCANTHRESH: u32 = 4;
const TTARMTHRESH: f64 = 29.0;
const TDARMTHRESH: f64 = 26.0;
const TASTEMTHRESH: f64 = 7.5;
const TASCANTHRESH: u32 = 8;
const SP1MAX: i64 = 3;
const SP2MIN: i64 = 0;
const SP2MAX: i64 = 2;
const MININTRONLEN: i64 = 0;

static BP: [[i32; 6]; 6] = [
    [0, 0, 0, 1, 1, 0],
    [0, 0, 1, 0, 1, 0],
    [0, 1, 0, 1, 1, 0],
    [1, 0, 1, 0, 1, 0],
    [1, 1, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];
static VBP: [[u32; 6]; 6] = [
    [0, 0, 1, 4, 4, 0],
    [0, 0, 4, 0, 4, 0],
    [1, 4, 0, 2, 4, 0],
    [4, 0, 2, 0, 4, 0],
    [4, 4, 4, 4, 4, 0],
    [0, 0, 0, 0, 0, 0],
];
static BEM: [[f64; 6]; 6] = [
    [-2.144, -0.428, -2.144, ATBOND, 0.0, 0.0],
    [-0.428, -2.144, 3.000, -2.144, 0.0, 0.0],
    [-2.144, 3.000, -2.144, 1.286, 0.0, 0.0],
    [ATBOND, -2.144, 1.286, -0.428, 0.0, 0.0],
    [0.0; 6],
    [0.0; 6],
];
static DFEM: [[f64; 6]; 6] = [
    [-4.000, -4.000, -4.000, ATBOND, 0.0, 0.0],
    [-4.000, -4.000, 3.000, -4.000, 0.0, 0.0],
    [-4.000, 3.000, -4.000, 1.286, 0.0, 0.0],
    [ATBOND, -4.000, 1.286, -4.000, 0.0, 0.0],
    [0.0; 6],
    [0.0; 6],
];
static CBEM: [[f64; 6]; 6] = [
    [-1.072, -0.214, -1.072, 2.0 * ATBOND, 0.0, 0.0],
    [-0.214, -1.072, 6.000, -1.072, 0.0, 0.0],
    [-1.072, 6.000, -1.072, 3.400, 0.0, 0.0],
    [2.0 * ATBOND, -1.072, 3.400, -0.214, 0.0, 0.0],
    [0.0; 6],
    [0.0; 6],
];

/// Amino acid names in ARAGORN's order (Phe=0 ... Pyl=22).
static AANAME: [&str; 23] = [
    "Phe", "Val", "Leu", "Ile", "Cys", "Gly", "Arg", "Ser", "Ala", "Pro", "Thr", "Tyr", "Asp", "His",
    "Asn", "Met", "Trp", "Glu", "Gln", "Lys", "Stop", "SeC", "Pyl",
];
/// Standard genetic code (ARAGORN `aamap[1]`), indexed by anticodon
/// (p1<<4) + (p2<<2) + p3.
static AAMAP_STD: [usize; 64] = [
    0, 1, 2, 3, 4, 5, 6, 7, 7, 8, 9, 10, 11, 12, 13, 14, 2, 1, 2, 15, 16, 5, 6, 6, 7, 8, 9, 10, 22, 17,
    18, 19, 0, 1, 2, 3, 4, 5, 6, 7, 7, 8, 9, 10, 11, 12, 13, 14, 2, 1, 2, 3, 21, 5, 6, 6, 7, 8, 9, 10,
    20, 17, 18, 19,
];

/// One tRNA gene call.  Coordinates are 1-based within the searched
/// sequence, `start <= end`, as ARAGORN prints them (`c[start,end]` on the
/// minus strand).
#[derive(Clone, Debug, PartialEq)]
pub struct TrnaCall {
    pub start: i64,
    pub end: i64,
    pub minus: bool,
    /// e.g. "tRNA-Met(cat)"; ARAGORN's ambiguous forms are kept as printed
    pub name: String,
    /// ARAGORN's raw score (not normalised)
    pub energy: f64,
    pub intron_len: i64,
    /// The gene as ARAGORN prints it (intron included), upper case, '.' for
    /// non-ACGT, on the gene's strand.
    pub seq: String,
    /// Bases at the end of `seq` (and of the coordinates) that lie after the
    /// discriminator: ARAGORN extends the gene over a genomic C, CC or CCA
    /// there (0 to 3).  `gene()` is the sequence without them.
    pub tail_ext: usize,
}

impl TrnaCall {
    /// The gene up to its discriminator base: `seq` without the C / CC / CCA
    /// that ARAGORN adds when the genome happens to continue that way.  The
    /// mature tRNA's CCA is added after transcription, so this is the part
    /// the 3' end is built from whether or not the genome encodes the tail.
    pub fn gene(&self) -> &str {
        &self.seq[..self.seq.len() - self.tail_ext.min(self.seq.len())]
    }

    /// Amino acid ("Met") when the name is a plain `tRNA-Xxx(nnn)`.
    pub fn amino_acid(&self) -> Option<&str> {
        let rest = self.name.strip_prefix("tRNA-")?;
        let aa = rest.split('(').next()?;
        if AANAME[..20].contains(&aa) { Some(aa) } else { None }
    }
    /// Short type label, e.g. "Met(cat)".
    pub fn short_name(&self) -> &str {
        self.name.strip_prefix("tRNA-").unwrap_or(&self.name)
    }
}

/// Map a byte to ARAGORN's base code.  ARAGORN skips characters outside the
/// IUPAC alphabet; here they become AMBIG so coordinates stay aligned with the
/// input text.
#[inline]
fn code(b: u8) -> u8 {
    match b {
        b'A' | b'a' => ADENINE,
        b'C' | b'c' => CYTOSINE,
        b'G' | b'g' => GUANINE,
        b'T' | b't' | b'U' | b'u' => THYMINE,
        b'.' | b'`' => NOBASE,
        _ => AMBIG,
    }
}

#[derive(Clone, Copy, Default)]
struct Loop {
    pos: i64,
    stem: i64,
    lp: i64,
    energy: f64,
}

#[derive(Clone, Copy)]
struct DLoop {
    pos: i64,
    end: Option<i64>,
    stem: i64,
    lp: i64,
    energy: f64,
}

/// ARAGORN's `gene` for a tRNA, positions as indices into the window.
#[derive(Clone, Default)]
struct Gene {
    ps: i64,
    nbase: i64,
    astem1: i64,
    astem2: i64,
    spacer1: i64,
    spacer2: i64,
    dstem: i64,
    dloop: i64,
    cstem: i64,
    cloop: i64,
    intron: i64,
    nintron: i64,
    anticodon: i64,
    var: i64,
    tstem: i64,
    tloop: i64,
    energy: f64,
}

/// A detection in absolute (1-based sequence) coordinates, ready for the
/// gene store.
#[derive(Clone)]
struct Detection {
    start: i64,
    stop: i64,
    minus: bool,
    energy: f64,
    nintron: i64,
    name: String,
    seq: String,
    /// bases of `seq` after the discriminator (`aatail`'s extension)
    tail_ext: usize,
}

/// Base code at window index `i`.  The search only reads inside the window
/// (the C code uses raw pointers there); debug builds check it.
#[inline(always)]
fn b(s: &[u8], i: i64) -> usize {
    debug_assert!(i >= 0 && (i as usize) < s.len(), "window index {} outside 0..{}", i, s.len());
    // SAFETY: indices stay inside the padded window, as in ARAGORN.
    let v = unsafe { *s.get_unchecked(i as usize) } as usize;
    // SAFETY: windows hold base codes 0..=5 only (see `code`).
    unsafe { std::hint::assert_unchecked(v < 6) };
    v
}

fn find_tstems(s: &[u8], ls: i64, off: i64, hit: &mut Vec<Loop>) {
    static TEM: [u32; 6] = [0x0100, 0x0002, 0x2000, 0x0220, 0x0000, 0x0000];
    static A: [f64; 6] = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    static C: [f64; 6] = [0.0, 2.0, 0.0, 0.0, 0.0, 0.0];
    static G: [f64; 6] = [0.0, 0.0, 2.0, 0.0, 0.0, 0.0];
    static T: [f64; 6] = [0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
    hit.clear();
    let ss = off;
    let mut si = ss + 4 - 1;
    let sl = ls - off + 5 + 3;
    let mut r: u32 = TEM[b(s, si)];
    si += 1;
    r = (r >> 4) + TEM[b(s, si)];
    si += 1;
    r = (r >> 4) + TEM[b(s, si)];
    si += 1;
    while si < sl {
        r = (r >> 4) + TEM[b(s, si)];
        si += 1;
        let c = r & 0xF;
        if c < TTSCANTHRESH {
            continue;
        }
        let mut sb = si - 7;
        let mut sf = sb + 13;
        let ec = (3 * c) as f64;
        for tstem in 4..=5i64 {
            if sb < sl - 8 {
                let mut sc = sf;
                let mut sx = si - 2;
                for tloop in 5..=9i64 {
                    let penalty = if tloop > 7 { 3.0 * (tloop - tstem - 2) as f64 } else { 3.0 * (12 - tloop - tstem) as f64 };
                    let mut s1 = sb;
                    let mut s2 = sc;
                    let se = s1 + tstem;
                    s2 -= 1;
                    let mut energy = ec + BEM[b(s, se)][b(s, se + 4)] + BEM[b(s, s1)][b(s, s2)] - penalty;
                    s1 += 1;
                    while s1 < se {
                        s2 -= 1;
                        energy += BEM[b(s, s1)][b(s, s2)];
                        s1 += 1;
                    }
                    energy += G[b(s, sx)] + A[b(s, sx + 1)] + T[b(s, sx + 3)] + C[b(s, sx + 4)] + C[b(s, sx + 5)];
                    if energy >= TTARMTHRESH {
                        if hit.len() >= NTH {
                            return; // "Too many tstem hits"
                        }
                        hit.push(Loop { pos: sb, stem: tstem, lp: tloop, energy });
                    }
                    sx += 1;
                    sc += 1;
                }
            }
            sb -= 1;
            if sb < ss {
                break;
            }
            sf += 1;
        }
    }
}

fn find_astem5(s: &[u8], mut si: i64, sl: i64, astem3: i64, n3: i64, hit: &mut Vec<Loop>, nh: usize) {
    static A: [u32; 6] = [0, 0, 0, 2, 0, 0];
    static C: [u32; 6] = [0, 0, 2, 0, 0, 0];
    static G: [u32; 6] = [0, 2, 0, 1, 0, 0];
    static T: [u32; 6] = [2, 0, 1, 0, 0, 0];
    hit.clear();
    let sl = sl + n3;
    let mut tem = [0u32; 6];
    let mut se = astem3 + n3 - 1;
    tem[0] = A[b(s, se)];
    tem[1] = C[b(s, se)];
    tem[2] = G[b(s, se)];
    tem[3] = T[b(s, se)];
    se -= 1;
    while se >= astem3 {
        tem[0] = (tem[0] << 4) + A[b(s, se)];
        tem[1] = (tem[1] << 4) + C[b(s, se)];
        tem[2] = (tem[2] << 4) + G[b(s, se)];
        tem[3] = (tem[3] << 4) + T[b(s, se)];
        se -= 1;
    }
    let mut r = tem[b(s, si)];
    si += 1;
    let mut k = 1;
    loop {
        k += 1;
        if k >= n3 {
            break;
        }
        r = (r >> 4) + tem[b(s, si)];
        si += 1;
    }
    while si < sl {
        r = (r >> 4) + tem[b(s, si)];
        si += 1;
        if (r & 15) >= TASCANTHRESH {
            let mut s1 = astem3;
            let mut s2 = si;
            let se = s1 + n3;
            s2 -= 1;
            let mut energy = BEM[b(s, s1)][b(s, s2)];
            s1 += 1;
            while s1 < se {
                s2 -= 1;
                energy += BEM[b(s, s1)][b(s, s2)];
                s1 += 1;
            }
            if energy >= TASTEMTHRESH {
                if hit.len() >= nh {
                    return; // "Too many astem5 hits"
                }
                hit.push(Loop { pos: si - n3, stem: 0, lp: 0, energy });
            }
        }
    }
}

fn vloop_stability(s: &[u8], sb: i64, var: i64) -> f64 {
    static A: [u32; 6] = [0, 0, 0x100, 0x400, 0, 0];
    static C: [u32; 6] = [0, 0, 0x400, 0, 0, 0];
    static G: [u32; 6] = [0x100, 0x400, 0, 0x200, 0, 0];
    static T: [u32; 6] = [0x400, 0, 0x200, 0, 0, 0];
    let mut te = [0u32; 6];
    let mut e: u32 = 0;
    let mut vstem: i64 = 0;
    let sc = sb + 3;
    let mut se = sb + var - 2;
    let sf = se - 2;
    te[0] = A[b(s, se)];
    te[1] = C[b(s, se)];
    te[2] = G[b(s, se)];
    te[3] = T[b(s, se)];
    loop {
        se -= 1;
        if se <= sf {
            break;
        }
        te[0] = (te[0] >> 4) | A[b(s, se)];
        te[1] = (te[1] >> 4) | C[b(s, se)];
        te[2] = (te[2] >> 4) | G[b(s, se)];
        te[3] = (te[3] >> 4) | T[b(s, se)];
    }
    while se >= sc {
        te[0] = (te[0] >> 4) | A[b(s, se)];
        te[1] = (te[1] >> 4) | C[b(s, se)];
        te[2] = (te[2] >> 4) | G[b(s, se)];
        te[3] = (te[3] >> 4) | T[b(s, se)];
        let mut p = se - 5;
        let sd = se - 7;
        let mut m = te[b(s, p)];
        loop {
            p -= 1;
            if p <= sd {
                break;
            }
            m = (m >> 4) + te[b(s, p)];
        }
        while p >= sb {
            m = (m >> 4) + te[b(s, p)];
            let mut c = m & 0xf;
            if c >= 9 {
                let mut stem = 3;
                let mut lp = (se - p) - 3;
                let mut sen = se;
                let mut sn = p + 2;
                while lp >= 6 {
                    let cn = VBP[b(s, sen - 1)][b(s, sn + 1)];
                    if cn == 0 {
                        break;
                    }
                    c += cn;
                    stem += 1;
                    lp -= 2;
                    sen -= 1;
                    sn += 1;
                }
                if c > e {
                    e = c;
                    vstem = stem;
                }
            }
            p -= 1;
        }
        se -= 1;
    }
    if e > 0 { 3.0 * (vstem - 4) as f64 } else { -12.0 }
}

fn aatail(s: &[u8], p: i64) -> (i32, i64) {
    let mut score = 1;
    let mut e = 1;
    if s[(p + 1) as usize] == CYTOSINE {
        score += 1;
        e = 2;
        if s[(p + 2) as usize] == CYTOSINE {
            score += 1;
            e = 3;
            if s[(p + 3) as usize] == ADENINE {
                score += 1;
                e = 4;
            }
        }
    }
    (score, e)
}

fn cbase(c: u8) -> char {
    match c {
        ADENINE => 'A',
        CYTOSINE => 'C',
        GUANINE => 'G',
        THYMINE => 'T',
        _ => '.',
    }
}

fn aa(s: &[u8]) -> &'static str {
    if s[0] >= AMBIG || s[1] >= AMBIG || s[2] >= AMBIG {
        return "???";
    }
    AANAME[AAMAP_STD[((s[0] as usize) << 4) + ((s[1] as usize) << 2) + s[2] as usize]]
}

/// ARAGORN's `name(t, .., proc=1, ..)` for a cytosolic tRNA, from the spliced
/// sequence.
fn trna_name(seq: &[u8], anticodon: i64, cloop: i64) -> String {
    let i = anticodon as usize;
    let s = [seq[i - 1], seq[i], seq[i + 1], seq[i + 2], seq[i + 3]];
    let lc = |c: u8| cbase(c).to_ascii_lowercase();
    if cloop == 8 {
        format!("tRNA-?({}|{})({}{}{}{})", aa(&s[1..]), aa(&s[2..]), lc(s[1]), lc(s[2]), lc(s[3]), lc(s[4]))
    } else if cloop == 6 {
        format!("tRNA-?({}|{})({}{})", aa(&s[0..]), aa(&s[1..]), lc(s[1]), lc(s[2]))
    } else {
        format!("tRNA-{}({}{}{})", aa(&s[1..]), lc(s[1]), lc(s[2]), lc(s[3]))
    }
}

/// ARAGORN's `ti_genedetected` up to `find_slot`: finish the gene and turn it
/// into an absolute detection.  `wstart` is the 1-based sequence coordinate
/// of window index 0 (forward) or `start + len` (minus strand, as `sw->start`).
fn ti_genedetected(s: &[u8], te: &mut Gene, wstart: i64, minus: bool, out: &mut Vec<Detection>) {
    te.nbase = te.astem1 + te.spacer1 + te.spacer2 + 2 * te.dstem + te.dloop + 2 * te.cstem + te.cloop + te.var
        + 2 * te.tstem
        + te.tloop
        + te.astem2;
    let p = te.ps + te.nbase + te.nintron;
    let (mut as_, mut aext) = aatail(s, p);
    if te.astem1 == 7 && BP[b(s, te.ps - 1)][b(s, p)] != 0 {
        let (as8, aext8) = aatail(s, p + 1);
        if as8 >= as_ {
            te.ps -= 1;
            te.nbase += 2;
            te.anticodon += 1;
            if te.nintron > 0 {
                te.intron += 1;
            }
            te.astem1 = 8;
            te.astem2 = 8;
            as_ = as8;
            aext = aext8;
        }
    }
    let _ = as_;
    let nbasefext = te.nbase + ASTEM2_EXT;
    te.nbase += aext;
    // spliced copy (t->seq) for the name; printed copy for the sequence
    let ps = te.ps as usize;
    let spliced: Vec<u8> = if te.nintron == 0 {
        s[ps..ps + nbasefext as usize].to_vec()
    } else {
        let i = te.intron as usize;
        let mut v = s[ps..ps + i].to_vec();
        v.extend_from_slice(&s[ps + i + te.nintron as usize..ps + nbasefext as usize + te.nintron as usize]);
        v
    };
    let printed: String = s[ps..ps + (te.nbase + te.nintron) as usize].iter().map(|&c| cbase(c)).collect();
    let (start, stop) = if minus {
        let stop = wstart - te.ps - 1;
        (stop - te.nbase - te.nintron + 1, stop)
    } else {
        let start = te.ps + wstart;
        (start, start + te.nbase + te.nintron - 1)
    };
    out.push(Detection {
        start,
        stop,
        minus,
        energy: te.energy,
        nintron: te.nintron,
        name: trna_name(&spliced, te.anticodon, te.cloop),
        seq: printed,
        tail_ext: (aext - 1).max(0) as usize,
    });
}

struct Scratch {
    thit: Vec<Loop>,
    ahit: Vec<Loop>,
    dhit: Vec<DLoop>,
    chit: Vec<Loop>,
}

impl Scratch {
    fn new() -> Self {
        Scratch {
            thit: Vec::with_capacity(NTH),
            ahit: Vec::with_capacity(NA),
            dhit: Vec::with_capacity(ND),
            chit: Vec::with_capacity(NC),
        }
    }
}

/// ARAGORN's `tmioptimise` (tRNA part) on one window.
#[allow(clippy::too_many_lines)]
fn tmioptimise(s: &[u8], ls: i64, off: i64, maxintronlen: i64, wstart: i64, minus: bool, sc_: &mut Scratch, out: &mut Vec<Detection>) {
    static TT: [u32; 6] = [0x00, 0x00, 0x00, 0x11, 0x00, 0x00];
    static GG: [u32; 6] = [0x00, 0x00, 0x11, 0x00, 0x00, 0x00];
    static CA: [u32; 6] = [0, 0, 0, 2, 0, 0];
    static CC: [u32; 6] = [0, 0, 2, 0, 0, 0];
    static CG: [u32; 6] = [0, 2, 0, 1, 0, 0];
    static CT: [u32; 6] = [2, 0, 1, 0, 0, 0];
    static YIC: [i64; 9] = [1, 0, 0, 0, 0, 0, 0, 0, 0];
    static TIC: [i64; 9] = [1, 1, 0, 0, 0, 0, 0, 0, 0];
    static A1IC: [i64; 9] = [1, 1, 1, 0, 0, 0, 0, 0, 0];
    static A2IC: [i64; 9] = [1, 1, 1, 1, 0, 0, 0, 0, 0];
    static A3IC: [i64; 9] = [1, 1, 1, 1, 1, 0, 0, 0, 0];
    static RIC: [i64; 9] = [1, 1, 1, 1, 1, 1, 0, 0, 0];
    // D-loop G-run windows: ige[GOFFB[dloop]..=GOFFE[dloop]] with
    // GOFFB = [.., 1, 2, 2, 2, 2, 2, 2, 2] and GOFFE = [.., 2, 3, 4, 4, 5, 6, 6, 6]
    // for dloop 4..=11 (unrolled in the D-stem search).
    static CY: [i32; 6] = [0, 1, 0, 1, 0, 0];
    static CR: [i32; 6] = [1, 0, 1, 0, 0, 0];
    const ILW: f64 = 0.002;
    static G: [f64; 6] = [0.0, 0.0, 6.0, 0.0, 0.0, 0.0];
    static T: [f64; 6] = [0.0, 0.0, 0.0, 7.0, 0.0, 0.0];
    static Y: [f64; 6] = [0.0, 3.0, 0.0, 3.0, 0.0, 0.0];
    static R: [f64; 6] = [2.0, 0.0, 2.0, 0.0, 0.0, 0.0];
    static YP: [f64; 6] = [0.0, 3.0, 0.0, 3.0, 0.0, 0.0];
    static RP: [f64; 6] = [2.0, 0.0, 2.0, 0.0, 0.0, 0.0];
    static RI: [f64; 6] = [0.1, 0.0, 0.05, 0.0, 0.0, 0.0];
    static GI: [f64; 6] = [0.0, 0.0, 0.1, 0.0, 0.0, 0.0];
    static YI: [f64; 6] = [0.0, 0.1, 0.0, 0.1, 0.0, 0.0];
    static AI: [f64; 6] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    static GC: [f64; 6] = [0.0, 1.5, 6.0, 0.0, 0.0, 0.0];
    static G3: [f64; 6] = [0.0, 6.0, 12.0, 12.0, 0.0, 0.0];
    static DR: [f64; 6] = [6.0, 0.0, 6.0, 0.0, 0.0, 0.0];
    static RH: [f64; 6] = [3.0, 0.0, 3.0, 0.0, 0.0, 0.0];
    static AGT: [f64; 6] = [6.0, 0.0, 6.0, 6.0, 0.0, 0.0];
    static DT: [f64; 6] = [0.0, 0.0, 0.0, 6.0, 0.0, 0.0];
    static DBEM: [[f64; 6]; 6] = BEM;

    let ethresh = TRNATHRESH;
    let tmindist = MINTRNALEN + MININTRONLEN - MAXTSTEM_DIST;
    let tmaxdist = MAXTRNALEN + maxintronlen - MINTSTEM_DIST;
    let mindist = tmindist;
    let maxdist = tmaxdist;
    let Scratch { thit, ahit, dhit, chit } = sc_;
    find_tstems(s, ls, off, thit);
    let mut t = Gene { astem1: 7, astem2: 7, spacer1: 2, spacer2: 1, dstem: 3, dloop: 9, cstem: 5, cloop: 7, var: 15, tstem: 5, tloop: 7, ..Default::default() };
    let mut te = Gene::default();
    for th in thit.iter() {
        let tpos = th.pos;
        t.tloop = th.lp;
        t.tstem = th.stem;
        let tfold = b(s, tpos - 1);
        let tloopfold = tpos + t.tstem + 1;
        let tarm = 2 * t.tstem + t.tloop;
        let tend = tpos + tarm;
        let tmv = tpos - VARMIN;
        let mut flag = false;
        te.energy = ethresh;
        let the = th.energy;
        find_astem5(s, tpos - maxdist, tpos - mindist, tend, 7, ahit, NA);
        for ah in ahit.iter() {
            let apos = ah.pos;
            if apos < tpos - tmaxdist {
                continue;
            }
            if apos > tpos - tmindist {
                break;
            }
            let he = the + ah.energy;

            // ---- find dstems ----
            dhit.clear();
            let mut sc = apos + 8;
            let mut energyf = DFEM[b(s, sc + 5)][tfold];
            let sl = sc + SP1MAX;
            'dstems: while sc < sl {
                let energy2 = DT[b(s, sc - 2)] + RH[b(s, sc - 1)] + GC[b(s, sc)] + DFEM[b(s, sc - 2)][b(s, sc + 4)];
                let energyf6 = DFEM[b(s, sc + 6)][tfold];
                for dstem in 3..=4i64 {
                    let sd = sc + dstem;
                    let dloop = 3;
                    let mut se = sd + dloop;
                    let mut energy = energy2 + 6.0 + DR[b(s, se - 1)] + energyf;
                    if dstem == 3 && energyf < 0.0 {
                        energyf = energyf6;
                    }
                    se += dstem;
                    let (mut s1, mut s2, sf) = (sc, se, sc + dstem);
                    while s1 < sf {
                        s2 -= 1;
                        energy += DBEM[b(s, s1)][b(s, s2)];
                        s1 += 1;
                    }
                    if energy >= TDARMTHRESH {
                        if dhit.len() >= ND {
                            break 'dstems; // "Too many D-stem hits"
                        }
                        dhit.push(DLoop { pos: sc, end: Some(se), lp: dloop, stem: dstem, energy });
                    }
                    let mut ige = [0u32; 7];
                    let mut sg1 = sd + 1;
                    let sg2 = sd + 6;
                    let mut q = GG[b(s, sg1)];
                    sg1 += 1;
                    ige[1] = q & 3;
                    let mut j = 2;
                    while sg1 <= sg2 {
                        q = (q >> 4) + GG[b(s, sg1)];
                        sg1 += 1;
                        ige[j] = q & 3;
                        j += 1;
                    }
                    // OR of ige[GOFFB[dloop]..=GOFFE[dloop]], unrolled
                    let o23 = ige[2] | ige[3];
                    let o234 = o23 | ige[4];
                    let o2345 = o234 | ige[5];
                    let gor = [ige[1] | ige[2], o23, o234, o234, o2345, o2345 | ige[6], o2345 | ige[6], o2345 | ige[6]];
                    for dloop in 4..=11i64 {
                        let genergy = G3[(gor[(dloop - 4) as usize] & 3) as usize];
                        let mut se = sd + dloop;
                        let mut energy = energy2 + genergy + DR[b(s, se - 1)] + energyf;
                        se += dstem;
                        let (mut s1, mut s2, sf) = (sc, se, sc + dstem);
                        while s1 < sf {
                            s2 -= 1;
                            energy += DBEM[b(s, s1)][b(s, s2)];
                            s1 += 1;
                        }
                        if energy >= TDARMTHRESH {
                            if dhit.len() >= ND {
                                break 'dstems;
                            }
                            dhit.push(DLoop { pos: sc, end: Some(se), lp: dloop, stem: dstem, energy });
                        }
                    }
                }
                // 6 bp D-stem
                {
                    let (mut s1, mut s2) = (sc, sc + 16);
                    let mut sd = sc + 6;
                    s2 -= 1;
                    let mut j = BP[b(s, s1)][b(s, s2)];
                    loop {
                        s1 += 1;
                        if s1 >= sd {
                            break;
                        }
                        s2 -= 1;
                        j += BP[b(s, s1)][b(s, s2)];
                    }
                    if j >= 6 {
                        let mut energy = DT[b(s, sc - 1)] + RH[b(s, sc)] + GC[b(s, sc + 1)] + energyf6;
                        sd += 1;
                        energy += G[b(s, sd)];
                        sd += 1;
                        energy += G[b(s, sd)];
                        sd += 1;
                        energy += AGT[b(s, sd)] + DFEM[b(s, sc - 1)][b(s, sc + 4)];
                        sd += 7;
                        let (mut s1, mut s2, sf) = (sc, sd, sc + 6);
                        while s1 < sf {
                            s2 -= 1;
                            energy += DBEM[b(s, s1)][b(s, s2)];
                            s1 += 1;
                        }
                        if energy >= TDARMTHRESH {
                            if dhit.len() >= ND {
                                break 'dstems;
                            }
                            dhit.push(DLoop { pos: sc, end: Some(sd), lp: 4, stem: 6, energy });
                        }
                    }
                }
                // 7 bp D-stem
                {
                    let (mut s1, mut s2) = (sc, sc + 18);
                    let mut sd = sc + 7;
                    s2 -= 1;
                    let mut j = BP[b(s, s1)][b(s, s2)];
                    loop {
                        s1 += 1;
                        if s1 >= sd {
                            break;
                        }
                        s2 -= 1;
                        j += BP[b(s, s1)][b(s, s2)];
                    }
                    if j >= 7 {
                        let mut energy = energy2 + DFEM[b(s, sc + 7)][tfold];
                        sd += 1;
                        energy += G[b(s, sd)];
                        sd += 1;
                        energy += G[b(s, sd)];
                        sd += 1;
                        energy += AGT[b(s, sd)];
                        sd += 8;
                        let (mut s1, mut s2, sf) = (sc, sd, sc + 7);
                        while s1 < sf {
                            s2 -= 1;
                            energy += DBEM[b(s, s1)][b(s, s2)];
                            s1 += 1;
                        }
                        if energy >= TDARMTHRESH {
                            if dhit.len() >= ND {
                                break 'dstems;
                            }
                            dhit.push(DLoop { pos: sc, end: Some(sd), lp: 4, stem: 7, energy });
                        }
                    }
                }
                energyf = energyf6;
                sc += 1;
            }
            // ---- end of find dstems ----

            let ndh = dhit.len();
            for nd1 in (0..ndh).rev() {
                let dstem = dhit[nd1].stem;
                let dpos = dhit[nd1].pos;
                if dpos - apos < 9 {
                    dhit[nd1].energy -= 3.0;
                }
                if s[tloopfold as usize] == GUANINE {
                    let mut sb = dpos + dstem + 2;
                    let mut sc = sb;
                    let se = sb + dhit[nd1].lp - 3;
                    let mut r = TT[b(s, sb)];
                    sb += 1;
                    while sb < se {
                        r = (r >> 4) + TT[b(s, sb)];
                        sb += 1;
                        if r & 2 != 0 {
                            dhit[nd1].energy += 10.0;
                            break;
                        }
                    }
                    r = GG[b(s, sc)];
                    sc += 1;
                    while sc < se {
                        r = (r >> 4) + GG[b(s, sc)];
                        sc += 1;
                        if r & 2 != 0 {
                            dhit[nd1].energy -= 12.0;
                            break;
                        }
                    }
                }
            }
            for nd1 in (0..ndh).rev() {
                let Some(cpos) = dhit[nd1].end else { continue };
                let mut denergy = dhit[nd1].energy;
                let mut ndx = nd1;
                for nd2 in (0..nd1).rev() {
                    if dhit[nd2].end != Some(cpos) {
                        continue;
                    }
                    let e = dhit[nd2].energy;
                    if e > denergy {
                        denergy = e;
                        dhit[ndx].end = None;
                        ndx = nd2;
                    }
                }
            }
            let mut cposmin: Option<i64> = None;
            let mut cposmax: i64 = 0;
            let mut nd1 = ndh as i64;
            loop {
                nd1 -= 1;
                if nd1 < 0 {
                    break;
                }
                if let Some(e) = dhit[nd1 as usize].end {
                    cposmin = Some(e);
                    cposmax = e;
                    break;
                }
            }
            let Some(mut cmin) = cposmin else { continue };
            let mut nd2 = nd1;
            loop {
                nd2 -= 1;
                if nd2 < 0 {
                    break;
                }
                if let Some(c) = dhit[nd2 as usize].end {
                    if c < cmin {
                        cmin = c;
                    }
                    if c > cposmax {
                        cposmax = c;
                    }
                }
            }
            let mut cpos = cmin + SP2MIN;
            while cpos <= cposmax + SP2MAX {
                let mut denergy = -INACTIVE;
                let mut ndx: i64 = -1;
                for nd1 in (0..ndh).rev() {
                    let Some(end) = dhit[nd1].end else { continue };
                    if end + SP2MAX < cpos {
                        continue;
                    }
                    if end + SP2MIN > cpos {
                        continue;
                    }
                    let e = dhit[nd1].energy;
                    if e > denergy {
                        denergy = e;
                        ndx = nd1 as i64;
                    }
                }
                if ndx < 0 {
                    cpos += 1;
                    continue;
                }
                let ndx = ndx as usize;
                denergy += he;
                if denergy < te.energy - 49.0 {
                    cpos += 1;
                    continue;
                }

                // ---- find cstems ----
                chit.clear();
                {
                    let mut ct = [0u32; 6];
                    let mut si = cpos;
                    let sc = cpos + 5;
                    let mut se = cpos + 4;
                    ct[0] = CA[b(s, se)];
                    ct[1] = CC[b(s, se)];
                    ct[2] = CG[b(s, se)];
                    ct[3] = CT[b(s, se)];
                    se -= 1;
                    while se >= cpos {
                        ct[0] = (ct[0] << 4) + CA[b(s, se)];
                        ct[1] = (ct[1] << 4) + CC[b(s, se)];
                        ct[2] = (ct[2] << 4) + CG[b(s, se)];
                        ct[3] = (ct[3] << 4) + CT[b(s, se)];
                        se -= 1;
                    }
                    si += 11;
                    let se = tmv - VARDIFF - 5;
                    if si < se {
                        si = se;
                    }
                    let mut r = ct[b(s, si)];
                    si += 1;
                    for _ in 0..3 {
                        r = (r >> 4) + ct[b(s, si)];
                        si += 1;
                    }
                    while si < tmv {
                        r = (r >> 4) + ct[b(s, si)];
                        si += 1;
                        if (r & 0xf) >= 5 {
                            if chit.len() >= NC {
                                break; // "Too many cstem hits"
                            }
                            let mut h = Loop { pos: si, stem: 5, lp: si - sc - 5, energy: 0.0 };
                            if h.lp == 9
                                && BP[b(s, sc)][b(s, si - 6)] != 0
                                && CY[b(s, sc + 2)] != 0
                                && CR[b(s, sc + 6)] != 0
                                && CY[b(s, sc + 1)] != 0
                            {
                                h.stem = 6;
                                h.lp = 7;
                            }
                            let (mut s1, mut s2) = (cpos, si);
                            let se = s1 + h.stem;
                            s2 -= 1;
                            h.energy = CBEM[b(s, s1)][b(s, s2)];
                            s1 += 1;
                            while s1 < se {
                                s2 -= 1;
                                h.energy += CBEM[b(s, s1)][b(s, s2)];
                                s1 += 1;
                            }
                            chit.push(h);
                        }
                    }
                }
                // ---- end of find cstems ----

                for ch in chit.iter() {
                    let mut energy = denergy + ch.energy;
                    if energy < te.energy - 19.0 {
                        continue;
                    }
                    let cend = ch.pos;
                    t.var = tpos - cend;
                    t.cloop = ch.lp;
                    t.cstem = ch.stem;
                    let mut intron: i64 = 0;
                    let sb;
                    if t.cloop < 9 {
                        if MININTRONLEN > 0 {
                            continue;
                        }
                        t.nintron = 0;
                        if t.var > 17 {
                            energy += vloop_stability(s, cend, t.var);
                        }
                        sb = cpos + t.cstem;
                        energy += T[b(s, sb + 1)] + Y[b(s, sb)] + R[b(s, sb + 5)] - 0.05 * t.var as f64
                            - if t.cloop == 7 { 0.0 } else { 6.0 };
                    } else {
                        t.nintron = t.cloop - 7;
                        if t.nintron > maxintronlen {
                            continue;
                        }
                        if t.nintron < MININTRONLEN {
                            continue;
                        }
                        if t.var > 17 {
                            energy += vloop_stability(s, cend, t.var);
                        }
                        if energy < te.energy - 9.0 {
                            continue;
                        }
                        t.cloop = 7;
                        sb = cpos + t.cstem;
                        let se = sb + t.nintron;
                        let mut cenergy = YP[b(s, se)] + T[b(s, se + 1)] + RP[b(s, se + 5)];
                        let mut ienergy = cenergy + RI[b(s, sb)] + GI[b(s, se - 1)] + AI[b(s, se - 2)] * YI[b(s, se - 1)];
                        for j in 1..=7i64 {
                            let si = se + j - 1;
                            let ju = j as usize;
                            let ec = YP[b(s, sb + YIC[ju] * t.nintron)]
                                + T[b(s, sb + TIC[ju] * t.nintron + 1)]
                                + RP[b(s, sb + RIC[ju] * t.nintron + 5)];
                            let mut e = ec + RI[b(s, sb + j)] + GI[b(s, si)] + AI[b(s, si - 1)] * YI[b(s, si)];
                            if j == 6 {
                                e += 0.01;
                            }
                            if e > ienergy {
                                ienergy = e;
                                cenergy = ec;
                                intron = j;
                            }
                        }
                        energy += cenergy - 10.0 - ILW * (t.nintron as f64 + 1.1 * t.var as f64);
                        if t.nintron >= 130 {
                            let si = se + intron;
                            let j = s[(si - 1) as usize];
                            if j != GUANINE {
                                if s[(si - 2) as usize] != ADENINE {
                                    energy -= 4.0;
                                }
                                if j != CYTOSINE && j != THYMINE {
                                    energy -= 8.0;
                                }
                            }
                        }
                    }
                    let dstem = dhit[ndx].stem;
                    let dpos = dhit[ndx].pos;
                    let iu = intron as usize;
                    if dstem >= 6 {
                        if s[(sb + 2 + A1IC[iu] * t.nintron) as usize] != THYMINE {
                            continue;
                        }
                        if s[(sb + 3 + A2IC[iu] * t.nintron) as usize] != CYTOSINE {
                            continue;
                        }
                        if s[(sb + 4 + A3IC[iu] * t.nintron) as usize] != ADENINE {
                            continue;
                        }
                        energy += 3.0;
                    } else if s[(dpos - 1) as usize] & 5 == 0 {
                        let mut i = 0;
                        let mut si = cend;
                        let se = cend + 4;
                        while si < se {
                            let v = s[si as usize];
                            si += 1;
                            if v & 5 == 0 {
                                i += 1;
                                if i >= 2 {
                                    energy += 3.0;
                                    break;
                                }
                            } else {
                                i = 0;
                            }
                        }
                    }
                    if t.cstem >= 6
                        && s[(sb + 2 + A1IC[iu] * t.nintron) as usize] == CYTOSINE
                        && s[(sb + 3 + A2IC[iu] * t.nintron) as usize] == THYMINE
                        && s[(sb + 4 + A3IC[iu] * t.nintron) as usize] == ADENINE
                    {
                        energy += 4.0;
                    }
                    if energy < ethresh {
                        continue;
                    }
                    t.energy = energy;
                    t.dstem = dstem;
                    t.astem1 = if t.dstem < 6 { 7 } else if t.tstem < 5 { 9 } else { 8 };
                    t.astem2 = t.astem1;
                    t.ps = apos + 7 - t.astem1;
                    t.nbase = (tend - t.ps) + t.astem2;
                    t.dloop = dhit[ndx].lp;
                    t.spacer1 = (dpos - apos) - 7;
                    t.spacer2 = cpos - dhit[ndx].end.unwrap_or(cpos);
                    let j = (cpos - t.ps) + t.cstem;
                    t.anticodon = j + 2;
                    if t.nintron > 0 {
                        t.intron = j + intron;
                        if t.nbase + t.nintron > MAXTRNALEN {
                            ti_genedetected(s, &mut t, wstart, minus, out);
                            continue;
                        }
                    }
                    if energy < te.energy {
                        continue;
                    }
                    flag = true;
                    te = t.clone();
                }
                cpos += 1;
            }
        }
        if flag {
            ti_genedetected(s, &mut te, wstart, minus, out);
        }
    }
}

/// Search one window (both strands), returning its detections in ARAGORN's
/// order: forward strand first, then the reverse complement.
fn search_window(w: &[u8], off: i64, maxintronlen: i64, wstart: i64) -> Vec<Detection> {
    let mut out = Vec::new();
    let mut scratch = Scratch::new();
    let len = w.len() as i64;
    tmioptimise(w, len, off, maxintronlen, wstart, false, &mut scratch, &mut out);
    let c: Vec<u8> = w.iter().rev().map(|&x| if x <= THYMINE { THYMINE - x } else { x }).collect();
    tmioptimise(&c, len, off, maxintronlen, wstart + len, true, &mut scratch, &mut out);
    out
}

/// The gene store with ARAGORN's `find_slot` / `nearest_trna_gene` logic.
struct Store {
    genes: Vec<Detection>,
    by_start: BTreeSet<(i64, usize)>,
    maxintronlen: i64,
    maxlen: i64,
}

impl Store {
    fn new(maxintronlen: i64) -> Self {
        Store { genes: Vec::new(), by_start: BTreeSet::new(), maxintronlen, maxlen: MAXTRNALEN + maxintronlen + 16 }
    }

    /// Index of the overlapping same-strand gene with the lowest energy
    /// (first in slot order on ties), as `nearest_trna_gene` picks it.
    fn nearest(&self, t: &Detection) -> Option<usize> {
        let (a, b) = (t.start, t.stop);
        let thresh = b - a;
        let mut best: Option<(f64, usize)> = None;
        for &(_, i) in self.by_start.range((a - self.maxlen, 0)..=(b, usize::MAX)) {
            let g = &self.genes[i];
            let (c, e) = (g.start, g.stop);
            if a > e || b < c || g.minus != t.minus {
                continue;
            }
            if self.maxintronlen > 0 {
                let ilength = e - c;
                if 2 * thresh > 5 * ilength || 2 * ilength > 5 * thresh {
                    continue;
                }
            }
            let score = if a >= c { if b >= e { e - a } else { thresh } } else if b >= e { e - c } else { b - c };
            if score >= PROXIMITY {
                let better = match best {
                    None => g.energy < INACTIVE,
                    Some((be, bi)) => g.energy < be || (g.energy == be && i < bi),
                };
                if better {
                    best = Some((g.energy, i));
                }
            }
        }
        best.map(|(_, i)| i)
    }

    fn add(&mut self, t: Detection) {
        match self.nearest(&t) {
            Some(i) => {
                if t.energy <= self.genes[i].energy {
                    return;
                }
                self.by_start.remove(&(self.genes[i].start, i));
                self.by_start.insert((t.start, i));
                self.genes[i] = t;
            }
            None => {
                let i = self.genes.len();
                self.by_start.insert((t.start, i));
                self.genes.push(t);
            }
        }
    }

    /// `remove_overlapping_trna` (run when introns are searched), then
    /// `gene_sort`: by start, longer first on equal starts.
    fn finish(mut self) -> Vec<Detection> {
        let n = self.genes.len();
        let mut dead = vec![false; n];
        let ov = |a: i64, b: i64, c: i64, e: i64| -> i64 {
            let len = b - a;
            if a >= c { if b >= e { e - a } else { len } } else if b >= e { len } else { b - c }
        };
        let near = |g: &Detection| -> Vec<usize> {
            self.by_start.range((g.start - self.maxlen, 0)..=(g.stop, usize::MAX)).map(|&(_, i)| i).collect()
        };
        if self.maxintronlen > 0 {
            for k in 0..n {
                if dead[k] || self.genes[k].nintron <= 0 {
                    continue;
                }
                let (a, b) = (self.genes[k].start, self.genes[k].stop);
                let mut cand = near(&self.genes[k]);
                cand.sort_unstable();
                for i in cand {
                    if i == k || dead[i] || self.genes[i].minus != self.genes[k].minus {
                        continue;
                    }
                    let (c, e) = (self.genes[i].start, self.genes[i].stop);
                    if ov(a, b, c, e) >= PROXIMITY && self.genes[k].energy < self.genes[i].energy {
                        dead[k] = true;
                        break;
                    }
                }
            }
            for k in 0..n.saturating_sub(1) {
                if dead[k] {
                    continue;
                }
                let (a, b) = (self.genes[k].start, self.genes[k].stop);
                let mut cand: Vec<usize> = near(&self.genes[k]).into_iter().filter(|&i| i > k).collect();
                cand.sort_unstable();
                for i in cand {
                    if dead[i] || self.genes[i].minus != self.genes[k].minus {
                        continue;
                    }
                    let (c, e) = (self.genes[i].start, self.genes[i].stop);
                    if ov(a, b, c, e) >= PROXIMITY {
                        if self.genes[k].energy < self.genes[i].energy {
                            dead[k] = true;
                            break;
                        } else if self.genes[i].energy < self.genes[k].energy {
                            dead[i] = true;
                        }
                    }
                }
            }
        }
        let mut keep: Vec<Detection> =
            self.genes.drain(..).zip(dead).filter(|(_, d)| !d).map(|(g, _)| g).collect();
        keep.sort_by(|x, y| x.start.cmp(&y.start).then(y.stop.cmp(&x.stop)));
        keep
    }
}

/// Window layout of `iopt_fastafile` for a linear sequence of `n` bases:
/// (offset into the padded sequence, length) per window.
fn windows(n: i64, rewind: i64) -> Vec<(i64, i64)> {
    let drewind = 2 * rewind;
    let real_end = rewind + n;
    let pend = real_end + rewind;
    let mut out = Vec::new();
    let mut o = 0;
    loop {
        if o + LSEQ <= real_end {
            out.push((o, LSEQ));
            o += LSEQ - drewind;
        } else {
            out.push((o, pend - o));
            break;
        }
    }
    out
}

/// Search one linear sequence for tRNA genes, as `aragorn -t -gc1 -l`
/// (`maxintronlen` 0) or `aragorn -t -gc1 -l -i` (`maxintronlen` 3000).
/// Windows are searched in parallel.
pub fn scan_sequence(seq: &[u8], maxintronlen: i64) -> Vec<TrnaCall> {
    let n = seq.len() as i64;
    if n == 0 {
        return Vec::new();
    }
    let rewind = (MAXTAGDIST + 20).max(MAXTRNALEN + maxintronlen);
    let mut padded = vec![NOBASE; (n + 2 * rewind) as usize];
    for (i, &c) in seq.iter().enumerate() {
        padded[rewind as usize + i] = code(c);
    }
    let start0 = 1 - rewind;
    let wins = windows(n, rewind);
    let dets: Vec<Vec<Detection>> = wins
        .par_iter()
        .map(|&(o, len)| search_window(&padded[o as usize..(o + len) as usize], rewind, maxintronlen, start0 + o))
        .collect();
    let mut store = Store::new(maxintronlen);
    for d in dets.into_iter().flatten() {
        store.add(d);
    }
    store
        .finish()
        .into_iter()
        .map(|d| TrnaCall {
            // ARAGORN prints a linear start <= 0 one lower
            start: if d.start <= 0 { d.start - 1 } else { d.start },
            end: d.stop,
            minus: d.minus,
            name: d.name,
            energy: d.energy,
            intron_len: d.nintron,
            seq: d.seq,
            tail_ext: d.tail_ext,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // Arabidopsis chloroplast (NC_000932.1:51996-52188) tRNA-Met(cat) and
    // chromosome 1 (NC_003070.9:18308843-18309047) tRNA-Tyr(gta) with a 12 bp
    // intron, each with 60 bp of flank; ARAGORN 1.2.41 -t -i -gc1 -l calls
    // tRNA-Met(cat) [61,133] and tRNA-Tyr(gta) [61,145].
    const MET: &str = "GATTTATAATCCTATCGATGTGATAGGTATCCCCTTTCCTTCTCTTTACAATGATAAATAACCTACTTAACTCAGTGGTTAGAGTATTGCTTTCATACGGCAGGAGTCATTGGTTCAAATCCAATAGTAGGTATAACTTATTAGACACCATGATCAATGGTGTCTAATAAGTTTTTGTAACCAGCTTTTTTTA";
    const TYR: &str = "TACACTTACAGCTCATCCCTCTTTTAGAATTTAAATAATTTTCTGAAATATTAAGCACAACCGACCTTAGCTCAGTTGGTAGAGCGGAAGACTGTAGTGTTTGCAGACAATCTTTAGGTCGCTGGTTCGATTCCGGCAGGTCGGATTTTTAATGTTTTGGTGGTTTTTAAGGATCCTTCTAAAAACAAAGCCCAAGAGTTCTCTT";

    fn rc(s: &str) -> String {
        s.bytes().rev().map(|b| match b { b'A' => 'T', b'C' => 'G', b'G' => 'C', _ => 'A' }).collect()
    }

    #[test]
    fn matches_aragorn_calls() {
        let c = scan_sequence(MET.as_bytes(), MAXINTRONLEN);
        assert_eq!(c.len(), 1);
        assert_eq!((c[0].start, c[0].end, c[0].minus, c[0].name.as_str(), c[0].intron_len), (61, 133, false, "tRNA-Met(cat)", 0));
        assert_eq!(c[0].seq, MET[60..133]);
        let c = scan_sequence(TYR.as_bytes(), MAXINTRONLEN);
        assert_eq!(c.len(), 1);
        assert_eq!((c[0].start, c[0].end, c[0].name.as_str(), c[0].intron_len), (61, 145, "tRNA-Tyr(gta)", 12));
        // without the intron search ARAGORN reads the intron as an 8-base
        // anticodon loop: tRNA-?(Tyr|Leu)(gtag), no intron
        let c = scan_sequence(TYR.as_bytes(), 0);
        assert_eq!((c[0].name.as_str(), c[0].intron_len), ("tRNA-?(Tyr|Leu)(gtag)", 0));
        assert_eq!(c[0].amino_acid(), None);
    }

    #[test]
    fn genomic_c_after_the_gene_is_not_gene() {
        // the three bases after the discriminator (TAA) changed to CAA, CCT and CCA: ARAGORN
        // reports the gene 1, 2 and 3 bases longer; gene() stays the 73 bases
        assert_eq!(&MET[133..136], "TAA");
        for (k, after) in ["CAA", "CCT", "CCA"].iter().enumerate() {
            let s = format!("{}{}{}", &MET[..133], after, &MET[136..]);
            let c = scan_sequence(s.as_bytes(), MAXINTRONLEN);
            assert_eq!(c.len(), 1);
            assert_eq!((c[0].start, c[0].end, c[0].tail_ext), (61, 134 + k as i64, k + 1), "{after}");
            assert_eq!(c[0].seq, format!("{}{}", &MET[60..133], &after[..k + 1]));
            assert_eq!(c[0].gene(), &MET[60..133]);
        }
        let c = scan_sequence(MET.as_bytes(), MAXINTRONLEN);
        assert_eq!((c[0].tail_ext, c[0].gene()), (0, &MET[60..133]));
    }

    #[test]
    fn minus_strand_coordinates() {
        let r = rc(MET);
        let c = scan_sequence(r.as_bytes(), MAXINTRONLEN);
        assert_eq!(c.len(), 1);
        let n = MET.len() as i64;
        assert_eq!((c[0].start, c[0].end, c[0].minus), (n - 133 + 1, n - 61 + 1, true));
        assert_eq!(c[0].seq, MET[60..133]);
    }

    #[test]
    fn windows_cover_sequence() {
        let w = windows(100_000, 3119);
        assert_eq!(w[0], (0, LSEQ));
        assert_eq!(w[1].0, LSEQ - 2 * 3119);
        let last = w.last().unwrap();
        assert_eq!(last.0 + last.1, 100_000 + 2 * 3119);
        assert_eq!(windows(10, 122), vec![(0, 10 + 2 * 122)]);
    }
}
