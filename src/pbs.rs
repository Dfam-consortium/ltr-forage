//! Primer binding site (PBS) annotation from the genome's own tRNA genes.
//!
//! Reverse transcription of an LTR retrotransposon is primed by a host tRNA
//! whose 3' end (with the CCA tail) pairs with the PBS just inside the 5' LTR.
//! Probes are built from the tRNA calls of `crate::trna`: for each gene of
//! 70-120 bp with a standard amino acid, the reverse complement of its last
//! 15 nt + CCA (18 nt, starting TGG).
//!
//! Primary test: the probe's first 12 nt (the tRNA's 3'-terminal 12 nt; the
//! TGG must match exactly) with at most one mismatch, starting within -5..+15
//! bp of the first base after the 5' LTR.  Many families pair over only these
//! 12-13 nt (rice iMet users read TGGTATCAGAGC and then diverge), which is why
//! the full 18-mer is a secondary column: its best edit distance (<= 2) over
//! the same window, seeded by its three 6-mers.  On rice, kept elements that
//! match a curated family have an exact 12-nt PBS 56% of the time and random
//! genome windows 0.1%.

use std::collections::{HashMap, HashSet};

use crate::trna::TrnaCall;

pub const K: usize = 12;
const PROBE: usize = 18;
const MAX_EDITS: usize = 2;
/// Window: probe start from PRE bp before to POST bp after the LTR boundary.
const PRE: usize = 5;
const POST: usize = 15;

/// Best PBS for one orientation or one element.
#[derive(Clone, Debug, PartialEq)]
pub struct PbsHit {
    /// Mismatches of the 12-nt probe (0 or 1)
    pub mm12: u8,
    /// tRNA types sharing the probe, e.g. ["Met(cat)"]
    pub trnas: Vec<String>,
    /// Probe start relative to the first base after the LTR
    pub offset: i64,
    /// The matched probe's 12 nt (its tRNA genes' 3'-terminal 12 nt with
    /// CCA, reverse complemented); identifies the probe group
    pub probe12: Vec<u8>,
    /// true when the PBS lies before the right LTR (element on the minus strand)
    pub minus: bool,
}

pub struct Probes {
    /// 18-nt probe -> tRNA types
    full: HashMap<Vec<u8>, Vec<String>>,
    /// 12-mer (exact or one-mismatch variant) -> (mismatches, tRNA types, probe 12-mer)
    short: HashMap<Vec<u8>, Vec<(u8, Vec<String>, Vec<u8>)>>,
    /// 6-mer seed -> (18-nt probe, offset in probe)
    seed: HashMap<Vec<u8>, Vec<(Vec<u8>, usize)>>,
}

fn rc(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

/// The 18-nt PBS probe of a tRNA call: reverse complement of its last 15 nt
/// + CCA, for genes of 70-120 bp with a standard amino acid; None otherwise.
pub fn probe18(c: &TrnaCall) -> Option<Vec<u8>> {
    let len = c.seq.len();
    if !(70..=120).contains(&len) || c.amino_acid().is_none() {
        return None;
    }
    let tail = &c.seq.as_bytes()[len - 15..];
    if !tail.iter().all(|b| b"ACGT".contains(b)) {
        return None;
    }
    let mut t = tail.to_vec();
    t.extend_from_slice(b"CCA");
    Some(rc(&t))
}

impl Probes {
    pub fn from_calls(calls: &[TrnaCall]) -> Self {
        let mut full: HashMap<Vec<u8>, HashSet<String>> = HashMap::new();
        for c in calls {
            if let Some(p) = probe18(c) {
                full.entry(p).or_default().insert(c.short_name().to_string());
            }
        }
        let full: HashMap<Vec<u8>, Vec<String>> = full
            .into_iter()
            .map(|(k, v)| {
                let mut v: Vec<String> = v.into_iter().collect();
                v.sort();
                (k, v)
            })
            .collect();
        let mut short: HashMap<Vec<u8>, Vec<(u8, Vec<String>, Vec<u8>)>> = HashMap::new();
        let mut seed: HashMap<Vec<u8>, Vec<(Vec<u8>, usize)>> = HashMap::new();
        for (p, names) in &full {
            let k = &p[..K];
            short.entry(k.to_vec()).or_default().push((0, names.clone(), k.to_vec()));
            for i in 3..K {
                for &b in b"ACGT" {
                    if b != k[i] {
                        let mut v = k.to_vec();
                        v[i] = b;
                        short.entry(v).or_default().push((1, names.clone(), k.to_vec()));
                    }
                }
            }
            for o in [0, 6, 12] {
                seed.entry(p[o..o + 6].to_vec()).or_default().push((p.clone(), o));
            }
        }
        Probes { full, short, seed }
    }

    pub fn len(&self) -> usize {
        self.full.len()
    }

    /// 12-nt test on a window whose index PRE is the first base after the
    /// LTR (the window runs on the element's strand).
    fn short_hit(&self, w: &[u8]) -> Option<PbsHit> {
        let mut best: Option<PbsHit> = None;
        if w.len() < K {
            return None;
        }
        for j in 0..=w.len() - K {
            if let Some(v) = self.short.get(&w[j..j + K]) {
                for (mm, names, p12) in v {
                    let better = match &best {
                        None => true,
                        Some(b) => {
                            *mm < b.mm12
                                || (*mm == b.mm12 && j as i64 - PRE as i64 == b.offset && (names, p12) < (&b.trnas, &b.probe12))
                        }
                    };
                    if better {
                        best = Some(PbsHit { mm12: *mm, trnas: names.clone(), offset: j as i64 - PRE as i64, probe12: p12.clone(), minus: false });
                    }
                }
            }
        }
        best
    }

    /// Best 18-nt edit distance (<= 2) on a window laid out like `short_hit`'s
    /// but `MAX_EDITS` bases longer.
    fn full_edist(&self, w: &[u8]) -> Option<usize> {
        let mut best: Option<usize> = None;
        let mut tried: HashSet<(&[u8], i64)> = HashSet::new();
        if w.len() < 6 {
            return None;
        }
        for j in 0..=w.len() - 6 {
            let Some(v) = self.seed.get(&w[j..j + 6]) else { continue };
            for (p, o) in v {
                let st = j as i64 - *o as i64;
                if st < -(MAX_EDITS as i64) || st > (PRE + POST) as i64 + MAX_EDITS as i64 {
                    continue;
                }
                if !tried.insert((p.as_slice(), st)) {
                    continue;
                }
                let lo = (st - MAX_EDITS as i64).max(0) as usize;
                let hi = ((st + (PROBE + MAX_EDITS) as i64).max(0) as usize).min(w.len());
                if lo >= hi {
                    continue;
                }
                let d = semiglobal(p, &w[lo..hi]);
                if d <= MAX_EDITS && best.map_or(true, |b| d < b) {
                    best = Some(d);
                }
            }
        }
        best
    }

    /// PBS after the 5' LTR of an element on the plus strand.  `text` is the
    /// genome text, `l1e` the last base of the left LTR and `[lo, hi]` the
    /// sequence bounds (inclusive).
    fn plus(&self, text: &[u8], l1e: usize, lo: usize, hi: usize) -> (Option<PbsHit>, Option<usize>) {
        let a = (l1e + 1).saturating_sub(PRE).max(lo);
        let pad = (l1e + 1 - a) as i64 - PRE as i64; // < 0 when clipped at the sequence start
        let e = (l1e + 1 + POST + K).min(hi + 1);
        let e18 = (l1e + 1 + POST + PROBE + MAX_EDITS).min(hi + 1);
        if a >= e {
            return (None, None);
        }
        let s = self.short_hit(&text[a..e]).map(|mut h| {
            h.offset -= pad;
            h
        });
        let d = self.full_edist(&text[a..e18.max(a)]);
        (s, d)
    }

    /// PBS for an element on the minus strand: it lies just before the right
    /// LTR, read on the reverse complement.  `l2s` is the right LTR's first base.
    fn minus(&self, text: &[u8], l2s: usize, lo: usize, hi: usize) -> (Option<PbsHit>, Option<usize>) {
        let e = (l2s + PRE).min(hi + 1);
        let a = l2s.saturating_sub(POST + K).max(lo);
        let a18 = l2s.saturating_sub(POST + PROBE + MAX_EDITS).max(lo);
        if a >= e {
            return (None, None);
        }
        let w = rc(&text[a..e]);
        let pad = (e - l2s) as i64 - PRE as i64;
        let s = self.short_hit(&w).map(|mut h| {
            h.offset -= pad;
            h.minus = true;
            h
        });
        let d = self.full_edist(&rc(&text[a18..e]));
        (s, d)
    }

    /// Annotate one element.  `strand` '+' or '-' searches that orientation
    /// only; anything else searches both and keeps the better (fewer
    /// mismatches, the plus orientation on ties).
    pub fn annotate(&self, text: &[u8], l1e: usize, l2s: usize, lo: usize, hi: usize, strand: &str) -> (Option<PbsHit>, Option<usize>) {
        match strand {
            "+" => self.plus(text, l1e, lo, hi),
            "-" => self.minus(text, l2s, lo, hi),
            _ => {
                let (s1, d1) = self.plus(text, l1e, lo, hi);
                let (s2, d2) = self.minus(text, l2s, lo, hi);
                let s = match (s1, s2) {
                    (None, x) | (x, None) => x,
                    (Some(a), Some(b)) => Some(if a.mm12 <= b.mm12 { a } else { b }),
                };
                let d = match (d1, d2) {
                    (None, x) | (x, None) => x,
                    (Some(a), Some(b)) => Some(a.min(b)),
                };
                (s, d)
            }
        }
    }
}

/// One entry of the tRNA library: a distinct tRNA gene sequence.
pub struct LibEntry {
    pub id: String,
    /// ARAGORN name without "tRNA-", e.g. "Met(cat)"
    pub kind: String,
    pub seq: String,
    /// (sequence name, start, end, minus) of every copy
    pub copies: Vec<(String, i64, i64, bool)>,
    pub intron: bool,
    pub probe12: Option<Vec<u8>>,
    /// candidates / kept elements whose best PBS is this gene's probe group:
    /// [exact, one mismatch]
    pub cand: [usize; 2],
    pub kept: [usize; 2],
}

/// Collapse tRNA calls into distinct gene sequences and credit each with the
/// elements whose best 12-nt PBS hit is its probe group (genes sharing the
/// same 3'-terminal 12 nt form one group, and an element counts for each of
/// them).  `hits` gives each candidate's best hit and whether it was kept.
/// Entries are sorted by kept-element usage, then candidates, then copies.
pub fn build_library(calls: &[(&str, &TrnaCall)], hits: &[(Option<&PbsHit>, bool)]) -> Vec<LibEntry> {
    let mut usage: HashMap<&[u8], [usize; 4]> = HashMap::new();
    for (h, kept) in hits {
        if let Some(h) = h {
            let u = usage.entry(h.probe12.as_slice()).or_default();
            let m = h.mm12.min(1) as usize;
            u[m] += 1;
            if *kept {
                u[2 + m] += 1;
            }
        }
    }
    let mut by_seq: HashMap<&str, usize> = HashMap::new();
    let mut lib: Vec<LibEntry> = Vec::new();
    for (name, c) in calls {
        let i = *by_seq.entry(c.seq.as_str()).or_insert_with(|| {
            let probe12 = probe18(c).map(|p| p[..K].to_vec());
            let u = probe12.as_ref().and_then(|p| usage.get(p.as_slice())).copied().unwrap_or_default();
            lib.push(LibEntry {
                id: String::new(),
                kind: c.short_name().to_string(),
                seq: c.seq.clone(),
                copies: Vec::new(),
                intron: c.intron_len > 0,
                probe12,
                cand: [u[0], u[1]],
                kept: [u[2], u[3]],
            });
            lib.len() - 1
        });
        lib[i].copies.push((name.to_string(), c.start, c.end, c.minus));
    }
    lib.sort_by(|a, b| {
        (b.kept[0] + b.kept[1]).cmp(&(a.kept[0] + a.kept[1]))
            .then((b.cand[0] + b.cand[1]).cmp(&(a.cand[0] + a.cand[1])))
            .then(b.copies.len().cmp(&a.copies.len()))
            .then(a.kind.cmp(&b.kind))
            .then(a.seq.cmp(&b.seq))
    });
    let mut per_kind: HashMap<String, usize> = HashMap::new();
    for e in lib.iter_mut() {
        let n = per_kind.entry(e.kind.clone()).or_default();
        *n += 1;
        e.id = format!("tRNA-{}_{}", e.kind, n);
    }
    lib
}

/// Write the library as FASTA; one header per distinct gene sequence.
pub fn write_library_fasta(out: &mut impl std::io::Write, lib: &[LibEntry]) -> std::io::Result<()> {
    for e in lib {
        let (sn, s, en, m) = &e.copies[0];
        writeln!(
            out,
            ">{} type={} copies={} intron={} len={} probe12={} pbs_kept={} pbs_kept_exact={} pbs_kept_mm1={} pbs_cand={} pbs_cand_exact={} pbs_cand_mm1={} first={}:{}-{}({})",
            e.id, e.kind, e.copies.len(), if e.intron { "yes" } else { "no" }, e.seq.len(),
            e.probe12.as_ref().map_or("none".to_string(), |p| String::from_utf8_lossy(p).into_owned()),
            e.kept[0] + e.kept[1], e.kept[0], e.kept[1], e.cand[0] + e.cand[1], e.cand[0], e.cand[1],
            sn, s, en, if *m { '-' } else { '+' },
        )?;
        for chunk in e.seq.as_bytes().chunks(60) {
            out.write_all(chunk)?;
            out.write_all(b"\n")?;
        }
    }
    Ok(())
}

/// Edit distance of all of `p` against its best-matching substring of `w`.
fn semiglobal(p: &[u8], w: &[u8]) -> usize {
    let mut prev = vec![0usize; w.len() + 1];
    let mut cur = vec![0usize; w.len() + 1];
    for i in 1..=p.len() {
        cur[0] = i;
        for j in 1..=w.len() {
            let sub = prev[j - 1] + (p[i - 1] != w[j - 1]) as usize;
            cur[j] = sub.min(prev[j] + 1).min(cur[j - 1] + 1);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    *prev.iter().min().unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn call(seq: &str, name: &str) -> TrnaCall {
        TrnaCall { start: 1, end: seq.len() as i64, minus: false, name: name.into(), energy: 150.0, intron_len: 0, seq: seq.into() }
    }

    #[test]
    fn imet_pbs_found_after_ltr() {
        // last 15 nt of a tRNA-Met gene; the PBS is rc(tail + CCA)
        let tail = "CCCGGGTTCGATCCC";
        let gene = format!("{}{}", "A".repeat(60), tail);
        let probes = Probes::from_calls(&[call(&gene, "tRNA-Met(cat)")]);
        assert_eq!(probes.len(), 1);
        let mut t = tail.as_bytes().to_vec();
        t.extend_from_slice(b"CCA");
        let pbs = rc(&t);
        assert!(pbs.starts_with(b"TGG"));
        // LTR (40 x C) then 2 bp then the PBS
        let mut text = vec![b'C'; 40];
        text.extend_from_slice(b"AT");
        text.extend_from_slice(&pbs);
        text.extend_from_slice(&[b'A'; 60]);
        let (s, d) = probes.annotate(&text, 39, 100, 0, text.len() - 1, "+");
        let s = s.expect("hit");
        assert_eq!((s.mm12, s.offset, d), (0, 2, Some(0)));
        assert_eq!(s.trnas, vec!["Met(cat)".to_string()]);
        // one mismatch after the TGG is still found, a mismatch in TGG is not
        let mut t1 = text.clone();
        t1[42 + 8] = if t1[42 + 8] == b'A' { b'C' } else { b'A' };
        assert_eq!(probes.annotate(&t1, 39, 100, 0, t1.len() - 1, "+").0.unwrap().mm12, 1);
        let mut t2 = text.clone();
        t2[42 + 1] = b'A';
        assert!(probes.annotate(&t2, 39, 100, 0, t2.len() - 1, "+").0.is_none());
        // the same element on the minus strand
        let rtext = rc(&text);
        let l2s = rtext.len() - 40; // the LTR is now the last 40 bases
        let (s, _) = probes.annotate(&rtext, 0, l2s, 0, rtext.len() - 1, "-");
        let s = s.unwrap();
        assert_eq!((s.offset, s.minus), (2, true));
        // unknown strand searches both
        assert!(probes.annotate(&rtext, 0, l2s, 0, rtext.len() - 1, "?").0.is_some());
    }

    #[test]
    fn library_credits_probe_groups() {
        let met = format!("{}CCCGGGTTCGATCCC", "A".repeat(60));
        let met2 = format!("{}CCCGGGTTCGATCCC", "C".repeat(60)); // same 3' end, other gene
        let lys = format!("{}GGGGTTCGAAACCCC", "A".repeat(60));
        let calls = [call(&met, "tRNA-Met(cat)"), call(&met, "tRNA-Met(cat)"), call(&met2, "tRNA-Met(cat)"), call(&lys, "tRNA-Lys(ttt)")];
        let p_met = probe18(&calls[0]).unwrap()[..K].to_vec();
        let hit = |mm, p: &Vec<u8>| PbsHit { mm12: mm, trnas: vec![], offset: 2, probe12: p.clone(), minus: false };
        let (h0, h1) = (hit(0, &p_met), hit(1, &p_met));
        let named: Vec<(&str, &TrnaCall)> = calls.iter().map(|c| ("chr1", c)).collect();
        let lib = build_library(&named, &[(Some(&h0), true), (Some(&h1), true), (Some(&h0), false), (None, true)]);
        assert_eq!(lib.len(), 3);
        // both Met genes share the probe group and each gets the credit
        assert_eq!((lib[0].kind.as_str(), lib[0].kept, lib[0].cand), ("Met(cat)", [1, 1], [2, 1]));
        assert_eq!(lib[1].kept, [1, 1]);
        assert_eq!(lib[0].copies.len(), 2); // the duplicated gene sorts first on copies
        assert_eq!((lib[2].kind.as_str(), lib[2].kept), ("Lys(ttt)", [0, 0]));
        assert_eq!(lib[0].id, "tRNA-Met(cat)_1");
    }

    #[test]
    fn semiglobal_distance() {
        assert_eq!(semiglobal(b"ACGT", b"TTACGTTT"), 0);
        assert_eq!(semiglobal(b"ACGT", b"TTACTTT"), 1);
    }
}
