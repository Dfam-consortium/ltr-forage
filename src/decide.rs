//! The keep-or-reject decision for each candidate, and the library FASTA.
//!
//! One rule with a strictness knob.  Every level rejects a candidate whose
//! protein hits are of not-LTR or mixed order; the levels differ in what
//! counts as positive evidence.  Measured against RepeatMasker (section 12
//! of the analysis document): on hg38 an LTR-order domain hit alone is 89%
//! precise at 68% recall; the TG..CA motif with a TSD of 4 bp or more adds
//! a few elements on every genome (human ERVs mostly make 4 bp TSDs, plant
//! Gypsy and Copia 5 bp); a TSD alone lets in many chance pairs.

use std::fmt::Write as _;
use std::io::Write;

use anyhow::Result;

use crate::domain::DomainCall;

#[derive(Debug, Clone, Copy, clap::ValueEnum, PartialEq)]
pub enum KeepLevel {
    /// An LTR-order domain hit
    Domain,
    /// A domain hit, or the TG..CA motif with a TSD of 4 bp or more
    Standard,
    /// A domain hit, a TSD of 5 bp or more, or TG..CA with a TSD of 4 bp or more
    Loose,
    /// Everything the veto does not reject
    None,
}

pub struct Evidence<'a> {
    pub domain: Option<&'a DomainCall>,
    pub tsd_len: usize,
    /// 1 when the TSD is the one-mismatch fallback (`features::tsd_or_mm1`).
    /// It counts only alongside TG..CA: elsewhere it is no commoner than
    /// in random flanks.
    pub tsd_mm: usize,
    pub motif: &'a str,
    /// Longest single-copy stretch of the internal region as a fraction of
    /// its tiles, when `--internal-stats` ran.
    pub gap_frac: Option<f64>,
    /// The larger of the two LTRs' self-repeat coverage (see
    /// `reextend::self_repeat_cov`).  The threshold has to sit near 1: real
    /// LTRs with tandem sub-repeats (LTR12 in human) reach 0.8 to 0.98,
    /// while satellite monomers and Alu dimers reach 1.0.
    pub self_cov: f64,
    /// LTR-order domain hits found inside either LTR.  An LTR does not
    /// encode the polyprotein; a hit there means the "LTR" is internal
    /// sequence (two adjacent copies of one family called as one element,
    /// or a boundary that landed in the coding region).
    pub ltr_domain_hits: usize,
    /// Both LTRs end in a poly-A tail (or both start with poly-T, the
    /// other strand); see `features::polya_ltrs`.
    pub polya: bool,
}

pub struct Decision {
    pub keep: bool,
    /// Short account of what decided it, for the table and the FASTA header.
    pub why: String,
}

pub fn decide(level: KeepLevel, gap_max: Option<f64>, tandem_max: f64, ev: &Evidence) -> Decision {
    let mut why = String::new();
    let order = ev.domain.map(|d| d.order.as_str()).unwrap_or(".");
    if order == "notLTR" || order == "mixture" {
        write!(why, "veto:{}", order).unwrap();
        return Decision { keep: false, why };
    }
    if ev.self_cov >= tandem_max {
        write!(why, "veto:tandem{:.2}", ev.self_cov).unwrap();
        return Decision { keep: false, why };
    }
    if ev.ltr_domain_hits > 0 {
        write!(why, "veto:ltrdomain{}", ev.ltr_domain_hits).unwrap();
        return Decision { keep: false, why };
    }
    let dom = order == "LTR";
    let tgca = ev.motif == "TGCA";
    // two LTRs that end in a poly-A tail without the TG..CA termini are the 3' ends of two copies of a
    // non-LTR retroposon flanking some other sequence (mm10: 22 kept pairs of B1/B2 SINE or L1 3' ends,
    // 115-185 and 700-1,150 bp; an L1 3' end called as an LTR masked 90 Mb of LINE sequence)
    if ev.polya && !tgca {
        why.push_str("veto:polyA");
        return Decision { keep: false, why };
    }
    if dom {
        write!(why, "domain").unwrap();
    }
    if ev.tsd_len > 0 {
        write!(why, "{}tsd{}{}", if why.is_empty() { "" } else { "," }, ev.tsd_len, if ev.tsd_mm > 0 { "mm1" } else { "" }).unwrap();
    }
    if tgca {
        write!(why, "{}tgca", if why.is_empty() { "" } else { "," }).unwrap();
    }
    let positive = match level {
        KeepLevel::Domain => dom,
        KeepLevel::Standard => dom || (tgca && ev.tsd_len >= 4),
        KeepLevel::Loose => dom || (ev.tsd_len >= 5 && ev.tsd_mm == 0) || (tgca && ev.tsd_len >= 4),
        KeepLevel::None => true,
    };
    if !positive {
        if why.is_empty() {
            why.push_str("none");
        }
        return Decision { keep: false, why };
    }
    if let (Some(max), Some(gap)) = (gap_max, ev.gap_frac) {
        if gap > max {
            write!(why, ",gap{:.2}", gap).unwrap();
            return Decision { keep: false, why };
        }
    }
    Decision { keep: true, why }
}

/// One kept element, ready to be written: left LTR, internal region and
/// right LTR as three FASTA records named the way LTR_retriever names its
/// redundant library (`seq:start..end_LTR#LTR/Gypsy`), coordinates 1-based
/// and inclusive, sequences in genome orientation.
pub fn write_element<W: Write>(
    out: &mut W,
    seq_name: &str,
    text: &[u8],
    seq_off: usize,
    ltrs: (usize, usize, usize, usize),
    superfamily: &str,
    why: &str,
) -> Result<()> {
    let (l1s, l1e, l2s, l2e) = ltrs;
    let rel1 = |abs: usize| abs - seq_off + 1;
    let mut rec = |s: usize, e: usize, kind: &str| -> Result<()> {
        writeln!(out, ">{}:{}..{}_{}#LTR/{} {}", seq_name, rel1(s), rel1(e), kind, superfamily, why)?;
        out.write_all(&text[s..=e])?;
        writeln!(out)?;
        Ok(())
    };
    rec(l1s, l1e, "LTR")?;
    if l2s > l1e + 1 {
        rec(l1e + 1, l2s - 1, "INT")?;
    }
    rec(l2s, l2e, "LTR")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn call(order: &str) -> DomainCall {
        DomainCall { order: order.into(), ..Default::default() }
    }

    #[test]
    fn levels() {
        let ltr = call("LTR");
        let not = call("notLTR");
        fn ev<'a>(d: Option<&'a DomainCall>, tsd: usize, motif: &'a str) -> Evidence<'a> {
            Evidence { domain: d, tsd_len: tsd, tsd_mm: 0, motif, gap_frac: None, self_cov: 0.0, ltr_domain_hits: 0, polya: false }
        }
        assert!(decide(KeepLevel::Domain, None, 0.5, &ev(Some(&ltr), 0, "AAAA")).keep);
        assert!(!decide(KeepLevel::None, None, 0.5, &ev(Some(&not), 6, "TGCA")).keep);
        assert!(!decide(KeepLevel::Domain, None, 0.5, &ev(None, 6, "TGCA")).keep);
        assert!(decide(KeepLevel::Standard, None, 0.5, &ev(None, 5, "TGCA")).keep);
        assert!(decide(KeepLevel::Standard, None, 0.5, &ev(None, 4, "TGCA")).keep);
        assert!(!decide(KeepLevel::Standard, None, 0.5, &ev(None, 3, "TGCA")).keep);
        assert!(!decide(KeepLevel::Standard, None, 0.5, &ev(None, 5, "TGAA")).keep);
        assert!(decide(KeepLevel::Loose, None, 0.5, &ev(None, 5, "TGAA")).keep);
        assert!(decide(KeepLevel::Loose, None, 0.5, &ev(None, 4, "TGCA")).keep);
        assert!(!decide(KeepLevel::Loose, None, 0.5, &ev(None, 4, "TGAA")).keep);
        let d = decide(KeepLevel::Standard, None, 0.5, &ev(Some(&ltr), 5, "TGCA"));
        assert!(d.keep);
        assert_eq!(d.why, "domain,tsd5,tgca");
        assert_eq!(decide(KeepLevel::Standard, None, 0.5, &ev(Some(&not), 0, "TGCA")).why, "veto:notLTR");
    }

    #[test]
    fn tandem_veto() {
        let ltr = call("LTR");
        let ev = Evidence { domain: Some(&ltr), tsd_len: 6, tsd_mm: 0, motif: "TGCA", gap_frac: None, self_cov: 0.8, ltr_domain_hits: 0, polya: false };
        let d = decide(KeepLevel::None, None, 0.5, &ev);
        assert!(!d.keep);
        assert_eq!(d.why, "veto:tandem0.80");
        assert!(decide(KeepLevel::None, None, 0.9, &ev).keep);
    }

    #[test]
    fn polya_veto() {
        let ltr = call("LTR");
        let mut ev = Evidence { domain: Some(&ltr), tsd_len: 0, tsd_mm: 0, motif: "ACAA", gap_frac: None, self_cov: 0.0, ltr_domain_hits: 0, polya: true };
        assert_eq!(decide(KeepLevel::Standard, None, 0.99, &ev).why, "veto:polyA");
        ev.motif = "TGCA";
        assert!(decide(KeepLevel::Standard, None, 0.99, &ev).keep);
    }

    #[test]
    fn ltr_domain_veto() {
        let ltr = call("LTR");
        let ev = Evidence { domain: Some(&ltr), tsd_len: 6, tsd_mm: 0, motif: "TGCA", gap_frac: None, self_cov: 0.0, ltr_domain_hits: 2, polya: false };
        let d = decide(KeepLevel::None, None, 0.99, &ev);
        assert!(!d.keep);
        assert_eq!(d.why, "veto:ltrdomain2");
    }

    #[test]
    fn gap_filter_applies_only_when_measured() {
        let ltr = call("LTR");
        let mut ev = Evidence { domain: Some(&ltr), tsd_len: 0, tsd_mm: 0, motif: "TGCA", gap_frac: Some(0.3), self_cov: 0.0, ltr_domain_hits: 0, polya: false };
        assert!(!decide(KeepLevel::Domain, Some(0.1), 0.5, &ev).keep);
        assert!(decide(KeepLevel::Domain, None, 0.5, &ev).keep);
        ev.gap_frac = None;
        assert!(decide(KeepLevel::Domain, Some(0.1), 0.5, &ev).keep);
    }

    #[test]
    fn fasta_records() {
        let text = b"NNAAAACCCCGGGGTTTTNN";
        let mut out = Vec::new();
        write_element(&mut out, "s", text, 2, (2, 5, 14, 17), "unknown", "domain").unwrap();
        let s = String::from_utf8(out).unwrap();
        assert_eq!(s, ">s:1..4_LTR#LTR/unknown domain\nAAAA\n>s:5..12_INT#LTR/unknown domain\nCCCCGGGG\n>s:13..16_LTR#LTR/unknown domain\nTTTT\n");
    }
}
