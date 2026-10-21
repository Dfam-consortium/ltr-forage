//! Tandem-array mask for the seed search (`--mask-tandem`).
//!
//! A satellite or simple-repeat array of unit p yields a seed pair between
//! almost every two units within `--maxdistltr`, so a telomere-to-telomere
//! assembly, with its centromeric and subtelomeric arrays assembled, floods
//! phases 1 and 2: gorilla (mGorGor1, 311 Mb of satellite) gave 953 M seed
//! pairs against hg38's 11 M at the same settings, and 4.5 h of the run.
//! HiRepEat collapses such arrays with FasTAN before its LTR search for the
//! same reason; here they are masked instead, so coordinates do not move.
//!
//! A position is periodic when its k-mer last occurred at most `max_period`
//! bases earlier.  An LTR element's two LTRs are at least `--mindistltr`
//! apart, so with `max_period` below that no element repeat makes a position
//! periodic, however the element is arranged; only shorter-period sequence
//! does.  Runs of at least `min_len` bases in which a window-smoothed share
//! of periodic positions reaches `min_frac` are masked: the seed search
//! treats them as wildcards (see `maxpairs::WildRuns`), and the text itself is
//! unchanged, so extension, TSDs and output still see the real sequence.

use rayon::prelude::*;

const K: usize = 10;
const MASK: u32 = (1 << (2 * K)) - 1;
const CHUNK: usize = 8 << 20;
const WINDOW: usize = 64;

#[derive(Clone, Copy, Debug)]
pub struct TandemParams {
    pub max_period: usize,
    pub min_len: usize,
    pub min_frac: f64,
}

#[inline]
fn code(c: u8) -> Option<u32> {
    match c {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Periodic flags for text[from..to], using text[from - warm..from] only to
/// fill the k-mer table.
fn flags(text: &[u8], from: usize, to: usize, warm: usize, max_period: usize, last: &mut [u32]) -> Vec<bool> {
    last.iter_mut().for_each(|x| *x = u32::MAX);
    let start = from.saturating_sub(warm);
    let mut out = vec![false; to - from];
    let (mut kmer, mut valid) = (0u32, 0usize);
    for i in start..to {
        match code(text[i]) {
            Some(c) => {
                kmer = ((kmer << 2) | c) & MASK;
                valid += 1;
            }
            None => {
                valid = 0;
                continue;
            }
        }
        if valid < K {
            continue;
        }
        let pos = (i - start) as u32; // k-mer ends at i
        let prev = last[kmer as usize];
        if i >= from && prev != u32::MAX && (pos - prev) as usize <= max_period {
            // the whole k-mer lies inside a repeat of period <= max_period
            for j in (i + 1 - K).max(from)..=i {
                out[j - from] = true;
            }
        }
        last[kmer as usize] = pos;
    }
    out
}

/// Masked runs `[start, end)`, sorted and disjoint.
pub fn tandem_runs(text: &[u8], p: TandemParams) -> Vec<(usize, usize)> {
    if p.max_period == 0 {
        return Vec::new();
    }
    let chunks: Vec<(usize, usize)> = (0..text.len()).step_by(CHUNK).map(|s| (s, (s + CHUNK).min(text.len()))).collect();
    let per: Vec<Vec<(usize, usize)>> = chunks
        .par_iter()
        .map_init(
            || vec![u32::MAX; 1 << (2 * K)],
            |last, &(from, to)| {
                // reach back so runs crossing a chunk boundary are judged on both sides
                let pad = p.min_len.max(WINDOW);
                let (f0, t0) = (from.saturating_sub(pad), (to + pad).min(text.len()));
                let fl = flags(text, f0, t0, p.max_period, p.max_period, last);
                // window share of periodic positions, via prefix sums
                let mut pre = vec![0u32; fl.len() + 1];
                for (i, &b) in fl.iter().enumerate() {
                    pre[i + 1] = pre[i] + b as u32;
                }
                let half = WINDOW / 2;
                let need = (p.min_frac * WINDOW as f64).ceil() as u32;
                let mut runs = Vec::new();
                let mut cur: Option<usize> = None;
                for i in 0..=fl.len() {
                    let hit = i < fl.len() && {
                        let (a, b) = (i.saturating_sub(half), (i + half).min(fl.len()));
                        pre[b] - pre[a] >= need.min((b - a) as u32)
                    };
                    match (hit, cur) {
                        (true, None) => cur = Some(i),
                        (false, Some(s)) => {
                            if i - s >= p.min_len {
                                runs.push((f0 + s, f0 + i));
                            }
                            cur = None;
                        }
                        _ => {}
                    }
                }
                // keep the parts inside this chunk; neighbours report theirs
                runs.into_iter()
                    .filter_map(|(a, b)| {
                        let (a, b) = (a.max(from), b.min(to));
                        (a < b).then_some((a, b))
                    })
                    .collect()
            },
        )
        .collect();
    let mut all: Vec<(usize, usize)> = per.into_iter().flatten().collect();
    all.sort_unstable();
    let mut out: Vec<(usize, usize)> = Vec::new();
    for (a, b) in all {
        match out.last_mut() {
            Some(l) if a <= l.1 => l.1 = l.1.max(b),
            _ => out.push((a, b)),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rand_seq(n: usize, mut x: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x % 4) as usize]
            })
            .collect()
    }

    #[test]
    fn masks_satellite_not_ltr_pair() {
        let p = TandemParams { max_period: 999, min_len: 1000, min_frac: 0.5 };
        // a 171-bp unit repeated 30 times between random flanks
        let unit = rand_seq(171, 7);
        let mut t = rand_seq(20_000, 1);
        let sat: Vec<u8> = unit.iter().cycle().take(171 * 30).cloned().collect();
        t.extend(&sat);
        t.extend(rand_seq(20_000, 2));
        let runs = tandem_runs(&t, p);
        let masked: usize = runs.iter().map(|(a, b)| b - a).sum();
        assert!(masked >= 171 * 28 && masked <= 171 * 30 + 200, "{runs:?}");
        // an "element": 500-bp LTRs 3 kb apart is not masked
        let ltr = rand_seq(500, 9);
        let mut e = rand_seq(10_000, 3);
        e.extend(&ltr);
        e.extend(rand_seq(2_500, 4));
        e.extend(&ltr);
        e.extend(rand_seq(10_000, 5));
        assert!(tandem_runs(&e, p).is_empty());
    }
}
