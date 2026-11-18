//! Internal-region statistics from k-mer counts against the genome index.
//!
//! An LTR element's internal region is repeat-derived end to end, so it
//! should not contain long stretches that occur once in the genome; the
//! region between two SINE copies that happen to sit a few kb apart is
//! ordinary sequence and usually does.  Measured on hg38, "longest single-copy
//! stretch at most 10% of the region" keeps 95% of RepeatMasker-annotated
//! elements and rejects 68% of the other candidates.  Its value is genome
//! dependent: in Arabidopsis, whose LTR families are low copy, the same test
//! rejects half the real elements, so the statistic is reported and the
//! filter is a separate, optional step.
//!
//! The region is tiled with `K`-mers every `STRIDE` bases; each tile is
//! counted on both strands with libsufr's suffix-array search, and a tile
//! seen fewer than `MULTI` times counts as single copy.

use std::collections::HashMap;

use anyhow::Result;
use libsufr::{
    sufr_file::SufrFile,
    types::{CountOptions, FromUsize, Int},
};

pub const K: usize = 20;
pub const STRIDE: usize = 10;
/// Occurrences (both strands) at or above which a tile counts as multi-copy.
pub const MULTI: usize = 3;

/// Coverage statistics for one internal region.
#[derive(Debug, Clone, Default)]
pub struct InternalStats {
    /// Tiles taken from the region (those containing N are skipped).
    pub tiles: usize,
    /// Fraction of tiles seen at least `MULTI` times.
    pub multi_frac: f64,
    /// Longest run of single-copy tiles, as a fraction of all tiles.
    pub gap_frac: f64,
    /// The same run in bases (tiles times stride).
    pub gap_bp: usize,
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

/// Tiles of one region, absolute inclusive coordinates.
pub fn tiles(text: &[u8], start: usize, end: usize) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    if end < start || end - start + 1 < K {
        return out;
    }
    let mut p = start;
    while p + K - 1 <= end {
        let t = &text[p..p + K];
        // only A/C/G/T tiles: N and the other IUPAC codes match nothing in the index
        if t.iter().all(|c| matches!(c, b'A' | b'C' | b'G' | b'T')) {
            out.push(t.to_vec());
        }
        p += STRIDE;
    }
    out
}

/// Count every distinct tile (and its reverse complement) once against the
/// index.  Returns occurrences per tile, summed over both strands.
pub fn count_tiles<T>(sf: &mut SufrFile<T>, regions: &[Vec<Vec<u8>>]) -> Result<HashMap<Vec<u8>, usize>>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    let mut index: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut queries: Vec<String> = Vec::new();
    for r in regions {
        for t in r {
            for s in [t.clone(), revcomp(t)] {
                if !index.contains_key(&s) {
                    index.insert(s.clone(), queries.len());
                    queries.push(String::from_utf8_lossy(&s).into_owned());
                }
            }
        }
    }
    let n_queries = queries.len();
    let mut raw = vec![0usize; n_queries];
    if n_queries > 0 {
        // No max_query_len: asking for one that differs from the value the
        // index was sorted with makes libsufr re-sort a copy of the suffix
        // array through a temporary file.  With none, the array is loaded as
        // built and each query is matched to its own length.
        let results = sf.count(CountOptions {
            queries,
            max_query_len: None,
            low_memory: false,
        })?;
        for r in results {
            raw[r.query_num] = r.count;
        }
    }
    Ok(both_strands(&index, &raw))
}

/// Occurrences of each queried tile on both strands: its own count plus its reverse
/// complement's.  A tile that is its own reverse complement (20 nt of AT repeats, say) has
/// one count, not two.
fn both_strands(index: &HashMap<Vec<u8>, usize>, raw: &[usize]) -> HashMap<Vec<u8>, usize> {
    let mut counts = HashMap::with_capacity(index.len());
    for (s, &i) in index {
        let both = raw[i] + index.get(&revcomp(s)).map_or(0, |&j| if j == i { 0 } else { raw[j] });
        counts.insert(s.clone(), both);
    }
    counts
}

/// Statistics for one region from the counted tiles.
pub fn stats(region: &[Vec<u8>], counts: &HashMap<Vec<u8>, usize>) -> InternalStats {
    let n = region.len();
    if n == 0 {
        return InternalStats { tiles: 0, multi_frac: 0.0, gap_frac: 1.0, gap_bp: 0 };
    }
    let multi: Vec<bool> = region.iter().map(|t| counts.get(t).copied().unwrap_or(0) >= MULTI).collect();
    let n_multi = multi.iter().filter(|&&m| m).count();
    let (mut run, mut best) = (0usize, 0usize);
    for &m in &multi {
        run = if m { 0 } else { run + 1 };
        best = best.max(run);
    }
    InternalStats {
        tiles: n,
        multi_frac: n_multi as f64 / n as f64,
        gap_frac: best as f64 / n as f64,
        gap_bp: best * STRIDE,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gap_statistics() {
        let region: Vec<Vec<u8>> = (0..10u8).map(|i| vec![i; K]).collect();
        let mut counts = HashMap::new();
        for (i, t) in region.iter().enumerate() {
            // tiles 3, 4, 5 are single copy; the rest multi-copy
            counts.insert(t.clone(), if (3..=5).contains(&i) { 1 } else { 5 });
        }
        let s = stats(&region, &counts);
        assert_eq!(s.tiles, 10);
        assert!((s.multi_frac - 0.7).abs() < 1e-9);
        assert!((s.gap_frac - 0.3).abs() < 1e-9);
        assert_eq!(s.gap_bp, 3 * STRIDE);
    }

    #[test]
    fn palindromic_tile_counted_once() {
        let pal = b"ATATATATATATATATATAT".to_vec();
        assert_eq!(revcomp(&pal), pal);
        let (fwd, rev) = (b"AAAAAAAAAAAAAAAAAAAC".to_vec(), b"GTTTTTTTTTTTTTTTTTTT".to_vec());
        let index: HashMap<Vec<u8>, usize> = [(pal.clone(), 0), (fwd.clone(), 1), (rev.clone(), 2)].into_iter().collect();
        let c = both_strands(&index, &[2, 1, 4]);
        assert_eq!(c[&pal], 2); // at two loci: below MULTI, where the doubled count was 4
        assert_eq!((c[&fwd], c[&rev]), (5, 5));
    }

    #[test]
    fn tiling_skips_n() {
        let mut text = vec![b'A'; 100];
        text[45] = b'N';
        let t = tiles(&text, 0, 99);
        // starts 0,10,...,80: those covering position 45 (starts 30 and 40) are skipped
        assert_eq!(t.len(), 7);
    }
}
