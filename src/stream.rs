/// Sequential, buffered reading of the suffix array and LCP array of a
/// `.sufr` file through libsufr's range reads.
///
/// The bottom-up traversal in `maxpairs.rs` consumes both arrays once, in
/// order, so there is no reason to hold them in memory: on hg38 they are
/// 24 GB each.  The index is opened with the arrays on disk, and
/// `ArrayReader` reads 262,144 values at a time (positioned reads, so the
/// workers of the parallel traversal share one opened index).

use anyhow::{Context, Result};
use libsufr::SuffixArray;

/// Values per read: a 1-2 MB read from the file into a 2 MB buffer, per reader (two per worker).
const CHUNK_ELEMS: usize = 1 << 18;

#[derive(Clone, Copy, Debug)]
pub enum Array {
    Suffixes,
    Lcps,
}

fn read(sa: &SuffixArray, which: Array, ranks: std::ops::Range<usize>) -> Result<Vec<usize>> {
    match which {
        Array::Suffixes => sa.suffixes(ranks.clone()),
        Array::Lcps => sa.lcps(ranks.clone()),
    }
    .with_context(|| format!("reading the index's {:?} at ranks {:?}", which, ranks))
}

pub struct ArrayReader<'a> {
    sa: &'a SuffixArray,
    which: Array,
    /// The next rank to read from the file, and the end of the range.
    next: usize,
    end: usize,
    buf: Vec<usize>,
    pos: usize,
}

impl<'a> ArrayReader<'a> {
    /// Yield the values of `which` at ranks `start..end`.
    pub fn new(sa: &'a SuffixArray, which: Array, start: usize, end: usize) -> Self {
        ArrayReader { sa, which, next: start, end, buf: Vec::new(), pos: 0 }
    }

    fn refill(&mut self) -> bool {
        if self.next >= self.end {
            return false;
        }
        let to = (self.next + CHUNK_ELEMS).min(self.end);
        // libsufr checked the file's length against its header when it opened it, so a failed read
        // here is an I/O error, which the traversal has no way to recover from
        self.buf = read(self.sa, self.which, self.next..to).unwrap_or_else(|e| panic!("{:#}", e));
        self.next = to;
        self.pos = 0;
        true
    }
}

impl Iterator for ArrayReader<'_> {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<usize> {
        if self.pos == self.buf.len() && !self.refill() {
            return None;
        }
        let v = self.buf[self.pos];
        self.pos += 1;
        Some(v)
    }
}

// --------------------------------------------------
/// Split the suffix array into ranges that can be traversed independently.
///
/// An lcp-interval of depth `d` covers SA positions `lb..=rb` with
/// `lcp[i] >= d` for every `i` in `lb+1..=rb`.  So any position `i` with
/// `lcp[i] < min_len` is never strictly inside an interval of depth
/// `min_len` or more, and cutting the array there keeps every interval the
/// maximal-repeat traversal cares about within one piece.  Cutting at such
/// positions and concatenating each piece's seeds in order therefore yields
/// exactly the sequential result.
///
/// `n_pieces` is a target: each cut is placed at the first qualifying
/// position at or after `k * num_suffixes / n_pieces`, so a long run of
/// `lcp >= min_len` (a large repeat family) stays in one piece and the
/// actual count may be smaller.  Returns half-open `(start, end)` ranges.
pub fn split_points(sa: &SuffixArray, num_suffixes: usize, min_len: usize, n_pieces: usize) -> Result<Vec<(usize, usize)>> {
    // the qualifying position is usually within a few values of the target: read small blocks
    const BLOCK: usize = 4096;
    let mut cuts = vec![0usize];
    for k in 1..n_pieces.max(1) {
        let target = k * num_suffixes / n_pieces;
        let last = *cuts.last().unwrap();
        if target <= last {
            continue;
        }
        let mut i = target;
        'scan: while i < num_suffixes {
            let to = (i + BLOCK).min(num_suffixes);
            for v in read(sa, Array::Lcps, i..to)? {
                if v < min_len {
                    break 'scan;
                }
                i += 1;
            }
        }
        if i < num_suffixes && i > last {
            cuts.push(i);
        }
    }
    cuts.push(num_suffixes);
    Ok(cuts.windows(2).map(|w| (w[0], w[1])).filter(|(s, e)| e > s).collect())
}
