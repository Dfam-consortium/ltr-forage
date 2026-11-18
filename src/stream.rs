/// Sequential, buffered reading of the suffix array and LCP array straight
/// from a `.sufr` file.
///
/// The bottom-up traversal in `maxpairs.rs` consumes both arrays once, in
/// order, so there is no reason to hold them in memory: on hg38 that was two
/// `Vec<usize>` of 24 GB each.  `ArrayReader` reads a few megabytes at a
/// time and yields each element as a `usize`.

use anyhow::{bail, Result};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};

/// Elements per read.  At 4 or 8 bytes each this is a 4-8 MB buffer.
const CHUNK_ELEMS: usize = 1 << 20;

pub struct ArrayReader {
    file: File,
    /// Bytes per element: 4 for a u32 index, 8 for u64.
    width: usize,
    remaining: usize,
    buf: Vec<u8>,
    pos: usize,
    filled: usize,
}

impl ArrayReader {
    /// Open `path` and position at `byte_offset`, ready to yield
    /// `num_elements` integers of `width` bytes each (native byte order,
    /// as written by sufr).
    pub fn open(path: &str, byte_offset: u64, num_elements: usize, width: usize) -> Result<Self> {
        if width != 4 && width != 8 {
            bail!("unsupported element width {width}");
        }
        let mut file = File::open(path)?;
        // a truncated index is an error here rather than a panic in a worker's read
        let (have, need) = (file.metadata()?.len(), byte_offset + (num_elements * width) as u64);
        if have < need {
            bail!("{path} is truncated: {have} bytes, but its header places array data up to byte {need}");
        }
        file.seek(SeekFrom::Start(byte_offset))?;
        Ok(ArrayReader {
            file,
            width,
            remaining: num_elements,
            buf: vec![0; CHUNK_ELEMS * width],
            pos: 0,
            filled: 0,
        })
    }

    fn refill(&mut self) -> bool {
        if self.remaining == 0 {
            return false;
        }
        let n = self.remaining.min(CHUNK_ELEMS);
        let bytes = n * self.width;
        self.file
            .read_exact(&mut self.buf[..bytes])
            .expect("read error in the .sufr file (its length was checked when it was opened)");
        self.remaining -= n;
        self.pos = 0;
        self.filled = bytes;
        true
    }
}

impl Iterator for ArrayReader {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<usize> {
        if self.pos == self.filled && !self.refill() {
            return None;
        }
        let b = &self.buf[self.pos..self.pos + self.width];
        self.pos += self.width;
        Some(if self.width == 4 {
            u32::from_ne_bytes([b[0], b[1], b[2], b[3]]) as usize
        } else {
            u64::from_ne_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]) as usize
        })
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
pub fn split_points(
    path: &str,
    lcp_pos: u64,
    num_suffixes: usize,
    width: usize,
    min_len: usize,
    n_pieces: usize,
) -> Result<Vec<(usize, usize)>> {
    let mut cuts = vec![0usize];
    for k in 1..n_pieces.max(1) {
        let target = k * num_suffixes / n_pieces;
        let last = *cuts.last().unwrap();
        if target <= last {
            continue;
        }
        let mut rd = ArrayReader::open(
            path,
            lcp_pos + (target * width) as u64,
            num_suffixes - target,
            width,
        )?;
        let mut i = target;
        while let Some(v) = rd.next() {
            if v < min_len {
                break;
            }
            i += 1;
        }
        if i < num_suffixes && i > last {
            cuts.push(i);
        }
    }
    cuts.push(num_suffixes);
    Ok(cuts.windows(2).map(|w| (w[0], w[1])).filter(|(s, e)| e > s).collect())
}
