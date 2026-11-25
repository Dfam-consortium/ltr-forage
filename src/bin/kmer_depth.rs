//! Prototype: genome copy depth along library records from exact k-mer counts.
//!
//! usage: kmer_depth <genome.sufr> <library.fa> <k> > depth.tsv
//!
//! For each record position i, depth[i] is the number of genome occurrences
//! (both strands) of the k-mer centred on i, read from the sufr suffix array.
//! A k-mer containing a non-ACGT base gets depth 0.  Output: one line per
//! record, `name<TAB>d0,d1,...`.  This is the alignment-free alternative to
//! RepeatMasking the genome with the library for depth cleanup (step 6).
use anyhow::{Context, Result};
use libsufr::{OpenOptions, SourceKind, SuffixArray};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};

fn rc(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&c| match c { b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', x => x }).collect()
}

fn read_fasta(path: &str) -> Result<Vec<(String, Vec<u8>)>> {
    let f = BufReader::new(std::fs::File::open(path).with_context(|| path.to_string())?);
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    for line in f.lines() {
        let line = line?;
        if let Some(h) = line.strip_prefix('>') {
            let name = h.split_whitespace().next().unwrap_or("").split('#').next().unwrap_or("").to_string();
            out.push((name, Vec::new()));
        } else if let Some(last) = out.last_mut() {
            last.1.extend(line.trim().bytes().map(|c| c.to_ascii_uppercase()));
        }
    }
    Ok(out)
}

fn run(sufr: &str, lib: &str, k: usize) -> Result<()> {
    let recs = read_fasta(lib)?;
    // canonical k-mer -> count (filled below)
    let mut uniq: HashMap<Vec<u8>, u32> = HashMap::new();
    let ok = |w: &[u8]| w.iter().all(|c| matches!(c, b'A' | b'C' | b'G' | b'T'));
    for (_, s) in &recs {
        if s.len() < k { continue; }
        for w in s.windows(k) {
            if ok(w) {
                let r = rc(w);
                let c = if r.as_slice() < w { r } else { w.to_vec() };
                uniq.entry(c).or_insert(0);
            }
        }
    }
    eprintln!("{} records, {} distinct canonical {}-mers", recs.len(), uniq.len(), k);
    let sf = SuffixArray::open_with(sufr, &OpenOptions::default().with_text(SourceKind::Memory).with_suffixes(SourceKind::Memory))?;
    let keys: Vec<Vec<u8>> = uniq.keys().cloned().collect();
    let t0 = std::time::Instant::now();
    for chunk in keys.chunks(2_000_000) {
        // each canonical k-mer and its reverse complement (the index holds one strand)
        let mut q: Vec<Vec<u8>> = Vec::with_capacity(chunk.len() * 2);
        for c in chunk {
            q.push(c.clone());
            q.push(rc(c));
        }
        let by_q: Vec<usize> = sf.par_search_all(&q)?.iter().map(|h| h.count()).collect();
        for (i, c) in chunk.iter().enumerate() {
            let fwd = by_q[2 * i];
            let rev = if rc(c) == *c { 0 } else { by_q[2 * i + 1] };
            *uniq.get_mut(c).unwrap() = (fwd + rev).min(u32::MAX as usize) as u32;
        }
    }
    eprintln!("counted in {:.1}s", t0.elapsed().as_secs_f64());
    let mut out = BufWriter::new(std::io::stdout().lock());
    let half = k / 2;
    for (name, s) in &recs {
        let n = s.len();
        let mut d = vec![0u32; n];
        if n >= k {
            let kc: Vec<u32> = s.windows(k).map(|w| {
                if !ok(w) { return 0; }
                let r = rc(w);
                let c = if r.as_slice() < w { r } else { w.to_vec() };
                uniq[&c]
            }).collect();
            for (i, v) in d.iter_mut().enumerate() {
                // k-mer centred on i, clamped at the record ends
                let st = i.saturating_sub(half).min(kc.len() - 1);
                *v = kc[st];
            }
        }
        write!(out, "{}\t", name)?;
        for (i, v) in d.iter().enumerate() {
            if i > 0 { out.write_all(b",")?; }
            write!(out, "{}", v)?;
        }
        out.write_all(b"\n")?;
    }
    Ok(())
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    if a.len() != 4 {
        anyhow::bail!("usage: kmer_depth <genome.sufr> <library.fa> <k>");
    }
    let k: usize = a[3].parse().with_context(|| format!("k: {}", a[3]))?;
    if k == 0 {
        anyhow::bail!("k must be at least 1");
    }
    run(&a[1], &a[2], k)
}
