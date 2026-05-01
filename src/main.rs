use anyhow::Result;
use clap::Parser;
use libsufr::{
    sufr_file::SufrFile,
    types::{FromUsize, Int},
    util::read_text_length,
};
use rayon::prelude::*;

mod ltr;
mod maxpairs;
mod tsd;
mod xdrop;

// --------------------------------------------------
#[derive(Debug, Clone, clap::ValueEnum, PartialEq)]
enum OutputFormat {
    /// GenomeTools LTR_harvest compatible output (default)
    Gt,
}

// --------------------------------------------------
#[derive(Parser, Debug)]
#[command(
    name = "ltr-sufrage",
    about = "Identify LTR retrotransposon candidates using a sufr suffix array index",
    allow_negative_numbers = true,
)]
struct Args {
    /// Path to a .sufr index file (created with `sufr create --dna`)
    sufr: String,

    /// Minimum seed length (minimum exact match length between LTR pair)
    #[arg(short = 'm', long, default_value_t = 30, value_name = "BP")]
    min_seed: usize,

    /// Minimum LTR length
    #[arg(short = 'l', long, default_value_t = 100, value_name = "BP")]
    ltr_min: usize,

    /// Maximum LTR length
    #[arg(short = 'L', long, default_value_t = 1000, value_name = "BP")]
    ltr_max: usize,

    /// Minimum distance between the 5′ ends of the two LTRs
    #[arg(short = 'd', long, default_value_t = 1000, value_name = "BP")]
    dist_min: usize,

    /// Maximum distance between the 5′ ends of the two LTRs
    #[arg(short = 'D', long, default_value_t = 15000, value_name = "BP")]
    dist_max: usize,

    /// Minimum percent similarity between LTR copies [0..100]
    #[arg(short = 's', long, default_value_t = 85.0, value_name = "PCT")]
    similar: f64,

    /// X-drop score threshold for extension alignment
    #[arg(short = 'x', long, default_value_t = 5, value_name = "SCORE")]
    xdrop: i32,

    /// Match score for extension alignment
    #[arg(long, default_value_t = 2)]
    mat: i32,

    /// Mismatch score for extension alignment (negative)
    #[arg(long, default_value_t = -2)]
    mis: i32,

    /// Insertion score for extension alignment (negative)
    #[arg(long, default_value_t = -3)]
    ins: i32,

    /// Deletion score for extension alignment (negative)
    #[arg(long, default_value_t = -3)]
    del: i32,

    /// Minimum TSD length (0 = disable TSD search)
    #[arg(long, default_value_t = 4, value_name = "BP")]
    mintsd: usize,

    /// Maximum TSD length
    #[arg(long, default_value_t = 20, value_name = "BP")]
    maxtsd: usize,

    /// Search window around LTR boundaries for TSD/motif (bp)
    #[arg(long, default_value_t = 60, value_name = "BP")]
    vicinity: usize,

    /// Enable TG...CA palindromic motif search
    #[arg(long, default_value_t = false)]
    motif: bool,

    /// Allowed mismatches in the motif (0–3)
    #[arg(long, default_value_t = 0, value_name = "N")]
    motifmis: u32,

    /// Overlap removal strategy: all (keep all), best (keep best per cluster), no (remove all overlapping)
    #[arg(long, default_value = "best", value_name = "STRATEGY")]
    overlaps: String,

    /// Number of threads for parallel seed extension (0 = use all available CPUs)
    #[arg(short = 'j', long, default_value_t = 0, value_name = "N")]
    threads: usize,

    /// Use tiered TSD selection (LTR_retriever style) instead of GT-compatible deviation-only
    #[arg(long, default_value_t = false)]
    tiered_tsd: bool,

    /// Output format
    #[arg(long, default_value = "gt", value_name = "FORMAT")]
    format: OutputFormat,
}

// --------------------------------------------------
/// Remove overlapping candidates, mirroring GT's gt_removeoverlapswithlowersimilarity.
///
/// Requires input to be sorted by (seq_num, left_ltr_start, right_ltr_end).
///
/// With `no_overlaps=false` ("best" mode): within each cluster of mutually
/// overlapping elements, keep the one with the highest similarity.
/// With `no_overlaps=true` ("no" mode): discard all elements in any cluster
/// that has more than one element (no overlaps at all).
///
/// Returns the indices of surviving candidates (into the sorted input slice).
fn overlapping_keep_indices(
    candidates: &[ltr::LtrCandidate],
    no_overlaps: bool,
) -> Vec<usize> {
    if candidates.is_empty() {
        return vec![];
    }

    let n = candidates.len();
    let mut skipped = vec![false; n];

    // Track the current overlap cluster.
    let mut ref_seq   = candidates[0].seq_num;
    let mut ref_start = candidates[0].left_ltr_start;
    let mut ref_end   = candidates[0].right_ltr_end;
    let mut max_sim_i = 0usize; // index of highest-similarity element in cluster

    for i in 1..n {
        if skipped[i] { continue; }
        let c = &candidates[i];
        let overlaps = c.seq_num == ref_seq
            && ref_start <= c.right_ltr_end
            && ref_end   >= c.left_ltr_start;

        if overlaps {
            // Extend cluster right boundary.
            ref_end = ref_end.max(c.right_ltr_end);
            if no_overlaps {
                // "no" mode: mark both the cluster best and this one.
                skipped[max_sim_i] = true;
                skipped[i] = true;
            } else {
                // "best" mode: keep higher similarity.
                if candidates[max_sim_i].similarity < c.similarity {
                    skipped[max_sim_i] = true;
                    max_sim_i = i;
                } else {
                    skipped[i] = true;
                }
            }
        } else {
            // New cluster.
            ref_seq   = c.seq_num;
            ref_start = c.left_ltr_start;
            ref_end   = c.right_ltr_end;
            max_sim_i = i;
        }
    }

    (0..n).filter(|&i| !skipped[i]).collect()
}

// --------------------------------------------------
fn run<T>(args: &Args) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok(); // ok() because build_global fails if already initialised
    }
    // threads == 0: leave Rayon at its default (all logical CPUs, or RAYON_NUM_THREADS)

    let mut sf: SufrFile<T> = SufrFile::read(&args.sufr, false)?;

    let sa: Vec<usize> = sf.suffix_array_file.iter().map(|v| v.to_usize()).collect();
    let lcp: Vec<usize> = sf.lcp_file.iter().map(|v| v.to_usize()).collect();
    let text = sf.text.clone();
    let num_suffixes = sf.len_suffixes.to_usize();
    let seq_starts: Vec<usize> = sf.sequence_starts.iter().map(|v| v.to_usize()).collect();
    let seq_names = sf.sequence_names.clone();

    eprintln!(
        "Index: {} sequence(s), {} suffixes, text length {}",
        sf.num_sequences.to_usize(),
        num_suffixes,
        sf.text_len.to_usize(),
    );
    eprintln!(
        "Parameters: min_seed={} ltr=[{},{}] dist=[{},{}] similar={:.0}% xdrop={}",
        args.min_seed,
        args.ltr_min,
        args.ltr_max,
        args.dist_min,
        args.dist_max,
        args.similar,
        args.xdrop,
    );

    let scores = xdrop::XdropScores {
        mat: args.mat,
        mis: args.mis,
        ins: args.ins,
        del: args.del,
        xdrop: args.xdrop,
    };

    // ---- Phase 1: collect seeds ----
    let mut seeds: Vec<ltr::Seed> = Vec::new();
    maxpairs::enumerate_maximal_repeats(
        &sa,
        &lcp,
        &text,
        num_suffixes,
        args.min_seed,
        &mut |pos1, pos2, len| {
            if let Some(seed) = ltr::filter_seed(
                pos1, pos2, len,
                args.ltr_min, args.ltr_max,
                args.dist_min, args.dist_max,
                &seq_starts,
            ) {
                seeds.push(seed);
            }
        },
    );
    eprintln!("Phase 1: {} seed pair(s) found", seeds.len());
    // ---- Phase 2 + TSD: extend seeds, annotate TSD, compute post-TSD similarity ----
    // GT runs TSD annotation per-seed inside gt_searchforLTRs, before the results are
    // collected into the arrayLTRboundaries.  The sort and dedup that follow therefore
    // operate on POST-TSD outer boundaries.  We mirror that by annotating TSD
    // immediately after extension, before collecting into the candidate pool.
    let motif_cfg = tsd::MotifConfig {
        firstleft:          b'T',
        secondleft:         b'G',
        firstright:         b'C',
        secondright:        b'A',
        allowed_mismatches: args.motifmis,
    };

    // Process seeds in parallel: extension + TSD annotation + post-TSD similarity.
    // Each seed is independent (all inputs are read-only), so this is embarrassingly
    // parallel — mirroring GT's gt_multithread approach.
    //
    // To guarantee deterministic output regardless of thread completion order, we
    // attach the seed's SA-enumeration index and use it as a tiebreaker in the sort.
    // Seeds that produce candidates with the same post-TSD outer boundaries always
    // yield the same survivor (the one from the lowest-index seed), matching what
    // a single-threaded stable sort on insertion order would produce.
    let mut paired: Vec<(usize, ltr::LtrCandidate, tsd::TsdMotifInfo)> = seeds
        .par_iter()
        .enumerate()
        .filter_map(|(seed_idx, seed)| {
            let mut cand = ltr::extend_seed(
                seed, &text, &seq_starts,
                args.ltr_min, args.ltr_max,
                args.dist_min, args.dist_max,
                &scores,
            )?;

            let ti = tsd::annotate_tsd_and_motif(
                &cand, &text, &seq_starts,
                args.mintsd, args.maxtsd, args.vicinity,
                &motif_cfg, args.motif, args.tiered_tsd,
            );

            // Compute post-TSD similarity and apply threshold (GT does this per-seed).
            let llen = ti.left_ltr_end  - ti.left_ltr_start  + 1;
            let rlen = ti.right_ltr_end - ti.right_ltr_start + 1;
            let edist = ltr::unit_edit_distance(
                &text[ti.left_ltr_start..=ti.left_ltr_end],
                &text[ti.right_ltr_start..=ti.right_ltr_end],
            );
            let post_sim = 100.0 * (1.0 - edist as f64 / llen.max(rlen) as f64);
            if post_sim < args.similar {
                return None;
            }

            cand.similarity      = post_sim;
            cand.left_ltr_start  = ti.left_ltr_start;
            cand.left_ltr_end    = ti.left_ltr_end;
            cand.right_ltr_start = ti.right_ltr_start;
            cand.right_ltr_end   = ti.right_ltr_end;
            Some((seed_idx, cand, ti))
        })
        .collect();
    eprintln!("Phase 2: {} candidate(s) after extension + TSD + similarity filter", paired.len());

    // ---- Phase 3a: remove duplicate predictions ----
    // Sort by (seq_num, post-TSD left_ltr_start, post-TSD right_ltr_end, seed_idx).
    // The seed_idx tiebreaker makes the sort fully deterministic regardless of the
    // order parallel threads deposited results, and ensures the lowest-index seed
    // (first in SA-enumeration order) survives dedup — the same choice a
    // single-threaded stable sort on insertion order would make.
    paired.sort_unstable_by_key(|(idx, c, _)| (c.seq_num, c.left_ltr_start, c.right_ltr_end, *idx));
    paired.dedup_by_key(|(_, c, _)| (c.seq_num, c.left_ltr_start, c.right_ltr_end));

    // Strip the seed index — no longer needed.
    let mut paired: Vec<(ltr::LtrCandidate, tsd::TsdMotifInfo)> =
        paired.into_iter().map(|(_, c, ti)| (c, ti)).collect();
    eprintln!("Phase 3a: {} candidate(s) after duplicate removal", paired.len());

    // ---- Phase 3c: remove overlapping predictions (if requested) ----
    // Already sorted by (seq_num, left_ltr_start, right_ltr_end) from Phase 3a.
    // No re-sort needed.

    let cands_for_overlap: Vec<ltr::LtrCandidate> = paired.iter().map(|(c, _)| c.clone()).collect();
    let keep_indices = match args.overlaps.as_str() {
        "no" => {
            let idx = overlapping_keep_indices(&cands_for_overlap, true);
            eprintln!("Phase 3c: {} candidate(s) after overlap removal (no)", idx.len());
            idx
        }
        "best" => {
            let idx = overlapping_keep_indices(&cands_for_overlap, false);
            eprintln!("Phase 3c: {} candidate(s) after overlap removal (best)", idx.len());
            idx
        }
        _ => {
            let idx: Vec<usize> = (0..paired.len()).collect();
            eprintln!("Phase 3c: {} candidate(s) (overlaps=all, no removal)", idx.len());
            idx
        }
    };

    let surviving: Vec<(ltr::LtrCandidate, tsd::TsdMotifInfo)> = keep_indices
        .into_iter()
        .map(|i| paired[i].clone())
        .collect();

    eprintln!(
        "Phase 4: {} candidate(s) with TSD, {} with outer motif, {} with inner motif",
        surviving.iter().filter(|(_, t)| t.has_tsd).count(),
        surviving.iter().filter(|(_, t)| t.motif_near_tsd).count(),
        surviving.iter().filter(|(_, t)| t.motif_far_tsd).count(),
    );

    // ---- Output ----
    match args.format {
        OutputFormat::Gt => {
            println!(
                "# args=-index {} -seed {} -minlenltr {} -maxlenltr {} -mindistltr {} -maxdistltr {} -similar {} -xdrop {} -mintsd {} -maxtsd {} -overlaps {}",
                args.sufr,
                args.min_seed, args.ltr_min, args.ltr_max,
                args.dist_min, args.dist_max,
                args.similar as u32, args.xdrop,
                args.mintsd, args.maxtsd,
                args.overlaps,
            );
            println!("# predictions are reported in the following way");
            println!("# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr ");
            println!("# where:");
            println!("# s = starting position");
            println!("# e = ending position");
            println!("# l = length");
            println!("# ret = LTR-retrotransposon");
            println!("# lLTR = left LTR");
            println!("# rLTR = right LTR");
            println!("# sim = similarity");
            println!("# seq-nr = sequence number");

            for (c, ti) in &surviving {
                let seq_off = seq_starts.get(c.seq_num).copied().unwrap_or(0);
                let rel1 = |abs: usize| abs - seq_off + 1;
                println!(
                    "{}  {}  {}  {}  {}  {}  {}  {}  {}  {:.2}  {}",
                    rel1(ti.left_ltr_start),
                    rel1(ti.right_ltr_end),
                    ti.right_ltr_end - ti.left_ltr_start + 1,
                    rel1(ti.left_ltr_start),
                    rel1(ti.left_ltr_end),
                    ti.left_ltr_end - ti.left_ltr_start + 1,
                    rel1(ti.right_ltr_start),
                    rel1(ti.right_ltr_end),
                    ti.right_ltr_end - ti.right_ltr_start + 1,
                    c.similarity,
                    c.seq_num,
                );
            }
        }
    }

    Ok(())
}

// --------------------------------------------------
fn main() -> Result<()> {
    let args = Args::parse();
    let text_len = read_text_length(&args.sufr)?;
    if text_len < u32::MAX as usize {
        run::<u32>(&args)
    } else {
        run::<u64>(&args)
    }
}
