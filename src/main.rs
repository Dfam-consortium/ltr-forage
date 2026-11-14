use anyhow::{Context, Result};
use clap::Parser;
use libsufr::{
    sufr_file::SufrFile,
    types::{FromUsize, Int},
    util::read_text_length,
};
use rayon::prelude::*;

mod decide;
mod domain;
mod features;
mod internal;
mod ltr;
mod maxpairs;
mod pbs;
mod reextend;
mod stream;
mod fastan;
mod trna;
mod tsd;
mod xdrop;

// --------------------------------------------------
#[derive(Debug, Clone, clap::ValueEnum, PartialEq)]
enum OutputFormat {
    /// GenomeTools LTR_harvest compatible output (default)
    Gt,
    /// One tab-separated row per candidate with boundaries before and after
    /// re-extension, alignment details, TSD, motif, PPT and PBS candidates
    Tsv,
}

// --------------------------------------------------
#[derive(Parser, Debug)]
#[command(
    version,
    name = "ltr-forage",
    about = "Identify LTR retrotransposon candidates using a sufr suffix array index",
    allow_negative_numbers = true,
)]
struct Args {
    /// Path to a .sufr index file (created with `sufr create --dna`)
    #[arg(required_unless_present = "trna_scan")]
    sufr: Option<String>,

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

    /// LTRharvest mode: GenomeTools LTRharvest's search and nothing more.
    /// The TSD search moves each candidate's outer ends (--mintsd, --maxtsd,
    /// --vicinity, --motif, --motifmis, --tiered-tsd apply), overlaps are
    /// resolved by seed similarity (--overlaps best), and the output is gt's
    /// table; no tandem-array mask, re-extension, annotation, decisions or
    /// library
    #[arg(long)]
    ltrharvest: bool,

    /// Minimum TSD length; 0 disables the TSD search [default: 4]
    #[arg(long, value_name = "BP", requires = "ltrharvest")]
    mintsd: Option<usize>,

    /// Maximum TSD length [default: 20]
    #[arg(long, value_name = "BP", requires = "ltrharvest")]
    maxtsd: Option<usize>,

    /// Search window around the LTR boundaries for the TSD and motif [default: 60]
    #[arg(long, value_name = "BP", requires = "ltrharvest")]
    vicinity: Option<usize>,

    /// Require the TG...CA motif at the LTR ends
    #[arg(long, requires = "ltrharvest")]
    motif: bool,

    /// Allowed mismatches in the motif (0-3) [default: 0]
    #[arg(long, value_name = "N", requires = "ltrharvest")]
    motifmis: Option<u32>,

    /// Choose among TSD candidates LTR_retriever's way (tiers) instead of
    /// gt's (smallest boundary shift)
    #[arg(long, requires = "ltrharvest")]
    tiered_tsd: bool,

    /// After re-extension: treat candidates whose final outer ends both lie
    /// within this many bases of a better candidate's as duplicates (higher
    /// final LTR identity wins); seeds of one element stop their X-drop runs
    /// hundreds of bases apart but re-extend to the same ends; 0 disables
    #[arg(long, default_value_t = 20, value_name = "BP")]
    dedup_tol: usize,

    /// Overlapping candidates: evidence (keep all through annotation, then
    /// among overlapping kept candidates keep the one with the most evidence:
    /// domain, TSD, TG..CA and PBS, then LTR identity), best (LTRharvest's:
    /// before re-extension, keep the highest seed similarity of each group of
    /// overlapping candidates), no (drop every overlapping candidate), all
    /// (keep them all) [default: evidence; best with --ltrharvest]
    #[arg(long, value_enum, value_name = "STRATEGY")]
    overlaps: Option<Overlaps>,

    /// Number of threads (0 = all available CPUs)
    #[arg(short = 'j', long, default_value_t = 0, value_name = "N")]
    threads: usize,

    /// Output format on stdout [default: tsv; gt with --ltrharvest]
    #[arg(long, value_enum, value_name = "FORMAT")]
    format: Option<OutputFormat>,

    /// Longest LTR the re-extension window reaches beyond each X-drop end
    /// (bp).  7000 covers the longest LTR in the curated libraries we tested
    /// (rice, 6,144 bp); 2500 is enough for vertebrates and Arabidopsis and
    /// is faster on genomes with many tandem arrays (see the README)
    #[arg(long, default_value_t = 7000, value_name = "BP")]
    reextend_ltr_max: usize,

    /// Extra window beyond that on the outer side of each LTR (bp)
    #[arg(long, default_value_t = 100, value_name = "BP")]
    reextend_vicinity: usize,

    /// Re-extend with a window of this many bp first (0: off).  A
    /// candidate is searched again with the full --reextend-ltr-max window
    /// only when the first alignment fails or ends within
    /// --reextend-first-margin of the window edge with both LTRs below
    /// --tandem-max self-repeat; satellite pairs that run to the edge are
    /// vetoed on the short window's LTRs (see the README)
    #[arg(long, default_value_t = 2500, value_name = "BP")]
    reextend_first: usize,

    /// With --reextend-first: a first-pass alignment ending within this many
    /// bp of a window edge is searched again with the full window
    #[arg(long, default_value_t = 500, value_name = "BP")]
    reextend_first_margin: usize,

    /// Terminal alignment window that must be clean for an end to stand (columns)
    #[arg(long, default_value_t = 6, value_name = "N")]
    trim_window: usize,

    /// Matches required inside that terminal window
    #[arg(long, default_value_t = 6, value_name = "N")]
    trim_matches: usize,

    /// Count 20-mers from each internal region against the genome index and
    /// report their coverage and the longest single-copy stretch (tsv output);
    /// loads the suffix array into memory
    #[arg(long, default_value_t = false)]
    internal_stats: bool,

    /// HMMER3 file of transposable-element protein profiles; each
    /// candidate's six-frame translation is searched against it with nail
    /// (needs nail and MMseqs2, see --nail and --mmseqs) and the order,
    /// superfamily and strand reported.  The default is LTR_retriever's
    /// TEfam.hmm, shipped in data/
    #[arg(long, value_name = "FILE", default_value = DEFAULT_HMM)]
    domain_hmm: String,

    /// Skip the protein-domain search (decisions then rest on the TSD,
    /// TG..CA and PBS only)
    #[arg(long)]
    no_domains: bool,

    /// E-value threshold for a domain hit
    #[arg(long, default_value_t = 0.05, value_name = "E")]
    domain_evalue: f64,

    /// Bases added on each side of the element before translation
    #[arg(long, default_value_t = 50, value_name = "BP")]
    domain_flank: usize,

    /// nail binary for the domain search (a path, or a name on PATH)
    #[arg(long, default_value = "nail", value_name = "PATH")]
    nail: String,

    /// MMseqs2 binary nail seeds with (a path, or a name on PATH)
    #[arg(long, default_value = "mmseqs", value_name = "PATH")]
    mmseqs: String,

    /// Directory for scratch files (the domain search's six-frame proteins
    /// and nail's databases: about 2 GB for a primate genome); a
    /// subdirectory is made and removed [default: $TMPDIR, else /tmp]
    #[arg(long, value_name = "DIR")]
    temp_dir: Option<String>,

    /// Keep candidates with a polyprotein domain hit inside an LTR (by
    /// default they are vetoed: an LTR does not encode the polyprotein, so
    /// the hit means the "LTR" is internal sequence)
    #[arg(long)]
    keep_coding_ltrs: bool,

    /// Write the kept elements as a FASTA library: left LTR, internal
    /// region and right LTR per element, named like LTR_retriever's
    /// redundant library
    #[arg(long, value_name = "FILE", conflicts_with = "ltrharvest")]
    library: Option<String>,

    /// Evidence a candidate needs to be kept in the library and marked kept
    /// in the table; every level rejects candidates whose protein hits are
    /// not LTR-order
    #[arg(long, value_enum, default_value_t = decide::KeepLevel::Standard, value_name = "LEVEL")]
    keep: decide::KeepLevel,

    /// Also reject a candidate whose internal region has a single-copy
    /// stretch longer than this fraction of it (needs --internal-stats)
    #[arg(long, value_name = "FRAC")]
    internal_gap_max: Option<f64>,

    /// Reject a candidate when either LTR aligns to itself off the diagonal
    /// over at least this fraction of its length (satellite monomers, Alu
    /// dimers, microsatellites); values above 1 disable the check
    #[arg(long, default_value_t = 0.99, value_name = "FRAC")]
    tandem_max: f64,

    /// Do not reject candidates whose two LTRs both end in a poly-A tail
    /// without the TG..CA termini (non-LTR retroposon 3' ends)
    #[arg(long, default_value_t = false)]
    no_polya_veto: bool,

    /// Run the domain searches on candidates the LTR self-repeat check has
    /// already rejected too (their table domain columns are "." otherwise);
    /// the decisions do not change
    #[arg(long)]
    domain_search_vetoed: bool,

    /// Skip the genome tRNA scan behind the PBS columns of --format tsv
    /// (pbs_offset is then -1 and pbs_trna, pbs_mm12, pbs_edist18 are ".")
    #[arg(long)]
    no_pbs_trna: bool,

    /// Search for the PBS only on the strand the protein domains give
    /// (default: both strands; the domain strand follows the majority of
    /// hits and is wrong for some rice Gypsy families)
    #[arg(long)]
    pbs_domain_strand: bool,

    /// Write the genome's tRNA library to FILE (FASTA, one entry per distinct
    /// gene sequence with its copy count and PBS usage) and every tRNA call
    /// to FILE.tsv
    #[arg(long, value_name = "FILE", conflicts_with = "ltrharvest")]
    trna_out: Option<String>,

    /// Longest tRNA intron searched for (ARAGORN -i default 3000; 0 = none)
    #[arg(long, default_value_t = trna::MAXINTRONLEN, value_name = "BP")]
    trna_intron_max: i64,

    /// Only scan FASTA for tRNA genes (ARAGORN -t -gc1 -l -i port), print
    /// the calls as a table and exit; no SUFR index is needed
    #[arg(long, value_name = "FASTA")]
    trna_scan: Option<String>,

    /// Also mask the intervals in FILE (BED: name, 0-based start, end) from
    /// the seed search, e.g. tandem repeats from ULTRA or TRF
    #[arg(long, value_name = "FILE")]
    mask_bed: Option<String>,

    /// Search tandem arrays for seeds too.  By default arrays found by a port
    /// of Gene Myers' FasTAN (src/fastan.rs: any unit up to 30 kb, about two
    /// copies or more at roughly 70% identity or better) are hidden from the
    /// seed search; the sequence itself is unchanged.  An array of unit p
    /// yields a seed pair between almost every two units, so without the mask
    /// centromeric and satellite arrays dominate the run time
    #[arg(long, conflicts_with = "ltrharvest")]
    no_fastan: bool,

    /// Write the FasTAN-port arrays to FILE (BED: name, start, end, period,
    /// 1000 x (1 - differences / length), as FasTAN's tanbed writes them)
    #[arg(long, value_name = "FILE", conflicts_with_all = ["no_fastan", "ltrharvest"])]
    fastan_out: Option<String>,

    /// Only run the FasTAN port (writes --fastan-out) and exit
    #[arg(long, requires = "fastan_out")]
    fastan_only: bool,
}

const DEFAULT_HMM: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/data/TEfam.hmm");

#[derive(Debug, Clone, Copy, clap::ValueEnum, PartialEq)]
enum Overlaps {
    Evidence,
    Best,
    No,
    All,
}

/// Elapsed time as 41.2s or 9m21s.
fn el(t: &std::time::Instant) -> String {
    let s = t.elapsed().as_secs_f64();
    if s < 60.0 { format!("{:.1}s", s) } else { format!("{}m{:02}s", (s / 60.0) as u64, (s % 60.0) as u64) }
}

/// The run's settings after mode defaults, and every problem the run would
/// hit late (missing programs, unreadable profiles, unwritable paths),
/// reported together before the index is loaded.
struct Setup {
    overlaps: Overlaps,
    format: OutputFormat,
    nail: Option<domain::Nail>,
}

fn version_of(prog: &std::path::Path, flag: &str) -> Option<(u32, u32, u32)> {
    let out = std::process::Command::new(prog).arg(flag).output().ok()?;
    let txt = String::from_utf8_lossy(if out.stdout.is_empty() { &out.stderr } else { &out.stdout }).to_string();
    txt.split(|c: char| !(c.is_ascii_digit() || c == '.')).find_map(|w| {
        let v: Vec<u32> = w.split('.').filter_map(|x| x.parse().ok()).collect();
        (v.len() >= 3).then(|| (v[0], v[1], v[2]))
    })
}

/// A file can be made in `dir` (made and removed again).
fn writable(dir: &std::path::Path) -> bool {
    let probe = dir.join(format!(".ltr-forage-probe-{}", std::process::id()));
    let ok = std::fs::File::create(&probe).is_ok();
    std::fs::remove_file(&probe).ok();
    ok
}

fn setup(args: &Args) -> Result<Setup> {
    let mut bad: Vec<String> = Vec::new();
    let overlaps = args.overlaps.unwrap_or(if args.ltrharvest { Overlaps::Best } else { Overlaps::Evidence });
    let format = args.format.clone().unwrap_or(if args.ltrharvest { OutputFormat::Gt } else { OutputFormat::Tsv });
    if args.ltrharvest && overlaps == Overlaps::Evidence {
        bad.push("--overlaps evidence needs the annotation, which --ltrharvest skips".into());
    }
    if args.ltrharvest && format == OutputFormat::Tsv {
        bad.push("--format tsv reports the annotation, which --ltrharvest skips".into());
    }
    if let Some(p) = &args.sufr {
        if !std::path::Path::new(p).is_file() { bad.push(format!("index not found: {}", p)); }
    }
    let mut nail = None;
    if !args.ltrharvest && !args.no_domains {
        match std::fs::read_to_string(&args.domain_hmm) {
            Err(e) => bad.push(format!("cannot read --domain-hmm {}: {}", args.domain_hmm, e)),
            Ok(t) if !t.trim_start().starts_with("HMMER3") => bad.push(format!("--domain-hmm {} is not a HMMER3 profile file", args.domain_hmm)),
            Ok(_) => {}
        }
        let n = domain::find_program(&args.nail);
        if let Ok(p) = &n {
            match version_of(p, "--version") {
                Some(v) if v >= domain::NAIL_MIN => {}
                v => bad.push(format!("nail {} is {}; the domain search needs {}.{}.{} or newer", p.display(),
                                      v.map_or("of unknown version".to_string(), |v| format!("{}.{}.{}", v.0, v.1, v.2)),
                                      domain::NAIL_MIN.0, domain::NAIL_MIN.1, domain::NAIL_MIN.2)),
            }
        }
        let m = domain::find_program(&args.mmseqs);
        if let Ok(p) = &m {
            if !std::process::Command::new(p).arg("version").output().map_or(false, |o| o.status.success()) {
                bad.push(format!("MMseqs2 {} does not run (`{} version` failed)", p.display(), p.display()));
            }
        }
        let tmp = args.temp_dir.as_ref().map(std::path::PathBuf::from).unwrap_or_else(std::env::temp_dir);
        if !tmp.is_dir() || !writable(&tmp) {
            bad.push(format!("temporary directory {} is missing or not writable (--temp-dir)", tmp.display()));
        }
        match (n, m) {
            (Ok(n), Ok(m)) => {
                nail = Some(domain::Nail { nail: n, mmseqs: m, threads: 0, work: tmp.join(format!("ltr-forage-nail-{}", std::process::id())) })
            }
            (n, m) => {
                for e in [n.err(), m.err()].into_iter().flatten() {
                    bad.push(format!("the domain search needs nail and MMseqs2 (--nail, --mmseqs, or --no-domains): {}", e));
                }
            }
        }
    }
    for (opt, p) in [("--library", &args.library), ("--fastan-out", &args.fastan_out), ("--trna-out", &args.trna_out)] {
        if let Some(p) = p {
            let dir = std::path::Path::new(p).parent().filter(|d| !d.as_os_str().is_empty()).unwrap_or(std::path::Path::new("."));
            if !dir.is_dir() || !writable(dir) {
                bad.push(format!("{} {}: directory {} is missing or not writable", opt, p, dir.display()));
            }
        }
    }
    if !bad.is_empty() {
        anyhow::bail!("cannot start:\n  {}", bad.join("\n  "));
    }
    Ok(Setup { overlaps, format, nail })
}

/// BED intervals as half-open text coordinates, clipped to their sequence;
/// returns the intervals, their total length and the lines naming unknown sequences
fn read_bed(path: &str, seq_names: &[String], seq_starts: &[usize], text_len: usize) -> Result<(Vec<(usize, usize)>, usize, usize)> {
    let idx: std::collections::HashMap<&str, usize> = seq_names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let (mut out, mut bp, mut unknown) = (Vec::new(), 0usize, 0usize);
    for line in std::io::BufRead::lines(std::io::BufReader::new(std::fs::File::open(path)?)) {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 || line.starts_with('#') || line.starts_with("track") { continue; }
        let (Ok(a), Ok(b)) = (f[1].parse::<usize>(), f[2].parse::<usize>()) else { continue };
        let Some(&k) = idx.get(f[0]) else { unknown += 1; continue };
        let end = if k + 1 < seq_starts.len() { seq_starts[k + 1] } else { text_len };
        let (a, b) = (seq_starts[k] + a, (seq_starts[k] + b).min(end));
        if a < b { out.push((a, b)); bp += b - a; }
    }
    Ok((out, bp, unknown))
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
fn run<T>(args: &Args, st: &Setup) -> Result<()>
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

    // Loads the text (needed for random access throughout); the suffix array
    // and LCP array stay on disk and are streamed once in phase 1.
    let t_start = std::time::Instant::now();
    let sufr = args.sufr.as_deref().context("a SUFR index is required")?;
    let mut sf: SufrFile<T> = SufrFile::read(sufr, false)?;

    // sufr stores whole FASTA header lines; output names are the first word
    let seq_names: Vec<String> = sf
        .sequence_names
        .iter()
        .map(|n| n.split_whitespace().next().unwrap_or("").to_string())
        .collect();
    let text: &[u8] = &sf.text;
    let num_suffixes = sf.len_suffixes.to_usize();
    let seq_starts: Vec<usize> = sf.sequence_starts.iter().map(|v| v.to_usize()).collect();
    let width = std::mem::size_of::<T>();
    let sa_pos = sf.suffix_array_pos as u64;
    let lcp_pos = sf.lcp_pos as u64;

    eprintln!(
        "[1 index] {} sequence(s), {} suffixes, text length {}; text loaded in {}",
        sf.num_sequences.to_usize(),
        num_suffixes,
        sf.text_len.to_usize(),
        el(&t_start),
    );
    eprintln!(
        "[1 index] {}parameters: min_seed={} ltr=[{},{}] dist=[{},{}] similar={:.0}% xdrop={}",
        if args.ltrharvest { "LTRharvest mode; " } else { "" },
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

    // ---- 2, 3: tandem mask, seed pairs ----
    // The suffix array is cut at positions whose LCP is below the seed length
    // (see stream::split_points); each piece is traversed by its own worker,
    // streaming its slice of the SA and LCP arrays from the file.  Seeds are
    // concatenated in piece order, which is the order a single sequential
    // traversal would emit them, so later tie-breaks are unaffected.
    let mut extra: Vec<(usize, usize)> = Vec::new();
    if let Some(path) = &args.mask_bed {
        let (iv, bp, unknown) = read_bed(path, &seq_names, &seq_starts, text.len())?;
        let n = iv.len();
        extra.extend(iv);
        eprintln!("[2 tandem mask] BED: {} interval(s), {:.1} Mb{}", n, bp as f64 / 1e6,
                  if unknown > 0 { format!(" ({} on sequences not in the index, ignored)", unknown) } else { String::new() });
    }
    if !args.ltrharvest && !args.no_fastan {
        let t0 = std::time::Instant::now();
        let arrays = fastan::tandem_arrays(text, &seq_starts);
        let bp: usize = arrays.iter().map(|a| a.aend - a.bbeg).sum();
        eprintln!("[2 tandem mask] FasTAN port: {} array(s), {:.1} Mb before merging ({})", arrays.len(), bp as f64 / 1e6, el(&t0));
        if let Some(path) = &args.fastan_out {
            use std::io::Write;
            let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
            for a in &arrays {
                // as tanbed reports it
                let id = (1000.0 * (1.0 - a.diffs as f64 / (a.aend - a.bbeg).max(1) as f64)).max(0.0) as i64;
                writeln!(w, "{}\t{}\t{}\t{}\t{}", seq_names[a.seq], a.bbeg, a.aend, a.period, id)?;
            }
            w.flush()?;
        }
        if args.fastan_only {
            std::process::exit(0);
        }
        extra.extend(arrays.iter().map(|a| (seq_starts[a.seq] + a.bbeg, seq_starts[a.seq] + a.aend)));
    }
    let wild = maxpairs::WildRuns::with_extra(text, &extra);
    let n_pieces = rayon::current_num_threads() * 4;
    let pieces = stream::split_points(sufr, lcp_pos, num_suffixes, width, args.min_seed, n_pieces)?;
    let per_piece: Vec<Result<Vec<ltr::Seed>>> = pieces
        .par_iter()
        .map(|&(start, end)| {
            let n = end - start;
            let sa = stream::ArrayReader::open(sufr, sa_pos + (start * width) as u64, n, width)?;
            let lcp = stream::ArrayReader::open(sufr, lcp_pos + (start * width) as u64, n, width)?;
            let mut seeds: Vec<ltr::Seed> = Vec::new();
            maxpairs::enumerate_maximal_repeats(
                sa,
                lcp,
                &text,
                n,
                args.min_seed,
                &wild,
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
            Ok(seeds)
        })
        .collect();
    let mut seeds: Vec<ltr::Seed> = Vec::new();
    for piece in per_piece {
        seeds.extend(piece?);
    }
    eprintln!("[3 seed pairs] {} seed pair(s) ({} piece(s); {} elapsed)", seeds.len(), pieces.len(), el(&t_start));
    // ---- 4: extend, TSD/motif, constraints, similarity (GT's order) ----
    // GT runs these per seed inside gt_searchforLTRs: x-drop extension, then the
    // TSD/motif search (which moves the outer boundaries), then the length and
    // distance check on the moved boundaries, then the edit distance.  The sort
    // and dedup that follow operate on the post-TSD outer boundaries.
    // Outside --ltrharvest the TSD search runs with mintsd 0: it records a TSD (the
    // xdrop_tsd_len column) but neither moves the ends nor rejects; the TSD the
    // decisions use is searched at the final ends (features::tsd_or_mm1).
    let (mintsd, maxtsd, vicinity) = if args.ltrharvest {
        (args.mintsd.unwrap_or(4), args.maxtsd.unwrap_or(20), args.vicinity.unwrap_or(60))
    } else {
        (0, 20, 60)
    };
    let motif_cfg = tsd::MotifConfig {
        firstleft:          b'T',
        secondleft:         b'G',
        firstright:         b'C',
        secondright:        b'A',
        allowed_mismatches: args.motifmis.unwrap_or(0),
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
            let mut cand = ltr::extend_seed(seed, &text, &seq_starts, args.ltr_max, &scores);

            // GT (gt_searchforLTRs): when TSD and/or motif search is on, a
            // candidate survives only if a TSD was found (and the motifs, when
            // requested), or, with TSD search off, if both motifs were found.
            let ti = tsd::annotate_tsd_and_motif(
                &cand, &text, &seq_starts,
                mintsd, maxtsd, vicinity,
                &motif_cfg, args.motif, args.tiered_tsd,
            );
            if mintsd >= 2 || args.motif {
                let full_with_tsd = ti.has_tsd
                    && (!args.motif || (ti.motif_near_tsd && ti.motif_far_tsd));
                let full_motif_only = mintsd < 2
                    && ti.motif_near_tsd && ti.motif_far_tsd;
                if !full_with_tsd && !full_motif_only {
                    return None;
                }
            }

            // Length and distance constraints on the adjusted boundaries.
            if !ltr::constraints_ok(
                ti.left_ltr_start, ti.left_ltr_end,
                ti.right_ltr_start, ti.right_ltr_end,
                args.ltr_min, args.ltr_max, args.dist_min, args.dist_max,
            ) {
                return None;
            }

            // Similarity on the adjusted boundaries, then the threshold.
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
    eprintln!(
        "[4 X-drop extension] {} candidate(s) {} ({} elapsed)",
        paired.len(),
        if args.ltrharvest { "after extension, TSD/motif search and similarity" } else { "at the similarity threshold" },
        el(&t_start),
    );

    // ---- 5: remove duplicate predictions ----
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
    eprintln!("[5 duplicates] {} candidate(s) after removing identical ends", paired.len());

    // ---- 10 (best, no): remove overlapping predictions before re-extension ----
    // Already sorted by (seq_num, left_ltr_start, right_ltr_end) by step 5.
    // No re-sort needed.

    let cands_for_overlap: Vec<ltr::LtrCandidate> = paired.iter().map(|(c, _)| c.clone()).collect();
    let keep_indices = match st.overlaps {
        Overlaps::No | Overlaps::Best => {
            let no = st.overlaps == Overlaps::No;
            let idx = overlapping_keep_indices(&cands_for_overlap, no);
            eprintln!("[10 overlaps] {} candidate(s) after overlap removal ({}, before re-extension)", idx.len(), if no { "no" } else { "best" });
            idx
        }
        Overlaps::All | Overlaps::Evidence => (0..paired.len()).collect(),
    };

    let surviving: Vec<(ltr::LtrCandidate, tsd::TsdMotifInfo)> = keep_indices
        .into_iter()
        .map(|i| paired[i].clone())
        .collect();

    if args.ltrharvest {
        print_gt(args, st, sufr, &seq_starts, &surviving, &vec![None; surviving.len()]);
        return Ok(());
    }

    // ---- 6: re-extension ----
    // A local alignment of the two LTR windows replaces the X-drop
    // boundaries, which stop at the first short indel between the copies.
    // Candidates whose alignment fails keep their X-drop boundaries.
    let reext: Vec<Option<reextend::Reextended>> = {
        let aligner = reextend::Aligner::new()?;
        let cfg = reextend::ReextendConfig {
            ltr_max: args.reextend_ltr_max,
            vicinity: args.reextend_vicinity,
            trim_window: args.trim_window,
            trim_matches: args.trim_matches,
            edge_margin: 5,
        };
        // a first-pass alignment that stops short of the edge may still be cut by it (X-drop and trim pull
        // the end back from where the window ends), so 'at the edge' gets a wide margin there
        let cfg1 = reextend::ReextendConfig { ltr_max: args.reextend_first, edge_margin: args.reextend_first_margin, ..cfg.clone() };
        let two_pass = args.reextend_first > 0 && args.reextend_first < args.reextend_ltr_max;
        let n_second = std::sync::atomic::AtomicUsize::new(0);
        let out: Vec<Option<reextend::Reextended>> = surviving
            .par_iter()
            .map(|(c, ti)| {
                let seq_start = seq_starts[c.seq_num];
                let seq_end = seq_starts
                    .get(c.seq_num + 1)
                    .map(|&s| s - 1)
                    .unwrap_or(text.len() - 1);
                let run = |cf: &reextend::ReextendConfig| reextend::reextend(
                    &aligner, cf, &text, seq_start, seq_end,
                    ti.left_ltr_start, ti.left_ltr_end, ti.right_ltr_start, ti.right_ltr_end,
                );
                if !two_pass {
                    return run(&cfg);
                }
                // Short window first.  An alignment that ends inside it stands; one that runs to the
                // edge is a long LTR or a satellite array.  Arrays are vetoed by 8a on these LTRs
                // (their self-repeat coverage is kept for it); everything else gets the full window.
                match run(&cfg1) {
                    Some(r) if !r.at_window_edge => Some(r),
                    Some(mut r) => {
                        let a = reextend::self_repeat_cov(&aligner, &text[r.left_ltr_start..=r.left_ltr_end]);
                        let b = reextend::self_repeat_cov(&aligner, &text[r.right_ltr_start..=r.right_ltr_end]);
                        if a.max(b) >= args.tandem_max {
                            r.self_cov = Some((a, b));
                            Some(r)
                        } else {
                            n_second.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                            run(&cfg)
                        }
                    }
                    None => {
                        n_second.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        run(&cfg)
                    }
                }
            })
            .collect();
        if two_pass {
            eprintln!(
                "[6 re-extension] {} bp window first; {} candidate(s) searched again with {} bp",
                args.reextend_first, n_second.into_inner(), args.reextend_ltr_max,
            );
        }
        let n_ok = out.iter().filter(|r| r.is_some()).count();
        let n_edge = out.iter().filter(|r| r.as_ref().map_or(false, |x| x.at_window_edge)).count();
        let n_moved = out
            .iter()
            .zip(surviving.iter())
            .filter(|(r, (_, ti))| r.as_ref().map_or(false, |x| {
                x.left_ltr_start != ti.left_ltr_start || x.left_ltr_end != ti.left_ltr_end
                    || x.right_ltr_start != ti.right_ltr_start || x.right_ltr_end != ti.right_ltr_end
            }))
            .count();
        eprintln!(
            "[6 re-extension] re-extended {} of {} candidate(s), {} moved, {} at a window edge ({} elapsed)",
            n_ok, surviving.len(), n_moved, n_edge, el(&t_start),
        );
        out
    };

    // ---- 7: near duplicates after re-extension ----
    let (surviving, reext) = if args.dedup_tol > 0 {
        let tol = args.dedup_tol;
        let ends: Vec<(usize, usize, usize, f64)> = surviving
            .iter()
            .zip(reext.iter())
            .map(|((c, ti), re)| match re {
                Some(r) => (c.seq_num, r.left_ltr_start, r.right_ltr_end, r.identity),
                None => (c.seq_num, ti.left_ltr_start, ti.right_ltr_end, c.similarity),
            })
            .collect();
        let mut order: Vec<usize> = (0..ends.len()).collect();
        order.sort_by(|&a, &b| ends[b].3.partial_cmp(&ends[a].3).unwrap().then(a.cmp(&b)));
        let mut claimed: std::collections::HashMap<usize, std::collections::BTreeMap<usize, Vec<usize>>> = Default::default();
        let mut keep = vec![false; ends.len()];
        for i in order {
            let (sq, ls, re_, _) = ends[i];
            let by_start = claimed.entry(sq).or_default();
            let near = by_start
                .range(ls.saturating_sub(tol)..=ls + tol)
                .any(|(_, v)| v.iter().any(|&e| e.abs_diff(re_) <= tol));
            if !near {
                by_start.entry(ls).or_default().push(re_);
                keep[i] = true;
            }
        }
        let before = surviving.len();
        let (s, r): (Vec<_>, Vec<_>) = surviving.into_iter().zip(reext).zip(keep).filter(|(_, k)| *k).map(|(p, _)| p).unzip();
        eprintln!("[7 near-duplicate merge] {} of {} candidate(s) left after merging final ends within {} bp", s.len(), before, tol);
        (s, r)
    } else {
        (surviving, reext)
    };

    // ---- 7b: internal-region statistics (optional) ----
    let final_coords: Vec<(usize, usize, usize, usize)> = surviving
        .iter()
        .zip(reext.iter())
        .map(|((_, ti), re)| match re {
            Some(r) => (r.left_ltr_start, r.left_ltr_end, r.right_ltr_start, r.right_ltr_end),
            None => (ti.left_ltr_start, ti.left_ltr_end, ti.right_ltr_start, ti.right_ltr_end),
        })
        .collect();
    let internal_stats: Vec<Option<internal::InternalStats>> = if args.internal_stats {
        let regions: Vec<Vec<Vec<u8>>> = final_coords
            .par_iter()
            .map(|&(_, l1e, l2s, _)| if l2s > l1e + 1 { internal::tiles(text, l1e + 1, l2s - 1) } else { Vec::new() })
            .collect();
        let n_tiles: usize = regions.iter().map(|r| r.len()).sum();
        let counts = internal::count_tiles(&mut sf, &regions)?;
        eprintln!(
            "[7b internal 20-mers] counted {} tile(s), {} distinct k-mer(s) with reverse complements ({} elapsed)",
            n_tiles, counts.len(), el(&t_start),
        );
        regions.iter().map(|r| Some(internal::stats(r, &counts))).collect()
    } else {
        vec![None; surviving.len()]
    };
    let text: &[u8] = &sf.text;

    // ---- 8a: LTR self-alignment (tandem-repeat check) ----
    // First, so the domain searches below can skip what it vetoes (decide()
    // rejects those whatever their domains).
    let self_cov: Vec<(f64, f64)> = {
        let aligner = reextend::Aligner::new()?;
        final_coords
            .par_iter()
            .zip(reext.par_iter())
            .map(|(&(l1s, l1e, l2s, l2e), re)| {
                // computed in step 6 for the same LTRs when the short window ran to its edge
                if let Some(c) = re.as_ref().and_then(|r| r.self_cov) {
                    return c;
                }
                (reextend::self_repeat_cov(&aligner, &text[l1s..=l1e]), reextend::self_repeat_cov(&aligner, &text[l2s..=l2e]))
            })
            .collect()
    };
    eprintln!(
        "[8a self-alignment] {} candidate(s) with an LTR at least {:.0}% self-repeat ({} elapsed)",
        self_cov.iter().filter(|(a, b)| a.max(*b) >= args.tandem_max).count(),
        100.0 * args.tandem_max,
        el(&t_start),
    );

    let skip_domains: Vec<bool> = self_cov
        .iter()
        .map(|(a, b)| !args.domain_search_vetoed && a.max(*b) >= args.tandem_max)
        .collect();

    // ---- 8b: protein-domain search ----
    // One nail run over each candidate (element plus --domain-flank) and,
    // unless --keep-coding-ltrs, each LTR alone: hits to the polyprotein
    // profiles inside an LTR mean the "LTR" is internal sequence.
    let n = final_coords.len();
    let (domain_calls, ltr_domain_hits): (Vec<Option<domain::DomainCall>>, Vec<Vec<String>>) = if let Some(nl) = &st.nail {
        let t0 = std::time::Instant::now();
        let path = &args.domain_hmm;
        let nail = domain::Nail { nail: nl.nail.clone(), mmseqs: nl.mmseqs.clone(), threads: rayon::current_num_threads(), work: nl.work.clone() };
        let origins: Vec<usize> = final_coords
            .iter()
            .zip(surviving.iter())
            .map(|(&(l1s, _, _, _), (c, _))| l1s.saturating_sub(args.domain_flank).max(seq_starts[c.seq_num]))
            .collect();
        let mut seqs: Vec<Option<&[u8]>> = Vec::with_capacity(3 * n);
        for ((&(_, _, _, l2e), (c, _)), (&s, &skip)) in final_coords.iter().zip(surviving.iter()).zip(origins.iter().zip(skip_domains.iter())) {
            let seq_end = seq_starts.get(c.seq_num + 1).map(|&s| s - 1).unwrap_or(text.len() - 1);
            seqs.push((!skip).then(|| &text[s..=(l2e + args.domain_flank).min(seq_end)]));
        }
        let ltrs = !args.keep_coding_ltrs;
        if ltrs {
            for (&(l1s, l1e, l2s, l2e), &skip) in final_coords.iter().zip(skip_domains.iter()) {
                seqs.push((!skip).then(|| &text[l1s..=l1e]));
                seqs.push((!skip).then(|| &text[l2s..=l2e]));
            }
        }
        let mut calls = domain::search(&nail, std::path::Path::new(path), &seqs, 6 * n.max(1), args.domain_evalue)?;
        let ltr_hits: Vec<Vec<String>> = if ltrs {
            calls[n..]
                .chunks(2)
                .map(|p| p.iter().flatten().flat_map(|c| c.ltr_order_hits()).collect())
                .collect()
        } else {
            vec![Vec::new(); n]
        };
        calls.truncate(n);
        for (c, &o) in calls.iter_mut().zip(origins.iter()) {
            if let Some(c) = c {
                c.origin = o;
            }
        }
        let n_skip = skip_domains.iter().filter(|&&x| x).count();
        if n_skip > 0 {
            eprintln!("[8b protein domains] skipping {} candidate(s) vetoed by LTR self-repeat (--domain-search-vetoed searches them)", n_skip);
        }
        let with_hits = calls.iter().filter(|d| d.as_ref().map_or(false, |d| !d.hits.is_empty())).count();
        let ltr = calls.iter().filter(|d| d.as_ref().map_or(false, |d| d.order == "LTR")).count();
        let not = calls.iter().filter(|d| d.as_ref().map_or(false, |d| d.order == "notLTR" || d.order == "mixture")).count();
        eprintln!(
            "[8b protein domains] nail search ({}): {} candidate(s) with a domain hit, {} LTR order, {} notLTR or mixture ({} elapsed)",
            el(&t0), with_hits, ltr, not, el(&t_start),
        );
        if ltrs {
            eprintln!(
                "[8b protein domains] {} candidate(s) with an LTR-list domain hit inside an LTR, {} with a core polyprotein domain",
                ltr_hits.iter().filter(|v| !v.is_empty()).count(),
                ltr_hits.iter().filter(|v| v.iter().any(|h| domain::is_core(h))).count(),
            );
        }
        (calls, ltr_hits)
    } else {
        (vec![None; n], vec![Vec::new(); n])
    };

    // ---- 8c: genome tRNA scan and PBS annotation ----
    let want_pbs = st.format == OutputFormat::Tsv && !args.no_pbs_trna;
    let (pbs_calls, trna_calls): (Vec<(Option<pbs::PbsHit>, Option<usize>)>, Option<(Vec<Vec<trna::TrnaCall>>, usize)>) = if want_pbs || args.trna_out.is_some() {
        let t0 = std::time::Instant::now();
        let ranges: Vec<(usize, usize)> = (0..seq_starts.len())
            .map(|i| (seq_starts[i], seq_starts.get(i + 1).map(|&s| s - 1).unwrap_or(text.len() - 1)))
            .collect();
        let calls: Vec<Vec<trna::TrnaCall>> = ranges
            .par_iter()
            .map(|&(a, e)| trna::scan_sequence(&text[a..e], args.trna_intron_max))
            .collect();
        let n_calls: usize = calls.iter().map(|c| c.len()).sum();
        let all: Vec<trna::TrnaCall> = calls.iter().flatten().cloned().collect();
        let probes = pbs::Probes::from_calls(&all);
        let out: Vec<(Option<pbs::PbsHit>, Option<usize>)> = if want_pbs || args.trna_out.is_some() {
            surviving
                .par_iter()
                .zip(final_coords.par_iter())
                .zip(domain_calls.par_iter())
                .map(|(((c, _), &(_, l1e, l2s, _)), dc)| {
                    let lo = seq_starts[c.seq_num];
                    let hi = seq_starts.get(c.seq_num + 1).map(|&s| s - 2).unwrap_or(text.len() - 2);
                    let strand = if args.pbs_domain_strand { dc.as_ref().map(|d| d.strand.as_str()).unwrap_or("?") } else { "?" };
                    probes.annotate(text, l1e, l2s, lo, hi, strand)
                })
                .collect()
        } else {
            vec![(None, None); surviving.len()]
        };
        eprintln!(
            "[8c tRNA → PBS] {} tRNA gene(s), {} PBS probe(s) ({}); {} candidate(s) with a 12-nt PBS, {} exact ({} elapsed)",
            n_calls, probes.len(), el(&t0),
            out.iter().filter(|(h, _)| h.is_some()).count(),
            out.iter().filter(|(h, _)| h.as_ref().map_or(false, |h| h.mm12 == 0)).count(),
            el(&t_start),
        );
        (out, Some((calls, n_calls)))
    } else {
        (vec![(None, None); surviving.len()], None)
    };


    // ---- 8d, 9, 10: TSD and TG..CA at the final ends, decisions, overlaps ----
    if args.internal_gap_max.is_some() && !args.internal_stats {
        eprintln!("Warning: --internal-gap-max has no effect without --internal-stats");
    }
    if st.nail.is_none() && args.keep == decide::KeepLevel::Domain {
        eprintln!("Warning: --keep domain keeps nothing without the domain search");
    }
    let identity: Vec<f64> = surviving
        .iter()
        .zip(reext.iter())
        .map(|((c, _), re)| re.as_ref().map_or(c.similarity, |r| r.identity))
        .collect();
    let ends: Vec<(features::Tsd, String, bool)> = surviving
        .par_iter()
        .zip(final_coords.par_iter())
        .map(|((c, _), &(l1s, l1e, l2s, l2e))| {
            let seq_off = seq_starts[c.seq_num];
            let seq_end = seq_starts.get(c.seq_num + 1).map(|&s| s - 1).unwrap_or(text.len() - 1);
            (features::tsd_or_mm1(text, l1s, l2e, seq_off, seq_end, 4, 20, 1), features::motif_at(text, l1s, l2e),
             features::polya_ltrs(text, (l1s, l1e, l2s, l2e)))
        })
        .collect();
    eprintln!(
        "[8d TSD, TG..CA] {} candidate(s) with a TSD of 4-20 bp ({} with one mismatch), {} with TG..CA ends",
        ends.iter().filter(|(t, _, _)| t.len > 0).count(),
        ends.iter().filter(|(t, _, _)| t.len > 0 && t.mismatches > 0).count(),
        ends.iter().filter(|(_, m, _)| features::motif_class(m) == "TGCA").count(),
    );
    let mut decisions: Vec<decide::Decision> = ends
        .iter()
        .zip(domain_calls.iter())
        .zip(internal_stats.iter())
        .zip(self_cov.iter())
        .zip(ltr_domain_hits.iter())
        .map(|(((((tsd, motif, polya), dc), ist), &(sc_l, sc_r)), ldh)| {
            decide::decide(
                args.keep,
                args.internal_gap_max,
                args.tandem_max,
                &decide::Evidence { domain: dc.as_ref(), tsd_len: tsd.len, tsd_mm: tsd.mismatches, motif, gap_frac: ist.as_ref().map(|s| s.gap_frac), self_cov: sc_l.max(sc_r), ltr_domain_hits: ldh.iter().filter(|h| domain::is_core(h)).count(), polya: *polya && !args.no_polya_veto },
            )
        })
        .collect();
    let kept = decisions.iter().filter(|d| d.keep).count();
    let vetoed = decisions.iter().filter(|d| d.why.starts_with("veto:notLTR") || d.why.starts_with("veto:mixture")).count();
    let tandem = decisions.iter().filter(|d| d.why.starts_with("veto:tandem")).count();
    let coding = decisions.iter().filter(|d| d.why.starts_with("veto:ltrdomain")).count();
    let polya = decisions.iter().filter(|d| d.why.starts_with("veto:polyA")).count();
    eprintln!(
        "[9 decision] kept {} of {} candidate(s) at --keep {:?}, {} vetoed by protein hits, {} by LTR self-repeat, {} by a domain inside an LTR, {} by poly-A-tailed LTRs ({} elapsed)",
        kept, decisions.len(), args.keep, vetoed, tandem, coding, polya, el(&t_start),
    );
    if st.overlaps == Overlaps::Evidence {
        // among overlapping kept candidates, the one with the most evidence; a candidate rejected on its own
        // evidence never blocks one that was kept
        let score = |i: usize| {
            let w = &decisions[i].why;
            (w.contains("domain") as u8) + (w.contains("tsd") as u8) + (w.contains("tgca") as u8) + (pbs_calls[i].0.is_some() as u8)
        };
        let mut order: Vec<usize> = (0..decisions.len()).filter(|&i| decisions[i].keep).collect();
        order.sort_by(|&a, &b| score(b).cmp(&score(a)).then(identity[b].partial_cmp(&identity[a]).unwrap()).then(a.cmp(&b)));
        let mut taken: std::collections::HashMap<usize, std::collections::BTreeMap<usize, usize>> = Default::default();
        let longest = final_coords.iter().map(|&(a, _, _, d)| d - a).max().unwrap_or(0);
        let mut dropped = 0;
        for i in order {
            let (sq, (s0, _, _, e0)) = (surviving[i].0.seq_num, final_coords[i]);
            let spans = taken.entry(sq).or_default();
            let hit = spans.range(s0.saturating_sub(longest)..=e0).any(|(&s, &e)| s <= e0 && s0 <= e);
            if hit {
                decisions[i].keep = false;
                decisions[i].why.push_str(",overlap");
                dropped += 1;
            } else {
                spans.insert(s0, e0);
            }
        }
        eprintln!("[10 overlaps] {} kept candidate(s) dropped as overlapping one with more evidence; {} left", dropped, decisions.iter().filter(|d| d.keep).count());
    }
    if let (Some(path), Some((calls, _))) = (&args.trna_out, &trna_calls) {
        write_trna_library(path, &seq_names, calls, &pbs_calls, &decisions)?;
    }
    if let Some(path) = &args.library {
        let mut out = std::io::BufWriter::new(std::fs::File::create(path).with_context(|| format!("creating {}", path))?);
        for (((c, _), &ltrs), (d, dc)) in surviving.iter().zip(final_coords.iter()).zip(decisions.iter().zip(domain_calls.iter())) {
            if !d.keep {
                continue;
            }
            let name = seq_names.get(c.seq_num).map(|s| s.as_str()).unwrap_or("?");
            let sf = dc.as_ref().map(|d| d.superfamily.as_str()).unwrap_or("unknown");
            decide::write_element(&mut out, name, text, seq_starts[c.seq_num], ltrs, sf, &d.why)?;
        }
        std::io::Write::flush(&mut out)?;
        eprintln!("[11 output] wrote {} element(s) to {}", decisions.iter().filter(|d| d.keep).count(), path);
    }

    // ---- Output ----
    match st.format {
        OutputFormat::Gt => print_gt(args, st, sufr, &seq_starts, &surviving, &reext),
        OutputFormat::Tsv => {
            println!(
                "seq\telem_start\telem_end\telem_len\tlltr_start\tlltr_end\tlltr_len\trltr_start\trltr_end\trltr_len\tidentity\t\
                 xdrop_lltr_start\txdrop_lltr_end\txdrop_rltr_start\txdrop_rltr_end\txdrop_similarity\t\
                 reextended\taln_score\taln_columns\ttrim_l5\ttrim_r5\ttrim_l3\ttrim_r3\twindow_edge\t\
                 tsd_len\ttsd_seq\ttsd_shift_l\ttsd_shift_r\tmotif\tmotif_class\tppt_purines15\tppt_offset\tpbs_offset\tn_count\txdrop_tsd_len\t\
                 int_tiles\tint_multi_frac\tint_gap_frac\tint_gap_bp\t\
                 dom_order\tdom_superfamily\tdom_strand\tdom_best_e\tdom_hits\tltr_dom_hits\tlltr_self_cov\trltr_self_cov\tkeep\tevidence\t\
                 pbs_trna\tpbs_mm12\tpbs_edist18\tpbs_strand\ttsd_mm\tdom_pos\tdom_nest"
            );
            for ((((((((c, ti), re), ist), dc), dec), &(sc_l, sc_r)), ldh), (pbs_hit, pbs_ed)) in surviving.iter().zip(reext.iter()).zip(internal_stats.iter()).zip(domain_calls.iter()).zip(decisions.iter()).zip(self_cov.iter()).zip(ltr_domain_hits.iter()).zip(pbs_calls.iter()) {
                let seq_off = seq_starts.get(c.seq_num).copied().unwrap_or(0);
                let seq_end = seq_starts.get(c.seq_num + 1).map(|&s| s - 1).unwrap_or(text.len() - 1);
                let rel1 = |abs: usize| abs - seq_off + 1;
                let (l1s, l1e, l2s, l2e, sim) = match re {
                    Some(r) => (r.left_ltr_start, r.left_ltr_end, r.right_ltr_start, r.right_ltr_end, r.identity),
                    None => (ti.left_ltr_start, ti.left_ltr_end, ti.right_ltr_start, ti.right_ltr_end, c.similarity),
                };
                let tsd = features::tsd_or_mm1(&text, l1s, l2e, seq_off, seq_end, 4, 20, 1);
                let motif = features::motif_at(&text, l1s, l2e);
                let (ppt_n, ppt_off) = features::ppt(&text, l2s, seq_off);
                let pbs_off = pbs_hit.as_ref().map_or(-1, |h| h.offset);
                let (pbs_trna, pbs_mm, pbs_strand) = match pbs_hit {
                    Some(h) => (h.trnas.join("/"), h.mm12.to_string(), if h.minus { "-" } else { "+" }),
                    None => (".".to_string(), ".".to_string(), "."),
                };
                let pbs_ed = pbs_ed.map_or(".".to_string(), |d| d.to_string());
                let ns = features::n_count(&text, l1s, l2e);
                let name = seq_names.get(c.seq_num).map(|s| s.as_str()).unwrap_or("?");
                let (reextended, score, cols, tr, edge) = match re {
                    Some(r) => (1, r.score, r.columns, r.trimmed, r.at_window_edge as u8),
                    None => (0, 0, 0, (0, 0, 0, 0), 0),
                };
                let (dord, dsf, dstr, dbe, dhits) = match dc {
                    Some(d) => {
                        let best = d.hits.iter().map(|h| h.e_value).fold(f64::INFINITY, f64::min);
                        let list = if d.hits.is_empty() { ".".to_string() } else {
                            d.hits.iter().map(|h| format!("{}:{:.1e}:{}", h.accession, h.e_value, domain::FRAME_NAMES[h.frame])).collect::<Vec<_>>().join(",")
                        };
                        (d.order.clone(), d.superfamily.clone(), d.strand.clone(), if best.is_finite() { format!("{:.1e}", best) } else { ".".into() }, list)
                    }
                    None => (".".into(), ".".into(), ".".into(), ".".into(), ".".into()),
                };
                // hit positions (genome coordinates, 1-based) and the nesting call
                let (dpos, dnest) = match dc {
                    Some(d) if !d.hits.is_empty() => (
                        d.hits.iter().map(|h| format!("{}:{}:{}-{}", h.accession, domain::FRAME_NAMES[h.frame], rel1(d.origin + h.start), rel1(d.origin + h.end - 1))).collect::<Vec<_>>().join(","),
                        domain::nesting(d).unwrap_or(".").to_string(),
                    ),
                    _ => (".".to_string(), ".".to_string()),
                };
                let (it, imf, igf, igb) = match ist {
                    Some(s) => (s.tiles.to_string(), format!("{:.3}", s.multi_frac), format!("{:.3}", s.gap_frac), s.gap_bp.to_string()),
                    None => (".".into(), ".".into(), ".".into(), ".".into()),
                };
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    name, rel1(l1s), rel1(l2e), l2e - l1s + 1,
                    rel1(l1s), rel1(l1e), l1e - l1s + 1,
                    rel1(l2s), rel1(l2e), l2e - l2s + 1,
                    sim,
                    rel1(ti.left_ltr_start), rel1(ti.left_ltr_end), rel1(ti.right_ltr_start), rel1(ti.right_ltr_end), c.similarity,
                    reextended, score, cols, tr.0, tr.1, tr.2, tr.3, edge,
                    tsd.len, if tsd.len > 0 { tsd.seq.as_str() } else { "." }, tsd.shift_left, tsd.shift_right,
                    motif, features::motif_class(&motif), ppt_n, ppt_off, pbs_off, ns, ti.tsd_len,
                    it, imf, igf, igb,
                    dord, dsf, dstr, dbe, dhits, if ldh.is_empty() { ".".to_string() } else { ldh.join(",") },
                    sc_l, sc_r,
                    dec.keep as u8, dec.why,
                    pbs_trna, pbs_mm, pbs_ed, pbs_strand, tsd.mismatches, dpos, dnest,
                );
            }
        }
    }

    Ok(())
}

// --------------------------------------------------
/// The candidates as GenomeTools LTRharvest's table (final ends when re-extended).
fn print_gt(
    args: &Args,
    st: &Setup,
    sufr: &str,
    seq_starts: &[usize],
    surviving: &[(ltr::LtrCandidate, tsd::TsdMotifInfo)],
    reext: &[Option<reextend::Reextended>],
) {
    println!(
        "# args=-index {} -seed {} -minlenltr {} -maxlenltr {} -mindistltr {} -maxdistltr {} -similar {} -xdrop {}{} -overlaps {}{}",
        sufr,
        args.min_seed, args.ltr_min, args.ltr_max,
        args.dist_min, args.dist_max,
        args.similar as u32, args.xdrop,
        if args.ltrharvest {
            format!(" -mintsd {} -maxtsd {}", args.mintsd.unwrap_or(4), args.maxtsd.unwrap_or(20))
        } else { String::new() },
        format!("{:?}", st.overlaps).to_lowercase(),
        if args.ltrharvest { String::new() } else {
            format!(" -reextend {}{} -reextend-vicinity {} -trim {}/{}",
                args.reextend_ltr_max,
                if args.reextend_first > 0 { format!(" -reextend-first {}", args.reextend_first) } else { String::new() },
                args.reextend_vicinity, args.trim_matches, args.trim_window)
        },
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
    for ((c, ti), re) in surviving.iter().zip(reext.iter()) {
        let seq_off = seq_starts.get(c.seq_num).copied().unwrap_or(0);
        let rel1 = |abs: usize| abs - seq_off + 1;
        let (l1s, l1e, l2s, l2e, sim) = match re {
            Some(r) => (r.left_ltr_start, r.left_ltr_end, r.right_ltr_start, r.right_ltr_end, r.identity),
            None => (ti.left_ltr_start, ti.left_ltr_end, ti.right_ltr_start, ti.right_ltr_end, c.similarity),
        };
        println!(
            "{}  {}  {}  {}  {}  {}  {}  {}  {}  {:.2}  {}",
            rel1(l1s),
            rel1(l2e),
            l2e - l1s + 1,
            rel1(l1s),
            rel1(l1e),
            l1e - l1s + 1,
            rel1(l2s),
            rel1(l2e),
            l2e - l2s + 1,
            sim,
            c.seq_num,
        );
    }
}

// --------------------------------------------------
/// Write tRNA calls as a table: 1-based coordinates within the sequence,
/// start <= end, strand, ARAGORN name, intron length, score, gene sequence.
fn write_trna_calls(
    out: &mut impl std::io::Write,
    names: &[String],
    calls: &[Vec<trna::TrnaCall>],
    ids: Option<&std::collections::HashMap<&str, &str>>,
) -> Result<()> {
    writeln!(out, "seq\tstart\tend\tstrand\tname\tlength\tintron\tscore\tsequence{}", if ids.is_some() { "\tlibrary_id" } else { "" })?;
    for (name, cs) in names.iter().zip(calls) {
        for c in cs {
            write!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}",
                name, c.start, c.end, if c.minus { '-' } else { '+' }, c.name, c.seq.len(), c.intron_len, c.energy, c.seq,
            )?;
            if let Some(ids) = ids {
                write!(out, "\t{}", ids.get(c.seq.as_str()).unwrap_or(&"."))?;
            }
            writeln!(out)?;
        }
    }
    Ok(())
}

// --------------------------------------------------
/// `--trna-out`: the tRNA library (FASTA, one entry per distinct gene
/// sequence, annotated with PBS usage), every call in `<FILE>.tsv`, and the
/// top primers on stderr.
fn write_trna_library(
    path: &str,
    names: &[String],
    calls: &[Vec<trna::TrnaCall>],
    hits: &[(Option<pbs::PbsHit>, Option<usize>)],
    decisions: &[decide::Decision],
) -> Result<()> {
    let flat: Vec<(&str, &trna::TrnaCall)> =
        names.iter().zip(calls).flat_map(|(n, cs)| cs.iter().map(move |c| (n.as_str(), c))).collect();
    let usage: Vec<(Option<&pbs::PbsHit>, bool)> = hits.iter().zip(decisions).map(|((h, _), d)| (h.as_ref(), d.keep)).collect();
    let lib = pbs::build_library(&flat, &usage);
    let mut out = std::io::BufWriter::new(std::fs::File::create(path).with_context(|| format!("creating {}", path))?);
    pbs::write_library_fasta(&mut out, &lib)?;
    std::io::Write::flush(&mut out)?;
    let ids: std::collections::HashMap<&str, &str> = lib.iter().map(|e| (e.seq.as_str(), e.id.as_str())).collect();
    let tsv = format!("{}.tsv", path);
    let mut out = std::io::BufWriter::new(std::fs::File::create(&tsv).with_context(|| format!("creating {}", tsv))?);
    write_trna_calls(&mut out, names, calls, Some(&ids))?;
    std::io::Write::flush(&mut out)?;
    // Top primers by tRNA type: kept elements whose best hit names the type.
    let mut by_type: std::collections::HashMap<&str, [usize; 2]> = std::collections::HashMap::new();
    for (h, keep) in &usage {
        if let (Some(h), true) = (h, keep) {
            for t in &h.trnas {
                by_type.entry(t.as_str()).or_default()[(h.mm12.min(1)) as usize] += 1;
            }
        }
    }
    let mut top: Vec<(&str, [usize; 2])> = by_type.into_iter().collect();
    top.sort_by(|a, b| (b.1[0] + b.1[1]).cmp(&(a.1[0] + a.1[1])).then(a.0.cmp(b.0)));
    eprintln!(
        "[11 output] tRNA library: {} distinct gene sequence(s) from {} call(s) written to {} (calls in {}); top primers of kept elements: {}",
        lib.len(), flat.len(), path, tsv,
        top.iter().take(8).map(|(t, c)| format!("{} {} ({} exact)", t, c[0] + c[1], c[0])).collect::<Vec<_>>().join(", "),
    );
    Ok(())
}

// --------------------------------------------------
/// `--trna-scan`: tRNA genes of every sequence in a FASTA file.
fn trna_scan_fasta(path: &str, args: &Args) -> Result<()> {
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().ok();
    }
    let t0 = std::time::Instant::now();
    let text = std::fs::read(path).with_context(|| format!("reading {}", path))?;
    let mut names: Vec<String> = Vec::new();
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    for line in text.split(|&b| b == b'\n') {
        if let Some(h) = line.strip_prefix(b">") {
            let h = String::from_utf8_lossy(h);
            names.push(h.split_whitespace().next().unwrap_or("").to_string());
            seqs.push(Vec::new());
        } else if let Some(s) = seqs.last_mut() {
            s.extend(line.iter().copied().filter(|b| !b.is_ascii_whitespace()));
        }
    }
    let calls: Vec<Vec<trna::TrnaCall>> = seqs.par_iter().map(|s| trna::scan_sequence(s, args.trna_intron_max)).collect();
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    write_trna_calls(&mut out, &names, &calls, None)?;
    std::io::Write::flush(&mut out)?;
    eprintln!(
        "tRNA scan: {} sequence(s), {} bp, {} gene(s) in {:.1}s",
        seqs.len(), seqs.iter().map(|s| s.len()).sum::<usize>(), calls.iter().map(|c| c.len()).sum::<usize>(), t0.elapsed().as_secs_f64(),
    );
    Ok(())
}

// --------------------------------------------------
fn main() -> Result<()> {
    let args = Args::parse();
    if let Some(fasta) = &args.trna_scan {
        return trna_scan_fasta(fasta, &args);
    }
    let st = setup(&args)?;
    let text_len = read_text_length(args.sufr.as_deref().context("a SUFR index is required")?)?;
    if text_len < u32::MAX as usize {
        run::<u32>(&args, &st)
    } else {
        run::<u64>(&args, &st)
    }
}
