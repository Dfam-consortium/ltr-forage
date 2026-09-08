# ltr-forage

**ltr-forage** is a drop-in replacement for GenomeTools
[LTRharvest](https://github.com/genometools/genometools) that identifies LTR
retrotransposon candidates in genomic sequences.  It implements the same core
search algorithm as LTRharvest, the method of Ellinghaus, Kurtz, and
Willhoeft, while replacing the underlying index with the 
[sufr](https://github.com/TravisWheelerLab/sufr) suffix-array format and 
running both the seed-enumeration and seed-extension phases in parallel.

Its output is not identical to LTRharvest's for two reasons.  First,
LTRharvest uses an unstable sort on extended boundaries, so the two tools
will make different choices when several seeds extend to the same 
outer boundaries. Second, LTRharvest treats N as a wildcard that
matches nothing, not even another N; sufr and ltr-forage compare
bytes, so a run of N in one LTR copy aligns to a run of N in the other.
On hg38 with matching parameters ltr-forage
reports 41,661 candidates against LTRharvest's 41,665: 41,245 agree in every
coordinate, 254 share outer boundaries but differ inside, 134 differ in
position, and 60 appear in only one of the two outputs.

`ltr-forage` holds neither the suffix array nor the LCP array in memory.  It
streams both from the file during phase 1, reading each once, in order, so
peak memory scales with the genome rather than with the index: searching the
3.2 Gb hg38 genome needs about 4.5 GB.

Both compute-heavy phases run in parallel with
[Rayon](https://github.com/rayon-rs/rayon).  Seed extension is easy to
parallelise because each seed is independent.  Seed enumeration is a bottom-up
suffix-tree traversal that LTRharvest runs on one core; ltr-forage cuts the
suffix array at positions whose LCP value falls below the seed length, which no
repeat of at least that length can span, and traverses the resulting pieces
independently.  A whole-genome hg38 search then takes about 45 seconds on 32
cores, against 56 minutes for single-threaded LTRharvest.  Indexing is a
separate one-time step.

Predictions do not depend on how many threads ran.  Where the surviving
candidate would otherwise follow the order in which seeds arrived,
`ltr-forage` breaks the tie on the seed's enumeration index, so every thread
count yields byte-identical output.  LTRharvest's own `-j` mode gives
different results at different job counts: on human chromosome 22 it returns
different boundaries for seven of its 783 predictions at 32 jobs, and takes
nearly four times as long as it does on one.

`--tiered-tsd` ranks candidate target site duplications by the length and
terminal motif combinations that characterise real insertions, following
LTR_retriever, instead of by distance from the X-drop boundary alone.

`--reextend` replaces the X-drop boundaries of each surviving candidate with
a local alignment of the two LTR regions, computed in process with the
RepeatMasker scoring scheme. X-drop with the LTRharvest defaults stops at any
2 bp indel between the two LTR copies, which leaves many real elements with
LTRs short by tens to hundreds of bases; on hg38, re-extension places 1,155
of the 1,929 RepeatMasker-annotated elements it finds within 20 bp of their
annotated boundaries, against 419 of 1,902 without it, and adds about ten
seconds to the run.

`--domain-hmm` searches each candidate's six-frame translation against a
HMMER3 profile file in process, using `libnail`, and reports whether the hits
belong to LTR retrotransposon proteins or to other transposon and host
proteins. On hg38 with LTR_retriever's TEfam profiles, a hit of LTR order
marks 68% of the RepeatMasker-annotated elements and 0.5% of the other
candidates, which is the evidence LTR_retriever gets from HMMER, found here
on more elements.

By default, `ltr-forage` produces output that is format-compatible with
LTRharvest, making it a transparent replacement in existing pipelines.

---

## Background and Credits

The search algorithm implemented here is a Claude (AI) generated port of the
original LTRharvest algorithm described in:

> Ellinghaus D, Kurtz S, Willhoeft U.
> **LTRharvest, an efficient and flexible software for de novo detection of LTR
> retrotransposons.**
> *BMC Bioinformatics.* 2008;**9**:18.
> DOI:[10.1186/1471-2105-9-18](https://doi.org/10.1186/1471-2105-9-18)

LTRharvest is distributed as part of
[GenomeTools](https://github.com/genometools/genometools):

> Gremme G, Steinbiss S, Kurtz S.
> **GenomeTools: A Comprehensive Software Library for Efficient Processing of
> Structured Genome Annotations.**
> *IEEE/ACM Trans Comput Biol Bioinform.* 2013;**10**(3):645–56.
> DOI:[10.1109/TCBB.2013.68](https://doi.org/10.1109/TCBB.2013.68)

The suffix-array index used by ltr-forage is provided by the
[`sufr`](https://github.com/TravisWheelerLab/sufr) project, developed by Ken
Youens-Clark, Jack Roddy, and Travis Wheeler at the University of Arizona.

---

## The LTRharvest Algorithm

LTR retrotransposons are mobile genetic elements bounded by a pair of nearly
identical Long Terminal Repeats (LTRs).  LTRharvest detects candidates by
exploiting the defining property that the two LTRs share a long, nearly-exact
sequence match.  The algorithm proceeds in four main phases.

### Phase 1 — Maximal Repeat Enumeration (Seed Finding)

The genome is indexed as a suffix array (SA) together with an LCP (Longest
Common Prefix) array.  These two arrays together encode the same information as
a suffix tree and enable efficient enumeration of all *maximal repeats*: exact
matches between two genomic positions that cannot be extended one character to
the left without introducing a mismatch (left-maximality).

Seed enumeration is performed by a bottom-up traversal of the implicit suffix
tree, based on the algorithm of Kurtz (1999).  The traversal visits each
internal node of the suffix tree exactly once.  At each node — representing a
set of suffixes sharing a common prefix of length *d* ≥ `--min-seed` — it
tracks the start positions of the child suffixes grouped by the character
immediately preceding the shared prefix (the *left-extension character*).  Any
two positions that differ in their left-extension character constitute a
left-maximal repeat of length *d* and are emitted as a seed pair.

Seeds are then filtered to retain only those satisfying the LTR structural
constraints:
- Distance between the two occurrences (5′-end to 5′-end) within
  [`--dist-min`, `--dist-max`]
- Repeat length ≤ `--ltr-max`
- Both occurrences within the same input sequence

### Phase 2 — X-Drop Gapped Extension

Each seed represents a short exact match that is the core of a putative LTR
pair.  It is extended in both directions using *X-drop gapped alignment*, a
heuristic that continues the alignment as long as the accumulated score remains
within a threshold *X* (`--xdrop`) of the running maximum.  The scoring scheme
(match award `--mat`, mismatch penalty `--mis`, indel penalties `--ins`/`--del`)
is configurable; the defaults match LTRharvest's defaults.

Extension is run twice per seed — leftward from the 5′ ends and rightward from
the 3′ ends — to determine the full LTR boundaries.

Extension yields provisional LTR boundaries.  ltr-forage applies the length,
distance and similarity filters after the TSD search described below, since
that search can move the outer boundaries.  GenomeTools uses the same order.

### Phase 3 — Post-Extension Refinement

#### Duplicate Removal

Multiple seeds can yield candidates with identical post-TSD outer LTR
boundaries.  These duplicates are removed by sorting on (sequence, left_ltr_start,
right_ltr_end) and collapsing identical records.

#### Target Site Duplication (TSD) Search

When TSD search is enabled (`--mintsd` ≥ 2), ltr-forage searches for *target
site duplications* — short exact direct repeats flanking the element at the
outermost LTR boundaries.  TSDs arise because the integrase creates a short
staggered cut in the host DNA during insertion; when the single-stranded gap is
repaired, the cut site is duplicated on both sides of the inserted element.

The search is performed by building a local suffix array of a window of width
`--vicinity` bp around the 5′ end of the left LTR, then querying each position
in a corresponding window around the 3′ end of the right LTR for an exact match
of length in [`--mintsd`, `--maxtsd`].  This mirrors GenomeTools'
`searchforbestTSDandormotifatborders` and `gt_sarrquerysubstringmatch`
procedures, including their tie-breaking order (right window position ascending,
then SA rank ascending), which is preserved for GT-compatible output.

For each candidate TSD the *total deviation* — the sum of the absolute
differences between the TSD-adjusted LTR boundaries and the original
post-extension boundaries — is computed.  By default (GT-compatible mode) the
TSD with minimum total deviation is selected.  The optional `--tiered-tsd` flag
enables a biologically motivated ranking that prefers TSDs of exactly 5 bp with
a TGCA motif (tier 1), then 5 bp with a known alternative motif (tier 2), then
shorter TSDs with TGCA (tier 3), and so on, with deviation as a tiebreaker
within each tier (inspired by the LTR_retriever project).

When TSD search is enabled, ltr-forage discards any candidate for which it
found no TSD, as GenomeTools does.  It then checks the adjusted boundaries
against the LTR length and 5′-to-5′ distance bounds, computes the similarity of
the adjusted LTRs with Ukkonen's greedy edit-distance algorithm (O(*n*·*d*)
time, O(*d*) space), and discards candidates below `--similar` percent
identity.

If `--motif` is specified, only TSD candidates where the adjusted LTR junctions
match the TG...CA palindromic motif (configurable via `--motifmis`) are kept.
The inner LTR borders are also checked for the complementary CA/TG motif.

#### Overlap Removal

Overlapping predictions on the same sequence are resolved according to
`--overlaps`:
- `best` (default) — within each cluster of mutually overlapping elements,
  keep the one with the highest similarity
- `no` — discard all elements in any overlapping cluster
- `all` — report all elements regardless of overlap

### Phase 4 — Re-extension (optional, `--reextend`)

The X-drop extension in Phase 2 uses LTRharvest's scores (match +2,
mismatch −2, gap −3, threshold 5), so any indel of 2 bp or more between the
two LTR copies ends the extension there, and the bases beyond it are left out
of both LTRs even when they are otherwise identical. Re-extension recovers
them. For each candidate that survives overlap removal:

1. A window is cut around each LTR, reaching `--reextend-ltr-max` (default
   2000) plus `--reextend-vicinity` (default 100) bases beyond the current
   LTR end on the outer side, and the same distance on the inner side up to
   half of the internal region, so the two windows never overlap.
2. The left window is aligned against the right window with `rmblast-lib`,
   using RepeatMasker's nucleotide comparison matrix (embedded in the binary),
   gap open 20 and extend 5, word size 14, and complexity-adjusted scoring.
3. The best-scoring plus-strand HSP that overlaps both X-drop LTRs becomes
   the new LTR pair. The overlap requirement matters: the windows can reach a
   neighbouring copy of the same family, and its alignment may score higher
   than the pair the seed came from.
4. Each end of the alignment is trimmed until its terminal `--trim-window`
   columns (default 6) are gap-free and contain at least `--trim-matches`
   (default 6) matches. This removes the few bases of chance matches a local
   alignment carries past the true end, without consulting the TSD, so the
   TSD found afterwards remains evidence rather than a consequence of the
   trim.
5. The reported similarity becomes the identity of the trimmed alignment,
   with gaps counted as differences.

A candidate whose alignment produces no qualifying HSP keeps its X-drop
boundaries. Larger windows are not better in genomes with nested elements:
in rice, 5 kb windows lost real elements to alignments that ran into
neighbours, and 2 kb windows gained them.

### Phase 5 — Protein-domain search (optional, `--domain-hmm`)

Most of the candidates a whole-genome run reports are not LTR elements: on
hg38, 95% are pairs of SINEs, simple repeats or other similar sequence a few
kb apart. Of the signals we measured on hg38, rice and Arabidopsis, a
retrotransposon protein domain in the internal region separates the real
elements best. `--domain-hmm FILE` searches each candidate's six-frame
translation against a HMMER3 file of profile HMMs in process, using
`libnail`'s implementation of HMMER's profile scoring. LTR_retriever's
`TEfam.hmm` (102 Pfam profiles) is the intended input; any HMMER3 amino-acid
profile file works.

ltr-forage translates each element plus `--domain-flank` bases (default 50)
on each side in six frames, with stop codons as X. For every frame and
profile it scores each diagonal that shares a 3-mer with the profile's
consensus sequence as an ungapped alignment against the profile's match
emissions; when the best diagonal reaches 18 bits, `libnail` runs the
Forward algorithm over a band of rows around it, applies its composition-bias
correction, and converts the score to an E-value with the profile's own
statistics, taking six frames per candidate as the database size, as
hmmsearch would. Hits with an E-value at or below `--domain-evalue`
(default 0.05) count. The seeding replaces the MMseqs2 stage nail's own
pipeline uses; on chr22 it keeps 97% of the pairs a full Forward pass would
report at that threshold while passing 7% of all pairs. On the same chr22
candidates ltr-forage finds 39 of the 43 domain hits hmmsearch reports, with
E-values within a few fold, and adds 16 more, mostly near the threshold. On
hg38 the search scores 25.5 million frame-profile pairs, seeds 1.7 million
of them, and takes 27 minutes on 64 threads, most of it in `libnail`'s
Forward pass. We did not time hmmsearch on the same input; its SIMD filters
make it likely faster, so a stronger prefilter is the next thing to try if
the run time matters.

The hits are classified as LTR_retriever's `annotate_TE.pl` does: profiles
belong to an LTR list or a not-LTR list (DNA transposons, LINEs, helitrons,
host genes); a candidate whose hits are all on the LTR list is `LTR`, all on
the other list `notLTR`, both `mixture`, and hits on neither list `NA`.
Superfamily is `Gypsy` when Copia-list hits are fewer than 0.3 of the
Gypsy-list hits, `Copia` when Gypsy-list hits are fewer than 0.9 of the
Copia-list hits, and strand `+` or `-` when hits on the other strand are
fewer than 0.4 of those on it. The HMMER3
reader is adapted from nail's, under its BSD-3-Clause licence.

### Phase 6 — Decision and library (`--keep`, `--library`)

Every candidate gets one decision, recorded in the table's `keep` and
`evidence` columns and used to select what `--library FILE` writes. The
rule has one knob, `--keep LEVEL`:

| Level | Kept when |
|---|---|
| `domain` | an LTR-order domain hit |
| `standard` (default) | an LTR-order domain hit, or TG..CA at the ends with a TSD of 5 bp or more |
| `loose` | an LTR-order domain hit, a TSD of 5 bp or more, or TG..CA with a 4 bp TSD |
| `none` | always |

At every level ltr-forage rejects a candidate whose protein hits are of
not-LTR or mixed order; those are LINE, DNA-transposon and host-gene pairs,
and on hg38 1.3% of them are annotated LTR elements. `--internal-gap-max FRAC` adds a
second rejection for candidates whose internal region has a single-copy
stretch longer than that fraction of it (needs `--internal-stats`); it is
off by default because it costs recall in genomes with low-copy families
(Arabidopsis).

ltr-forage also rejects satellite arrays and Alu dimers with an ERV
inserted between two monomers; those pass the evidence rule as
domain-bearing elements with a satellite or Alu "LTR". It aligns each LTR
against itself with
`rmblast-lib` and records the fraction of it covered by plus-strand
alignments off the main diagonal (`lltr_self_cov`, `rltr_self_cov`); a
satellite monomer, an Alu dimer or a microsatellite covers itself
completely. `--tandem-max` (default 0.99) rejects a candidate when either
LTR reaches that coverage. The threshold has to sit near 1: real LTRs with
tandem sub-repeats, the LTR12 family in human among them, reach 0.5 to 0.98,
and a threshold of 0.5 would have removed 71 annotated elements on hg38. At
0.99 the check removes 17 kept candidates on hg38 (10 satellite cases, 5
Alu-dimer pairs, 2 on unannotated contigs), 52 on rice and none on
Arabidopsis, and no annotated element on any of the three. The alignment
adds 40 s to an hg38 run.

Measured against RepeatMasker annotation (strict or loose harness), the
levels keep:

| Genome | Candidates | Annotated | `domain` | `standard` | `loose` |
|---|---|---|---|---|---|
| hg38 | 41,661 | 2,267 | 68.5% of annotated, 88.7% precise | 68.7%, 86.7% | 71.6%, 70.5% |
| Arabidopsis | 1,789 | 415 | 78.1%, 62.2% | 81.4%, 62.8% | 82.7%, 61.5% |
| rice | 7,092 | 3,216 | 84.3%, 69.5% | 86.4%, 69.9% | 87.1%, 69.3% |

The plant precisions are floors set by the annotation: when we inspected
the kept candidates the harness calls false there (analysis document,
sections 9 and 10), most were elements whose two LTRs the library names
differently, or nested copies. On rice the motif-and-TSD route alone is 98%
precise and adds 2% of recall over `domain`. There is no separate
tandem-repeat finder: of the 1,797 candidates `standard` keeps on hg38
before the self-alignment check, 11 have an LTR that is mostly satellite,
against 1,091 of all 41,661 candidates, and the check above removes 10 of
the 11.

`--library FILE` writes the kept elements as three FASTA records each, left
LTR, internal region and right LTR, named the way LTR_retriever names its
redundant library and in genome orientation:

```
>chr22:7291483..7292450_LTR#LTR/unknown domain,tsd6,tgca
>chr22:7292451..7299693_INT#LTR/unknown domain,tsd6,tgca
>chr22:7299694..7300657_LTR#LTR/unknown domain,tsd6,tgca
```

The superfamily comes from the domain hits (`Gypsy`, `Copia` or `unknown`)
and the description repeats the evidence. RepeatModeler's LTRPipeline reads
this file in place of LTR_retriever's output.

---

## Differences from GenomeTools LTRharvest

| Feature | LTRharvest (GenomeTools) | ltr-forage |
|---|---|---|
| Index format | GenomeTools ESA (`suffixerator`) | sufr suffix array (`.sufr` file) |
| Seed enumeration | Bottom-up ESA traversal (C), single-threaded | Same algorithm, streaming the sufr SA/LCP from disk, parallel over LCP-delimited pieces (`-j`) |
| Seed extension | Single-threaded | Parallelised with Rayon (`-j`) |
| TSD selection | Minimum total deviation | Minimum deviation (default); tiered biological ranking (`--tiered-tsd`) |
| Output | Space-separated, GT-compatible | GT-compatible by default (`--format gt`) |

The seed enumeration (`maxpairs.rs`) and X-drop extension (`xdrop.rs`) are
faithful Rust ports of the corresponding GenomeTools routines
(`esa-maxpairs.c`, `esa-bottomup-maxpairs.inc`,
`gt_evalxdroparbitscoresextend`).  The TSD search (`tsd.rs`) mirrors
`searchforbestTSDandormotifatborders` and `gt_sarrquerysubstringmatch`,
including the tie-breaking emission order required to reproduce GT-compatible
boundary choices.

The first structural departure is the use of the `sufr` index rather than a
GenomeTools ESA.  Because `sufr` stores a combined suffix array over all input
sequences with separators excluded from the active suffix count, the bottom-up
traversal reads directly from `sufr`'s pre-built SA and LCP arrays without the
intermediate ESA query layer present in GenomeTools.  An LCP shift
(look-ahead vs. look-back convention) is applied during traversal to match the
GenomeTools array layout.

The suffix array and LCP array are never loaded into memory.  Phase 1 consumes
each of them once, in order, so `ltr-forage` streams both from the `.sufr`
file through a small buffer (`stream.rs`) and holds only the text in memory.
Peak memory is therefore close to the genome size rather than to the index
size.

Phase 1 also runs in parallel.  A suffix array position whose LCP value is
below `--min-seed` cannot lie inside any repeat of seed length or more, so
cutting the array at such positions gives pieces that the traversal can
process independently.  `ltr-forage` cuts the array into about four pieces
per thread (`stream::split_points`), traverses each piece on its own worker
with its own file readers, and concatenates the seeds in piece order.  That is
the order a single sequential traversal would emit them, so the seed-index
tie-break below and the output are unchanged by the thread count.

The parallelisation of Phase 2 (X-drop extension + TSD annotation + similarity
computation) is embarrassingly parallel: each seed is independent.  A
deterministic seed-index tiebreaker is attached before the parallel scatter and
used during the subsequent sort so that the surviving candidate from any set of
seeds with identical post-TSD boundaries is always the one with the lowest
SA-enumeration index — the same choice a single-threaded stable sort on
insertion order would make.

---

## Installation

### Prerequisites

- Rust toolchain ≥ 1.80 ([rustup](https://rustup.rs) recommended)
- The `sufr` tool for building `.sufr` index files
  ([TravisWheelerLab/sufr](https://github.com/TravisWheelerLab/sufr))

### Build from source

```bash
git clone <repository-url>
cd ltr-forage
cargo build --release
# Binary: target/release/ltr-forage
```

### Create a sufr index

Before running `ltr-forage`, build a `.sufr` index of your genome with the
`sufr` tool.  The `--dna` flag is required to restrict the alphabet to ACGT:

```bash
sufr create --dna genome.fa -o genome.sufr
```

`sufr create --dna` stores the sequence upper-cased, so a soft-masked FASTA
(lowercase runs, as NCBI distributes them) needs no preprocessing; the index
and everything `ltr-forage` reads from it are case-free. Query strings given
to `sufr count` are not normalised, so pass them in upper case.

---

## Usage

```
ltr-forage [OPTIONS] <SUFR>

Arguments:
  <SUFR>  Path to a .sufr index file (created with `sufr create --dna`)

Seed and extension options:
  -m, --min-seed <BP>       Minimum seed length (min exact match) [default: 30]
  -l, --ltr-min <BP>        Minimum LTR length [default: 100]
  -L, --ltr-max <BP>        Maximum LTR length [default: 1000]
  -d, --dist-min <BP>       Minimum 5'-to-5' distance between LTRs [default: 1000]
  -D, --dist-max <BP>       Maximum 5'-to-5' distance between LTRs [default: 15000]
  -s, --similar <PCT>       Minimum percent similarity between LTR copies [default: 85]
  -x, --xdrop <SCORE>       X-drop threshold for extension alignment [default: 5]
      --mat <MAT>           Match score (positive) [default: 2]
      --mis <MIS>           Mismatch score (negative) [default: -2]
      --ins <INS>           Insertion score (negative) [default: -3]
      --del <DEL>           Deletion score (negative) [default: -3]

TSD and motif options:
      --mintsd <BP>         Minimum TSD length; 0 or 1 disables TSD search [default: 4]
      --maxtsd <BP>         Maximum TSD length [default: 20]
      --vicinity <BP>       Search window around LTR boundaries for TSD/motif [default: 60]
      --motif               Enable TG...CA palindromic motif search
      --motifmis <N>        Allowed mismatches in the motif [default: 0]
      --tiered-tsd          Use tiered TSD selection instead of GT-compatible deviation-only

Re-extension options:
      --reextend            Re-extend each LTR pair with a local alignment after overlap removal
      --reextend-ltr-max <BP>   Longest LTR the alignment windows allow for [default: 2000]
      --reextend-vicinity <BP>  Extra window beyond that on the outer side of each LTR [default: 100]
      --trim-window <N>     Terminal alignment columns that must be clean for an end to stand [default: 6]
      --trim-matches <N>    Matches required inside that terminal window [default: 6]
      --internal-stats      Count 20-mers of each internal region against the index and
                            report coverage and the longest single-copy stretch (tsv output)

Domain search options:
      --domain-hmm <FILE>   HMMER3 file of protein profiles; search each candidate's six-frame
                            translation against it and report order, superfamily and strand
      --domain-evalue <E>   E-value threshold for a domain hit [default: 0.05]
      --domain-flank <BP>   Bases added on each side of the element before translation [default: 50]

Decision and library options:
      --keep <LEVEL>        Evidence a candidate needs: domain | standard | loose | none [default: standard]
      --internal-gap-max <FRAC>  Also reject candidates whose internal region has a single-copy
                            stretch longer than this fraction of it (needs --internal-stats)
      --library <FILE>      Write the kept elements as a FASTA library (left LTR, internal
                            region, right LTR per element)
      --tandem-max <FRAC>   Reject a candidate when either LTR aligns to itself off the diagonal
                            over at least this fraction of its length [default: 0.99]

Output and performance options:
      --overlaps <STRATEGY> Overlap removal strategy: all | best | no [default: best]
      --format <FORMAT>     Output format [default: gt] [possible values: gt, tsv]
  -j, --threads <N>         Threads for parallel extension (0 = all CPUs) [default: 0]
  -h, --help                Print help
```

### Example — matching GenomeTools LTRharvest defaults

```bash
sufr create --dna hg38.fa -o hg38.sufr

ltr-forage hg38.sufr \
    --min-seed 30 \
    --ltr-min 100 --ltr-max 1000 \
    --dist-min 1000 --dist-max 15000 \
    --similar 85 \
    --xdrop 5 \
    --mat 2 --mis -2 --ins -3 --del -3 \
    --mintsd 4 --maxtsd 20 \
    --vicinity 60 \
    --overlaps best \
    > hg38_ltrs.out
```

---

## Output Format

The default output (`--format gt`) is compatible with GenomeTools LTRharvest.
The file begins with a comment block describing the parameters and column
layout, followed by one space-separated data line per prediction.

```
# args=-index hg38.sufr -seed 30 -minlenltr 100 -maxlenltr 1000 ...
# predictions are reported in the following way
# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr
# where:
# s = starting position
# e = ending position
# l = length
# ret = LTR-retrotransposon
# lLTR = left LTR
# rLTR = right LTR
# sim = similarity
# seq-nr = sequence number
135253  138333  3081  135253  135376  124  138219  138333  115  87.10  0
```

All positions are **1-based and inclusive**.  The `seq-nr` column is a
**0-based** index into the sequences in the `.sufr` index.

| Column | Field | Description |
|---|---|---|
| 1 | s(ret) | Element start (= left LTR start) |
| 2 | e(ret) | Element end (= right LTR end) |
| 3 | l(ret) | Element length |
| 4 | s(lLTR) | Left LTR start |
| 5 | e(lLTR) | Left LTR end |
| 6 | l(lLTR) | Left LTR length |
| 7 | s(rLTR) | Right LTR start |
| 8 | e(rLTR) | Right LTR end |
| 9 | l(rLTR) | Right LTR length |
| 10 | sim(LTRs) | Percent similarity between LTR copies |
| 11 | seq-nr | Sequence number (0-based) |

---

### Per-candidate table (`--format tsv`)

`--format tsv` writes a header line and one tab-separated row per candidate,
with every value the later stages decide on, so a decision can be inspected
rather than inferred. Coordinates are 1-based and inclusive, relative to the
sequence named in the first column.

| Column | Description |
|---|---|
| `seq` | Sequence name from the index |
| `elem_start`, `elem_end`, `elem_len` | Element span |
| `lltr_start`, `lltr_end`, `lltr_len`, `rltr_start`, `rltr_end`, `rltr_len` | LTR spans after re-extension (or the X-drop spans when `--reextend` is off or no alignment qualified) |
| `identity` | Identity of the trimmed alignment in percent, gaps counted as differences; the X-drop similarity when not re-extended |
| `xdrop_lltr_start` … `xdrop_rltr_end`, `xdrop_similarity` | The boundaries and similarity the X-drop stage produced |
| `reextended` | 1 when the row's boundaries come from the alignment |
| `aln_score`, `aln_columns` | Raw score and aligned columns of the trimmed alignment |
| `trim_l5`, `trim_r5`, `trim_l3`, `trim_r3` | Bases removed by the end trim at the 5′ end of the left and right LTR and at their 3′ ends |
| `window_edge` | 1 when the alignment reached within 5 bp of a window edge that was not a sequence end: the LTR may be longer than the window, or the pair sits in a larger duplication |
| `tsd_len`, `tsd_seq`, `tsd_shift_l`, `tsd_shift_r` | Longest exact direct repeat (4 to 20 bp) flanking the element, allowing 0 or 1 bp of play at each end; `.` and 0 when none |
| `motif` | First two bases of the left LTR and last two of the right LTR (TGCA for the canonical termini) |
| `ppt_purines15`, `ppt_offset` | Purine count in the most purine-rich 15 bp window within 30 bp before the right LTR, and that window's start relative to the right LTR |
| `pbs_offset` | Offset of the first TGG within 20 bp after the left LTR, or −1 |
| `n_count` | Ns inside the element |
| `xdrop_tsd_len` | TSD length the X-drop stage's own search reported |
| `int_tiles`, `int_multi_frac`, `int_gap_frac`, `int_gap_bp` | With `--internal-stats`: number of 20-mer tiles taken every 10 bp across the internal region (tiles containing N skipped), the fraction seen at least three times in the genome on either strand, and the longest run of single-copy tiles as a fraction of the tiles and in bases; `.` otherwise |
| `dom_order`, `dom_superfamily`, `dom_strand` | With `--domain-hmm`: `LTR`, `notLTR`, `mixture`, `NA` or `none`; `Gypsy`, `Copia` or `unknown`; `+`, `-` or `?`, as described under Phase 5; `.` otherwise |
| `dom_best_e`, `dom_hits` | Best E-value, and every hit as `accession:E-value:frame` (frames `aa1`..`aa3` forward, `rev_aa1`..`rev_aa3` on the reverse complement), comma separated; `.` when none |
| `lltr_self_cov`, `rltr_self_cov` | Fraction of each LTR covered by plus-strand alignments of the LTR against itself off the main diagonal; 1.00 for a satellite monomer or an Alu dimer, 0 for unique sequence |
| `keep`, `evidence` | 1 when the candidate passes `--keep`, and what decided it: `domain`, `tsdN`, `tgca` (comma separated), `veto:notLTR`, `veto:mixture`, `veto:tandemF`, `gapF` when `--internal-gap-max` rejected it, or `none` |

The TSD and motif columns are computed at the final boundaries with at most
1 bp of play, so they keep their discriminating power: on hg38 a 5 bp TSD
appears on 17% of RepeatMasker-annotated elements (41% of those over 95%
identical) and on 1.3% of the other candidates. A wider search would find a
short repeat almost everywhere.

The internal-region columns describe whether the sequence between the LTRs
looks repeat-derived end to end. An LTR element's internal region rarely
contains a long single-copy stretch; the sequence between two SINE copies a
few kb apart usually does. On hg38 a longest single-copy stretch of at most
10% of the region holds for 96% of annotated elements and 32% of the other
candidates; in Arabidopsis, whose LTR families are low copy, the same cut
keeps only half the real elements, so the statistic is reported here and any
filter on it is a separate choice. The count loads the suffix array (4 bytes
per base for genomes under 4 Gb) and the tile table: on hg38, 25.5 million
tiles, 140 s and 45 GB at peak; on rice, 30 s.

## Comparing with GenomeTools LTRharvest on hg38

The numbers quoted in this README come from one comparison on hg38 (455
sequences, 3.2 Gb, all uppercase).  Everything needed to repeat it is in
`eval/`: three Python 3 scripts with no dependencies beyond the standard
library, and `eval/run_hg38_comparison.sh`, which runs the steps below in
order and skips any step whose output already exists.  The steps need
GenomeTools 1.6.2 (`gt`), `sufr` 0.7.12, `samtools` for the FASTA index,
`/usr/bin/time` for the timings, and a RepeatMasker `.out` annotation of the
same assembly for the last two steps.

### Indexes

```bash
samtools faidx hg38.fa
gt suffixerator -dna -pl -tis -suf -lcp -sds -des -db hg38.fa -indexname hg38_esa
sufr create --dna hg38.fa -o hg38.sufr
```

Both indexes are built from the same FASTA, so the 0-based sequence numbers in
the two outputs refer to the same sequences, and the `.fai` maps them to
names.  `gt suffixerator` is single-threaded; `sufr create` splits the work
into 16 partitions by default.  Do not pass `--ignore-softmask` to
`sufr create` on a soft-masked assembly: it replaces lowercase bases with N
instead of uppercasing them, which removes the repeats you are looking for.

### Searches

```bash
gt ltrharvest -index hg38_esa -seed 30 -minlenltr 100 -maxlenltr 1000 \
    -mindistltr 1000 -maxdistltr 15000 -similar 85 -xdrop 5 \
    -mat 2 -mis -2 -ins -3 -del -3 -mintsd 4 -maxtsd 20 -vic 60 \
    -overlaps best > hg38_gt.out

ltr-forage hg38.sufr --min-seed 30 --ltr-min 100 --ltr-max 1000 \
    --dist-min 1000 --dist-max 15000 --similar 85 --xdrop 5 \
    --mat 2 --mis -2 --ins -3 --del -3 --mintsd 4 --maxtsd 20 --vicinity 60 \
    --overlaps best -j 32 > hg38_ltr-forage.out
```

These are LTRharvest's default parameters written out in full.  Wrap each
command in `/usr/bin/time -v` for wall time and peak memory; `ltr-forage`
also prints the elapsed time at the end of each phase on stderr.  `gt`
accepts a global `-j N` flag, but on chr22 that ran slower than
single-threaded at 16 and 32 jobs and changed seven of 783 predictions, so
LTRharvest was timed on one thread.

Measured on two 16-core Xeon Silver 4314 sockets with a warm page cache:

| | LTRharvest | ltr-forage |
|---|---|---|
| Index build | 2 h 08 min, 25 GB | 18 min, 50 GB |
| Search, 1 thread | 56 min, 1.35 GB | 19 min 50 s, 3.9 GB |
| Search, 32 threads | | 45 s, 4.5 GB |

### Comparing the predictions

```bash
eval/compare_predictions.py hg38_gt.out hg38_ltr-forage.out --fasta hg38.fa --out-prefix hg38cmp
```

The script pairs predictions within each sequence by best reciprocal overlap
and classifies each pair by how far its coordinates differ.  Unmatched
predictions go to `hg38cmp.a_only.tsv` and `hg38cmp.b_only.tsv`.

| Category | Count |
|---|---|
| LTRharvest predictions | 41,665 |
| ltr-forage predictions | 41,661 |
| Identical in all six coordinates and similarity | 41,245 |
| Same outer boundaries, different inner boundaries | 254 |
| Differ by 1-20 bp / 21-200 bp / over 200 bp | 2 / 61 / 71 |
| Only in LTRharvest / only in ltr-forage | 32 / 28 |

The "same outer boundaries" class is the duplicate tie-break described in the
introduction.  Eight of the 28 predictions found only by ltr-forage contain N
inside an LTR, which is the N-handling difference.

### Scoring against RepeatMasker

Neither tool's output is ground truth, so the scripts score both against the
intact elements in a RepeatMasker annotation of hg38: an LTR-class annotation
followed on the same strand by an internal region and then the nearest
annotation of the same family, with both LTRs 100-1000 bp long and their 5'
ends 1000-15000 bp apart, the window the tools were run with.  The annotation
used here has 7,178 such elements.  A prediction hits one when each of its LTRs overlaps the
corresponding annotated LTR by at least half the length of the shorter.

```bash
eval/rm_truth_recall.py hg38.fa.out hg38.fa.fai hg38_gt.out hg38_ltr-forage.out \
    hg38cmp.a_only.tsv hg38cmp.b_only.tsv
```

Recall by the larger of the two LTRs' divergence from consensus:

| Max LTR divergence | Elements | LTRharvest | ltr-forage |
|---|---|---|---|
| 0-5% | 586 | 90.8% | 90.8% |
| 5-10% | 1,320 | 63.6% | 63.6% |
| 10-15% | 1,794 | 21.2% | 21.2% |
| 15-20% | 1,301 | 7.5% | 7.5% |
| 20% and up | 2,177 | 2.3% | 2.3% |

Both tools hit the same 1,902 elements, and none of the 60 predictions unique
to one tool hits any of them.  The differences between the two outputs fall
entirely on predictions that do not correspond to intact annotated elements.

```bash
eval/rm_classify.py hg38.fa.out hg38.fa.fai hg38_gt.out
```

This labels each prediction by the repeat class under its two LTRs.  For
either tool, 62% of hg38 predictions are pairs of SINE elements, 14% pairs of
LINE fragments, 8.6% LTR class (2,777 of LTRharvest's with the same family on
both sides), and 11% mostly unannotated.  That is a property of the LTRharvest search at these
parameters on a SINE-rich genome rather than of either implementation, and it
is why a RepeatMasker-based check is worth running on any change to the
search.

---

## References

1. Ellinghaus D, Kurtz S, Willhoeft U. LTRharvest, an efficient and flexible
   software for de novo detection of LTR retrotransposons. *BMC Bioinformatics.*
   2008;**9**:18.
   DOI:[10.1186/1471-2105-9-18](https://doi.org/10.1186/1471-2105-9-18)

2. Gremme G, Steinbiss S, Kurtz S. GenomeTools: A Comprehensive Software
   Library for Efficient Processing of Structured Genome Annotations.
   *IEEE/ACM Trans Comput Biol Bioinform.* 2013;**10**(3):645–56.
   DOI:[10.1109/TCBB.2013.68](https://doi.org/10.1109/TCBB.2013.68)

3. Youens-Clark K, Roddy J, Wheeler T. sufr: Parallel Construction of Suffix
   Arrays in Rust.
   [github.com/TravisWheelerLab/sufr](https://github.com/TravisWheelerLab/sufr)

4. Kurtz S. Reducing the space requirement of suffix trees.
   *Software — Practice & Experience.* 1999;**29**(13):1149–1171.

---

## License

This software is released to the public domain under the
[CC0 1.0 Universal](LICENSE) public domain dedication.

The LTRharvest algorithm is implemented based on the published description and
the open-source GenomeTools reference implementation (GPL v2+).  The `sufr`
library is used under its BSD-3-Clause license.
