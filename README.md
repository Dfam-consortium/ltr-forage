# ltr-forage

**ltr-forage** finds LTR retrotransposons in genome assemblies.  It starts
from the search algorithm of GenomeTools
[LTRharvest](https://github.com/genometools/genometools) (Ellinghaus, Kurtz
and Willhoeft), reimplemented on the
[sufr](https://github.com/TravisWheelerLab/sufr) suffix-array index and run
in parallel, and adds the steps LTRharvest leaves to later tools: a mask for
tandem arrays, LTR boundaries from a local alignment of the two copies, and
evidence from each element itself (protein domains, target site duplication,
TG..CA ends, a primer binding site matched to the genome's own tRNA genes)
behind a keep-or-reject decision for every candidate.  It writes a table of
every candidate with that evidence and a FASTA library of the kept elements.

`--ltrharvest` runs LTRharvest's search alone and writes LTRharvest's output
format, as a drop-in replacement in existing pipelines.  The output is not
quite identical: LTRharvest uses an unstable sort on extended boundaries, so
the two tools make different choices when several seeds extend to the same
outer boundaries.  On hg38 with LTRharvest's default parameters both report
41,665 candidates: 41,264 agree in every coordinate, 254 share outer
boundaries but differ inside, 133 differ in position, and 14 on each side
appear in only one output (see
[Comparing with GenomeTools LTRharvest on hg38](#comparing-with-genometools-ltrharvest-on-hg38)).

`ltr-forage` holds neither the suffix array nor the LCP array in memory.  It
streams both from the file during the seed search, reading each once, in
order, so peak memory scales with the genome rather than with the index:
searching the 3.2 Gb hg38 genome needs about 4.5 GB.

The seed search and the per-candidate steps run in parallel with
[Rayon](https://github.com/rayon-rs/rayon).  Seed extension is easy to
parallelise because each seed is independent.  Seed enumeration is a bottom-up
suffix-tree traversal that LTRharvest runs on one core; ltr-forage cuts the
suffix array at positions whose LCP value falls below the seed length, which no
repeat of at least that length can span, and traverses the resulting pieces
independently.  LTRharvest's search of hg38 takes about two minutes on 16
cores, against 56 minutes for single-threaded LTRharvest.  Indexing is a
separate one-time step.

Results do not depend on how many threads ran.  Where the surviving
candidate would otherwise follow the order in which seeds arrived,
`ltr-forage` breaks the tie on the seed's enumeration index, so every thread
count yields byte-identical output.  LTRharvest's own `-j` mode gives
different results at different job counts: on human chromosome 22 it returns
different boundaries for seven of its 783 predictions at 32 jobs, and takes
nearly four times as long as it does on one.

`pipeline/ltr_library.py` builds on the kept elements to produce an LTR
family library from a genome assembly: one LTR and one internal model per
family, cleaned of nested insertions and segmental duplications, in place of
LTRharvest followed by LTR_retriever (see
[Building a family library](#building-a-family-library)).

![Pipeline](docs/pipeline_figure.svg)

The numbered boxes are the steps `ltr-forage` (1-11) and the library driver
(12-24) name in their logs.

---

## Quick start

See [Prerequisites](#prerequisites) for sufr, nail and MMseqs2.

### 1. Index the genome (once)

```bash
sufr create --dna genome.fa -o genome.sufr
```

Do not pass `--ignore-softmask` on a soft-masked assembly (see
[Create a sufr index](#create-a-sufr-index)).

### 2. Detect LTR retrotransposon copies

```bash
ltr-forage genome.sufr --library elements.fa --fastan-out tandem.bed -j 24 \
    > cand.tsv 2> ltr-forage.log
```

The defaults are the settings the library driver uses: tandem arrays
masked, LTR pairs re-aligned in a 7 kb window, calls whose final ends agree
within 20 bp merged, protein domains searched with LTR_retriever's
`TEfam.hmm` (shipped in `data/`), overlapping calls resolved by evidence,
and the per-candidate table on stdout.  The domain search writes about 2 GB
of scratch files for a primate genome to `$TMPDIR` or `/tmp`;
`--temp-dir DIR` puts them elsewhere.  `ltr-forage` checks its programs,
the profile file and the output directories before it loads the index, and
reports every problem at once.

| File | Contents |
|---|---|
| `cand.tsv` | every candidate, kept or not, with the evidence behind the call: boundaries, LTR identity, TSD, TG..CA, PBS and its tRNA, protein domains, LTR self-repeat; `keep` is 1 for a kept element and `evidence` says which rule kept or rejected it (see [Per-candidate table](#per-candidate-table)) |
| `elements.fa` | each kept element's 5′ LTR, internal region and 3′ LTR |
| `tandem.bed` | the masked tandem arrays (name, start, end, period, 1000 × identity) |
| `ltr-forage.log` | one line per step, numbered as in the figure and in [How ltr-forage finds elements](#how-ltr-forage-finds-elements) |

On rice (MSU7, 373 Mb) with 16 threads the log reads:

```
[1 index] 12 sequence(s), 373128866 suffixes, text length 373245531; text loaded in 0.7s
[1 index] parameters: min_seed=30 ltr=[100,1000] dist=[1000,15000] similar=85% xdrop=5
[2 tandem mask] FasTAN port: 115154 array(s), 21.6 Mb before merging (29.3s)
[3 seed pairs] 61599 seed pair(s) (64 piece(s); 37.9s elapsed)
[4 X-drop extension] 40218 candidate(s) at the similarity threshold (38.1s elapsed)
[5 duplicates] 17569 candidate(s) after removing identical ends
[6 re-extension] re-extended 17565 of 17569 candidate(s), 16583 moved, 983 at a window edge (2m07s elapsed)
[7 near-duplicate merge] 11246 of 17569 candidate(s) left after merging final ends within 20 bp
[8a self-alignment] 146 candidate(s) with an LTR at least 99% self-repeat (3m52s elapsed)
[8b protein domains] skipping 146 candidate(s) vetoed by LTR self-repeat (--domain-search-vetoed searches them)
[8b protein domains] nail search (2m05s): 6356 candidate(s) with a domain hit, 5597 LTR order, 759 notLTR or mixture (5m57s elapsed)
[8b protein domains] 1657 candidate(s) with an LTR-list domain hit inside an LTR, 123 with a core polyprotein domain
[8c tRNA → PBS] 1410 tRNA gene(s), 361 PBS probe(s) (1m13s); 3440 candidate(s) with a 12-nt PBS, 2807 exact (7m11s elapsed)
[8d TSD, TG..CA] 2152 candidate(s) with a TSD of 4-20 bp (167 with one mismatch), 1352 with TG..CA ends
[9 decision] kept 5568 of 11246 candidate(s) at --keep Standard, 759 vetoed by protein hits, 146 by LTR self-repeat, 121 by a domain inside an LTR, 14 by poly-A-tailed LTRs (7m11s elapsed)
[10 overlaps] 1229 kept candidate(s) dropped as overlapping one with more evidence; 4339 left
[11 output] wrote 4339 element(s) to elements.fa
```

The kept copies as BED:

```bash
awk -F'\t' 'NR==1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
            $c["keep"] == 1 {print $1"\t"$2-1"\t"$3"\t"$c["evidence"]"\t"$c["identity"]}' cand.tsv > copies.bed
```

`--format gt` writes LTRharvest's format instead of the table.  The
evidence for each copy is in `cand.tsv`; there is no per-copy HTML report
yet (the evidence page below is per family and needs the library run).

### 3. Build a family library, with the evidence page

```bash
pipeline/ltr_library.py genome.fa -o out -t 24 --evidence-page
```

This runs step 2 itself (with the same settings; `out/cand.tsv` holds the
detected copies) and then clusters the kept elements into families. The
library is `out/library.fa`, one LTR and one internal model per family;
`out/families.tsv` describes each model, and `out/steps.log` times each
step, numbered on from ltr-forage's (12-24).

`out/evidence.html` has one dossier per family, drawn from the element it
was founded on: a track with LTRs, gaps, TSD, PBS, polypurine tract and
protein domains, the 5′ against 3′ LTR alignment, and the family's
full-length copies. Each family gets a verdict with its reasons:

| Verdict | Means |
|---|---|
| Well supported (Likely novel family with a reference) | at least three lines of evidence (two or more domains, TSD, TG..CA, PBS, copies on several sequences, many genome copies), or two domains and two more |
| Known family | with `--evidence-reference`: an LTR or internal model at least 80% identical over 80% of a curated family |
| Weak evidence | little beyond the LTR pair: the families to check by hand |
| Possible artifact | domains on both strands (a nested pair), copies sharing both flanks (a segmental duplication), or all copies on one sequence and few elsewhere |

`--evidence-reference curated_LTR.fa` measures each family against a curated
library. To add the page to a finished run, run the same command with
`--evidence-page --from evidence`. With a curated annotation of the genome,
`pipeline/evidence_page.py` can also select the families that mask most
outside it (see [Evidence page](#evidence-page)).

### LTRharvest's search only

```bash
ltr-forage genome.sufr --ltrharvest -j 24 > ltrharvest.out
```

This is LTRharvest's default search and output format (see
[LTRharvest mode](#ltrharvest-mode---ltrharvest)).

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

The tandem-array mask (`src/fastan.rs`) is a port of Gene Myers' FasTAN
([github.com/thegenemyers/FASTAN](https://github.com/thegenemyers/FASTAN)),
including the wave aligner it takes from daligner (Myers 2014).

The default protein profiles, `data/TEfam.hmm`, are LTR_retriever's
`database/TEfam.hmm` (102 Pfam profiles of transposable-element proteins;
Ou and Jiang 2018), distributed under GPL v3.

---

## How ltr-forage finds elements

LTR retrotransposons are mobile genetic elements bounded by a pair of nearly
identical Long Terminal Repeats (LTRs).  LTRharvest detects candidates by
exploiting the defining property that the two LTRs share a long, nearly-exact
sequence match; ltr-forage runs that search (steps 1, 3, 4, 5) and then
judges each candidate on evidence from its own sequence.  The step numbers
are the ones `ltr-forage` writes in its log and the figure shows.

| Step | Log label | What it does | In LTRharvest |
|---|---|---|---|
| 1 | `[1 index]` | load the sufr index (text only; the arrays are streamed) | yes (ESA) |
| 2 | `[2 tandem mask]` | hide tandem arrays from the seed search | no |
| 3 | `[3 seed pairs]` | maximal repeats at LTR spacing | yes |
| 4 | `[4 X-drop extension]` | extend each seed; length, spacing and similarity limits | yes (plus the TSD boundary move) |
| 5 | `[5 duplicates]` | drop candidates with identical outer ends | yes |
| 6 | `[6 re-extension]` | local alignment of the LTR pair sets the final ends and identity | no |
| 7 | `[7 near-duplicate merge]` | merge candidates whose final ends agree within 20 bp | no |
| 8a | `[8a self-alignment]` | LTRs that are tandem repeats themselves | no |
| 8b | `[8b protein domains]` | domains in the element and in each LTR (nail) | no |
| 8c | `[8c tRNA → PBS]` | the genome's tRNA genes, primer binding sites | no |
| 8d | `[8d TSD, TG..CA]` | TSD and terminal motif at the final ends | TSD only, at the X-drop ends |
| 9 | `[9 decision]` | keep or reject each candidate, with the reason | no |
| 10 | `[10 overlaps]` | resolve overlapping candidates | yes, by seed similarity before step 6 |
| 11 | `[11 output]` | table, element library, tRNA library | gt table only |

### 1. Index

`sufr create --dna` builds the suffix array, LCP array and text of all
sequences (see [Create a sufr index](#create-a-sufr-index)).  ltr-forage
loads only the text; the seed search streams the two arrays from the file.

### 2. Tandem-array mask (not in LTRharvest)

A tandem array of unit *p* yields a seed pair between almost every two of
its units within `--dist-max`, so satellite and centromeric arrays flood the
search: on a telomere-to-telomere gorilla assembly they gave 953 M seed
pairs against hg38's 11 M. Worse, a pair of units that passes extension is
re-extended across the array to the window limit and only then rejected by
the LTR self-alignment check, which is where most of the run time went.
ltr-forage finds tandem arrays with a port of Gene Myers' FasTAN
(`src/fastan.rs`) and hides them from the seed search; the sequence itself
is unchanged, so extension, TSDs and output still see it.

FasTAN scans each contig (a maximal A/C/G/T run) in 32 kb blocks 24 kb
apart. In a block, consecutive occurrences of the same 8-mer give seeds at
offsets *d* < 30,000; offsets holding more than max(1, *d*/64) seeds are
tried, most seeds first, by a local alignment of the sequence against
itself along diagonal *d*. The aligner is daligner's wave aligner: unit-cost
edit distance, furthest-reaching points per diagonal (so its cost grows with
the differences, not the area), confined to *d* ± 0.2 *d*, and extended while
the last 60 columns stay near 70% identity. An alignment spanning at least
1.95 *d* is an array of period about *d* (the mean diagonal over its trace
points); its bases are set to N so the next offsets and blocks skip it. An
LTR element's two LTRs align over far less than twice their spacing, so a
single element is never masked; the cost falls on genuine tandem copies of
an element and on satellites grown from an LTR, which keep their first
copies.

The port reproduces FasTAN v0.8's first pass exactly: the same 115,154 arrays
on rice (MSU7) and 613,069 on gorilla (mGorGor1 v2.1), base for base, with
the period and identity `tanbed` reports. It takes 47 s on rice and 386 s on
gorilla with 8 threads (FasTAN: 32 s and 542 s). `--fastan-out` writes the
arrays as BED (name, start, end, period, 1000 × (1 − differences / length)),
and `--fastan-only` writes them and exits. Arrays cover 21.6 Mb of rice
(5.8%) and 713.7 Mb of gorilla (20.1%).

With the mask, gorilla's candidate detection (no TSD adjustment and no
overlap removal before re-extension) dropped from 320,350 candidates
and 7.3 h to 24,283 candidates and 40 min at 8 threads, and the kept set
grew from 1,805 to 2,137 elements while their curated LTR support rose from
96% to 98%. `--mask-bed` adds intervals of one's own (name, 0-based start,
end) to the mask.
`--no-fastan` turns the mask off.

### 3. Seed pairs

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

### 4. X-drop extension

Each seed represents a short exact match that is the core of a putative LTR
pair.  It is extended in both directions using *X-drop gapped alignment*, a
heuristic that continues the alignment as long as the accumulated score remains
within a threshold *X* (`--xdrop`) of the running maximum.  The scoring scheme
(match award `--mat`, mismatch penalty `--mis`, indel penalties `--ins`/`--del`)
is configurable; the defaults match LTRharvest's defaults.

Extension is run twice per seed — leftward from the 5′ ends and rightward from
the 3′ ends — to determine the full LTR boundaries.

The X-drop ends then face LTRharvest's checks: both LTRs within
[`--ltr-min`, `--ltr-max`], their 5′ ends [`--dist-min`, `--dist-max`]
apart, and at least `--similar` percent similarity between the two LTRs
(Ukkonen's greedy edit distance, O(*n*·*d*) time, O(*d*) space).
LTRharvest first moves the outer ends to a target site duplication and
checks the moved ends; ltr-forage does that only in
[LTRharvest mode](#ltrharvest-mode---ltrharvest), and looks for the TSD at
the final ends instead (step 8d).

### 5. Duplicate removal

Multiple seeds can yield candidates with identical outer LTR boundaries.
These duplicates are removed by sorting on (sequence, left_ltr_start,
right_ltr_end) and collapsing identical records.

### 6. Re-extension (not in LTRharvest)

The X-drop extension (step 4) uses LTRharvest's scores (match +2,
mismatch −2, gap −3, threshold 5), so any indel of 2 bp or more between the
two LTR copies ends the extension there, and the bases beyond it are left out
of both LTRs even when they are otherwise identical. Re-extension recovers
them. For each candidate:

1. A window is cut around each LTR, reaching `--reextend-ltr-max` (default
   7000) plus `--reextend-vicinity` (default 100) bases beyond the current
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
boundaries.  On hg38, re-extension places 1,155 of the 1,929
RepeatMasker-annotated elements it finds within 20 bp of their annotated
boundaries, against 419 of 1,902 without it.

#### Window length (`--reextend-ltr-max`)

The window has to hold the longest LTR in the genome: an LTR longer than
the window is cut at the window edge (`window_edge` is 1 in the table).
The longest LTRs in the curated libraries we score against:

| Library | Longest LTR | Note |
|---|---|---|
| Rice (Jiang 6.9) | 6,144 bp | 99th percentile 5,825; 16% of LTRs over 2 kb |
| Lotus (Jiang 3.3) | 4,183 bp | |
| Human chr22 | 2,123 bp | |
| Arabidopsis (RepBase) | 2,033 bp | |
| Gorilla (curated) | 1,579 bp | |
| Zebra finch (curated) | 1,410 bp | |

The default, 7,000, covers all of them.  A shorter window is faster where
re-extension runs into masked tandem arrays: a pair of satellite units
re-extends across the array to the window edge before the self-alignment
check (step 8a) rejects it, and that alignment grows with the window.
At 2,500 bp, gorilla's candidate detection took 31 min instead of 106 (8
threads) and kept the same elements (2,135 against 2,136, 98.0% on curated
LTR annotation either way).  On rice the time barely changed (6.3 min
against 8.4), but 512 kept elements ended at a window edge instead of 68:
their LTRs were cut short, and 20 of the 23 elements lost had LTRs over
2.5 kb.

Keep the default for plants with large Gypsy LTRs, such as rice and
lotus; 2,500 is enough for the vertebrates and Arabidopsis above, and much
faster on assemblies rich in tandem arrays.  Set it from the longest LTR you expect, not
from the median.

### 7. Near-duplicate merge (not in LTRharvest)

Seeds of one element stop their X-drop runs a few to hundreds of bases
apart, so exact duplicates are rare; after re-extension, though, they end at
the same boundaries.  `--dedup-tol N` (default 20; 0 disables) merges
candidates whose final outer ends both lie within N bases of a candidate
with higher final LTR identity.  On rice with no overlap removal, 13,394 of
22,873 candidates have another within 20 bp after re-extension and none
before it.

### 8. Evidence annotation (not in LTRharvest)

Each candidate is annotated at its final ends.  The four annotations are
independent; 8a runs first so that 8b can skip what it rejects.

#### 8a. LTR self-alignment

Satellite arrays and Alu dimers with an ERV inserted between two monomers
pass the evidence rule as domain-bearing elements with a satellite or Alu
"LTR".  ltr-forage aligns each LTR
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
adds 40 s to an hg38 run.  It runs before the domain search, which skips the candidates it rejects
(`--domain-search-vetoed` searches them too; the decisions do not change).

#### 8b. Protein domains

Most of the candidates a whole-genome run reports are not LTR elements: on
hg38, 95% are pairs of SINEs, simple repeats or other similar sequence a few
kb apart. Of the signals we measured on hg38, rice and Arabidopsis, a
retrotransposon protein domain in the internal region separates the real
elements best.  ltr-forage searches each candidate's six-frame
translation against a HMMER3 file of profile HMMs with
[nail](https://github.com/TravisWheelerLab/nail): LTR_retriever's
`TEfam.hmm` (102 Pfam profiles, shipped as `data/TEfam.hmm`) by default,
or any HMMER3 amino-acid profile file given with `--domain-hmm`.
`--no-domains` skips the search.  On hg38 a hit of LTR order marks 68% of
the RepeatMasker-annotated elements and 0.5% of the other candidates.

ltr-forage translates each element plus `--domain-flank` bases (default 50)
on each side in six frames, with stop codons as X, writes the frames as a
protein FASTA, and searches it with nail's full pipeline (MMseqs2 seeds, cloud
search, Forward, composition-bias correction), taking six frames per
candidate as the database size for the E-value, as hmmsearch would, and
with nail's seed filter at P = 0.01 instead of its default 1e-4 (twice the
time, half as many weak hits lost against hmmsearch). Each LTR
is searched alone in the same run (see the coding-LTR veto below). Hits with
an E-value at or below `--domain-evalue` (default 0.05) count. `nail` and
`mmseqs` are found on PATH or given with `--nail` and `--mmseqs`; the frames
and nail's databases go to a subdirectory of `--temp-dir` (default:
`$TMPDIR`, else `/tmp`), removed afterwards.

On rice (39,924 frames, TEfam.hmm's 102 profiles, 8 threads) nail takes 36 s
and 206 CPU seconds at its default seed filter (70 s and 467 at 0.01),
hmmsearch 49 s and 250 CPU seconds, and the in-process
`libnail` search ltr-forage used before (its own 3-mer diagonal seeds, then
Forward on a band of rows) about 18,000 CPU seconds. At E ≤ 1e-5, 4,043
candidates have a nail hit and 4,044 a hmmsearch hit; each finds about 40
the other misses, against 95 hmmsearch candidates the in-process search
missed. The whole ltr-forage run on rice drops from 825 s to 106 s on 16
threads; it keeps 3,898 elements instead of 3,878.

The hits are classified as LTR_retriever's `annotate_TE.pl` does: profiles
belong to an LTR list or a not-LTR list (DNA transposons, LINEs, helitrons,
host genes); a candidate whose hits are all on the LTR list is `LTR`, all on
the other list `notLTR`, both `mixture`, and hits on neither list `NA`.
Superfamily is `Gypsy` when Copia-list hits are fewer than 0.3 of the
Gypsy-list hits, `Copia` when Gypsy-list hits are fewer than 0.9 of the
Copia-list hits, and strand `+` or `-` when hits on the other strand are
fewer than 0.4 of those on it.

#### 8c. tRNA genes and primer binding sites

Reverse transcription of an LTR retrotransposon is primed by a host tRNA
whose 3′ end, with its CCA tail, pairs with the primer binding site (PBS)
just inside the 5′ LTR. ltr-forage finds the genome's own tRNA genes and
looks for their 3′ ends at that position, so the probes come from the
genome being searched and no curated tRNA set is needed.

**tRNA genes.** `src/trna.rs` is a Rust port of the cytosolic tRNA search
in ARAGORN 1.2.41 (Laslett and Canback 2004) as run with
`aragorn -t -gc1 -l -i`: T-arm, acceptor-stem, D-arm and anticodon-stem
searches with ARAGORN's scores and thresholds, introns up to 3,000 bp in
the anticodon loop (`--trna-intron-max`), and its gene store and overlap
rules. ARAGORN's 20 kb windows are searched in parallel and the detections
replayed through its gene store in order, so the result equals a
single-threaded run. Against ARAGORN compiled with default floating point
(`-O2`), the calls are identical: 789 of 789 on Arabidopsis (whole
chromosomes and 2 Mb chunks) and 1,411 of 1,411 on rice chunks. A build with
`-Ofast` reorders floating-point sums and finds one more borderline gene on
rice (a tRNA-Trp with a 220 bp intron). The scan takes 39 s on Arabidopsis
and about 100 s on rice with 32 threads, under 5% of a rice run with the
domain search. `--trna-scan FASTA` runs it alone; `--no-pbs-trna` skips it.

**Probes.** Each tRNA gene of 70 to 120 bp with a standard amino acid gives
an 18-nt probe: the reverse complement of its last 15 nt plus CCA, which
starts with TGG. The test is the probe's first 12 nt, the tRNA's
3′-terminal 12 nt, with the TGG exact and at most one mismatch elsewhere,
starting 5 bp before to 15 bp after the first base past the 5′ LTR. Many
families pair over only these 12 to 13 nt: rice elements primed by the
initiator Met read TGGTATCAGAGC and then diverge, and an 18-nt test misses
them. The full probe's best edit distance (at most 2) is reported as a
secondary column. Both orientations are searched by default and the better
hit kept: the strand ltr-forage infers from the protein domains follows the
majority of hits, and for some rice Gypsy families those are accessory ORFs
(DUF828, DUF1263, Transposase_28) on the strand opposite the integrase,
while the PBS agrees with the integrase. Restricting the search to the
domain strand (`--pbs-domain-strand`) loses 306 PBS calls on rice, 280 of
them exact matches. `pbs_strand` gives the orientation the PBS was found in.

On rice, an exact 12-nt PBS sits after 56% of the kept elements that match
a curated family and after 0.1% of random genome windows (2.5% with one
mismatch allowed); unannotated singletons and those carrying non-LTR
transposon sequence have one at the random-window rate. A PBS is therefore
strong evidence that a 5′ LTR end is real, but its absence is not evidence
against an element: some families carry no recognizable PBS.

**tRNA library (`--trna-out FILE`).** FILE is a FASTA library with one entry
per distinct tRNA gene sequence, sorted by the number of kept elements it
primes:

```
>tRNA-Met(cat)_1 type=Met(cat) copies=7 intron=no len=72 probe12=TGGTATCAGAGC pbs_kept=1294 pbs_kept_exact=1161 pbs_kept_mm1=133 pbs_cand=1455 pbs_cand_exact=1296 pbs_cand_mm1=159 first=Chr2:780010-780081(+)
ATCAGAGTGGCGCAGCGGAAGCGTGGTGGGCCCATAACCCACAGGTCCCAGGATCGAAAC
CTGGCTCTGATA
```

That is the top entry on rice (MSU7), where the scan finds 1,410 tRNA genes,
1,057 distinct sequences and 361 probes, and the kept elements are primed
mostly by Met(cat) (1,294 elements), Arg(cct) (817), Arg(acg) (73),
Lys(ttt) (72), Arg(tct) (49) and Asp(gtc) (43); ltr-forage prints these
counts to stderr.

| Field | Meaning |
|---|---|
| ID | `tRNA-<type>_<n>`, numbered within each type in library order |
| `type` | Amino acid and anticodon as ARAGORN names them; ambiguous calls keep ARAGORN's `?(…)` form |
| `copies` | Genomic copies with this exact sequence (strand-independent: a copy is reported in its gene orientation) |
| `intron` | `yes` when the gene has an intron in the anticodon loop (the sequence includes it) |
| `len` | Gene length in bp, intron included |
| `probe12` | The 12-nt PBS probe: the reverse complement of the tRNA's 3′-terminal 12 nt (its last 9 nt plus CCA); `none` for genes outside 70 to 120 bp or without a standard amino acid, which make no probe |
| `pbs_kept`, `pbs_kept_exact`, `pbs_kept_mm1` | Kept elements whose best 12-nt PBS hit is this probe, in total, exact and with one mismatch |
| `pbs_cand`, … | The same over all candidates |
| `first` | Location of the first copy |

Usage is credited to probe groups: distinct genes that share their
3′-terminal 12 nt share one probe, and an element whose best hit is that
probe counts for each of them, so the per-entry counts add up to more than
the number of elements. The stderr summary counts each element once per
tRNA type. Every call, with the ID of its library entry, is written to
`FILE.tsv` (`seq`, `start`, `end`, `strand`, `name`, `length`, `intron`,
`score`, `sequence`, `library_id`; coordinates 1-based within the sequence,
start ≤ end).

`--no-pbs-trna` skips the tRNA scan; the PBS columns are then empty and no
candidate counts a PBS as evidence.

#### 8d. TSD and TG..CA

The TSD and the TG..CA motif are looked for at the final, re-extended
ends.  The TSD is the longest exact direct repeat of 4 to 20 bp flanking
the element with at most 1 bp of play at each end; failing that, a 5 or
6 bp repeat with one mismatch directly against both LTR ends.  The motif is
the first two bases of the 5′ LTR and the last two of the 3′ LTR.  Held to
the final ends, both keep their power to tell elements from chance: on
hg38 a 5 bp TSD appears on 17% of RepeatMasker-annotated elements (41% of
those over 95% identical) and on 1.3% of the other candidates.  A wider
search would find a short repeat almost everywhere.

That is why ltr-forage does not search for the TSD as LTRharvest does, at
the X-drop ends with a 60 bp window, moving the ends to it.  A 4 bp exact
match between the two 120 bp windows is almost always there by chance: of
1,235,994 gorilla seed pairs 7 lack one, of 93,050 on rice 5.  What the
search does is move the ends, and the similarity check on the moved ends
then rejects: 163,540 gorilla and 5,584 rice candidates pass at their X-drop
ends but fail once their ends are moved to the nearest 4 bp match.  On rice
those are mostly real elements (dropping the adjustment gains elements, 97%
of them on curated LTR annotation); on gorilla it halved the candidates,
yet removed the expensive satellite pairs no more often than any others
(52% of the candidates, 28% of the downstream time).  The tandem-array mask
(step 2) removes those directly.

### 9. Decision

Every candidate gets one decision, recorded in the table's `keep` and
`evidence` columns and used to select what `--library FILE` writes. The
rule has one knob, `--keep LEVEL`:

| Level | Kept when |
|---|---|
| `domain` | an LTR-order domain hit |
| `standard` (default) | an LTR-order domain hit, or TG..CA at the ends with a TSD of 4 bp or more |
| `loose` | an LTR-order domain hit, a TSD of 5 bp or more, or TG..CA with a TSD of 4 bp or more |
| `none` | always |

At every level ltr-forage rejects a candidate whose protein hits are of
not-LTR or mixed order; those are LINE, DNA-transposon and host-gene pairs,
and on hg38 1.3% of them are annotated LTR elements. `--internal-gap-max FRAC` adds a
second rejection for candidates whose internal region has a single-copy
stretch longer than that fraction of it (needs `--internal-stats`); it is
off by default because it costs recall in genomes with low-copy families
(Arabidopsis).

A third veto uses the domain search on each LTR alone. An LTR does not
encode the polyprotein, so a hit to a core polyprotein domain (integrase,
gag, reverse transcriptase, the aspartyl proteases, chromo; the list is
`CORE` in `src/domain.rs`) inside an LTR means the record is internal
sequence. The case it catches on rice is two adjacent copies of one
family called as a single element, with the repeated internal-plus-LTR
unit as the "LTRs": Chr4:6791971..6808083 is a 16 kb call whose two 7.2
kb LTRs are each a Truncator internal followed by 1.6 kb of Truncator
LTR, and its LTRs carry gag and protease hits at E-values down to 1e-69.
The self-alignment check cannot see these. On rice the veto removes 39
candidates: 18 whose two LTRs the Jiang library labels internal
sequence, 16 whose LTRs it labels LTR (10 of them annotated elements,
0.3% of the annotated set), and 5 unannotated. It leaves 22 of the 40
both-LTRs-internal candidates, whose short LTRs (106 to 2,400 bp) carry
no core domain. The veto counts only the core domains because other profiles on
the LTR-order list occur inside real LTRs: DUF834 (PF05754) and DUF1263
(PF06882) hit 748 candidate pairs of annotated LTRs on rice, and a veto
on the whole list removed 10% of the annotated elements. The two extra searches
add about two minutes to the ten-minute rice run; `--keep-coding-ltrs`
turns the veto off.

A fourth veto rejects a candidate whose two LTRs both end in a poly-A tail
(at least 12 A in the last 15 bp) or both start with a poly-T run (the same
on the other strand), unless its ends are TG..CA.  Real LTRs end in CA; two
copies of a SINE or a LINE 3′ end on either side of some sequence look like
an LTR pair, and on mouse they carry the domains of whatever lies between
them.  On mm10 the veto removes 878 of 29,032 candidates, and the library's
masking of SINEs outside the curated LTR annotation falls from 46 Mb to
0.4 Mb (library precision 42% before this veto and the model depth check,
step 22b of the family library, and 76% after).  It removes 14 candidates on rice, all of them
overlapping a call with more evidence, and 1 on Arabidopsis.

Measured against RepeatMasker annotation (strict or loose harness), with
the self-repeat check applied, the levels keep:

| Genome | Candidates | Annotated | `domain` | `standard` | `loose` |
|---|---|---|---|---|---|
| hg38 | 41,661 | 2,267 | 68.5% of annotated, 89.4% precise | 70.0%, 87.1% | 71.6%, 74.1% |
| Arabidopsis | 1,789 | 415 | 78.1%, 62.2% | 82.2%, 63.0% | 82.7%, 62.0% |
| rice | 7,092 | 3,216 | 84.3%, 70.4% | 86.6%, 70.9% | 87.1%, 70.6% |

The plant precisions are floors set by the annotation: when we inspected
the kept candidates the harness calls false there (analysis document,
sections 9 and 10), most were elements whose two LTRs the library names
differently, or nested copies. The motif route accepts a 4 bp TSD because
retroviruses make 4 to 6 bp target site duplications while plant Gypsy and
Copia elements mostly make 5: among candidates with no domain hit, TG..CA with a
5 bp TSD marks 6 annotated elements and 41 others on hg38, TG..CA with a
4 bp TSD 30 and 13; on rice the 5 bp case dominates (68 and 2 against 7 and
2). LTR_retriever ranks a 5 bp TSD first.  Of the 1,797 candidates the
earlier 5 bp rule kept on hg38 before the self-alignment check, 11 have an
LTR that is mostly satellite, against 1,091 of all 41,661 candidates, and
the self-alignment check (8a) removes 10 of the 11.

### 10. Overlap resolution

Overlapping candidates on the same sequence are resolved according to
`--overlaps`:

- `evidence` (default) — every candidate goes through re-extension,
  annotation and the decision; then, among overlapping kept candidates, the
  one with the most evidence (an LTR-order domain hit, a TSD, TG..CA, a PBS;
  then final LTR identity) stays.  A candidate rejected on its own evidence
  never displaces a kept one.  The dropped candidates stay in the table
  with `keep` 0 and `,overlap` added to `evidence`.
- `best` (LTRharvest's rule) — before re-extension, keep the candidate
  with the highest X-drop similarity in each group of overlapping
  candidates.
- `no` — before re-extension, discard every candidate in an overlapping
  group.
- `all` — keep them all.

`best` judges an overlapping group by the X-drop seed alignment, before
anything else is known.  On Arabidopsis and rice, `evidence` keeps 10% and
5% more elements than `best` at the same curated support, at 1.8–11× the
candidate-detection time (most of it spent re-extending and rejecting pairs
next to masked arrays).  In the library driver it gave higher sensitivity
on four of five genomes (see [Accuracy](#accuracy)).

### 11. Output

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

`--trna-out FILE` also writes the genome's tRNA library (see step 8c).
Without `--library` the table is the only output.

### LTRharvest mode (`--ltrharvest`)

`--ltrharvest` runs LTRharvest's search and nothing more: steps 1, 3, 4
(with the TSD search below), 5 and 10 (`best` by default; `no` and `all`
are accepted), and LTRharvest's output format.  No tandem-array mask,
re-extension, annotation or decision; `--library`, `--trna-out`,
`--fastan-out`, `--format tsv` and `--overlaps evidence` are refused.
The TSD and motif options `--mintsd`, `--maxtsd`, `--vicinity`, `--motif`,
`--motifmis` and `--tiered-tsd` exist only in this mode.

With `--mintsd` 2 or more (default 4), ltr-forage searches for *target
site duplications* — short exact direct repeats flanking the element at the
outermost LTR boundaries.  TSDs arise because the integrase creates a short
staggered cut in the host DNA during insertion; when the single-stranded gap is
repaired, the cut site is duplicated on both sides of the inserted element.

The search is performed by building a local suffix array of a window of width
`--vicinity` bp around the 5′ end of the left LTR, then querying each position
in a corresponding window around the 3′ end of the right LTR for an exact match
of length in [`--mintsd`, `--maxtsd`] (default 20), with `--vicinity` 60 by default.  This mirrors GenomeTools'
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

---

## Differences from GenomeTools LTRharvest

| Feature | LTRharvest (GenomeTools) | ltr-forage |
|---|---|---|
| Index format | GenomeTools ESA (`suffixerator`) | sufr suffix array (`.sufr` file) |
| Seed enumeration | Bottom-up ESA traversal (C), single-threaded | Same algorithm, streaming the sufr SA/LCP from disk, parallel over LCP-delimited pieces (`-j`) |
| Seed extension | Single-threaded | Parallelised with Rayon (`-j`) |
| Tandem arrays | Seeded like any repeat | Masked from the seed search with a FasTAN port (`--no-fastan` to seed them) |
| TSD | Searched at the X-drop ends, which move to it; required | Searched at the final ends and used as evidence (`--ltrharvest`: LTRharvest's search, or LTR_retriever-style ranking with `--tiered-tsd`) |
| LTR boundaries | X-drop ends | Local alignment of the LTR pair (re-extension) |
| Near duplicates | Exact outer boundaries only | Also final ends within `--dedup-tol` bp |
| Annotation and decision | None | LTR self-alignment, protein domains, PBS, TSD and TG..CA; keep or reject with the reason |
| Overlaps | Best seed similarity per group, before anything else | Most evidence among kept candidates, after annotation (`--overlaps best` for LTRharvest's rule) |
| Output | Space-separated table | Per-candidate table and element library (`--format gt`, or `--ltrharvest`, for LTRharvest's table) |

The seed enumeration (`maxpairs.rs`) and X-drop extension (`xdrop.rs`) are
faithful Rust ports of the corresponding GenomeTools routines
(`esa-maxpairs.c`, `esa-bottomup-maxpairs.inc`,
`gt_evalxdroparbitscoresextend`), including GenomeTools' treatment of N as a
wildcard that matches nothing.  The TSD search of LTRharvest mode (`tsd.rs`)
mirrors `searchforbestTSDandormotifatborders` and
`gt_sarrquerysubstringmatch`, including the tie-breaking emission order
required to reproduce GT-compatible boundary choices.

The first structural departure is the use of the `sufr` index rather than a
GenomeTools ESA.  Because `sufr` stores a combined suffix array over all input
sequences with separators excluded from the active suffix count, the bottom-up
traversal reads directly from `sufr`'s pre-built SA and LCP arrays without the
intermediate ESA query layer present in GenomeTools.  An LCP shift
(look-ahead vs. look-back convention) is applied during traversal to match the
GenomeTools array layout.

The suffix array and LCP array are never loaded into memory.  The seed
search (step 3) consumes each of them once, in order, so `ltr-forage`
streams both from the `.sufr` file through a small buffer (`stream.rs`) and
holds only the text in memory.  Peak memory is therefore close to the genome
size rather than to the index size.

The seed search also runs in parallel.  A suffix array position whose LCP
value is below `--min-seed` cannot lie inside any repeat of seed length or
more, so cutting the array at such positions gives pieces that the traversal
can process independently.  `ltr-forage` cuts the array into about four
pieces per thread (`stream::split_points`), traverses each piece on its own
worker with its own file readers, and concatenates the seeds in piece order.
That is the order a single sequential traversal would emit them, so the
seed-index tie-break below and the output are unchanged by the thread count.

X-drop extension (step 4) is embarrassingly parallel: each seed is
independent.  A deterministic seed-index tiebreaker is attached before the
parallel scatter and used during the subsequent sort so that the surviving
candidate from any set of seeds with identical outer boundaries is always
the one with the lowest SA-enumeration index — the same choice a
single-threaded stable sort on insertion order would make.  The later
per-candidate steps keep the candidates' order.

---

## Installation

### Prerequisites

- Rust toolchain ≥ 1.80 ([rustup](https://rustup.rs) recommended)
- The `sufr` tool for building `.sufr` index files
  ([TravisWheelerLab/sufr](https://github.com/TravisWheelerLab/sufr))
- For the protein-domain search (step 8b): `nail` ≥ 0.7.1
  ([TravisWheelerLab/nail](https://github.com/TravisWheelerLab/nail)) and
  MMseqs2 (`mmseqs`), which nail seeds with; not needed with `--no-domains`
  or `--ltrharvest`.  The protein profiles ship in `data/TEfam.hmm`; the
  binary looks for them in the checkout it was built from, so give
  `--domain-hmm` if you move it.

For the family-library driver (`pipeline/ltr_library.py`) also:

- Python ≥ 3.8 with numpy
- The Rust port of RMBlast, `rmblastn` ≥ 3.0.9
  ([Dfam-consortium/rmblast](https://github.com/Dfam-consortium/rmblast));
  NCBI's C++ `rmblastn` does not take the options the driver uses; 3.0.11 or later
  is recommended (3.0.10's streaming `--mt-mode 1` can deadlock, and the driver
  falls back to `--mt-mode 0` with it)
- `te-composer` from dfam-curator ≥ 0.6.6 (consensus models)

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
ltr-forage --ltrharvest [LTRharvest options] <SUFR>
ltr-forage --trna-scan <FASTA> [-j N] [--trna-intron-max BP]

Arguments:
  <SUFR>  Path to a .sufr index file (created with `sufr create --dna`)

Seed and extension options:
  -m, --min-seed <BP>       Minimum seed length (min exact match) [default: 30]
  -l, --ltr-min <BP>        Minimum LTR length [default: 100]
  -L, --ltr-max <BP>        Maximum LTR length at the X-drop stage [default: 1000]
  -d, --dist-min <BP>       Minimum 5'-to-5' distance between LTRs [default: 1000]
  -D, --dist-max <BP>       Maximum 5'-to-5' distance between LTRs [default: 15000]
  -s, --similar <PCT>       Minimum percent similarity between LTR copies [default: 85]
  -x, --xdrop <SCORE>       X-drop threshold for extension alignment [default: 5]
      --mat <MAT>           Match score (positive) [default: 2]
      --mis <MIS>           Mismatch score (negative) [default: -2]
      --ins <INS>           Insertion score (negative) [default: -3]
      --del <DEL>           Deletion score (negative) [default: -3]

Tandem-array mask (step 2):
      --no-fastan           Search tandem arrays for seeds too
      --fastan-out <FILE>   Write the masked arrays as BED (name, start, end, period, 1000 x identity)
      --fastan-only         Only find the arrays (writes --fastan-out) and exit
      --mask-bed <FILE>     Also mask these intervals (BED: name, 0-based start, end)

Re-extension and near duplicates (steps 6, 7):
      --reextend-ltr-max <BP>   Longest LTR the alignment window reaches beyond each X-drop
                            end [default: 7000]
      --reextend-vicinity <BP>  Extra window beyond that on the outer side of each LTR [default: 100]
      --trim-window <N>     Terminal alignment columns that must be clean for an end to stand [default: 6]
      --trim-matches <N>    Matches required inside that terminal window [default: 6]
      --dedup-tol <BP>      Merge candidates whose final outer ends both lie within BP of a
                            better one's; 0 disables [default: 20]
      --internal-stats      Count 20-mers of each internal region against the index and
                            report coverage and the longest single-copy stretch

Protein domains (step 8b):
      --domain-hmm <FILE>   HMMER3 file of protein profiles [default: data/TEfam.hmm]
      --no-domains          Skip the domain search
      --domain-evalue <E>   E-value threshold for a domain hit [default: 0.05]
      --domain-flank <BP>   Bases added on each side of the element before translation [default: 50]
      --nail <PATH>         nail binary (a path, or a name on PATH) [default: nail]
      --mmseqs <PATH>       MMseqs2 binary nail seeds with [default: mmseqs]
      --temp-dir <DIR>      Directory for scratch files [default: $TMPDIR, else /tmp]
      --domain-search-vetoed  Search candidates the self-alignment check rejected too
      --keep-coding-ltrs    Keep candidates with a polyprotein domain hit inside an LTR

Primer binding site (step 8c):
      --no-pbs-trna         Skip the genome tRNA scan behind the PBS
      --pbs-domain-strand   Search for the PBS only on the strand the protein domains give (default: both)
      --trna-out <FILE>     Write the genome's tRNA library (FASTA) and every call (FILE.tsv)
      --trna-intron-max <BP>  Longest tRNA intron searched for [default: 3000]
      --trna-scan <FASTA>   Only scan FASTA for tRNA genes, print the calls and exit

Decision and overlaps (steps 9, 10):
      --keep <LEVEL>        Evidence a candidate needs: domain | standard | loose | none [default: standard]
      --tandem-max <FRAC>   Reject a candidate when either LTR aligns to itself off the diagonal
                            over at least this fraction of its length [default: 0.99]
      --internal-gap-max <FRAC>  Also reject candidates whose internal region has a single-copy
                            stretch longer than this fraction of it (needs --internal-stats)
      --overlaps <STRATEGY> evidence | best | no | all [default: evidence; best with --ltrharvest]

Output (step 11):
      --format <FORMAT>     tsv | gt [default: tsv; gt with --ltrharvest]
      --library <FILE>      Write the kept elements as a FASTA library (left LTR, internal
                            region, right LTR per element)
  -j, --threads <N>         Threads (0 = all CPUs) [default: 0]
  -h, --help                Print help

LTRharvest mode:
      --ltrharvest          LTRharvest's search only: TSD boundary move, overlaps best, gt output
      --mintsd <BP>         Minimum TSD length; 0 or 1 disables the TSD search [default: 4]
      --maxtsd <BP>         Maximum TSD length [default: 20]
      --vicinity <BP>       Search window around the LTR boundaries for TSD and motif [default: 60]
      --motif               Require the TG...CA motif at the LTR ends
      --motifmis <N>        Allowed mismatches in the motif [default: 0]
      --tiered-tsd          Rank TSD candidates LTR_retriever's way instead of by smallest shift
```

### Example — LTRharvest's defaults

```bash
ltr-forage hg38.sufr --ltrharvest -j 32 > hg38_ltrs.out
```

which is, with every LTRharvest default written out,

```bash
ltr-forage hg38.sufr --ltrharvest \
    --min-seed 30 \
    --ltr-min 100 --ltr-max 1000 \
    --dist-min 1000 --dist-max 15000 \
    --similar 85 \
    --xdrop 5 \
    --mat 2 --mis -2 --ins -3 --del -3 \
    --mintsd 4 --maxtsd 20 \
    --vicinity 60 \
    --overlaps best \
    -j 32 > hg38_ltrs.out
```

---

## Output Format

The default output is the per-candidate table below.  `--format gt` (the
default with `--ltrharvest`) is compatible with GenomeTools LTRharvest.
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

### Per-candidate table

The table (`--format tsv`, the default) has a header line and one tab-separated row per candidate,
with every value the later stages decide on, so a decision can be inspected
rather than inferred. Coordinates are 1-based and inclusive, relative to the
sequence named in the first column.

| Column | Description |
|---|---|
| `seq` | Sequence name from the index |
| `elem_start`, `elem_end`, `elem_len` | Element span |
| `lltr_start`, `lltr_end`, `lltr_len`, `rltr_start`, `rltr_end`, `rltr_len` | LTR spans after re-extension (the X-drop spans when no alignment qualified) |
| `identity` | Identity of the trimmed alignment in percent, gaps counted as differences; the X-drop similarity when not re-extended |
| `xdrop_lltr_start` … `xdrop_rltr_end`, `xdrop_similarity` | The boundaries and similarity the X-drop stage produced |
| `reextended` | 1 when the row's boundaries come from the alignment |
| `aln_score`, `aln_columns` | Raw score and aligned columns of the trimmed alignment |
| `trim_l5`, `trim_r5`, `trim_l3`, `trim_r3` | Bases removed by the end trim at the 5′ end of the left and right LTR and at their 3′ ends |
| `window_edge` | 1 when the alignment reached within 5 bp of a window edge that was not a sequence end: the LTR may be longer than the window, or the pair sits in a larger duplication |
| `tsd_len`, `tsd_seq`, `tsd_shift_l`, `tsd_shift_r` | Longest exact direct repeat (4 to 20 bp) flanking the element, allowing 0 or 1 bp of play at each end; failing that, a 6 or 5 bp repeat with one mismatch directly against both LTR ends (see `tsd_mm`). `tsd_seq` is the left copy; `.` and 0 when none |
| `motif`, `motif_class` | First two bases of the left LTR and last two of the right LTR, and LTR_retriever's reading of them: `TGCA` for the canonical termini, `known` for the other pairs its rules accept (TGCT, TACA, TACT, TGGA, TGGT, TATA, TGTA, TGCC), `other` |
| `ppt_purines15`, `ppt_offset` | Purine count in the most purine-rich 15 bp window within 30 bp before the right LTR, and that window's start relative to the right LTR |
| `pbs_offset` | Start of the primer binding site relative to the first base after the 5′ LTR (−5 to +15), from the 12-nt tRNA match described under step 8c; −1 when none or with `--no-pbs-trna` |
| `n_count` | Ns inside the element |
| `xdrop_tsd_len` | Always 0: the TSD search at the X-drop ends runs only in LTRharvest mode, which writes no table (the column is kept so the others do not move) |
| `int_tiles`, `int_multi_frac`, `int_gap_frac`, `int_gap_bp` | With `--internal-stats`: number of 20-mer tiles taken every 10 bp across the internal region (tiles containing N skipped), the fraction seen at least three times in the genome on either strand, and the longest run of single-copy tiles as a fraction of the tiles and in bases; `.` otherwise |
| `dom_order`, `dom_superfamily`, `dom_strand` | `LTR`, `notLTR`, `mixture`, `NA` or `none`; `Gypsy`, `Copia` or `unknown`; `+`, `-` or `?`, as described under step 8b; `.` with `--no-domains` or for a candidate the self-alignment check rejected |
| `dom_best_e`, `dom_hits` | Best E-value, and every hit as `accession:E-value:frame` (frames `aa1`..`aa3` forward, `rev_aa1`..`rev_aa3` on the reverse complement), comma separated; `.` when none |
| `ltr_dom_hits` | Hits to the LTR-order (polyprotein) profiles when each LTR is searched on its own, summed over the two LTRs; 0 when `--keep-coding-ltrs` skips the search |
| `lltr_self_cov`, `rltr_self_cov` | Fraction of each LTR covered by plus-strand alignments of the LTR against itself off the main diagonal; 1.00 for a satellite monomer or an Alu dimer, 0 for unique sequence |
| `keep`, `evidence` | 1 when the candidate is kept, and what decided it: `domain`, `tsdN`, `tgca` (comma separated), `veto:notLTR`, `veto:mixture`, `veto:tandemF`, `veto:ltrdomainN`, `veto:polyA`, `gapF` when `--internal-gap-max` rejected it, or `none`; `,overlap` is appended when step 10 dropped a kept candidate for an overlapping one with more evidence |
| `pbs_trna`, `pbs_mm12` | The tRNA type(s) whose 3′-terminal 12 nt match the PBS, e.g. `Met(cat)` (several joined by `/` when they share the probe), and the mismatches (0 or 1); `.` when none |
| `pbs_edist18` | Best edit distance (0 to 2) of a full 18-nt tRNA probe over the same window; `.` when none |
| `pbs_strand` | `+` when the PBS was found after the left LTR, `-` when before the right LTR on the reverse strand; `.` when none |
| `tsd_mm` | 1 when the TSD is the one-mismatch fallback, else 0. It is enriched about 20-fold over shifted flanks only for TG..CA elements (lotus), so it counts toward `--keep standard`/`loose` only alongside TG..CA; the evidence tag reads `tsd5mm1` |

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

## Building a family library

`pipeline/ltr_library.py` runs ltr-forage on a genome and turns its
candidates into a family library for RepeatMasker or RepeatModeler's
LTRPipeline.  It finds integrase-type LTR retrotransposons (Gypsy, Copia,
Pao/BEL, ERVs).  Tyrosine-recombinase elements (DIRS, Ngaro), which
RepeatMasker files under the LTR class, have no pair of direct LTRs and are
not found.

```bash
pipeline/ltr_library.py genome.fa -o out -t 32
```

It finds `ltr-forage`, `sufr`, `rmblastn`, `te-composer`, `nail` and `mmseqs`
from the `--ltr-forage`, `--sufr-bin`, `--rmblastn`, `--te-composer`, `--nail`
and `--mmseqs` options, then the `LTRFORAGE`, `SUFR`, `RMBLASTN`,
`TE_COMPOSER`, `NAIL` and `MMSEQS` environment variables, then `PATH`, then
the usual install locations and checkouts next to this one. It skips a
binary older than it needs (rmblastn 3.0.9, te-composer 0.6.6, nail 0.7.1),
so an old install on `PATH` does not shadow a newer one. The protein
profiles default to the shipped `data/TEfam.hmm` (`--domain-hmm`). A sufr
index beside the genome (`genome.sufr` for `genome.fa`) is used when it
exists; `--sufr` names another. Each run writes the command and every
program's path and version to `out/versions.txt`, and its step times to
`out/steps.log`. The step scripts live beside the driver in `pipeline/`.

### Steps

The numbers continue ltr-forage's and appear in `out/steps.log` and the
figure.

0. **Depth sample.**  A genome sample for copy depth: the whole genome up
   to 200 Mb, else 10% of it in 100 kb windows.  Depth thresholds scale
   with genome size and the sample fraction.
1. **Index.**  A sufr index (skipped when one exists).
1–11. **ltr-forage** with its defaults (see
   [How ltr-forage finds elements](#how-ltr-forage-finds-elements);
   `--overlaps best` for LTRharvest's overlap rule, `--no-fastan` to seed
   tandem arrays), giving the candidate table (`cand.tsv`), the masked
   arrays (`tandem.bed`) and a raw library of each kept element's LTR and
   internal region.
12. **Pre-clustering.**  The raw library is clustered all-vs-all with
    rmblastn, LTRs against LTRs and internals against internals (the
    clustering reads only same-kind hits of 80% identity or more, and drops
    the rest as they are written).
13. **Consensus representatives.**  Each cluster of at least `--min-cons`
    (3) elements is represented from here on by a te-composer consensus,
    which drops alignment columns held by fewer than half the median number
    of copies.  This keeps an insertion carried by a few copies out of the
    model.  A cluster with more than 150 LTR or internal records is built
    from a seeded random sample of 150, the same on every run: a random
    sample keeps the mix of full-length and truncated copies the column rule
    depends on, and the consensus matches the full one to within 31 bases
    on mm10's 22 such clusters.  A 2,006-copy mm10 internal cluster took
    107 minutes, and now 40 s.
14. **Depth cleanup.**  Three rounds (`--depth-rounds`) of rmblastn search
    of the sample against the library, one single-thread batch per thread.  Stretches far
    above their record's own depth are excised as nested passengers: 3x for
    an element, 10x for a consensus representative, whose family is
    high-copy itself.  High-copy LTRs are dropped, but never a consensus
    representative.  The search uses mask level 80, so a passenger that
    many records share splits its genome copies among them, and each record
    may see too little depth to cut it (step 22b catches these).
15. **Nested-LTR vote.**  A first element clustering excises LTRs nested in
    internals by vote.
16. **Element clustering.**  Clusters and one exemplar per cluster; each big
    cluster's members rejoin the cluster of their representative.
17. **Flank filter.**  A cluster is dropped when most of its pairs, or
    every one of its copies, share both flanks (300 bp at 80% identity)
    with another copy: the copies are one duplicated locus, not independent
    insertions (`--no-flank` skips this).
18. **Containment.**  Exemplars contained in another at 90% or more are
    dropped.
19. **Internal nested excision.**  Two rounds of nested-LTR removal from
    exemplar internals.
20. **PBS** of each exemplar, from ltr-forage's tRNA search.
21. **Tiering.**  Multi-element clusters are kept; singletons only with LTR
    identity ≥ 93% or a primer binding site (`--tier strict`: identity ≥ 90%
    and a PBS).  `--early-tier` applies the singleton rule to the raw
    clusters before step 14 instead, which is faster but not the default
    (see [Rule ablation](#rule-ablation)).
22. **Consensus models.**  A consensus model for each cluster of at least
    `--min-cons` elements and the exemplar element for smaller ones (at
    most 150 sampled copies, as in step 13).
22b. **Model depth check.**  One more depth round on the finished models:
    a consensus model is cut where its depth reaches 10x its median, an
    exemplar at 5x (`--model-depth-ratio`); no model is dropped
    (`--no-model-depth` skips the step).  It removes passengers step 14
    left shared among many records, and passengers a consensus took back
    from copies its founder had lost them from: on mm10, L1 and B1/B2
    pieces in a few models masked 191 Mb outside the curated annotation.
    On Arabidopsis it cuts a complete 2 kb ATGP2N solo LTR (39x the
    model's depth) out of an ATCopia76 internal.  At 3x it also cut
    conserved stretches of old, 30%-diverged Gypsy families at 3–4x, with
    no gain in precision anywhere.  The same search drops an internal model
    whose genome copies exceed 2,000 per element of its cluster
    (`--int-copies-max`; 0 turns it off).  Internals are never solo, so
    such a model is mostly another repeat: with `--early-tier`, an mm10
    exemplar internal of L1 sequence had 5,999 copies per element and
    masked 133 Mb outside the curated LTRs.  No other model in the six
    genomes' default runs went above 707, so the check changes none of
    them.
23. **80–80 families**, LTR and internal models grouped separately (and
    searched separately), with Gypsy and Copia classes from the domain hits.
23b. **Seed alignments** (`library.stk`): each model's best
    `--seed-copies` (50) non-overlapping copies in the depth sample, or its
    origin element when the sample holds none, aligned to the model by
    `te-composer --fixed-cons`.  The model is not refined, so each
    alignment's `#=GC RF` line is the model itself.  Rows are named by
    genome range.  On rice: 865 alignments, median 23 rows, 51 s.  A
    te-composer without `--fixed-cons` (before 0.6.8) gives one-row
    alignments.
24. **Evidence page** (`--evidence-page`).

`--stage-a` skips steps 12–13 and depth-cleans and re-clusters every
element, as the driver did before.  It is slower on large libraries and no
more accurate.  `--min-cons 0` gives exemplar models only.  Each step's
command lines and output go to `out/logs/`, and timings to
`out/steps.log`.  A rerun skips finished steps; `--from STEP` reruns from a
step on.

### Output

| File | Contents |
|------|----------|
| `library.fa` | The family library: `fam<N>_v<k>_LTR` and `fam<N>_v<k>_INT` records, classified `#LTR/Gypsy`, `#LTR/Copia` or `#LTR` |
| `library.stk` | A Stockholm seed alignment per model, with the model as `#=GC RF` (step 23b) |
| `families.tsv` | Each model's family, origin element, kind, length, status and class |
| `cand.tsv` | ltr-forage's per-candidate table, with every value behind its decision (see [Per-candidate table](#per-candidate-table)) |
| `sd_clusters.tsv` | Clusters tested for segmental duplication, and the call |
| `evidence.html` | With `--evidence-page`: one dossier per family (see below) |

### Evidence page

`--evidence-page` writes `evidence.html`, a self-contained page with one
dossier per family, drawn from the family's founding element: LTRs, assembly
gaps, TSD, primer binding site, polypurine tract and protein domains on a
genome-browser style track, the 5′ and 3′ LTR alignment, and the family's
full-length copies and the sequences they lie on.  Each family gets a
triage verdict (likely novel or well supported, known family, weak
evidence, possible artifact) with the evidence behind it.
`--evidence-reference LIB` adds the nearest family in a curated library.

Run by hand, `pipeline/evidence_page.py` also takes a RepeatMasker run of
the genome with the library (`--mask`), the curated annotation
(`--standard`) and another tool's mask (`--other NAME=OUT`).  It then shows
each family's masking inside and outside the curated LTR annotation and how
much of the rest the other tool finds too, and `--families top-fp:N`
selects the families with the most masking outside it:

```bash
pipeline/evidence_page.py out genome.fa -o fp.html --families top-fp:30 \
    --reference curated.fa --mask genome.fa.out --standard curated.fa.out \
    --other HiRepEat=hirepeat.fa.out
```

### Accuracy

Library masking was scored against curated libraries: sensitivity and
precision over the bases the curated LTR annotation covers.  Recovered
families are the curated LTR and internal models covered at 80% or more.
All ltr_library.py rows are runs of the current driver (Arabidopsis, rice,
mouse and gorilla 2026-10-25; zebra finch and lotus 2026-10-28), rerun
from ltr-forage on 2026-10-30 without the domain-only identity floor (see
[Rule ablation](#rule-ablation)); gorilla was unchanged.  Times are from the
full runs, which the floor barely affects.  Driver
times are with 16 threads (Arabidopsis, rice, zebra finch, lotus) or 24
(mouse, gorilla), HiRepEat's with 16.

| Genome | Library | Sensitivity | Precision | Recovered LTR / INT | Time |
|--------|---------|-------------|-----------|---------------------|------|
| Arabidopsis (TAIR10.1) | ltr_library.py | 91.89% | 95.61% | 120 / 122 | 50 min |
| Rice (MSU7) | ltr_library.py | 93.44% | 98.43% | 309 / 308 | 85 min |
| Rice (MSU7) | HiRepEat | 83.10% | 99.61% | 80 / 76 | 19 min |
| Mouse (mm10) | ltr_library.py | 42.93% | 76.20% | 105 / 41 | 124 min |
| Mouse (mm10) | HiRepEat | 24.04% | 65.46% | 61 / 30 | 31 min |
| Zebra finch (bTaeGut7) | ltr_library.py | 75.21% | 86.71% | 25 / 29 | 59 min |
| Zebra finch (bTaeGut7) | HiRepEat | 70.58% | 90.86% | 16 / 27 | 20 min |
| Lotus (ChineseLotus 1.1) | ltr_library.py | 92.43% | 93.78% | 283 / 298 | 303 min |
| Gorilla (mGorGor1 v2.1) | ltr_library.py | 36.78% | 97.96% | 2 / 13 | 107 min |
| Gorilla (mGorGor1 v2.1) | HiRepEat | 17.68% | 95.97% | 1 / 3 | 72 min |

On mouse the library was at 44.29% and 42.13% before the poly-A veto and
the model depth check: 141 Mb of its false-positive masking was L1 and 46 Mb
B1/B2 SINEs, now 35 Mb and 0.4 Mb.  Most of the L1 that remains comes from
two LTR models made of L1 sequence that carry no poly-A tail.  The copy cap
took the mouse run from 266 to 210 minutes; splitting the all-vs-all
searches by record kind, dropping hits the clustering does not read and
running the depth search on every thread took it to 124 minutes, with
identical clusters.  On Arabidopsis and rice the new rules cost
0.9 and 0.3 points of sensitivity at equal precision; on zebra finch 1.5
points at equal precision (76.64% / 86.59% before), and on lotus 0.9 points
for 2.4 points of precision (93.36% / 91.66%, 288 / 296 families before).  On Arabidopsis that
is mostly the ATGP2N solo LTR that step 22b cuts out of an ATCopia76
internal: it had masked 40 ATGP2N copies (81 kb) as ATCopia76.

Against the previous driver (LTRharvest's TSD boundary move, a k-mer
tandem mask and `best` overlaps), sensitivity rose on rice (+1.3 points),
zebra finch (+2.2), lotus (+0.7) and gorilla (+2.9) and fell 0.3 on
Arabidopsis; precision rose on zebra finch (+2.3) and fell 1.7 points on
lotus.  Most of lotus's new false-positive masking comes from families with
retrotransposon domains, TSDs or primer binding sites, and many genome
copies, that the curated library does not annotate.

Much of the zebra finch false-positive masking is in regions the curated
library leaves unannotated.  About 95% of HiRepEat's falls there too, and
several of these models are dispersed elements with retrotransposon domains,
which suggests families missing from the curated library rather than
artifacts.

#### Rule ablation

Several of the driver's rules were added while tuning on one genome.  To
check whether each one generalizes, every rule was switched off on its own
on all six genomes, with the rest of the run unchanged (2026-10-29/31, the
runs in the table above).  The table gives the change in sensitivity /
precision, in points, with the rule off:

| Rule off (added for) | Arabidopsis | Zebra finch | Rice | Gorilla | Lotus | Mouse |
|---|---|---|---|---|---|---|
| Step 22b model depth check (mouse) | +0.87 / +0.04 | +1.52 / −0.10 | +0.25 / −0.05 | +0.61 / −0.16 | +0.93 / −2.43 | +1.50 / −32.47 |
| Step 17 flank filter (zebra finch) | +0.16 / −0.82 | +0.21 / −2.71 | +0.04 / −0.09 | 0 / −0.06 | 0 / −0.01 | +1.59 / −38.21 |
| Step 21 tier (`--tier all`, rice) | +3.97 / −7.15 | +5.45 / −8.08 | +2.83 / −5.85 | +2.09 / −1.75 | +5.40 / −8.52 | +2.65 / −42.39 |
| Poly-A veto (`--no-polya-veto`, mouse) | 0 / 0 | 0 / 0 | 0 / 0 | 0 / −0.17 | 0 / 0 | −0.04 / −20.32 |
| Domain-only identity floor (gorilla) | −0.05 / +0.05 | +0.10 / +0.02 | −0.19 / +0.01 | 0 / 0 | −0.01 / −0.32 | −0.08 / −0.05 |

- **Flank filter and poly-A veto.** These cost no sensitivity anywhere, and
  without them mouse loses 38 and 20 points of precision.  Mouse stands for
  genomes rich in SINEs, LINEs and segmental duplications, most mammals
  among them, not for a quirk of one assembly.
- **Step 22b.** It costs 0.25–1.5 points of sensitivity on every genome.
  It keeps mouse precision from falling to 44%, and is worth 2.4 points on
  lotus.
- **Tier.** It trades 2–5.5 points of sensitivity for 2–8 points of
  precision on every genome, and 42 on mouse.  `--tier all` keeps every
  singleton, for users who want sensitivity first.
- **Domain-only identity floor.** It changed no genome by more than 0.3
  points and was removed.

Two simplifications were tested the same way (2026-10-31), rerunning each
genome from the changed step:

| Change (sensitivity / precision) | Arabidopsis | Zebra finch | Rice | Gorilla | Lotus | Mouse |
|---|---|---|---|---|---|---|
| Step 14 with one round | 92.02 / 95.57 | 75.87 / 84.08 | 93.53 / 98.24 | 39.71 / 64.37 | 92.51 / 89.55 | 46.01 / 46.98 |
| No step 14 (step 22b only) | 92.68 / 94.65 | 79.59 / 83.17 | 93.65 / 97.81 | 40.23 / 51.83 | 92.48 / 85.15 | 46.20 / 44.68 |
| `--early-tier` | 93.05 / 95.37 | 75.10 / 85.24 | 92.78 / 98.35 | 35.18 / 97.90 | 92.41 / 95.26 | 42.71 / 45.69 |
| Base (with the identity floor) | 91.94 / 95.56 | 75.11 / 86.69 | 93.63 / 98.42 | 36.78 / 97.96 | 92.44 / 94.10 | 43.01 / 76.25 |

- **Depth rounds.**  Rounds 2 and 3 of step 14 are needed: each round's
  cuts expose pieces the next round can cut, and step 22b cannot take
  their place.  Without them gorilla and mouse lose 34–46 and 29–32 points
  of precision.
- **Early tier.**  It drops the singletons step 21 would reject before
  the depth cleanup, which halves the time from step 13 on (rice 21
  instead of 45 minutes, lotus 67).  The results are mixed: rice and
  gorilla lose 0.9 and 1.6 points of sensitivity and zebra finch 1.45
  points of precision; Arabidopsis gains 1.1 points of sensitivity and
  lotus 1.2 of precision.  On mouse it first cost 30 points of precision:
  it changes which elements cluster together, and an 8-element cluster
  that had been part of a 578-element one became its own family, whose
  exemplar internal, mostly L1 sequence, masked 133 Mb outside the curated
  annotation.  Step 22b's internal-copies check now drops that model, and
  mouse reaches 42.65% / 79.40% (default 42.93% / 76.20%).  It stays off
  by default; it is worth trying on large genomes, where the time matters
  most.

#### Held-out genomes

The driver was then run with its defaults, unchanged, on genomes that were
never used for tuning.  The standards are RepeatMasker runs:

- fruit fly, *C. elegans* and zebrafish: the curated Dfam 3.9 families of
  the species and its ancestors;
- soybean: RepBase 31.09's entries for the species.

All runs used 8 threads.

| Genome | LTR share | Sensitivity | Precision | Recovered LTR / INT | Time |
|--------|-----------|-------------|-----------|---------------------|------|
| *D. melanogaster* (Release 6) | 11.8% | 92.71% | 96.44% | 94/142 / 85/130 | 86 min |
| *C. elegans* (WBcel235) | 0.6% | 23.74% | 53.85% | 9/23 / 8/18 | 7 min |
| Zebrafish (GRCz13ab) | 5.6% | 46.58% | 95.01% | 266/390 / 267/366 | 565 min |
| Zebrafish, without DIRS and Ngaro | 3.5% | 75.45% | 94.98% | | |
| Soybean (Glycine_max_v4.0) | 28.2% | 99.27% | 79.34% | 235/251 / 241/250 | 881 min |

- **Zebrafish: DIRS and Ngaro.** Dfam files DIRS and Ngaro elements under
  the LTR class.  They make up 38% of zebrafish's LTR-class bases, and the
  library finds none of them.  These elements integrate with a tyrosine
  recombinase: DIRS has inverted terminal repeats and Ngaro split direct
  repeats, so neither has a pair of direct LTRs, TG..CA ends or a TSD.
  ltr-forage looks for integrase-type LTR elements only.  Without DIRS and
  Ngaro, zebrafish scores like zebra finch.
- ***C. elegans*.** The genome has little LTR sequence, and almost all of it
  is old: families at about 17% divergence with few intact copies.  Half of
  the standard's LTR bases are one family, CER15-1, and nearly all of those
  are hits to a 500 bp tandem repeat inside its consensus.  The tandem-array
  mask hides that repeat on purpose.
- **Soybean.** The library masks 35% of the genome against the standard's
  28%.  91.5% of the false-positive bases are unannotated by the standard;
  only 5% come from families with no curated match, and the rest are extra
  copies of curated families.  RepBase's 251 soybean LTR models seem to
  miss copies more than the library invents them.
- **Maize.** The run did not finish.  ltr-forage kept 77,485 elements
  (155,000 LTR records), and step 12's all-vs-all search, at about 10,750
  hits per query, would have taken about 10 days for the LTRs alone.
  Genomes this rich in LTR elements need a different pre-clustering step.

---

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

ltr-forage hg38.sufr --ltrharvest --min-seed 30 --ltr-min 100 --ltr-max 1000 \
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

The current version, in LTRharvest mode, searched hg38 in 2 min with 4.4 GB
on 16 threads of a shared machine.

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
| ltr-forage predictions | 41,665 |
| Identical in all six coordinates and similarity | 41,264 |
| Same outer boundaries, different inner boundaries | 254 |
| Differ by 1-20 bp / 21-200 bp / over 200 bp | 2 / 61 / 70 |
| Only in LTRharvest / only in ltr-forage | 14 / 14 |

The "same outer boundaries" class is the duplicate tie-break described in
the introduction, and the rest follow from it: a different survivor of one
duplicate group changes which candidates the overlap rule later compares.
None of the 28 unique predictions contains N inside an LTR.

### Scoring against RepeatMasker

Neither tool's output is ground truth, so the scripts score both against the
intact elements in a RepeatMasker annotation of hg38: an LTR-class annotation
followed on the same strand by an internal region and then the nearest
annotation of the same family, with both LTRs 100-1000 bp long and their 5'
ends 1000-15000 bp apart, the window the tools were run with.  The annotation
used here has 7,170 such elements.  A prediction hits one when each of its LTRs overlaps the
corresponding annotated LTR by at least half the length of the shorter.

```bash
eval/rm_truth_recall.py hg38.fa.out hg38.fa.fai hg38_gt.out hg38_ltr-forage.out \
    hg38cmp.a_only.tsv hg38cmp.b_only.tsv
```

Recall by the larger of the two LTRs' divergence from consensus:

| Max LTR divergence | Elements | LTRharvest | ltr-forage |
|---|---|---|---|
| 0-5% | 585 | 90.9% | 90.9% |
| 5-10% | 1,323 | 63.6% | 63.6% |
| 10-15% | 1,785 | 21.2% | 21.2% |
| 15-20% | 1,301 | 7.6% | 7.6% |
| 20% and up | 2,176 | 2.3% | 2.3% |

Both tools hit the same 1,902 elements, and none of the 28 predictions
unique to one tool hits any of them.  The differences between the two
outputs fall entirely on predictions that do not correspond to intact
annotated elements.

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

5. Laslett D, Canback B. ARAGORN, a program to detect tRNA genes and tmRNA
   genes in nucleotide sequences. *Nucleic Acids Res.* 2004;**32**(1):11–16.
   DOI:[10.1093/nar/gkh152](https://doi.org/10.1093/nar/gkh152)

6. Myers G. FasTAN: a fast tandem repeat finder.
   [github.com/thegenemyers/FASTAN](https://github.com/thegenemyers/FASTAN)

7. Myers G. Efficient local alignment discovery amongst noisy long reads.
   In: *Algorithms in Bioinformatics (WABI 2014)*, LNCS 8701:52–67. Springer; 2014.

8. Ou S, Jiang N. LTR_retriever: a highly accurate and sensitive program for
   identification of long terminal repeat retrotransposons. *Plant Physiol.*
   2018;**176**(2):1410–1422.
   DOI:[10.1104/pp.17.01310](https://doi.org/10.1104/pp.17.01310)

---

## License

ltr-forage is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.  See [LICENSE](LICENSE).

`src/trna.rs` is a Rust port of the tRNA gene search in ARAGORN 1.2.41
(Copyright (C) 2003-present Dean Laslett; Laslett and Canback 2004, Nucleic
Acids Res. 32:11-16), which is distributed under GPL v3 or later; ltr-forage is
licensed under GPL v3 or later to match.

`src/fastan.rs` is a port of FasTAN and of daligner's local aligner
(Copyright (c) 2020, Dr. Eugene W. Myers), distributed under a BSD-style
license whose notice the file reproduces.

`data/TEfam.hmm` is LTR_retriever's `database/TEfam.hmm` (Ou and Jiang
2018), distributed under GPL v3.

The LTRharvest algorithm is implemented based on the published description and
the open-source GenomeTools reference implementation (ISC license).  The `sufr`
library is used under its BSD-3-Clause license.
