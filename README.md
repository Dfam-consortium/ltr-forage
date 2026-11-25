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

Install a release or build from source, with sufr, nail, MMseqs2,
rmblastn and te-composer findable (see [Installation](#installation) and
[Finding the programs](#finding-the-programs)).  To try it on a small
genome first, run the
[example on Arabidopsis chromosome 4](#example-arabidopsis-chromosome-4).

### 1. Index the genome (once)

```bash
sufr create --dna -o genome.sufr genome.fa
```

Do not pass `--ignore-softmask` on a soft-masked assembly (see
[Create a sufr index](#create-a-sufr-index)).  sufr writes its partition files to `$TMPDIR` (`/tmp` when
unset) and needs up to the size of the finished index there, about 9 bytes
per genome base (see [Create a sufr index](#create-a-sufr-index)).

### 2. Detect LTR retrotransposon copies

```bash
ltr-forage genome.sufr --library elements.fa --fastan-out tandem.bed -t 24 \
    --log ltr-forage.log > cand.tsv
```

The defaults are the settings the library driver uses: tandem arrays
masked, LTR pairs re-aligned in a 7 kb window, calls whose final ends agree
within 20 bp merged, protein domains searched with `data/TEdomains.hmm`
(122 Pfam profiles), overlapping calls resolved by evidence,
and the per-candidate table on stdout.  The domain search writes about 2 GB
of scratch files for a primate genome to `$TMPDIR` or `/tmp`;
`--temp-dir DIR` puts them elsewhere.  nail and MMseqs2 come from `--nail`
and `--mmseqs`, else the `NAIL` and `MMSEQS` environment variables, else
`PATH`.  `ltr-forage` checks its programs, the profile file and the output
directories before it loads the index, and reports every problem at once.
`-h` lists the options most runs need; `--help-all` adds the development and
tuning ones.

| File | Contents |
|---|---|
| `cand.tsv` | every candidate, kept or not, with the evidence behind the call: boundaries, LTR identity, TSD, TG..CA, PBS and its tRNA, protein domains, LTR self-repeat; `keep` is 1 for a kept element and `evidence` says which rule kept or rejected it (see [Per-candidate table](#per-candidate-table)) |
| `elements.fa` | each kept element's 5′ LTR, internal region and 3′ LTR, on the genome's plus strand; `strand=` in the header is the element's coding strand (`+`, `-`, `?`) |
| `tandem.bed` | the masked tandem arrays (name, start, end, period, 1000 × identity) |
| `ltr-forage.log` | `--log`: the screen's lines and every step's own lines, numbered as in the figure and in [How ltr-forage finds elements](#how-ltr-forage-finds-elements) |

The screen (stderr) shows each stage as it starts and a summary at the end;
on Arabidopsis chromosome 4 (18.6 Mb, 8 threads):

```
[14:23:01] ltr-forage ath_chr4.sufr --library elements.fa --fastan-out tandem.bed -t 8 --log ltr-forage.log
[14:23:01] Index and Tandem-Array Mask (stages 1-2)
[14:23:01]   1 sequence(s), 18.6 Mb
[14:23:09]   2,729 tandem array(s) masked (0.8 Mb)
[14:23:09] Candidate Detection (stages 3-7)
[14:23:27]   823 candidate element(s) from 3,470 seed pair(s)
[14:23:27] Evidence (stage 8)
[14:23:48]   230 with a protein domain hit (154 LTR order), 34 with a PBS, 82 with a TSD, 30 with TG..CA ends
[14:23:48] Decision and Output (stages 9-11)
[14:23:48] Results
[14:23:48]   116 LTR retrotransposon element(s) kept of 823 candidate(s)
[14:23:48]   rejected: 76 non-LTR protein hits, 54 LTR self-repeat, 9 domain inside an LTR, 0 poly-A-tailed LTRs, 35 overlapping a call with more evidence, 533 too little evidence
[14:23:48]   output: candidate table on stdout, elements.fa, tandem.bed
[14:23:48] done: 0.8 min
```

`-v` adds every step's own lines to the screen, and `--log FILE` always has
them.  On rice (MSU7, 373 Mb) with 16 threads they read (each line also
starts with the time):

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

`--format gt` writes the kept elements in LTRharvest's format instead of the
table.  The evidence for each copy is in `cand.tsv`; there is no per-copy HTML report
yet (the evidence page below is per family and needs the library run).

### 3. Build a family library, with the evidence page

```bash
pipeline/ltr_library.py genome.fa -o out -t 24 --evidence-page
```

This runs step 2 itself (with the same settings) and then clusters the
kept elements into families.  The screen shows each stage as it starts and
a summary at the end; on Arabidopsis chromosome 4 (18.6 Mb, 8 threads):

```
[13:47:58] ltr_library.py ath_chr4.fa -o out -t 8 --evidence-page
[13:47:58] Sequence & Suffix Array Preparation (stages 0-1)
[13:47:58]   genome 18.6 Mb non-N; depth sample 100.0% of the genome; floors 5 / 50
[13:48:01] LTR Detection: ltr-forage (stages 1-11)
[13:48:35]   117 LTR retrotransposon copies kept of 823 candidates (cand.tsv)
[13:48:35] Library Building
[13:48:35]   stages 12-13: clustering and consensus building
[13:48:41]   stages 14-22: artifact cleanup
[13:51:37]   stage    23: families and seed alignments
[13:51:41]   stage    24: evidence page
[13:51:43] Results
[13:51:43]   families: 48 (48 LTR and 48 internal models)
[13:51:43]   evidence: Well supported 23, Weak evidence 21, Possible artifact 4
[13:51:43]   output in out: library.fa, library.stk, families.tsv, cand.tsv, evidence.html, run.log
[13:51:43] done: 3.7 min
```

| File in `out/` | Contents |
|---|---|
| `library.fa` | the family library: one LTR and one internal model per family |
| `library.stk` | a seed alignment per model |
| `families.tsv` | each model's family, origin element, class, orientation (what set its strand; `flipped` when it was reverse-complemented to it) and estimated share of the genome |
| `cand.tsv` | ltr-forage's per-candidate table: the detected copies and the evidence behind each call (see [Per-candidate table](#per-candidate-table)) |
| `evidence.html` | with `--evidence-page`: one dossier per family |
| `run.log` | every step (numbered on from ltr-forage's, 12-24) with its time and results, ltr-forage's own step lines, and the programs, versions and paths used |

`-v` prints all of `run.log` on the screen as the run goes.  The steps work
in `out/work/`, which a successful run removes; `--keep-work` keeps it (each
step's files, their command logs in `out/work/logs/`, and the sufr index the
run built), and `--from STEP` needs it.  A run that stops keeps
`out/work/`, and the same command resumes it.

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
library. To add the page to a finished run made with `--keep-work`, run the
same command with `--evidence-page --from evidence`. With a curated annotation of the genome,
`pipeline/evidence_page.py` can also select the families that mask most
outside it (see [Evidence page](#evidence-page)).

### LTRharvest's search only

```bash
ltr-forage genome.sufr --ltrharvest -t 24 > ltrharvest.out
```

This is LTRharvest's default search and output format (see
[LTRharvest mode](#ltrharvest-mode---ltrharvest)).

---

## Example: Arabidopsis chromosome 4

`examples/ath_chr4.fa.gz` is chromosome 4 of *Arabidopsis thaliana*
(TAIR10.1, NCBI RefSeq NC_003075.7): 18.6 Mb, with a pericentromere full of
Gypsy and Copia elements.  The whole example takes about five minutes on 8
threads.  The commands below run in an empty directory, with `ltr-forage`,
`ltr_library.py` and the other programs findable (from a source checkout,
use `target/release/ltr-forage` and `pipeline/ltr_library.py`, or set
`LTRFORAGE`):

```bash
gzip -dc /path/to/ltr-forage/examples/ath_chr4.fa.gz > ath_chr4.fa
sufr create --dna -o ath_chr4.sufr ath_chr4.fa                  # 2 s, 167 MB
```

**Copies only.**  ltr-forage alone writes the candidate table, the kept
elements and the masked tandem arrays (under a minute):

```bash
ltr-forage ath_chr4.sufr -t 8 --library elements.fa --fastan-out tandem.bed \
    --log ltr-forage.log > cand.tsv
```

```
[15:03:21] ltr-forage ath_chr4.sufr -t 8 --library elements.fa --fastan-out tandem.bed --log ltr-forage.log
[15:03:21] Index and Tandem-Array Mask (stages 1-2)
[15:03:21]   1 sequence(s), 18.6 Mb
[15:03:28]   2,729 tandem array(s) masked (0.8 Mb)
[15:03:28] Candidate Detection (stages 3-7)
[15:03:46]   823 candidate element(s) from 3,470 seed pair(s)
[15:03:46] Evidence (stage 8)
[15:04:08]   230 with a protein domain hit (154 LTR order), 34 with a PBS, 82 with a TSD, 30 with TG..CA ends
[15:04:08] Decision and Output (stages 9-11)
[15:04:08] Results
[15:04:08]   116 LTR retrotransposon element(s) kept of 823 candidate(s)
[15:04:08]   rejected: 76 non-LTR protein hits, 54 LTR self-repeat, 9 domain inside an LTR, 0 poly-A-tailed LTRs, 35 overlapping a call with more evidence, 533 too little evidence
[15:04:08]   output: candidate table on stdout, elements.fa, tandem.bed
[15:04:08] done: 0.8 min
```

`cand.tsv` has one row per candidate (all 823), with `keep` and the
`evidence` that decided it; `elements.fa` has each kept element as its 5′
LTR, internal region and 3′ LTR, with the class, evidence and coding strand
in the header.  To tally why candidates were kept or rejected:

```bash
awk -F'\t' 'NR==1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
            {print $c["keep"], $c["evidence"]}' cand.tsv | sort | uniq -c | sort -nr | head
```

**A family library.**  The driver runs ltr-forage itself (it finds
`ath_chr4.sufr` beside the genome) and builds the families:

```bash
ltr_library.py ath_chr4.fa -o out -t 8 --evidence-page
```

```
[15:04:08] ltr_library.py ath_chr4.fa -o out -t 8 --evidence-page
[15:04:08] Sequence & Suffix Array Preparation (stages 0-1)
[15:04:08]   genome 18.6 Mb non-N; depth sample 100.0% of the genome; floors 5 / 50
[15:04:08] LTR Detection: ltr-forage (stages 1-11)
[15:04:54]   116 LTR retrotransposon copies kept of 823 candidates (cand.tsv)
[15:04:54] Library Building
[15:04:54]   stages 12-13: clustering and consensus building
[15:05:00]   stages 14-22: artifact cleanup
[15:07:53]   stage    23: families and seed alignments
[15:07:57]   stage    24: evidence page
[15:07:58] Results
[15:07:58]   families: 47 (47 LTR and 47 internal models)
[15:07:58]   evidence: Well supported 22, Weak evidence 21, Possible artifact 4
[15:07:58]   output in out: library.fa, library.stk, families.tsv, cand.tsv, evidence.html, run.log
[15:07:58] done: 3.8 min
```

`out/library.fa` is the library, one LTR and one internal model per family,
each in its element's coding strand where the domains or PBS give it;
`out/library.stk` holds a seed alignment per model; `out/families.tsv`
describes each model; `out/evidence.html` has a dossier per family, with a
verdict (open it in a browser); `out/run.log` has every step's numbers and
the programs and versions used.  The intermediate files are removed after a
successful run; `--keep-work` keeps them in `out/work/`, and `-v` prints
every step as it runs.  `ltr_library.py -h` lists the options most runs
need, `--help-all` all of them.

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

The default protein profiles, `data/TEdomains.hmm`, are 122 Pfam profiles
of transposable-element proteins: 116 from Pfam 39.0 and six taken from
LTR_retriever's `database/TEfam.hmm` (Ou and Jiang 2018).  The choice of
profiles starts from LTR_retriever's list.  `data/TEfam.hmm` is
LTR_retriever's file itself (102 Pfam profiles), distributed under GPL v3.

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
end; tabs or spaces) to the mask.  A line that does not parse stops the
run, since skipping it would leave that interval unmasked.
`--no-fastan` turns the mask off.

**Known defect.** A repeat that runs into the mask on its right is cut at
the mask.  A repeat that starts under the mask in either copy is dropped,
even when `--min-seed` or more unmasked matching bases follow: the
remainder is not left-maximal, and the test for that does not know about
the mask.  A pair of identical LTRs is one long match, so a masked array
over the first bases of either LTR removes its only seed.  A patched seed
search (candidate step only, no library scored) kept 4 more elements on
Arabidopsis, 36 on rice, 12 on zebra finch, 83 on lotus, 7 on gorilla and
81 on mouse, 0.3% to 0.9% of the kept elements; 201 of the 223 have
curated LTR annotation on both LTRs.  Not fixed yet.

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

#### Two-pass window (`--reextend-first`)

By default each candidate is aligned with a 2,500 bp window first
(`--reextend-first`), which costs a fraction of the full window's
alignment.  The result stands unless it ends within
`--reextend-first-margin` (500) bases of a window edge, or no alignment
qualifies.  Of those, a candidate with an LTR at least `--tandem-max`
self-repeat is a satellite array running to the edge: it keeps the short
window's LTRs, and step 8a vetoes it on them without aligning a 7 kb
array against itself.  Every other one is aligned again with the full
`--reextend-ltr-max` window, so long LTRs keep their full length.

The two passes are not identical to one full-window pass: a short window
occasionally holds a different best alignment, even well inside it.

| Genome | Candidate detection (8 threads) | Kept elements | Library sens. / prec. |
|---|---|---|---|
| Rice | 9.8 -> 6.9 min | 4,339 -> 4,340 (4 lost, 5 gained) | 93.44 / 98.43 -> 93.43 / 98.43 |
| Gorilla | 103 -> 45 min | 2,142 -> 2,144 (none lost) | 36.78 / 97.96 either way |

`--reextend-first 0` gives the single full-window pass.

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
[nail](https://github.com/TravisWheelerLab/nail): `data/TEdomains.hmm` by
default, or any HMMER3 amino-acid profile file given with `--domain-hmm`.
`--no-domains` skips the search.  On hg38, with LTR_retriever's
`TEfam.hmm`, a hit of LTR order marks 68% of the RepeatMasker-annotated
elements and 0.5% of the other candidates.

`data/TEdomains.hmm` holds 122 Pfam profiles (67 LTR, 55 not LTR).  Each
carries its RepeatMasker class in its name (`rve#LTR`,
`RVT_2#LTR/Copia`), and polyprotein domains are marked `[role=core]`.
It is LTR_retriever's list rebuilt from Pfam 39.0, with retrotransposon
core domains added, Atrophin-1 removed, and six profiles kept as their 2015
builds because the Pfam 39.0 rebuilds find far fewer elements.
`data/domains/README.md` gives the list, the evidence for each class and
the build script.  `--domain-hmm data/TEfam.hmm` uses LTR_retriever's own
file (102 Pfam profiles), the default until 2026-11-17.

ltr-forage translates each element plus `--domain-flank` bases (default 50)
on each side in six frames, with stop codons as X, writes the frames as a
protein FASTA, and searches it with nail's full pipeline (MMseqs2 seeds, cloud
search, Forward, composition-bias correction), with a fixed database
size for the E-value (`--domain-z`, default 120,000: six frames of 20,000
candidates), with nail's seed filter at P = 0.01 instead of its default
1e-4 (twice the time, half as many weak hits lost against hmmsearch), and
with nail's progressive seeding never stopping early (`--prog-f 0`). Each LTR
is searched alone in the same run (see the coding-LTR veto below). Hits with
an E-value at or below `--domain-evalue` (default 0.05) count. `nail` and
`mmseqs` are given with `--nail` and `--mmseqs`, or the `NAIL` and `MMSEQS`
environment variables, else found on PATH; the frames
and nail's databases go to a subdirectory of `--temp-dir` (default:
`$TMPDIR`, else `/tmp`), removed afterwards.

The database size is a constant so that a hit's E-value, and with it
whether an element is kept, does not depend on how many other candidates
the genome has. Until 2026-11-18 it was six frames per candidate of the
run (`--domain-z 0` gives that rule): gorilla's 157,206 candidates made
every E-value 7.9 times what it is now, Arabidopsis' 5,573 made it 0.28
times. The constant sits in the range the candidate counts of the tuning
genomes had when the thresholds were set (13,000 to 29,000). With the
constant and `--prog-f 0` (2026-11-18/19, sensitivity / precision, before
then after): Arabidopsis 91.09 / 95.50 to 91.05 / 95.46, zebra finch
76.13 / 86.28 to 76.08 / 86.28, rice 93.57 / 98.42 to 93.56 / 98.42, lotus
92.64 / 94.32 to 92.63 / 94.20.  Lotus ran with rmblastn 3.1.0, the other
three with 3.0.12.

**Known limitation.** On a large target set nail does not report some
strong hits that it reports when the same sequences are searched in a
smaller file.  On gorilla (2.7 M frames, `TEfam.hmm`) 51 of about 11,500
candidate and profile hits were missing, the best at E = 2.1e-14; 50 of
them came back when their candidates were searched alone.  Two causes:
nail's progressive seeding stops early for a profile with many targets
(ltr-forage now runs it with `--prog-f 0`, which brings back 21), and
MMseqs2 gives no seed for the rest.
The nail authors know of it; the MMseqs2 part stays until nail drops or
integrates MMseqs2.  Splitting the search is not a remedy (it loses more
hits than it brings back).  In that run 15 gorilla and 5 mouse elements
with LTR-annotated LTRs lost their evidence this way, 12 and 4 of them
kept again with `TEdomains.hmm` through other hits.

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
The lists come from the classes in the profile names; a file without
classes, such as `TEfam.hmm`, gets `annotate_TE.pl`'s own lists.
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
starts with TGG. The gene ends at the discriminator base: where the genome
happens to continue with C, CC or CCA, ARAGORN reports the gene that much
longer, and those bases are left out before the CCA is added (counting them
shifted the probe by one to three bases). The test is the probe's first 12 nt, the tRNA's
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
per distinct tRNA gene sequence (up to the discriminator base), sorted by the number of kept elements it
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
`--fastan-out`, `--format tsv` and `--overlaps evidence` are refused, and so
is any option that only the skipped steps read (the re-extension, domain,
decision and PBS options), which would otherwise be ignored without a word.
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
| Seed enumeration | Bottom-up ESA traversal (C), single-threaded | Same algorithm, streaming the sufr SA/LCP from disk, parallel over LCP-delimited pieces (`-t`) |
| Seed extension | Single-threaded | Parallelised with Rayon (`-t`) |
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

### Installing a release

Each [release](https://github.com/Dfam-consortium/ltr-forage/releases) has an
archive for Linux (x86_64) and macOS (Apple Silicon):

```bash
tar -xzf ltr-forage-0.1.0-x86_64-unknown-linux-gnu.tar.gz
export PATH=$PWD/ltr-forage-0.1.0/bin:$PATH          # tcsh: set path = ($cwd/ltr-forage-0.1.0/bin $path)
ltr-forage --version
```

The archive holds `bin/ltr-forage`, `bin/ltr_library.py` (a link to
`pipeline/ltr_library.py`), the driver's step scripts in `pipeline/`, the
protein profiles and scoring matrices in `data/`, and the example in
`examples/`.  The tree can be moved: `ltr-forage` finds
`data/TEdomains.hmm` beside its `bin/`, and the driver finds its scripts and
`data/` beside itself.  The other programs below are not included.

### Prerequisites

- The `sufr` tool for building `.sufr` index files
  ([TravisWheelerLab/sufr](https://github.com/TravisWheelerLab/sufr)),
  version 0.7.12.  sufr 0.8.0 writes a new index format that this release
  cannot read.
- For the protein-domain search (step 8b): `nail` ≥ 0.7.1
  ([TravisWheelerLab/nail](https://github.com/TravisWheelerLab/nail)) and
  MMseqs2 (`mmseqs`), which nail seeds with; not needed with `--no-domains`
  or `--ltrharvest`.  The protein profiles ship in `data/TEdomains.hmm`.
- To build from source: a recent stable Rust toolchain
  ([rustup](https://rustup.rs) recommended).

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
git clone https://github.com/Dfam-consortium/ltr-forage.git
cd ltr-forage
make                               # cargo build --release; the binary is target/release/ltr-forage
make install PREFIX=$HOME/opt      # optional: the release layout in $HOME/opt/ltr-forage-<version>/
```

The build fetches libsufr and rmblast-lib from GitHub (pinned in
`Cargo.toml`).  A binary built in the checkout finds `data/TEdomains.hmm`
two levels up, and `--version` names the commit it was built from, with
`-dirty` when the source had uncommitted changes.

### Finding the programs

Neither program has built-in install locations.  `ltr-forage` takes nail and
MMseqs2 from `--nail` and `--mmseqs`, else the `NAIL` and `MMSEQS`
environment variables, else `PATH`.  `ltr_library.py` finds each program the
same way:

| Program | Option | Environment variable |
|---|---|---|
| ltr-forage | `--ltr-forage` | `LTRFORAGE` |
| sufr (only when the driver builds the index) | `--sufr-bin` | `SUFR` |
| rmblastn (the Rust port) | `--rmblastn` | `RMBLASTN` |
| te-composer | `--te-composer` | `TE_COMPOSER` |
| nail | `--nail` | `NAIL` |
| MMseqs2 | `--mmseqs` | `MMSEQS` |

It skips a program older than it needs and says which one it skipped when
none qualifies.  Each library run writes the path and version of every
program it used to `out/run.log`.  For example, with the programs in
`/opt/tools`:

```bash
export NAIL=/opt/tools/nail/nail MMSEQS=/opt/tools/mmseqs/bin/mmseqs        # bash
setenv NAIL /opt/tools/nail/nail; setenv MMSEQS /opt/tools/mmseqs/bin/mmseqs  # tcsh
```

### Create a sufr index

Before running `ltr-forage`, build a `.sufr` index of your genome with the
`sufr` tool.  The `--dna` flag is required to restrict the alphabet to ACGT:

```bash
sufr create --dna -o genome.sufr genome.fa
```

sufr 0.7.12 takes its thread count before the subcommand
(`sufr -t 8 create ...`).

`sufr create --dna` stores the sequence upper-cased, so a soft-masked FASTA
(lowercase runs, as NCBI distributes them) needs no preprocessing; the index
and everything `ltr-forage` reads from it are case-free. Query strings given
to `sufr count` are not normalised, so pass them in upper case.

sufr writes its partition files to `$TMPDIR`, or `/tmp` when it is unset,
and needs up to the size of the finished index there (about 9 bytes per
genome base: 0.58 GB peak for the 100-Mb *C. elegans* genome, whose index is
0.90 GB).  Point `TMPDIR` at a disk with room.  sufr 0.7.12 does not check
these writes: when the disk fills, it builds an index of the expected size
with corrupt LCP values and reports no error, and ltr-forage then stalls in
step 3 (a 3.4-Gb genome on a machine whose `/tmp` had 20 GB free).

---

## Usage

```
ltr-forage [OPTIONS] <SUFR>
ltr-forage --ltrharvest [LTRharvest options] <SUFR>
ltr-forage --trna-scan <FASTA> [-t N] [--trna-intron-max BP]

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
      --mask-bed <FILE>     Also mask these intervals (BED: name, 0-based start, end);
                            a line that does not parse is an error

Re-extension and near duplicates (steps 6, 7):
      --reextend-ltr-max <BP>   Longest LTR the alignment window reaches beyond each X-drop
                            end [default: 7000]
      --reextend-vicinity <BP>  Extra window beyond that on the outer side of each LTR [default: 100]
      --reextend-first <BP>     Re-extend with this window first, the full one only near its
                            edge; 0 disables [default: 2500]
      --reextend-first-margin <BP>  A first-pass end within BP of the window edge counts as
                            reaching it [default: 500]
      --trim-window <N>     Terminal alignment columns that must be clean for an end to stand [default: 6]
      --trim-matches <N>    Matches required inside that terminal window [default: 6]
      --dedup-tol <BP>      Merge candidates whose final outer ends both lie within BP of a
                            better one's; 0 disables [default: 20]
      --internal-stats      Count 20-mers of each internal region against the index and
                            report coverage and the longest single-copy stretch

Protein domains (step 8b):
      --domain-hmm <FILE>   HMMER3 file of protein profiles [default: data/TEdomains.hmm]
      --no-domains          Skip the domain search
      --domain-evalue <E>   E-value threshold for a domain hit [default: 0.05]
      --domain-z <N>        Database size for the domain E-values (nail's -Z); 0: six frames per candidate of this run [default: 120000]
      --domain-flank <BP>   Bases added on each side of the element before translation [default: 50]
      --nail <PATH>         nail binary (a path, or a name on PATH) [env: NAIL] [default: nail]
      --mmseqs <PATH>       MMseqs2 binary nail seeds with [env: MMSEQS] [default: mmseqs]
      --temp-dir <DIR>      Directory for scratch files [default: $TMPDIR, else /tmp]
      --domain-search-vetoed  Search candidates the self-alignment check rejected too
      --keep-coding-ltrs    Keep candidates with a polyprotein domain hit inside an LTR

Primer binding site (step 8c):
      --no-pbs-trna         Skip the genome tRNA scan behind the PBS
      --pbs-domain-strand   Search for the PBS only on the strand the protein domains give (default: both)
      --trna-out <FILE>     Write the genome's tRNA library (FASTA) and every call (FILE.tsv)
      --trna-intron-max <BP>  Longest tRNA intron searched for, 0 to 3000 [default: 3000]
      --trna-scan <FASTA>   Only scan FASTA for tRNA genes, print the calls and exit

Decision and overlaps (steps 9, 10):
      --keep <LEVEL>        Evidence a candidate needs: domain | standard | loose | none [default: standard]
      --tandem-max <FRAC>   Reject a candidate when either LTR aligns to itself off the diagonal
                            over at least this fraction of its length [default: 0.99]
      --internal-gap-max <FRAC>  Also reject candidates whose internal region has a single-copy
                            stretch longer than this fraction of it (needs --internal-stats)
      --overlaps <STRATEGY> evidence | best | no | all [default: evidence; best with --ltrharvest]

Output (step 11):
      --format <FORMAT>     tsv | gt [default: tsv; gt with --ltrharvest]; gt lists the kept
                            elements (every prediction with --ltrharvest)
      --library <FILE>      Write the kept elements as a FASTA library (left LTR, internal
                            region, right LTR per element)
  -t, --threads <N>         Threads (0 = all CPUs; -j also works) [default: 0]
  -v, --verbose             Print every step's own lines (the screen otherwise shows stages and results)
      --log <FILE>          Also write the whole log to FILE
  -h, --help                Print help (the options most runs need)
      --help-all            Print help with the development and tuning options too

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
ltr-forage hg38.sufr --ltrharvest -t 32 > hg38_ltrs.out
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
    -t 32 > hg38_ltrs.out
```

---

## Output Format

The default output is the per-candidate table below.  `--format gt` (the
default with `--ltrharvest`) is compatible with GenomeTools LTRharvest.
LTRharvest's table has no keep column, so outside `--ltrharvest` it lists
the kept elements only, with their final (re-extended) coordinates.
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
| `pbs_offset` | Start of the primer binding site relative to the first base after the 5′ LTR (−5 to +15), from the 12-nt tRNA match described under step 8c; `.` when none or with `--no-pbs-trna` (it was −1 before 2026-11-15, which is also a valid offset) |
| `n_count` | Ns inside the element |
| `xdrop_tsd_len` | Always 0: the TSD search at the X-drop ends runs only in LTRharvest mode, which writes no table (the column is kept so the others do not move) |
| `int_tiles`, `int_multi_frac`, `int_gap_frac`, `int_gap_bp` | With `--internal-stats`: number of 20-mer tiles taken every 10 bp across the internal region (tiles containing N skipped), the fraction seen at least three times in the genome on either strand, and the longest run of single-copy tiles as a fraction of the tiles and in bases; `.` otherwise |
| `dom_order`, `dom_superfamily`, `dom_strand` | `LTR`, `notLTR`, `mixture`, `NA` or `none`; `Gypsy`, `Copia` or `unknown`; `+`, `-` or `?`, as described under step 8b; `.` with `--no-domains` or for a candidate the self-alignment check rejected |
| `dom_best_e`, `dom_hits` | Best E-value, and every hit as `accession:E-value:frame` (frames `aa1`..`aa3` forward, `rev_aa1`..`rev_aa3` on the reverse complement), comma separated; `.` when none |
| `ltr_dom_hits` | Hits to the LTR-order (polyprotein) profiles when each LTR is searched on its own, summed over the two LTRs; 0 when `--keep-coding-ltrs` skips the search |
| `lltr_self_cov`, `rltr_self_cov` | Fraction of each LTR covered by plus-strand alignments of the LTR against itself off the main diagonal; 1.00 for a satellite monomer or an Alu dimer, 0 for unique sequence |
| `keep`, `evidence` | 1 when the candidate is kept, and what decided it: `domain`, `tsdN`, `tgca` (comma separated), `veto:notLTR`, `veto:mixture`, `veto:tandemF`, `veto:ltrdomainN`, `veto:polyA`, `gapF` when `--internal-gap-max` rejected it, or `none`; `,overlap` is appended when step 10 dropped a kept candidate for an overlapping one with more evidence.  In the library driver's `cand.tsv`, step 11b adds `veto:polyAfamily` to the elements it vetoes and writes `rescued` in place of `overlap` on those it keeps again |
| `pbs_trna`, `pbs_mm12` | The tRNA type(s) whose 3′-terminal 12 nt match the PBS, e.g. `Met(cat)` (several joined by `/` when they share the 12 nt), and the mismatches (0 or 1); `.` when none |
| `pbs_edist18` | Best edit distance (0 to 2) of a full 18-nt tRNA probe over the same window; `.` when none |
| `pbs_strand` | `+` when the PBS was found after the left LTR, `-` when before the right LTR on the reverse strand; `.` when none |
| `tsd_mm` | 1 when the TSD is the one-mismatch fallback, else 0. It is enriched about 20-fold over shifted flanks only for TG..CA elements (lotus), so it counts toward `--keep standard`/`loose` only alongside TG..CA; the evidence tag reads `tsd5mm1` |
| `strand` | The element's coding strand: `dom_strand` when it is `+` or `-`, else `pbs_strand`; `?` when neither is known or the two disagree.  Coordinates and `elements.fa` stay on the genome's plus strand; `elements.fa` headers carry it as `strand=` |

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

`-h` lists the options a user sets; `--help-all` adds the development and
tuning options (depth sample, clustering, cleanup rules, ltr-forage
pass-through) and the list of steps.  Every option works either way.

It finds `ltr-forage`, `sufr`, `rmblastn`, `te-composer`, `nail` and `mmseqs`
from the `--ltr-forage`, `--sufr-bin`, `--rmblastn`, `--te-composer`, `--nail`
and `--mmseqs` options, then the `LTRFORAGE`, `SUFR`, `RMBLASTN`,
`TE_COMPOSER`, `NAIL` and `MMSEQS` environment variables, then `PATH`; it
has no built-in install locations.  It skips a binary older than it needs
(rmblastn 3.0.9, te-composer 0.6.6, nail 0.7.1) and names it when no
other qualifies. The protein
profiles default to the shipped `data/TEdomains.hmm` (`--domain-hmm`). A sufr
index beside the genome (`genome.sufr` for `genome.fa`) is used when it
exists and is not older than the genome; `--sufr` names another. Each run
appends the command and every program's path and version to
`out/run.log`, with each step's time and results.  Every step's temporary files go to
`out/work/tmp` (`--tmpdir` names another place), not the system `/tmp`; the
index step needs up to the size of the finished index there.  The genome must
be an uncompressed FASTA file whose sequence names are unique and hold no
`#` (library records use it to separate name and class). The step scripts live beside the driver in `pipeline/`.

### Steps

The numbers continue ltr-forage's and appear in `out/run.log` (on the
screen with `-v`) and the figure.

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

11b. **Poly-A family veto** (`pipeline/polya_family.py`; `--polya-family 0`
    turns it off).  ltr-forage vetoes a candidate whose two LTRs both end
    in a poly-A tail: a pair of SINE or LINE copies.  The same repeats also
    flank old ERV fragments, and such a pair is kept on the domain hit
    between them.  This step compares both LTRs of every kept element with
    the left LTRs of 2,000 poly-A-vetoed candidates (a seeded sample,
    `--seed`) and vetoes the element when at least half of each LTR is
    covered by at least `--polya-family` (10) of them.
    - **It runs only on genomes with many such pairs.**  The pool of
      poly-A-vetoed candidates, not counting those with an LTR domain hit of
      their own, must hold at least `--polya-pool` (2,000).  Gorilla has
      about 43,000, mouse about 17,400 and hg38 about 47,600; Arabidopsis, zebra finch, rice
      and lotus have 1 to 42, so the step does nothing there.  In a small
      pool many of the candidates are real elements the poly-A rule vetoed
      wrongly (soybean 32 of 80, maize 18 of 24), and kept members of their
      families match them.
    - **Not measured: a pool between 80 and 17,000.**  No genome we ran
      falls there, and the limit of 2,000 sits inside that range.  What the
      step does on a genome just above the limit, or what such a genome
      loses just below it, is unknown.  The step's line in `run.log`
      gives the pool size.
    - On mouse and gorilla no kept element with LTR annotation on both
      LTRs matched more than 5 of the 2,000; SINE-flanked ones matched
      1,400 to 1,950 and LINE-flanked ones 10 to 380.  The step vetoes 148
      to 298 kept elements on gorilla and 758 to 1,154 on mouse, depending
      on the profile file.
    - The step was developed on mouse and gorilla.  On hg38 (candidate
      step only, no library) it vetoes 348 of 2,903 kept elements: 302
      flanked by SINEs, 43 by LINEs and 3 by Alu and L1 pieces, none with
      LTR annotation on both LTRs.  No LTR-annotated element matched more
      than 1 of the 2,000, and the rescue kept 38 again, 37 of them
      LTR-annotated.  Human and gorilla share their SINE and LINE
      families; a mammal with other SINE families has not been run.
    - **Overlap rescue** (`--no-polya-rescue` turns it off).  ltr-forage
      resolves overlaps before this step.  A candidate that lost the
      overlap choice only to elements vetoed here is kept again when no
      part of either of its LTRs is covered by 10 of the sampled LTRs: 8 to
      59 elements a run, 109 of 122 with LTR annotation on both LTRs.
    - The vetoed rows of `cand.tsv` get `keep` 0 and `veto:polyAfamily`;
      a rescued row reads `rescued` in place of `overlap`.
    - The sample is seeded and the search does not depend on the thread
      count, so a rerun gives the same library.  Another seed moves the
      vetoed set by a few elements near the count of 10.
12. **Pre-clustering.**  The raw library is clustered all-vs-all with
    rmblastn, LTRs against LTRs and internals against internals (the
    clustering reads only same-kind hits of 80% identity or more, and drops
    the rest as they are written).  A probe of 150 records of each kind first
    projects the hits of that search.  Above `--greedy-above` million (100)
    the clustering is greedy instead: batches of `--greedy-batch` (2,000)
    elements are clustered all-vs-all, and the remaining elements are searched
    against at most `--greedy-reps` (100) members of each cluster, by the same
    linking rules.  None of the six tuning genomes reaches 100 million (mouse
    projects 45); maize projects 2,365.  The searches seed with 16-base
    words (`--word`); 12 gave the same scores within 0.05 points on rice,
    lotus and mouse, and took this step 1.3-2x as long.
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
14. **Depth cleanup.**  Three rounds (`--depth-rounds`) of sample depth
    over the library, from an rmblastn search of the sample against it, one
    single-thread batch per thread.  Stretches far
    above their record's own depth are excised as nested passengers: 3x for
    an element, 10x for a consensus representative, whose family is
    high-copy itself.  High-copy LTRs are dropped, but never a consensus
    representative.  The search uses mask level 80, so a passenger that
    many records share splits its genome copies among them, and each record
    may see too little depth to cut it (step 22b catches these).  Only
    the first round is a search, run without the mask-level filter and
    keeping each hit's alignment.  A later round applies the previous
    round's cuts to those hits (`masklevel.py remap`): a hit on a dropped
    record goes; a hit on a cut record keeps the alignment columns outside
    the cuts, rescored as rmblastn scores them, with each cut end trimmed
    or extended across the junction; two pieces either side of a cut are
    joined when they are collinear in the genome (a copy without the
    passenger).  rmblastn's cross-record mask-level filter is then rebuilt
    over the moved hits.  This is close to searching the cut records again
    but not identical (see [Accuracy](#accuracy)); `--no-depth-remap`
    searches the records the previous round trimmed instead, which gives
    the same cuts as three full searches (identical libraries on rice,
    lotus and mouse).  Step 22b's second depth search reuses the first
    one's hits in that exact way.  The unfiltered hits (gigabytes) are
    deleted when steps 14 and 22b finish.
15. **Nested-LTR vote.**  A first element clustering excises LTRs nested in
    internals by vote.
16. **Element clustering.**  Clusters and one exemplar per cluster; each big
    cluster's members rejoin the cluster of their representative.  The
    clustering links elements that cover each other, then attaches each
    remaining single element to the cluster that covers it.  A consensus
    representative is never attached that way: it stands for a whole
    cluster, and attaching it would merge two clusters into one set of
    exemplars (until 2026-11-18 it was attached, 4 to 22 times a tuning
    genome).
17. **Flank filter.**  A cluster is dropped when most of its pairs, or
    every one of its copies, share both flanks (300 bp at 80% identity)
    with another copy: the copies are one duplicated locus, not independent
    insertions (`--no-flank` skips this).  Only pairs within a cluster
    count, so the flanks are searched in bins of whole clusters rather
    than all against all: the same calls on seven genomes, and maize's
    search falls from 689 to 48 minutes (44 GB of hits to 1 GB).
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
22b. **Model depth check.**  When step 11b ran, the models are first cut
    where at least `--polya-family` (10) of its sampled poly-A-terminated
    repeats cover them (`pipeline/polya_cut.py`), and a model with under 50
    bases left is dropped.  Until step 11b the elements those repeats flank
    told the nested-LTR vote (step 15) where an L1 or SINE lay inside a
    real element's internal region.  With them vetoed that evidence is
    gone, and the depth check cannot replace it where the insertion is most
    of the record: one mouse MTA internal of 3,594 bp, 2,807 of them L1,
    masked 184 Mb outside the curated LTRs.
    Then one more depth round on the finished models:
    a consensus model is cut where its depth reaches 10x its median, an
    exemplar at 5x (`--model-depth-ratio`); the cuts drop no model
    (`--no-model-depth` skips the step).  It removes passengers step 14
    left shared among many records, and passengers a consensus took back
    from copies its founder had lost them from: on mm10, L1 and B1/B2
    pieces in a few models masked 191 Mb outside the curated annotation.
    On Arabidopsis it cuts a complete 2 kb ATGP2N solo LTR (39x the
    model's depth) out of an ATCopia76 internal.  At 3x it also cut
    conserved stretches of old, 30%-diverged Gypsy families at 3–4x, with
    no gain in precision anywhere.  With `--early-tier`, the same search
    drops an internal model whose genome copies exceed 2,000 per element
    of its cluster (`--int-copies-max`; 0 turns it off).  Internals are
    never solo, so such a model is mostly another repeat: with
    `--early-tier`, an mm10 exemplar internal of L1 sequence had 5,999
    copies per element and masked 133 Mb outside the curated LTRs.
    Without `--early-tier` the check is off unless `--int-copies-max` is
    given: no model in the six genomes' default runs went above 707, so it
    changed none of them there.  The depth search is then repeated on the cut models, because a
    cut segment's copies go to whichever model matches them next.  An LTR
    model is dropped when its median depth rises at least 20-fold
    (`--depth-jump-fold`; 0 turns it off) and its element has neither a TSD
    nor TG..CA ends.  Depth alone would not do: a real zebra finch LTR nested
    in another family's internal rose 170-fold, and has both.  On mouse, two
    LTR models of L1 sequence rose 101- and 60-fold when an internal's L1
    segment was cut; no model in the other five genomes meets the rule.  A
    dropped model's copies move on in the same way, so after a drop the
    depth is taken once more without the dropped models (no search: their
    hits are removed and the mask-level filter rebuilt, 5 s on mouse), and
    a model this pass cuts differently from the second takes its cuts.  On
    mouse an internal's 87-bp L1 piece stood at depth 38 beside the two
    dropped models and at 7,447 without them; left in, it masked 5.2 Mb
    outside the curated LTRs (precision 89.51 instead of 92.46).  The
    second search adds 3 to 16 minutes.
23. **80–80 families**, LTR and internal models grouped separately (and
    searched separately), with Gypsy and Copia classes from the domain hits.
    Then each model is turned to its element's coding strand
    (`pipeline/orient_models.py`).  A model is in the genome strand of the
    element it was founded on, and that element's strand comes from
    ltr-forage's table: the protein domains' strand, else the PBS's.  A
    model of a minus-strand element is reverse-complemented.  Models whose
    element has neither, or whose domains and PBS disagree, stay as found.
    `families.tsv` records it (`orientation`, `flipped`).  On Arabidopsis
    chromosome 4: 42 of 94 records turned, 8 unknown; after turning, every
    domain hit of the known-orientation founding elements lies on the plus
    strand.
23b. **Seed alignments** (`library.stk`): each model's best
    `--seed-copies` (50) non-overlapping copies in the depth sample, or its
    origin element when the sample holds none, aligned to the model by
    `te-composer --fixed-cons`.  The model is not refined, so each
    alignment's `#=GC RF` line is the model itself.  Rows are named by
    genome range.  On rice: 865 alignments, median 23 rows, 51 s.  A
    te-composer without `--fixed-cons` (before 0.6.8) gives one-row
    alignments.

    **Large-record warning (report only).**  The same search is read as a
    masking run (each sample base goes to the best-scoring record over it),
    and each record's estimated share of the non-N genome goes into
    `families.tsv` (`est_genome_pct`).  A record at or above `--mask-warn`
    (2%) is named in the log line.  Nothing is dropped or changed.  The
    failure it is for is one record of L1 or SINE sequence called LTR that
    masks a large part of the genome: in three mouse runs with a rule
    switched off, one record each was estimated at 3.5% to 5.5% of the
    genome and made 60% to 90% of the run's false-positive bases.  On the
    six tuning libraries the largest record is 1.2% (lotus).  The 2% was set
    once from these nine libraries and is not tuned.  **A warning means
    "look at this record", not "this record is wrong"**: in a genome made
    largely of a few LTR families, real records pass 2% (maize: 3 of 2,268
    records, all one Gypsy family, up to 6.2%).
24. **Evidence page** (`--evidence-page`).

`--stage-a` skips steps 12–13 and depth-cleans and re-clusters every
element, as the driver did before.  It is slower on large libraries and no
more accurate.  `--min-cons 0` gives exemplar models only.  The steps
work in `out/work/`; each step's command lines and output go to
`out/work/logs/`, and its time and results to `out/run.log` (the screen
shows a summary, `-v` all of it).  A successful run removes `out/work/`
unless `--keep-work`.  A rerun skips finished steps and `--from STEP` reruns
from a step on; both need `out/work/` (a `--keep-work` run, or one that
stopped).  Each finished step records the options it read, and the outdir
records the genome (path, size and modification time): a rerun whose options
or genome differ from a finished step's stops and names the step to rerun
from, instead of mixing old and new results.  Programs and profile files
are not part of that record; `run.log` shows what each run used.  When
ltr-forage keeps no element the run ends with an empty `library.fa`,
`library.stk` and a header-only `families.tsv`.

### Output

| File | Contents |
|------|----------|
| `library.fa` | The family library: `fam<N>_v<k>_LTR` and `fam<N>_v<k>_INT` records, classified `#LTR/Gypsy`, `#LTR/Copia` or `#LTR` |
| `library.stk` | A Stockholm seed alignment per model, with the model as `#=GC RF` (step 23b) |
| `families.tsv` | Each model's family, origin element, kind, length, status and class; `orientation` (what set the model's strand: `domains+PBS`, `domains`, `PBS`, or `conflict` / `unknown` when it was left as found) and `flipped` (`yes` when it was reverse-complemented to the coding strand); `est_genome_pct` |
| `cand.tsv` | ltr-forage's per-candidate table, with every value behind its decision (see [Per-candidate table](#per-candidate-table)) |
| `run.log` | Each run's command, programs, versions and paths; every step's time and results; ltr-forage's step lines |
| `evidence.html` | With `--evidence-page`: one dossier per family (see below) |

With `--keep-work`, `out/work/` also holds every step's files, among them
`sd_clusters.tsv` (the clusters tested for segmental duplication, and the
call) and the sufr index the run built.

### Evidence page

`--evidence-page` writes `evidence.html`, a self-contained page with one
dossier per family, drawn from the family's founding element: LTRs, assembly
gaps, TSD, primer binding site, polypurine tract and protein domains on a
genome-browser style track, the 5′ and 3′ LTR alignment, and the family's
full-length copies and the sequences they lie on.  The element is drawn in
the orientation of the family's models: a minus-strand element is shown
reverse-complemented, with its location and the orientation evidence in the
dossier.  Each family gets a
triage verdict (likely novel or well supported, known family, weak
evidence, possible artifact) with the evidence behind it.
`--evidence-reference LIB` adds the nearest family in a curated library.

Run by hand on the `work/` directory of a `--keep-work` run, `pipeline/evidence_page.py` also takes a RepeatMasker run of
the genome with the library (`--mask`), the curated annotation
(`--standard`) and another tool's mask (`--other NAME=OUT`).  It then shows
each family's masking inside and outside the curated LTR annotation and how
much of the rest the other tool finds too, and `--families top-fp:N`
selects the families with the most masking outside it:

```bash
pipeline/evidence_page.py out/work genome.fa -o fp.html --families top-fp:30 \
    --reference curated.fa --mask genome.fa.out --standard curated.fa.out \
    --other HiRepEat=hirepeat.fa.out
```

### Accuracy

Each library was scored by masking the genome with RepeatMasker and
comparing with a curated annotation of the same genome: sensitivity and
precision over the bases the curated LTR annotation covers.  Recovered
families are the curated LTR and internal models that one library record
covers over at least 80% of its length.  Every row was scored with the same
recipe (`eval/benchmarks/recipe/compare_libs.sh`) against the same
standards.  The first six genomes were used for tuning; the others were not.

| Genome | Library | Sensitivity | Precision | Recovered LTR / INT | Time | Peak memory |
|--------|---------|-------------|-----------|---------------------|------|-------------|
| Arabidopsis (TAIR10.1) | ltr_library.py | 91.05% | 95.46% | 121/150 / 123/150 | 28 min | 1.3 GB |
| | LTRPipeline | 72.85% | 92.08% | 62/150 / 50/150 | 10 min | 13.9 GB |
| | HiRepEat | 64.24% | 99.49% | 16/150 / 14/150 | 3 min | 0.5 GB |
| Rice (MSU7) | ltr_library.py | 93.56% | 98.42% | 309/488 / 306/409 | 76 min | 12.0 GB |
| | LTRPipeline | 92.07% | 87.58% | 204/488 / 206/409 | 134 min | 13.7 GB |
| | HiRepEat | 83.10% | 99.61% | 87/488 / 75/409 | 19 min | 1.0 GB |
| Zebra finch (bTaeGut7) | ltr_library.py | 76.08% | 86.28% | 130/258 / 32/34 | 50 min | 2.0 GB |
| | LTRPipeline | 77.56% | 68.73% | 98/258 / 27/34 | 298 min | 14.8 GB |
| | HiRepEat | 70.58% | 90.86% | 92/258 / 29/34 | 20 min | 1.4 GB |
| Lotus (ChineseLotus 1.1) | ltr_library.py | 92.63% | 94.20% | 288/433 / 296/367 | 276 min | 14.3 GB |
| | LTRPipeline | 87.96% | 87.32% | 202/433 / 173/367 | 68 min | 17.8 GB |
| | HiRepEat | 72.55% | 99.34% | 63/433 / 89/367 | 6 min | 0.8 GB |
| Mouse (mm10) | ltr_library.py | 44.21% | 96.38% | 117/525 / 43/116 | 258 min | 38.5 GB |
| | LTRPipeline | 48.64% | 25.43% | 72/525 / 35/116 | 552 min | 21.7 GB |
| | HiRepEat | 24.04% | 65.46% | 61/525 / 30/116 | 31 min | 3.1 GB |
| Gorilla (mGorGor1 v2.1) | ltr_library.py | 40.25% | 97.59% | 104/485 / 55/135 | 269 min | 15.8 GB |
| | LTRPipeline | did not finish | | | | |
| | HiRepEat | 17.68% | 95.97% | 48/485 / 23/135 | 72 min | 5.0 GB |
| *C. elegans* (held out) | ltr_library.py | 31.14% | 66.17% | 11/23 / 9/18 | 5 min | 0.9 GB |
| | LTRPipeline | 55.24% | 53.28% | 13/23 / 10/18 | 7 min | 7.4 GB |
| | HiRepEat | no library (0 records) | | | 2 min | |
| Fruit fly (held out) | ltr_library.py | 92.67% | 91.05% | 95/142 / 86/130 | 34 min | 1.5 GB |
| | LTRPipeline | 64.21% | 80.09% | 28/142 / 41/130 | 40 min | 16.4 GB |
| | HiRepEat | 82.47% | 98.60% | 62/142 / 64/130 | 3 min | 0.3 GB |
| Soybean (held out) | ltr_library.py | 98.73% | 79.72% | 235/251 / 241/250 | 520 min | 33.3 GB |
| | LTRPipeline | 92.42% | 81.74% | 56/251 / 149/250 | 677 min | 23.8 GB |
| | HiRepEat | 97.75% | 84.28% | 119/251 / 125/250 | 57 min | 2.7 GB |
| Horse (held out) | ltr_library.py | 23.86% | 92.75% | 31/170 / 14/39 | 66 min | 24.1 GB |
| | LTRPipeline | 41.22% | 14.67% | 24/170 / 5/39 | 81 min | 19.9 GB |
| | HiRepEat | 13.79% | 96.38% | 13/170 / 3/39 | 10 min | 2.6 GB |
| Opossum (held out) | ltr_library.py | 79.13% | 89.30% | 53/116 / 47/60 | 451 min | 10.7 GB |
| | LTRPipeline | 76.59% | 22.27% | 35/116 / 26/60 | 765 min | 28.6 GB |
| | HiRepEat | 48.63% | 70.73% | 19/116 / 16/60 | 32 min | 4.2 GB |

How the rows were made:

- **ltr_library.py**: this release's code and defaults (frozen 2026-11-19,
  commit d80fbdb), 8 threads, run 2026-11-20/21.  The changes made after
  that commit (output layout, logging, models turned to their coding
  strand) are not expected to change what a library masks; they were not
  rescored.  The time does not count
  the sufr index where it was built beforehand (10 to 28 min on the large
  genomes); horse's includes it.
- **LTRPipeline**: RepeatModeler 2.0.9's, `-threads 8`, from the Dfam
  TETools image with digest
  `af5f0716b49fa6421548f389f2d7930a76de413ba85e5fe3d4ab5972aa294603`
  (GenomeTools 1.6.4, LTR_retriever 2.9.0, RepeatMasker 4.2.6).  It runs
  LTRharvest with its defaults, LTR_retriever with `-noanno`, then MAFFT,
  NINJA and Refiner; its time includes building the index.  Suffixerator and
  LTRharvest run on one thread whatever `-threads` says.  On gorilla it did
  not finish: the first attempt stopped after 168 hours in LTRharvest, and a
  second was stopped after 54 hours, still in LTRharvest.  Given the gorilla genome with ltr-forage's tandem arrays
  (713.7 Mb) replaced by N, it finished in 82 min (23.3 GB, 128 models).
  Scored on the original genome, that library gives 36.32 / 17.44 (40/485
  LTR and 17/135 INT families) and masks 15.9% of the genome; four records
  that are Alu or L1 sequence make 71% of its false-positive bases.  Its
  input differed from every other row's, so it is not in the table.
- **HiRepEat** (Gene Myers et al.): 16 threads, its published setting.  It
  produced an empty library on *C. elegans*.
- **Times and memory show scale, not a benchmark.**  The runs used two
  machines (64 and 56 cores), each shared with other jobs for part of most
  runs, and HiRepEat ran on twice the threads.  Fair timings are for a later
  release.

What the numbers do not show:

- **Base-level sensitivity rewards masking broadly.**  On mouse and horse
  LTRPipeline has the higher sensitivity, and it masks 23.6% and 15.8% of
  the genome against 5.7% and 1.5% for ltr_library.py.  Most of its lead
  comes from records that are mostly not LTR sequence (mouse: one record
  masks 157 Mb of L1) or from single-copy "internal" records that carry a
  whole LTR of another family and so mask that family's solo LTRs genome
  wide.  Restricted to its records whose own masking is at least half LTR,
  LTRPipeline scores 21.92 / 93.18 on horse (ltr_library.py: 23.86 / 92.75)
  and 35.88 / 97.16 on mouse (44.21 / 96.38).  Curated-family recovery,
  which a broad record does not raise, favours ltr_library.py on every
  mammal.  Neither tool models the old mammalian MaLR families, which
  survive mostly as solo LTRs.
- **Scope.**  The release was validated on plants, birds, insects,
  primates and rodents.  Other SINE- and LINE-rich genomes are untested:
  on opossum one record (`fam0251_v1_LTR`, an L1 fragment between two
  direct copies read as LTRs) makes 37% of the false-positive bases.  Read
  the large-record line in `run.log` (`--mask-warn`) before using a library.
- **Fly precision.**  The fly's precision is 5.3 points below the earlier
  code's (96.39%).  The change of profile file to `data/TEdomains.hmm`
  accounts for most of it (with `data/TEfam.hmm`: 95.03%); see
  [Held-out genomes](#held-out-genomes).

Earlier comparisons, from the development of the driver, are in
[docs/accuracy-history.md](docs/accuracy-history.md).

#### Rule ablation

Several of the driver's rules were added while tuning on one genome.  To
check whether each one generalizes, every rule was switched off on its own
on all six genomes, with the rest of the run unchanged (2026-10-29/31, with
the code of that date; see [docs/accuracy-history.md](docs/accuracy-history.md)).  The table gives the change in sensitivity /
precision, in points, with the rule off:

Mouse values in this section predate the depth-jump check.

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
  mouse reaches 42.65% / 79.40% (default at the time 42.93% / 76.20%).  It stays off
  by default; it is worth trying on large genomes, where the time matters
  most.

Greedy step 12 was forced on (`--greedy-above 0`, 2026-11-01; Arabidopsis
with `--greedy-batch 300`, as it has too few elements for 2,000) and each
genome rerun from step 12 to compare it with the all-vs-all search:

| Sensitivity / precision (families) | Arabidopsis | Zebra finch | Rice | Gorilla | Lotus | Mouse |
|---|---|---|---|---|---|---|
| Greedy step 12 | 91.89 / 95.61 (220) | 75.22 / 86.71 (112) | 93.42 / 98.43 (440) | 36.53 / 97.92 (124) | 92.39 / 93.60 (641) | 42.94 / 92.90 (165) |
| All-vs-all | 91.89 / 95.61 (220) | 75.21 / 86.71 (112) | 93.44 / 98.43 (442) | 36.78 / 97.96 (127) | 92.43 / 93.78 (645) | 42.92 / 92.90 (167) |

In stand-alone tests of step 12, greedy clusters are nearly pure (about 99.9%
or more of the element pairs it puts together are together in the all-vs-all
clusters), and it errs toward splitting; lotus keeps 76% of the all-vs-all pairs with 100 representatives
per cluster against 49% with 20, for 3 more minutes.  The later steps absorb
most of the splitting.  The search itself took 0.5, 4.6, 15.5, 8.5 and 70
minutes with 8 threads, against 0.5, 3, 20–27, 5.5 and 87–105 with 16.  It
helps most where hits per element are high: 27 minutes on mouse (8 threads)
against 60–110 (16 threads).  On mouse the depth-jump check drops the same two
L1 models as after all-vs-all clustering.

Greedy holds one bucket of a round's hits at a time (`--bucket-hits`, 2 M by
default, split by query element), and its singleton scores are dropped as
elements join clusters.  Lotus peaks at 2.0 GB, and its clusters are identical
with the buckets on or off.

#### Held-out genomes

The driver was run with its defaults, unchanged, on genomes that were
never used for tuning.  The table below is the code of 2026-11-05/06,
before the release; zebrafish and maize have not been rerun with the
release.  The release's own held-out runs are in the
[Accuracy](#accuracy) table and under Release validation below.  The
standards are RepeatMasker runs:

- fruit fly, *C. elegans* and zebrafish: the curated Dfam 3.9 families of
  the species and its ancestors;
- soybean and maize: RepBase 31.09's entries for the species;
- horse and opossum (release only): RepBase 31.10's entries for the species.

All runs used 8 threads.  Fly, *C. elegans*, zebrafish and soybean were rerun
(2026-11-05/06) with the code maize used, greedy step 12 and the depth-jump
check included; their times are from a busier machine than the first runs
(86, 7, 565 and 881 min).

| Genome | LTR share | Sensitivity | Precision | Recovered LTR / INT | Time |
|--------|-----------|-------------|-----------|---------------------|------|
| *D. melanogaster* (Release 6) | 11.8% | 92.66% | 96.39% | 94/142 / 85/130 | 113 min |
| *C. elegans* (WBcel235) | 0.6% | 23.74% | 53.85% | 9/23 / 8/18 | 11 min |
| Zebrafish (GRCz13ab) | 5.6% | 46.53% | 95.00% | 266/390 / 267/366 | 823 min |
| Zebrafish, without DIRS and Ngaro | 3.5% | 75.37% | 94.97% | | |
| Soybean (Glycine_max_v4.0) | 28.2% | 99.20% | 79.49% | 235/251 / 241/250 | 1,094 min |
| Maize (Zm-B73-REFERENCE-NAM-5.0) | 75.0% | 98.07% | 98.75% | 521/788 / 537/667 | 3,012 min |

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
  miss copies more than the library invents them.  The depth-jump check
  dropped one LTR model here, and it is a real LTR: RepeatMasker calls its
  origin GmGYPSY10_LTR (9% divergence).  Its median depth rose 45-fold
  after the cuts and its element has no TSD or TG..CA ends.  Left as is:
  held-out results are not used for tuning.
- **Maize.** ltr-forage kept 77,485 elements (155,000 LTR records).  Step
  12's all-vs-all search was first estimated at about 10 days for the LTRs
  alone, so this run (2026-11-02/04) used the later greedy step 12, which
  the probe chose by itself (2,365 M projected hits), and the depth-jump
  check.  The library has 839
  families.  The genome is 75% LTR by the standard, so precision alone says
  little here: specificity is 96.26%, with 20.4 Mb of false-positive
  masking.  Most of the time went to depth cleanup (1,360 min), the flank
  filter (690) and greedy clustering (477, with its probe).  The flank
  search grows much faster than its input: 151 MB of flanks took 690 min
  against 33 min for soybean's 22 MB.  The check on internal copies dropped
  four models here, and all four are LTR sequence, not another repeat: they
  match curated Copia and Huck internals at 92-95% identity, and their
  copies probably belong to close relatives with many more elements.  The
  rule assumes otherwise.  It now applies only with `--early-tier`, the
  one setting where it changed a tuning genome, so a default run keeps
  these four models; the maize row predates that.

**Release validation (2026-11-19/20).**  The frozen release (commit
d80fbdb, rmblastn 3.1.0, nail 0.7.1), run with all defaults and 8 threads
on a 56-core machine of our own, scored against the same standards.  For
part of the time two RepeatMasker jobs (the horse and opossum standards)
shared that machine:

| Genome | Sensitivity | Precision | Earlier row | Time |
|--------|-------------|-----------|-------------|------|
| *D. melanogaster* | 92.67% | 91.05% | 92.66% / 96.39% | 34 min |
| *C. elegans* | 31.14% | 66.17% | 23.74% / 53.85% | 5 min |
| Soybean | 98.73% | 79.72% | 99.20% / 79.49% | 520 min |

- **Fly lost 5.3 points of precision, and it is a known limitation of this
  release.**  The release set itself a limit of 3 points on held-out
  precision, and fly fails it.  Sensitivity did not move, but
  false-positive masking went from 0.59 to 1.55 Mb and kept elements from
  947 to 1,111.  The cause is the profile set, `data/TEdomains.hmm`, which
  replaced LTR_retriever's `TEfam.hmm` on 2026-11-17.  Rerunning the release
  on fly with one change reverted at a time gives:

  | Fly, frozen release with | Sensitivity | Precision | Kept elements |
  |--------------------------|-------------|-----------|---------------|
  | nothing changed | 92.67% | 91.05% | 1,111 |
  | `--domain-hmm data/TEfam.hmm` | 93.08% | 95.03% | 968 |
  | `--domain-z 0` (E-values from the run's own count) | 92.63% | 88.92% | 1,113 |
  | nail `--prog-f 0.01` (nail's default) | 92.87% | 90.99% | 1,110 |
  | all three reverted | 93.22% | 94.09% | 996 |

  The new profiles' Pfam domains are not in `TEfam.hmm`: the two reverse
  transcriptase RNase H domains (PF17917, PF17919), the Gypsy protein
  (PF07253) and the integrase zinc-binding domain (PF17921).  Those
  profiles now give strong domain calls (most under E = 1e-10) to
  candidates whose two "LTRs" are young copies of a LINE (R1, CR1, Jockey)
  or a DNA transposon (ProtoP).  The copies sit inside or beside a real
  Gypsy element, so the candidate's "internal" region carries that
  element's coding sequence.  The domain calls are right; the elements are
  nesting artifacts.  Fly is full of young, active, densely nested
  elements, so it shows this most.  The tuning genomes chose these profiles,
  and held-out results are not used for tuning, so the profile set stays.
  The remaining 1.4 points against the earlier row come from the other
  changes since 2026-11-05 and were not traced.  The fix that would address
  this, a test of whether the coding sequence belongs to another element
  nested between the two copies, is future work.

### Known limitations and future work

- **Nested elements.**  A candidate whose two "LTRs" are copies of another
  repeat beside or inside a real element carries that element's coding
  sequence and is kept on its domains (fly precision, above).  A test of
  whether the coding sequence belongs to an element nested between the two
  copies is the fix.
- **SINE- and LINE-rich genomes outside primates and rodents** are untested;
  on opossum one L1 record makes 37% of the false-positive bases.  Its
  "LTRs" grew from 218 bp to 2.25 kb in re-extension over a direct
  duplication containing L1.
- **Benchmarks.**  Base-level sensitivity rewards libraries that mask
  broadly (see [Accuracy](#accuracy)); measures that resist this, such as
  each record's own precision or curated-family recovery, are to come.
- **Timings** in this README are from shared machines; fair timings are for
  a later release.
- **Tandem mask left cut.**  The mask misses 0.3% to 0.9% of kept elements
  at an array's left edge; a patch exists.
- **Overlap choice.**  A PBS call outranks LTR identity outright, so a
  one-mismatch PBS can beat a candidate several points more identical.
- **The poly-A rule vetoes some real elements** (soybean, maize), and the
  poly-A family veto's pool guard is untested between 80 and 17,000.
- **Gorilla `fam0118`**, a SATR1 satellite pair, masks 0.52 Mb.
- **The depth-jump check drops a real soybean LTR** (held out, so left).
- **nail and MMseqs2 miss some strong hits** on large target sets.
- **A general "is this model an LTR element" test** to replace the stack of
  rules.
- **The merged RepeatModeler library** (this pipeline in place of
  LTRPipeline) is to be scored.
- **Run time and memory:** step 14 is the longest step on large genomes, and
  mouse peaks at 38 GB.
- **sufr 0.8.0** changes the index format; supporting it needs ltr-forage's
  index reading ported.
- **ltr-forage's PPT columns** are computed on the plus strand only, so for
  minus-strand elements they describe the wrong side (the evidence page
  recomputes the PPT on the coding strand).

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
    --overlaps best -t 32 > hg38_ltr-forage.out
```

These are LTRharvest's default parameters written out in full.  Wrap each
command in `/usr/bin/time -v` for wall time and peak memory; `ltr-forage`
also prints the elapsed time at the end of each phase with `-v`.  `gt`
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
used here has 7,002 such elements (internal regions are the `-int` annotations; an
earlier version of the script also took any name ending in `I`, such as the LTR
MLT1I, for one, and counted 7,170).  A prediction hits one when each of its LTRs overlaps the
corresponding annotated LTR by at least half the length of the shorter.

```bash
eval/rm_truth_recall.py hg38.fa.out hg38.fa.fai hg38_gt.out hg38_ltr-forage.out \
    hg38cmp.a_only.tsv hg38cmp.b_only.tsv
```

Recall by the larger of the two LTRs' divergence from consensus:

| Max LTR divergence | Elements | LTRharvest | ltr-forage |
|---|---|---|---|
| 0-5% | 585 | 90.9% | 90.9% |
| 5-10% | 1,321 | 63.7% | 63.7% |
| 10-15% | 1,754 | 21.6% | 21.6% |
| 15-20% | 1,264 | 7.8% | 7.8% |
| 20% and up | 2,078 | 2.5% | 2.5% |

Both tools hit the same 1,903 elements, and none of the 28 predictions
unique to one tool hits any of them.  The differences between the two
outputs fall entirely on predictions that do not correspond to intact
annotated elements.

```bash
eval/rm_classify.py hg38.fa.out hg38.fa.fai hg38_gt.out
```

This labels each prediction by the repeat class under its two LTRs.  For
either tool, 62% of hg38 predictions are pairs of SINE elements, 14% pairs of
LINE fragments, 8.6% LTR class (3,229 of LTRharvest's 3,595 with the same family
covering most of both sides), and 11% mostly unannotated.  That is a property of the LTRharvest search at these
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
2018), distributed under GPL v3.  `data/TEdomains.hmm` is built from Pfam
39.0 profiles (Pfam is released under CC0) and six profiles copied from
`TEfam.hmm`.

The LTRharvest algorithm is implemented based on the published description and
the open-source GenomeTools reference implementation (ISC license).  The `sufr`
library is used under its BSD-3-Clause license.
