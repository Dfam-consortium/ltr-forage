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
HMMER3 profile file with nail, and reports whether the hits
belong to LTR retrotransposon proteins or to other transposon and host
proteins. On hg38 with LTR_retriever's TEfam profiles, a hit of LTR order
marks 68% of the RepeatMasker-annotated elements and 0.5% of the other
candidates, which is the evidence LTR_retriever gets from HMMER, found here
on more elements.

By default, `ltr-forage` produces output that is format-compatible with
LTRharvest, making it a transparent replacement in existing pipelines.

`pipeline/ltr_library.py` builds on the candidates to produce a curated LTR
family library from a genome assembly: one LTR and one internal model per
family, cleaned of nested insertions and segmental duplications, in place of
LTRharvest followed by LTR_retriever (see
[Building a family library](#building-a-family-library)).

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
translation against a HMMER3 file of profile HMMs with
[nail](https://github.com/TravisWheelerLab/nail). LTR_retriever's
`TEfam.hmm` (102 Pfam profiles) is the intended input; any HMMER3 amino-acid
profile file works.

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
and nail's databases go to a subdirectory of `--domain-tmp` (default: the
system temp directory), removed afterwards.

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

### Phase 6 — Decision and library (`--keep`, `--library`)

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

`--domain-only-min-identity` (default 70) rejects a candidate whose only
evidence is a domain hit (no TSD, no TG..CA and no PBS) when its final LTR
identity, measured after re-extension, is below that percentage. `--similar`
is checked on the X-drop seed alignment, before re-extension grows the LTRs
to their boundaries, so a pair of duplicated segments with a stray integrase
hit can pass it at 89% and end at 70%. Such pairs gave the gorilla library
1.2 kb "LTRs" in a 440 bp family. Among kept elements below 70% that rest on a
domain hit alone, 41% of gorilla's lie in curated LTR annotation, against
83-84% for rice and lotus and 98-100% for those with a PBS; the rule removes
17 gorilla, 39 rice, 115 lotus, 2 zebra finch and 1 Arabidopsis element.

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
2). LTR_retriever ranks a 5 bp TSD first. There is no separate
tandem-repeat finder: of the 1,797 candidates the earlier 5 bp rule kept
on hg38 before the self-alignment check, 11 have an LTR that is mostly satellite,
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

### Phase 7 — Primer binding site (`--format tsv`, `--trna-out`)

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
- For the domain search (`--domain-hmm`): `nail` ≥ 0.7.1
  ([TravisWheelerLab/nail](https://github.com/TravisWheelerLab/nail)) and
  MMseqs2 (`mmseqs`), which nail seeds with

For the family-library driver (`pipeline/ltr_library.py`) also:

- Python ≥ 3.8 with numpy
- The Rust port of RMBlast, `rmblastn` ≥ 3.0.9
  ([Dfam-consortium/rmblast](https://github.com/Dfam-consortium/rmblast));
  NCBI's C++ `rmblastn` does not take the options the driver uses; 3.0.11 or later
  is recommended (3.0.10's streaming `--mt-mode 1` can deadlock, and the driver
  falls back to `--mt-mode 0` with it)
- `te-composer` from dfam-curator ≥ 0.6.6 (consensus models)
- An HMMER3 profile file of TE proteins for the domain search, such as
  LTR_retriever's `TEfam.hmm`

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
ltr-forage --trna-scan <FASTA> [-j N] [--trna-intron-max BP]

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
      --nail <PATH>         nail binary (a path, or a name on PATH) [default: nail]
      --mmseqs <PATH>       MMseqs2 binary nail seeds with [default: mmseqs]
      --domain-tmp <DIR>    Directory for the domain search's scratch files [default: system temp dir]
      --keep-coding-ltrs    Keep candidates with a polyprotein domain hit inside an LTR (by default
                            they are vetoed: an LTR does not encode the polyprotein, so the hit
                            means the "LTR" is internal sequence)

Decision and library options:
      --keep <LEVEL>        Evidence a candidate needs: domain | standard | loose | none [default: standard]
      --internal-gap-max <FRAC>  Also reject candidates whose internal region has a single-copy
                            stretch longer than this fraction of it (needs --internal-stats)
      --library <FILE>      Write the kept elements as a FASTA library (left LTR, internal
                            region, right LTR per element)
      --tandem-max <FRAC>   Reject a candidate when either LTR aligns to itself off the diagonal
                            over at least this fraction of its length [default: 0.99]
      --domain-only-min-identity <PCT>  Reject a candidate kept on a domain hit alone (no TSD,
                            TG..CA or PBS) below this final LTR identity; 0 disables [default: 70]

Primer binding site options (tsv output, or for `--domain-only-min-identity`):
      --no-pbs-trna         Skip the genome tRNA scan behind the PBS columns
      --pbs-domain-strand   Search for the PBS only on the strand the protein domains give (default: both)
      --trna-out <FILE>     Write the genome's tRNA gene calls (the PBS probe source)
      --trna-intron-max <BP>  Longest tRNA intron searched for [default: 3000]
      --trna-scan <FASTA>   Only scan FASTA for tRNA genes, print the calls and exit

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
| `tsd_len`, `tsd_seq`, `tsd_shift_l`, `tsd_shift_r` | Longest exact direct repeat (4 to 20 bp) flanking the element, allowing 0 or 1 bp of play at each end; failing that, a 6 or 5 bp repeat with one mismatch directly against both LTR ends (see `tsd_mm`). `tsd_seq` is the left copy; `.` and 0 when none |
| `motif`, `motif_class` | First two bases of the left LTR and last two of the right LTR, and LTR_retriever's reading of them: `TGCA` for the canonical termini, `known` for the other pairs its rules accept (TGCT, TACA, TACT, TGGA, TGGT, TATA, TGTA, TGCC), `other` |
| `ppt_purines15`, `ppt_offset` | Purine count in the most purine-rich 15 bp window within 30 bp before the right LTR, and that window's start relative to the right LTR |
| `pbs_offset` | Start of the primer binding site relative to the first base after the 5′ LTR (−5 to +15), from the 12-nt tRNA match described under Phase 7; −1 when none or with `--no-pbs-trna` |
| `n_count` | Ns inside the element |
| `xdrop_tsd_len` | TSD length the X-drop stage's own search reported |
| `int_tiles`, `int_multi_frac`, `int_gap_frac`, `int_gap_bp` | With `--internal-stats`: number of 20-mer tiles taken every 10 bp across the internal region (tiles containing N skipped), the fraction seen at least three times in the genome on either strand, and the longest run of single-copy tiles as a fraction of the tiles and in bases; `.` otherwise |
| `dom_order`, `dom_superfamily`, `dom_strand` | With `--domain-hmm`: `LTR`, `notLTR`, `mixture`, `NA` or `none`; `Gypsy`, `Copia` or `unknown`; `+`, `-` or `?`, as described under Phase 5; `.` otherwise |
| `dom_best_e`, `dom_hits` | Best E-value, and every hit as `accession:E-value:frame` (frames `aa1`..`aa3` forward, `rev_aa1`..`rev_aa3` on the reverse complement), comma separated; `.` when none |
| `ltr_dom_hits` | With `--domain-hmm`: hits to the LTR-order (polyprotein) profiles when each LTR is searched on its own, summed over the two LTRs; 0 when `--keep-coding-ltrs` skips the search |
| `lltr_self_cov`, `rltr_self_cov` | Fraction of each LTR covered by plus-strand alignments of the LTR against itself off the main diagonal; 1.00 for a satellite monomer or an Alu dimer, 0 for unique sequence |
| `keep`, `evidence` | 1 when the candidate passes `--keep`, and what decided it: `domain`, `tsdN`, `tgca` (comma separated), `veto:notLTR`, `veto:mixture`, `veto:tandemF`, `veto:ltrdomainN`, `veto:lowidP` (`--domain-only-min-identity`), `gapF` when `--internal-gap-max` rejected it, or `none` |
| `pbs_trna`, `pbs_mm12` | The tRNA type(s) whose 3′-terminal 12 nt match the PBS, e.g. `Met(cat)` (several joined by `/` when they share the probe), and the mismatches (0 or 1); `.` when none |
| `pbs_edist18` | Best edit distance (0 to 2) of a full 18-nt tRNA probe over the same window; `.` when none |
| `pbs_strand` | `+` when the PBS was found after the left LTR, `-` when before the right LTR on the reverse strand; `.` when none |
| `tsd_mm` | 1 when the TSD is the one-mismatch fallback, else 0. It is enriched about 20-fold over shifted flanks only for TG..CA elements (lotus), so it counts toward `--keep standard`/`loose` only alongside TG..CA; the evidence tag reads `tsd5mm1` |

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

## Building a family library

`pipeline/ltr_library.py` runs ltr-forage on a genome and turns its
candidates into a family library for RepeatMasker or RepeatModeler's
LTRPipeline:

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
profiles default to LTR_retriever's `TEfam.hmm` (`--domain-hmm`). A sufr
index beside the genome (`genome.sufr` for `genome.fa`) is used when it
exists; `--sufr` names another. Each run writes the command and every
program's path and version to `out/versions.txt`, and its step times to
`out/steps.log`. The step scripts live beside the driver in `pipeline/`.

### Steps

1. **Index and depth sample.**  A sufr index, and a genome sample for copy
   depth: the whole genome up to 200 Mb, else 10% of it in 100 kb windows.
   Depth thresholds scale with genome size and the sample fraction.
2. **ltr-forage** with re-extension, tiered TSDs, protein domains and
   the decision step, giving the candidate table (`cand.tsv`) and a raw
   library of each kept element's LTR and internal region.
3. **Raw clustering and consensus representatives.**  The raw library is
   clustered all-vs-all with rmblastn.  Each cluster of at least
   `--min-cons` (3) elements is represented from here on by a te-composer
   consensus, which drops alignment columns held by fewer than half the
   median number of copies.  This keeps an insertion carried by a few copies
   out of the model.
4. **Depth cleanup.**  Three rounds of rmblastn search of the sample against
   the library.  Stretches far above their record's own depth are excised
   as nested passengers: 3x for an element, 10x for a consensus
   representative, whose family is high-copy itself.  High-copy LTRs are
   dropped, but never a consensus representative.
5. **Nested-insertion vote and element clusters.**  A second clustering
   excises LTRs nested in internals by vote.  Each big cluster's members
   rejoin the cluster of their representative.
6. **Segmental duplications.**  A cluster is dropped when most of its
   pairs, or every one of its copies, share both flanks (300 bp at 80%
   identity) with another copy: the copies are one duplicated locus, not
   independent insertions (`--no-flank` skips this).
7. **Containment, nested excision, PBS and tiering.**  Singletons are kept
   only with LTR identity ≥ 93% or a primer binding site (`--tier strict`:
   identity ≥ 90% and a PBS).
8. **Models and families.**  A consensus model for each cluster of at least
   `--min-cons` elements and the exemplar element for smaller ones, grouped
   into 80-80 families (LTR and internal alike) with Gypsy and Copia
   classes from the domain hits.

`--stage-a` skips step 3 and depth-cleans and re-clusters every element, as
the driver did before.  It is slower on large libraries and no more
accurate.  `--min-cons 0` gives exemplar models only.  Each step's command
lines and output go to `out/logs/`, and timings to `out/steps.log`.  A rerun
skips finished steps; `--from STEP` reruns from a step on.

### Output

| File | Contents |
|------|----------|
| `library.fa` | The family library: `fam<N>_v<k>_LTR` and `fam<N>_v<k>_INT` records, classified `#LTR/Gypsy`, `#LTR/Copia` or `#LTR` |
| `families.tsv` | Each model's family, origin element, kind, length, status and class |
| `cand.tsv` | ltr-forage's per-candidate table, with every value behind its decision (see [Per-candidate table](#per-candidate-table---format-tsv)) |
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

| Genome | Library | Sensitivity | Precision | Recovered LTR / INT |
|--------|---------|-------------|-----------|---------------------|
| Rice (MSU7) | ltr_library.py | 92.70% | 98.53% | 312 / 306 |
| Rice (MSU7) | HiRepEat | 83.10% | 99.61% | 80 / 76 |
| Zebra finch (bTaeGut7) | ltr_library.py | 74.76% | 83.03% | 24 / 30 |
| Zebra finch (bTaeGut7) | HiRepEat | 70.58% | 90.86% | 16 / 27 |
| Gorilla (mGorGor1 v2.1) | ltr_library.py | 34.57% | 97.65% | 2 / 8 |
| Gorilla (mGorGor1 v2.1) | HiRepEat | 17.68% | 95.97% | 1 / 3 |

Much of the zebra finch false-positive masking is in regions the curated
library leaves unannotated.  About 95% of HiRepEat's falls there too, and
several of these models are dispersed elements with retrotransposon domains,
which suggests families missing from the curated library rather than
artifacts.

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

5. Laslett D, Canback B. ARAGORN, a program to detect tRNA genes and tmRNA
   genes in nucleotide sequences. *Nucleic Acids Res.* 2004;**32**(1):11–16.
   DOI:[10.1093/nar/gkh152](https://doi.org/10.1093/nar/gkh152)

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

The LTRharvest algorithm is implemented based on the published description and
the open-source GenomeTools reference implementation (ISC license).  The `sufr`
library is used under its BSD-3-Clause license.
