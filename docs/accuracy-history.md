# Accuracy history

Dated comparisons from the development of `pipeline/ltr_library.py`, moved here from the README's Accuracy section for release 0.1.0.  They describe earlier code and defaults; the README has the release's scores.  Section links point to the README.

## Before release 0.1.0

Library masking was scored against curated libraries: sensitivity and
precision over the bases the curated LTR annotation covers.  Recovered
families are the curated LTR and internal models covered at 80% or more.
Every library in the table was scored with the same recipe
(`eval/benchmarks/recipe/compare_libs.sh`).

The three tools were run under different conditions, so the time column
shows scale, not a benchmark.  All times are wall-clock on a shared 64-core
machine with other jobs running.

- **ltr_library.py**: the driver as it was on 2026-11-13, with its defaults
  of that date (`data/TEfam.hmm`, no step 11b), 8 threads, one genome at a
  time.  The current code and defaults are in the second table below.
  Steps 1-13 are from end-to-end runs of 2026-11-11/12
  and steps 14 on from a rerun on 2026-11-13 with the remapped depth rounds;
  the time is the two parts added.  The sufr index was already built and is
  not counted (it took 10 min for zebra finch and 15 min for gorilla).
- **LTRPipeline**: RepeatModeler 2.0.9's, run with `-threads 8` from the Dfam
  TETools image with digest
  `af5f0716b49fa6421548f389f2d7930a76de413ba85e5fe3d4ab5972aa294603`
  (GenomeTools 1.6.4, LTR_retriever 2.9.0, RepeatMasker 4.2.6), starting
  2026-11-14.  It runs LTRharvest with its defaults, LTR_retriever with
  `-noanno` (which uses RepeatMasker, CD-HIT and TRF), then MAFFT, NINJA and
  Refiner.  Its time includes building the index.  In this release
  suffixerator and LTRharvest run on one thread whatever `-threads` says;
  they took 2, 6, 244 and 8 min of the Arabidopsis, rice, zebra finch and
  lotus totals.  A 16-thread RepeatModeler job shared the machine during the
  rice, zebra finch and lotus runs and an 8-thread job during Arabidopsis.
  The mouse and gorilla runs had not finished when this was written.
- **HiRepEat**: 16 threads.

| Genome | Library | Sensitivity | Precision | Recovered LTR / INT | Time |
|--------|---------|-------------|-----------|---------------------|------|
| Arabidopsis (TAIR10.1) | ltr_library.py | 91.12% | 95.54% | 120/150 / 120/150 | 22 min |
| Arabidopsis (TAIR10.1) | LTRPipeline | 72.85% | 92.08% | 62/150 / 50/150 | 10 min |
| Rice (MSU7) | ltr_library.py | 93.38% | 98.43% | 308/488 / 307/409 | 69 min |
| Rice (MSU7) | LTRPipeline | 92.07% | 87.58% | 204/488 / 206/409 | 134 min |
| Rice (MSU7) | HiRepEat | 83.10% | 99.61% | 87/488 / 75/409 | 19 min |
| Mouse (mm10) | ltr_library.py | 42.88% | 92.86% | 106/525 / 41/116 | 225 min |
| Mouse (mm10) | HiRepEat | 24.04% | 65.46% | 61/525 / 30/116 | 31 min |
| Zebra finch (bTaeGut7) | ltr_library.py | 75.18% | 86.70% | 113/258 / 30/34 | 43 min |
| Zebra finch (bTaeGut7) | LTRPipeline | 77.56% | 68.73% | 98/258 / 27/34 | 298 min |
| Zebra finch (bTaeGut7) | HiRepEat | 70.58% | 90.86% | 92/258 / 29/34 | 20 min |
| Lotus (ChineseLotus 1.1) | ltr_library.py | 92.46% | 94.14% | 287/433 / 298/367 | 301 min |
| Lotus (ChineseLotus 1.1) | LTRPipeline | 87.96% | 87.32% | 202/433 / 173/367 | 68 min |
| Gorilla (mGorGor1 v2.1) | ltr_library.py | 36.90% | 98.00% | 94/485 / 49/135 | 179 min |
| Gorilla (mGorGor1 v2.1) | HiRepEat | 17.68% | 95.97% | 48/485 / 23/135 | 72 min |

The recovered counts are out of each curated library's LTR and internal
records.  Zebra finch and gorilla use RepeatMasker's names, where most LTRs
carry no `LTR` suffix (THE1B, MLT1A0) and an internal may have
no suffix in the library (HERVH, annotated as HERVH-int); `family_recovery.py`
takes the part from the name where it has one and otherwise from the record's
LTR class and the standard annotation (`--standard`).  Before 2026-11-15 it
skipped every record without `LTR` or a trailing `I` in its name, and the
ltr_library.py rows of the two read 25 / 29 and 2 / 13, out of 44 / 33 and
11 / 79.  The 45 LINE,
SINE and other non-LTR records in the Arabidopsis file are left out.

**Current code and defaults.**  The `ltr_library.py` rows above are from
before the code review of 2026-11-15.  The table below is the same six
genomes with that review's fixes and the current defaults
(`data/TEdomains.hmm`, step 11b with its guards and rescue, the poly-A
repeat cut in step 22b), 8 threads, one genome at a time, scored with the
same recipe:

| Genome | Records / families | Sensitivity | Precision | Recovered LTR / INT | Time | Before the review (sensitivity / precision) |
|--------|--------------------|-------------|-----------|---------------------|------|------|
| Arabidopsis (TAIR10.1) | 412 / 205 | 91.09% | 95.50% | 121/150 / 121/150 | 33 min | 91.12% / 95.54% |
| Rice (MSU7) | 929 / 443 | 93.57% | 98.42% | 309/488 / 306/409 | 101 min | 93.38% / 98.43% |
| Zebra finch (bTaeGut7) | 216 / 105 | 76.13% | 86.28% | 130/258 / 32/34 | 65 min | 75.18% / 86.70% |
| Lotus (ChineseLotus 1.1) | 3,405 / 799 | 92.64% | 94.32% | 283/433 / 293/367 | 512 min | 92.46% / 94.14% |
| Mouse (mm10) | 458 / 202 | 44.05% | 96.47% | 116/525 / 44/116 | 300 min | 42.88% / 92.86% |
| Gorilla (mGorGor1 v2.1) | 347 / 165 | 40.36% | 97.65% | 102/485 / 54/135 | 261 min | 36.90% / 98.00% |

- Mouse and gorilla were run end to end from the repository code on
  2026-11-17/18.  Step 11b vetoed 873 of 9,444 kept elements on mouse and
  274 of 2,808 on gorilla and kept 55 and 18 again; the cut removed 7,210
  bases of 7 mouse models and nothing on gorilla.  Peak memory was 38.0 GB
  (mouse) and 15.7 GB (gorilla).
- Gorilla's library holds one known false family, `fam0118_v1_LTR`: two
  141-bp repeats 1.7 kb apart inside a SATR1 satellite array
  (`CM055456.2:70120079..70122095`), kept on a weak RNase H hit
  (PF17919, E = 1.4e-3) to sequence the satellite's repeat unit does
  carry.  It masks 0.52 Mb of satellite; without it gorilla's precision
  would be 98.10%.  Two tandem arrays cover 54% of its internal region but
  leave its LTRs free.  Widening every masked array by one unit would hide
  the pair and 77 other kept gorilla elements, 75 of them with curated LTR
  annotation on both LTRs, so the family stays.
- Arabidopsis, zebra finch, rice and lotus were run from the copy of the
  pipeline in which the two steps were developed.  Step
  11b does not run on these genomes (the pool guard), and the first
  version of it, which did run on rice and lotus, vetoed nothing there.
  They have not been rerun from the repository code.
- The times are not comparable with the first table's: the machine was
  busier (its load average was 35 to 105 during the lotus run).
- The zebra finch recovered counts rise with the profile file: the fixed
  code with `data/TEfam.hmm` recovers 114 / 28.

Step 14's remapped rounds (the default since 2026-11-13) against the same
runs with rounds 2 and 3 searched (`--no-depth-remap`), which is the
reference: the remap only saves time.

| Genome | Rounds searched | Rounds remapped | TP | FP | Step 14, min | Whole run, min |
|---|---|---|---|---|---|---|
| Arabidopsis | 91.12 / 95.52 | 91.12 / 95.54 | 0 | -1 kb | 18 -> 13 | 27 -> 22 |
| Rice | 93.38 / 98.43 | 93.38 / 98.43 | +1 kb | -1 kb | 26 -> 20 | 75 -> 69 |
| Mouse | 42.92 / 92.46 | 42.88 / 92.86 | -154 kb | -683 kb | 146 -> 100 | 271 -> 225 |
| Zebra finch | 75.18 / 86.71 | 75.18 / 86.70 | +1 kb | +1 kb | 24 -> 14 | 55 -> 43 |
| Lotus | 92.46 / 93.83 | 92.46 / 94.14 | +3 kb | -695 kb | 100 -> 86 | 315 -> 301 |
| Gorilla | 36.87 / 97.96 | 36.90 / 98.00 | +101 kb | -35 kb | 110 -> 106 | 183 -> 179 |

Rounds 2 and 3 take 0.2 to 11 minutes remapped against 5 to 52 searched;
what is left of step 14 is the round 1 search (gorilla's ran on a busier
machine than its reference).  The precision gained on lotus and mouse is
not a property of the remap.  Depth over a small record is thin, so one
hit more or less moves a segment across the 3x rule, and the cut then
moves that segment's copies to other records:

- **Lotus.**  One internal is a hit deeper in round 2 (51 against 50, at a
  threshold near 60) and loses 264 bp instead of 50.  About 80 sample loci
  then go to another internal, whose depth doubles; it is dropped where
  the searched rounds trim it to a 737-bp model that masks 0.58 Mb outside
  the curated annotation and 4 kb inside.
- **Mouse.**  One internal (chrX at 86.61 Mb) comes out 43 bp shorter and
  masks 0.65 Mb less outside the curated LTRs: the 80-bp segment noted
  below.  The 154 kb lost inside them is spread over other models (7 of
  354 models differ in origin, 35 in length).

With the searched rounds the mouse row is the end-to-end run repeated from
step 22b with the depth pass after the depth-jump drop, which changes no
other genome (without it: 89.51% precision).

The rows before the runtime changes of 2026-11-10 (incremental depth
rounds, word 16, two-pass re-extension) were Arabidopsis 91.89 / 95.61,
rice 93.44 / 98.43, mouse 42.92 / 92.90, zebra finch 75.21 / 86.71, lotus
92.43 / 93.78 and gorilla 36.78 / 97.96, in 50, 85, 124, 59, 303 and 107
minutes with 16 threads (24 on mouse and gorilla).  Incremental depth
rounds give byte-identical libraries.  Word 16 and two-pass re-extension
are not identical, and were first scored on their own, each run alongside
the previous code (8 threads each):

| Change | Rice | Lotus | Mouse | Gorilla |
|---|---|---|---|---|
| `--word 16` (step 12 on) | 93.40 / 98.43; 72 -> 60 min | 92.46 / 93.83; 539 -> 346 min | 42.93 / 92.90; 248 -> 225 min | |
| `--reextend-first 2500` | 93.43 / 98.43; detection 9.8 -> 6.9 min | | | 36.78 / 97.96; detection 103 -> 45 min |

Two genomes moved by more than 0.1 points when all three ran together
(searched rounds):

- **Arabidopsis, sensitivity 91.89 -> 91.12.**  Two-pass re-extension gives
  one old Gypsy element without a TSD (NC_003070.9 at 14.24 Mb) LTRs of
  2.9 kb at 79% identity where the full window gave 4.0 and 4.3 kb at 70%.
  The shorter LTR aligns over its whole length to two other families'
  internals, so the nested-LTR vote (step 15) cuts 1,450 bp out of each;
  the longer one aligned over 73% of its length, under the vote's 80%.
  Those two models lose 210 kb of masking inside the curated annotation,
  60 kb net of what other models take up.  A third internal keeps 950 bp
  more; the other 402 records are identical.
- **Mouse, precision 92.90 -> 92.46.**  Detection reports 113 more
  candidates, and in round 2 of step 14 one internal's median depth was 4
  instead of 3, which kept a segment under the 3x rule.  The L1 piece left
  in it is the one the pass after the depth-jump drop now removes; an
  80-bp segment of another internal (0.65 Mb outside the curated LTRs)
  remains with the searched rounds.

On mouse, depth cleanup is the longest step (100 of the 225 minutes, 94 of
them the round 1 search).

On mouse the library was at 44.29% and 42.13% before the poly-A veto and
the model depth check: 141 Mb of its false-positive masking was L1 and 46 Mb
B1/B2 SINEs, 35 Mb and 0.4 Mb after those two rules.  Most of that L1 came
from two LTR models made of L1 sequence that carry no poly-A tail; the
depth-jump check in step 22b drops them, which takes all false-positive
masking from 45 Mb to 11 Mb (76.20% to 92.90% precision).  The copy cap
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

