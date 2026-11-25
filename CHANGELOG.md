# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## Unreleased
### Changed
- sufr 0.8.0: ltr-forage reads its index format (version 7) through libsufr
  0.8.0's API, and no longer reads indexes made by sufr 0.7 or older; rebuild
  them with `sufr create --dna --write-lcp` (ltr-forage reads the LCP array,
  which sufr 0.8.0 writes only when asked).  ltr-forage and `ltr_library.py`
  say so when given an old index or one without the LCP array, and the driver
  builds its index that way, with sufr 0.8.0 or newer.  The results are
  unchanged: on the Arabidopsis example every output is byte-identical, and on
  rice and hg38 LTRharvest mode's output too.

### Fixed
- The evidence page named the run after the driver's work directory ("work")
  instead of the output directory.
- No step script falls back to a built-in program or data path any more
  (`rmb_ava.sh`, `clean_library_rmb.sh`, `evidence_page.py`).

### Documentation
- README reorganised and shortened: option tables in place of the full help,
  one output section, an example evidence-page dossier
  (`docs/evidence_example.png`, made by `docs/evidence_figure.py`); the
  algorithm notes left it.  The hg38 comparison with LTRharvest was rerun
  with this code: the same counts.

## 0.1.0 - 11/23/26
The first release: `ltr-forage`, which finds LTR retrotransposon copies in a
genome from a sufr suffix-array index, and `pipeline/ltr_library.py`, which
turns them into a family library for RepeatMasker or RepeatModeler's
LTRPipeline.

### ltr-forage
- LTRharvest's seed search (maximal repeats from the suffix array, X-drop
  extension), streaming the index from disk and parallel over threads, with
  output that does not depend on the thread count.  `--ltrharvest` runs
  LTRharvest's search alone and writes its table.
- Tandem arrays found by a port of Gene Myers' FasTAN are masked from the
  seed search, so satellite pairs no longer dominate the run.
- Each candidate's LTR pair is re-aligned in a window of up to 7 kb, and
  calls whose final ends agree within 20 bp are merged.
- Evidence for every candidate: LTR self-repeat, protein domains (nail against
  `data/TEdomains.hmm`, 122 Pfam profiles), the primer binding site from the
  genome's own tRNA genes (an ARAGORN port), TSD and TG..CA ends.  A candidate
  is kept on that evidence and vetoed on non-LTR protein hits, LTRs that
  repeat themselves, polyprotein domains inside an LTR, or poly-A-tailed LTRs.
  Overlapping calls are resolved by the evidence.
- A per-candidate table with every value behind each decision, the element's
  coding strand (`strand`: protein domains, else the PBS), and a FASTA of the
  kept elements with the strand in each header.
- The screen shows each stage and a summary; `-v` adds every step's numbers,
  `--log FILE` keeps all of it.  `-h` lists the options most runs need,
  `--help-all` the development and tuning ones too.  `-t` sets the threads.
  nail and MMseqs2 come from `--nail`/`--mmseqs`, `NAIL`/`MMSEQS`, or `PATH`.
- `--version` names the commit the binary was built from, with `-dirty` for a
  build from uncommitted changes.

### ltr_library.py
- Families from ltr-forage's kept elements: clustering, consensus models for
  families of three or more copies, depth cleanup of nested insertions and
  high-copy passengers, segmental-duplication and containment filters, 80-80
  families with Gypsy/Copia classes, and a seed alignment per model
  (`library.stk`).
- Each model is turned to its element's coding strand (`families.tsv`
  records what set it, and whether it was flipped).
- `--evidence-page` writes `evidence.html`, one dossier per family with a
  triage verdict.
- The output directory holds the library, its seed alignments, `families.tsv`,
  ltr-forage's table, the evidence page and `run.log`; the intermediate files
  are removed after a successful run unless `--keep-work`.  A stopped run
  resumes; `--from STEP` reruns from a step.
- Programs are found from their option, environment variable, or `PATH`, and
  their paths and versions are written to `run.log`.  Every step's temporary
  files go under the output directory (`--tmpdir` to move them).

### Validation
Scored against curated annotations of six tuning genomes (Arabidopsis, rice,
zebra finch, lotus, gorilla, mouse) and five held out (C. elegans, fly,
soybean, horse, opossum); see the README's Accuracy section.  Validated on
plants, birds, insects, primates and rodents; other SINE- and LINE-rich
genomes are untested, and the large-record check in the run log should be
read before a library is used.
