# Protein-domain profiles for ltr-forage

`../TEdomains.hmm` is the profile file for ltr-forage's domain search
(`--domain-hmm`).  It is built from a Pfam release by `build_domains.py`
using the profile list in `profiles.tsv`:

    python3 build_domains.py /path/to/Pfam-A.hmm.gz        # writes ../TEdomains.hmm

The file holds 122 profiles: 67 LTR (32 of them core) and 55 not LTR.
116 are the Pfam 39.0 builds.  Six are the 2015 builds from
`../TEfam.hmm`, because the Pfam 39.0 rebuilds of those six find far fewer
elements (see "Kept as the 2015 builds"); the build reads them from
`../TEfam.hmm`, which has to be present.  One of LTR_retriever's profiles,
Atrophin-1, has been removed (see "Removed").

## How a profile carries its class and role

Each profile names its RepeatMasker class after a `#`, and an LTR
polyprotein domain ends its DESC line with `[role=core]`:

    NAME  rve#LTR
    ACC   PF00665.33
    DESC  Integrase core domain [role=core]

ltr-forage (`src/domain.rs`, `ProfileTable`) reads them as follows:

- **Type `LTR`** is evidence for an LTR candidate.  Any other type (`DNA`,
  `LINE`, `RC`, `Other` ...) is evidence against it.  A candidate whose hits
  disagree is a "mixture" and is vetoed, as in LTR_retriever.  A profile
  with no class in an annotated file also counts as disagreement.
- **Subtype `Gypsy` or `Copia`** counts toward the superfamily call.  Other
  LTR subtypes (`ERV1`, `ERVK`, `Pao`, `Caulimovirus` ...) record what the
  domain belongs to, but the call does not name them yet.
- **`[role=core]`** marks polyprotein domains: gag, protease, RT, RNase H,
  integrase and chromodomain.  A core hit inside an LTR vetoes it as coding
  sequence, and core hits decide the nesting test.  Other LTR profiles are
  accessory evidence.

HMMER and nail read the file unchanged.  A profile file with no `#` class
anywhere (LTR_retriever's TEfam.hmm) gets LTR_retriever's own lists from
`annotate_TE.pl`, which are built into the code.

## profiles.tsv

| column | meaning |
|---|---|
| accession | Pfam accession without its version |
| class | RepeatMasker class (`LTR`, `LTR/Gypsy`, `DNA/hAT`, `LINE`, `Other` ...) |
| role | `core` or `accessory` for LTR profiles, `-` otherwise |
| status | `kept` (in LTR_retriever's TEfam.hmm) or `added` |
| evidence | where the class comes from |
| source | optional: the profile file to take this profile from instead of the Pfam release, as a path relative to `data/` (`TEfam.hmm` for the six 2015 builds) |

A line that starts with `#` is skipped by the build.  A removed profile
stays in the table as such a line, with the status `removed` and the
reason.

## How the classes were set (2026-11)

- **Kept profiles** are TEfam.hmm's 102 Pfam accessions.  100 are in the
  file: 94 rebuilt from Pfam 39.0 and 6 as the 2015 builds; DUF390 was
  retired by Pfam and Atrophin-1 was removed.  Their order (LTR or not)
  and core role are LTR_retriever's.
- **Subtypes come from RepeatPeps.**  Each profile was searched with
  hmmsearch (domain E <= 1e-5) against RepeatPeps, the class-labelled TE
  proteins in `RepeatModeler/Libraries/PROT-Dfam_4.0`.  A subtype is used
  when at least 90% of the proteins hit agree (at least 3 proteins), and
  the type alone when they don't.  This is a labelling aid only; no
  RepeatPeps sequence is in the profiles.  It changes LTR_retriever's
  superfamily lists in a few places:
  - `rve` (integrase) is split Gypsy 1,816 / Copia 1,787, so it is now
    plain `LTR`.  LTR_retriever counted it as Gypsy.
  - `RVP` and `zf-H2C2` are plain `LTR` (mostly ERV, or under 90% Gypsy).
  - `zf-CCHC_5` is `LTR/ERVK`, `DUF1758` is `LTR/Pao` and `Peptidase_A3`
    is `LTR/Caulimovirus`.  They were counted as Copia, Gypsy and Gypsy.
  - Not-LTR profiles take their RepeatPeps class (`DNA/MULE-MuDR`, `LINE`,
    `RC/Helitron` ...).  `Other` means no TE class: a host or unknown
    domain, still counted as evidence against.
  - Profiles with no RepeatPeps hits keep LTR_retriever's class.  Several
    LTR accessories are in this group: DUF834, HGWP, DUF1263, eIF3g and
    others.  Some of them (DUF834, HGWP) are the only domain evidence for
    many rice elements, so they stay until a benchmark says otherwise.
    Atrophin-1 was in this group until it was removed (see "Removed").
- **Added profiles** are retrotransposon core domains missing from
  TEfam.hmm, whose RepeatPeps hits are LTR only:
  - Gypsy and retroviral RNase H (`RT_RNaseH`, `RT_RNaseH_2`)
  - integrase zinc-binding and C-terminal domains (`Integrase_H2C2`,
    `Integrase_Zn`, `IN_DBD_C`, `Integrase_p58_C`)
  - RT connection and thumb
  - Ty3, Ty5 and retroviral capsid, matrix and gag (`Ty3_capsid`,
    `Retrotran_gag_4`, `Gag_p24`, `Gag_p30`, `Gag_MA`, `Gag_p10` ...)
  - Pao RNase H (`Pao_retrotransp`, and `DUF1759`, which LTR_retriever
    listed as Copia)
  - env proteins, as accessory evidence only
- **Retired by Pfam:** DUF390 (PF04094) is in neither Pfam 38.2 nor 39.0,
  and nothing replaces it: its consensus is mostly low-complexity, and the
  best Pfam 39 match is E = 6e-4.  It was the only domain hit for 2
  candidates across the six tuning genomes, so it is dropped.
- **Left out** because they are shared with LINEs, DNA transposons or host
  genes:
  - RT `RVT_1` (46% of its hits are LTR) and `RVT_3` (35%)
  - `RNase_H` (43%)
  - `PMD` (MULE and Gypsy)
  - `Phage_integrase` (DIRS and Crypton)
  - `zf-CCHC` (a generic zinc knuckle)

## Kept as the 2015 builds

Six of the kept profiles are less sensitive in Pfam 39.0 than in the 2015
builds that LTR_retriever's TEfam.hmm holds.  Pfam re-seeded them from
fewer sequences, and the new Chromo no longer finds the chromodomains of
chromoviruses.  A file built from Pfam 39.0 alone left 344 kept rice
elements with no LTR domain evidence that had some with TEfam.hmm.  So
these six are copied from `../TEfam.hmm` (the `source` column), with the
class and role written into them like the others:

| profile | 2015 build (seeds) | Pfam 39.0 build (seeds) | candidates hit, 2015 build | candidates hit, Pfam 39.0 build |
|---|---|---|---|---|
| Retrotrans_gag | PF03732.14 (56) | PF03732.26 (9) | rice 2,976 | rice 1,174 |
| DUF834 | PF05754.11 (78) | PF05754.20 (13) | rice 1,830 | rice 363 |
| Chromo | PF00385.21 (121) | PF00385.32 (80) | rice 1,160 | rice 17 |
| Asp_protease_2 | PF13650.3 (159) | PF13650.13 (103) | rice 1,026 | rice 308 |
| DUF4219 | PF13961.3 (67) | PF13961.14 (21) | rice 384 | rice 68 |
| zf-CCHC_5 | PF14787.3 (5) | PF14787.13 (4) | mouse 2,103 | mouse 1,829 |

The rule for taking the 2015 build: it hits at least 10% more candidates
than the Pfam 39.0 build, and at least 20 more, in rice or in mouse.  Every
other kept profile is the Pfam 39.0 build.

With the six, 11 kept rice elements lose their LTR domain evidence instead
of 344, and 3,552 kept rice elements have a core hit (3,537 with TEfam.hmm,
3,374 with Pfam 39.0 alone).  Library scores against the curated
annotation, as sensitivity / precision and library records:

| genome | TEfam.hmm | Pfam 39.0 alone | with the six 2015 builds |
|---|---|---|---|
| Arabidopsis | 90.91% / 96.06%, 408 | 89.59% / 95.71%, 395 | 91.09% / 95.50%, 412 |
| rice | 93.47% / 98.45%, 912 | 93.76% / 98.41%, 920 | 93.57% / 98.42%, 929 |
| zebra finch | 74.61% / 86.51%, 198 | 76.13% / 86.28%, 216 | 76.13% / 86.28%, 216 (same library) |

These numbers were measured on 2026-11-15 (`runs/domains26`), with
Atrophin-1 still in all three files.

When a later Pfam release rebuilds one of the six, compare the two builds
by the same rule before switching.

## Removed

**Atrophin-1 (PF03154), removed 2026-11-17.**  LTR_retriever lists it as
an LTR domain.  It is a host protein family (Pfam 39.0 builds it from 3
seed sequences, 992 positions), it has no hits in RepeatPeps, and its
hits on LTR candidates are weak.  As an LTR accessory it was enough to
keep an element by itself:

| genome | candidates with a hit | kept elements with a hit | kept, Atrophin-1 their only domain hit | of those, kept on it alone | E value of the only-hit kept elements (best / median / worst) |
|---|---|---|---|---|---|
| Arabidopsis | 2 | 1 | 1 | 1 | 5.2e-4 |
| rice | 3 | 3 | 2 | 2 | 9.4e-3 / 4.7e-2 / 4.7e-2 |
| zebra finch | 1 | 1 | 1 | 1 | 1.3e-3 |
| mouse (mm10) | 702 | 78 | 71 | 64 | 4.3e-7 / 7.3e-3 / 4.8e-2 |
| gorilla | 63 | 6 | 5 | 5 | 3.1e-3 / 3.8e-3 / 2.0e-2 |

"Kept on it alone" means no TSD or TG..CA evidence either.  In the curated
annotation, 46 of the 71 mouse elements have no LTR annotation on either
LTR, 23 have LTR annotation on both and 2 on one; of the 5 gorilla
elements, 4 have none and 1 has both.  One of the mouse elements is a pair
of L1 copies (chr1:21148651..21151309, hit E = 7.3e-3) whose internal
region became a library model that masked 10.5 Mb of L1 sequence.

The cost: the 23 mouse elements with LTR annotation on both LTRs lose
their only domain hit, and 19 of them had no other evidence.

The row is kept in `profiles.tsv` as a comment.  The lists built into
`src/domain.rs` for a profile file without classes (TEfam.hmm) are
LTR_retriever's and still count PF03154 as LTR.

## Future upgrades

- **Profiles built from RepeatPeps.**  Per-domain, per-class profiles would
  give RT and RNase H profiles that tell LTR from LINE, and name ERV1, ERVK,
  ERVL, Pao, DIRS, Ngaro and Caulimovirus elements directly.  The labels
  would match RepeatMasker's.
- **REXdb** (Neumann et al. 2019) has domain sequences split by domain
  and by lineage: Ty1/copia and Ty3/gypsy GAG, PROT, RT, RH and INT, the
  plant chromodomains (CHD, CHDCR), and LINE, DIRS, Penelope, Bel-Pao and
  DNA-transposon domains.  Its RT and RH are per class, which Pfam's are
  not.  Licence:
  - The GitHub repository (github.com/repeatexplorer/rexdb) has no licence
    file.
  - Its one tagged release, v1.0 (2023-11-20), is archived on Zenodo
    (doi:10.5281/zenodo.10160279) under CC BY 4.0.  That release holds
    Viridiplantae v3.0 (80,446 sequences) and Metazoa v3.1 (11,192); a
    copy is in `LTRHarvest/data/rexdb/zenodo-10160279/`.
  - CC BY 4.0 allows profiles built from it, with attribution (cite the
    paper and the Zenodo record, and say what was changed).
  - Viridiplantae v4.0 is on the repository's main branch but in no
    release, so it has no stated licence; ask the authors before using it.
  Nothing from REXdb is in TEdomains.hmm yet.
- **Name more superfamilies** in the call (ERV, Pao) once the profiles
  separate them.
