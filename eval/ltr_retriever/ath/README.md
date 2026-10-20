Arabidopsis run (2026-09-06) against the RepeatMasker run now in `data/annotations/ath-TAIR10.1`
(TAIR10.1, RepeatMasker 4.2.4, Dfam 4.0 RBRM).  Steps: upper-case the
soft-masked NCBI FASTA, `samtools faidx`, `sufr create --dna`, `ltr-forage`
at defaults and at LTR_retriever's recommended settings, the re-extension
scripts in ../reextend (windows, pairwise rmblastn, HSP, 6-of-6 trim),
LTR_retriever 2.9.0 on the default candidates with the parallel TRF wrapper,
and ath_features.py for TSD/motif/internal-gap features (it imports
../truth_loose.py as truth2 for the loose harness).  Results: section 9 of
ltrstruct_without_ltr_retriever.md.
