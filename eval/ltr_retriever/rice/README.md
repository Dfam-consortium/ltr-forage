Rice run (2026-09-06) against the RepeatMasker run now in `data/annotations/rice-ASM3414082v1` (ASM3414082v1,
RepeatMasker 4.2.4, rice 7.0.0 TE library).  Same steps as ../ath, with two
changes: the internal-name test is "contains _INT" (truth_rice.py is truth.py
with that rule; truth_loose.py takes the same patch), and the re-extension
windows are sized for a 2 kb LTR (LTR_MAX in cut_windows.py), because 5 kb
windows ran into neighbouring nested elements.  Results: sections 10 and 11
of ltrstruct_without_ltr_retriever.md.
