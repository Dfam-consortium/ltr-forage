# Nucleotide scoring matrices

`20p<GC>g.matrix`: RepeatMasker's 20% divergence matrices for background GC 35-53%
(copied from RepeatMasker 4.2.6, `Matrices/ncbi/nt/`).  The recipe picks one by GC
(`pipeline/pick_matrix.sh`) for the rmblastn searches in step 6 (depth
cleanup, against the genome sample) and steps 7-14 (library all-vs-all, `rmb_ava.sh`).
