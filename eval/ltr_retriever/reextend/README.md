Re-extension experiment (2026-09-06). Takes ltr-forage's default predictions,
cuts a window around each LTR wide enough for a 2 kb LTR plus 100 bp, aligns the
left window against the right window with the rmblastn port (one pairwise run
per candidate, RepeatMasker matrix, gap 20/5, word 14), and uses the best HSP as
the re-extended LTR pair. Optional trims to a TSD or TG..CA motif are in
reextend.py (RANK and VIC environment variables); none beat the raw HSP ends.

    python3 cut_windows.py preds.gt genome.fa.fai genome.fa 50
    for i in ...; do rmblastn --db pairs/$i.s.fa --query pairs/$i.q.fa ... > pairs/$i.hits; done
    cat pairs/*.hits > all.hits
    python3 reextend.py all.hits windows.tsv genome.fa.fai genome.fa hsp out.gt

The rmblastn invocation and the results are in ltrstruct_without_ltr_retriever.md.

trim_ends.py: alignment-intrinsic end trimming.  Needs the pairwise hits with
qseq/sseq in the outfmt; walks in from each end until W aligned columns are
gap-free with >= M matches.  `python3 trim_ends.py all_aln.hits windows.tsv out.gt 6 6`
