Signal analysis (2026-09-06): which per-candidate signals separate real
LTR/int/LTR elements from the rest of the ltr-forage candidate set on hg38.

    K=20 python3 features.py            # features20.tsv, kmers_int20.txt, kmers_ltr20.txt
    # add reverse complements, then:
    sufr count hg38.sufr kmers_int20_both.txt > counts_int20.txt
    sufr count hg38.sufr kmers_ltr20_both.txt > counts_ltr20.txt
    python3 analyze.py

features.py expects the re-extended predictions (../reext/reext_hsp.gt) and the
LTR_retriever .defalse from the run on those coordinates, for the domain calls.
Results are in section 8 of ltrstruct_without_ltr_retriever.md.

Coverage along the internal region (coverage.py): tile each internal region with
20-mers at a 10 bp stride, count the tiles with `sufr count` against the genome
index and against an index of the candidate internals (`sufr create --dna` on
internals.fa), then write per-candidate coverage, per-tenth minimum, longest
unique gap and count uniformity to coverage_stats.tsv. The longest unique gap
(<= 10% of the region) is the useful statistic; see section 8 of the doc.

truth_loose.py (one directory up) is the loosened truth harness: LTR blocks of any
subfamily name, adjacent LTR fragments merged (gap <= 100 bp), the span between two
blocks >= 50% covered by LTR-internal rows on the same strand, 5'-5' distance 1-20 kb,
mutual-nearest pairing.  `python3 truth_loose.py hg38.out hg38.fa.fai preds.gt ...`
