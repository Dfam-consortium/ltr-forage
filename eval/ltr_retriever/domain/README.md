# Domain search validation

`score_tsv.py` scores an `ltr-forage --format tsv` table against RepeatMasker
truth (the strict and loose harnesses from the parent directory) and prints the
keep rates of the domain columns and the candidate rules. With a
`signals/features.py` table as a third argument it also joins LTR_retriever's
hmmsearch order call to each candidate and prints both callers side by side.

```
ltr-forage hg38.sufr --reextend --internal-stats --format tsv \
    --domain-hmm /usr/local/LTR_retriever-2.9.0/database/TEfam.hmm > hg38.tsv
python3 score_tsv.py ~/hg38.out hg38.tsv signals/features.tsv
```

`cmp_domains.py` compares the `dom_hits` column with an hmmsearch `--tblout`
file on the same candidates (LTR_retriever's `*.extend.fa.aa.tbl`), matching
candidates by 50% overlap and listing accessions found by one caller only.

`self_eval.py` scores the LTR self-alignment columns: it takes the
RepeatMasker `.out`, a table with the `lltr_self_cov` columns, and the
domain-search table of the same run (joined by element coordinates), and
prints how many real and non-real kept candidates each `--tandem-max`
threshold would remove, with the annotation of the flagged ones.
