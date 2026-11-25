# Example: Arabidopsis thaliana chromosome 4

`ath_chr4.fa.gz` is chromosome 4 of the TAIR10.1 assembly (NCBI RefSeq
NC_003075.7, assembly GCF_000001735.4): 18.6 Mb, with a pericentromere full
of Gypsy and Copia elements.  It is small enough to run in a few minutes and
large enough to give a few dozen families.

The commands are in the main README's
[Example: Arabidopsis chromosome 4](../README.md#example-arabidopsis-chromosome-4)
section.  In short:

```bash
gzip -dc examples/ath_chr4.fa.gz > ath_chr4.fa
sufr create -t 8 --dna --write-lcp -o ath_chr4.sufr ath_chr4.fa
ltr-forage ath_chr4.sufr -t 8 --library elements.fa --fastan-out tandem.bed --log ltr-forage.log > cand.tsv
pipeline/ltr_library.py ath_chr4.fa -o out -t 8 --evidence-page
```

With ltr-forage 0.1.0 at 8 threads (2026-11-23): the index takes a few
seconds; ltr-forage finds 823 candidates and keeps 116 elements in under a
minute; the library run gives 47 families (47 LTR and 47 internal models) in
about 4 minutes, 22 of them "Well supported" on the evidence page.
`ltr-forage --trna-scan ath_chr4.fa` finds 100 tRNA genes; the release
workflow checks that number.
