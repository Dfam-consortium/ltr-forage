#!/bin/bash
# Drop-in for trf as called by LTR_retriever's cleanup.pl: split the FASTA into 48
# chunks, run trf on each in parallel, concatenate the -ngs output to stdout.
f=$1; shift
d=$(mktemp -d "${TMPDIR:-/tmp}/trfpar.XXXXXX")
awk -v d="$d" 'BEGIN{n=0} /^>/{if (n % 900 == 0) {c=sprintf("%s/c%04d.fa", d, int(n/900))} n++} {print > c}' "$f"
ls "$d"/c*.fa | xargs -P 48 -I{} sh -c '/usr/local/bin/trf {} '"$*"' > {}.out 2>/dev/null; true'
cat "$d"/c*.fa.out
rm -rf "$d"
