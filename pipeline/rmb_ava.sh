#!/bin/bash
# Library-vs-library search with the rmblast port, written as the BLAST+ tabular columns the
# recipe scripts read (drop-in for blastn -task dc-megablast -outfmt "6 <fields>").
# usage: rmb_ava.sh <query.fa> <db.fa> <max_hsps> <out.tsv> <fields...>
#   fields: any of qseqid sseqid pident length qstart qend sstart send qlen slen bitscore
# Env: THREADS (32), RMBLASTN (port binary), WORD (seed word size, 12; 9-16 give the same 80-80 pairs and scores, see runs/aligner_rmb_word), MATRIX (default 20p<GC>g for the db's GC, pick_matrix.sh),
#      CHUNK (2000 queries per rmblastn call: --mt-mode 1 holds a call's results in memory until it ends),
#      MT_MODE (1: split by queries, streamed in input order, the default; rmblastn 3.0.10's --mt-mode 1 can
#      deadlock, use 3.0.9 or >= 3.0.11, or MT_MODE=0).
# Search: RepeatMasker's -lib settings except the word (gap 24/6, X-drops 450/112/225, minscore 225),
# no mask level (a self-comparison would hide every hit behind the self-hit), E <= 1e-10.
# pident/length are BLAST's: length = alignment columns = query span + query gaps,
# pident = matches / columns (matches = columns - substitutions - all gaps).
# At most max_hsps hits per query-subject pair, best first (the port reports them by score).
set -eo pipefail
Q=$1; DB=$2; MAXH=$3; OUT=$4; shift 4
# nothing to search (e.g. no multi-element clusters for the flank filter): no hits
if ! grep -q "^>" "$Q" || ! grep -q "^>" "$DB"; then : > "$OUT"; : > "$OUT.log"; exit 0; fi
RMB=${RMBLASTN:-$(dirname "$0")/../../../rmblast-port/rmblast/target/release/rmblastn}   # rmblast-port beside LTRHarvest
MATRIX=${MATRIX:-$(bash "$(dirname "$0")/pick_matrix.sh" "$DB")}
CH=$(mktemp -d "$OUT.chunks.XXXX")
awk -v n=${CHUNK:-2000} -v d="$CH" '/^>/{if (k++ % n == 0) f = sprintf("%s/q%05d.fa", d, k / n)} {print > f}' "$Q"
: > "$OUT"; : > "$OUT.log"
for qf in "$CH"/q*.fa; do
/usr/bin/time -f "chunk $(basename $qf) %e s %M KB" -a -o "$OUT.log" \
"$RMB" --query "$qf" --db "$DB" --matrix "$MATRIX" --gapopen 24 --gapextend 6 --word-size ${WORD:-12} \
  --xdrop-ungap 450 --xdrop-gap 112 --xdrop-gap-final 225 --min-raw-gapped-score 225 --evalue 1e-10 \
  --complexity-adjust --dust no --mask-level -1 --num-threads ${THREADS:-32} --mt-mode ${MT_MODE:-1} \
  --outfmt "6 qseqid sseqid qstart qend sstart send qlen slen perc_sub perc_query_gap perc_db_gap bitscore" 2>> "$OUT.log" |
awk -F'\t' -v OFS='\t' -v maxh="$MAXH" -v fields="$*" '
BEGIN { nf = split(fields, F, " ") }
$1 != cq { cq = $1; delete seen }
++seen[$2] <= maxh {
  qspan = $4 - $3 + 1; sspan = ($5 > $6 ? $5 - $6 : $6 - $5) + 1
  qg = int($10 * qspan / 100 + 0.5); sg = int($11 * sspan / 100 + 0.5); x = int($9 * qspan / 100 + 0.5)
  len = qspan + qg; v["length"] = len; v["pident"] = sprintf("%.3f", 100 * (len - x - qg - sg) / len)
  v["qseqid"] = $1; v["sseqid"] = $2; v["qstart"] = $3; v["qend"] = $4; v["sstart"] = $5; v["send"] = $6
  v["qlen"] = $7; v["slen"] = $8; v["bitscore"] = $12
  line = v[F[1]]; for (i = 2; i <= nf; i++) line = line OFS v[F[i]]
  print line
}' >> "$OUT"
done
rm -r "$CH"
