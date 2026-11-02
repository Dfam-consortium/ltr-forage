#!/bin/bash
# Reproduce the hg38 comparison between GenomeTools LTRharvest and ltr-forage
# described in the README.  Every step is skipped when its output exists, so
# the script can be re-run after a failure.
#
# Requires: gt (GenomeTools 1.6.2), sufr (0.7.12), samtools, python3,
# /usr/bin/time, and a RepeatMasker .out annotation of the same assembly
# for the last two steps.
set -euo pipefail

GENOME=${GENOME:-hg38.fa}          # input FASTA
RM_OUT=${RM_OUT:-hg38.fa.out}      # RepeatMasker annotation of GENOME
THREADS=${THREADS:-32}             # ltr-forage -j
GT=${GT:-gt}
SUFR=${SUFR:-sufr}
LTR_FORAGE=${LTR_FORAGE:-ltr-forage}
EVAL=$(cd "$(dirname "$0")" && pwd)
TIME="/usr/bin/time -v"

PARAMS_GT="-seed 30 -minlenltr 100 -maxlenltr 1000 -mindistltr 1000 -maxdistltr 15000 \
  -similar 85 -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -mintsd 4 -maxtsd 20 -vic 60 -overlaps best"
PARAMS_PORT="--min-seed 30 --ltr-min 100 --ltr-max 1000 --dist-min 1000 --dist-max 15000 \
  --similar 85 --xdrop 5 --mat 2 --mis -2 --ins -3 --del -3 --mintsd 4 --maxtsd 20 --vicinity 60 --overlaps best --ltrharvest"

base=$(basename "$GENOME" .fa)

echo "== FASTA index"
[ -s "$GENOME.fai" ] || samtools faidx "$GENOME"

echo "== GenomeTools index (single-threaded)"
[ -s "${base}_esa.suf" ] || $TIME $GT suffixerator -dna -pl -tis -suf -lcp -sds -des \
    -db "$GENOME" -indexname "${base}_esa" 2> "${base}_gt_index.time"

echo "== sufr index (do not add --ignore-softmask: it turns lowercase into N)"
[ -s "$base.sufr" ] || $TIME $SUFR create --dna "$GENOME" -o "$base.sufr" 2> "${base}_sufr_index.time"

echo "== LTRharvest search (single-threaded; gt -j is slower and not deterministic)"
[ -s "${base}_gt.out" ] || $TIME $GT ltrharvest -index "${base}_esa" $PARAMS_GT \
    > "${base}_gt.out" 2> "${base}_gt_search.time"

echo "== ltr-forage search, 1 thread and $THREADS threads"
[ -s "${base}_ltr-forage_j1.out" ] || $TIME $LTR_FORAGE "$base.sufr" $PARAMS_PORT -j 1 \
    > "${base}_ltr-forage_j1.out" 2> "${base}_ltr-forage_j1.time"
[ -s "${base}_ltr-forage.out" ] || $TIME $LTR_FORAGE "$base.sufr" $PARAMS_PORT -j "$THREADS" \
    > "${base}_ltr-forage.out" 2> "${base}_ltr-forage_j$THREADS.time"
cmp "${base}_ltr-forage_j1.out" "${base}_ltr-forage.out" && echo "ltr-forage output identical at 1 and $THREADS threads"

echo "== compare predictions"
python3 "$EVAL/compare_predictions.py" "${base}_gt.out" "${base}_ltr-forage.out" \
    --fasta "$GENOME" --out-prefix "${base}cmp" | tee "${base}_compare.txt"

if [ -s "$RM_OUT" ]; then
    echo "== recall against intact RepeatMasker elements"
    python3 "$EVAL/rm_truth_recall.py" "$RM_OUT" "$GENOME.fai" "${base}_gt.out" "${base}_ltr-forage.out" \
        "${base}cmp.a_only.tsv" "${base}cmp.b_only.tsv" | tee "${base}_recall.txt"
    echo "== RepeatMasker class under each prediction's LTRs"
    for f in "${base}_gt.out" "${base}_ltr-forage.out"; do
        python3 "$EVAL/rm_classify.py" "$RM_OUT" "$GENOME.fai" "$f"
    done | tee "${base}_classes.txt"
else
    echo "no RepeatMasker annotation at $RM_OUT; skipping recall and classification"
fi

echo "== timings (wall clock, peak RSS)"
for t in "${base}"_*.time; do
    [ -e "$t" ] || continue
    printf '%-32s %s  %s\n' "$t" "$(grep 'Elapsed' "$t" | sed 's/.*: //')" "$(grep 'Maximum resident' "$t" | sed 's/.*: //') KB"
done
