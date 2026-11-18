#!/bin/bash
# Depth cleanup with the rmblast port (ltr_library.py step 6; eval/benchmarks/clean_library.sh is the RepeatMasker version).
# Usage: clean_library_rmb.sh <lib.fa> <genome_sample.fa> <outdir> <rounds> <scale> <matrix> [excise options...]
# Each round: rmblastn of the sample with RepeatMasker's -lib settings (NCBIBlastSearchEngine.pm:
# word 9, gap 24/6, X-drops 450/112/225, minscore 225; matrix = 20p<GC>g for the genome's GC),
# except that the port applies mask level 80 itself instead of RepeatMasker's 101 + Perl filter; split into
# BATCHES query batches run in parallel, against the current library; then excise by hit depth.
# Env: RMBLASTN (the Rust rmblastn; default /usr/local/rmblast/bin/rmblastn), BATCHES (16), FILTER_JOBS (mask-level filters at a time, 4), THREADS per batch (4; BATCHES x THREADS = cores, as RepeatMasker -pa 16), ROUND1_HITS (reuse for round 1; round 2 then searches the whole library), INCREMENTAL (1), REMAP (0; see below), PREV_RUN (see below), DEPTH_WORD (seed word size, 9 as RepeatMasker).
set -euo pipefail
lib=$1; sample=$2; out=$3; rounds=$4; scale=$5; matrix=$6; shift 6
HERE=$(cd "$(dirname "$0")" && pwd)
# default: the installed release, never a development checkout, which a rebuild swaps mid-run
# (ltr_library.py always sets RMBLASTN)
RMB=${RMBLASTN:-/usr/local/rmblast/bin/rmblastn}
B=${BATCHES:-16}; T=${THREADS:-4}
# a failed batch ends the script: stop the batches still running with it
trap 'p=$(jobs -p); [ -z "$p" ] || kill $p 2>/dev/null || true' EXIT
# a fresh batch directory: a rerun with fewer batches would otherwise concatenate the earlier run's extra h*.tsv
# (mm10, 2026-10-25: a 16-batch rerun of a 24-batch step 22b counted a third of the sample twice)
rm -rf "$out/batches"; mkdir -p "$out/batches"; cp "$lib" "$out/round0.fa"
# split the sample into B batches of whole sequences, balanced by length
python3 - "$sample" "$out/batches" "$B" <<'PY'
import sys
src, d, B = sys.argv[1], sys.argv[2], int(sys.argv[3])
recs = []; 
for l in open(src):
    if l.startswith('>'): recs.append([l, []])
    else: recs[-1][1].append(l)
size = [0] * B; fh = [open(f'{d}/q{i}.fa', 'w') for i in range(B)]
for h, s in sorted(recs, key=lambda r: -sum(len(x) for x in r[1])):
    i = size.index(min(size)); fh[i].write(h + ''.join(s)); size[i] += sum(len(x) for x in s)
PY
# Incremental rounds (INCREMENTAL=1, the default): every round keeps its unfiltered hits (--mask-level 100,
# batches/raw<i>.tsv).  A record's hits do not depend on the other library records except through the mask-level
# filter, so from round 2 on only the records the last round changed are searched; the unchanged records' hits
# carry over and masklevel.py applies mask level 80 to the union - the same hits as a search of the whole library
# (rice, 2026-11-09: identical; that round 1,121 -> 220 s).  INCREMENTAL=0 searches the whole library every round.
INC=${INCREMENTAL:-1}; raw=0
# PREV_RUN=<an earlier run's outdir, same sample and BATCHES>: round 1 starts from that run's last unfiltered hits
# (its round<last-1>.fa library) as a later round would (step 22b's second depth search of the trimmed models)
prev_db=""
if [ "$INC" = 1 ] && [ -n "${PREV_RUN:-}" ]; then
  prev_db=$(ls "$PREV_RUN"/round*.fa | sed 's/.*round\([0-9]*\)\.fa/\1/' | sort -n | tail -2 | head -1)
  prev_db="$PREV_RUN/round$prev_db.fa"
  for i in $(seq 0 $((B-1))); do
    cmp -s "$PREV_RUN/batches/q$i.fa" "$out/batches/q$i.fa" && [ -f "$PREV_RUN/batches/raw$i.tsv" ] || { prev_db=""; break; }
  done
  if [ -n "$prev_db" ]; then for i in $(seq 0 $((B-1))); do cp "$PREV_RUN/batches/raw$i.tsv" "$out/batches/raw$i.tsv"; done; raw=1
  else echo "clean_library_rmb.sh: $PREV_RUN has other batches or no unfiltered hits; searching the whole library" >&2; fi
fi
# REMAP=1 (with INCREMENTAL=1; ltr_library.py sets it for step 14 unless --no-depth-remap): rounds after the first are not searched.  A round's cuts are applied to the last
# round's unfiltered hits instead (masklevel.py remap): a hit on a dropped record goes, a hit on a cut record keeps
# the alignment columns outside the cuts, rescored as rmblastn scores them, each cut end trimmed or extended across
# the junction, and two pieces either side of a cut joined when they are collinear on the query (a copy without the
# insertion that was cut).  For that the unfiltered hits carry their alignment rows (about 4 times the size).
# Not exact: rice 9 of 862 final records differ (531 bp), 93.42 / 98.42 against 93.40 / 98.43; mm10 31 of 353
# (1,826 bp), 42.88 / 92.89 against 42.93 / 92.90; rounds 2 and 3 take 0.5 min instead of 6.8 and 6.0 instead of 51.8.
ALN=""; if [ "$INC" = 1 ] && [ "${REMAP:-0}" = 1 ] && [ "$rounds" -gt 1 ]; then ALN=" qseq sseq"; fi
run_batches() {  # db mask_level outfile-prefix: rmblastn of every query batch against db, in parallel
  local db=$1 ml=$2 pre=$3 pids=() i p
  for i in $(seq 0 $((B-1))); do
    "$RMB" --query "$out/batches/q$i.fa" --db "$db" --matrix "$matrix" --gapopen 24 --gapextend 6 \
      --word-size ${DEPTH_WORD:-9} --xdrop-ungap 450 --xdrop-gap 112 --xdrop-gap-final 225 --min-raw-gapped-score 225 \
      --complexity-adjust --dust no --mask-level "$ml" --num-threads "$T" --mt-mode 0 \
      --outfmt "6 score qseqid qstart qend sseqid sstart send sstrand$ALN" > "$out/batches/$pre$i.tsv" &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p" || { echo "clean_library_rmb.sh: a rmblastn batch failed" >&2; exit 1; }; done
}
for r in $(seq 1 "$rounds"); do
  hits="$out/hits$r.tsv"; cur="$out/round$((r-1)).fa"
  if [ "$r" = 1 ] && [ -n "${ROUND1_HITS:-}" ]; then cp "$ROUND1_HITS" "$hits"
  elif [ "$INC" = 1 ]; then
    if [ "$raw" = 1 ] && [ "$r" -gt 1 ] && [ -n "$ALN" ]; then  # no search: the last round's hits moved onto the cut records
      pids=()
      for i in $(seq 0 $((B-1))); do
        python3 "$HERE/masklevel.py" remap "$out/cuts$((r-1)).tsv" "$out/batches/raw$i.tsv" 30 225 1 "$matrix" "$out/batches/q$i.fa" "$cur" \
          > "$out/batches/raw$i.next" 2>> "$out/remap$r.log" &
        pids+=($!)
      done
      for p in "${pids[@]}"; do wait "$p" || { echo "clean_library_rmb.sh: masklevel.py remap failed" >&2; exit 1; }; done
      for i in $(seq 0 $((B-1))); do mv "$out/batches/raw$i.next" "$out/batches/raw$i.tsv"; done
    elif [ "$raw" = 1 ]; then  # search the changed records only, carry the rest
      prev="$out/round$((r-2)).fa"; [ "$r" = 1 ] && prev=$prev_db
      python3 "$HERE/masklevel.py" changed "$prev" "$cur" "$out/batches/changed.fa" 2>&1 | sed "s/^/round $r: /"
      if grep -q '^>' "$out/batches/changed.fa"; then run_batches "$out/batches/changed.fa" 100 n
      else for i in $(seq 0 $((B-1))); do : > "$out/batches/n$i.tsv"; done; fi
      for i in $(seq 0 $((B-1))); do
        { python3 "$HERE/masklevel.py" carry "$prev" "$cur" "$out/batches/raw$i.tsv"; cat "$out/batches/n$i.tsv"; } > "$out/batches/raw$i.next"
        mv "$out/batches/raw$i.next" "$out/batches/raw$i.tsv"
      done
    else
      run_batches "$cur" 100 raw
    fi
    # FILTER_JOBS filters at a time (default 4): each reads a whole raw batch (gorilla 2026-11-16: 28 GB each)
    raw=1; FJ=${FILTER_JOBS:-4}
    for w in $(seq 0 "$FJ" $((B-1))); do
      pids=()
      for i in $(seq "$w" $((w+FJ-1 < B-1 ? w+FJ-1 : B-1))); do
        python3 "$HERE/masklevel.py" filter "$out/batches/raw$i.tsv" "$cur" "$out/batches/q$i.fa" 80 > "$out/batches/h$i.tsv" &
        pids+=($!)
      done
      for p in "${pids[@]}"; do wait "$p" || { echo "clean_library_rmb.sh: masklevel.py failed" >&2; exit 1; }; done
    done
    for i in $(seq 0 $((B-1))); do cat "$out/batches/h$i.tsv"; done > "$hits"
  else
    run_batches "$cur" 80 h
    for i in $(seq 0 $((B-1))); do cat "$out/batches/h$i.tsv"; done > "$hits"
  fi
  python3 "$HERE/hits_to_profile.py" "$out/round$((r-1)).fa" "$hits" > "$out/depth$r.tsv"
  python3 "$HERE/excise_by_depth.py" "$out/round$((r-1)).fa" "$out/depth$r.tsv" "$out/round$r.fa" --profile --scale "$scale" --cuts "$out/cuts$r.tsv" "$@" 2>&1 | tee "$out/round$r.log"
done
cp "$out/round$rounds.fa" "$out/clean.fa"
