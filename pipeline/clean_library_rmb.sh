#!/bin/bash
# Depth cleanup with the rmblast port (ltr_library.py step 6; eval/benchmarks/clean_library.sh is the RepeatMasker version).
# Usage: clean_library_rmb.sh <lib.fa> <genome_sample.fa> <outdir> <rounds> <scale> <matrix> [excise options...]
# Each round: rmblastn of the sample with RepeatMasker's -lib settings (NCBIBlastSearchEngine.pm:
# word 9, gap 24/6, X-drops 450/112/225, minscore 225; matrix = 20p<GC>g for the genome's GC),
# except that the port applies mask level 80 itself instead of RepeatMasker's 101 + Perl filter; split into
# BATCHES query batches run in parallel, against the current library; then excise by hit depth.
# Env: RMBLASTN (port binary), BATCHES (16), THREADS per batch (4; BATCHES x THREADS = cores, as RepeatMasker -pa 16), ROUND1_HITS (reuse for round 1).
set -euo pipefail
lib=$1; sample=$2; out=$3; rounds=$4; scale=$5; matrix=$6; shift 6
HERE=$(cd "$(dirname "$0")" && pwd)
# default: a rmblast-port checkout beside the LTRHarvest directory (ltr_library.py always sets RMBLASTN)
RMB=${RMBLASTN:-$HERE/../../../rmblast-port/rmblast/target/release/rmblastn}
B=${BATCHES:-16}; T=${THREADS:-4}
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
pids=()
for r in $(seq 1 "$rounds"); do
  hits="$out/hits$r.tsv"
  if [ "$r" = 1 ] && [ -n "${ROUND1_HITS:-}" ]; then cp "$ROUND1_HITS" "$hits"
  else
    for i in $(seq 0 $((B-1))); do
      "$RMB" --query "$out/batches/q$i.fa" --db "$out/round$((r-1)).fa" --matrix "$matrix" --gapopen 24 --gapextend 6 \
        --word-size 9 --xdrop-ungap 450 --xdrop-gap 112 --xdrop-gap-final 225 --min-raw-gapped-score 225 \
        --complexity-adjust --dust no --mask-level 80 --num-threads "$T" --mt-mode 0 \
        --outfmt "6 score qseqid qstart qend sseqid sstart send sstrand" > "$out/batches/h$i.tsv" &
      pids+=($!)
    done
    for p in "${pids[@]}"; do wait "$p" || { echo "clean_library_rmb.sh: a rmblastn batch failed" >&2; exit 1; }; done
    pids=()
    for i in $(seq 0 $((B-1))); do cat "$out/batches/h$i.tsv"; done > "$hits"
  fi
  python3 "$HERE/hits_to_profile.py" "$out/round$((r-1)).fa" "$hits" > "$out/depth$r.tsv"
  python3 "$HERE/excise_by_depth.py" "$out/round$((r-1)).fa" "$out/depth$r.tsv" "$out/round$r.fa" --profile --scale "$scale" "$@" 2>&1 | tee "$out/round$r.log"
done
cp "$out/round$rounds.fa" "$out/clean.fa"
