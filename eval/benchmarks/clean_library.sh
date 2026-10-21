#!/bin/bash
# Iterative nested-insertion cleanup of a library against a genome sample.
# Usage: clean_library.sh <lib.fa> <genome_sample.fa> <outdir> <rounds> <scale> [excise options...]
# Each round: RepeatMask the sample with the library, excise by hit depth, repeat.
set -euo pipefail
lib=$1; sample=$2; out=$3; rounds=$4; scale=$5; shift 5
RM=${RM:-/usr/local/RepeatMasker/RepeatMasker}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out"; cp "$lib" "$out/round0.fa"
for r in $(seq 1 "$rounds"); do
  d="$out/mask$r"; mkdir -p "$d"
  $RM -pa 32 -e rmblast -no_is -norna -nolow -lib "$out/round$((r-1)).fa" -dir "$d" "$sample" > "$d/rmsk.log" 2>&1
  python3 "$HERE/../../pipeline/excise_by_depth.py" "$out/round$((r-1)).fa" "$d/$(basename "$sample").out" "$out/round$r.fa" --scale "$scale" "$@" 2>&1 | tee "$out/round$r.log"
done
cp "$out/round$rounds.fa" "$out/clean.fa"
