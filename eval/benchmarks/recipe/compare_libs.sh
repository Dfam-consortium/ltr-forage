#!/bin/bash
# Score candidate libraries: masking vs a standard .out and family recovery vs a curated LTR library.
# usage: compare_libs.sh <genome.fa> <standard.out> <curated_LTR.fa> <outdir> name=lib.fa [name=lib.fa ...]
# Env: PA (RepeatMasker -pa per library, default 12).
set -eo pipefail
G=$1; STD=$2; CUR=$3; OUT=$4; shift 4
B=$(cd "$(dirname "$0")/.." && pwd)
RM="/usr/local/RepeatMasker/RepeatMasker -e rmblast -no_is -norna -nolow"
mkdir -p $OUT; cd $OUT
for kv in "$@"; do
  n=${kv%%=*}; src=$(realpath ${kv#*=}); [ "$src" = "$(realpath -m $n.fa)" ] || cp $src $n.fa
  ( mkdir -p mask_$n; $RM -pa ${PA:-12} -lib $n.fa -dir mask_$n $G > mask_$n/log 2>&1
    python3 $B/mask_benchmark.py $G.fai $STD mask_$n/$(basename $G).out --test-class all > mask_$n/metric.txt
    python3 $B/family_recovery.py $CUR $n.fa | sed -n 3p > famrec_$n.txt ) &
done
wait
for kv in "$@"; do n=${kv%%=*}; echo "== $n $(grep -c '>' $n.fa) seqs"; cat famrec_$n.txt; grep -E "sensitivity|^genome|^TP" mask_$n/metric.txt; done
