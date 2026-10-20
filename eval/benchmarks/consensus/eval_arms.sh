#!/bin/bash
# Experiment 2 (2026-10-19): score the cleanup arms written by build_arms.py.
# usage: eval_arms.sh <arm_dir> <run_dir> <genome.fa> <standard.out> <curated_LTR.fa> <artifact_lib> <arm...>
# For each arm X (X.fa = one LTR + one INT record per core cluster):
#   families  rmblastn all-vs-all + family.py (--rule and --gap-aware, classes from cand.tsv) -> X/library.fa
#   masking   compare_libs.sh of X/library.fa (masking vs standard, family recovery vs curated)
#   artifacts RepeatMasker of X.fa against <artifact_lib> (rice 7.0.0), classify2.py with the run's clusters
# Env: THREADS (default 32), PA (RepeatMasker -pa, default 8).
set -eo pipefail
AD=$(realpath $1); RUN=$(realpath $2); G=$3; STD=$4; CUR=$5; ART=$6; shift 6
B=$(cd "$(dirname "$0")/.." && pwd); R=$B/recipe
RM="/usr/local/RepeatMasker/RepeatMasker -e rmblast -no_is -norna -nolow"
for X in "$@"; do (
  mkdir -p $AD/$X; cd $AD/$X
  THREADS=${THREADS:-32} bash $R/rmb_ava.sh ../$X.fa ../$X.fa 5 ava.tsv qseqid sseqid pident length qstart qend sstart send qlen slen bitscore > /dev/null 2>&1
  python3 $R/family.py ../$X.fa ava.tsv library.fa families.tsv --rule and --cand $RUN/cand.tsv --gap-aware > family.log 2>&1
  mkdir -p art; $RM -pa ${PA:-8} -lib $ART -dir art ../$X.fa > art/log 2>&1
  python3 $R/classify2.py ../$X.fa art/$X.fa.out $RUN/elem/clusters_nested.json > classes.txt
  PA=${PA:-8} bash $R/compare_libs.sh $G $STD $CUR $AD/$X/score $X=$AD/$X/library.fa > score.txt 2>&1
) & done
wait
for X in "$@"; do echo "== $X"; tail -2 $AD/$X/family.log; cat $AD/$X/score.txt; cat $AD/$X/classes.txt; done
