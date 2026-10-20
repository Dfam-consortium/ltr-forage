#!/bin/bash
# ltr-forage library recipe, steps 1-14 (see docs/benchmarks.md, 2026-10-15).
# usage: library_recipe.sh <genome.fa> <genome.sufr> <depth_ref.fa> <depth_scale> <min_depth> <abs_depth> <outdir>
#   depth_ref.fa: genome (or a sample of it) RepeatMasked in each depth-cleanup round;
#   depth_scale: 0.1 for a 10% sample; min_depth/abs_depth: excise_by_depth floors before scaling.
# Env: THREADS (default 32), HMM (TEfam.hmm), LTRFORAGE (binary),
#      TIER (pbs | strict, default pbs), GAPAWARE (1 for gappy assemblies such as lotus),
#      START (7: reuse cand.tsv and clean/clean.fa; 12: also final.fa and elem/clusters_nested.json),
#      ALIGNER (rmb: the rmblast port via rmb_ava.sh, default; blast: BLAST+ dc-megablast).
set -eo pipefail
G=$1; SUFR=$2; DREF=$3; SCALE=$4; MIND=$5; ABSD=$6; OUT=$7
HERE=$(cd "$(dirname "$0")" && pwd); B=$(dirname "$HERE")
T=${THREADS:-32}; HMM=${HMM:-/usr/local/LTR_retriever-2.9.0/database/TEfam.hmm}
LF=${LTRFORAGE:-$B/../../target/release/ltr-forage}   # this checkout's release build
AVA="qseqid sseqid pident length qstart qend sstart send qlen slen bitscore"
mkdir -p $OUT/elem; cd $OUT
# search <query.fa> <db.fa> <max_hsps> <out.tsv> <fields...>: BLAST+ tabular columns from either aligner
search() {
  local q=$1 db=$2 mh=$3 o=$4; shift 4
  if [ "${ALIGNER:-rmb}" = rmb ]; then THREADS=$T bash $HERE/rmb_ava.sh $q $db $mh $o "$@"
  else makeblastdb -in $db -dbtype nucl -out $db.db > /dev/null
       blastn -task dc-megablast -query $q -db $db.db -evalue 1e-10 -num_threads $T -max_target_seqs 20000 -max_hsps $mh -outfmt "6 $*" > $o 2> /dev/null; fi; }
ava() { search $1 $1 5 $2 $AVA; }
step() { echo "[$(date +%T)] $*" | tee -a steps.log; }
export TSV=$OUT/cand.tsv JIANG=/dev/null
if [ "${START:-1}" -le 6 ]; then
# 1-5: detect, re-extend, annotate/decide, library
step "1-5 ltr-forage"
$LF $SUFR --reextend --reextend-ltr-max 7000 --format tsv --domain-hmm $HMM --library lib.raw.fa -j $T > cand.tsv 2> ltr-forage.err
awk '/^>/{print $1} !/^>/' lib.raw.fa > lib.fa
# 6: depth cleanup, three rounds (rmblastn of the depth sample against the library, RepeatMasker -lib settings)
step "6 depth cleanup"
BATCHES=16 THREADS=$(( T / 16 > 0 ? T / 16 : 1 )) bash $B/clean_library_rmb.sh lib.fa $DREF clean 3 $SCALE $(bash $HERE/pick_matrix.sh $DREF) --ratio 3 --min-depth $MIND --min-len 50 --abs-depth $ABSD --drop-frac 0.7 > clean.log 2>&1
fi
if [ "${START:-1}" -le 11 ]; then
# 7: nested LTR excision with voting (needs a first clustering of the cleaned records)
step "7 nested vote"
ava clean/clean.fa clean_ava.tsv
python3 $HERE/step4.py clean_ava.tsv clean_ava.tsv elem/clusters_clean.json > elem/clean.txt 2>&1
python3 $HERE/nested_vote.py clean_ava.tsv elem/clusters_clean.json clean/clean.fa clean_nested.fa > nested_cut.txt
# 8-9: element clustering and exemplars
step "8-9 cluster + exemplar"
ava clean_nested.fa nested_ava.tsv
python3 $HERE/step4.py nested_ava.tsv nested_ava.tsv elem/clusters_nested.json > elem/nested.txt 2>&1
python3 $HERE/exemplar_lib.py elem/clusters_nested.json clean_nested.fa exemplar.fa > exemplar.txt
awk '/^>/{print $1} !/^>/' exemplar.fa > exemplar.ids.fa
# 10: containment at 0.9
step "10 containment"
search exemplar.ids.fa exemplar.ids.fa 20 ex_self.tsv qseqid sseqid pident length qstart qend qlen
python3 $HERE/containment.py elem/clusters_nested.json exemplar.ids.fa ex_self.tsv contained 0.9 | tee containment.txt
# 11: nested LTR excision from internals (idea 2), drop heavily nested singletons, two rounds
step "11 exemplar nested excision"
cp contained_90.fa r0.fa
for i in 1 2; do
  p=r$((i-1))
  search $p.fa $p.fa 20 $p.hits qseqid sseqid pident length qstart qend sstart send qlen slen
  python3 $HERE/excise_contained.py $p.fa elem/clusters_nested.json $p.hits r$i.fa cuts$i.tsv --int-cov 2 --drop-multi 2 2>&1 | tee -a excise.txt
done
cp r2.fa final.fa
fi
# 12: primer binding site per element, from ltr-forage's own tRNA search (pbs_* columns;
# matches eval/benchmarks/pbs/pbs_annotate.py + ARAGORN on 99.5% of rice elements)
step "12 PBS"
python3 - cand.tsv > pbs.tsv <<'PY'
import csv, sys
print('element\tpbs_mm12\tpbs_trna\tpbs_offset')
for r in csv.DictReader(open(sys.argv[1]), delimiter='\t'):
    print(f"{r['seq']}:{r['elem_start']}..{r['elem_end']}", r['pbs_mm12'], r['pbs_trna'], r['pbs_offset'], sep='\t')
PY
# 13: tier.  pbs (default): a singleton is kept on LTR identity >= 93 OR a PBS;
# strict: identity >= 90 AND a PBS (fewer families, less masking; see tier.py)
TIER=${TIER:-pbs}
step "13 tier ($TIER)"
if [ "$TIER" = strict ]; then TOPT="--min-id 90 --pbs pbs.tsv --strict"; else TOPT="--min-id 93 --pbs pbs.tsv"; fi
python3 $HERE/tier.py final.fa elem/clusters_nested.json cand.tsv core.fa candidates.fa $TOPT | tee tier.txt
# 14: 80-80 families, LTRs and internals must both match (--gap-aware when GAPAWARE=1)
step "14 families"
ava core.fa core_ava.tsv
python3 $HERE/family.py core.fa core_ava.tsv library.fa families.tsv --rule and --cand cand.tsv ${GAPAWARE:+--gap-aware} 2>&1 | tee family.txt
step "done: $(grep -c '>' library.fa) records, $(cut -f1 families.tsv | tail -n +2 | sort -u | wc -l) families"
