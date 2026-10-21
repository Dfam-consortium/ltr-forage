#!/usr/bin/env python3
"""Depth profiles along library records from a non-competitive rmblastn search.

usage: hits_to_profile.py <lib.fa> <hits.tsv> [--min-score 225] > depth.tsv

hits.tsv: rmblastn --outfmt "6 score qseqid qstart qend sseqid sstart send sstrand"
with the genome (or a sample) as query and the library as the database, run
without mask level, so every record collects every genomic copy it aligns to.
Each hit adds 1 along its record span (sstart..send).  Output: record<TAB>d0,d1,...
(the --profile input of excise_by_depth.py)."""
import sys, argparse
ap = argparse.ArgumentParser(); ap.add_argument('lib'); ap.add_argument('hits')
ap.add_argument('--min-score', type=int, default=225); a = ap.parse_args()
L = {}; order = []; n = None
for l in open(a.lib):
    if l.startswith('>'): n = l[1:].split()[0].split('#')[0]; L[n] = 0; order.append(n)
    else: L[n] += len(l.strip())
diff = {k: [0] * (v + 2) for k, v in L.items()}
for line in open(a.hits):
    f = line.split()
    if len(f) < 7 or not f[0].isdigit() or int(f[0]) < a.min_score: continue
    r = f[4].split('#')[0]                       # library names may carry a #class suffix
    if r not in diff: continue
    s, e = sorted((int(f[5]), int(f[6]))); n = L[r]
    s = max(1, min(s, n)); e = max(1, min(e, n))
    diff[r][s] += 1; diff[r][e + 1] -= 1
out = sys.stdout
for k in order:
    d = diff[k]; r = 0; p = []
    for i in range(1, L[k] + 1): r += d[i]; p.append(r)
    out.write(k + '\t' + ','.join(map(str, p)) + '\n')
