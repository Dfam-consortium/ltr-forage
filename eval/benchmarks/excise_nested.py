#!/usr/bin/env python3
"""Prototype nested-insertion cleanup of a library from an all-versus-all
alignment of its records.

Usage: excise_nested.py <lib.fa> <clusters.json> <allvsall.tsv> <out.fa> [--min-clusters 3] [--max-frac 0.7] [--min-len 100]

A position of a record is "shared" when records from at least --min-clusters
other cd-hit clusters align to it.  Shared stretches of --min-len or more
are cut out of the record, unless they cover more than --max-frac of it, in
which case the record is a member of a family the clustering split and is
left alone.  Records shorter than --min-len after cutting are dropped.
"""
import sys, json, collections, argparse
ap = argparse.ArgumentParser()
ap.add_argument('lib'); ap.add_argument('clusters'); ap.add_argument('hits'); ap.add_argument('out')
ap.add_argument('--min-clusters', type=int, default=3)
ap.add_argument('--max-frac', type=float, default=0.7)
ap.add_argument('--min-len', type=int, default=100)
a = ap.parse_args()
clu = json.load(open(a.clusters))
seqs = {}; order = []; name = None
for l in open(a.lib):
    if l.startswith('>'): name = l[1:].split()[0]; seqs[name] = []; order.append(name)
    else: seqs[name].append(l.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}
# per query position: set of other clusters aligning
cov = {k: [set() for _ in range(len(v) + 1)] for k, v in seqs.items()}
for l in open(a.hits):
    q, s, pid, ln, qs, qe, ss, se, ev = l.split('\t')
    if q == s or clu.get(q) == clu.get(s): continue
    cs = clu.get(s)
    for i in range(int(qs), int(qe) + 1): cov[q][i].add(cs)
out = open(a.out, 'w'); cut_bases = 0; dropped = 0; touched = 0
for k in order:
    s = seqs[k]; L = len(s)
    shared = [len(cov[k][i]) >= a.min_clusters for i in range(1, L + 1)]
    # runs of shared positions
    runs = []; i = 0
    while i < L:
        if shared[i]:
            j = i
            while j < L and shared[j]: j += 1
            if j - i >= a.min_len: runs.append((i, j))
            i = j
        else: i += 1
    total = sum(j - i for i, j in runs)
    if runs and total <= a.max_frac * L:
        keep = []; pos = 0
        for i, j in runs: keep.append(s[pos:i]); pos = j
        keep.append(s[pos:]); new = ''.join(keep); cut_bases += total; touched += 1
    else:
        new = s
    if len(new) < a.min_len: dropped += 1; continue
    out.write(f'>{k}\n{new}\n')
out.close()
print(f'{touched} records trimmed, {cut_bases:,} bases removed, {dropped} records dropped', file=sys.stderr)
