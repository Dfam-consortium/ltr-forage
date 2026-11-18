#!/usr/bin/env python3
"""Containment test: drop a single-member cluster's exemplar when the rest of the
library already covers it.

usage: containment.py <clusters.json> <exemplar.fa> <self_hits.tsv> <out_prefix> [threshold ...]
self_hits.tsv: blastn -outfmt '6 qseqid sseqid pident length qstart qend qlen' of the
exemplar singletons' records against the whole exemplar library.
Coverage is counted over the element (all its records together), from HSPs at
>= 80% identity to records of other elements that are still kept.  Singletons are
visited shortest first, so a dropped one never served as another's only cover.
"""
import sys, os, json, collections
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from records import element_lookup
clusfile, fa, hits, prefix = sys.argv[1:5]
ths = [float(x) for x in sys.argv[5:]] or [0.9]
clus = json.load(open(clusfile))
size = collections.Counter(clus.values())
recs, order = {}, []
n = None
for l in open(fa):
    if l.startswith('>'):
        n = l[1:].split()[0]; recs[n] = []; order.append(n)
    else: recs[n].append(l.strip())
recs = {k: ''.join(v) for k, v in recs.items()}
# record -> cluster, through the record's own element (env TSV: the ltr-forage table)
elem = element_lookup(clus)
rec_elem = {r: elem(r) for r in order}
lost = [r for r in order if rec_elem[r] is None]
if lost: sys.exit(f'containment.py: {len(lost)} record(s) belong to no clustered element, e.g. {lost[0]}')
rec_clus = {r: clus[e] for r, e in rec_elem.items()}
clus_recs = collections.defaultdict(list)
for r in order: clus_recs[rec_clus[r]].append(r)
single = [c for c in clus_recs if size[c] == 1]
iv = collections.defaultdict(lambda: collections.defaultdict(list))   # q record -> subject cluster -> intervals
for l in open(hits):
    q, s, pid, ln, qs, qe, ql = l.split()[:7]
    if float(pid) < 80 or rec_clus.get(q) == rec_clus.get(s): continue
    iv[q][rec_clus[s]].append((min(int(qs), int(qe)), max(int(qs), int(qe))))
def covered(c, kept):
    tot = cov = 0
    for r in clus_recs[c]:
        L = len(recs[r]); tot += L
        v = sorted(x for sc, xs in iv[r].items() if sc in kept for x in xs)
        last = 0
        for a, b in v:
            if b > last: cov += b - max(a, last + 1) + 1; last = b
    return cov / tot
for th in ths:
    kept = set(clus_recs)
    for c in sorted(single, key=lambda c: sum(len(recs[r]) for r in clus_recs[c])):
        if covered(c, kept) >= th: kept.discard(c)   # iv holds no hit within a cluster, so c itself never counts
    ns = sum(1 for c in kept if size[c] == 1)
    out = f'{prefix}_{round(th * 100)}.fa'
    with open(out, 'w') as f:
        for r in order:
            if rec_clus[r] in kept: f.write(f'>{r}\n{recs[r]}\n')
    print(f'threshold {th:.2f}: {len(kept)} clusters kept ({len(kept)-ns} multi-element, {ns} singletons; '
          f'{len(single)-ns} of {len(single)} singletons dropped) -> {out}')
