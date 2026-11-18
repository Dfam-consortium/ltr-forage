#!/usr/bin/env python3
"""Drop internal models with far more genome copies than their cluster has elements (ltr_library.py step 22b).
usage: int_copies.py <models.fa> <depth.tsv> <clusters.json> <scale> <max> <dropped.tsv>

An LTR family's internal region can't outnumber its intact elements by much: internals are not solo, while
LTRs are (solo LTRs, two per element).  An internal model whose genome copies (median sample depth / scale)
exceed <max> per element of its cluster is mostly some other repeat (mm10: an exemplar internal of L1 sequence,
5,999 copies per element, masked 133 Mb outside the curated LTRs; no other model in six genomes went above 707).
models.fa is rewritten without them; the LTR models are kept."""
import sys, os, json, collections, statistics
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from records import element_lookup
mf, depf, clf, scale, mx, outf = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5]), sys.argv[6]

dep = {}
for l in open(depf):
    n, v = l.rstrip('\n').split('\t')[:2]
    x = [int(t) for t in v.split(',') if t]
    dep[n] = statistics.median(x) if x else 0
cl = json.load(open(clf)); size = collections.Counter(cl.values())
element = element_lookup(cl)   # the element a model's origin record is of (env TSV: the ltr-forage table)

recs, drop = [], []
for l in open(mf):
    if l.startswith('>'): recs.append([l, ''])
    else: recs[-1][1] += l
for h, _ in recs:
    n = h[1:].split()[0].split('#')[0]
    e = element(n) if n.endswith('_INT') else None
    if e is None: continue
    per = dep.get(n, 0) / scale / size[cl[e]]
    if per > mx: drop.append((n, size[cl[e]], round(dep.get(n, 0) / scale), round(per)))
gone = {d[0] for d in drop}
with open(mf, 'w') as fh:
    for h, s in recs:
        if h[1:].split()[0].split('#')[0] not in gone: fh.write(h + s)
with open(outf, 'w') as fh:
    fh.write('model\telements\tgenome_copies\tcopies_per_element\n')
    for d in drop: fh.write('\t'.join(map(str, d)) + '\n')
print(f'{len(drop)} internal model(s) with more than {mx:g} genome copies per element dropped')
