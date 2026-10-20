#!/usr/bin/env python3
"""What did the cuts remove?  Each cut base is labelled by the rice 7.0.0 annotation of the
original record: the host's own family (main LTR family of the record), its own family's LTR,
another LTR family, a non-LTR TE, or unannotated.  usage: cut_quality.py <r0.fa.out> <cuts.tsv>..."""
import sys, re, collections
rm, cutfiles = sys.argv[1], sys.argv[2:]
segs = collections.defaultdict(list)
for l in open(rm):
    c = l.split()
    if len(c) < 14 or not c[0].isdigit(): continue
    segs[c[4]].append((int(c[5]), int(c[6]), re.match(r'(Os\d+|[^_]+)', c[9])[1], '_LTR' in c[9], c[10]))
prim = {}
for r, v in segs.items():
    cnt = collections.Counter()
    for a, b, f, isl, cl in v:
        if cl.startswith('LTR') and not isl: cnt[f] += b - a + 1
    prim[r] = cnt.most_common(1)[0][0] if cnt else None
tot = collections.Counter()
for cf in cutfiles[:1]:   # round-1 cuts are in r0 coordinates
    for l in open(cf):
        r, x, y = l.split()[:3]; x, y = int(x), int(y)
        lab = {}
        for a, b, f, isl, cl in segs[r]:
            if not cl.startswith('LTR'): k = 'non-LTR TE'
            elif f == prim[r]: k = "host family's own LTR" if isl else "host family's own internal"
            else: k = 'other LTR family'
            for p in range(max(a, x), min(b, y) + 1): lab.setdefault(p, k)
        for p in range(x, y + 1): tot[lab.get(p, 'unannotated')] += 1
s = sum(tot.values())
for k, v in tot.most_common(): print(f'{k:30s} {v:9,d} bp  {100*v/s:5.1f}%')
