#!/usr/bin/env python3
"""Characterize singleton elements of an exemplar library from a RepeatMasker run of the
whole library against the full rice 7.0.0 library.
usage: classify2.py <lib.fa> <lib.fa.out> <clusters.json>"""
import sys, os, re, json, collections
def famkey(n):
    if os.environ.get('FAMKEY') == 'strip': return re.sub(r'([-_]?LTR|[-_]?I|_INT)$', '', n)
    return re.match(r'(Os\d+|[^_]+)', n)[1]
fa, rmout, clusfile = sys.argv[1:4]
L = {}; n = None
for l in open(fa):
    if l.startswith('>'): n = l[1:].split()[0]; L[n] = 0
    else: L[n] += len(l.strip())
clus = json.load(open(clusfile)); size = collections.Counter(clus.values())
by = collections.defaultdict(list)
for e, c in clus.items():
    m = re.match(r'(\S+):(\d+)\.\.(\d+)', e); by[m[1]].append((int(m[2]), int(m[3]), e))
def elem(r):
    m = re.match(r'(\S+):(\d+)\.\.(\d+)_', r)
    for a, b, e in by[m[1]]:
        if a <= int(m[2]) and int(m[3]) <= b: return e
recs = collections.defaultdict(list)
for r in L: recs[elem(r)].append(r)
hits = collections.defaultdict(list)
for l in open(rmout):
    c = l.split()
    if len(c) < 14 or not c[0].isdigit(): continue
    hits[c[4]].append((int(c[5]), int(c[6]), famkey(c[9]), c[10], float(c[1])))
def segs(r):
    v = sorted(hits[r], key=lambda h: h[1]-h[0], reverse=True); taken = []; out = []
    for h in v:
        if sum(max(0, min(h[1], b) - max(h[0], a) + 1) for a, b in taken) < 0.5*(h[1]-h[0]+1):
            taken.append((h[0], h[1])); out.append(h)
    return sorted(out)
info = {}
for e, rs in recs.items():
    ints = [r for r in rs if '_INT' in r]; ltrs = [r for r in rs if '_LTR' in r]
    S = {r: segs(r) for r in rs}
    def famb(rr):
        c = collections.Counter()
        for r in rr:
            for a, b, f, cl, d in S[r]:
                if cl.startswith('LTR'): c[f] += b - a + 1
        return c
    fi, fl = famb(ints), famb(ltrs)
    prim = (fi or fl).most_common(1)[0][0] if (fi or fl) else None
    tot = sum(L[r] for r in rs); ann = sum(b-a+1 for ss in S.values() for a, b, *_ in ss)
    nonltr = sum(b-a+1 for ss in S.values() for a, b, f, cl, d in ss if not cl.startswith('LTR'))
    pat = []
    for r in ints:
        for a, b, f, cl, d in S[r]:
            if b - a + 1 >= 300:
                lab = f if cl.startswith('LTR') else 'nonLTR'
                if not pat or pat[-1] != lab: pat.append(lab)
    ltrpat = [x for x in pat if x != 'nonLTR']
    ins = any(ltrpat[i] != ltrpat[0] and ltrpat[0] in ltrpat[i+1:] for i in range(1, len(ltrpat)))
    fams_int = set(ltrpat)
    divs = sorted(d for ss in S.values() for a, b, f, cl, d in ss if f == prim)
    mdiv = divs[len(divs)//2] if divs else None
    ltr_main = {fl_.most_common(1)[0][0] for r in ltrs for fl_ in [famb([r])] if fl_}
    unann = 1 - ann / tot
    if unann >= 0.5: c = 'A mostly unannotated by rice 7.0.0'
    elif nonltr >= 200: c = 'B carries non-LTR TE sequence (>=200 bp)'
    elif ins: c = 'C another LTR family inserted in the internal (A-B-A)'
    elif len(fams_int) >= 2: c = 'D internal is two LTR families end to end (A-B)'
    elif ltr_main and prim and prim not in ltr_main: c = 'E LTRs are another family than the internal'
    elif mdiv is not None and mdiv > 15: c = 'F one family, >15% from its consensus'
    else: c = 'G one family, close to its consensus'
    info[e] = dict(prim=prim, cat=c, single=size[clus[e]] == 1, pat=pat, unann=unann, mdiv=mdiv, len=tot)
famcount = collections.Counter(v['prim'] for v in info.values())
tab = collections.defaultdict(lambda: [0, 0])
for e, v in info.items():
    if v['single']: tab[v['cat']][famcount[v['prim']] == 1 and v['prim'] is not None] += 1
print('singleton class                                            | family also elsewhere in library | sole representative')
for c in sorted(tab): print(f'{c:58s} | {tab[c][0]:5d} | {tab[c][1]:5d}')
m = [v for v in info.values() if not v['single']]
print('multi-element clusters:', len(m), collections.Counter(v['cat'][0] for v in m))
json.dump(info, open('singleton_classes.json', 'w'), indent=0)
