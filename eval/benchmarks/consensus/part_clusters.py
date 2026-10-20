#!/usr/bin/env python3
"""Experiment 1b (2026-10-19): cluster LTRs and internals as separate units, from the existing all-vs-all.

usage: part_clusters.py <run_dir> <standard.out> [out_prefix]

Units: one INT unit and one LTR unit per element (the element's two LTRs pooled: a hit from either counts).
Links (rmblastn HSPs >= 80% id, coverage = union of HSPs, single linkage):
  INT  covers 80% of both internals                                   (as step4.py)
  LTR  variants: 'shorter' covers 80% of the shorter LTR (step4.py), 'both' 80% of both LTRs,
       'idcov' identity-weighted coverage (matched bases / length) >= 0.8 of both
Scored against curated part labels (Os####_LTR / Os####_INT models, as diag_parts.py): B-cubed precision and
recall of LTR units vs LTR labels and INT units vs INT labels, for the part clusterings and for the run's
element clusters (clusters_nested.json; each element's parts inherit its cluster). Then the pairing graph:
LTR cluster -- INT cluster edges weighted by elements, and whether curated shared parts show up as one LTR
cluster paired with several INT clusters (or the reverse). Writes <out_prefix>.parts.json (unit -> cluster).
"""
import os, sys, re, csv, json, bisect, collections

run, std = sys.argv[1:3]
pre = sys.argv[3] if len(sys.argv) > 3 else None

elems = []; rec2elem = {}; parts = {}
for r in csv.DictReader(open(os.path.join(run, 'cand.tsv')), delimiter='\t'):
    if r['keep'] != '1': continue
    s = r['seq']; e = f"{s}:{r['elem_start']}..{r['elem_end']}"; elems.append(e)
    ls, le, rs, re_ = int(r['lltr_start']), int(r['lltr_end']), int(r['rltr_start']), int(r['rltr_end'])
    rec2elem[f'{s}:{ls}..{le}_LTR'] = (e, 'L'); rec2elem[f'{s}:{rs}..{re_}_LTR'] = (e, 'L')
    rec2elem[f'{s}:{le+1}..{rs-1}_INT'] = (e, 'I')
    parts[e] = (s, (ls, le), (rs, re_), (le + 1, rs - 1))


def union(iv):
    iv = sorted(iv); tot = 0; cs = ce = None
    for s, e in iv:
        if cs is None: cs, ce = s, e
        elif s <= ce + 1: ce = max(ce, e)
        else: tot += ce - cs + 1; cs, ce = s, e
    if cs is not None: tot += ce - cs + 1
    return tot


# best record pair per (element, element, part): coverage of each side and matched bases
pair = collections.defaultdict(lambda: {'q': [], 's': [], 'ql': 0, 'sl': 0, 'm': 0.0})
for l in open(os.path.join(run, 'nested_ava.tsv')):
    q, s, pid, ln, qs, qe, ss, se, ql, sl, bits = l.split('\t'); q = q.split('#')[0]; s = s.split('#')[0]
    if q == s or float(pid) < 80 or q.endswith('_INT') != s.endswith('_INT'): continue
    p = pair[(q, s)]
    p['q'].append((min(int(qs), int(qe)), max(int(qs), int(qe)))); p['s'].append((min(int(ss), int(se)), max(int(ss), int(se))))
    p['ql'] = int(ql); p['sl'] = int(sl); p['m'] += float(pid) * int(ln) / 100
best = {}  # (a, b, kind) with a < b -> (cov_a, cov_b, idcov_a, idcov_b)
for (q, s), p in pair.items():
    if q not in rec2elem or s not in rec2elem: continue
    (a, k), (b, _) = rec2elem[q], rec2elem[s]
    if a == b: continue
    ca, cb = union(p['q']) / p['ql'], union(p['s']) / p['sl']
    ia, ib = min(1, p['m'] / p['ql']), min(1, p['m'] / p['sl'])
    if a > b: a, b, ca, cb, ia, ib = b, a, cb, ca, ib, ia
    v = (ca, cb, ia, ib); o = best.get((a, b, k))
    if o is None or min(v[:2]) > min(o[:2]): best[(a, b, k)] = v

RULES = {
    'shorter': lambda v: max(v[0], v[1]) >= 0.8,
    'both': lambda v: min(v[0], v[1]) >= 0.8,
    'idcov': lambda v: min(v[2], v[3]) >= 0.8,
}
INT_RULE = lambda v: min(v[0], v[1]) >= 0.8


def cluster(kind, rule):
    parent = {e: e for e in elems}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for (a, b, k), v in best.items():
        if k == kind and rule(v):
            ra, rb = find(a), find(b)
            if ra != rb: parent[ra] = rb
    return {e: find(e) for e in elems}


# curated part labels
iv = collections.defaultdict(list)
for l in open(std):
    c = l.split()
    if len(c) < 15 or not c[0].isdigit() or not c[10].startswith('LTR'): continue
    m = re.match(r'(Os\d+)_(LTR|INT)', c[9])
    if m: iv[c[4]].append((int(c[5]), int(c[6]), m[1], m[2]))
for s in iv: iv[s].sort()
st = {s: [x[0] for x in v] for s, v in iv.items()}
mx = max(b - a for v in iv.values() for a, b, *_ in v)
def lab(seq, s, e, part):
    v = iv.get(seq, []); i = bisect.bisect_left(st.get(seq, []), s - mx); cov = collections.Counter()
    for a, b, f, p in v[i:]:
        if a > e: break
        o = min(b, e) - max(a, s) + 1
        if o > 0 and p == part: cov[f] += o
    if not cov: return None
    f, n = cov.most_common(1)[0]
    return f if n >= 0.5 * (e - s + 1) else None
labL = {e: lab(s, *l, 'LTR') or lab(s, *r, 'LTR') for e, (s, l, r, i) in parts.items()}
labI = {e: lab(s, *i, 'INT') for e, (s, l, r, i) in parts.items()}


def b3(cl, labels):
    el = {e: c for e, c in cl.items() if labels[e]}
    byc = collections.defaultdict(list); byf = collections.defaultdict(list)
    for e, c in el.items(): byc[c].append(e); byf[labels[e]].append(e)
    P = R = 0.0
    for e, c in el.items():
        both = sum(1 for x in byc[c] if labels[x] == labels[e])
        P += both / len(byc[c]); R += both / len(byf[labels[e]])
    n = len(el); multi = [f for f, v in byf.items() if len(v) >= 2]
    split = sum(1 for f in multi if len({el[x] for x in byf[f]}) > 1)
    mixed = sum(1 for c, v in byc.items() if len(v) >= 2 and len({labels[x] for x in v}) > 1)
    return f'n {n:5d} clusters {len(byc):5d} (multi {sum(len(v) > 1 for v in byc.values()):4d})  B3 P {P/n:.3f} R {R/n:.3f}  split {split:3d}/{len(multi)}  mixed {mixed}'


ref = json.load(open(os.path.join(run, 'elem', 'clusters_nested.json')))
elemcl = {e: ref[e] for e in elems}
print(f'{len(elems)} elements; labelled LTR {sum(1 for e in elems if labL[e])}, INT {sum(1 for e in elems if labI[e])}')
print(f"{'element clusters':22s} LTR {b3(elemcl, labL)}")
print(f"{'':22s} INT {b3(elemcl, labI)}")
ci = cluster('I', INT_RULE)
print(f"{'INT units (80/80)':22s} INT {b3(ci, labI)}")
cls = {}
for name, rule in RULES.items():
    cls[name] = cluster('L', rule)
    print(f"{'LTR units (' + name + ')':22s} LTR {b3(cls[name], labL)}")

# pairing graph for the 'both' LTR rule: curated LTR/INT co-occurrence vs cluster co-occurrence
cl = cls['both']
size = lambda c: collections.Counter(c.values())
sL, sI = size(cl), size(ci)
edges = collections.Counter((cl[e], ci[e]) for e in elems)
ml = collections.Counter(a for (a, b) in edges if sL[a] > 1 and sI[b] > 1)
mi = collections.Counter(b for (a, b) in edges if sL[a] > 1 and sI[b] > 1)
print(f"pairing graph (LTR 'both' x INT, multi-member clusters only): {sum(1 for a, b in edges if sL[a] > 1 and sI[b] > 1)} edges; "
      f"LTR clusters with >1 INT partner {sum(1 for v in ml.values() if v > 1)}/{len(ml)}; INT clusters with >1 LTR partner {sum(1 for v in mi.values() if v > 1)}/{len(mi)}")
# curated shared parts: does each curated LTR family land in one LTR cluster, and its INT partners in distinct INT clusters?
def major(c, labels, f):
    k = collections.Counter(c[e] for e in elems if labels[e] == f)
    return k.most_common(1)[0] if k else (None, 0)
for L, Is in (('Os0173', ('Os0173', 'Os2463', 'Os3061', 'Os3063')), ('Os0028', ('Os0028', 'Os0592')),
              ('Os0060', ('Os0060', 'Os2125')), ('Os0032', ('Os0032', 'Os3344'))):
    lc, ln = major(cl, labL, L)
    lp = collections.Counter(labL[e] for e in elems if cl[e] == lc and labL[e])
    row = []
    for I in Is:
        icl, n = major(ci, labI, I)
        ip = collections.Counter(labI[e] for e in elems if ci[e] == icl and labI[e])
        row.append(f"{I}_INT->{n}/{sum(1 for e in elems if labI[e] == I)} in cluster of {sum(ip.values())} ({','.join(f'{k}:{v}' for k, v in ip.most_common(3))})")
    print(f"  {L}_LTR->{ln}/{sum(1 for e in elems if labL[e] == L)} in cluster of {sum(lp.values())} ({','.join(f'{k}:{v}' for k, v in lp.most_common(3))}) | " + '; '.join(row))
if pre:
    json.dump({'LTR': {k: {e: c for e, c in v.items()} for k, v in cls.items()}, 'INT': ci}, open(f'{pre}.parts.json', 'w'))
