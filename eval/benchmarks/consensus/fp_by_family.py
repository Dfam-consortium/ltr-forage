#!/usr/bin/env python3
"""Where a driver library's false-positive masking comes from (2026-10-19, zebra finch segmental duplications).

usage: fp_by_family.py <run_dir> <test.out> <standard.out> [top]

Per library family: bases masked inside / outside the standard's LTR annotation (TP / FP), what the standard
calls the FP bases (unannotated, or another TE class), the size of the element cluster the model came from
(<run_dir>/elem/clusters_nested.json via families.tsv origins), and the number of test hits.  Then FP totals
by cluster size, which says whether the 1-2 copy models (which the flank test cannot screen) carry it.
"""
import sys, os, re, json, bisect, collections, csv

run, test, std = sys.argv[1:4]
top = int(sys.argv[4]) if len(sys.argv) > 4 else 25


def rows(p):
    for l in open(p):
        c = l.split()
        if len(c) >= 15 and c[0].isdigit(): yield c


# standard: merged LTR intervals, and every annotation for labelling FP bases
ltr = collections.defaultdict(list); ann = collections.defaultdict(list)
for c in rows(std):
    s, a, b = c[4], int(c[5]), int(c[6])
    if c[10].startswith('LTR'): ltr[s].append((a, b))
    else: ann[s].append((a, b, c[10].split('/')[0]))
def merge(v):
    o = []
    for a, b in sorted(v):
        if o and a <= o[-1][1] + 1: o[-1][1] = max(o[-1][1], b)
        else: o.append([a, b])
    return o
ltr = {s: merge(v) for s, v in ltr.items()}; lst = {s: [x[0] for x in v] for s, v in ltr.items()}
for s in ann: ann[s].sort()
ast = {s: [x[0] for x in v] for s, v in ann.items()}
amx = max((b - a for v in ann.values() for a, b, _ in v), default=0)


def overlap(iv, starts, a, b, mx=0):
    i = max(0, bisect.bisect_right(starts, a - mx) - 1)
    for x in iv[i:]:
        if x[0] > b: break
        o = min(b, x[1]) - max(a, x[0]) + 1
        if o > 0: yield x, o


fam2o = {r['new_name']: r['origin'] for r in csv.DictReader(open(os.path.join(run, 'families.tsv')), delimiter='\t')}
clus = json.load(open(os.path.join(run, 'elem', 'clusters_nested.json')))
size = collections.Counter(clus.values())
by_seq = collections.defaultdict(list)
for e in clus:
    m = re.match(r'(\S+):(\d+)\.\.(\d+)$', e); by_seq[m[1]].append((int(m[2]), int(m[3]), e))
def csize(rec):
    m = re.match(r'(\S+):(\d+)\.\.(\d+)_', rec)
    for s, t, e in by_seq.get(m[1], ()):
        if s <= int(m[2]) and int(m[3]) <= t: return size[clus[e]]
    return 0

tp = collections.Counter(); fp = collections.Counter(); hits = collections.Counter()
fpcls = collections.defaultdict(collections.Counter)
for c in rows(test):
    f = c[9]; s, a, b = c[4], int(c[5]), int(c[6]); hits[f] += 1
    x = sum(o for _, o in overlap(ltr.get(s, []), lst.get(s, []), a, b))
    tp[f] += x; fp[f] += b - a + 1 - x
    if b - a + 1 - x > 0:   # label the FP part roughly: bases under other annotations (may double count overlaps)
        other = collections.Counter()
        for (p, q, k), o in overlap(ann.get(s, []), ast.get(s, []), a, b, amx): other[k] += o
        rest = b - a + 1 - x
        for k, o in other.items():
            o = min(o, rest); fpcls[f][k] += o; rest -= o
        fpcls[f]['unannotated'] += max(0, rest)

fams = sorted(fp, key=lambda f: -fp[f])
tot = sum(fp.values())
print(f'FP {tot/1e6:.2f} Mb, TP {sum(tp.values())/1e6:.2f} Mb over {len(fp)} models')
print(f"{'model':22s} {'clu':>4s} {'hits':>6s} {'TP Mb':>7s} {'FP Mb':>7s}  FP called by the standard")
for f in fams[:top]:
    cs = csize(fam2o.get(f, ''))
    print(f'{f:22s} {cs:4d} {hits[f]:6d} {tp[f]/1e6:7.2f} {fp[f]/1e6:7.2f}  ' +
          ', '.join(f'{k} {v/1e6:.2f}' for k, v in fpcls[f].most_common(3)))
bucket = collections.defaultdict(lambda: [0, 0, 0])
for f in fp:
    cs = csize(fam2o.get(f, '')); k = '1' if cs == 1 else '2' if cs == 2 else '3-9' if cs < 10 else '10+'
    bucket[k][0] += 1; bucket[k][1] += tp[f]; bucket[k][2] += fp[f]
print('by element-cluster size of the model:')
for k in ('1', '2', '3-9', '10+'):
    if k in bucket:
        n, t, x = bucket[k]; print(f'  {k:4s} models {n:4d}  TP {t/1e6:6.2f} Mb  FP {x/1e6:6.2f} Mb  precision {t/max(1,t+x):.3f}')
allc = collections.Counter()
for v in fpcls.values(): allc.update(v)
print('FP bases by what the standard calls them: ' + ', '.join(f'{k} {v/1e6:.2f} Mb' for k, v in allc.most_common()))
