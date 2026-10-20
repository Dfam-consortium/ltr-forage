#!/usr/bin/env python3
"""Experiment 1 (2026-10-19): score element clusterings against curated families.

usage: cluster_vs_curated.py <genome.fa.fai> <standard.out> <famkey-regex> <name=kind:path[:path2]> ...

Each element is labelled with the curated LTR family (key = first group of <famkey-regex> on the
RepeatMasker family name, e.g. '(Os\\d+)' for the Jiang rice library) that covers most of it in the
standard annotation; elements less than half covered by curated LTR families are unlabelled and left out.
Clusterings:
  ours:<clusters.json>                       ltr-forage element clusters (element -> cluster)
  ourfam:<clusters.json>:<families.tsv>      ltr-forage families (clusters grouped by family.py)
  hirepeat:<original.bed>:<membership.tsv>   HiRepEat consensus membership (only elements in a family)
Reported per clustering: elements clustered and labelled, B-cubed precision (how pure an element's
cluster is) and recall (how much of its curated family is in its cluster), curated families touched,
and families split over >1 cluster / clusters mixing >1 family (counting clusters with >= 2 labelled members).
"""
import os, sys, re, json, csv, bisect, collections

fai, std, famre = sys.argv[1:4]
seq_order = [l.split('\t')[0] for l in open(fai)]
fre = re.compile(famre)

# curated LTR intervals per sequence
iv = collections.defaultdict(list)
for l in open(std):
    c = l.split()
    if len(c) < 15 or not c[0].isdigit() or not c[10].startswith('LTR'): continue
    m = fre.match(c[9]); key = m[1] if m else c[9]
    iv[c[4]].append((int(c[5]), int(c[6]), key))
for s in iv: iv[s].sort()
starts = {s: [x[0] for x in v] for s, v in iv.items()}
maxlen = max((b - a for v in iv.values() for a, b, _ in v), default=0)


def label(seq, s, e):
    v = iv.get(seq)
    if not v: return None
    i = bisect.bisect_left(starts[seq], s - maxlen); cov = collections.Counter()
    for a, b, k in v[i:]:
        if a > e: break
        o = min(b, e) - max(a, s) + 1
        if o > 0: cov[k] += o
    if not cov: return None
    k, n = cov.most_common(1)[0]
    return k if sum(cov.values()) >= 0.5 * (e - s + 1) else None


def load(kind, paths):
    """element (seq, start, end) -> cluster id"""
    out = {}
    if kind in ('ours', 'ourfam'):
        cl = json.load(open(paths[0]))
        fam = {}
        if kind == 'ourfam':
            for r in csv.DictReader(open(paths[1]), delimiter='\t'):
                m = re.match(r'(\S+):(\d+)\.\.(\d+)_', r['origin'])
                fam[(m[1], int(m[2]), int(m[3]))] = r['family']
            # record -> element cluster -> family
            ce = collections.defaultdict(set)
            for e, c in cl.items():
                m = re.match(r'(\S+):(\d+)\.\.(\d+)$', e); ce[m[1]].add((int(m[2]), int(m[3]), c))
            cfam = {}
            for (s, a, b), f in fam.items():
                for x, y, c in ce[s]:
                    if x <= a and b <= y: cfam[c] = f
        for e, c in cl.items():
            m = re.match(r'(\S+):(\d+)\.\.(\d+)$', e)
            if kind == 'ourfam':
                if c not in cfam: continue
                c = cfam[c]
            out[(m[1], int(m[2]), int(m[3]))] = c
    else:
        loc = {}
        for l in open(paths[0]):
            if l.startswith('#'): continue
            si, s, e, eid, rnd = l.split()[:5]
            loc[f'r{rnd}_{eid}'] = (seq_order[int(si)], int(s) + 1, int(e))
        for r in csv.DictReader(open(paths[1]), delimiter='\t'):
            for m in r['members'].split(','):
                if m in loc: out[loc[m]] = r['consensus']
    return out


def score(name, el):
    lab = {e: label(*e) for e in el}
    el = {e: c for e, c in el.items() if lab[e]}
    byc = collections.defaultdict(list); byf = collections.defaultdict(list)
    for e, c in el.items(): byc[c].append(e); byf[lab[e]].append(e)
    P = R = 0.0
    for e, c in el.items():
        same_c = byc[c]; f = lab[e]
        both = sum(1 for x in same_c if lab[x] == f)
        P += both / len(same_c); R += both / len(byf[f])
    n = len(el)
    multi = [c for c, v in byc.items() if len(v) >= 2]
    mixed = sum(1 for c in multi if len({lab[x] for x in byc[c]}) > 1)
    fam_multi = [f for f, v in byf.items() if len(v) >= 2]
    split = sum(1 for f in fam_multi if len({el[x] for x in byf[f]}) > 1)
    if os.environ.get('MIXED'):
        with open(f"{os.environ['MIXED']}.{name}.mixed.tsv", 'w') as fh:
            for c in multi:
                k = collections.Counter(lab[x] for x in byc[c])
                if len(k) > 1: fh.write(f"{c}\t{len(byc[c])}\t" + ','.join(f'{f}:{n}' for f, n in k.most_common()) + '\n')
    print(f'{name:10s} elements {n:5d}  clusters {len(byc):5d}  B3 precision {P/n:.3f} recall {R/n:.3f}  '
          f'curated families {len(byf):4d}  split {split}/{len(fam_multi)}  mixed {mixed}/{len(multi)}')


def common(sets):
    """elements present (>= 80% reciprocal overlap) in every clustering, keyed by the first one's element"""
    def idx(el):
        d = collections.defaultdict(list)
        for e in el: d[e[0]].append(e)
        for v in d.values(): v.sort()
        return d
    def match(e, d):
        for x in d.get(e[0], ()):
            if x[1] > e[2]: break
            o = min(x[2], e[2]) - max(x[1], e[1]) + 1
            if o >= 0.8 * (e[2] - e[1] + 1) and o >= 0.8 * (x[2] - x[1] + 1): return x
    idxs = [idx(el) for el in sets[1:]]
    keep = {}
    for e in sets[0]:
        m = [match(e, d) for d in idxs]
        if all(m): keep[e] = m
    out = [{e: sets[0][e] for e in keep}]
    for i, el in enumerate(sets[1:]):
        out.append({e: el[m[i]] for e, m in keep.items()})
    return out


specs = []
for spec in sys.argv[4:]:
    name, rest = spec.split('=', 1)
    kind, *paths = rest.split(':')
    specs.append((name, load(kind, paths)))
print('all elements each tool clustered:')
for name, el in specs: score(name, el)
print('elements clustered by every tool (>= 80% reciprocal overlap):')
for (name, _), el in zip(specs, common([el for _, el in specs])): score(name, el)
