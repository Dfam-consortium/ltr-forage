#!/usr/bin/env python3
"""Experiment 1 follow-up (2026-10-19): why do element clusters mix curated families?

usage: diag_links.py <run_dir> <standard.out> <famkey-regex> [out_prefix]

Rebuilds step4.py's phase-1 links (internals 80% both, LTRs 80% of the shorter, HSPs >= 80% id) and
phase-2 singleton attachments from <run_dir>/nested_ava.tsv and cand.tsv, labels elements with the
curated family covering most of them (as cluster_vs_curated.py), and for every cluster mixing families
says how each minority-family element got in:
  direct    phase-1 link straight to an element of another curated family (link rule too loose)
  chained   no such link; joined through unlabelled elements (single-linkage chaining)
  attached  phase-2 attachment of a singleton
For links it also reports LTR/INT coverage and aligned identity, within-family vs cross-family, to see
whether a stricter LTR rule separates them. Writes <out_prefix>.links.tsv (all labelled links).
"""
import os, sys, re, csv, json, bisect, collections

run, std, famre = sys.argv[1:4]
pre = sys.argv[4] if len(sys.argv) > 4 else None
fre = re.compile(famre)

# elements and their records (as step4.py)
elems = []; rec2elem = {}; ltrlen = {}
for r in csv.DictReader(open(os.path.join(run, 'cand.tsv')), delimiter='\t'):
    if r['keep'] != '1': continue
    s = r['seq']; e = (s, int(r['elem_start']), int(r['elem_end'])); elems.append(e)
    ls, le, rs, re_ = int(r['lltr_start']), int(r['lltr_end']), int(r['rltr_start']), int(r['rltr_end'])
    for k, n in (('L', f'{s}:{ls}..{le}_LTR'), ('I', f'{s}:{le+1}..{rs-1}_INT'), ('R', f'{s}:{rs}..{re_}_LTR')):
        rec2elem[n] = (e, k)


def union(iv):
    iv = sorted(iv); tot = 0; cs = ce = None
    for s, e in iv:
        if cs is None: cs, ce = s, e
        elif s <= ce + 1: ce = max(ce, e)
        else: tot += ce - cs + 1; cs, ce = s, e
    if cs is not None: tot += ce - cs + 1
    return tot


# per record pair: union coverage of each side, and identity over the aligned columns (pid-weighted)
cov = collections.defaultdict(dict)
pair = collections.defaultdict(lambda: {'q': [], 's': [], 'ql': 0, 'sl': 0, 'bits': 0.0, 'm': 0.0, 'n': 0})
for l in open(os.path.join(run, 'nested_ava.tsv')):
    q, s, pid, ln, qs, qe, ss, se, ql, sl, bits = l.split('\t'); q = q.split('#')[0]; s = s.split('#')[0]
    if q == s or float(pid) < 80 or q.endswith('_INT') != s.endswith('_INT'): continue
    p = pair[(q, s)]
    p['q'].append((min(int(qs), int(qe)), max(int(qs), int(qe)))); p['s'].append((min(int(ss), int(se)), max(int(ss), int(se))))
    p['ql'] = int(ql); p['sl'] = int(sl); p['bits'] += float(bits); p['m'] += float(pid) * int(ln) / 100; p['n'] += int(ln)
for (q, s), p in pair.items():
    if q not in rec2elem or s not in rec2elem: continue
    (eq, kq), (es, _) = rec2elem[q], rec2elem[s]
    if eq == es: continue
    kind = 'I' if kq == 'I' else 'L'
    cq = union(p['q']) / p['ql']; cs = union(p['s']) / p['sl']; pid = 100 * p['m'] / p['n']
    # identity-weighted coverage of the shorter record: matched columns / shorter length
    idc = min(1.0, p['m'] / min(p['ql'], p['sl']))
    for a, b, ca, cb in ((eq, es, cq, cs), (es, eq, cs, cq)):
        d = cov[(a, b)]; old = d.get(kind)
        if old is None or ca > old[0] or (ca == old[0] and p['bits'] > old[2]): d[kind] = (ca, cb, p['bits'], pid, idc)

parent = {e: e for e in elems}
def find(x):
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
links = []
for (a, b), d in cov.items():
    if a < b and 'I' in d and 'L' in d:
        ci, cj = d['I'][:2]; li, lj = d['L'][:2]
        if ci >= 0.8 and cj >= 0.8 and max(li, lj) >= 0.8:
            links.append((a, b, d)); ra, rb = find(a), find(b)
            if ra != rb: parent[ra] = rb
comp = collections.defaultdict(list)
for e in elems: comp[find(e)].append(e)
attach = {}
cand = collections.defaultdict(lambda: collections.defaultdict(float))
for (a, b), d in cov.items():
    if len(comp[find(a)]) != 1 or len(comp[find(b)]) == 1: continue
    if 'I' in d and 'L' in d and d['I'][0] >= 0.8 and d['L'][0] >= 0.8:
        cand[a][find(b)] = max(cand[a][find(b)], d['I'][2] + d['L'][2])
for a, c in cand.items(): attach[a] = max(c.items(), key=lambda x: x[1])[0]
clus = {e: attach.get(e, find(e)) for e in elems}

# sanity: same clustering as the run's clusters_nested.json
ref = json.load(open(os.path.join(run, 'elem', 'clusters_nested.json')))
key = lambda e: f'{e[0]}:{e[1]}..{e[2]}'
grp = collections.defaultdict(set)
for e in elems: grp[clus[e]].add(key(e))
rgrp = collections.defaultdict(set)
for k, c in ref.items(): rgrp[c].add(k)
same = len({frozenset(v) for v in grp.values()} & {frozenset(v) for v in rgrp.values()})
print(f'rebuilt {len(grp)} clusters, run has {len(rgrp)}, identical {same}')

# curated labels (as cluster_vs_curated.py)
iv = collections.defaultdict(list)
for l in open(std):
    c = l.split()
    if len(c) < 15 or not c[0].isdigit() or not c[10].startswith('LTR'): continue
    m = fre.match(c[9]); iv[c[4]].append((int(c[5]), int(c[6]), m[1] if m else c[9]))
for s in iv: iv[s].sort()
starts = {s: [x[0] for x in v] for s, v in iv.items()}
maxlen = max(b - a for v in iv.values() for a, b, _ in v)
def label(seq, s, e):
    v = iv.get(seq)
    if not v: return None
    i = bisect.bisect_left(starts[seq], s - maxlen); cov_ = collections.Counter()
    for a, b, k in v[i:]:
        if a > e: break
        o = min(b, e) - max(a, s) + 1
        if o > 0: cov_[k] += o
    if not cov_: return None
    k, _ = cov_.most_common(1)[0]
    return k if sum(cov_.values()) >= 0.5 * (e - s + 1) else None
lab = {e: label(*e) for e in elems}

# how each minority-family element joined its cluster
adj = collections.defaultdict(set)
for a, b, _ in links: adj[a].add(b); adj[b].add(a)
how = collections.Counter(); mixed = 0; rows = []
for c, mem in collections.defaultdict(list, {k: [e for e in elems if clus[e] == k] for k in set(clus.values())}).items():
    fams = collections.Counter(lab[e] for e in mem if lab[e])
    if len(fams) < 2: continue
    mixed += 1; major = fams.most_common(1)[0][0]
    for e in mem:
        f = lab[e]
        if not f or f == major: continue
        if e in attach: h = 'attached'
        elif any(lab[x] and lab[x] != f for x in adj[e]): h = 'direct'
        else: h = 'chained'
        how[h] += 1
    rows.append((len(mem), key(c), ','.join(f'{f}:{n}' for f, n in fams.most_common())))
print(f'{mixed} mixed clusters; minority-family elements by how they joined: ' + ', '.join(f'{k} {v}' for k, v in how.most_common()))
for r in sorted(rows, reverse=True)[:20]: print(f'  {r[0]:4d} {r[1]:28s} {r[2]}')

# link statistics: within- vs cross-family
def q(v, p): v = sorted(v); return v[min(len(v) - 1, int(p * len(v)))] if v else float('nan')
st = {'within': collections.defaultdict(list), 'cross': collections.defaultdict(list)}
out = open(f'{pre}.links.tsv', 'w') if pre else None
if out: out.write('a\tb\tfam_a\tfam_b\tint_cov_a\tint_cov_b\tint_pid\tltr_cov_a\tltr_cov_b\tltr_pid\tltr_idcov_shorter\n')
for a, b, d in links:
    if not lab[a] or not lab[b]: continue
    t = 'within' if lab[a] == lab[b] else 'cross'
    I, L = d['I'], d['L']
    for k, v in (('LTR cov longer', min(L[0], L[1])), ('LTR pid', L[3]), ('LTR idcov', L[4]), ('INT cov min', min(I[0], I[1])), ('INT pid', I[3])):
        st[t][k].append(v)
    if out: out.write(f'{key(a)}\t{key(b)}\t{lab[a]}\t{lab[b]}\t{I[0]:.3f}\t{I[1]:.3f}\t{I[3]:.1f}\t{L[0]:.3f}\t{L[1]:.3f}\t{L[3]:.1f}\t{L[4]:.3f}\n')
print(f"links between labelled elements: within-family {len(st['within']['LTR pid'])}, cross-family {len(st['cross']['LTR pid'])}")
print(f"{'':16s} {'within: p5 p25 p50':>26s}   {'cross: p50 p75 p95':>26s}")
for k in ('LTR cov longer', 'LTR pid', 'LTR idcov', 'INT cov min', 'INT pid'):
    w, x = st['within'][k], st['cross'][k]
    print(f'{k:16s} {q(w,.05):8.3f} {q(w,.25):8.3f} {q(w,.5):8.3f}   {q(x,.5):8.3f} {q(x,.75):8.3f} {q(x,.95):8.3f}')
# what a stricter LTR rule would cut: identity-weighted coverage of the shorter LTR >= t
for t in (0.70, 0.75, 0.80, 0.85):
    kw = sum(v >= t for v in st['within']['LTR idcov']); kx = sum(v >= t for v in st['cross']['LTR idcov'])
    print(f'LTR idcov >= {t:.2f}: keeps within {kw}/{len(st["within"]["LTR idcov"])}, cross {kx}/{len(st["cross"]["LTR idcov"])}')
