#!/usr/bin/env python3
"""2026-10-19: do a library's own LTR (and INT) models satisfy the 80/80/80 rule among themselves?
usage: lib_808080.py <ava.tsv> [pairs.tsv]
ava.tsv: rmblastn all-vs-all of the library (qseqid sseqid pident length qstart qend sstart send qlen slen bitscore).
Models are named <family>_LTR / <family>_INT, optionally with a _<name> suffix, (RepeatMasker '#class' suffix ignored). Two models of the same part
pass 80/80/80 when HSPs at >= 80% identity cover >= 80% of both ('both') or of the shorter ('shorter') and >= 80 bp.
Reports passing pairs, single-linkage groups, and the LTR-group x INT-group association graph (a family links its
LTR model's group to its INT model's group)."""
import sys, re, collections
ava = sys.argv[1]
def union(iv):
    iv = sorted(iv); tot = 0; cs = ce = None
    for s, e in iv:
        if cs is None: cs, ce = s, e
        elif s <= ce + 1: ce = max(ce, e)
        else: tot += ce - cs + 1; cs, ce = s, e
    return tot + (ce - cs + 1 if cs is not None else 0)
pair = collections.defaultdict(lambda: [[], [], 0, 0]); models = {}
for l in open(ava):
    q, s, pid, ln, qs, qe, ss, se, ql, sl, _ = l.split('\t'); q = q.split('#')[0]; s = s.split('#')[0]
    for n, L in ((q, ql), (s, sl)):
        m = re.match(r'(.+?)_(LTR|INT)(_[^_]+)?$', n)
        if m: models[n] = (m[1], m[2], int(L))
    if q == s or float(pid) < 80 or q not in models or s not in models or models[q][1] != models[s][1]: continue
    p = pair[tuple(sorted((q, s)))]; a, b = (q, s) if q < s else (s, q)
    qa = (min(int(qs), int(qe)), max(int(qs), int(qe))); sa = (min(int(ss), int(se)), max(int(ss), int(se)))
    if q == a: p[0].append(qa); p[1].append(sa)
    else: p[0].append(sa); p[1].append(qa)
cov = {k: (union(v[0]) / models[k[0]][2], union(v[1]) / models[k[1]][2], union(v[0]), union(v[1])) for k, v in pair.items()}
out = open(sys.argv[2], 'w') if len(sys.argv) > 2 else None
for rule in ('both', 'shorter'):
    ok = lambda c: (min(c[0], c[1]) if rule == 'both' else max(c[0], c[1])) >= 0.8 and min(c[2], c[3]) >= 80
    print(f'== rule: 80% identity over 80% of {rule}')
    grp = {}
    for part in ('LTR', 'INT'):
        ms = [m for m, v in models.items() if v[1] == part]
        parent = {m: m for m in ms}
        def find(x):
            while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
            return x
        npass = 0
        for (a, b), c in cov.items():
            if models[a][1] == part and ok(c):
                npass += 1; parent[find(a)] = find(b)
                if out and rule == 'both': out.write(f'{a}\t{b}\t{c[0]:.3f}\t{c[1]:.3f}\n')
        g = collections.Counter(find(m) for m in ms)
        for m in ms: grp[m] = find(m)
        sizes = sorted(g.values(), reverse=True)
        print(f'{part}: {len(ms)} models, {npass} passing pairs, {len(g)} groups ({sum(1 for v in sizes if v > 1)} multi-model holding '
              f'{sum(v for v in sizes if v > 1)} models; largest {sizes[:6]})')
    fam = collections.defaultdict(dict)
    for m, (f, part, _) in models.items(): fam[f][part] = m
    full = [(grp[d['LTR']], grp[d['INT']]) for d in fam.values() if 'LTR' in d and 'INT' in d]
    e = set(full); lg = collections.Counter(a for a, b in e); ig = collections.Counter(b for a, b in e)
    print(f'families with both models {len(full)} -> distinct LTR-group/INT-group associations {len(e)}; '
          f'LTR groups {len(lg)} ({sum(v > 1 for v in lg.values())} with >1 INT group, max {max(lg.values())}); '
          f'INT groups {len(ig)} ({sum(v > 1 for v in ig.values())} with >1 LTR group, max {max(ig.values())})')
    for f in ('Os0365', 'Os0366', 'Os0119', 'Os0329', 'Os0850', 'Os0173', 'Os2463', 'Os0028', 'Os0592', 'Os0016', 'Os0419'):
        if 'LTR' in fam[f]: print(f'   {f}_LTR group size {sum(1 for m in grp if grp[m] == grp[fam[f]["LTR"]])}', end='')
    print()
