import sys, json, re, collections
clus = json.load(open(sys.argv[1])); jiang = sys.argv[2]; libs = sys.argv[3:]
size = collections.Counter(clus.values())
best = {}
for l in open(jiang):
    q, s, pid, ln, qs, qe, ql, bs = l.split()[:8]
    if float(pid) >= 80 and int(ln) >= 0.3*int(ql) and (q not in best or float(bs) > best[q][1]):
        best[q] = (s.split('_')[0], float(bs))
by_chr = collections.defaultdict(list)
for e, c in clus.items():
    m = re.match(r'(\S+):(\d+)\.\.(\d+)', e); by_chr[m[1]].append((int(m[2]), int(m[3]), c))
def cl(r):
    m = re.match(r'(\S+):(\d+)\.\.(\d+)_', r)
    for a, b, c in by_chr[m[1]]:
        if a <= int(m[2]) and int(m[3]) <= b: return c
for lib in libs:
    recs = [l[1:].split()[0] for l in open(lib) if l.startswith('>')]
    fam = collections.defaultdict(set); cs = set()
    for r in recs:
        c = cl(r); cs.add(c)
        if r in best: fam[c].add((0 if '_INT' in r else 1, best[r][0]))
    def f(c): return min(fam[c])[1] if fam[c] else None
    multi_f = {f(c) for c in cs if size[c] > 1} - {None}
    sing = [c for c in cs if size[c] == 1]
    a = sum(1 for c in sing if f(c) in multi_f); b = len({f(c) for c in sing if f(c) and f(c) not in multi_f})
    b_n = sum(1 for c in sing if f(c) and f(c) not in multi_f); u = sum(1 for c in sing if f(c) is None)
    print(f'{lib}: {len(cs)} clusters; multi-element clusters name {len(multi_f)} curated families; '
          f'singletons {len(sing)}: {a} of a family that has a multi-element cluster, '
          f'{b_n} from {b} families seen only as singletons, {u} matching no curated family')
