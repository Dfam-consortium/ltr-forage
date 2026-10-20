import statistics, math, collections
comp = str.maketrans('ACGT','TGCA'); K, STRIDE = 20, 10
def load(fn):
    d = {}
    for l in open(fn):
        c = l.split()
        if len(c) == 2: d[c[0]] = int(c[1])
    return d
G = load('tiles_genome.txt'); Sx = load('tiles_self.txt')
def cnt(d, q): return d.get(q, 0) + d.get(q.translate(comp)[::-1], 0)
seqs = {}; name = None
for l in open('internals.fa'):
    if l.startswith('>'): name = l[1:].strip()
    else: seqs[name] = l.strip()
feat = {}
for l in open('features20.tsv'):
    c = l.rstrip('\n').split('\t')
    if c[0] == 'idx': continue
    feat[c[0]] = dict(real=int(c[6]), domain=int(c[10]), order=c[11], motif=int(c[14]), tsd=int(c[13]), intl=int(c[8]))
def stats(counts, thr):
    n = len(counts); multi = [x >= thr for x in counts]
    frac = sum(multi) / n
    nb = 10; bins = [multi[i*n//nb:(i+1)*n//nb] for i in range(nb)]
    binfrac = sorted(sum(b)/max(1,len(b)) for b in bins if b)
    minbin = binfrac[0]; p20 = binfrac[max(0, len(binfrac)//5)]
    run = best = 0
    for m in multi:
        run = 0 if m else run + 1; best = max(best, run)
    gap = best * STRIDE / (n * STRIDE)
    logs = [math.log10(x) for x in counts if x >= thr]
    unif = (sum(1 for v in logs if abs(v - statistics.median(logs)) <= 0.5) / len(logs)) if logs else 0
    return frac, minbin, p20, gap, unif
rows = []
for i, s in seqs.items():
    tiles = [s[p:p+K] for p in range(0, len(s)-K+1, STRIDE)]
    tiles = [t for t in tiles if 'N' not in t]
    if len(tiles) < 20: continue
    cg = [cnt(G, t) for t in tiles]; cs = [cnt(Sx, t) for t in tiles]
    f = feat[i]
    gaps_bp = None
    rows.append((i, f, stats(cg, 3), stats(cs, 2), len(tiles)))
out = open('coverage_stats.tsv', 'w')
out.write('idx\treal\tdomain\torder\tmotif\ttsd\tintl\tntiles\tg_frac\tg_minbin\tg_p20\tg_gap\tg_unif\ts_frac\ts_minbin\ts_p20\ts_gap\ts_unif\n')
for i, f, g, s, nt in rows:
    out.write('\t'.join(map(str, [i, f['real'], f['domain'], f['order'], f['motif'], f['tsd'], f['intl'], nt] + list(g) + list(s))) + '\n')
out.close(); print('wrote', len(rows))
