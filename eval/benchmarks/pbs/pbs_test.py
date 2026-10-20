#!/usr/bin/env python3
"""Does a tRNA-derived PBS separate real LTR elements from junk?

usage: pbs_test.py <genome.fa> <trna_calls.json> <features.tsv> <lib_vs_curated.tsv> <singleton_classes.json> [--max-edits 2]

Probes: for each ARAGORN tRNA call (70-120 bp, standard amino acid), the
reverse complement of its last 15 nt + CCA (the 18 nt a tRNA primer pairs
with).  A PBS hit = a probe aligned (semi-global, edit distance <= max-edits)
starting within -5..+15 bp of the 5' LTR end.  Matching uses pigeonhole
seeding (three 6-mers; with <= 2 edits one is exact) and then a full DP.

Groups: kept elements whose LTR matches a curated LTR (>=80% id over >=80% of
the LTR), kept elements with no curated LTR match, singleton classes A-G,
and two controls: the same window 300 bp further into the internal region,
and random genome windows."""
import sys, json, csv, random, collections, argparse, re

ap = argparse.ArgumentParser()
for p in ('genome', 'trna', 'tsv', 'curated_hits', 'classes'): ap.add_argument(p)
ap.add_argument('--max-edits', type=int, default=2)
a = ap.parse_args()
E = a.max_edits

seqs = {}; n = None
for l in open(a.genome):
    if l.startswith('>'): n = l[1:].split()[0]; seqs[n] = []
    else: seqs[n].append(l.strip().upper())
seqs = {k: ''.join(v) for k, v in seqs.items()}
RC = str.maketrans('ACGTN', 'TGCAN')
rc = lambda s: s.translate(RC)[::-1]

AA = {'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'}
probes = collections.defaultdict(set)   # probe -> {tRNA-Xxx(anticodon)}
for chrom, s, e, strand, name, gene in json.load(open(a.trna)):
    if not 70 <= len(gene) <= 120 or name.split('(')[0][5:] not in AA: continue
    probes[rc(gene[-15:] + 'CCA')].add(name)
print(f'{len(probes)} distinct probes from tRNA calls of 70-120 bp', file=sys.stderr)
seed = collections.defaultdict(list)    # 6-mer -> [(probe, offset in probe)]
for p in probes:
    for o in (0, 6, 12): seed[p[o:o + 6]].append((p, o))

def edist_semiglobal(p, w):
    """min edit distance of all of p against any substring of w; returns (d, end)."""
    prev = [0] * (len(w) + 1)
    for i in range(1, len(p) + 1):
        cur = [i] + [0] * len(w)
        for j in range(1, len(w) + 1):
            cur[j] = min(prev[j - 1] + (p[i - 1] != w[j - 1]), prev[j] + 1, cur[j - 1] + 1)
        prev = cur
    d = min(prev); return d, prev.index(d)

def best_pbs(chrom, ltr_end):
    """ltr_end: 1-based last base of the 5' LTR on the + strand."""
    s = seqs[chrom]; lo = max(0, ltr_end - 5); w = s[lo:ltr_end + 15 + 18 + E]
    best = None
    tried = set()
    for j in range(len(w) - 5):
        for p, o in seed.get(w[j:j + 6], ()):
            st = j - o                       # implied probe start in w
            if (p, st) in tried or not (-E <= st <= 20 + E): continue
            tried.add((p, st))
            sub = w[max(0, st - E): st + 18 + E]
            d, _ = edist_semiglobal(p, sub)
            if d <= E and (best is None or d < best[0]): best = (d, p, lo + st - ltr_end)
    return best

def best_pbs_minus(chrom, ltr_start):
    """element on the - strand: PBS lies before the right-hand LTR's start; search the reverse complement."""
    s = seqs[chrom]; lo = max(0, ltr_start - 1 - 15 - 18 - E); hi = ltr_start - 1 + 5
    w = rc(s[lo:hi]); best = None; tried = set()
    for j in range(len(w) - 5):
        for p, o in seed.get(w[j:j + 6], ()):
            st = j - o
            if (p, st) in tried or not (-E <= st <= 20 + E): continue
            tried.add((p, st))
            d, _ = edist_semiglobal(p, w[max(0, st - E): st + 18 + E])
            if d <= E and (best is None or d < best[0]): best = (d, p, st - 5)
    return best

def better(r1, r2):
    if r1 is None: return r2
    if r2 is None: return r1
    return r1 if r1[0] <= r2[0] else r2


# Short test: the tRNA's 3'-terminal 12 nt (probe prefix, starts TGG), <= 1 mismatch.
# Many families pair only over these 12-13 nt (e.g. rice iMet users: TGGTATCAGAGC then diverge).
K = 12
short = collections.defaultdict(set)    # 12-mer (incl. 1-mismatch variants) -> {(mm, tRNA names)}
for p, names in probes.items():
    k = p[:K]
    short[k].add((0, tuple(sorted(names))))
    for i in range(3, K):                # keep the TGG exact
        for b in 'ACGT':
            if b != k[i]: short[k[:i] + b + k[i + 1:]].add((1, tuple(sorted(names))))

def short_hit(w, base):
    """best (mm, names, offset) for a 12-mer starting at w[j]; offset = base + j."""
    best = None
    for j in range(len(w) - K + 1):
        for mm, names in short.get(w[j:j + K], ()):
            if best is None or mm < best[0]: best = (mm, names, base + j)
    return best

def short_both(c, lltr_end, rltr_start):
    s = seqs[c]
    a = short_hit(s[lltr_end - 5:lltr_end + 15 + K], -5)
    b = short_hit(rc(s[rltr_start - 1 - 15 - K:rltr_start - 1 + 5]), -5)
    return better(a, b)

# curated matches per LTR record (query coverage >= 80%, identity >= 80%)
cur = collections.defaultdict(float)
for l in open(a.curated_hits):
    q, s, pid, ln, qs, qe, ql, bits = l.split('\t')[:8]
    if '_LTR' in q.split('#')[0] and '_LTR' in s and float(pid) >= 80:
        cur[q.split('#')[0]] = max(cur[q.split('#')[0]], int(ln) / int(ql))
cls = json.load(open(a.classes))
prim = collections.Counter(v['prim'] for v in cls.values())

rows = list(csv.DictReader(open(a.tsv), delimiter='\t'))
groups = collections.defaultdict(list); g12 = collections.defaultdict(list); per_el = []
for x in rows:
    if x['keep'] != '1': continue
    c = x['seq']; k = f"{c}:{x['elem_start']}..{x['elem_end']}"
    l5 = f"{c}:{x['lltr_start']}..{x['lltr_end']}_LTR"
    # orientation is often unknown, so both strands are searched everywhere (controls too)
    res = better(best_pbs(c, int(x['lltr_end'])), best_pbs_minus(c, int(x['rltr_start'])))
    ctrl = better(best_pbs(c, int(x['lltr_end']) + 300), best_pbs_minus(c, int(x['rltr_start']) - 300))
    g = 'kept, LTR matches curated' if cur.get(l5, 0) >= 0.8 else 'kept, no curated LTR match'
    r12 = short_both(c, int(x['lltr_end']), int(x['rltr_start']))
    c12 = short_both(c, int(x['lltr_end']) + 300, int(x['rltr_start']) - 300)
    groups[g].append(res); groups['control: +300 bp into internal'].append(ctrl)
    g12[g].append(r12); g12['control: +300 bp into internal'].append(c12)
    per_el.append((k, g, res, r12))
    if k in cls and cls[k]['single']:
        c0 = cls[k]['cat'][0]
        gg = 'singleton ' + ('G sole' if c0 == 'G' and prim[cls[k]['prim']] == 1 else c0)
        groups[gg].append(res); g12[gg].append(r12)
    elif k in cls:
        groups['multi-element cluster exemplar'].append(res); g12['multi-element cluster exemplar'].append(r12)
random.seed(1)
names = [k for k in seqs if len(seqs[k]) > 100000]
for _ in range(5000):
    c = random.choice(names); p = random.randrange(1000, len(seqs[c]) - 1000)
    if 'N' in seqs[c][p - 10:p + 50]: continue
    q = random.randrange(1000, len(seqs[c]) - 1000)
    if 'N' in seqs[c][q - 50:q + 10]: continue
    groups['control: random genome window'].append(better(best_pbs(c, p), best_pbs_minus(c, q)))
    g12['control: random genome window'].append(better(short_hit(seqs[c][p - 5:p + 15 + K], -5), short_hit(rc(seqs[c][q - 1 - 15 - K:q - 1 + 5]), -5)))

order = ['kept, LTR matches curated', 'kept, no curated LTR match', 'multi-element cluster exemplar',
         'singleton G sole', 'singleton G', 'singleton F', 'singleton E', 'singleton A', 'singleton B', 'singleton C', 'singleton D',
         'control: +300 bp into internal', 'control: random genome window']
print(f"{'group':34}{'n':>6}{'18nt<=1':>9}{'18nt<=' + str(E):>9}{'12nt=0':>9}{'12nt<=1':>9}{'either':>9}   top tRNA (12nt)")
for g in order:
    v = groups.get(g, []); v12 = g12.get(g, [])
    if not v: continue
    n = len(v)
    o = sum(1 for r in v if r and r[0] <= 1); t = sum(1 for r in v if r)
    z12 = sum(1 for r in v12 if r and r[0] == 0); o12 = sum(1 for r in v12 if r)
    ei = sum(1 for r, r2 in zip(v, v12) if r or r2)
    top = collections.Counter(r[1][0][5:] for r in v12 if r).most_common(3)
    print(f'{g:34}{n:>6}{100*o/n:>8.1f}%{100*t/n:>8.1f}%{100*z12/n:>8.1f}%{100*o12/n:>8.1f}%{100*ei/n:>8.1f}%   ' + ', '.join(f'{k} {c}' for k, c in top))

with open('pbs_per_element.tsv', 'w') as o:
    o.write('element\tgroup\tedits\tprobe\ttrnas\toffset\tmm12\ttrnas12\toffset12\n')
    for k, g, r, r2 in per_el:
        o.write(f"{k}\t{g}\t" + (f"{r[0]}\t{r[1]}\t{','.join(sorted(probes[r[1]]))}\t{r[2]}" if r else '.\t.\t.\t.')
                + '\t' + (f"{r2[0]}\t{','.join(r2[1])}\t{r2[2]}" if r2 else '.\t.\t.') + '\n')
