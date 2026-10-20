#!/usr/bin/env python3
"""Add tRNA / primer binding site evidence to cases.json.

usage: add_trna.py <cases.json> <genome.fa> <aragorn_out_dir> <pbs.tsv>

aragorn_out_dir: *.trna files (ARAGORN -fasta/-fon, as for pbs_annotate.py).
pbs.tsv: pbs_annotate.py output for every kept element (for genome-wide primer usage).
Per case adds 'trna': the best 12-nt PBS hit (TGG exact, <= 1 mismatch elsewhere, starting
-5..+15 from the 5' LTR end, either strand), the site, the primer tRNA's 3' end, pairing over
18 nt, the tRNA's genome copy number and how many kept elements use the same primer."""
import sys, glob, re, json, collections, os

CASES, G, OUTDIR, PBS = sys.argv[1:5]
K, K18 = 12, 18
fai = {l.split()[0]: list(map(int, l.split()[1:5])) for l in open(G + '.fai')}
gf = open(G, 'rb')
def get(c, s, e):  # 1-based inclusive
    L, off, lb, lB = fai[c]; s = max(1, s); e = min(L, e); s -= 1
    a = off + s // lb * lB + s % lb; b = off + (e - 1) // lb * lB + (e - 1) % lb
    gf.seek(a); return gf.read(b - a + 1).replace(b'\n', b'').decode().upper()
RC = str.maketrans('ACGTN', 'TGCAN'); rc = lambda s: s.translate(RC)[::-1]
AA = {'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'}

# tRNA genes: probe12 -> types, copies, a representative 3' end (last 15 nt + CCA)
genes = []
for f in glob.glob(os.path.join(OUTDIR, '*.trna')):
    cur = None
    for l in open(f):
        l = l.rstrip()
        if l.startswith('>'):
            m = re.match(r'>(?:\d+-\d+\s+)?tRNA-(\S+)', l); cur = [m.group(1) if m else None, '']; genes.append(cur)
        elif cur and re.fullmatch(r'[acgtnACGTN]+', l): cur[1] += l.upper()
types = collections.defaultdict(collections.Counter); tail = {}; ncopy = collections.Counter()
total_genes = 0
for name, g in genes:
    if not name or not (70 <= len(g) <= 120) or name.split('(')[0] not in AA: continue
    total_genes += 1
    t3 = g[-15:] + 'CCA'; p = rc(t3)[:K]
    types[p][name] += 1; ncopy[p] += 1; tail.setdefault(p, collections.Counter())[t3] += 1
usage = collections.Counter(); n_elem = 0; n_pbs = 0
for l in list(open(PBS))[1:]:
    e, mm, t, off = l.rstrip('\n').split('\t'); n_elem += 1
    if mm != '.': n_pbs += 1; usage[t] += 1
# one-mismatch index (TGG fixed)
idx = collections.defaultdict(set)
for p in types:
    idx[p].add((0, p))
    for i in range(3, K):
        for b in 'ACGT':
            if b != p[i]: idx[p[:i] + b + p[i + 1:]].add((1, p))

def scan(w, base):
    best = None
    for j in range(len(w) - K + 1):
        for mm, p in idx.get(w[j:j + K], ()):
            if best is None or mm < best[0]: best = (mm, p, base + j)
    return best

cases = json.load(open(CASES))
for c in cases:
    le, rs = c['lltr'][1], c['rltr'][0]
    plus_w = get(c['seq'], le + 1 - 5, le + 15 + K + 10)              # internal side of the 5' LTR, + strand
    minus_w = rc(get(c['seq'], rs - 15 - K - 10, rs - 1 + 5))          # same on the - strand (3' LTR is the 5' LTR there)
    hits = [(h, s, w) for h, s, w in ((scan(plus_w[:5 + 15 + K], -5), '+', plus_w), (scan(minus_w[:5 + 15 + K], -5), '-', minus_w)) if h]
    if not hits:
        c['trna'] = None; continue
    (mm, p, off), strand, w = min(hits, key=lambda x: x[0][0])
    site = w[off + 5: off + 5 + K18]
    t3 = tail[p].most_common(1)[0][0]                 # tRNA 3' end, 18 nt incl. CCA
    comp = rc(t3)                                     # what a perfect PBS reads
    pair = ''.join('|' if a == b else ' ' for a, b in zip(site, comp))
    names = types[p]
    label = ','.join(sorted({n for n in names}))
    c['trna'] = {
        'strand': strand, 'offset': off, 'mm12': mm, 'match18': pair.count('|'),
        'site': site, 'context': w[:5 + 15 + K + 10], 'site_at': off + 5,
        'trna_3end': t3, 'trna_types': dict(names), 'probe12': p,
        'copies': ncopy[p], 'uses': usage.get(label, 0),
    }
json.dump({'cases': cases, 'trna_summary': {'genes': total_genes, 'probes': len(types), 'elements': n_elem, 'with_pbs': n_pbs,
           'top': usage.most_common(8)}}, open(CASES.replace('.json', '_trna.json'), 'w'), indent=1)
for c in cases:
    t = c['trna']; print(c['id'], t and (t['strand'], t['offset'], t['mm12'], t['match18'], t['trna_types'], t['copies'], t['uses']))
