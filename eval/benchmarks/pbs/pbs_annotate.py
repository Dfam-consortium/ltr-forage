#!/usr/bin/env python3
"""Per-element PBS annotation from ARAGORN tRNA calls (stand-in until ltr-forage does it).

usage: pbs_annotate.py <genome.fa> <aragorn_out_dir> <cand.tsv> > pbs.tsv

aragorn_out_dir holds <seq>__<offset>.fa.trna files (aragorn -t -i -gc1 -l -fasta
on 2 Mb chunks).  Probe = the tRNA's 3'-terminal 12 nt incl. CCA, reverse
complemented (starts TGG; TGG exact, <= 1 mismatch elsewhere), starting within
-5..+15 bp of the 5' LTR end, searched on both strands.  Output columns: element, pbs_mm12 (0/1, or . when none), pbs_trna, pbs_offset."""
import sys, glob, re, csv, collections, os

G, OUTDIR, TSV = sys.argv[1:4]
seqs = {}; n = None
for l in open(G):
    if l.startswith('>'): n = l[1:].split()[0]; seqs[n] = []
    else: seqs[n].append(l.strip().upper())
seqs = {k: ''.join(v) for k, v in seqs.items()}
RC = str.maketrans('ACGTN', 'TGCAN'); rc = lambda s: s.translate(RC)[::-1]
AA = {'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'}

genes = {}
for f in glob.glob(os.path.join(OUTDIR, '*.trna')):
    name = None
    for l in open(f):
        l = l.rstrip()
        if l.startswith('>'):
            m = re.match(r'>(?:\d+-\d+\s+)?(tRNA-\S+)\s+(c?)\[(\d+),(\d+)\]', l); name = m.group(1) if m else None; key = (f, l); genes[key] = [name, '']
        elif name and re.fullmatch(r'[acgtnACGTN]+', l): genes[key][1] += l.upper()
probes = collections.defaultdict(set)
for name, g in genes.values():
    if name and 70 <= len(g) <= 120 and name.split('(')[0][5:] in AA:
        probes[rc(g[-15:] + 'CCA')[:12]].add(name[5:])
K = 12
idx = collections.defaultdict(set)
for p, names in probes.items():
    t = ','.join(sorted(names)); idx[p].add((0, t))
    for i in range(3, K):
        for b in 'ACGT':
            if b != p[i]: idx[p[:i] + b + p[i + 1:]].add((1, t))
print(f'{len(probes)} distinct 12-nt probes', file=sys.stderr)

def hit(w, base):
    best = None
    for j in range(len(w) - K + 1):
        for mm, t in idx.get(w[j:j + K], ()):
            if best is None or mm < best[0]: best = (mm, t, base + j)
    return best

print('element\tpbs_mm12\tpbs_trna\tpbs_offset')
for x in csv.DictReader(open(TSV), delimiter='\t'):
    c = x['seq']; s = seqs[c]; le, rs = int(x['lltr_end']), int(x['rltr_start'])
    # both strands always: dom_strand follows the majority of domain hits and is wrong for some
    # Gypsy families whose accessory ORFs sit on the other strand (ltr-forage does the same)
    plus = hit(s[le - 5:le + 15 + K], -5)
    minus = hit(rc(s[rs - 1 - 15 - K:rs - 1 + 5]), -5)
    r = min((h for h in (plus, minus) if h), key=lambda h: h[0], default=None)
    print(f"{c}:{x['elem_start']}..{x['elem_end']}\t" + (f'{r[0]}\t{r[1]}\t{r[2]}' if r else '.\t.\t.'))
