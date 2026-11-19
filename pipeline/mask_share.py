#!/usr/bin/env python3
"""Report-only outlier check (end of step 23b): how much of the genome each library record would mask.

usage: mask_share.py <seeds/hits.tsv> <families.tsv> <scale> <genome_mb> [--warn PCT]

Step 23b's search of library.fa against the depth sample (seeds/hits.tsv: qseqid sseqid sstart send bitscore,
at most 50 HSPs per record and sample window) is read as a masking run: each sample base goes to the
best-scoring hit over it, as RepeatMasker's competition would give it, and each record's bases are divided by
the sample's fraction of the genome (<scale>).  families.tsv gets the column est_genome_pct; records at or
above --warn (2.0) percent of the non-N genome are listed in the log line.  Nothing is dropped or changed.

Why: the known failure on SINE- and LINE-rich genomes is one record made of L1 or SINE sequence that masks a
large share of the genome (mouse ablation runs, 2026-10-31: one record each masking 92-147 Mb, 3.5-5.5% of
the genome by this estimate, 60-90% of the run's false-positive bases).  On the six tuning libraries the
largest record is at most 1.2% (lotus); 2% was set from these nine libraries once (2026-11-19) and is not
tuned further.  A warning means "look at this record", not "this record is wrong": a genome made largely of a
few LTR families has real records above it (maize: 3 of 2,268 records, all one Gypsy family, up to 6.2%).
"""
import argparse, collections, csv, os
import numpy as np

ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
ap.add_argument('hits'); ap.add_argument('families'); ap.add_argument('scale', type=float); ap.add_argument('genome_mb', type=float)
ap.add_argument('--warn', type=float, default=2.0, help='warn at this percent of the non-N genome (2.0; 0: no warning)')
a = ap.parse_args()

by = collections.defaultdict(list)
if os.path.exists(a.hits):
    for l in open(a.hits):
        q, s, ss, se, b = l.rstrip('\n').split('\t')[:5]; ss, se = int(ss), int(se)
        by[s].append((float(b), min(ss, se), max(ss, se), q.split('#')[0]))
mass = collections.Counter()
for hs in by.values():
    own = np.full(max(h[2] for h in hs) + 1, -1, dtype=np.int32); names = []
    for b, lo, hi, q in sorted(hs, reverse=True):   # best score first; a base keeps the first hit that claims it
        seg = own[lo:hi + 1]; free = seg < 0
        if free.any(): names.append(q); seg[free] = len(names) - 1
    for i, c in zip(*np.unique(own[own >= 0], return_counts=True)): mass[names[i]] += int(c)
pct = {r: 100 * m / a.scale / (a.genome_mb * 1e6) for r, m in mass.items()}

rows = list(csv.DictReader(open(a.families), delimiter='\t'))
fields = list(rows[0].keys()) if rows else ['family', 'new_name', 'origin', 'kind', 'len', 'status', 'class']
if 'est_genome_pct' not in fields: fields.append('est_genome_pct')
tmp = a.families + '.tmp'
with open(tmp, 'w', newline='') as fo:
    w = csv.DictWriter(fo, fieldnames=fields, delimiter='\t', lineterminator='\n'); w.writeheader()
    for r in rows:
        r['est_genome_pct'] = f"{pct.get(r['new_name'], 0.0):.3f}" if r.get('status') == 'kept' else ''
        w.writerow(r)
os.replace(tmp, a.families)

kept = [r['new_name'] for r in rows if r.get('status') == 'kept']
top = max(kept, key=lambda r: pct.get(r, 0.0)) if kept else None
flag = sorted((r for r in kept if a.warn > 0 and pct.get(r, 0.0) >= a.warn), key=lambda r: -pct[r])
msg = (f'largest record {top} est. {pct.get(top, 0.0):.2f}% of the genome' if top else 'no kept records')
if flag:
    msg = (f'WARNING: {len(flag)} record(s) estimated to mask >= {a.warn:g}% of the genome, look at them before use '
           f'(L1 or SINE sequence called LTR is the known failure; very abundant real LTR families also reach it): '
           + ', '.join(f'{r} {pct[r]:.2f}%' for r in flag[:10]) + (' ...' if len(flag) > 10 else '') + f'; {msg}')
else:
    msg = f'no record estimated to mask >= {a.warn:g}% of the genome; {msg}' if a.warn > 0 else msg
print(msg)
