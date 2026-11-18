#!/usr/bin/env python3
"""Cut out of the models what the run's poly-A-terminated repeats cover (ltr_library.py step 22b, before the depth check).
usage: polya_cut.py <models.fa> <polya.fa> <out.fa> <workdir> [--min-depth 10] [--min-len 50] [--cuts cuts.tsv]

polya.fa: the LTRs of the poly-A-vetoed candidates sampled by step 11b (polya_family.py), copies of the genome's SINE
and LINE families.  Until step 11b the elements those repeats flank were kept and their "LTR" records told the
nested-LTR vote (step 15) where an L1 or SINE lay inside a real element's internal region.  With those elements
vetoed that evidence is gone, and the depth check cannot replace it where the insertion is most of the record: it
cuts runs at a multiple of the record's median depth (mm10: an MTA internal of 3,594 bp, 2,807 of them L1, masked
184 Mb outside the curated LTRs).  Here every run of >= --min-len positions covered by >= --min-depth of the sampled
LTRs (rmb_ava.sh, word 12, the search of step 11b) is cut out and the flanks joined; a record with less than --min-len
left is dropped.  cuts.tsv: record<TAB>i-j,i-j,... (1-based, inclusive) or record<TAB>drop."""
import sys, os, argparse, subprocess
import numpy as np
ap = argparse.ArgumentParser()
for p in ('models', 'polya', 'out', 'work'): ap.add_argument(p)
ap.add_argument('--min-depth', type=int, default=10); ap.add_argument('--min-len', type=int, default=50); ap.add_argument('--cuts')
a = ap.parse_args()
HERE = os.path.dirname(os.path.realpath(__file__))
recs = []; n = None
for l in open(a.models):
    if l.startswith('>'): recs.append([l[1:].rstrip('\n'), []])
    else: recs[-1][1].append(l.strip())
recs = [(h, ''.join(s)) for h, s in recs]
key = lambda h: h.split()[0].split('#')[0]
os.makedirs(a.work, exist_ok=True)
with open(os.path.join(a.work, 'models.fa'), 'w') as fo:
    for h, s in recs: fo.write(f'>{key(h)}\n{s}\n')
env = dict(os.environ, WORD='12'); env.pop('MATRIX', None); env.pop('MIN_PID', None)
subprocess.run(['bash', os.path.join(HERE, 'rmb_ava.sh'), 'models.fa', os.path.abspath(a.polya), '1', 'cut_hits.tsv',
                *'qseqid sseqid pident length qstart qend sstart send qlen slen bitscore'.split()],
               cwd=a.work, env=env, check=True, stdout=sys.stderr)
depth = {}
for l in open(os.path.join(a.work, 'cut_hits.tsv')):
    q, s, pid, ln, qs, qe, ss, se, ql, sl, b = l.split('\t')
    d = depth.get(q)
    if d is None: d = depth[q] = np.zeros(int(ql) + 2, dtype=np.int32)
    lo, hi = sorted((int(qs), int(qe))); d[lo] += 1; d[hi + 1] -= 1
ncut = bp = gone = 0
cuts = open(a.cuts, 'w') if a.cuts else None
with open(a.out + '.tmp', 'w') as fo:
    for h, s in recs:
        d = depth.get(key(h))
        if d is not None:
            c = np.cumsum(d)[1:-1] >= a.min_depth
            e = np.flatnonzero(np.diff(np.concatenate(([0], c.astype(np.int8), [0]))))
            runs = [(int(x), int(y)) for x, y in zip(e[::2], e[1::2]) if y - x >= a.min_len]   # 0-based half-open
            if runs:
                keep = np.ones(len(s), dtype=bool)
                for x, y in runs: keep[x:y] = False
                t = ''.join(ch for ch, k in zip(s, keep) if k)
                if len(t) < a.min_len:
                    gone += 1; bp += len(s)
                    if cuts: cuts.write(f'{key(h)}\tdrop\n')
                    continue
                ncut += 1; bp += len(s) - len(t); s = t
                if cuts: cuts.write(f'{key(h)}\t' + ','.join(f'{x + 1}-{y}' for x, y in runs) + '\n')
        fo.write(f'>{h}\n{s}\n')
os.replace(a.out + '.tmp', a.out)
if cuts: cuts.close()
print(f'{ncut} model(s) cut where >= {a.min_depth} sampled poly-A-terminated repeats cover them, {bp:,} bases removed, {gone} dropped')
