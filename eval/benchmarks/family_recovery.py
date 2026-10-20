#!/usr/bin/env python3
"""Family-level recovery of a curated LTR library by a candidate library.

Usage: family_recovery.py <curated.fa> <library.fa> [--min-identity 80] [--out table.tsv]

Every curated sequence (named <family>_LTR... or <family>_INT...) is aligned
against the candidate library with blastn (dc-megablast).  For each curated
sequence the script takes the union of the HSPs at or above the identity cut
and reports the fraction of it covered, then per family the best LTR and
best INT coverage.  The summary counts families whose LTR, INT, or both reach
50%, 80% and 95% coverage.  Weighting is per family, not per base.
"""
import sys, os, re, subprocess, tempfile, collections, argparse

ap = argparse.ArgumentParser()
ap.add_argument('curated'); ap.add_argument('library')
ap.add_argument('--min-identity', type=float, default=80.0)
ap.add_argument('--out', default=None)
ap.add_argument('--task', default='dc-megablast')
a = ap.parse_args()

def read_fasta(p):
    seqs, name = {}, None
    for l in open(p):
        if l.startswith('>'):
            name = l[1:].split()[0]; seqs[name] = []
        elif name: seqs[name].append(l.strip())
    return {k: ''.join(v) for k, v in seqs.items()}

def family_kind(name):
    """Family and part (LTR or INT) from a curated name.  Handles the Jiang
    libraries (Os0008_LTR, Os0020_INT_Bajie) and Repbase (ATCOPIA44LTR,
    ATGP8I, ATLANTYS2_I, Copia-3_AT-LTR)."""
    base = name.split('#')[0]
    parts = base.split('_')
    if 'LTR' in parts[1:]: return parts[0], 'LTR'
    if 'INT' in parts[1:]: return parts[0], 'INT'
    m = re.match(r'^(.+?)[-_]?LTR$', base)
    if m: return m.group(1), 'LTR'
    m = re.match(r'^(.+?)[-_]?(INT|int|I)$', base)
    if m: return m.group(1), 'INT'
    return base, 'other'

cur = read_fasta(a.curated)
with tempfile.TemporaryDirectory() as td:
    db = os.path.join(td, 'lib')
    subprocess.run(['makeblastdb', '-in', a.library, '-dbtype', 'nucl', '-out', db], check=True, capture_output=True)
    res = subprocess.run(['blastn', '-task', a.task, '-query', a.curated, '-db', db, '-evalue', '1e-10',
                          '-max_target_seqs', '100', '-max_hsps', '50', '-num_threads', '8',
                          '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send bitscore'],
                         check=True, capture_output=True, text=True).stdout
hsps = collections.defaultdict(list)
for l in res.splitlines():
    q, s, pid, ln, qs, qe, ss, se, bits = l.split('\t')
    if float(pid) < a.min_identity: continue
    hsps[q].append((int(qs), int(qe), float(pid)))

def coverage(name):
    L = len(cur[name]); iv = sorted((min(s, e), max(s, e)) for s, e, _ in hsps.get(name, []))
    cov, cur_s, cur_e = 0, None, None
    for s, e in iv:
        if cur_s is None: cur_s, cur_e = s, e
        elif s <= cur_e + 1: cur_e = max(cur_e, e)
        else: cov += cur_e - cur_s + 1; cur_s, cur_e = s, e
    if cur_s is not None: cov += cur_e - cur_s + 1
    return cov / L if L else 0.0

fam = collections.defaultdict(lambda: {'LTR': 0.0, 'INT': 0.0, 'has': set()})
for name in cur:
    f, k = family_kind(name)
    if k == 'other': continue
    fam[f]['has'].add(k)
    fam[f][k] = max(fam[f][k], coverage(name))
rows = sorted(fam.items())
if a.out:
    with open(a.out, 'w') as o:
        o.write('family\tltr_cov\tint_cov\n')
        for f, d in rows: o.write(f"{f}\t{d['LTR']:.3f}\t{d['INT']:.3f}\n")
n = len(rows); n_ltr = sum(1 for _, d in rows if 'LTR' in d['has']); n_int = sum(1 for _, d in rows if 'INT' in d['has'])
print(f"{os.path.basename(a.library)} vs {os.path.basename(a.curated)}: {n} families ({n_ltr} with an LTR, {n_int} with an INT), identity >= {a.min_identity:.0f}%")
for t in (0.5, 0.8, 0.95):
    l = sum(1 for _, d in rows if 'LTR' in d['has'] and d['LTR'] >= t)
    i = sum(1 for _, d in rows if 'INT' in d['has'] and d['INT'] >= t)
    b = sum(1 for _, d in rows if 'LTR' in d['has'] and 'INT' in d['has'] and d['LTR'] >= t and d['INT'] >= t)
    print(f"  coverage >= {t:.2f}: LTR {l}/{n_ltr} ({100*l/max(1,n_ltr):.1f}%)  INT {i}/{n_int} ({100*i/max(1,n_int):.1f}%)  both {b}")
miss = [f for f, d in rows if d['LTR'] < 0.5 and d['INT'] < 0.5]
print(f"  families with neither LTR nor INT half covered: {len(miss)}")
