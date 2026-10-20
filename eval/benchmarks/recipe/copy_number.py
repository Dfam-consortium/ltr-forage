#!/usr/bin/env python3
"""Genome copy number per library element (LTR/INT/LTR record set).
usage: copy_number.py <library.fa> <genome_mask.out> [--min-frac 0.5] > table.tsv
A copy of a record = one RepeatMasker hit (fragments sharing an ID merged) covering
>= min-frac of the record.  Reports copies per record."""
import sys, collections, argparse
ap = argparse.ArgumentParser(); ap.add_argument('lib'); ap.add_argument('out')
ap.add_argument('--min-frac', type=float, default=0.5); a = ap.parse_args()
L = {}; n = None
for l in open(a.lib):
    if l.startswith('>'): n = l[1:].split()[0].split('#')[0]; L[n] = 0
    else: L[n] += len(l.strip())
cov = collections.Counter()   # (record, rm id) -> consensus bases covered
for l in open(a.out):
    f = l.split()
    if len(f) < 14 or not f[0].isdigit(): continue
    rec = f[9]
    if f[8] == '+': b, e = int(f[11]), int(f[12])
    else: b, e = int(f[13]), int(f[12])
    cov[(rec, f[-1] if f[-1] != '*' else f[-2])] += abs(e - b) + 1
cop = collections.Counter()
for (rec, _), c in cov.items():
    if rec in L and c >= a.min_frac * L[rec]: cop[rec] += 1
print('record\tlen\tcopies')
for r in L: print(f'{r}\t{L[r]}\t{cop[r]}')
