#!/usr/bin/env python3
"""Turn LTR_FINDER -w2 table output into a library FASTA of LTR, INT, LTR per
element, named like ltr-forage's library.

Usage: ltrfinder_to_library.py <genome.fa> <table.w2> [table.w2 ...] > library.fa
"""
import sys, subprocess

genome = sys.argv[1]
def fetch(chrom, s, e):
    out = subprocess.run(['samtools', 'faidx', genome, f'{chrom}:{s}-{e}'], capture_output=True, text=True, check=True).stdout
    return ''.join(out.split('\n')[1:]).upper()

n = 0
for path in sys.argv[2:]:
    for line in open(path):
        if not line.startswith('['): continue
        c = line.rstrip('\n').split('\t')
        chrom = c[1]; s, e = map(int, c[2].split('-')); l1, l2 = map(int, c[3].split(','))
        recs = [(s, s + l1 - 1, 'LTR'), (s + l1, e - l2, 'INT'), (e - l2 + 1, e, 'LTR')]
        for a, b, kind in recs:
            if b < a: continue
            print(f'>{chrom}:{a}..{b}_{kind}#LTR/unknown ltr_finder score={c[13]} strand={c[12]}')
            print(fetch(chrom, a, b))
        n += 1
print(f'{n} elements', file=sys.stderr)
