#!/usr/bin/env python3
"""Base-level masking benchmark, as in Ou and Jiang 2018.

Usage: mask_benchmark.py <genome.fa.fai> <standard.out> <test.out> [--standard-class LTR] [--test-class LTR]

The standard annotation's rows of the given class (default: class/family
starting with LTR) define the LTR bases of the genome; the test annotation's
rows (same rule, or all rows with --test-class all) define the predicted
bases.  Both are unions of intervals per sequence.  Prints TP, FP, FN, TN in
bases and sensitivity TP/(TP+FN), specificity TN/(FP+TN), accuracy and
precision TP/(TP+FP).
"""
import sys, argparse, collections

ap = argparse.ArgumentParser()
ap.add_argument('fai'); ap.add_argument('standard'); ap.add_argument('test')
ap.add_argument('--standard-class', default='LTR')
ap.add_argument('--test-class', default='LTR')
a = ap.parse_args()

lengths = {l.split()[0]: int(l.split()[1]) for l in open(a.fai)}

def intervals(path, cls):
    iv = collections.defaultdict(list)
    for line in open(path):
        c = line.split()
        if len(c) < 15 or not c[0].isdigit(): continue
        if cls != 'all' and not c[10].startswith(cls): continue
        iv[c[4]].append((int(c[5]), int(c[6])))
    out = {}
    for chrom, v in iv.items():
        v.sort(); merged = []
        for s, e in v:
            if merged and s <= merged[-1][1] + 1: merged[-1][1] = max(merged[-1][1], e)
            else: merged.append([s, e])
        out[chrom] = merged
    return out

std = intervals(a.standard, a.standard_class)
tst = intervals(a.test, a.test_class)

def overlap_bases(A, B):
    i = j = 0; tot = 0
    while i < len(A) and j < len(B):
        s = max(A[i][0], B[j][0]); e = min(A[i][1], B[j][1])
        if s <= e: tot += e - s + 1
        if A[i][1] < B[j][1]: i += 1
        else: j += 1
    return tot

TP = FP = FN = TN = 0
for chrom, L in lengths.items():
    s = std.get(chrom, []); t = tst.get(chrom, [])
    sb = sum(e - b + 1 for b, e in s); tb = sum(e - b + 1 for b, e in t)
    tp = overlap_bases(s, t)
    TP += tp; FP += tb - tp; FN += sb - tp; TN += L - sb - (tb - tp)
G = sum(lengths.values())
print(f"genome {G:,} bp; standard LTR {TP+FN:,} bp ({100*(TP+FN)/G:.2f}%); test masked {TP+FP:,} bp ({100*(TP+FP)/G:.2f}%)")
print(f"TP {TP:,}  FP {FP:,}  FN {FN:,}  TN {TN:,}")
print(f"sensitivity {100*TP/max(1,TP+FN):.2f}%  specificity {100*TN/max(1,FP+TN):.2f}%  accuracy {100*(TP+TN)/max(1,G):.2f}%  precision {100*TP/max(1,TP+FP):.2f}%")
