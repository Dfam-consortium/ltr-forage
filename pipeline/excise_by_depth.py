#!/usr/bin/env python3
"""Nested-insertion cleanup from genome hit depth.

Usage: excise_by_depth.py <lib.fa> <genome_mask.out> <out.fa> [--ratio 10] [--min-depth 50] [--min-len 100] [--max-frac 0.7]

genome_mask.out is a RepeatMasker run of the genome with lib.fa as the
library.  For each record, the number of genomic hits covering each record
position (from the "position in repeat" columns) is its copy-depth profile.
Positions whose depth is at least --ratio times the record's median depth
and at least --min-depth are nested high-copy repeats; runs of --min-len or
longer are cut out unless they cover more than --max-frac of the record.
"""
import sys, argparse, collections, statistics
ap = argparse.ArgumentParser()
ap.add_argument('lib'); ap.add_argument('mask'); ap.add_argument('out')
ap.add_argument('--ratio', type=float, default=10); ap.add_argument('--min-depth', type=int, default=50)
ap.add_argument('--min-len', type=int, default=100); ap.add_argument('--max-frac', type=float, default=0.7)
ap.add_argument('--abs-depth', type=int, default=0, help='a position at or above this depth is nested whatever the record median (0: off)')
ap.add_argument('--drop-frac', type=float, default=0, help='drop a record whose high-depth runs cover more than this fraction of it (0: off)')
ap.add_argument('--scale', type=float, default=1.0, help='divide depth thresholds by this (a genome sample of 10%% needs 0.1)')
ap.add_argument('--profile', action='store_true', help='mask is a depth TSV (record<TAB>d0,d1,...; e.g. from kmer_depth) instead of a RepeatMasker .out')
ap.add_argument('--consensus', help='fasta of consensus records (stage B representatives): cut only runs at --consensus-ratio '
                'times their median depth, never by --abs-depth or --drop-frac; a big family\'s consensus is high-copy itself')
ap.add_argument('--consensus-ratio', type=float, default=10)
ap.add_argument('--cuts', help='write the changes here: record<TAB>drop, or record<TAB>i-j,i-j,... (0-based half-open runs cut out)')
a = ap.parse_args()
a.min_depth = a.min_depth * a.scale; a.abs_depth = a.abs_depth * a.scale
seqs = {}; order = []; name = None
for l in open(a.lib):
    if l.startswith('>'): name = l[1:].split()[0].split('#')[0]; seqs[name] = []; order.append(name)
    else: seqs[name].append(l.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}
cons = {l[1:].split()[0].split('#')[0] for l in open(a.consensus) if l.startswith('>')} if a.consensus else set()
depth = {k: [0] * (len(v) + 2) for k, v in seqs.items()}
prof_in = {}
if a.profile:
    for line in open(a.mask):
        n, v = line.rstrip('\n').split('\t')
        if n in seqs: prof_in[n] = list(map(int, v.split(','))) if v else []
for line in ([] if a.profile else open(a.mask)):
    c = line.split()
    if len(c) < 15 or not c[0].isdigit(): continue
    name = c[9]
    if name not in depth: continue
    if c[8] == '+': rs, re_ = int(c[11]), int(c[12])
    else: rs, re_ = int(c[13]), int(c[12])
    d = depth[name]; L = len(seqs[name])
    rs = max(1, min(rs, L)); re_ = max(1, min(re_, L))
    if rs > re_: rs, re_ = re_, rs
    d[rs] += 1; d[re_ + 1] -= 1
cuts = open(a.cuts, 'w') if a.cuts else None
out = open(a.out, 'w'); trimmed = 0; removed = 0; dropped = 0; profile = {}
for k in order:
    s = seqs[k]; L = len(s); d = depth[k]; prof = []; run = 0
    if a.profile: prof = prof_in.get(k, [0] * L)[:L]
    else:
        for i in range(1, L + 1): run += d[i]; prof.append(run)
    med = statistics.median(prof) if prof else 0
    # representatives (2026-10-20): at --ratio 3 depth cut 3-9x family variation from them (zebra finch: -1.5 points
    # sensitivity) and --abs-depth dropped high-copy LTR families whole; the insertions worth cutting were 10x+
    # (gorilla: a 5' Alu at 19x, which uncut masked 294 Mb)
    rep = k in cons
    thr = max(a.min_depth, (a.consensus_ratio if rep else a.ratio) * max(med, 1))
    hi = [p >= thr or (not rep and a.abs_depth > 0 and p >= a.abs_depth) for p in prof]
    runs = []; i = 0
    while i < L:
        if hi[i]:
            j = i
            while j < L and hi[j]: j += 1
            if j - i >= a.min_len: runs.append((i, j))
            i = j
        else: i += 1
    total = sum(j - i for i, j in runs)
    if runs and a.drop_frac > 0 and total > a.drop_frac * L and not rep:
        dropped += 1
        if cuts: cuts.write(f'{k}\tdrop\n')
        continue
    if runs and total <= a.max_frac * L:
        keep = []; pos = 0
        for i, j in runs: keep.append(s[pos:i]); pos = j
        keep.append(s[pos:]); new = ''.join(keep); trimmed += 1; removed += total
        cut = True
    else: new = s; cut = False
    if len(new) < a.min_len:
        dropped += 1
        if cuts: cuts.write(f'{k}\tdrop\n')
        continue
    if cut and cuts: cuts.write(k + '\t' + ','.join(f'{i}-{j}' for i, j in runs) + '\n')
    out.write(f'>{k}\n{new}\n')
out.close()
print(f'{trimmed} records trimmed, {removed:,} bases removed, {dropped} dropped', file=sys.stderr)
