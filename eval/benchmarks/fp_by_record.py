#!/usr/bin/env python3
"""Where the false-positive masked bases come from: for a RepeatMasker run of a
library made by ltr-forage (records named seq:start..end_LTR or _INT), sum the
masked bases outside the standard mask per record kind and, with the
ltr-forage table, per evidence class of the element the record came from.

Usage: fp_by_record.py <standard.out> <test.out> [table.tsv]
"""
import sys, collections, bisect, csv, re

std_path, tst_path = sys.argv[1], sys.argv[2]
tsv = sys.argv[3] if len(sys.argv) > 3 else None

def load(path, cls='LTR'):
    iv = collections.defaultdict(list)
    for line in open(path):
        c = line.split()
        if len(c) < 15 or not c[0].isdigit(): continue
        if cls and not c[10].startswith(cls): continue
        iv[c[4]].append((int(c[5]), int(c[6])))
    out = {}
    for k, v in iv.items():
        v.sort(); m = []
        for s, e in v:
            if m and s <= m[-1][1] + 1: m[-1][1] = max(m[-1][1], e)
            else: m.append([s, e])
        out[k] = m
    return out
std = load(std_path)
starts = {k: [a for a, b in v] for k, v in std.items()}
def covered(chrom, s, e):
    v = std.get(chrom, []); i = bisect.bisect_left(starts.get(chrom, []), s) - 1; tot = 0
    i = max(i, 0)
    while i < len(v) and v[i][0] <= e:
        o = min(e, v[i][1]) - max(s, v[i][0]) + 1
        if o > 0: tot += o
        i += 1
    return tot

evidence = {}
if tsv:
    for r in csv.DictReader(open(tsv), delimiter='\t'):
        key = (r['seq'], int(r['elem_start']), int(r['elem_end']))
        ev = r['evidence']
        cls = 'domain' if ev.startswith('domain') else 'motif+tsd' if 'tgca' in ev else 'other'
        evidence[key] = cls
by_kind = collections.Counter(); by_ev = collections.Counter(); tot_kind = collections.Counter()
per_record = collections.Counter()
# element lookup: record coords lie inside the element
elems = collections.defaultdict(list)
for (seq, s, e), cls in evidence.items(): elems[seq].append((s, e, cls))
for v in elems.values(): v.sort()
estarts = {k: [a for a, b, c in v] for k, v in elems.items()}
def elem_class(seq, s, e):
    v = elems.get(seq, []); i = bisect.bisect_right(estarts.get(seq, []), s) - 1
    while i >= 0 and v[i][1] >= e:
        if v[i][0] <= s: return v[i][2]
        i -= 1
    return '?'
for line in open(tst_path):
    c = line.split()
    if len(c) < 15 or not c[0].isdigit(): continue
    chrom, s, e, name = c[4], int(c[5]), int(c[6]), c[9]
    fp = (e - s + 1) - covered(chrom, s, e)
    kind = 'INT' if name.endswith('_INT') else 'LTR' if name.endswith('_LTR') else 'other'
    by_kind[kind] += fp; tot_kind[kind] += e - s + 1; per_record[name] += fp
    m = re.match(r'(.+):(\d+)\.\.(\d+)_(LTR|INT)$', name)
    if m and evidence:
        by_ev[elem_class(m.group(1), int(m.group(2)), int(m.group(3)))] += fp
print("masked bases outside the standard, by record kind (and total masked by that kind):")
for k in by_kind: print(f"  {k}: {by_kind[k]:,} of {tot_kind[k]:,} ({100*by_kind[k]/max(1,tot_kind[k]):.1f}%)")
if evidence:
    print("by evidence class of the source element:")
    for k, v in by_ev.most_common(): print(f"  {k}: {v:,}")
print("top 15 records by false-positive bases:")
for n, v in per_record.most_common(15): print(f"  {n}: {v:,}")
