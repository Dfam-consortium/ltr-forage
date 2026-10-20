#!/usr/bin/env python3
"""Build a full-length LTR element truth set from a RepeatMasker .out (LTR - int - LTR, same family, same strand)
and score LTRharvest-style predictions against it."""
import sys, collections, bisect
rm, fai = sys.argv[1], sys.argv[2]
preds = sys.argv[3:]   # GT-format outputs
names = [l.split()[0] for l in open(fai)]
ann = collections.defaultdict(list)
with open(rm) as f:
    for line in f:
        c = line.split()
        if len(c) < 15 or not c[0].isdigit(): continue
        if not c[10].startswith('LTR'): continue
        ann[c[4]].append((int(c[5]), int(c[6]), c[8], c[9], float(c[1]), c[10]))
truth = []  # (chrom, l1s, l1e, l2s, l2e, maxdiv, name)
for chrom, lst in ann.items():
    lst.sort()
    n = len(lst)
    for i in range(n):
        s1, e1, st1, nm1, dv1, cl1 = lst[i]
        if '_INT' in nm1.upper() or nm1.endswith('-int') or nm1.endswith('_I') or nm1.endswith('-I') or nm1.endswith('_int'): continue
        if not (100 <= e1 - s1 + 1 <= 1000): continue
        seen_int = False
        for j in range(i + 1, n):
            s2, e2, st2, nm2, dv2, cl2 = lst[j]
            if s2 - s1 > 15000: break
            if nm2 == nm1 and st2 == st1:
                if seen_int and s2 - s1 >= 1000 and 100 <= e2 - s2 + 1 <= 1000:
                    truth.append((chrom, s1, e1, s2, e2, max(dv1, dv2), nm1))
                break  # nearest same-name copy only (whether or not it qualified)
            if nm2 != nm1 and ('_INT' in nm2.upper() or nm2.endswith('-int') or nm2.endswith('_I') or nm2.endswith('-I') or nm2.endswith('_int') or nm2.endswith('I')) and st2 == st1:
                seen_int = True
print(f"truth elements (LTR/int/LTR, same family+strand, LTR 100-1000, 5'-5' dist 1000-15000): {len(truth)}")
bins = [(0,5),(5,10),(10,15),(15,20),(20,100)]
def binof(d): 
    for lo,hi in bins:
        if lo <= d < hi: return (lo,hi)
    return bins[-1]
tot = collections.Counter(binof(t[5]) for t in truth)
by_chrom = collections.defaultdict(list)
for t in truth: by_chrom[t[0]].append(t)
for v in by_chrom.values(): v.sort(key=lambda t: t[1])
starts = {c: [t[1] for t in v] for c, v in by_chrom.items()}

def ovl(a1, a2, b1, b2):
    o = min(a2, b2) - max(a1, b1) + 1
    return o / min(a2 - a1 + 1, b2 - b1 + 1) if o > 0 else 0.0

def load(path):
    out = []
    for line in open(path):
        c = line.split()
        if len(c) >= 11 and c[0].isdigit():
            out.append((names[int(c[10])], int(c[3]), int(c[4]), int(c[6]), int(c[7]), float(c[9])))
        elif len(c) == 8 and c[0].isdigit():   # a_only/b_only dumps: seq s e l1s l1e l2s l2e sim
            out.append((names[int(c[0])], int(c[3]), int(c[4]), int(c[5]), int(c[6]), float(c[7])))
    return out

def score(pred, label):
    hit = set(); supported = 0
    for chrom, l1s, l1e, l2s, l2e, sim in pred:
        v = by_chrom.get(chrom, [])
        k = bisect.bisect_left(starts.get(chrom, []), l1s - 2000)
        found = False
        while k < len(v) and v[k][1] <= l1e + 2000:
            t = v[k]
            if ovl(l1s, l1e, t[1], t[2]) >= 0.5 and ovl(l2s, l2e, t[3], t[4]) >= 0.5:
                hit.add((chrom, t[1], t[3])); found = True
            k += 1
        # weaker support: both predicted LTRs overlap some LTR-class annotation (any family)
        if not found:
            pass
        supported += found
    got = collections.Counter(binof(t[5]) for t in truth if (t[0], t[1], t[3]) in hit)
    print(f"\n{label}: {len(pred)} predictions; {supported} ({100*supported/max(1,len(pred)):.1f}%) match a truth element (both LTRs >=50% overlap)")
    print("  recall by max RM divergence of the two LTRs:")
    for b in bins:
        print(f"    {b[0]:>2}-{b[1]:<3}%: {got[b]:>6} / {tot[b]:<6} = {100*got[b]/max(1,tot[b]):5.1f}%")
    return hit

hits = {}
for p in preds:
    hits[p] = score(load(p), p)
if len(preds) >= 2:
    a, b = hits[preds[0]], hits[preds[1]]
    print(f"\ntruth hit by {preds[0]} only: {len(a-b)}; by {preds[1]} only: {len(b-a)}; both: {len(a&b)}")
