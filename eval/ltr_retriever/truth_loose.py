"""Loosened truth harness.  An element is LTR block / internal / LTR block on one strand:
LTR blocks are runs of non-internal LTR-class rows (any names) with gaps <= MERGE bp,
100..3000 bp long; the span between two blocks must be >= 50% covered by LTR-internal rows
on the same strand; 5'-5' distance 1..20 kb; the partner is the nearest qualifying block."""
import sys, collections, bisect
MERGE = 100
def isint(n): return n.endswith(('-int','_I','-I','_int','I'))
def build(rm):
    ann = collections.defaultdict(list)
    for line in open(rm):
        c = line.split()
        if len(c) < 15 or not c[0].isdigit() or not c[10].startswith('LTR'): continue
        ann[c[4]].append((int(c[5]), int(c[6]), c[8], c[9]))
    truth = collections.defaultdict(list)
    for chrom, lst in ann.items():
        lst.sort()
        # LTR blocks per strand
        for strand in ('+', 'C'):
            rows = [r for r in lst if r[2] == strand]
            blocks = []; ints = []
            for s, e, st, nm in rows:
                if isint(nm): ints.append((s, e)); continue
                if blocks and s - blocks[-1][1] <= MERGE: blocks[-1][1] = max(blocks[-1][1], e)
                else: blocks.append([s, e])
            blocks = [b for b in blocks if 100 <= b[1]-b[0]+1 <= 3000]
            istart = [i[0] for i in ints]
            def intcov(a, b):
                k = bisect.bisect_left(istart, a - 20000); cov = 0
                while k < len(ints) and ints[k][0] <= b:
                    o = min(b, ints[k][1]) - max(a, ints[k][0]) + 1
                    if o > 0: cov += o
                    k += 1
                return cov / max(1, b - a + 1)
            fwd = {}; back = {}
            for i, b1 in enumerate(blocks):
                for j in range(i+1, len(blocks)):
                    b2 = blocks[j]; d = b2[0] - b1[0]
                    if d > 20000: break
                    if d < 1000 or b2[0] - b1[1] - 1 < 100: continue
                    if intcov(b1[1]+1, b2[0]-1) >= 0.5:
                        fwd[i] = j
                        back[j] = i
                        break
            for i, j in fwd.items():
                if back.get(j) == i:
                    b1, b2 = blocks[i], blocks[j]
                    truth[chrom].append((b1[0], b1[1], b2[0], b2[1]))
    for v in truth.values(): v.sort()
    return truth
def scorer(truth):
    ts = {c: [t[0] for t in v] for c, v in truth.items()}
    def ovl(a1,a2,b1,b2):
        o = min(a2,b2)-max(a1,b1)+1; return o/min(a2-a1+1,b2-b1+1) if o > 0 else 0
    def hit(chrom, l1s, l1e, l2s, l2e):
        v = truth.get(chrom, []); k = bisect.bisect_left(ts.get(chrom, []), l1s-3000)
        while k < len(v) and v[k][0] <= l1e+3000:
            t = v[k]
            if ovl(l1s,l1e,t[0],t[1]) >= 0.5 and ovl(l2s,l2e,t[2],t[3]) >= 0.5: return t
            k += 1
        return None
    return hit
if __name__ == '__main__':
    rm, fai = sys.argv[1], sys.argv[2]
    names = [l.split()[0] for l in open(fai)]
    truth = build(rm); hit = scorer(truth)
    print(f"loose truth elements: {sum(len(v) for v in truth.values())}")
    for p in sys.argv[3:]:
        n = 0; h = 0; hits = set()
        for line in open(p):
            c = line.split()
            if len(c) < 11 or not c[0].isdigit(): continue
            n += 1; t = hit(names[int(c[10])], int(c[3]), int(c[4]), int(c[6]), int(c[7]))
            if t: h += 1; hits.add((names[int(c[10])],) + t)
        print(f"  {p}: {n} predictions, {h} hit a loose truth element ({100*h/max(1,n):.1f}%), {len(hits)} distinct elements")
