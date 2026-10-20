"""Alignment-intrinsic end trimming.  For the best plus-strand HSP of each candidate, walk in
from each end of the aligned strings and trim until the terminal window of W aligned columns
has at least M matches and no gap; the trimmed columns map back to query/subject offsets.
Writes GT-format predictions with the trimmed boundaries."""
import sys
import os
ROOT = os.environ.get('LTR_WORKSPACE') or os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..'))
hits, windows, out = sys.argv[1], sys.argv[2], sys.argv[3]
W, M = int(sys.argv[4]), int(sys.argv[5])
names=[l.split()[0] for l in open(ROOT+'/data/genomes/hg38/hg38.fa.fai')]; cidx={n:i for i,n in enumerate(names)}
best = {}
for line in open(hits):
    c = line.rstrip('\n').split('\t')
    if len(c) < 12 or c[5] != 'plus': continue
    sc, i = int(c[0]), int(c[1])
    if i not in best or sc > best[i][0]:
        best[i] = (sc, int(c[2]), int(c[3]), int(c[7]), int(c[8]), c[10], c[11])
def trim(qs, ss):
    """return (q_trim_left, s_trim_left, q_trim_right, s_trim_right): bases to drop from each end"""
    n = len(qs)
    def ok(a, b):  # columns a..b-1 form a clean window
        seg_q, seg_s = qs[a:b], ss[a:b]
        if '-' in seg_q or '-' in seg_s: return False
        return sum(1 for x, y in zip(seg_q, seg_s) if x == y) >= M
    left = 0
    while left + W <= n and not ok(left, left + W): left += 1
    right = n
    while right - W >= left and not ok(right - W, right): right -= 1
    if left + W > n or right - W < left: return None
    ql = sum(1 for x in qs[:left] if x != '-'); sl = sum(1 for x in ss[:left] if x != '-')
    qr = sum(1 for x in qs[right:] if x != '-'); sr = sum(1 for x in ss[right:] if x != '-')
    return ql, sl, qr, sr
n = moved = 0
with open(out, 'w') as o:
    for line in open(windows):
        i, chrom, l1s, l1e, l2s, l2e, sim, qs0, qe0, ss0, se0, k = line.split()
        i = int(i); qs0, ss0 = int(qs0), int(ss0)
        if i not in best:
            L1, E1, L2, E2 = map(int, (l1s, l1e, l2s, l2e))
        else:
            sc, qstart, qend, sstart, send, qseq, sseq = best[i]
            t = trim(qseq, sseq)
            if t: ql, sl, qr, sr = t; moved += (ql + sl + qr + sr) > 0
            else: ql = sl = qr = sr = 0
            L1, E1 = qs0 + qstart - 1 + ql, qs0 + qend - 1 - qr
            L2, E2 = ss0 + sstart - 1 + sl, ss0 + send - 1 - sr
        o.write(f'{L1} {E2} {E2-L1+1} {L1} {E1} {E1-L1+1} {L2} {E2} {E2-L2+1} {sim} {cidx[chrom]}\n'); n += 1
print(f'{n} written, {moved} trimmed', file=sys.stderr)
