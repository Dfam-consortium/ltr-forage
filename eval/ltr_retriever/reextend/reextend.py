"""Map the best pairwise rmblastn HSP back to genome coordinates, optionally trim to
the best TSD/motif placement, and write GT-format predictions.
Variants: hsp  = HSP ends as LTR boundaries
          tsd  = HSP ends, outer boundaries moved to the best TSD/motif within +-VIC
          keep = original ltr-forage boundaries (control)"""
import sys, collections
hits, windows, fai_path, genome, variant, out = sys.argv[1:7]
import os
VIC, MIN_TSD, MAX_TSD = int(os.environ.get("VIC", 60)), 4, 20
RANK = os.environ.get("RANK", "tier")
fai = {}; names = []
for l in open(fai_path):
    c = l.split(); fai[c[0]] = tuple(map(int, c[1:5])); names.append(c[0])
idx = {n: i for i, n in enumerate(names)}
fh = open(genome, 'rb')
def fetch(chrom, s, e):
    ln, off, lb, lw = fai[chrom]; s = max(1, s); e = min(ln, e)
    out = b''; pos = s - 1
    while pos < e:
        line, col = divmod(pos, lb); fh.seek(off + line*lw + col); n = min(lb-col, e-pos)
        out += fh.read(n); pos += n
    return out.decode().upper()
best = {}
for line in open(hits):
    c = line.split()
    if c[6] != 'plus': continue
    score, i = int(c[0]), int(c[2])
    if i not in best or score > best[i][0]:
        best[i] = (score, int(c[3]), int(c[4]), int(c[8]), int(c[9]), float(c[1]))
KNOWN = {'TGCT','TACA','TACT','TGGA','TATA','TGTA','TGCA'}
def tier(motif, tl):
    if motif == 'TGCA' and tl == 5: return 1
    if motif in KNOWN and tl == 5: return 2
    if motif == 'TGCA' and tl < 5: return 3
    if tl == 5: return 4
    if tl < 5: return 5
    return 6
def trim(chrom, L, R):
    """RANK=motif: move the outer ends to the nearest placement (within VIC) where the LTR
    starts with TG and ends with CA, minimal total deviation; else keep.
    RANK=tsdmotif: first try TSD >= 5 bp with a TG..CA motif within VIC (min deviation),
    then fall back to the motif-only search; else keep."""
    lo = max(1, L - VIC - MAX_TSD); left = fetch(chrom, lo, L + VIC + 1)
    rlo = max(1, R - VIC - 1);      right = fetch(chrom, rlo, R + VIC + MAX_TSD)
    def motif_at(newL, newR):
        p = newL - lo; q = newR - rlo
        if p < 0 or p + 2 > len(left) or q - 1 < 0 or q + 1 > len(right): return ''
        return left[p:p+2] + right[q-1:q+1]
    bestc = None
    if RANK == 'tsdmotif':
        for k in range(MAX_TSD, 4, -1):
            kmers = collections.defaultdict(list)
            for p in range(0, len(left) - k + 1): kmers[left[p:p+k]].append(p)
            for q in range(0, len(right) - k + 1):
                s = right[q:q+k]
                if s not in kmers: continue
                newR = rlo + q - 1
                if abs(newR - R) > VIC: continue
                for p in kmers[s]:
                    newL = lo + p + k
                    if abs(newL - L) > VIC or motif_at(newL, newR) != 'TGCA': continue
                    dev = abs(newL - L) + abs(newR - R)
                    if bestc is None or dev < bestc[0]: bestc = (dev, newL, newR)
        if bestc: return bestc[1], bestc[2]
    for dl in range(-VIC, VIC + 1):
        for dr in range(-VIC, VIC + 1):
            newL, newR = L + dl, R + dr
            if motif_at(newL, newR) == 'TGCA':
                dev = abs(dl) + abs(dr)
                if bestc is None or dev < bestc[0]: bestc = (dev, newL, newR)
    return (bestc[1], bestc[2]) if bestc else (L, R)
n = 0; nohit = 0; trimmed = 0
with open(out, 'w') as o:
    for line in open(windows):
        i, chrom, l1s, l1e, l2s, l2e, sim, qs, qe, ss, se, k = line.split()
        i = int(i); l1s, l1e, l2s, l2e, qs, ss = map(int, (l1s, l1e, l2s, l2e, qs, ss))
        if variant == 'keep' or i not in best:
            nohit += (i not in best); L1, E1, L2, E2 = l1s, l1e, l2s, l2e
        else:
            sc, qstart, qend, sstart, send, psub = best[i]
            L1, E1 = qs + qstart - 1, qs + qend - 1
            L2, E2 = ss + sstart - 1, ss + send - 1
            if variant == 'tsd':
                nL, nR = trim(chrom, L1, E2)
                trimmed += (nL, nR) != (L1, E2); L1, E2 = nL, nR
            sim = f'{100 - psub:.2f}'
        o.write(f'{L1} {E2} {E2-L1+1} {L1} {E1} {E1-L1+1} {L2} {E2} {E2-L2+1} {sim} {idx[chrom]}\n'); n += 1
print(f'{variant}: {n} written, {nohit} without an HSP, {trimmed} moved by TSD trim', file=sys.stderr)
