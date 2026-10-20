"""Cut left/right windows around each ltr-forage candidate for re-extension.
Window: outward by E = (LTR_MAX - ltr_len) + VIC, inward by min(E, half the internal span).
Writes chunk_NNN.q.fa (left windows), chunk_NNN.s.fa (right windows), and windows.tsv."""
import sys
pred, fai_path, genome, nchunks = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
LTR_MAX, VIC = 2000, 100
fai = {}; names = []
for l in open(fai_path):
    c = l.split(); fai[c[0]] = tuple(map(int, c[1:5])); names.append(c[0])
fh = open(genome, 'rb')
def fetch(chrom, s, e):  # 1-based inclusive, clamped
    ln, off, lb, lw = fai[chrom]; s = max(1, s); e = min(ln, e)
    out = b''; pos = s - 1
    while pos < e:
        line, col = divmod(pos, lb); fh.seek(off + line*lw + col); n = min(lb-col, e-pos)
        out += fh.read(n); pos += n
    return out.decode().upper(), s, e
cands = []
for line in open(pred):
    c = line.split()
    if len(c) < 11 or not c[0].isdigit(): continue
    cands.append((names[int(c[10])], int(c[3]), int(c[4]), int(c[6]), int(c[7]), c[9]))
per = -(-len(cands) // nchunks)
meta = open('windows.tsv', 'w')
files = {}
for i, (chrom, l1s, l1e, l2s, l2e, sim) in enumerate(cands):
    k = i // per
    if k not in files: files[k] = (open(f'chunk_{k:03d}.q.fa', 'w'), open(f'chunk_{k:03d}.s.fa', 'w'))
    ltr = max(l1e-l1s+1, l2e-l2s+1)
    E = (LTR_MAX - ltr) + VIC if ltr < LTR_MAX else VIC
    Ein = min(E, (l2s - l1e - 1) // 2)
    qseq, qs, qe = fetch(chrom, l1s - E, l1e + Ein)
    sseq, ss, se = fetch(chrom, l2s - Ein, l2e + E)
    files[k][0].write(f'>{i}\n{qseq}\n'); files[k][1].write(f'>{i}\n{sseq}\n')
    meta.write(f'{i}\t{chrom}\t{l1s}\t{l1e}\t{l2s}\t{l2e}\t{sim}\t{qs}\t{qe}\t{ss}\t{se}\t{k}\n')
for q, s in files.values(): q.close(); s.close()
print(f'{len(cands)} candidates, {len(files)} chunks', file=sys.stderr)
