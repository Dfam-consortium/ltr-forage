"""Per-candidate features on the re-extended (HSP) boundaries, plus truth label and
LTR_retriever's domain call.  Also writes k-mer query files for sufr count."""
import sys, re, collections, bisect
import os
ROOT = os.environ.get('LTR_WORKSPACE') or os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..'))
rm, fai_path, genome = ROOT+'/data/annotations/hg38/hg38.out', ROOT+'/data/genomes/hg38/hg38.fa.fai', ROOT+'/data/genomes/hg38/hg38.fa'
pred, defalse = '../reext/reext_hsp.gt', '../lr_reext/seq.fa.defalse'
K, NK = int(os.environ.get("K", 30)), 8
fai = {}; names = []
for l in open(fai_path):
    c = l.split(); fai[c[0]] = tuple(map(int, c[1:5])); names.append(c[0])
fh = open(genome, 'rb')
def fetch(chrom, s, e):
    ln, off, lb, lw = fai[chrom]; s = max(1, s); e = min(ln, e)
    out = b''; pos = s - 1
    while pos < e:
        line, col = divmod(pos, lb); fh.seek(off + line*lw + col); n = min(lb-col, e-pos)
        out += fh.read(n); pos += n
    return out.decode().upper()
# truth
ann = collections.defaultdict(list)
for line in open(rm):
    c = line.split()
    if len(c) < 15 or not c[0].isdigit() or not c[10].startswith('LTR'): continue
    ann[c[4]].append((int(c[5]), int(c[6]), c[8], c[9]))
# a bare trailing I names an internal only in Repbase-named annotations (ATCOPIA19I); where internals are
# marked -int (RepeatMasker's libraries) it is part of an LTR's name (MLT1I)
BARE_I = not any(a[3].endswith('-int') for lst in ann.values() for a in lst)
def isint(n): return n.endswith(('-int','_I','-I','_int')) or (BARE_I and n.endswith('I'))
truth = collections.defaultdict(list)
for chrom, lst in ann.items():
    lst.sort()
    for i, (s1,e1,st1,nm1) in enumerate(lst):
        if isint(nm1) or not (100 <= e1-s1+1 <= 1000): continue
        seen = False
        for j in range(i+1, len(lst)):
            s2,e2,st2,nm2 = lst[j]
            if s2-s1 > 15000: break
            if nm2 == nm1 and st2 == st1:
                if seen and s2-s1 >= 1000 and 100 <= e2-s2+1 <= 1000: truth[chrom].append((s1,e1,s2,e2))
                break
            if nm2 != nm1 and isint(nm2) and st2 == st1: seen = True
ts = {c: [t[0] for t in v] for c, v in truth.items()}
def ovl(a1,a2,b1,b2):
    o = min(a2,b2)-max(a1,b1)+1; return o/min(a2-a1+1,b2-b1+1) if o > 0 else 0
# LTR_retriever domain evidence per locus (from the re-extended run)
dom = collections.defaultdict(list)
for blk in open(defalse).read().split('\n\n'):
    m = re.match(r'seq(\d+):(\d+)\.\.(\d+)\t(\w+)\tmotif:(\S+)\tTSD:(\S+)\t\S+\t\S+\tIN:\S+\t(\S+)\t(\S+)\t(\S+)\t(\S+)', blk)
    if not m: continue
    sn,s,e,dec,motif,tsd,sim,strand,sf,order = m.groups(); s,e = int(s),int(e)
    dom[names[int(sn)-1]].append((min(s,e), max(s,e), order.upper() != 'NA' or sf.lower() != 'unknown' or strand != '?', order.upper(), dec))
for c in dom: dom[c].sort()
ds = {c: [x[0] for x in v] for c, v in dom.items()}
def domain(chrom, s, e):
    v = dom.get(chrom, []); k = bisect.bisect_left(ds.get(chrom, []), s-20000)
    while k < len(v) and v[k][0] <= e:
        if ovl(s,e,v[k][0],v[k][1]) >= 0.5: return v[k][2], v[k][3], v[k][4]
        k += 1
    return None, 'none', 'none'
def tsd_len(chrom, L, R):
    """Longest exact direct repeat immediately flanking the element (allowing 0-2 bp slack at each end)."""
    left = fetch(chrom, L-22, L+1); right = fetch(chrom, R-2, R+21)   # left[20]=L-2 .. ; right[2]=R
    best = 0
    for dl in range(0, 3):          # element start slack: TSD ends at L-1-dl
        for dr in range(0, 3):
            for k in range(20, 3, -1):
                a = left[22-dl-k:22-dl]; b = right[3+dr:3+dr+k]
                if len(a) == k and a == b: best = max(best, k); break
    return best
def motif_ok(chrom, L, R):
    s = fetch(chrom, L-2, L+3) ; e = fetch(chrom, R-4, R+1)   # allow +-2 slack
    return any(s[i:i+2] == 'TG' for i in range(0, 5)) and any(e[i:i+2] == 'CA' for i in range(0, 5))
def ppt(chrom, l2s):
    """Max purine fraction (as percent) over 15 bp windows in the 30 bp before the 3' LTR."""
    w = fetch(chrom, l2s-30, l2s-1); best = 0
    for i in range(0, len(w)-14):
        s = w[i:i+15]; best = max(best, sum(1 for ch in s if ch in 'AG'))
    return best   # out of 15
def pbs(chrom, l1e):
    w = fetch(chrom, l1e+1, l1e+20)
    return 'TGG' in w
def kmers(seq, n):
    seq = seq.replace('N', '')
    if len(seq) < K: return []
    step = max(1, (len(seq) - K) // max(1, n-1))
    return [seq[i:i+K] for i in range(0, len(seq)-K+1, step)][:n]
out = open(f'features{K}.tsv', 'w'); qint = open(f'kmers_int{K}.txt', 'w'); qltr = open(f'kmers_ltr{K}.txt', 'w')
out.write('idx\tchrom\tl1s\tl1e\tl2s\tl2e\treal\tltr_len\tint_len\tsim\tdomain\torder\tdecision\ttsd\tmotif\tppt\tpbs\tnk_int\tnk_ltr\n')
for i, line in enumerate(open(pred)):
    c = line.split(); chrom = names[int(c[10])]
    l1s,l1e,l2s,l2e = map(int, (c[3],c[4],c[6],c[7])); sim = c[9]
    v = truth.get(chrom, []); k = bisect.bisect_left(ts.get(chrom, []), l1s-2000); real = 0
    while k < len(v) and v[k][0] <= l1e+2000:
        t = v[k]
        if ovl(l1s,l1e,t[0],t[1]) >= 0.5 and ovl(l2s,l2e,t[2],t[3]) >= 0.5: real = 1; break
        k += 1
    d, order, dec = domain(chrom, l1s, l2e)
    ltr = max(l1e-l1s+1, l2e-l2s+1); intl = l2s - l1e - 1
    ki = kmers(fetch(chrom, l1e+1, l2s-1), NK) if intl >= K else []
    kl = kmers(fetch(chrom, l1s, l1e), NK)
    for q in ki: qint.write(f'{q}\n')
    for q in kl: qltr.write(f'{q}\n')
    out.write(f'{i}\t{chrom}\t{l1s}\t{l1e}\t{l2s}\t{l2e}\t{real}\t{ltr}\t{intl}\t{sim}\t{int(bool(d))}\t{order}\t{dec}\t{tsd_len(chrom,l1s,l2e)}\t{int(motif_ok(chrom,l1s,l2e))}\t{ppt(chrom,l2s)}\t{int(pbs(chrom,l1e))}\t{len(ki)}\t{len(kl)}\n')
print('done', file=sys.stderr)
