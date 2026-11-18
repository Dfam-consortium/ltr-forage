import sys, collections, bisect, re
rm, fai, pred, defalse = sys.argv[1:5]
names = [l.split()[0] for l in open(fai)]
ann = collections.defaultdict(list)
for line in open(rm):
    c = line.split()
    if len(c) < 15 or not c[0].isdigit() or not c[10].startswith('LTR'): continue
    ann[c[4]].append((int(c[5]), int(c[6]), c[8], c[9], float(c[1])))
# a bare trailing I names an internal only in Repbase-named annotations (ATCOPIA19I); where internals are
# marked -int (RepeatMasker's libraries) it is part of an LTR's name (MLT1I)
BARE_I = not any(a[3].endswith('-int') for lst in ann.values() for a in lst)
def isint(n): return n.endswith(('-int','_I','-I','_int')) or (BARE_I and n.endswith('I'))
truth = collections.defaultdict(list)
for chrom, lst in ann.items():
    lst.sort()
    for i, (s1,e1,st1,nm1,dv1) in enumerate(lst):
        if isint(nm1) or not (100 <= e1-s1+1 <= 1000): continue
        seen = False
        for j in range(i+1, len(lst)):
            s2,e2,st2,nm2,dv2 = lst[j]
            if s2-s1 > 15000: break
            if nm2 == nm1 and st2 == st1:
                if seen and s2-s1 >= 1000 and 100 <= e2-s2+1 <= 1000: truth[chrom].append((s1,e1,s2,e2,max(dv1,dv2),nm1))
                break
            if nm2 != nm1 and isint(nm2) and st2 == st1: seen = True
tstarts = {c: [t[0] for t in v] for c, v in truth.items()}
def ovl(a1,a2,b1,b2):
    o = min(a2,b2)-max(a1,b1)+1
    return o/min(a2-a1+1,b2-b1+1) if o > 0 else 0
# defalse records keyed by chrom
recs = collections.defaultdict(list)
for blk in open(defalse).read().split('\n\n'):
    m = re.match(r'seq(\d+):(\d+)\.\.(\d+)\t(\w+)\tmotif:(\S+)\tTSD:(\S+)\t\S+\t\S+\tIN:\S+\t(\S+)\t(\S+)\t(\S+)\t(\S+)', blk)
    if not m: continue
    sn,s,e,dec,motif,tsd,sim,strand,sf,order = m.groups(); s,e = int(s),int(e)
    ll = re.search(r'lLTR: (\S+)', blk).group(1); ov = int(re.search(r'TSD-LTR overlap: (\S+)', blk).group(1))
    noprot = strand=='?' and sf.lower()=='unknown' and order.upper()=='NA'
    known = re.search(r'TGCT|TACA|TACT|TGGA|TGGT|TATA|TGTA|TGCC|TGCA', motif, re.I) is not None
    r = ('A: self-blastn pair / LTR-length mismatch' if ll=='NA' else 'C: HMM notLTR' if order.upper()=='NOTLTR' else 'D: no protein evidence, TSD-LTR overlap' if noprot and ov>0 else 'E: no protein evidence, motif unknown' if noprot and not known else 'F: no protein evidence, nonTGCA, TSD!=5' if noprot and 'TGCA' not in motif.upper() and len(tsd)!=5 else 'H: flank windows aligned' if dec=='false' else 'ok')
    rescued = (not noprot) and (ov>0 or not known or ('TGCA' not in motif.upper() and len(tsd)!=5))
    recs[names[int(sn)-1]].append((min(s,e), max(s,e), dec, r, noprot, rescued))
rstarts = {c: [x[0] for x in sorted(v)] for c, v in recs.items()}
for c in recs: recs[c].sort()
def find_rec(chrom, s, e):
    v = recs.get(chrom, []); k = bisect.bisect_left(rstarts.get(chrom, []), s-20000)
    while k < len(v) and v[k][0] <= e:
        if ovl(s,e,v[k][0],v[k][1]) >= 0.5: return v[k]
        k += 1
    return None
by_reason = collections.Counter(); offs = collections.defaultdict(list); kept_truth = collections.Counter()
n_truth = 0
for line in open(pred):
    c = line.split()
    if len(c) < 11 or not c[0].isdigit(): continue
    chrom = names[int(c[10])]; s,e,l1s,l1e,l2s,l2e = map(int, (c[0],c[1],c[3],c[4],c[6],c[7]))
    v = truth.get(chrom, []); k = bisect.bisect_left(tstarts.get(chrom, []), l1s-2000); hit = None
    while k < len(v) and v[k][0] <= l1e+2000:
        t = v[k]
        if ovl(l1s,l1e,t[0],t[1]) >= 0.5 and ovl(l2s,l2e,t[2],t[3]) >= 0.5: hit = t; break
        k += 1
    if not hit: continue
    n_truth += 1
    r = find_rec(chrom, s, e)
    key = f"{r[2]}\t{r[3]}" if r else "no defalse record (TRF/N stage)"
    by_reason[key] += 1
    off = abs(l1s - hit[0]) + abs(l1e - hit[1]) + abs(l2s - hit[2]) + abs(l2e - hit[3])
    offs[r[2] if r else 'trf'].append(off)
print(f"GT candidates hitting a truth element: {n_truth}")
print("== decision / reason for those real elements")
for k, v in sorted(by_reason.items(), key=lambda x: -x[1]): print(f"{v:6d}  {k}")
print("== sum of |LTR boundary offsets| vs RepeatMasker (bp), by decision: median / mean / fraction <= 20 bp")
import statistics
for d, L in offs.items(): print(f"  {d:10s} n={len(L):5d} median={statistics.median(L):7.1f} mean={statistics.mean(L):8.1f} exact(<=20bp)={100*sum(1 for x in L if x<=20)/len(L):5.1f}%")
# among kept records, truth-hit rate for protein-rescued vs others
kept = [(c, x) for c, v in recs.items() for x in v if x[2] != 'false']
hitset = set()
for line in open(pred):
    c = line.split()
    if len(c) < 11 or not c[0].isdigit(): continue
    chrom = names[int(c[10])]; s,e,l1s,l1e,l2s,l2e = map(int, (c[0],c[1],c[3],c[4],c[6],c[7]))
    v = truth.get(chrom, []); k = bisect.bisect_left(tstarts.get(chrom, []), l1s-2000)
    while k < len(v) and v[k][0] <= l1e+2000:
        t = v[k]
        if ovl(l1s,l1e,t[0],t[1]) >= 0.5 and ovl(l2s,l2e,t[2],t[3]) >= 0.5: hitset.add((chrom, s, e)); break
        k += 1
tab = collections.Counter()
for chrom, x in kept:
    is_hit = any(ovl(x[0],x[1],hs,he) >= 0.5 for (hc,hs,he) in [h for h in hitset if h[0]==chrom and abs(h[1]-x[0])<20000])
    tab[(x[2], 'rescued by protein evidence' if x[5] else ('protein evidence' if not x[4] else 'no protein evidence'), is_hit)] += 1
print("== kept elements: truth-hit rate by evidence class")
for (dec, cls) in sorted({(k[0],k[1]) for k in tab}):
    h, n = tab[(dec,cls,True)], tab[(dec,cls,True)] + tab[(dec,cls,False)]
    print(f"  {dec:10s} {cls:28s} n={n:4d} truth-hit={h:4d} ({100*h/max(1,n):5.1f}%)")
