import sys, collections, bisect, re, statistics
rm, fai, pred, defalse = sys.argv[1:5]
names = [l.split()[0] for l in open(fai)]
ann = collections.defaultdict(list)
for line in open(rm):
    c = line.split()
    if len(c) < 15 or not c[0].isdigit() or not c[10].startswith('LTR'): continue
    ann[c[4]].append((int(c[5]), int(c[6]), c[8], c[9], float(c[1])))
def isint(n): return n.endswith(('-int','_I','-I','_int','I'))
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
                if seen and s2-s1 >= 1000 and 100 <= e2-s2+1 <= 1000: truth[chrom].append((s1,e1,s2,e2))
                break
            if nm2 != nm1 and isint(nm2) and st2 == st1: seen = True
tstarts = {c: [t[0] for t in v] for c, v in truth.items()}
def ovl(a1,a2,b1,b2):
    o = min(a2,b2)-max(a1,b1)+1
    return o/min(a2-a1+1,b2-b1+1) if o > 0 else 0
dec = {}
if defalse != '-':
    for blk in open(defalse).read().split('\n\n'):
        m = re.match(r'seq(\d+):(\d+)\.\.(\d+)\t(\w+)\t', blk)
        if m: dec[(names[int(m.group(1))-1], min(int(m.group(2)),int(m.group(3))), max(int(m.group(2)),int(m.group(3))))] = m.group(4)
    dk = collections.defaultdict(list)
    for (c,s,e),d in dec.items(): dk[c].append((s,e,d))
    for c in dk: dk[c].sort()
    dks = {c: [x[0] for x in v] for c,v in dk.items()}
def decision(chrom, s, e):
    if defalse == '-': return 'all'
    v = dk.get(chrom, []); k = bisect.bisect_left(dks.get(chrom, []), s-20000)
    while k < len(v) and v[k][0] <= e:
        if ovl(s,e,v[k][0],v[k][1]) >= 0.5: return v[k][2]
        k += 1
    return 'none'
sides = collections.defaultdict(lambda: collections.defaultdict(list)); n_hit = 0; n_pred = 0; exact = collections.Counter()
for line in open(pred):
    c = line.split()
    if len(c) < 11 or not c[0].isdigit(): continue
    n_pred += 1
    chrom = names[int(c[10])]; s,e,l1s,l1e,l2s,l2e = map(int, (c[0],c[1],c[3],c[4],c[6],c[7]))
    v = truth.get(chrom, []); k = bisect.bisect_left(tstarts.get(chrom, []), l1s-2000)
    while k < len(v) and v[k][0] <= l1e+2000:
        t = v[k]
        if ovl(l1s,l1e,t[0],t[1]) >= 0.5 and ovl(l2s,l2e,t[2],t[3]) >= 0.5:
            n_hit += 1; d = decision(chrom, s, e)
            offs = (l1s-t[0], l1e-t[1], l2s-t[2], l2e-t[3])
            for name, o in zip(('outer-left','inner-left','inner-right','outer-right'), offs): sides[d][name].append(o)
            exact[d] += sum(abs(o) for o in offs) <= 20
            break
        k += 1
print(f"{pred}: {n_pred} predictions, {n_hit} truth hits, exact (sum offsets <= 20 bp): {sum(exact.values())}")
for d in sides:
    n = len(sides[d]['outer-left'])
    print(f"  decision={d} n={n} exact={exact[d]}")
    for name in ('outer-left','inner-left','inner-right','outer-right'):
        L = sides[d][name]; off = [x for x in L if abs(x) > 20]
        short = sum(1 for x in off if (x > 0) == (name in ('outer-left','inner-right')))  # inside the true LTR
        print(f"    {name:12s} off by >20 bp: {100*len(off)/n:5.1f}%  (of those, {100*short/max(1,len(off)):5.1f}% short, i.e. inside the RM LTR); median |off| {statistics.median([abs(x) for x in L]):.0f}")
