import statistics, collections, itertools
comp = str.maketrans('ACGT', 'TGCA')
def load_counts(fn):
    d = {}
    for l in open(fn):
        c = l.split()
        if len(c) == 2: d[c[0]] = int(c[1])
    return d
ci, cl = load_counts('counts_int20.txt'), load_counts('counts_ltr20.txt')
def both(d, q): return d.get(q, 0) + d.get(q.translate(comp)[::-1], 0)
rows = [l.rstrip('\n').split('\t') for l in open('features20.tsv')][1:]
ki = [l.strip() for l in open('kmers_int20.txt')]; kl = [l.strip() for l in open('kmers_ltr20.txt')]
pi = pl = 0
feats = []
for r in rows:
    nki, nkl = int(r[17]), int(r[18])
    qi = ki[pi:pi+nki]; pi += nki; ql = kl[pl:pl+nkl]; pl += nkl
    ci_l = [both(ci, q) for q in qi]; cint = statistics.median(ci_l) if qi else 0; fint = sum(1 for x in ci_l if x >= 3) / len(ci_l) if ci_l else 0
    cltr = statistics.median([both(cl, q) for q in ql]) if ql else 0
    feats.append(dict(real=int(r[6]), ltr=int(r[7]), intl=int(r[8]), sim=float(r[9]), domain=int(r[10]), order=r[11],
                      tsd=int(r[13]), motif=int(r[14]), ppt=int(r[15]), pbs=int(r[16]), cint=cint, cltr=cltr, fint=fint))
real = [f for f in feats if f['real']]; other = [f for f in feats if not f['real']]
print(f"candidates {len(feats)}: real {len(real)}, other {len(other)}")
def rate(pred, name):
    tp = sum(1 for f in real if pred(f)); fp = sum(1 for f in other if pred(f))
    prec = tp / max(1, tp+fp)
    print(f"  {name:52s} real {100*tp/len(real):5.1f}%  other {100*fp/len(other):5.1f}%  precision {100*prec:5.1f}%  (keeps {tp+fp})")
print("== single signals (share of real / other candidates showing it)")
rate(lambda f: f['domain'], 'domain hit (any)')
rate(lambda f: f['order'] == 'LTR', 'domain hit, LTR order')
rate(lambda f: f['motif'], 'TG..CA at boundaries (+-2 bp)')
rate(lambda f: f['tsd'] >= 4, 'TSD >= 4 bp flanking')
rate(lambda f: f['tsd'] >= 5, 'TSD >= 5 bp flanking')
rate(lambda f: f['motif'] and f['tsd'] >= 4, 'motif and TSD >= 4')
rate(lambda f: f['ppt'] >= 10, 'PPT: purine run >= 10 within 40 bp of 3prime LTR')
rate(lambda f: f['pbs'], 'PBS-like TGG within 20 bp after 5prime LTR')
rate(lambda f: f['cint'] >= 3, 'internal region multi-copy (median 30-mer count >= 3)')
rate(lambda f: f['cint'] >= 10, 'internal region multi-copy (>= 10)')
rate(lambda f: f['cltr'] >= 3, 'LTR multi-copy (>= 3)')
rate(lambda f: 200 <= f['ltr'] <= 1500, 'LTR length 200-1500')
rate(lambda f: f['intl'] >= 1000, 'internal >= 1000 bp')
rate(lambda f: f['fint'] >= 0.25, 'internal: >=25% of 20-mers seen >=3 times')
rate(lambda f: f['fint'] >= 0.5, 'internal: >=50% of 20-mers seen >=3 times')
rate(lambda f: f['ppt'] >= 13, 'PPT: >=13/15 purines in a window within 30 bp of 3prime LTR')
rate(lambda f: f['ppt'] >= 14, 'PPT: >=14/15 purines')
rate(lambda f: f['domain'] or (f['motif'] and f['tsd'] >= 4) or (f['fint'] >= 0.25 and f['motif'] and f['ppt'] >= 13), 'domain or motif+TSD or (multi-copy and motif and PPT13)')
print("== combinations")
rate(lambda f: f['cint'] >= 3 and f['motif'] and f['tsd'] >= 4, 'multi-copy internal and motif and TSD')
rate(lambda f: f['cint'] >= 3 and (f['domain'] or (f['motif'] and f['tsd'] >= 4)), 'multi-copy internal and (domain or motif+TSD)')
rate(lambda f: f['domain'] or (f['motif'] and f['tsd'] >= 4), 'domain or motif+TSD')
rate(lambda f: f['cint'] >= 3 and f['ppt'] >= 10, 'multi-copy internal and PPT')
rate(lambda f: f['cint'] >= 3 and f['motif'] and f['ppt'] >= 10, 'multi-copy internal and motif and PPT')
rate(lambda f: f['domain'] or (f['cint'] >= 3 and f['motif'] and (f['tsd'] >= 4 or f['ppt'] >= 10)), 'domain or (multi-copy and motif and (TSD or PPT))')
print("== internal-region copy number distribution (median 30-mer count)")
bins = [(0,1),(1,2),(2,3),(3,10),(10,100),(100,10**9)]
for lo, hi in bins:
    r = sum(1 for f in real if lo <= f['cint'] < hi); o = sum(1 for f in other if lo <= f['cint'] < hi)
    print(f"  {lo:>4}-{hi if hi<10**9 else 'inf':<4} real {100*r/len(real):5.1f}%  other {100*o/len(other):5.1f}%")
print("== real elements without a domain hit: what else do they show?")
nd = [f for f in real if not f['domain']]
print(f"  n={len(nd)}: multi-copy internal {100*sum(1 for f in nd if f['cint']>=3)/len(nd):.0f}%, motif {100*sum(1 for f in nd if f['motif'])/len(nd):.0f}%, TSD>=4 {100*sum(1 for f in nd if f['tsd']>=4)/len(nd):.0f}%, PPT {100*sum(1 for f in nd if f['ppt']>=10)/len(nd):.0f}%")
