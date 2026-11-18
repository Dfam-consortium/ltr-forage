import sys, collections, bisect, re
rm, pred, defalse = sys.argv[1:4]
ann=[]
for line in open(rm):
    c=line.split()
    if len(c)<15 or not c[0].isdigit() or not c[10].startswith('LTR'): continue
    ann.append((int(c[5]),int(c[6]),c[8],c[9],float(c[1]),c[10]))
ann.sort()
# a bare trailing I names an internal only in Repbase-named annotations (ATCOPIA19I); where internals are
# marked -int (RepeatMasker's libraries) it is part of an LTR's name (MLT1I)
BARE_I = not any(a[3].endswith('-int') for a in ann)
def isint(n): return n.endswith(('-int','_I','-I','_int')) or (BARE_I and n.endswith('I'))
truth=[]
for i,(s1,e1,st1,nm1,dv1,cl1) in enumerate(ann):
    if isint(nm1) or not (100<=e1-s1+1<=1000): continue
    seen=False
    for j in range(i+1,len(ann)):
        s2,e2,st2,nm2,dv2,cl2=ann[j]
        if s2-s1>15000: break
        if nm2==nm1 and st2==st1:
            if seen and s2-s1>=1000 and 100<=e2-s2+1<=1000: truth.append((s1,e1,s2,e2,max(dv1,dv2),nm1))
            break
        if nm2!=nm1 and isint(nm2) and st2==st1: seen=True
def ovl(a1,a2,b1,b2):
    o=min(a2,b2)-max(a1,b1)+1
    return o/min(a2-a1+1,b2-b1+1) if o>0 else 0
# defalse records
recs=[]
for blk in open(defalse).read().split('\n\n'):
    m=re.match(r'\S+:(\d+)\.\.(\d+)\t(\w+)\tmotif:(\S+)\tTSD:(\S+)\t\S+\t\S+\tIN:\S+\t(\S+)\t(\S+)\t(\S+)\t(\S+)',blk)
    if not m: continue
    s,e,dec,motif,tsd,sim,strand,sf,order=m.groups(); s,e=int(s),int(e)
    ll=re.search(r'lLTR: (\S+)',blk).group(1); ov=int(re.search(r'TSD-LTR overlap: (\S+)',blk).group(1))
    noprot = strand=='?' and sf.lower()=='unknown' and order.upper()=='NA'
    known = re.search(r'TGCT|TACA|TACT|TGGA|TGGT|TATA|TGTA|TGCC|TGCA',motif,re.I)
    r = ('A:blastn-len' if ll=='NA' else 'C:notLTR' if order.upper()=='NOTLTR' else 'D:overlap' if noprot and ov>0 else 'E:motif' if noprot and not known else 'F:tsd5' if noprot and 'TGCA' not in motif.upper() and len(tsd)!=5 else 'H:flank' if dec=='false' else 'ok')
    recs.append((s,e,dec,motif,tsd,r,sf,order,strand,sim))
print(f"truth elements: {len(truth)}")
for line in open(pred):
    c=line.split()
    if len(c)<11 or not c[0].isdigit(): continue
    s,e,l1s,l1e,l2s,l2e,sim=int(c[0]),int(c[1]),int(c[3]),int(c[4]),int(c[6]),int(c[7]),c[9]
    hits=[t for t in truth if ovl(l1s,l1e,t[0],t[1])>=0.5 and ovl(l2s,l2e,t[2],t[3])>=0.5]
    if not hits: continue
    t=hits[0]
    d=[r for r in recs if ovl(s,e,r[0],r[1])>=0.5]
    d=d[0] if d else None
    print(f"GT {s}..{e} sim={sim} LTRlen={l1e-l1s+1}/{l2e-l2s+1}  truth={t[5]} div={t[4]}  ->  " + (f"{d[2]} [{d[5]}] motif={d[3]} TSD={d[4]} sf={d[6]} order={d[7]} strand={d[8]} id={d[9]}" if d else "no defalse record (TRF/N filtered)"))
