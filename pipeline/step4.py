# Two-phase element clustering. Phase 1: link elements whose internals cover each other (80% of both)
# and whose LTRs cover the shorter LTR. Phase 2: attach each remaining singleton element to the
# cluster whose members cover the singleton's own internal and LTR records; attachments never
# merge clusters.
import os, csv, collections, sys, json
int_hits, ltr_hits = sys.argv[1], sys.argv[2]
elem_recs={}; elems=[]
if 'TSV' not in os.environ: sys.exit('step4.py: set TSV to the ltr-forage table (cand.tsv)')
for r in csv.DictReader(open(os.environ['TSV']),delimiter='\t'):
    if r['keep']!='1': continue
    s=r['seq']; e=(s,int(r['elem_start']),int(r['elem_end'])); elems.append(e)
    ls,le,rs,re_=int(r['lltr_start']),int(r['lltr_end']),int(r['rltr_start']),int(r['rltr_end'])
    elem_recs[e]={'L':f"{s}:{ls}..{le}_LTR",'I':f"{s}:{le+1}..{rs-1}_INT",'R':f"{s}:{rs}..{re_}_LTR"}
rec2elem={n:(e,k) for e,d in elem_recs.items() for k,n in d.items()}
def union(iv):
    iv=sorted(iv); tot=0; cs=ce=None
    for s,e in iv:
        if cs is None: cs,ce=s,e
        elif s<=ce+1: ce=max(ce,e)
        else: tot+=ce-cs+1; cs,ce=s,e
    if cs is not None: tot+=ce-cs+1
    return tot
# directed coverage: cov[(a,b)][kind] = fraction of a's record covered by b's record, plus bits
cov=collections.defaultdict(dict)
for path,kind in ((int_hits,'I'),(ltr_hits,'L')):
    pair=collections.defaultdict(lambda: {'q':[], 's':[], 'ql':0, 'sl':0, 'bits':0.0})
    for l in open(path):
        q,s,pid,ln,qs,qe,ss,se,ql,sl,bits=l.split('\t'); q=q.split('#')[0]; s=s.split('#')[0]
        if q==s or float(pid)<80: continue
        if not (q.endswith('_INT')==(kind=='I') and s.endswith('_INT')==(kind=='I')): continue
        p=pair[(q,s)]; p['q'].append((min(int(qs),int(qe)),max(int(qs),int(qe)))); p['s'].append((min(int(ss),int(se)),max(int(ss),int(se)))); p['ql']=int(ql); p['sl']=int(sl); p['bits']+=float(bits)
    for (q,s),p in pair.items():
        if q not in rec2elem or s not in rec2elem: continue
        eq,es=rec2elem[q][0],rec2elem[s][0]
        if eq==es: continue
        cq=union(p['q'])/p['ql']; cs=union(p['s'])/p['sl']
        for a,b,ca,cb in ((eq,es,cq,cs),(es,eq,cs,cq)):
            d=cov[(a,b)]; old=d.get(kind,(0,0,0))
            if ca>old[0] or (ca==old[0] and p['bits']>old[2]): d[kind]=(ca,cb,p['bits'])
parent={e:e for e in elems}
def find(x):
    while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
    return x
for (a,b),d in cov.items():
    if a<b and 'I' in d and 'L' in d:
        ci,cj,_=d['I']; li,lj,_=d['L']
        if ci>=0.8 and cj>=0.8 and max(li,lj)>=0.8:
            ra,rb=find(a),find(b)
            if ra!=rb: parent[ra]=rb
def report(label, member):
    cores=collections.defaultdict(list)
    for e in elems: cores[member(e)].append(e)
    print(f"{label}: {len(cores)} clusters ({sum(1 for v in cores.values() if len(v)==1)} singletons, largest {max(len(v) for v in cores.values())})")
    return cores
cores=report('phase 1 (INT both 80%, LTR shorter 80%)', find)
# phase 2: attach singletons
attach={}
members={r:set(v) for r,v in cores.items() if len(v)>1}
cand=collections.defaultdict(lambda: collections.defaultdict(float))
for (a,b),d in cov.items():
    if len(cores[find(a)])!=1: continue
    rb=find(b)
    if len(cores[rb])==1: continue
    if 'I' in d and 'L' in d and d['I'][0]>=0.8 and d['L'][0]>=0.8:
        cand[a][rb]=max(cand[a][rb], d['I'][2]+d['L'][2])
for a,c in cand.items(): attach[a]=max(c.items(),key=lambda x:x[1])[0]
print(f"phase 2 attached {len(attach)} of {sum(1 for v in cores.values() if len(v)==1)} singletons")
cores2=report('phase 1 + phase 2 attachment', lambda e: attach.get(e, find(e)))
json.dump({f"{e[0]}:{e[1]}..{e[2]}": f"{r[0]}:{r[1]}..{r[2]}" for e in elems for r in [attach.get(e,find(e))]}, open(sys.argv[3],'w')) if len(sys.argv)>3 else None
