# One exemplar element per cluster: the member with the longest internal record after cleanup,
# written as its three records (names kept, so the RepeatMasker output stays traceable).
# usage: exemplar_lib.py <clusters.json> <records.fa> <out.fa>   (env TSV as in step4.py)
import os, sys, json, csv, collections
clusfile, recfa, outfa = sys.argv[1:4]
clus=json.load(open(clusfile))
seqs={}; names={}; n=None
for l in open(recfa):
    if l.startswith('>'): full=l[1:].strip().split()[0]; n=full.split('#')[0]; names[n]=full; seqs[n]=[]
    else: seqs[n].append(l.strip())
seqs={k:''.join(v) for k,v in seqs.items()}
elem_recs={}
for r in csv.DictReader(open(os.environ.get('TSV','rice.m5.tsv')),delimiter='\t'):
    if r['keep']!='1': continue
    s=r['seq']; e=f"{s}:{r['elem_start']}..{r['elem_end']}"
    ls,le,rs,re_=int(r['lltr_start']),int(r['lltr_end']),int(r['rltr_start']),int(r['rltr_end'])
    elem_recs[e]=[f"{s}:{ls}..{le}_LTR",f"{s}:{le+1}..{rs-1}_INT",f"{s}:{rs}..{re_}_LTR"]
members=collections.defaultdict(list)
for e,c in clus.items(): members[c].append(e)
n_out=0
with open(outfa,'w') as out:
    for c,v in members.items():
        best=max(v,key=lambda e: len(seqs.get(elem_recs[e][1],'')))
        for rec in elem_recs[best]:
            if rec in seqs and len(seqs[rec])>=50:
                out.write(f">{names[rec]}\n{seqs[rec]}\n"); n_out+=1
print(f"{len(members)} clusters, {n_out} records written to {outfa}")
