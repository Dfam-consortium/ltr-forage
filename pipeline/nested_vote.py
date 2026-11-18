# Cut nested LTR copies out of internal records, judged by a voting cluster.  For an internal
# record the voting cluster is its own element cluster when that has two or more members,
# otherwise the multi-element cluster whose internal records cover most of it.  A 50 bp bin is
# cut when an LTR record aligns across it as a whole copy (80% of the LTR record, at least 100 bp)
# and no internal record of the voting cluster covers the bin.  Runs shorter than 100 bp stay.
# usage: nested_vote.py <allvsall.tsv> <clusters.json> <records.fa> <out.fa>   (env TSV as in step4.py)
import sys, os, json, collections
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from records import element_lookup
ava, clusfile, recfa, outfa = sys.argv[1:5]
BIN=50; MIN_ID=80.0
clus=json.load(open(clusfile))
seqs={}; names={}; n=None
for l in open(recfa):
    if l.startswith('>'): full=l[1:].strip().split()[0]; n=full.split('#')[0]; names[n]=full; seqs[n]=[]
    else: seqs[n].append(l.strip())
seqs={k:''.join(v) for k,v in seqs.items()}
elem=element_lookup(clus)
size=collections.Counter(clus.values())
rec_clu={r:clus.get(elem(r)) for r in seqs}
int_bins=collections.defaultdict(lambda: collections.defaultdict(set))   # q -> cluster -> bins covered by that cluster's INTs
ltr_bins=collections.defaultdict(set)                                     # q -> bins covered by whole-copy LTR hits
for l in open(ava):
    q,s,pid,ln,qs,qe,ss,se,ql,sl,bits=l.split('\t'); q=q.split('#')[0]; s=s.split('#')[0]
    if q==s or not q.endswith('_INT') or float(pid)<MIN_ID: continue
    a,b=min(int(qs),int(qe)),max(int(qs),int(qe)); bins=range((a-1)//BIN,(b-1)//BIN+1)
    if s.endswith('_INT'):
        c=rec_clu.get(s)
        if c is not None and size[c]>=2 and c!=rec_clu.get(q): int_bins[q][c].update(bins)
        elif c is not None and c==rec_clu.get(q): int_bins[q][c].update(bins)
    else:
        sl=int(sl)
        if sl>=100 and abs(int(se)-int(ss))+1>=0.8*sl: ltr_bins[q].update(bins)
cuts={}; n_bp=0
for q,sq in seqs.items():
    if not q.endswith('_INT') or q not in ltr_bins: continue
    own=rec_clu.get(q); nb=(len(sq)-1)//BIN+1
    if own is not None and size[own]>=2: voter=int_bins[q].get(own,set())
    else:
        cands=int_bins[q]
        voter=max(cands.values(),key=len) if cands else set()
    if not voter: continue
    cut=[bi for bi in sorted(ltr_bins[q]) if bi not in voter]
    runs=[]
    for bi in cut:
        if runs and bi==runs[-1][1]+1: runs[-1][1]=bi
        else: runs.append([bi,bi])
    iv=[(r[0]*BIN+1,min(len(sq),(r[1]+1)*BIN)) for r in runs if (r[1]-r[0]+1)*BIN>=100]
    if iv: cuts[q]=iv
with open(outfa,'w') as out:
    for q,sq in seqs.items():
        iv=cuts.get(q,[]); pos=0; keep=[]
        for a,b in iv: keep.append(sq[pos:a-1]); pos=b
        keep.append(sq[pos:]); new=''.join(keep); n_bp+=len(sq)-len(new)
        if len(new)>=50: out.write(f">{names[q]}\n{new}\n")
print(f"records cut: {len(cuts)}, bases removed: {n_bp:,}")
