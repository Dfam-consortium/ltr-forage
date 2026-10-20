"""What the motif route adds at each TSD length: counts of real and other candidates
without an LTR-order domain hit.  Usage: tsd_motif.py <rm.out> <table.tsv>"""
import csv, sys, collections, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); import truth_loose
def isint(n): return n.endswith(('-int','_I','-I','_int','I')) or '_INT' in n
truth_loose.isint=isint
ns={'collections':collections,'isint':isint}
src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'score_tsv.py')).read(); exec(src[src.index('def strict'):src.index('hit_s = ')], ns)
rm, tsv = sys.argv[1], sys.argv[2]
hs=truth_loose.scorer(ns['strict'](rm)); hl=truth_loose.scorer(truth_loose.build(rm))
rows=list(csv.DictReader(open(tsv),delimiter='\t'))
for r in rows:
    a=(r['seq'],int(r['lltr_start']),int(r['lltr_end']),int(r['rltr_start']),int(r['rltr_end'])); r['real']=hs(*a) is not None or hl(*a) is not None
    r['tsd']=int(r['tsd_len'])
real=[r for r in rows if r['real']]; other=[r for r in rows if not r['real']]
dom=lambda r: r['dom_order']=='LTR'; veto=lambda r: r['dom_order'] in ('notLTR','mixture')
nodom=lambda rs: [r for r in rs if not dom(r) and not veto(r)]
print(f"{len(real)} real, {len(other)} other; rows below count candidates WITHOUT an LTR domain hit (what the motif route adds)")
for L in range(4,8):
    print(f"  TSD {L}: real {sum(1 for r in nodom(real) if r['tsd']==L)} (TGCA {sum(1 for r in nodom(real) if r['tsd']==L and r['motif']=='TGCA')})  other {sum(1 for r in nodom(other) if r['tsd']==L)} (TGCA {sum(1 for r in nodom(other) if r['tsd']==L and r['motif']=='TGCA')})")
for name,pred in [("TGCA & TSD>=5",lambda r: r['motif']=='TGCA' and r['tsd']>=5),("TGCA & TSD==4",lambda r: r['motif']=='TGCA' and r['tsd']==4),("TGCA & TSD>=4",lambda r: r['motif']=='TGCA' and r['tsd']>=4),("TGCA alone",lambda r: r['motif']=='TGCA'),("TSD==5 alone",lambda r: r['tsd']==5)]:
    tp=sum(1 for r in nodom(real) if pred(r)); fp=sum(1 for r in nodom(other) if pred(r))
    print(f"  {name:16s} adds real {tp:4d} other {fp:4d} precision {100*tp/max(1,tp+fp):5.1f}%")
