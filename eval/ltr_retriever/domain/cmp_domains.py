"""Compare the dom_hits column of an ltr-forage table with an hmmsearch --tblout on the
same candidates.  Usage: cmp_domains.py table.tsv hits.tbl"""
import csv, sys, re, collections

tsv, tbl = sys.argv[1], sys.argv[2]
ours=list(csv.DictReader(open(tsv),delimiter='\t'))
# hmmsearch tblout: target name like chr22:1255845..1260930|chr22:1255895..1260880|aa1 ?
hm=collections.defaultdict(set)
names=set()
for l in open(tbl):
    if l.startswith('#'): continue
    f=l.split()
    tname, acc, ev = f[0], f[3].split('.')[0], float(f[4])
    m=re.match(r'(\S+?):(\d+)\.\.(\d+)\|', tname)
    if not m: m=re.match(r'(\S+?):(\d+)\.\.(\d+)', tname)
    key=(m.group(1),int(m.group(2)),int(m.group(3)))
    hm[key].add(acc); names.add(tname)
print('hmmsearch: %d hits on %d candidates' % (sum(len(v) for v in hm.values()), len(hm)))
def ov(a,b): return max(0, min(a[2],b[2])-max(a[1],b[1])+1)
ourhits={}
for r in ours:
    if r['dom_order']=='none': continue
    key=(r['seq'],int(r['elem_start']),int(r['elem_end']))
    ourhits[key]=set(h.split(':')[0] for h in r['dom_hits'].split(','))
print('ours: %d hits on %d candidates' % (sum(len(v) for v in ourhits.values()), len(ourhits)))
matched=0; shared=0; hm_only=0; our_only=0
used=set()
for k,accs in hm.items():
    best=None
    for k2 in ourhits:
        if k2[0]==k[0] and ov(k,k2)>0.5*min(k[2]-k[1],k2[2]-k2[1]):
            best=k2; break
    if best is None:
        # is there any candidate of ours overlapping at all?
        anyc=[r for r in ours if r['seq']==k[0] and ov(k,(r['seq'],int(r['elem_start']),int(r['elem_end'])))>0]
        print('hmmsearch-only', k, sorted(accs), 'our candidate present' if anyc else 'no overlapping candidate')
        hm_only+=len(accs); continue
    used.add(best); matched+=1
    shared+=len(accs & ourhits[best]); hm_only+=len(accs-ourhits[best]); our_only+=len(ourhits[best]-accs)
    if accs-ourhits[best] or ourhits[best]-accs: print('diff', k, 'hmm-only', sorted(accs-ourhits[best]), 'ours-only', sorted(ourhits[best]-accs))
for k2 in ourhits:
    if k2 not in used:
        our_only+=len(ourhits[k2]); print('ours-only candidate', k2, sorted(ourhits[k2]))
print(f'candidates matched {matched}; accessions shared {shared}, hmmsearch-only {hm_only}, ours-only {our_only}')
