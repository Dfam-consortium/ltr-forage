"""Per-candidate features on the trimmed Arabidopsis boundaries: truth (strict|loose), TSD at ends (+-1),
TG..CA (+-1), LTR/internal lengths, identity, and internal-region 20-mer coverage/gap against the genome index."""
import sys, collections, bisect, subprocess, statistics
import truth2
import os
ROOT = os.environ.get('LTR_WORKSPACE') or os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..'))
fai={}; names=[]
for l in open('rice.fa.fai'):
    c=l.split(); fai[c[0]]=tuple(map(int,c[1:5])); names.append(c[0])
fh=open('rice.fa','rb')
def fetch(chrom,s,e):
    ln,off,lb,lw=fai[chrom]; s=max(1,s); e=min(ln,e); out=b''; pos=s-1
    while pos<e:
        line,col=divmod(pos,lb); fh.seek(off+line*lw+col); n=min(lb-col,e-pos); out+=fh.read(n); pos+=n
    return out.decode().upper()
# strict truth (truth.py rules)
ann=collections.defaultdict(list)
for line in open('rice.out'):
    c=line.split()
    if len(c)<15 or not c[0].isdigit() or not c[10].startswith('LTR'): continue
    ann[c[4]].append((int(c[5]),int(c[6]),c[8],c[9]))
def isint(n): return '_INT' in n.upper() or n.endswith(('-int','_I','-I','_int','I'))
strict=collections.defaultdict(list)
for chrom,lst in ann.items():
    lst.sort()
    for i,(s1,e1,st1,nm1) in enumerate(lst):
        if isint(nm1) or not (100<=e1-s1+1<=1000): continue
        seen=False
        for j in range(i+1,len(lst)):
            s2,e2,st2,nm2=lst[j]
            if s2-s1>15000: break
            if nm2==nm1 and st2==st1:
                if seen and s2-s1>=1000 and 100<=e2-s2+1<=1000: strict[chrom].append((s1,e1,s2,e2))
                break
            if nm2!=nm1 and isint(nm2) and st2==st1: seen=True
for v in strict.values(): v.sort()
hs=truth2.scorer(strict); hl=truth2.scorer(truth2.build('rice.out'))
def tsd_at(chrom,L,R,slack=1):
    left=fetch(chrom,L-22,L+1); right=fetch(chrom,R-2,R+21); best=0
    for dl in range(0,slack+1):
        for dr in range(0,slack+1):
            for k in range(20,3,-1):
                a=left[22-dl-k:22-dl]; b=right[3+dr:3+dr+k]
                if len(a)==k and a==b: best=max(best,k); break
    return best
def motif_at(chrom,L,R,slack=1):
    s=fetch(chrom,L-slack,L+slack+1); e=fetch(chrom,R-slack-1,R+slack)
    return int(any(s[i:i+2]=='TG' for i in range(len(s)-1)) and any(e[i:i+2]=='CA' for i in range(len(e)-1)))
K,STRIDE=20,10; comp=str.maketrans('ACGT','TGCA')
cands=[]; tiles=set(); tl=open('tiles.txt','w')
for i,l in enumerate(open(sys.argv[1])):
    c=l.split(); chrom=names[int(c[10])]; l1s,l1e,l2s,l2e=map(int,(c[3],c[4],c[6],c[7]))
    internal=fetch(chrom,l1e+1,l2s-1); tt=[internal[p:p+K] for p in range(0,len(internal)-K+1,STRIDE)]; tt=[t for t in tt if 'N' not in t]
    for t in tt:
        for x in (t,t.translate(comp)[::-1]):
            if x not in tiles: tiles.add(x); tl.write(x+'\n')
    cands.append((i,chrom,l1s,l1e,l2s,l2e,float(c[9]),tt))
tl.close()
subprocess.run([ROOT+'/sufr/target/release/sufr','count','rice.sufr','tiles.txt'],stdout=open('tile_counts.txt','w'),stderr=subprocess.DEVNULL,check=True)
cnt={}
for l in open('tile_counts.txt'):
    c=l.split()
    if len(c)==2: cnt[c[0]]=int(c[1])
def both(t): return cnt.get(t,0)+cnt.get(t.translate(comp)[::-1],0)
out=open(sys.argv[2],'w'); out.write('idx\tchrom\tl1s\tl1e\tl2s\tl2e\tsim\tstrict\tloose\tltr\tintl\ttsd\tmotif\tntiles\tfrac\tgap\n')
for i,chrom,l1s,l1e,l2s,l2e,sim,tt in cands:
    st=int(bool(hs(chrom,l1s,l1e,l2s,l2e))); lo=int(bool(hl(chrom,l1s,l1e,l2s,l2e)))
    cs=[both(t) for t in tt]; multi=[x>=3 for x in cs]
    frac=sum(multi)/len(multi) if multi else 0
    run=best=0
    for m in multi:
        run=0 if m else run+1; best=max(best,run)
    gap=best/len(multi) if multi else 1.0
    out.write(f"{i}\t{chrom}\t{l1s}\t{l1e}\t{l2s}\t{l2e}\t{sim}\t{st}\t{lo}\t{max(l1e-l1s+1,l2e-l2s+1)}\t{l2s-l1e-1}\t{tsd_at(chrom,l1s,l2e)}\t{motif_at(chrom,l1s,l2e)}\t{len(tt)}\t{frac:.3f}\t{gap:.3f}\n")
out.close(); print('done', len(cands))
