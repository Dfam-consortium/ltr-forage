#!/usr/bin/env python3
"""Cut nested sequence out of internal records, from an all-vs-all of the library.

Idea 1 (LTR_retriever's cleanup_nestedIN.pl): an internal record of another cluster that
is >= --int-cov contained in an internal record is cut out of it.
Idea 2 (Module 6 analog): an LTR record of another cluster that is >= --ltr-cov contained
in an internal record (a solo LTR or a nested element's LTR) is cut out of it.
Only internal records are cut, and a cut that would remove more than --max-frac of a
record is skipped (that is redundancy, which containment handles, not nesting).
Records shorter than --min-len after cutting are dropped.

usage: excise_contained.py <lib.fa> <clusters.json> <hits.tsv> <out.fa> <cuts.tsv>
hits.tsv: blastn -outfmt '6 qseqid sseqid pident length qstart qend sstart send qlen slen'
"""
import sys, re, json, collections, argparse
ap = argparse.ArgumentParser()
for x in ('lib', 'clusters', 'hits', 'out', 'cuts'): ap.add_argument(x)
ap.add_argument('--int-cov', type=float, default=0.95)
ap.add_argument('--ltr-cov', type=float, default=0.90)
ap.add_argument('--min-id', type=float, default=80)
ap.add_argument('--max-frac', type=float, default=0.7)
ap.add_argument('--min-len', type=int, default=100)
ap.add_argument('--drop-multi', type=int, default=0, help='drop a single-element cluster whose internal has this many or more nested stretches, or nested sequence over --max-frac (0: off)')
a = ap.parse_args()
seqs, order, n = {}, [], None
for l in open(a.lib):
    if l.startswith('>'): n = l[1:].split()[0]; seqs[n] = []; order.append(n)
    else: seqs[n].append(l.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}
clus = json.load(open(a.clusters))
by = collections.defaultdict(list)
for e, c in clus.items():
    m = re.match(r'(\S+):(\d+)\.\.(\d+)', e); by[m[1]].append((int(m[2]), int(m[3]), c))
def cl(r):
    m = re.match(r'(\S+):(\d+)\.\.(\d+)_', r)
    for x, y, c in by[m[1]]:
        if x <= int(m[2]) and int(m[3]) <= y: return c
rc = {r: cl(r) for r in order}
size = collections.Counter(clus.values())
# same-family guard: clusters whose LTR records align to each other (>= min-id over half the shorter)
same = set()
for l in open(a.hits):
    q, s, pid, ln, qs, qe, ss, se, ql, sl = l.split()[:10]
    if '_LTR' in q and '_LTR' in s and q in rc and s in rc and rc[q] != rc[s] and float(pid) >= a.min_id \
       and int(ln) >= 0.5 * min(int(ql), int(sl)):
        same.add((rc[q], rc[s])); same.add((rc[s], rc[q]))
nguard = 0
cuts = collections.defaultdict(list)
for l in open(a.hits):
    q, s, pid, ln, qs, qe, ss, se, ql, sl = l.split()[:10]
    if '_INT' not in s or q not in seqs or s not in seqs or rc[q] == rc[s] or float(pid) < a.min_id: continue
    cov = (abs(int(qe) - int(qs)) + 1) / int(ql)
    need = a.int_cov if '_INT' in q else a.ltr_cov
    if cov >= need and '_INT' in q and (rc[q], rc[s]) in same:
        nguard += 1; continue
    if cov >= need:
        cuts[s].append((min(int(ss), int(se)), max(int(ss), int(se)), q))
def merged(r):
    iv = []
    for x, y, q in sorted(cuts[r]):
        if iv and x <= iv[-1][1] + 1: iv[-1][1] = max(iv[-1][1], y); iv[-1][2].add(q)
        else: iv.append([x, y, {q}])
    return iv
dropc = set()
if a.drop_multi:
    for r in cuts:
        iv = merged(r)
        if size[rc[r]] == 1 and (len(iv) >= a.drop_multi or sum(y - x + 1 for x, y, _ in iv) > a.max_frac * len(seqs[r])):
            dropc.add(rc[r])
nrec = ncut = nbase = ndrop = nskip = 0
with open(a.out, 'w') as out, open(a.cuts, 'w') as ct:
    for r in order:
        if rc[r] in dropc:
            ct.write(f'{r}\tDROP\t\t\n'); continue
        sq = seqs[r]
        if r in cuts:
            iv = []
            for x, y, q in sorted(cuts[r]):
                if iv and x <= iv[-1][1] + 1: iv[-1][1] = max(iv[-1][1], y); iv[-1][2].add(q)
                else: iv.append([x, y, {q}])
            tot = sum(y - x + 1 for x, y, _ in iv)
            if tot > a.max_frac * len(sq): nskip += 1
            else:
                keep, last = [], 0
                for x, y, qs_ in iv:
                    keep.append(sq[last:x-1]); last = y
                    ct.write(f'{r}\t{x}\t{y}\t{",".join(sorted(qs_))}\n')
                keep.append(sq[last:]); sq = ''.join(keep)
                nrec += 1; ncut += len(iv); nbase += tot
        if len(sq) >= a.min_len: out.write(f'>{r}\n{sq}\n')
        else: ndrop += 1
print(f'{len(dropc)} single-element clusters dropped; {nrec} internal records cut ({ncut} stretches, {nbase:,} bp); {nskip} skipped as over {a.max_frac}; {ndrop} records dropped; {nguard} same-family hits skipped', file=sys.stderr)
