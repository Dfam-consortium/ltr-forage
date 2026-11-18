#!/usr/bin/env python3
"""Group a tiered library (one LTR + one INT record per element cluster) into 80-80-80 families.

usage: family.py <lib.fa> <ava.tsv> <out.fa> <families.tsv> [--rule or|and] [--gap-aware] [--cand cand.tsv]
                 [--dedup-id 90] [--dedup-cov 0.9]

ava.tsv: dc-megablast all-vs-all of lib.fa, outfmt
  "6 qseqid sseqid pident length qstart qend sstart send qlen slen bitscore".
Two elements are in the same family when their LTRs, or their internal regions, are
>= 80% identical over >= 80% of the shorter record (Wicker et al. 2007), single linkage;
with --rule and, both the LTRs and the internal regions must match.

--gap-aware (for gappy assemblies; lotus): internal regions are compared only over the
positions where neither record is N (one record's gaps projected onto the other through
the pair's main alignment diagonal), and under --rule and an internal comparison that
cannot be judged -- under --min-comparable shared gap-free bases, or a record more than
--max-int-gap N -- does not block a link the LTRs make.  Lotus strict: --rule and gives
710 families, --gap-aware 568, LTR-only grouping 542; rice and Arabidopsis are unchanged.
The recipe's own clustering needs both to match, which splits lotus families with
variable internals ~3-fold (1,445 strict "families" -> 426).  Coverage is of the
record's non-N length, so assembly gaps do not split families.

Within a family every distinct variant is kept for masking; a record is dropped only
when a longer record of the same kind in the family covers >= dedup-cov of it at
>= dedup-id identity.  Records are renamed fam<N>_v<k>_LTR / _INT (origin=<old name>).
With --cand (the ltr-forage table), names carry a RepeatMasker class, fam<N>_v<k>_LTR#LTR/Gypsy:
the dom_superfamily of the family's elements when >= --class-majority of those with a Gypsy or
Copia call agree, otherwise plain LTR."""
import sys, argparse, collections
ap = argparse.ArgumentParser()
for p in ('lib', 'ava', 'out', 'tsv'): ap.add_argument(p)
ap.add_argument('--rule', choices=('or', 'and'), default='or')
ap.add_argument('--gap-aware', action='store_true', help='measure coverage only where neither record has N (for gappy assemblies such as lotus)')
ap.add_argument('--min-comparable', type=int, default=500)
ap.add_argument('--max-int-gap', type=float, default=0.2, help='with --gap-aware: an internal record more than this fraction N cannot block a link the LTRs make')
ap.add_argument('--cand', help='ltr-forage table: class names by the elements\' dom_superfamily')
ap.add_argument('--class-majority', type=float, default=0.8)
ap.add_argument('--dedup-id', type=float, default=90); ap.add_argument('--dedup-cov', type=float, default=0.9)
args = a = ap.parse_args()

seq = {}; order = []; n = None
for l in open(a.lib):
    if l.startswith('>'): n = l[1:].split()[0]; seq[n] = []; order.append(n)
    else: seq[n].append(l.strip())
seq = {k: ''.join(v) for k, v in seq.items()}
kind = lambda r: 'INT' if r.endswith('_INT') else 'LTR'
# element of each record.  With --cand: the kept element of the ltr-forage table one of whose records (5' LTR,
# internal region, 3' LTR) has this name, which also gives the element's superfamily.  Otherwise, and for a
# record the table does not name: an INT named <seq>:<s>..<e>_INT belongs to the LTR record that ends at s-1
# on the same sequence (tier.py keeps the 5' LTR), and any other record stands alone.
import re, csv
pos = lambda r: re.match(r'(.+):(\d+)\.\.(\d+)_', r).groups()
rec_elem = {}; sf = {}
if a.cand:
    for x in csv.DictReader(open(a.cand), delimiter='\t'):
        if x['keep'] != '1': continue
        c = x['seq']; e = f"{c}:{x['elem_start']}..{x['elem_end']}"; sf[e] = x['dom_superfamily']
        ls, le, rs, re_ = int(x['lltr_start']), int(x['lltr_end']), int(x['rltr_start']), int(x['rltr_end'])
        for nm in (f'{c}:{ls}..{le}_LTR', f'{c}:{le + 1}..{rs - 1}_INT', f'{c}:{rs}..{re_}_LTR'): rec_elem.setdefault(nm, e)
ltr_end = {}
for r in order:
    if kind(r) == 'LTR': c, s, e = pos(r); ltr_end[(c, int(e))] = r
e_of = {}
for r in order:
    if r in rec_elem: e_of[r] = rec_elem[r]
    elif kind(r) == 'LTR': e_of[r] = r
    else: c, s, e = pos(r); l = ltr_end.get((c, int(s) - 1), r); e_of[r] = rec_elem.get(l, l)

def union(iv):
    iv = sorted(iv); t = 0; cs = ce = None
    for s, e in iv:
        if cs is None: cs, ce = s, e
        elif s <= ce + 1: ce = max(ce, e)
        else: t += ce - cs + 1; cs, ce = s, e
    return t + (ce - cs + 1 if cs is not None else 0)
hsp80 = collections.defaultdict(list)
cov80 = collections.defaultdict(lambda: [[], []]); covd = collections.defaultdict(lambda: [[], []])
for l in open(a.ava):
    q, s, pid, ln, qs, qe, ss, se, ql, sl, b = l.split('\t')
    if q == s or q not in seq or s not in seq or kind(q) != kind(s): continue
    qi = (min(int(qs), int(qe)), max(int(qs), int(qe))); si = (min(int(ss), int(se)), max(int(ss), int(se)))
    if float(pid) >= 80: c = cov80[(q, s)]; c[0].append(qi); c[1].append(si); hsp80[(q, s)].append((int(qs), int(qe), int(ss), int(se)))
    if float(pid) >= a.dedup_id: c = covd[(q, s)]; c[0].append(qi); c[1].append(si)
# coverage is measured on the non-N length: lotus internal records are a median 26% assembly
# gap, which cannot align, and counting it split one internal family into dozens
L = {r: max(1, len(x) - x.upper().count('N')) for r, x in seq.items()}
par = {e: e for e in set(e_of.values())}
def find(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
# --gap-aware: project each record's N runs onto the other through the pair's dominant
# alignment diagonal and count only positions that are non-N in both.  Two copies of one
# internal region that each lost different parts to assembly gaps then compare over what
# they share.  The projection uses one diagonal, so it is approximate across indels.
import numpy as np
okmask = {r: np.frombuffer(x.upper().encode(), dtype=np.uint8) != ord('N') for r, x in seq.items()} if a.gap_aware else {}
def comparable(q, s):
    hs = hsp80[(q, s)]
    plus = [h for h in hs if h[2] <= h[3]]; minus = [h for h in hs if h[2] > h[3]]
    w = lambda v: sum(abs(h[1] - h[0]) for h in v)
    mq, ms = okmask[q], okmask[s]
    if w(plus) >= w(minus):
        d = int(np.median([h[2] - h[0] for h in plus]))           # s = q + d (1-based)
        lo, hi = max(1, 1 - d), min(len(mq), len(ms) - d)
        if hi < lo: return 0, 0
        both = mq[lo - 1:hi] & ms[lo - 1 + d:hi + d]
    else:
        t = int(np.median([h[2] + h[0] for h in minus]))           # s = t - q
        lo, hi = max(1, t - len(ms)), min(len(mq), t - 1)
        if hi < lo: return 0, 0
        both = mq[lo - 1:hi] & ms[t - hi - 1:t - lo][::-1]
    n = int(both.sum())
    return n, n                                                    # same count on each side
passed = collections.defaultdict(set)   # element pair -> record kinds that pass 80-80
uninformative = set()
gapfrac = {r: x.upper().count('N') / max(1, len(x)) for r, x in seq.items()}
for (q, s), (qi, si) in cov80.items():
    if e_of[q] == e_of[s]: continue
    if a.gap_aware:
        # the shared-region test needs enough comparable sequence; short records (LTRs)
        # and pairs that barely overlap fall back to the plain test
        nq, ns = comparable(q, s)
        ok = (union(qi) / L[q] if L[q] <= L[s] else union(si) / L[s]) >= 0.8
        if not ok and nq >= a.min_comparable:
            ok = min(union(qi) / nq, union(si) / ns) >= 0.8
        # internals too gappy, or sharing too little gap-free sequence, to judge: with
        # --rule and they then do not block a link the LTRs make
        if not ok and kind(q) == 'INT' and nq < a.min_comparable:
            uninformative.add(frozenset((e_of[q], e_of[s])))
    else:
        ok = (union(qi) / L[q] if L[q] <= L[s] else union(si) / L[s]) >= 0.8
    if ok: passed[frozenset((e_of[q], e_of[s]))].add(kind(q))
has = collections.defaultdict(set)
for r in order: has[e_of[r]].add(kind(r))
intrec = {e_of[r]: r for r in order if kind(r) == 'INT'}
for pair, kinds in passed.items():
    x, y = tuple(pair)
    # --rule or: the LTRs or the internals match (Wicker); --rule and: every record kind
    # both elements have must match (an element with no INT record links on its LTR)
    need = has[x] & has[y]
    if pair in uninformative or (a.gap_aware and 'INT' in need and max(gapfrac[intrec[x]], gapfrac[intrec[y]]) > a.max_int_gap):
        need = need - {'INT'}
    if args.rule == 'or' or kinds >= need:
        par[find(x)] = find(y)
fam = collections.defaultdict(list)
for r in order: fam[find(e_of[r])].append(r)
# family class from the domain superfamily of the family's elements (kept rows of the table only)
def fam_class(recs):
    votes = collections.Counter(sf.get(e) for e in {e_of[r] for r in recs})
    gc = {k: v for k, v in votes.items() if k in ('Gypsy', 'Copia')}
    if gc:
        k, v = max(gc.items(), key=lambda x: x[1])
        if v >= a.class_majority * sum(gc.values()): return 'LTR/' + k
    return 'LTR'
# within-family near-duplicate removal, longest first
covered = lambda r, by: union((covd.get((r, by)) or [[], []])[0] + (covd.get((by, r)) or [[], []])[1]) / L[r]
out = open(a.out, 'w'); tsv = open(a.tsv, 'w'); tsv.write('family\tnew_name\torigin\tkind\tlen\tstatus\tclass\n')
rank = {r: i for i, r in enumerate(order)}
fams = sorted(fam.values(), key=lambda v: rank[v[0]])
nkept = 0
ncls = collections.Counter()
for i, recs in enumerate(fams, 1):
    cl = fam_class(recs); ncls[cl] += 1
    for k in ('LTR', 'INT'):
        kept = []
        for r in sorted((r for r in recs if kind(r) == k), key=lambda r: -L[r]):
            by = next((x for x in kept if covered(r, x) >= a.dedup_cov), None)
            if by:
                tsv.write(f'fam{i:04d}\t.\t{r}\t{k}\t{L[r]}\tdup_of:{by}\t{cl}\n'); continue
            kept.append(r); nm = f'fam{i:04d}_v{len(kept)}_{k}'; nkept += 1
            out.write(f'>{nm}#{cl} origin={r}\n{seq[r]}\n'); tsv.write(f'fam{i:04d}\t{nm}\t{r}\t{k}\t{L[r]}\tkept\t{cl}\n')
print(f'{len(set(e_of.values()))} element clusters -> {len(fams)} families; {len(seq)} records -> {nkept} kept', file=sys.stderr)
if a.cand: print('classes:', dict(ncls.most_common()), file=sys.stderr)
