#!/usr/bin/env python3
"""Split a recipe library into a core library and a candidates pile.
usage: tier.py <final.fa> <clusters.json> <cand.tsv> <core.fa> <candidates.fa>
               [--min-id 90] [--max-len 9000] [--pbs pbs.tsv] [--strict] [--both-ltrs]

Multi-element clusters go to core.  A singleton cluster goes to core when its
element is <= max-len long and:
  without --pbs:  LTR identity >= min-id
  with --pbs:     LTR identity >= min-id OR a tRNA PBS (pbs_mm12 <= 1)
  with --strict:  LTR identity >= min-id AND a tRNA PBS
Rice singletons (classes from classify2.py): id>=90 keeps 87 of 108 clean
sole-representative families and 50 of 213 nested/chimeric/unannotated ones;
--min-id 93 --pbs keeps 92 and 42; --min-id 90 --pbs --strict keeps 51 and 23.

Recommended: --min-id 93 --pbs (the "PBS tier"; library_recipe.sh step 13).  Over
strict it adds +357 kb TP / +207 kb FP in rice, +84/+77 kb in Arabidopsis and
+3.9/+4.5 Mb in lotus, where ~1/4 of the extra FP is non-LTR (mostly DNA TE)
passenger sequence and the rest looks like LTR material missing from the reference.

Each family is written as one LTR record and one internal record (the 3' LTR
record is dropped; it is a near-identical copy of the 5' one), unless --both-ltrs."""
import sys, re, json, csv, bisect, collections, argparse
ap = argparse.ArgumentParser()
for p in ('fa', 'clusters', 'tsv', 'core', 'cand'): ap.add_argument(p)
ap.add_argument('--min-id', type=float, default=90); ap.add_argument('--max-len', type=int, default=9000)
ap.add_argument('--pbs'); ap.add_argument('--strict', action='store_true'); ap.add_argument('--both-ltrs', action='store_true')
a = ap.parse_args()
cl = json.load(open(a.clusters)); size = collections.Counter(cl.values())
feat = {f"{x['seq']}:{x['elem_start']}..{x['elem_end']}": x for x in csv.DictReader(open(a.tsv), delimiter='\t')}
pbs = {x['element']: x['pbs_mm12'] != '.' for x in csv.DictReader(open(a.pbs), delimiter='\t')} if a.pbs else {}
iv = collections.defaultdict(list)
for k in cl:
    c, s, e = re.match(r'(.+):(\d+)\.\.(\d+)$', k).groups(); iv[c].append((int(s), int(e), k))
for v in iv.values(): v.sort()
def element(name):
    c, s = re.match(r'(.+):(\d+)\.\.', name).groups(); s = int(s); v = iv[c]
    i = bisect.bisect_right(v, (s, float('inf'), '')) - 1
    while i >= 0 and v[i][1] < s: i -= 1
    return v[i][2]
def admit(k):
    x = feat[k]
    if int(x['elem_len']) > a.max_len: return False
    young = float(x['identity']) >= a.min_id
    if not a.pbs: return young
    return (young and pbs.get(k, False)) if a.strict else (young or pbs.get(k, False))
recs = []
for l in open(a.fa):
    if l.startswith('>'): recs.append([l, []])
    else: recs[-1][1].append(l)
out = {True: open(a.core, 'w'), False: open(a.cand, 'w')}; why = collections.Counter(); seen = {}; nltr = collections.Counter()
for h, s in recs:
    name = h[1:].split()[0]; k = element(name)
    if k not in seen:
        multi = size[cl[k]] > 1; ok = multi or admit(k)
        seen[k] = ok; why['multi' if multi else 'singleton kept' if ok else 'singleton -> candidates'] += 1
    if '_LTR' in name.split('#')[0] and not a.both_ltrs:
        nltr[k] += 1
        if nltr[k] > 1: continue
    out[seen[k]].write(h + ''.join(s))
print(dict(why), file=sys.stderr)
