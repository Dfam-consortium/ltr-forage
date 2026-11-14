#!/usr/bin/env python3
"""Drop LTR models that took over a cut segment's copies and have no LTR termini (ltr_library.py step 22b).
usage: depth_jump.py <models.fa> <depth_before.tsv> <depth_after.tsv> <cand.tsv> <ratio> <min_depth> <dropped.tsv>

Step 22b cuts high-depth segments out of models; the sample loci a cut segment held go to whichever model
matches them next.  A real LTR nested in another family's internal gains its solo copies that way (zebra finch
fam0093, median depth 3 -> 509, TSD and TG..CA); so does an "LTR" model that is really part of another repeat
(mm10: two LTR models of L1 sequence took over an internal's cut L1 segment, 99 -> 10,031 and 46 -> 2,730, and
masked 44 Mb outside the curated LTRs).  The depth search is repeated on the cut models; an LTR model whose
median depth rose <ratio>-fold, to at least <min_depth>, is dropped when its element has neither a TSD nor
TG..CA ends.  models.fa is rewritten without them.  (ltr_library.py then takes the depth once more without the
dropped models, whose copies move on to the next match in turn.)"""
import sys, csv, statistics
mf, bf, af, tsv, ratio, mind, outf = sys.argv[1:8]; ratio, mind = float(ratio), float(mind)

def medians(p):
    d = {}
    for l in open(p):
        n, v = l.rstrip('\n').split('\t')[:2]
        x = [int(t) for t in v.split(',') if t]
        d[n] = statistics.median(x) if x else 0
    return d
before, after = medians(bf), medians(af)
evid = {}   # LTR record -> the element has a TSD or TG..CA ends
for r in csv.DictReader(open(tsv), delimiter='\t'):
    if r['keep'] != '1': continue
    ok = int(r['tsd_len'] or 0) > 0 or r['motif'] == 'TGCA'
    for a, b in (('lltr_start', 'lltr_end'), ('rltr_start', 'rltr_end')):
        evid[f"{r['seq']}:{r[a]}..{r[b]}_LTR"] = ok

recs, drop = [], []
for l in open(mf):
    if l.startswith('>'): recs.append([l, ''])
    else: recs[-1][1] += l
for h, _ in recs:
    n = h[1:].split()[0].split('#')[0]
    if not n.endswith('_LTR') or evid.get(n, True): continue
    b, a = before.get(n, 0), after.get(n, 0)
    if a >= mind and a >= ratio * max(b, 1): drop.append((n, b, a, round(a / max(b, 1), 1)))
gone = {d[0] for d in drop}
with open(mf, 'w') as fh:
    for h, s in recs:
        if h[1:].split()[0].split('#')[0] not in gone: fh.write(h + s)
with open(outf, 'w') as fh:
    fh.write('model\tdepth_before\tdepth_after\tfold\n')
    for d in drop: fh.write('\t'.join(map(str, d)) + '\n')
print(f'{len(drop)} LTR model(s) without TSD or TG..CA whose depth rose {ratio:g}-fold after the cuts dropped')
