#!/usr/bin/env python3
"""Experiment 1 follow-up (2026-10-19): label each element's LTRs and internal separately with curated
model names (e.g. Os0028_LTR / Os0028_INT) and classify the links of diag_links.py.
usage: diag_parts.py <run_dir> <standard.out> <links.tsv>
Each part gets the curated model covering most of it (>= 50% covered by that family's models).
Link classes: same (LTR and INT families agree), LTR-swap (same INT family, different LTR family),
INT-swap (same LTR family, different INT family), both-differ; links with an unlabelled part are 'n/a'."""
import sys, csv, re, bisect, collections, statistics as S
run, std, links = sys.argv[1:4]
iv = collections.defaultdict(list)
for l in open(std):
    c = l.split()
    if len(c) < 15 or not c[0].isdigit() or not c[10].startswith('LTR'): continue
    m = re.match(r'(Os\d+)_(LTR|INT)', c[9])
    if m: iv[c[4]].append((int(c[5]), int(c[6]), m[1], m[2]))
for s in iv: iv[s].sort()
st = {s: [x[0] for x in v] for s, v in iv.items()}
mx = max(b - a for v in iv.values() for a, b, *_ in v)
def lab(seq, s, e, part):
    v = iv.get(seq, []); i = bisect.bisect_left(st.get(seq, []), s - mx); cov = collections.Counter()
    for a, b, f, p in v[i:]:
        if a > e: break
        o = min(b, e) - max(a, s) + 1
        if o > 0 and p == part: cov[f] += o
    if not cov: return None
    f, n = cov.most_common(1)[0]
    return f if n >= 0.5 * (e - s + 1) else None
part = {}
for r in csv.DictReader(open(f'{run}/cand.tsv'), delimiter='\t'):
    if r['keep'] != '1': continue
    s = r['seq']; k = f"{s}:{r['elem_start']}..{r['elem_end']}"
    L = lab(s, int(r['lltr_start']), int(r['lltr_end']), 'LTR') or lab(s, int(r['rltr_start']), int(r['rltr_end']), 'LTR')
    I = lab(s, int(r['lltr_end']) + 1, int(r['rltr_start']) - 1, 'INT')
    part[k] = (L, I)
cls = collections.Counter(); per = collections.defaultdict(collections.Counter); stat = collections.defaultdict(list)
for r in csv.DictReader(open(links), delimiter='\t'):
    (La, Ia), (Lb, Ib) = part[r['a']], part[r['b']]
    if not (La and Ia and Lb and Ib): c = 'n/a'
    elif La == Lb and Ia == Ib: c = 'same'
    elif Ia == Ib: c = 'LTR-swap'
    elif La == Lb: c = 'INT-swap'
    else: c = 'both-differ'
    cls[c] += 1; stat[c].append((float(r['ltr_pid']), float(r['ltr_idcov_shorter']), float(r['int_pid'])))
    if c in ('LTR-swap', 'INT-swap', 'both-differ'): per[c][tuple(sorted(((La, Ia), (Lb, Ib))))] += 1
print('links by part labels: ' + ', '.join(f'{k} {v}' for k, v in cls.most_common()))
for c in ('same', 'LTR-swap', 'INT-swap', 'both-differ'):
    v = stat[c]
    if not v: continue
    print(f'{c:12s} n={len(v):6d}  median LTR pid {S.median(x[0] for x in v):.1f}  LTR idcov {S.median(x[1] for x in v):.3f}  '
          f'INT pid {S.median(x[2] for x in v):.1f}  | LTR idcov<0.85 {sum(x[1] < 0.85 for x in v)}')
for c in ('LTR-swap', 'INT-swap', 'both-differ'):
    print(f'top {c}:', '; '.join(f'{a[0]}L+{a[1]}I ~ {b[0]}L+{b[1]}I x{n}' for (a, b), n in per[c].most_common(8)))
