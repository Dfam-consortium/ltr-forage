import csv, sys, collections, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import truth_loose
def isint(n): return n.endswith(('-int','_I','-I','_int')) or '_INT' in n or (BARE_I and n.endswith('I'))
truth_loose.isint = isint
ns = {'collections': collections, 'isint': isint}
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'score_tsv.py')).read()
exec(src[src.index('def strict'):src.index('hit_s = ')], ns)
rm, selftsv, m4tsv = sys.argv[1], sys.argv[2], sys.argv[3]
BARE_I = truth_loose.bare_i(rm)   # False where the annotation marks internals with -int (MLT1I is an LTR)
hit_s = truth_loose.scorer(ns['strict'](rm)); hit_l = truth_loose.scorer(truth_loose.build(rm))
cov = {}
for r in csv.DictReader(open(selftsv), delimiter='\t'):
    cov[(r['seq'], r['elem_start'], r['elem_end'])] = max(float(r['lltr_self_cov']), float(r['rltr_self_cov']))
rows = list(csv.DictReader(open(m4tsv), delimiter='\t'))
dom = lambda r: r['dom_order'] == 'LTR'; veto = lambda r: r['dom_order'] in ('notLTR', 'mixture')
n = 0
for r in rows:
    a = (r['seq'], int(r['lltr_start']), int(r['lltr_end']), int(r['rltr_start']), int(r['rltr_end']))
    r['real'] = hit_s(*a) is not None or hit_l(*a) is not None
    r['cov'] = cov.get((r['seq'], r['elem_start'], r['elem_end']))
    r['std'] = not veto(r) and (dom(r) or (r['motif'] == 'TGCA' and int(r['tsd_len']) >= 5))
    n += r['cov'] is None
print('unmatched', n)
rows = [r for r in rows if r['cov'] is not None]
bins = [(0, 0.0001), (0.0001, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1.01)]
print("self-cov bin: all real / all other | kept-at-standard real / other")
for lo, hi in bins:
    f = lambda rs: sum(1 for r in rs if lo <= r['cov'] < hi)
    print(f"  [{lo:.2f},{hi:.2f}): {f([r for r in rows if r['real']]):5d} / {f([r for r in rows if not r['real']]):5d} | {f([r for r in rows if r['real'] and r['std']]):5d} / {f([r for r in rows if not r['real'] and r['std']]):5d}")
for thr in (0.5, 0.7, 0.9, 0.95, 0.98, 0.99):
    lost = [r for r in rows if r['real'] and r['std'] and r['cov'] > thr]
    junk = [r for r in rows if not r['real'] and r['std'] and r['cov'] > thr]
    print(f"threshold {thr}: removes {len(lost)} real and {len(junk)} other kept candidates")
lost = [r for r in rows if r['real'] and r['std'] and r['cov'] > 0.5]
print("real kept candidates removed at 0.5:")
for r in lost[:40]:
    print(f"  {r['seq']}:{r['elem_start']}-{r['elem_end']} ltr {r['lltr_len']}/{r['rltr_len']} id {r['identity']} cov {r['cov']:.2f} {r['dom_hits'][:40]}")

# finer bins near the top, with RepeatMasker annotation of the flagged kept candidates
import bisect
ann = collections.defaultdict(list)
for line in open(rm):
    c = line.split()
    if len(c) < 15 or not c[0].isdigit(): continue
    ann[c[4]].append((int(c[5]), int(c[6]), c[9], c[10]))
for v in ann.values(): v.sort()
starts = {k: [a[0] for a in v] for k, v in ann.items()}
def desc(chrom, s, e):
    v = ann.get(chrom, []); k = bisect.bisect_left(starts.get(chrom, []), s - 200000); d = collections.Counter()
    while k < len(v) and v[k][0] <= e:
        o = min(e, v[k][1]) - max(s, v[k][0]) + 1
        if o > 0: d[f"{v[k][2]}({v[k][3]})"] += o
        k += 1
    return ', '.join(f"{a}:{b}" for a, b in d.most_common(2))
print("fine bins, kept-at-standard real / other:")
for lo in (0.8, 0.85, 0.9, 0.95, 0.98, 0.99):
    hi = {0.95: 0.98, 0.98: 0.99, 0.99: 1.01}.get(lo, lo + 0.05)
    print(f"  [{lo:.2f},{hi:.2f}): {sum(1 for r in rows if r['real'] and r['std'] and lo <= r['cov'] < hi)} / {sum(1 for r in rows if not r['real'] and r['std'] and lo <= r['cov'] < hi)}")
print("kept candidates with cov > 0.8, with LTR annotation:")
for r in sorted([r for r in rows if r['std'] and r['cov'] >= 0.95], key=lambda r: -r['cov']):
    print(f"  {'REAL ' if r['real'] else 'other'} {r['seq']}:{r['elem_start']}-{r['elem_end']} ltr {r['lltr_len']} cov {r['cov']:.2f} L[{desc(r['seq'], int(r['lltr_start']), int(r['lltr_end']))}] R[{desc(r['seq'], int(r['rltr_start']), int(r['rltr_end']))}]")
