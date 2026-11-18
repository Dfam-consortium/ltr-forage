"""Score an `ltr-forage --format tsv` table against RepeatMasker truth.

Usage: score_tsv.py <rm.out> <table.tsv> [features.tsv|-] [self.tsv]

Truth is the union of the strict harness (truth.py rules: LTR/internal/LTR of
one family and strand) and the loose one (truth_loose.py: any LTR names,
merged fragments, >= 50% internal coverage).  A candidate is real when both
of its LTRs overlap a truth element's LTRs by >= 50%.  The optional
features.tsv (signals/features.py output) supplies LTR_retriever's hmmsearch
order call per candidate for a side-by-side comparison (pass - to skip it).
The optional self.tsv is a table from the same run with the lltr_self_cov
columns (any run of the same candidates); when given, the self-repeat veto at
0.99 is applied to the --keep levels.
"""
import sys, csv, bisect, collections, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import truth_loose

rm, tsv = sys.argv[1], sys.argv[2]
feat = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] != '-' else None
selftsv = sys.argv[4] if len(sys.argv) > 4 else None

# a bare trailing I names an internal in Repbase libraries (ATCOPIA19I), but not where the annotation marks
# internals with -int (RepeatMasker's libraries: MLT1I is an LTR)
BARE_I = truth_loose.bare_i(rm)
def isint(n): return n.endswith(('-int', '_I', '-I', '_int')) or '_INT' in n or (BARE_I and n.endswith('I'))
truth_loose.isint = isint  # rice names its internals *_INT_*

def strict(rm):
    ann = collections.defaultdict(list)
    for line in open(rm):
        c = line.split()
        if len(c) < 15 or not c[0].isdigit() or not c[10].startswith('LTR'): continue
        ann[c[4]].append((int(c[5]), int(c[6]), c[8], c[9]))
    truth = collections.defaultdict(list)
    for chrom, lst in ann.items():
        lst.sort()
        for i, (s1, e1, st1, nm1) in enumerate(lst):
            if isint(nm1) or not (100 <= e1 - s1 + 1 <= 1000): continue
            seen = False
            for j in range(i + 1, len(lst)):
                s2, e2, st2, nm2 = lst[j]
                if s2 - s1 > 15000: break
                if nm2 == nm1 and st2 == st1:
                    if seen and s2 - s1 >= 1000 and 100 <= e2 - s2 + 1 <= 1000: truth[chrom].append((s1, e1, s2, e2))
                    break
                if nm2 != nm1 and isint(nm2) and st2 == st1: seen = True
    for v in truth.values(): v.sort()
    return truth

hit_s = truth_loose.scorer(strict(rm))
hit_l = truth_loose.scorer(truth_loose.build(rm))

rows = list(csv.DictReader(open(tsv), delimiter='\t'))
for r in rows:
    a = (r['seq'], int(r['lltr_start']), int(r['lltr_end']), int(r['rltr_start']), int(r['rltr_end']))
    r['real_s'] = hit_s(*a) is not None
    r['real_l'] = hit_l(*a) is not None
    r['real'] = r['real_s'] or r['real_l']
    r['gap'] = float(r['int_gap_frac']) if r['int_gap_frac'] not in ('.', '') else None
    r['tsd'] = int(r['tsd_len'])
    r['self_cov'] = max(float(r.get('lltr_self_cov', 0) or 0), float(r.get('rltr_self_cov', 0) or 0))
    core = ('PF00665', 'PF13683', 'PF03732', 'PF14223', 'PF14244', 'PF13976', 'PF07727', 'PF13650', 'PF13975', 'PF09668', 'PF00077', 'PF08284', 'PF02160', 'PF12384', 'PF00385')
    r['ltr_dom'] = sum(1 for h in (r.get('ltr_dom_hits') or '.').split(',') if h.startswith(core))
if selftsv:
    sc = {}
    for x in csv.DictReader(open(selftsv), delimiter='\t'):
        sc[(x['seq'], x['elem_start'], x['elem_end'])] = max(float(x['lltr_self_cov']), float(x['rltr_self_cov']))
    for r in rows:
        r['self_cov'] = sc.get((r['seq'], r['elem_start'], r['elem_end']), r['self_cov'])

# LTR_retriever's hmmsearch order per candidate (join by overlap)
if feat:
    fx = collections.defaultdict(list)
    for f in csv.DictReader(open(feat), delimiter='\t'):
        fx[f['chrom']].append((int(f['l1s']), int(f['l1e']), int(f['l2s']), int(f['l2e']), f['order'], int(f['real'])))
    for v in fx.values(): v.sort()
    fs = {c: [t[0] for t in v] for c, v in fx.items()}
    def ovl(a1, a2, b1, b2):
        o = min(a2, b2) - max(a1, b1) + 1; return o / min(a2 - a1 + 1, b2 - b1 + 1) if o > 0 else 0
    for r in rows:
        v = fx.get(r['seq'], []); k = bisect.bisect_left(fs.get(r['seq'], []), int(r['lltr_start']) - 3000)
        r['lr_order'] = None
        while k < len(v) and v[k][0] <= int(r['lltr_end']) + 3000:
            t = v[k]
            if ovl(int(r['lltr_start']), int(r['lltr_end']), t[0], t[1]) >= 0.5 and ovl(int(r['rltr_start']), int(r['rltr_end']), t[2], t[3]) >= 0.5:
                r['lr_order'] = t[4]; break
            k += 1

real = [r for r in rows if r['real']]; other = [r for r in rows if not r['real']]
print(f"candidates {len(rows)}: real {len(real)} (strict {sum(r['real_s'] for r in rows)}, loose-only {sum(r['real_l'] and not r['real_s'] for r in rows)}), other {len(other)}")

def rate(pred, name):
    tp = sum(1 for r in real if pred(r)); fp = sum(1 for r in other if pred(r))
    print(f"  {name:58s} real {100*tp/max(1,len(real)):5.1f}%  other {100*fp/max(1,len(other)):5.1f}%  precision {100*tp/max(1,tp+fp):5.1f}%  keeps {tp+fp}")

dom_ltr = lambda r: r['dom_order'] == 'LTR'
dom_any = lambda r: r['dom_order'] not in ('none', '.')
veto = lambda r: r['dom_order'] in ('notLTR', 'mixture') or r['self_cov'] >= 0.99 or r['ltr_dom'] > 0
mt = lambda r: r['motif'] == 'TGCA' and r['tsd'] >= 4
mt5 = lambda r: r['motif'] == 'TGCA' and r['tsd'] >= 5
gap = lambda r: r['gap'] is not None and r['gap'] <= 0.10
print("== domain evidence (ltr-forage, in process)")
rate(dom_any, 'domain hit, any order')
rate(dom_ltr, 'domain hit, LTR order')
rate(lambda r: r['dom_order'] in ('notLTR', 'mixture'), 'notLTR or mixture (veto)')
rate(lambda r: r['self_cov'] >= 0.99, 'LTR self-repeat >= 0.99 (veto)')
rate(lambda r: r['ltr_dom'] > 0, 'polyprotein domain inside an LTR (veto)')
rate(lambda r: r['dom_order'] == 'NA', 'hits only outside both lists')
print("== rules")
rate(mt5, 'TG..CA at ends and TSD >= 5')
rate(mt, 'TG..CA at ends and TSD >= 4')
rate(lambda r: dom_ltr(r) or mt(r), 'LTR domain or motif+TSD')
rate(lambda r: not veto(r) and (dom_ltr(r) or mt(r)), 'no veto and (LTR domain or motif+TSD)')
print("== --keep levels (all reject notLTR/mixture)")
rate(lambda r: not veto(r) and dom_ltr(r), 'domain: LTR-order domain hit')
rate(lambda r: not veto(r) and (dom_ltr(r) or mt(r)), 'standard: domain, or TG..CA with TSD >= 4')
rate(lambda r: not veto(r) and (dom_ltr(r) or mt5(r)), '(old standard: domain, or TG..CA with TSD >= 5)')
rate(lambda r: not veto(r) and (dom_ltr(r) or r['tsd'] >= 5 or mt(r)), 'loose: domain, or TSD >= 5, or TG..CA with TSD >= 4')
rate(lambda r: not veto(r), 'none: veto only')
if any(r['gap'] is not None for r in rows):
    rate(gap, 'internal gap <= 10%')
    rate(lambda r: gap(r) and not veto(r) and (dom_ltr(r) or mt(r)), 'gap <= 10%, no veto, (LTR domain or motif+TSD)')
if feat:
    print("== LTR_retriever hmmsearch order on the same candidates (joined by overlap)")
    joined = [r for r in rows if r['lr_order'] is not None]
    print(f"  joined {len(joined)} of {len(rows)}")
    jr = [r for r in joined if r['real']]; jo = [r for r in joined if not r['real']]
    def rate2(pred, name):
        tp = sum(1 for r in jr if pred(r)); fp = sum(1 for r in jo if pred(r))
        print(f"  {name:58s} real {100*tp/max(1,len(jr)):5.1f}%  other {100*fp/max(1,len(jo)):5.1f}%  precision {100*tp/max(1,tp+fp):5.1f}%  keeps {tp+fp}")
    rate2(lambda r: r['lr_order'] == 'LTR', 'LTR_retriever: LTR order')
    rate2(dom_ltr, 'ltr-forage: LTR order')
    rate2(lambda r: r['lr_order'].lower() in ('notltr', 'mixture'), 'LTR_retriever: veto')
    rate2(veto, 'ltr-forage: veto')
    conf = collections.Counter((r['lr_order'], r['dom_order']) for r in joined)
    print("  confusion (LTR_retriever order -> ltr-forage order): count [real]")
    for (a, b), n in sorted(conf.items(), key=lambda x: -x[1])[:12]:
        nr = sum(1 for r in joined if r['lr_order'] == a and r['dom_order'] == b and r['real'])
        print(f"    {a:8s} -> {b:8s} {n:6d} [{nr}]")
