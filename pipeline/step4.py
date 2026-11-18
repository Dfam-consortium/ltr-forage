# Two-phase element clustering. Phase 1: link elements whose internals cover each other (80% of both)
# and whose LTRs cover the shorter LTR. Phase 2: attach each remaining singleton element to the
# cluster whose members cover the singleton's own internal and LTR records; attachments never
# merge clusters.
# usage: step4.py <int_hits.tsv> <ltr_hits.tsv> [clusters.json]   (env TSV: the ltr-forage table, cand.tsv)
# The functions are shared with greedy_cluster.py (step 12 on very large raw libraries).
import os, csv, collections, sys, json

def load_elements(tsv):
    """Kept elements in table order and their records: L (left LTR), I (internal), R (right LTR)."""
    elem_recs = {}; elems = []
    for r in csv.DictReader(open(tsv), delimiter='\t'):
        if r['keep'] != '1': continue
        s = r['seq']; e = (s, int(r['elem_start']), int(r['elem_end'])); elems.append(e)
        ls, le, rs, re_ = int(r['lltr_start']), int(r['lltr_end']), int(r['rltr_start']), int(r['rltr_end'])
        elem_recs[e] = {'L': f"{s}:{ls}..{le}_LTR", 'I': f"{s}:{le+1}..{rs-1}_INT", 'R': f"{s}:{rs}..{re_}_LTR"}
    return elems, elem_recs

def union(iv):
    iv = sorted(iv); tot = 0; cs = ce = None
    for s, e in iv:
        if cs is None: cs, ce = s, e
        elif s <= ce + 1: ce = max(ce, e)
        else: tot += ce - cs + 1; cs, ce = s, e
    if cs is not None: tot += ce - cs + 1
    return tot

def add_coverage(cov, path, kind, rec2elem):
    """directed coverage: cov[(a,b)][kind] = (fraction of a's record covered by b's, of b's by a's, bits);
    the record pair of the kind that covers most of a's record is kept (phase 2 asks for that).
    cov[(a,b)][kind + 'max'] is the highest coverage of either record over every record pair of the two
    elements, the same for (a,b) and (b,a): phase 1's LTR test, which is of the pair and has no direction."""
    pair = collections.defaultdict(lambda: {'q': [], 's': [], 'ql': 0, 'sl': 0, 'bits': 0.0})
    for l in open(path):
        q, s, pid, ln, qs, qe, ss, se, ql, sl, bits = l.split('\t'); q = q.split('#')[0]; s = s.split('#')[0]
        if q == s or float(pid) < 80: continue
        if not (q.endswith('_INT') == (kind == 'I') and s.endswith('_INT') == (kind == 'I')): continue
        p = pair[(q, s)]; p['q'].append((min(int(qs), int(qe)), max(int(qs), int(qe)))); p['s'].append((min(int(ss), int(se)), max(int(ss), int(se)))); p['ql'] = int(ql); p['sl'] = int(sl); p['bits'] += float(bits)
    for (q, s), p in pair.items():
        if q not in rec2elem or s not in rec2elem: continue
        eq, es = rec2elem[q][0], rec2elem[s][0]
        if eq == es: continue
        cq = union(p['q']) / p['ql']; cs = union(p['s']) / p['sl']
        for a, b, ca, cb in ((eq, es, cq, cs), (es, eq, cs, cq)):
            d = cov[(a, b)]; old = d.get(kind, (0, 0, 0))
            if ca > old[0] or (ca == old[0] and p['bits'] > old[2]): d[kind] = (ca, cb, p['bits'])
            d[kind + 'max'] = max(d.get(kind + 'max', 0), cq, cs)

def linked(d):
    """phase 1: internals cover each other at 80%, and some LTR of one covers 80% of an LTR of the other.
    The same answer for (a,b) and (b,a): an element has one internal record, and the LTR test takes the
    best of the four LTR pairs.  (It was read off the pair covering most of a's LTR, so two elements whose
    LTR records differ in length could link in one direction only, and the clustering tests one.)"""
    if 'I' not in d or 'L' not in d: return False
    ci, cj, _ = d['I']
    return ci >= 0.8 and cj >= 0.8 and d['Lmax'] >= 0.8

def covered(d):
    """phase 2: a's own internal and LTR are each covered at 80% by b's; the score of the cover"""
    if 'I' in d and 'L' in d and d['I'][0] >= 0.8 and d['L'][0] >= 0.8: return d['I'][2] + d['L'][2]
    return None

def report(label, elems, member):
    cores = collections.defaultdict(list)
    for e in elems: cores[member(e)].append(e)
    print(f"{label}: {len(cores)} clusters ({sum(1 for v in cores.values() if len(v)==1)} singletons, largest {max((len(v) for v in cores.values()), default=0)})")
    return cores

def main():
    int_hits, ltr_hits = sys.argv[1], sys.argv[2]
    if 'TSV' not in os.environ: sys.exit('step4.py: set TSV to the ltr-forage table (cand.tsv)')
    elems, elem_recs = load_elements(os.environ['TSV'])
    rec2elem = {n: (e, k) for e, d in elem_recs.items() for k, n in d.items()}
    cov = collections.defaultdict(dict)
    for path, kind in ((int_hits, 'I'), (ltr_hits, 'L')): add_coverage(cov, path, kind, rec2elem)
    parent = {e: e for e in elems}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for (a, b), d in cov.items():
        if a < b and linked(d):
            ra, rb = find(a), find(b)
            if ra != rb: parent[ra] = rb
    cores = report('phase 1 (INT both 80%, LTR shorter 80%)', elems, find)
    # phase 2: attach singletons
    attach = {}
    cand = collections.defaultdict(lambda: collections.defaultdict(float))
    for (a, b), d in cov.items():
        if len(cores[find(a)]) != 1: continue
        rb = find(b)
        if len(cores[rb]) == 1: continue
        sc = covered(d)
        if sc is not None: cand[a][rb] = max(cand[a][rb], sc)
    for a, c in cand.items(): attach[a] = max(c.items(), key=lambda x: x[1])[0]
    print(f"phase 2 attached {len(attach)} of {sum(1 for v in cores.values() if len(v)==1)} singletons")
    report('phase 1 + phase 2 attachment', elems, lambda e: attach.get(e, find(e)))
    if len(sys.argv) > 3:
        json.dump({f"{e[0]}:{e[1]}..{e[2]}": f"{r[0]}:{r[1]}..{r[2]}" for e in elems for r in [attach.get(e, find(e))]},
                  open(sys.argv[3], 'w'))

if __name__ == '__main__':
    main()
