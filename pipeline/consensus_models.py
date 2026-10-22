#!/usr/bin/env python3
"""Family models for the tiered clusters: a te-composer consensus for clusters of >= --min-cons elements,
the cleaned exemplar for smaller ones; clusters flagged as segmental duplications are left out.

usage: consensus_models.py <core.fa> <lib.fa> <clusters.json> <models.fa> [--sd sd_clusters.tsv]
                           [--min-cons 3] [--jobs N] [--te-composer PATH] [--work DIR] [--big-only] [--premade FA]

core.fa: one LTR and one INT record per tiered cluster (the exemplar, after depth cleanup and nested
excision).  lib.fa: every element's raw records.  For a cluster of >= --min-cons elements the LTR model is
the consensus of both LTRs of every member and the INT model the consensus of the internals, from the raw
records: te-composer's insertion packing keeps spans present in most copies and drops the minority
insertions that depth cleanup exists to remove.  Records are oriented to the longest one first (te-composer
before 0.6.6 misplaced minus-strand rows), consensus columns supported by fewer than max(2, half the median)
aligned copies are dropped, at the ends and inside (te-composer keeps a span carried by a single copy: on
lotus a 3.3 kb insertion in 1 of 4 copies cost 7.2 Mb of false positives), and a cluster whose consensus
cannot be built keeps its exemplar.  LTR inputs longer than twice the cluster's median length are left out.  Models keep the exemplar's record names, so
family.py and the reports trace them as before.  Validated on rice, 2026-10-19 (experiment 2, arm T): 273
curated families recovered in full against 260 for exemplars, precision 98.26 vs 98.24.  With the column
gate (arm G): rice 273, precision 98.27; lotus 250 vs 225 for exemplars, precision 92.50 vs 92.36 (ends-only
trimming: 83.90).
"""
import argparse, collections, json, os, re, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def read_fa(path):
    seqs, order, n = {}, [], None
    for l in open(path):
        if l.startswith('>'):
            n = l[1:].split()[0].split('#')[0]; seqs[n] = []; order.append(n)
        else:
            seqs[n].append(l.strip())
    return {k: ''.join(seqs[k]) for k in order}


def orient(names, seqs, k=12):
    """Turn each record to share more k-mers with the longest record."""
    km = lambda x: {x[i:i + k] for i in range(len(x) - k + 1)}
    ref = km(max((seqs[n].upper() for n in names), key=len))
    out = []
    for n in names:
        x = seqs[n].upper(); r = x[::-1].translate(RC)
        out.append((n, x if len(km(x) & ref) >= len(km(r) & ref) else r))
    return out


def gate(stk, seq):
    """Drop every consensus column (ends and interior) where fewer than max(2, half the median) copies align
    (Stockholm RF columns).  te-composer keeps a span carried by a single copy, at the ends or inside."""
    rows = collections.defaultdict(str); rf = ''
    for l in open(stk):
        if l.startswith('#=GC RF'): rf += l.split()[-1]
        elif l.strip() and not l.startswith(('#', '//')):
            n, x = l.split()[:2]; rows[n] += x
    cols = [i for i, c in enumerate(rf) if c not in '.-']
    if len(cols) != len(seq): return seq
    occ = [sum(1 for x in rows.values() if x[i] not in '.-') for i in cols]
    t = max(2, 0.5 * sorted(occ)[len(occ) // 2])
    out = ''.join(b for b, o in zip(seq, occ) if o >= t)
    return out if len(out) >= 50 else seq


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('core'); ap.add_argument('lib'); ap.add_argument('clusters'); ap.add_argument('models')
    ap.add_argument('--sd', help='flank_filter.py output; flagged clusters are left out')
    ap.add_argument('--min-cons', type=int, default=3)
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--te-composer', required=True)
    ap.add_argument('--work', default='cons')
    ap.add_argument('--big-only', action='store_true',
                    help='write only clusters of >= --min-cons elements (stage B representatives)')
    ap.add_argument('--premade', help='models already built (stage B reps.fa): core records of the same name get no new consensus')
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)

    core, raw = read_fa(a.core), read_fa(a.lib)
    clus = json.load(open(a.clusters))
    members = collections.defaultdict(list)
    for e, c in clus.items(): members[c].append(e)
    by_seq = collections.defaultdict(list)
    for e in clus:
        m = re.match(r'(\S+):(\d+)\.\.(\d+)$', e); by_seq[m[1]].append((int(m[2]), int(m[3]), e))

    def elem_of(rec):
        m = re.match(r'(\S+):(\d+)\.\.(\d+)_(LTR|INT)$', rec)
        for s, t, e in by_seq.get(m[1], ()):
            if s <= int(m[2]) and int(m[3]) <= t: return e

    recs = collections.defaultdict(lambda: {'LTR': [], 'INT': []})
    for r in raw:
        e = elem_of(r)
        if e: recs[e][r[-3:]].append(r)
    sd = set()
    if a.sd:
        for l in open(a.sd):
            f = l.rstrip('\n').split('\t')
            if f[0] != 'cluster' and f[-1] == 'sd': sd.add(f[0])

    core_recs = collections.defaultdict(dict)
    for r in core: core_recs[clus[elem_of(r)]][r[-3:]] = r

    def inputs(c, kind):
        names = [r for e in members[c] for r in recs[e][kind] if len(raw[r]) >= 50]
        if kind == 'LTR' and len(names) >= 3:
            # an LTR more than twice the cluster's median length is a lookalike (a duplicated segment called as an
            # element joins a real family's cluster); te-composer can pick it as the reference, and its extra span
            # then reaches the consensus: zebra finch 904 -> 4,433 bp (+10.8 Mb false positives), gorilla 439 ->
            # 1,229 bp (+8.3 Mb)
            med = sorted(len(raw[r]) for r in names)[len(names) // 2]
            names = [r for r in names if len(raw[r]) <= 2 * med]
        return names

    def consensus(c, kind):
        names = inputs(c, kind)
        if len(names) < 2: return None
        tag = re.sub(r'[^\w.-]', '_', f'c{c}_{kind}')
        fa, stk, out = (os.path.join(a.work, tag + x) for x in ('.fa', '.stk', '.cons.fa'))
        text = ''.join(f'>{n}\n{x}\n' for n, x in orient(names, raw))
        # reuse a finished consensus built from the same input (reruns with --from models)
        cached = all(os.path.exists(f) for f in (fa, stk, out)) and open(fa).read() == text \
            and os.path.getmtime(out) >= os.path.getmtime(fa)
        if not cached:
            for f in (stk, out):
                if os.path.exists(f): os.remove(f)
            with open(fa, 'w') as fh: fh.write(text)
            # a few large clusters dominate the run time (rice: one 247-copy INT took 1h54m on 2 threads), so
            # they get more threads; output does not depend on the thread count
            th = 8 if sum(len(raw[n]) for n in names) > 500_000 else 2
            r = subprocess.run([a.te_composer, '--silent', '--threads', str(th), '--name', tag, '--consensus', out,
                                '--format', 'fasta', fa, stk], stderr=open(os.path.join(a.work, tag + '.log'), 'w'))
            if r.returncode: return None
        if not os.path.exists(out): return None   # e.g. nothing aligns above te-composer's minimum score
        s = ''.join(l.strip() for l in open(out) if not l.startswith('>'))
        return gate(stk, s) if s else None

    premade = read_fa(a.premade) if a.premade else {}
    jobs = [(c, k) for c in core_recs if c not in sd and len(members[c]) >= a.min_cons for k in ('LTR', 'INT')
            if k in core_recs[c] and core_recs[c][k] not in premade]
    jobs.sort(key=lambda j: -sum(len(raw[n]) for n in inputs(*j)))   # largest first, so they don't finish last
    with ThreadPoolExecutor(a.jobs) as ex:
        cons = dict(zip(jobs, ex.map(lambda j: consensus(*j), jobs)))

    st = collections.Counter()
    with open(a.models, 'w') as fh:
        for c, cr in core_recs.items():
            if c in sd: st['segmental-duplication clusters left out'] += 1; continue
            if a.big_only and len(members[c]) < a.min_cons: continue
            for kind, name in sorted(cr.items()):
                # a stage B representative's core record is its consensus after depth, nested and exemplar
                # excision; premade only marks it as built (its raw consensus can still hold a cut insertion)
                s = cons.get((c, kind))
                if name in premade: st[f'{kind} consensus (stage B representative)'] += 1
                elif s: st[f'{kind} consensus'] += 1
                elif (c, kind) in cons: st[f'{kind} consensus failed, exemplar kept'] += 1
                else: st[f'{kind} exemplar (cluster < {a.min_cons})'] += 1
                fh.write(f'>{name}\n{s or core[name]}\n')
    for k, v in sorted(st.items()): print(f'{k}: {v}')


if __name__ == '__main__':
    main()
