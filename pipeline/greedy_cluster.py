#!/usr/bin/env python3
"""Step 12 for very large raw libraries: greedy element clustering with step4.py's rules (ltr_library.py).
usage: greedy_cluster.py probe <lib.fa> [--n N] [--seed S]
       greedy_cluster.py cluster <lib.fa> <clusters.json> [--batch B] [--reps C] [--seed S] [--work DIR]
Env as rmb_ava.sh (THREADS, RMBLASTN, WORD), plus TSV (the ltr-forage table, cand.tsv).

probe: searches N random records of each kind (LTR, internal) against all records of the kind and prints the
projected number of >= 80% hits of the all-vs-all search, its cost.  Step 12's all-vs-all costs about
1.5-2.5 min per million hits at 16 threads (mm10: 45 M hits, 60-110 min; maize would need 2.4 G).

cluster: elements in a seeded random order.  The first B unassigned elements are clustered all-vs-all;
every unassigned element is then searched against the new representatives only, at most C members per
cluster, and joins a cluster whose representative it links to by step4.py's phase-1 rule (clusters it links
join each other).  Joined elements become representatives while their cluster has fewer than C, and the
unassigned are searched against them too; when nothing more joins, the next B form the next batch.  Cost then
grows with elements x clusters, not elements^2.  At the end, singletons are attached by step4.py's phase-2
rule to the cluster that covers them best.  Links between two non-representative members are not seen, so
with C below a cluster's size two clusters a member would bridge can stay apart.

Memory: a round's hits are split on disk by query element, at most about --bucket-hits per bucket, and read one
bucket at a time; a pair's hits all come from one query element in those searches.  Phase-2 scores are kept per
element and cluster, and dropped when the element joins a cluster."""
import os, sys, re, json, math, random, argparse, subprocess, collections, shutil
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import step4

AVA = 'qseqid sseqid pident length qstart qend sstart send qlen slen bitscore'.split()

def read_fasta(p):
    name = None; seq = []
    for l in open(p):
        if l.startswith('>'):
            if name: yield name, ''.join(seq)
            name = l[1:].strip(); seq = []
        else: seq.append(l.strip())
    if name: yield name, ''.join(seq)

class Search:
    def __init__(self, lib, work):
        self.work = work; os.makedirs(work, exist_ok=True)
        self.recs = {n.split()[0].split('#')[0]: (n, s) for n, s in read_fasta(lib)}
        self.env = dict(os.environ, MIN_PID='80')
        if 'MATRIX' not in self.env:   # the whole library's matrix, as the all-vs-all search picks
            self.env['MATRIX'] = subprocess.run(['bash', os.path.join(HERE, 'pick_matrix.sh'), os.path.abspath(lib)],
                                                capture_output=True, text=True, check=True).stdout.strip()
        self.n = 0; self.hits = 0

    def write(self, path, names):
        k = 0
        with open(path, 'w') as fh:
            for r in names:
                if r in self.recs: n, s = self.recs[r]; fh.write(f'>{n}\n{s}\n'); k += 1
        return k

    def run(self, q, db, kind):
        """q, db: record names of one kind; returns the hit file"""
        self.n += 1; b = os.path.join(self.work, f's{self.n}.{kind}')
        if not self.write(b + '.q.fa', q) or not self.write(b + '.db.fa', db):
            open(b + '.tsv', 'w').close(); return b + '.tsv'
        with open(os.path.join(self.work, 'search.log'), 'a') as log:
            subprocess.run(['bash', os.path.join(HERE, 'rmb_ava.sh'), b + '.q.fa', b + '.db.fa', '5', b + '.tsv', *AVA],
                           env=self.env, stdout=log, stderr=log, check=True)
        with open(b + '.tsv') as fh: self.hits += sum(1 for _ in fh)
        return b + '.tsv'

    def clean(self):
        for f in os.listdir(self.work):
            if f.startswith('s') and (f.endswith('.fa') or f.endswith('.tsv')): os.remove(os.path.join(self.work, f))

def recs_of(elems, elem_recs, kind):
    return [elem_recs[e][k] for e in elems for k in (('I',) if kind == 'I' else ('L', 'R'))]

def probe(a):
    elems, elem_recs = step4.load_elements(os.environ['TSV'])
    S = Search(a.lib, a.work)
    rng = random.Random(a.seed); total = 0
    for kind in ('L', 'I'):
        allr = [r for r in recs_of(elems, elem_recs, kind) if r in S.recs]
        if not allr: continue
        q = rng.sample(allr, min(a.n, len(allr)))
        with open(S.run(q, allr, kind)) as fh: h = sum(1 for _ in fh)
        proj = h / len(q) * len(allr); total += proj
        print(f'{"LTR" if kind == "L" else "INT"}: {len(allr)} records, {h / len(q):.0f} hits per query, {proj / 1e6:.1f} M projected')
    S.clean()
    print(f'projected all-vs-all hits: {total / 1e6:.1f} M')

def cluster(a):
    elems, elem_recs = step4.load_elements(os.environ['TSV'])
    rec2elem = {n: (e, k) for e, d in elem_recs.items() for k, n in d.items()}
    S = Search(a.lib, a.work)
    order = [e for e in elems if any(r in S.recs for r in elem_recs[e].values())]
    random.Random(a.seed).shuffle(order)
    parent = {e: e for e in elems}; nrep = collections.Counter(); size = {e: 1 for e in elems}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def join(x, y):
        rx, ry = find(x), find(y)
        if rx != ry: parent[rx] = ry; nrep[ry] += nrep.pop(rx, 0); size[ry] += size.pop(rx)
    cover = collections.defaultdict(dict)   # phase 2, singletons only: cover[a][cluster] = best score of a member covering a's records
    def score(cov):
        for (x, y), d in cov.items():
            if size[find(x)] > 1: continue
            sc = step4.covered(d); r = find(y)
            if sc is not None and sc > cover[x].get(r, 0): cover[x][r] = sc
    def prune():
        for x in [x for x in cover if size[find(x)] > 1]: del cover[x]
    def files(q, db): return [(S.run(recs_of(q, elem_recs, k), recs_of(db, elem_recs, k), k), k) for k in ('I', 'L')]
    def read(fs):
        cov = collections.defaultdict(dict)
        for p, k in fs: step4.add_coverage(cov, p, k, rec2elem)
        return cov
    def buckets(fs):
        """the hits of files(), one bucket's coverage at a time; q and db are disjoint"""
        n = 0
        for p, _ in fs:
            with open(p) as fh: n += sum(1 for _ in fh)
        nb = math.ceil(n / a.bucket_hits)
        if nb <= 1:
            yield read(fs); return
        idx = {e: i for i, e in enumerate(order)}; out = []
        for p, k in fs:
            bf = [open(f'{p}.{b}.tsv', 'w') for b in range(nb)]
            with open(p) as fh:
                for l in fh:
                    q = l[:l.index('\t')].split('#')[0]
                    if q in rec2elem: bf[idx[rec2elem[q][0]] % nb].write(l)
            for f in bf: f.close()
            os.remove(p); out.append((p, k))
        for b in range(nb):
            bfs = [(f'{p}.{b}.tsv', k) for p, k in out]
            yield read(bfs)
            for p, _ in bfs: os.remove(p)
    def make_reps(cands):
        out = []
        for e in cands:
            r = find(e)
            if nrep[r] < a.reps: nrep[r] += 1; out.append(e)
        return out

    U = order; rnd = 0; nb = 0
    reps = []
    while U:
        # every unassigned element against the newest representatives, until none joins
        while reps and U:
            rnd += 1; Us = set(U); rs = set(reps)
            joined = set()
            for cov in buckets(files(U, reps)):
                score(cov)
                for (x, y), d in cov.items():
                    if x in Us and y in rs and step4.linked(d): join(x, y); joined.add(x)
            prune()
            U = [e for e in U if e not in joined]
            reps = make_reps([e for e in order if e in joined])
            print(f'round {rnd}: {len(joined)} joined, {len(U)} unassigned, {len(reps)} new representatives', flush=True)
            S.clean()
        if not U: break
        batch, U = U[:a.batch], U[a.batch:]
        nb += 1; rnd += 1
        cov = read(files(batch, batch)); score(cov)
        for (x, y), d in cov.items():
            if x < y and step4.linked(d): join(x, y)
        del cov; prune()
        reps = make_reps(batch)
        print(f'round {rnd}: batch {nb} of {len(batch)} elements, {len(set(find(e) for e in batch))} clusters, '
              f'{len(U)} unassigned', flush=True)
        S.clean()
    print(f'{rnd} rounds, {nb} batches, {S.hits / 1e6:.1f} M hits searched')
    cores = step4.report('phase 1 (INT both 80%, LTR shorter 80%)', elems, find)
    attach = {}
    for x, c in cover.items():
        if len(cores[find(x)]) != 1: continue
        best = collections.defaultdict(float)
        for r, sc in c.items():
            ry = find(r)
            if len(cores[ry]) > 1: best[ry] = max(best[ry], sc)
        if best: attach[x] = max(best.items(), key=lambda t: t[1])[0]
    print(f"phase 2 attached {len(attach)} of {sum(1 for v in cores.values() if len(v)==1)} singletons")
    step4.report('phase 1 + phase 2 attachment', elems, lambda e: attach.get(e, find(e)))
    json.dump({f"{e[0]}:{e[1]}..{e[2]}": f"{r[0]}:{r[1]}..{r[2]}" for e in elems for r in [attach.get(e, find(e))]},
              open(a.out, 'w'))
    shutil.rmtree(a.work, ignore_errors=True)

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('probe'); p.add_argument('lib'); p.add_argument('--n', type=int, default=150)
    c = sub.add_parser('cluster'); c.add_argument('lib'); c.add_argument('out')
    c.add_argument('--batch', type=int, default=2000); c.add_argument('--reps', type=int, default=100)
    c.add_argument('--bucket-hits', type=int, default=2000000)
    for x in (p, c):
        x.add_argument('--seed', type=int, default=1); x.add_argument('--work', default='greedy_work')
    a = ap.parse_args()
    if 'TSV' not in os.environ: sys.exit('greedy_cluster.py: set TSV to the ltr-forage table (cand.tsv)')
    probe(a) if a.cmd == 'probe' else cluster(a)
