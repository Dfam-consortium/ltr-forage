#!/usr/bin/env python3
"""Experiment 2 (2026-10-19): cleanup arms for the same element clusters.

usage: build_arms.py <run_dir> <out_dir> [--jobs 24] [--min-cons 3] [--te-composer PATH]

<run_dir> is a finished ltr_library.py run (lib.fa, clean/clean.fa, core.fa, elem/clusters_nested.json).
The clusters and the set of clusters in core.fa (the tiered library) are held fixed; only the records that
stand for each cluster change.  Each arm writes <out_dir>/<arm>.fa with one LTR and one INT record per core
cluster, named after the cluster's exemplar records so classify2.py and family.py can trace them:

  E  exemplar     core.fa as is (depth cleanup + nested vote + excision; today's recipe)
  A  raw-cons     clusters >= --min-cons: te-composer consensus of the members' raw records (lib.fa);
                  smaller clusters: the exemplar's raw records (no cleanup at all)
  B  clean-cons   clusters >= --min-cons: te-composer consensus of the members' depth-cleaned records
                  (clean/clean.fa); smaller clusters: core.fa records
  D  split        clusters >= --min-cons: as A; smaller clusters: core.fa records (the proposed design)
  T  split+trim   as D, with consensus ends trimmed where fewer than max(2, half the median) copies align
  G  split+gate   as D, with every consensus column dropped where fewer than max(2, half the median) copies
                  align (lotus, 2026-10-19: a 3.3 kb insertion carried by 1 of 4 copies sat inside the
                  consensus, so end trimming missed it; 7.2 Mb of false positives)

LTR consensi are built from both LTRs of every member, INT consensi from the internals.  te-composer runs
without --genome, so no extension.
"""
import argparse, collections, json, os, re, subprocess
from concurrent.futures import ThreadPoolExecutor

ap = argparse.ArgumentParser()
ap.add_argument('run'); ap.add_argument('out')
ap.add_argument('--jobs', type=int, default=24)
ap.add_argument('--min-cons', type=int, default=3)
ap.add_argument('--cons-dir', default='cons', help='consensus work directory under <out_dir>')
ap.add_argument('--no-orient', action='store_true', help='pass records to te-composer in genome orientation')
ap.add_argument('--te-composer', default=os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),   # dfam-curator beside LTRHarvest
                 '..', '..', '..', '..', '..', 'dfam-curator', 'target', 'release', 'te-composer')))
a = ap.parse_args()
os.makedirs(os.path.join(a.out, a.cons_dir), exist_ok=True)


def read_fa(path):
    seqs, n = {}, None
    for l in open(path):
        if l.startswith('>'):
            n = l[1:].split()[0].split('#')[0]; seqs[n] = []
        else:
            seqs[n].append(l.strip())
    return {k: ''.join(v) for k, v in seqs.items()}


raw = read_fa(os.path.join(a.run, 'lib.fa'))
clean = read_fa(os.path.join(a.run, 'clean', 'clean.fa'))
core = read_fa(os.path.join(a.run, 'core.fa'))
clus = json.load(open(os.path.join(a.run, 'elem', 'clusters_nested.json')))
members = collections.defaultdict(list)
for e, c in clus.items(): members[c].append(e)

# records of each element, by coordinate containment (record names are <seq>:<s>..<e>_LTR|_INT)
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

# core clusters and their exemplar records (one LTR + one INT each)
core_recs = collections.defaultdict(dict)
for r in core:
    core_recs[clus[elem_of(r)]][r[-3:]] = r
print(f'{len(core_recs)} core clusters; sizes', collections.Counter(min(len(members[c]), 10) for c in core_recs))


RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def orient(names, seqs, k=12):
    """Records are in genome orientation, so members on the minus strand come reverse-complemented.
    te-composer does not reconcile mixed orientations (2026-10-19: 160 of 258 consensi were > 2% N, all but 4
    from clusters with minus-strand members; orienting first took a test cluster from 18% N to 0%). Each record is
    turned to share more k-mers with the longest record."""
    if a.no_orient: return [(n, seqs[n]) for n in names]
    km = lambda x: {x[i:i + k] for i in range(len(x) - k + 1)}
    ref = km(max((seqs[n].upper() for n in names), key=len))
    out = []
    for n in names:
        x = seqs[n].upper(); r = x[::-1].translate(RC)
        out.append((n, x if len(km(x) & ref) >= len(km(r) & ref) else r))
    return out


def consensus(tag, kind, names, seqs):
    """te-composer consensus of the given records; None when fewer than 2 usable sequences."""
    fa = os.path.join(a.out, a.cons_dir, f'{tag}_{kind}.fa')
    use = [n for n in names if n in seqs and len(seqs[n]) >= 50]
    if len(use) < 2: return None
    with open(fa, 'w') as fh:
        for n, x in orient(use, seqs): fh.write(f'>{n}\n{x}\n')
    out = fa[:-3] + '.cons.fa'
    if not os.path.exists(out):
        r = subprocess.run([a.te_composer, '--silent', '--threads', '2', '--name', f'{tag}_{kind}', '--consensus', out,
                            '--format', 'fasta', fa, fa[:-3] + '.stk'], stderr=open(fa[:-3] + '.log', 'w'))
        # te-composer exits 1 when nothing aligns above its minimum score (e.g. short, diverged LTRs);
        # the cluster then keeps its exemplar record.
        if r.returncode != 0: return None
    s = read_fa(out)
    return next(iter(s.values())) if s else None


def occupancy(tag, kind, seq):
    """Copies aligned at each consensus column (Stockholm RF line); None if the columns do not match seq."""
    rows = collections.defaultdict(str); rf = ''
    for l in open(os.path.join(a.out, a.cons_dir, f'{tag}_{kind}.stk')):
        if l.startswith('#=GC RF'): rf += l.split()[-1]
        elif l.strip() and not l.startswith(('#', '//')):
            n, x = l.split()[:2]; rows[n] += x
    cols = [i for i, c in enumerate(rf) if c not in '.-']
    if len(cols) != len(seq): return None
    return [sum(1 for x in rows.values() if x[i] not in '.-') for i in cols]


def gate(tag, kind, seq):
    """Drop every consensus column where fewer than max(2, half the median) copies align."""
    if seq is None: return None
    occ = occupancy(tag, kind, seq)
    if occ is None: return seq
    t = max(2, 0.5 * sorted(occ)[len(occ) // 2])
    out = ''.join(b for b, o in zip(seq, occ) if o >= t)
    return out if len(out) >= 50 else seq


def trim_ends(tag, kind, seq):
    """Trim consensus ends supported by few copies: te-composer's insertion packing drops minority insertions
    inside the alignment but keeps a terminal overhang even when a single copy supplies it (2026-10-19: a CACTA
    + LTR tail from 1 of 15 copies added 1.6 Mb of rice false positives). Occupancy per consensus column comes
    from the Stockholm RF line; ends below max(2, half the median occupancy) are cut."""
    if seq is None: return None
    rows = collections.defaultdict(str); rf = ''
    for l in open(os.path.join(a.out, a.cons_dir, f'{tag}_{kind}.stk')):
        if l.startswith('#=GC RF'): rf += l.split()[-1]
        elif l.strip() and not l.startswith(('#', '//')):
            n, x = l.split()[:2]; rows[n] += x
    cols = [i for i, c in enumerate(rf) if c not in '.-']
    if len(cols) != len(seq): return seq
    occ = [sum(1 for x in rows.values() if x[i] not in '.-') for i in cols]
    t = max(2, 0.5 * sorted(occ)[len(occ) // 2])
    i, j = 0, len(occ)
    while i < j and occ[i] < t: i += 1
    while j > i and occ[j - 1] < t: j -= 1
    return seq[i:j] if j - i >= 50 else seq


jobs = []
for c, cr in core_recs.items():
    if len(members[c]) < a.min_cons: continue
    for kind in ('LTR', 'INT'):
        names = [r for e in members[c] for r in recs[e][kind]]
        jobs.append((c, kind, 'raw', names, raw)); jobs.append((c, kind, 'clean', names, clean))
print(f'{len(jobs)} te-composer runs')
with ThreadPoolExecutor(a.jobs) as ex:
    res = list(ex.map(lambda j: consensus(f'c{j[0]}_{j[2]}', j[1], j[3], j[4]), jobs))
cons = {(j[0], j[1], j[2]): s for j, s in zip(jobs, res)}
trim = {(j[0], j[1], j[2]): trim_ends(f'c{j[0]}_{j[2]}', j[1], s) for j, s in zip(jobs, res)}
gated = {(j[0], j[1], j[2]): gate(f'c{j[0]}_{j[2]}', j[1], s) for j, s in zip(jobs, res)}
gcut = [(k, len(cons[k]) - len(v)) for k, v in gated.items() if v and cons[k] and len(v) < len(cons[k])]
print(f'column gate shortened {len(gcut)} consensus sequences by {sum(x for _, x in gcut):,} bp')
cut = [(k, len(cons[k]) - len(v)) for k, v in trim.items() if v and cons[k] and len(v) < len(cons[k])]
print(f'end trimming shortened {len(cut)} consensus sequences by {sum(x for _, x in cut):,} bp')

arms = {k: open(os.path.join(a.out, f'{k}.fa'), 'w') for k in 'EABDTG'}
stats = collections.Counter()
for c, cr in sorted(core_recs.items()):
    big = len(members[c]) >= a.min_cons
    for kind, name in sorted(cr.items()):
        pick = {'E': core[name], 'A': raw.get(name, core[name]), 'B': core[name], 'D': core[name], 'T': core[name], 'G': core[name]}
        if big:
            for arm, src in (('A', 'raw'), ('B', 'clean'), ('D', 'raw')):
                s = cons.get((c, kind, src))
                if s: pick[arm] = s; stats[f'{arm} consensus {kind}'] += 1
                else: stats[f'{arm} no consensus {kind} (kept record)'] += 1
            if trim.get((c, kind, 'raw')): pick['T'] = trim[(c, kind, 'raw')]
            if gated.get((c, kind, 'raw')): pick['G'] = gated[(c, kind, 'raw')]
        for arm, s in pick.items():
            arms[arm].write(f'>{name}\n{s}\n')
for f in arms.values(): f.close()
for k, v in sorted(stats.items()): print(f'{k}: {v}')
