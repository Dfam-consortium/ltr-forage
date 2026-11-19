#!/usr/bin/env python3
"""Veto kept elements whose LTRs are copies of the run's poly-A-terminated repeats (ltr_library.py step 11b).
usage: polya_family.py <outdir> <genome.fa> [--min-depth 10] [--min-cov 0.5] [--sample 2000] [--min-pool 2000]
                       [--seed 1] [--no-rescue] [--off]

ltr-forage vetoes a candidate whose two LTRs both end in a poly-A tail (veto:polyA): a pair of SINE or LINE
copies, not an LTR element.  The same repeats also flank old ERV fragments, and such a pair is kept on the
domain hit between them even though its "LTRs" are the Alu, B1 or L1 copies that thousands of vetoed pairs
are made of.  This step compares the LTRs of every kept element with the left LTRs of up to --sample
poly-A-vetoed candidates (rmb_ava.sh, word 12) and vetoes an element when at least --min-cov of each of its
LTRs is covered by at least --min-depth of the sampled LTRs.

The pool the sample is drawn from:
  - A poly-A-vetoed candidate with an LTR domain hit of its own (dom_order LTR) is left out.  The poly-A rule
    also vetoes real LTR elements whose LTRs end A-rich, and their family's kept members match them.
  - The step does nothing when the pool is smaller than --min-pool.  The pairs it removes exist where the
    genome has tens of thousands of poly-A-tailed repeat pairs (gorilla 43,000, mouse 18,000 vetoed
    candidates); where there are a few dozen (Arabidopsis 1, zebra finch 3, rice 12, lotus 42; held out: worm 10,
    fly 12, maize 24, zebrafish 78, soybean 80) a large share of them can be real elements (maize 18 of 24,
    soybean 32 of 80), and 10 matches among them is no evidence of a SINE or LINE.  At the default the sample
    is always full, so --min-depth is a count out of --sample.
    NOT MEASURED: a pool between 80 and 17,263 (before the domain-hit candidates are left out).  The default
    is a cliff inside that range; a genome with a moderate SINE or LINE load can land near it, and what the
    step does just above it, or what such a genome loses just below it, is unknown.  The summary line always
    gives the pool size.

The sample is drawn with --seed (the driver's --seed) from the pool in cand.tsv order, and the search result does
not depend on the thread count: a rerun, and a run with 3 threads instead of 8, give the same cand.tsv and lib.fa.
Another seed moves the vetoed set by a few elements whose match count is near --min-depth (gorilla L1 pairs match
10 to 32 of 2,000).

Overlap rescue (--no-rescue turns it off): ltr-forage resolves overlaps before this step, so an element that
lost the overlap choice to a pair vetoed here would stay dropped.  A candidate dropped as overlapping
(evidence "overlap") that overlaps a vetoed element and no element still kept is searched against the same
sample and kept again when no part of either LTR is covered by --min-depth of the sampled LTRs.  That is
stricter than the veto: a candidate that shares ground with a vetoed pair and has the same repeats in its own
LTRs most likely has its boundaries from them (on gorilla and mm10 the LTR-annotated candidates are at 0 on
both LTRs, 114 of 125; the SINE- and LINE-flanked ones are not, 32 of 38).  Among candidates that overlap each
other the one with the most evidence is taken, in ltr-forage's order (domain, TSD, TG..CA and PBS count, then
LTR identity).  Its evidence reads "rescued" in place of "overlap", and its three records are added to lib.fa
from the genome.

Developed on gorilla and mm10 (2026-11-16): no kept element with LTR-annotated LTRs matched more than 5 of
2,000 sampled LTRs; SINE-flanked ones matched ~1,940 / ~1,420, LINE-flanked ones ~25 / ~150.
Run on hg38 afterwards (2026-11-18, pool 47,632): 348 vetoed, none LTR-annotated (302 SINE-, 43 LINE-flanked,
3 Alu and L1 pieces); LTR-annotated elements matched at most 1.  Human shares its SINE and LINE families with
gorilla; a mammal with other SINE families has not been run.

cand.tsv is rewritten with keep 0 and the evidence tag veto:polyAfamily on the vetoed rows (a rerun first takes
the tag off again and puts rescued rows back to "overlap", so the step can be redone or turned off with --off);
lib.fa is rewritten from lib.raw.fa without their records.  polya_family/ holds the sequences, the hits,
vetoed.tsv and rescued.tsv."""
import sys, os, csv, argparse, bisect, random, subprocess, collections
import numpy as np
ap = argparse.ArgumentParser()
ap.add_argument('outdir'); ap.add_argument('genome')
ap.add_argument('--min-depth', type=int, default=10); ap.add_argument('--min-cov', type=float, default=0.5)
ap.add_argument('--sample', type=int, default=2000); ap.add_argument('--min-pool', type=int, default=2000)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--no-rescue', action='store_true'); ap.add_argument('--off', action='store_true')
a = ap.parse_args()
HERE = os.path.dirname(os.path.realpath(__file__))
TAG, RESCUED = 'veto:polyAfamily', 'rescued'
P = lambda *p: os.path.join(a.outdir, *p)

with open(P('cand.tsv')) as fh:
    rd = csv.reader(fh, delimiter='\t'); head = next(rd); rows = list(rd)
col = {n: i for i, n in enumerate(head)}
K, EV = col['keep'], col['evidence']
for r in rows:   # a rerun: the earlier veto and rescue off again
    ev = r[EV].split(',')
    if TAG in ev: r[EV] = ','.join(e for e in ev if e != TAG) or 'none'; r[K] = '1'
    elif RESCUED in ev: r[EV] = ','.join('overlap' if e == RESCUED else e for e in ev); r[K] = '0'
get = lambda r, n: r[col[n]]
span = lambda r: (int(get(r, 'lltr_start')), int(get(r, 'rltr_end')))
regions = lambda r: ((int(get(r, 'lltr_start')), int(get(r, 'lltr_end')), 'LTR'), (int(get(r, 'lltr_end')) + 1, int(get(r, 'rltr_start')) - 1, 'INT'),
                     (int(get(r, 'rltr_start')), int(get(r, 'rltr_end')), 'LTR'))
records = lambda r: tuple(f"{get(r, 'seq')}:{s}..{e}_{k}" for s, e, k in regions(r))
in_lib = [i for i, r in enumerate(rows) if r[K] == '1']   # the elements lib.raw.fa holds, in its order

def fetch(want):
    """{(seq, start, end): sequence} for 1-based inclusive regions of the genome, upper case."""
    out = {}; by = collections.defaultdict(list)
    for c, s, e in want: by[c].append((s, e))
    if os.path.exists(a.genome + '.fai'):
        fai = {l.split('\t')[0]: [int(x) for x in l.split('\t')[1:5]] for l in open(a.genome + '.fai')}
        with open(a.genome, 'rb') as fa:
            for c, iv in by.items():
                ln, off, lb, lw = fai[c]
                for s, e in iv:
                    b0 = off + (s - 1) // lb * lw + (s - 1) % lb; b1 = off + (e - 1) // lb * lw + (e - 1) % lb
                    fa.seek(b0); out[c, s, e] = fa.read(b1 - b0 + 1).replace(b'\n', b'').decode().upper()
    else:   # no index: one pass over the genome, keeping only the contigs wanted
        def flush(c, parts):
            if c in by:
                x = ''.join(parts)
                for s, e in by[c]: out[c, s, e] = x[s - 1:e].upper()
        c = None; parts = []
        with open(a.genome) as fa:
            for l in fa:
                if l.startswith('>'):
                    flush(c, parts); c = l[1:].split()[0]; parts = []
                elif c in by: parts.append(l.rstrip('\n'))
        flush(c, parts)
    return out

def write(vetoed, rescued=()):
    """cand.tsv (written beside it and renamed) and lib.fa: lib.raw.fa's records of the kept rows other than
    `vetoed`, with the records of the `rescued` rows (from the genome) at their place in the row order."""
    for i in vetoed: rows[i][K] = '0'; rows[i][EV] = (rows[i][EV] + ',' if rows[i][EV] not in ('', 'none') else '') + TAG
    for i in rescued: rows[i][K] = '1'; rows[i][EV] = ','.join(RESCUED if e == 'overlap' else e for e in rows[i][EV].split(','))
    with open(P('cand.tsv.tmp'), 'w') as fo:
        fo.write('\t'.join(head) + '\n')
        for r in rows: fo.write('\t'.join(r) + '\n')
    os.replace(P('cand.tsv.tmp'), P('cand.tsv'))
    keep = {n for r in rows if r[K] == '1' for n in records(r)}
    before = collections.defaultdict(list); last = []   # first record of the next lib.raw.fa element -> rescued rows to write before it
    if rescued:
        seqs = fetch([(get(rows[i], 'seq'), s, e) for i in rescued for s, e, _ in regions(rows[i])])
        nxt = iter(in_lib); j = next(nxt, None)
        for i in sorted(rescued):
            while j is not None and j < i: j = next(nxt, None)
            (before[records(rows[j])[0]] if j is not None else last).append(i)
    def put(fo, i):
        r = rows[i]
        for (s, e, k), name in zip(regions(r), records(r)): fo.write(f">{name}#LTR/{get(r, 'dom_superfamily')}\n{seqs[get(r, 'seq'), s, e]}\n")
    n = kept = 0; on = True
    with open(P('lib.raw.fa')) as fi, open(P('lib.fa.tmp'), 'w') as fo:
        for l in fi:
            if l.startswith('>'):
                name = l.split()[0]; base = name[1:].split('#')[0]
                for i in before.get(base, ()): put(fo, i); n += 3; kept += 3
                on = base in keep; n += 1; kept += on
                if on: fo.write(name + '\n')
            elif on: fo.write(l)
        for i in last: put(fo, i); n += 3; kept += 3
    os.replace(P('lib.fa.tmp'), P('lib.fa'))
    return n, kept

nkept = len(in_lib)
polya_all = [i for i, r in enumerate(rows) if 'veto:polyA' in r[EV].split(',')]
polya = [i for i in polya_all if get(rows[i], 'dom_order') != 'LTR']
if a.off or a.min_depth <= 0 or len(polya) < max(a.min_pool, a.min_depth):
    write([])
    for f in ('polya.fa', 'vetoed.tsv', 'rescued.tsv'):   # an earlier run's files: step 22b would cut with the sample
        if os.path.exists(P('polya_family', f)): os.remove(P('polya_family', f))
    why = 'off' if a.off or a.min_depth <= 0 else (f'not run: {len(polya)} poly-A-vetoed candidates without an LTR domain hit '
                                                  f'({len(polya_all)} in all), fewer than {max(a.min_pool, a.min_depth)}')
    print(f'poly-A family veto: {why}; {nkept} kept elements')
    sys.exit()

os.makedirs(P('polya_family'), exist_ok=True)
sample = sorted(random.Random(a.seed).sample(polya, min(a.sample, len(polya))))
want = [(get(rows[i], 'seq'), int(get(rows[i], 'lltr_start')), int(get(rows[i], 'lltr_end'))) for i in sample]
got = fetch(want)
# two sampled candidates can share their left LTR (gorilla: 1 in 2,000); each copy keeps its own record, as before,
# under a numbered name, since rmblastn >= 3.1 refuses a database with a repeated sequence name
seen = collections.Counter()
with open(P('polya_family', 'polya.fa'), 'w') as fo:
    for c, s, e in want:
        seen[c, s, e] += 1; k = seen[c, s, e]
        fo.write(f'>{c}:{s}..{e}{"" if k == 1 else f"_copy{k}"}\n{got[c, s, e]}\n')

def coverage(fasta, hits):
    """{LTR record: share of it covered by >= --min-depth of the sampled LTRs}, from a search of `fasta`."""
    env = dict(os.environ, WORD='12'); env.pop('MATRIX', None); env.pop('MIN_PID', None)
    subprocess.run(['bash', os.path.join(HERE, 'rmb_ava.sh'), fasta, 'polya.fa', '1', hits,
                    *'qseqid sseqid pident length qstart qend sstart send qlen slen bitscore'.split()],
                   cwd=P('polya_family'), env=env, check=True, stdout=sys.stderr)
    depth = {}
    for l in open(P('polya_family', hits)):
        q, s, pid, ln, qs, qe, ss, se, ql, sl, b = l.split('\t')
        d = depth.get(q)
        if d is None: d = depth[q] = np.zeros(int(ql) + 2, dtype=np.int32)
        lo, hi = sorted((int(qs), int(qe))); d[lo] += 1; d[hi + 1] -= 1
    return {q: float((np.cumsum(d)[1:-1] >= a.min_depth).mean()) for q, d in depth.items()}
def covs(cov, r):
    l5, _, l3 = records(r); return cov.get(l5, 0.), cov.get(l3, 0.)

ltr_names = {n for i in in_lib for n in records(rows[i]) if n.endswith('_LTR')}
with open(P('lib.raw.fa')) as fi, open(P('polya_family', 'ltrs.fa'), 'w') as fo:
    on = False
    for l in fi:
        if l.startswith('>'): n = l[1:].split()[0].split('#')[0]; on = n in ltr_names; l = f'>{n}\n'
        if on: fo.write(l)
cov = coverage('ltrs.fa', 'hits.tsv')
vetoed = []
with open(P('polya_family', 'vetoed.tsv'), 'w') as fo:
    fo.write('elem\tevidence\telem_len\tlltr_len\tidentity\tlltr_cov\trltr_cov\n')
    for i in in_lib:
        r = rows[i]; c5, c3 = covs(cov, r)
        if min(c5, c3) >= a.min_cov:
            vetoed.append(i)
            fo.write(f"{get(r, 'seq')}:{get(r, 'elem_start')}..{get(r, 'elem_end')}\t{r[EV]}\t{get(r, 'elem_len')}\t{get(r, 'lltr_len')}\t{get(r, 'identity')}\t{c5:.3f}\t{c3:.3f}\n")

# ---- overlap rescue ----
rescued = []; n_cand = n_family = n_part = 0
if os.path.exists(P('polya_family', 'rescued.tsv')): os.remove(P('polya_family', 'rescued.tsv'))
if vetoed and not a.no_rescue:
    vset = set(vetoed)
    def spans_of(idx):
        d = collections.defaultdict(list)
        for i in idx: d[get(rows[i], 'seq')].append(span(rows[i]))
        for v in d.values(): v.sort()
        return d
    taken = spans_of(i for i in in_lib if i not in vset); vspans = spans_of(vetoed)
    longest = max(e - s for d in (taken, vspans) for v in d.values() for s, e in v)
    def hits_any(r, spans):
        v = spans.get(get(r, 'seq'), ()); s, e = span(r)
        j = bisect.bisect_left(v, (s - longest,))
        while j < len(v) and v[j][0] <= e:
            if v[j][1] >= s: return True
            j += 1
        return False
    # dropped as overlapping, overlaps a vetoed element and no element still kept: every element it lost to is vetoed
    cands = [i for i, r in enumerate(rows) if 'overlap' in r[EV].split(',') and hits_any(r, vspans) and not hits_any(r, taken)]
    n_cand = len(cands)
    if cands:
        got = fetch([(get(rows[i], 'seq'), s, e) for i in cands for s, e, k in regions(rows[i]) if k == 'LTR'])
        with open(P('polya_family', 'rescue_ltrs.fa'), 'w') as fo:
            for i in cands:
                for (s, e, k), name in zip(regions(rows[i]), records(rows[i])):
                    if k == 'LTR': fo.write(f">{name}\n{got[get(rows[i], 'seq'), s, e]}\n")
        rcov = coverage('rescue_ltrs.fa', 'rescue_hits.tsv')
        def rank(i):   # ltr-forage's --overlaps evidence order
            ev = rows[i][EV]
            n = ('domain' in ev) + ('tsd' in ev) + ('tgca' in ev) + (get(rows[i], 'pbs_trna') not in ('.', ''))
            return (-n, -float(get(rows[i], 'identity')), i)
        new = collections.defaultdict(list)
        with open(P('polya_family', 'rescued.tsv'), 'w') as fo:
            fo.write('elem\tevidence\telem_len\tlltr_len\tidentity\tlltr_cov\trltr_cov\toutcome\n')
            for i in sorted(cands, key=rank):
                r = rows[i]; c5, c3 = covs(rcov, r); c = get(r, 'seq'); s, e = span(r)
                if min(c5, c3) >= a.min_cov: out = 'poly-A family'; n_family += 1
                elif max(c5, c3) > 0: out = 'poly-A repeat in an LTR'; n_part += 1
                elif any(s <= e2 and s2 <= e for s2, e2 in new[c]): out = 'overlaps a rescued element with more evidence'
                else: out = 'rescued'; new[c].append((s, e)); rescued.append(i)
                fo.write(f"{c}:{get(r, 'elem_start')}..{get(r, 'elem_end')}\t{r[EV]}\t{get(r, 'elem_len')}\t{get(r, 'lltr_len')}\t{get(r, 'identity')}\t{c5:.3f}\t{c3:.3f}\t{out}\n")
n, kept = write(vetoed, rescued)
print(f'poly-A family veto: {len(vetoed)} of {nkept} kept elements have LTRs that are copies of poly-A-terminated repeats '
      f'(>= {a.min_cov:.0%} of each LTR covered by >= {a.min_depth} of {len(sample)} sampled LTRs of {len(polya)} poly-A-vetoed candidates; '
      f'{len(polya_all) - len(polya)} more left out for an LTR domain hit); '
      + ('' if a.no_rescue else f'{len(rescued)} element(s) that had lost the overlap choice only to vetoed ones kept again '
         f'({n_cand} such candidates, {n_family} of them poly-A family, {n_part} with such a repeat in an LTR); ')
      + f'{kept} of {n} records left')
