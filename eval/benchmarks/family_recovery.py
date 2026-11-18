#!/usr/bin/env python3
"""Family-level recovery of a curated LTR library by a candidate library.

Usage: family_recovery.py <curated.fa> <library.fa> [--min-identity 80] [--out table.tsv]
                          [--standard genome.fa.out]

Every curated sequence is an LTR or an internal region by its name (see
family_kind; --standard names the internals of a RepeatMasker library) and is aligned
against the candidate library with the Rust rmblastn (pipeline/rmb_ava.sh:
RepeatMasker -lib settings, E <= 1e-10; RMBLASTN, THREADS from the environment).  For each curated
sequence the script takes the union of the HSPs at or above the identity cut
and reports the fraction of it covered, then per family the best LTR and
best INT coverage.  The summary counts families whose LTR, INT, or both reach
50%, 80% and 95% coverage.  Weighting is per family, not per base.  "both"
needs a library that names a family's LTR and internal alike (Os0008_LTR /
Os0008_INT); RepeatMasker's names do not (LTR7 / HERVH).  The last line
gives the sequences left out, with no LTR or INT in the name and no LTR class.
"""
import sys, os, re, subprocess, tempfile, collections, argparse

HERE = os.path.dirname(os.path.abspath(__file__))

ap = argparse.ArgumentParser()
ap.add_argument('curated'); ap.add_argument('library')
ap.add_argument('--min-identity', type=float, default=80.0)
ap.add_argument('--out', default=None)
ap.add_argument('--standard', help='RepeatMasker .out of the curated library on the genome: a curated <name> '
                                   'annotated there as <name>-int is an internal region (HERVH)')
a = ap.parse_args()

def read_fasta(p):
    seqs, name = {}, None
    for l in open(p):
        if l.startswith('>'):
            name = l[1:].split()[0]; seqs[name] = []
        elif name: seqs[name].append(l.strip())
    return {k: ''.join(v) for k, v in seqs.items()}

def family_kind(name):
    """Family and part (LTR or INT) from a curated name.
    Named parts: the Jiang libraries (Os0008_LTR, Os0020_INT_Bajie), Repbase (ATCOPIA44LTR, ATLANTYS2_I,
    Copia-3_AT-LTR), RepeatMasker (MER50-int).  A bare trailing I is an internal only next to an LTR record
    of the same family (ATGP8I with ATGP8LTR; not MLT1I, an LTR, nor ATLINEIII).  Any other name is a
    RepeatMasker-library one and needs an LTR class in the header (LTR7#LTR/ERV1): an internal when the
    --standard annotation shows it as <name>-int and the library has no <name>-int of its own (HERVH), an
    LTR otherwise.  A name with neither a part nor an LTR class is left out ('other': the LINEs and SINEs
    of the Arabidopsis Repbase file)."""
    base, _, cls = name.partition('#')
    parts = base.split('_')
    if 'LTR' in parts[1:]: return parts[0], 'LTR'
    if 'INT' in parts[1:]: return parts[0], 'INT'
    m = re.match(r'^(.+?)[-_]?LTR$', base)
    if m: return m.group(1), 'LTR'
    m = re.match(r'^(.+?)(?:[-_](?:INT|int|I))$', base)
    if m: return m.group(1), 'INT'
    if base.endswith('I') and any(base[:-1] + x in cur_bases for x in ('LTR', '_LTR', '-LTR')): return base[:-1], 'INT'
    if not cls.startswith('LTR'): return base, 'other'
    if base + '-int' in std_int and base + '-int' not in cur_bases: return base, 'INT'
    return base, 'LTR'

cur = read_fasta(a.curated)
cur_bases = {n.split('#')[0] for n in cur}
std_int = set()
if a.standard:
    for l in open(a.standard):
        c = l.split()
        if len(c) >= 11 and c[0].isdigit() and c[10].startswith('LTR') and c[9].endswith('-int'): std_int.add(c[9])
with tempfile.TemporaryDirectory() as td:
    out = os.path.join(td, 'hits.tsv')
    env = dict(os.environ, THREADS=os.environ.get('THREADS', '8'))
    subprocess.run(['bash', os.path.join(HERE, '..', '..', 'pipeline', 'rmb_ava.sh'), a.curated, a.library, '50', out,
                    'qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'bitscore'],
                   check=True, capture_output=True, env=env)
    res = open(out).read()
hsps = collections.defaultdict(list)
for l in res.splitlines():
    q, s, pid, ln, qs, qe, ss, se, bits = l.split('\t')
    if float(pid) < a.min_identity: continue
    hsps[q].append((int(qs), int(qe), float(pid)))

def coverage(name):
    L = len(cur[name]); iv = sorted((min(s, e), max(s, e)) for s, e, _ in hsps.get(name, []))
    cov, cur_s, cur_e = 0, None, None
    for s, e in iv:
        if cur_s is None: cur_s, cur_e = s, e
        elif s <= cur_e + 1: cur_e = max(cur_e, e)
        else: cov += cur_e - cur_s + 1; cur_s, cur_e = s, e
    if cur_s is not None: cov += cur_e - cur_s + 1
    return cov / L if L else 0.0

fam = collections.defaultdict(lambda: {'LTR': 0.0, 'INT': 0.0, 'has': set()})
skipped = []
for name in cur:
    f, k = family_kind(name)
    if k == 'other': skipped.append(name); continue
    fam[f]['has'].add(k)
    fam[f][k] = max(fam[f][k], coverage(name))
rows = sorted(fam.items())
if a.out:
    with open(a.out, 'w') as o:
        o.write('family\tltr_cov\tint_cov\n')
        for f, d in rows: o.write(f"{f}\t{d['LTR']:.3f}\t{d['INT']:.3f}\n")
n = len(rows); n_ltr = sum(1 for _, d in rows if 'LTR' in d['has']); n_int = sum(1 for _, d in rows if 'INT' in d['has'])
print(f"{os.path.basename(a.library)} vs {os.path.basename(a.curated)}: {n} families ({n_ltr} with an LTR, {n_int} with an INT), identity >= {a.min_identity:.0f}%")
for t in (0.5, 0.8, 0.95):
    l = sum(1 for _, d in rows if 'LTR' in d['has'] and d['LTR'] >= t)
    i = sum(1 for _, d in rows if 'INT' in d['has'] and d['INT'] >= t)
    b = sum(1 for _, d in rows if 'LTR' in d['has'] and 'INT' in d['has'] and d['LTR'] >= t and d['INT'] >= t)
    print(f"  coverage >= {t:.2f}: LTR {l}/{n_ltr} ({100*l/max(1,n_ltr):.1f}%)  INT {i}/{n_int} ({100*i/max(1,n_int):.1f}%)  both {b}")
miss = [f for f, d in rows if d['LTR'] < 0.5 and d['INT'] < 0.5]
print(f"  families with neither LTR nor INT half covered: {len(miss)}")
if skipped: print(f"  {len(skipped)} of {len(cur)} curated sequences left out (no LTR / INT in the name, no LTR class), e.g. {', '.join(skipped[:4])}")
