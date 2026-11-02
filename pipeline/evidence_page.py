#!/usr/bin/env python3
"""Per-family evidence page for an ltr_library.py run: one dossier per family, built from its founding element.

usage: evidence_page.py <run_dir> <genome.fa> -o evidence.html [--families SEL] [--reference LIB]
                        [--mask OUT] [--standard OUT] [--other NAME=OUT] [--title T] [--threads N]

For each family (library.fa / families.tsv) the founding element (the origin of its first internal model,
else of its first LTR model) is drawn as in a genome browser: LTRs, assembly gaps, TSD, primer binding site
(ltr-forage's tRNA search), polypurine tract and protein domains (nail search of the six-frame ORFs against
--domain-hmm, as hmmscan would).  The element cluster behind the family gives its full-length copies and how many sequences
they lie on.  Optional evidence:
  --reference LIB   a curated library: nearest family by rmblastn (LTR, internal) and by nail on the six-frame
                    translations of the internal regions (the tblastx comparison)
  --mask OUT        RepeatMasker .out of the genome masked with this run's library: hits and copies
  --standard OUT    RepeatMasker .out of the curated annotation: bases inside / outside its LTR
                    annotation (TP / FP) and what it calls the rest
  --other NAME=OUT  another library's mask (e.g. HiRepEat): share of this family's FP bases it masks too

--families: 'all' (default), 'top-fp:N' (the N families with the most FP masking; needs --mask and
--standard), or a comma list of family names (fam0042,fam0007).  A verdict is computed per family from
the evidence; it is a triage aid, not a call.  External programs: the Rust rmblastn (through
pipeline/rmb_ava.sh, and for the 5' vs 3' LTR alignment), nail and MMseqs2.
"""
import argparse, bisect, collections, csv, json, os, re, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
RMB_AVA = os.path.join(HERE, 'rmb_ava.sh')
PICK_MATRIX = os.path.join(HERE, 'pick_matrix.sh')
FLANK = 30
LTR_ORDER = ('Gypsy', 'Copia', 'ERV', 'LTR', 'BEL', 'Pao', 'DIRS')


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout


def read_fa(p):
    d, n = {}, None
    for l in open(p):
        l = l.rstrip('\n')
        if l.startswith('>'):
            n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l.strip())
    return {k: ''.join(v) for k, v in d.items()}


class Genome:
    """Random access by a .fai index (written if missing)."""
    def __init__(self, path):
        self.f = open(path, 'rb'); fai = path + '.fai'
        if not os.path.exists(fai): fai = self.index(path)
        self.fai = {l.split()[0]: list(map(int, l.split()[1:5])) for l in open(fai)}

    @staticmethod
    def index(path):
        out = os.path.join(tempfile.mkdtemp(), 'g.fai'); rows = []; name = None
        with open(path, 'rb') as f:
            off = 0
            for line in f:
                if line.startswith(b'>'):
                    name = line[1:].split()[0].decode(); rows.append([name, 0, off + len(line), 0, 0])
                elif name:
                    r = rows[-1]
                    if r[3] == 0: r[3], r[4] = len(line.rstrip(b'\r\n')), len(line)
                    r[1] += len(line.rstrip(b'\r\n'))
                off += len(line)
        with open(out, 'w') as fh:
            for r in rows: fh.write('\t'.join(map(str, r)) + '\n')
        return out

    def get(self, c, s, e):   # 1-based inclusive, clipped
        L, off, lb, lB = self.fai[c]; s = max(1, s); e = min(L, e)
        if e < s: return ''
        s -= 1
        a = off + s // lb * lB + s % lb; b = off + (e - 1) // lb * lB + (e - 1) % lb
        self.f.seek(a); return self.f.read(b - a + 1).replace(b'\n', b'').replace(b'\r', b'').decode().upper()


CODON = {}
_b = 'TCAG'; _aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for i, (x, y, z) in enumerate((x, y, z) for x in _b for y in _b for z in _b): CODON[x + y + z] = _aa[i]
def rc(s): return s[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))
def orfs(seq):
    """Six-frame stop-to-stop ORFs of >= 60 aa: (id, aa, strand, frame, aa offset)."""
    for strand, s in (('+', seq), ('-', rc(seq))):
        for fr in range(3):
            aa = ''.join(CODON.get(s[j:j + 3], 'X') for j in range(fr, len(s) - 2, 3)); pos = 0
            for piece in aa.split('*'):
                if len(piece) >= 60: yield f'{strand}{fr}_{pos}', piece, strand, fr, pos
                pos += len(piece) + 1


def frames(seq):
    """The six frames of seq as proteins, stops as X (nail's alphabet has no '*'): (frame 0-5, aa)."""
    for k, x in enumerate((seq, rc(seq))):
        for fr in range(3):
            yield 3 * k + fr, ''.join(CODON.get(x[j:j + 3], 'X') for j in range(fr, len(x) - 2, 3)).replace('*', 'X')


def n_runs(s):
    return [[m.start() + 1, m.end()] for m in re.finditer(r'N+', s)]


def merge(v):
    o = []
    for a, b in sorted(v):
        if o and a <= o[-1][1] + 1: o[-1][1] = max(o[-1][1], b)
        else: o.append([a, b])
    return o


class Intervals:
    def __init__(self, d): self.iv = {s: merge(v) for s, v in d.items()}; self.st = {s: [x[0] for x in v] for s, v in self.iv.items()}
    def parts(self, s, a, b):
        """Pieces of [a, b] covered by the intervals."""
        iv = self.iv.get(s)
        if not iv: return []
        i = max(0, bisect.bisect_right(self.st[s], a) - 1); out = []
        for x, y in iv[i:]:
            if x > b: break
            lo, hi = max(a, x), min(b, y)
            if lo <= hi: out.append((lo, hi))
        return out


def out_rows(p):
    for l in open(p):
        c = l.split()
        if len(c) >= 15 and c[0].isdigit(): yield c


def best_family(rows, qlen):
    """Per target family: summed bitscore, weighted identity, union query coverage, HSPs; best three.
    rows: (sseqid, pident, length, qstart, qend, bitscore)."""
    fam = collections.defaultdict(list)
    for s, pid, ln, qs, qe, bits in rows:
        fam[re.sub(r'_(LTR|INT)$', '', s.split('#')[0])].append((int(qs), int(qe), float(pid), int(ln), float(bits)))
    res = []
    for s, h in fam.items():
        u = sum(b - a + 1 for a, b in merge([(min(a, b), max(a, b)) for a, b, *_ in h]))
        wid = sum(p * n for _, _, p, n, _ in h) / sum(n for _, _, p, n, _ in h)
        res.append({'target': s, 'bits': round(sum(x[4] for x in h), 1), 'ident': round(wid, 1), 'qcov': round(min(1, u / max(1, qlen)), 3),
                    'hsps': [{'qs': min(a, b), 'qe': max(a, b), 'pid': p} for a, b, p, _, _ in h]})
    res.sort(key=lambda r: -r['bits'])
    return res[:3]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('run'); ap.add_argument('genome'); ap.add_argument('-o', '--out', required=True)
    ap.add_argument('--families', default='all'); ap.add_argument('--reference'); ap.add_argument('--mask')
    ap.add_argument('--standard'); ap.add_argument('--other', help='NAME=mask.out')
    ap.add_argument('--title'); ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--domain-hmm', default=os.environ.get('HMM', '/usr/local/LTR_retriever-2.9.0/database/TEfam.hmm'))
    ap.add_argument('--work', help='scratch directory (default: <run>/evidence)')
    ap.add_argument('--rmblastn', default=os.environ.get('RMBLASTN', 'rmblastn'), help='the Rust rmblastn')
    ap.add_argument('--nail', default=os.environ.get('NAIL', 'nail'))
    ap.add_argument('--mmseqs', default=os.environ.get('MMSEQS', 'mmseqs'))
    a = ap.parse_args()
    R = os.path.abspath(a.run); work = os.path.abspath(a.work or os.path.join(R, 'evidence')); os.makedirs(work, exist_ok=True)
    G = Genome(a.genome)

    lib = read_fa(os.path.join(R, 'library.fa'))
    libname = {n.split('#')[0]: n for n in lib}
    models = collections.defaultdict(list)          # family -> rows of families.tsv
    for r in csv.DictReader(open(os.path.join(R, 'families.tsv')), delimiter='\t'):
        if r['new_name'] in libname: models[r['family']].append(r)
    cand = {f"{x['seq']}:{x['elem_start']}..{x['elem_end']}": x for x in csv.DictReader(open(os.path.join(R, 'cand.tsv')), delimiter='\t')}
    clus = json.load(open(os.path.join(R, 'elem', 'clusters_nested.json')))
    members = collections.defaultdict(list)
    for e, c in clus.items(): members[c].append(e)
    sd = set()
    p = os.path.join(R, 'sd_clusters.tsv')
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter='\t'):
            if r.get('call') == 'sd': sd.add(r['cluster'])
    by_seq = collections.defaultdict(list)
    for e in cand:
        m = re.match(r'(\S+):(\d+)\.\.(\d+)$', e); by_seq[m[1]].append((int(m[2]), int(m[3]), e))
    def element_of(rec):    # record name (elem_LTR / elem_INT, possibly a sub-range) -> candidate element
        # library records come from kept candidates; with --overlaps evidence cand.tsv also holds the dropped
        # calls overlapping them, so prefer a kept one, then the shortest that contains the record
        m = re.match(r'(\S+):(\d+)\.\.(\d+)_(LTR|INT)$', rec)
        if not m: return None
        s, x, y = m[1], int(m[2]), int(m[3])
        hit = [(cand[e].get('keep') != '1', q_ - p_, e) for p_, q_, e in by_seq.get(s, ()) if p_ <= x and y <= q_]
        return min(hit)[2] if hit else None

    # ---- genome masking, per model ------------------------------------------------------------
    mstat = collections.defaultdict(lambda: {'hits': 0, 'bp': 0, 'tp': 0, 'fp': 0, 'full': 0, 'seqs': set(), 'fpcls': collections.Counter(),
                                             'fp_other': 0})
    totals = {}
    if a.mask:
        std = oth = ann = None
        if a.standard:
            ltr = collections.defaultdict(list); other = collections.defaultdict(list)
            for c in out_rows(a.standard):
                (ltr if c[10].startswith('LTR') else other)[c[4]].append((int(c[5]), int(c[6]), c[10].split('/')[0]))
            std = Intervals({s: [(x, y) for x, y, _ in v] for s, v in ltr.items()})
            ann = {}
            for s, v in other.items():
                byc = collections.defaultdict(list)
                for x, y, k in v: byc[k].append((x, y))
                ann[s] = {k: Intervals({s: w}) for k, w in byc.items()}
            totals['standard_ltr_bp'] = sum(y - x + 1 for v in std.iv.values() for x, y in v)
        if a.other:
            oname, opath = a.other.split('=', 1)
            d = collections.defaultdict(list)
            for c in out_rows(opath): d[c[4]].append((int(c[5]), int(c[6])))
            oth = Intervals(d); totals['other'] = oname
        mlen = {n.split('#')[0]: len(s) for n, s in lib.items()}
        allmask = collections.defaultdict(list)
        for c in out_rows(a.mask):
            f, s, x, y = c[9], c[4], int(c[5]), int(c[6]); st = mstat[f]; allmask[s].append((x, y))
            st['hits'] += 1; st['bp'] += y - x + 1; st['seqs'].add(s)
            rs, re_ = (int(c[11]), int(c[12])) if c[8] == '+' else (int(c[13]), int(c[12]))
            if abs(re_ - rs) + 1 >= 0.5 * mlen.get(f, 1e18): st['full'] += 1
            if std is None: continue
            inside = std.parts(s, x, y); tp = sum(q - p_ + 1 for p_, q in inside); st['tp'] += tp; st['fp'] += y - x + 1 - tp
            if y - x + 1 - tp <= 0: continue
            # FP pieces: [x, y] minus the LTR annotation
            pieces = []; cur = x
            for p_, q in inside:
                if p_ > cur: pieces.append((cur, p_ - 1))
                cur = q + 1
            if cur <= y: pieces.append((cur, y))
            for p_, q in pieces:
                rest = q - p_ + 1
                for k, iv in (ann or {}).get(s, {}).items():
                    o = min(rest, sum(v - u + 1 for u, v in iv.parts(s, p_, q)))
                    if o: st['fpcls'][k] += o; rest -= o
                st['fpcls']['unannotated'] += rest
                if oth: st['fp_other'] += sum(v - u + 1 for u, v in oth.parts(s, p_, q))
        if std is not None:
            # genome totals over merged masking (per-family numbers can overlap one another; these do not)
            tp = fp = fpo = fpu = 0
            for sq, v in allmask.items():
                for x, y in merge(v):
                    inside = std.parts(sq, x, y); t_ = sum(q - p_ + 1 for p_, q in inside); tp += t_; fp += y - x + 1 - t_
                    cur = x; pieces = []
                    for p_, q in inside:
                        if p_ > cur: pieces.append((cur, p_ - 1))
                        cur = q + 1
                    if cur <= y: pieces.append((cur, y))
                    for p_, q in pieces:
                        o = sum(min(q, v2) - max(p_, u) + 1 for iv in (ann or {}).get(sq, {}).values() for u, v2 in merge(iv.parts(sq, p_, q)))
                        fpu += max(0, q - p_ + 1 - o)
                        if oth: fpo += sum(v2 - u + 1 for u, v2 in oth.parts(sq, p_, q))
            totals.update(tp=tp, fp=fp, fp_unannotated=fpu)
            if oth: totals['fp_other'] = fpo

    def fam_mask(fam):
        out = {'hits': 0, 'bp': 0, 'tp': 0, 'fp': 0, 'fp_other': 0, 'fpcls': collections.Counter(), 'seqs': set(), 'ltr_full': 0, 'int_full': 0}
        for r in models[fam]:
            st = mstat.get(r['new_name'])
            if not st: continue
            for k in ('hits', 'bp', 'tp', 'fp', 'fp_other'): out[k] += st[k]
            out['fpcls'].update(st['fpcls']); out['seqs'] |= st['seqs']
            out['ltr_full' if r['kind'] == 'LTR' else 'int_full'] += st['full']
        out['seqs'] = len(out['seqs']); out['fpcls'] = out['fpcls'].most_common(4)
        return out

    # ---- which families ---------------------------------------------------------------------------
    fams = sorted(models)
    if a.families.startswith('top-fp:'):
        if not (a.mask and a.standard): sys.exit('evidence_page: top-fp needs --mask and --standard')
        n = int(a.families.split(':')[1]); fams = sorted(fams, key=lambda f: -fam_mask(f)['fp'])[:n]
    elif a.families != 'all':
        want = a.families.split(','); fams = [f for f in want if f in models]
        missing = set(want) - set(fams)
        if missing: print('evidence_page: not in this library:', ' '.join(sorted(missing)), file=sys.stderr)

    # ---- profile names and reference libraries ---------------------------------------------------
    hmm_desc, name = {}, None      # accession (no version) -> name, description
    for l in open(a.domain_hmm):
        if l.startswith('NAME'): name = l.split()[1]
        elif l.startswith('ACC'): hmm_desc[l.split()[1].split('.')[0]] = {'name': name, 'desc': ''}; acc = l.split()[1].split('.')[0]
        elif l.startswith('DESC'): hmm_desc[acc]['desc'] = l[6:].strip()
    ref = read_fa(a.reference) if a.reference else {}
    ref_int = {n: s_ for n, s_ in ref.items() if not re.search(r'(_LTR|-LTR|LTR$)', n.split('#')[0])}
    if a.reference:
        with open(os.path.join(work, 'ref.fa'), 'w') as fh:
            for n, s_ in ref.items(): fh.write(f'>{n}\n{s_}\n')
    with open(os.path.join(work, 'own_ltr.fa'), 'w') as fh:
        for n, s_ in lib.items():
            if '_LTR' in n: fh.write(f'>{n.split("#")[0]}\n{s_}\n')
    env = dict(os.environ, RMBLASTN=a.rmblastn, THREADS=str(a.threads))

    def rmb(queries, db, tag):
        """rmblastn (RepeatMasker -lib settings, E <= 1e-10) of {name: seq} against a FASTA file:
        {query: [(sseqid, pident, length, qstart, qend, bitscore)]}."""
        out = collections.defaultdict(list)
        if not queries: return out
        q = os.path.join(work, f'{tag}.fa'); o = os.path.join(work, f'{tag}.tsv')
        with open(q, 'w') as fh:
            for n, s_ in queries.items(): fh.write(f'>{n}\n{s_}\n')
        subprocess.run(['bash', RMB_AVA, q, db, '50', o, 'qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 'bitscore'],
                       check=True, capture_output=True, env=env)
        for l in open(o):
            qq, ss, pid, ln, qs, qe, bits = l.rstrip('\n').split('\t')
            out[qq].append((ss, float(pid), int(ln), int(qs), int(qe), float(bits)))
        return out

    def domains(cases):
        """Protein domains of each founding element: nail of the profiles against all the elements' ORFs in
        one run, E-values as hmmscan gives them (Z = number of profiles), E <= 1e-3, the best of overlapping
        hits kept.  Sets each case's 'domains'."""
        d = os.path.join(work, 'nail_dom'); os.makedirs(d, exist_ok=True); meta = {}
        with open(f'{d}/orfs.fa', 'w') as o:
            for c in cases:
                for oid, aa, strand, fr, pos in orfs(c['_elem']):
                    o.write(f'>{c["family"]}|{oid}\n{aa}\n'); meta[f'{c["family"]}|{oid}'] = (strand, fr, pos)
        # nail 0.7.1 panics on a blank line in the profile file (TEfam.hmm starts with one)
        with open(a.domain_hmm) as fi, open(f'{d}/profiles.hmm', 'w') as fo:
            for l in fi:
                if l.strip(): fo.write(l)
        found = collections.defaultdict(list)
        if meta:
            subprocess.run([a.nail, 'search', '-t', str(a.threads), '-E', '1e-3', '-Z', str(max(1, len(hmm_desc))), '-S', '0.01', '-X',
                            '--use-accession', '--mmseqs-path', a.mmseqs, '--tmp-dir', f'{d}/tmp', '--tbl-out', f'{d}/hits.tbl',
                            f'{d}/profiles.hmm', f'{d}/orfs.fa'], check=True, capture_output=True)
            # target query tstart tend qstart qend score bias evalue cellfrac
            for l in open(f'{d}/hits.tbl'):
                if l.startswith('#') or not l.strip(): continue
                f = l.split(); fam, _ = f[0].split('|', 1); acc = f[1].split('.')[0]; ev = float(f[8])
                strand, fr, pos = meta[f[0]]; L = len(by_fam[fam]['_elem'])
                n1, n2 = fr + 3 * (pos + int(f[2]) - 1), fr + 3 * (pos + int(f[3]))
                if strand == '-': n1, n2 = L - n2, L - n1
                hd = hmm_desc.get(acc, {})
                found[fam].append({'name': hd.get('name', acc), 'acc': acc, 'desc': hd.get('desc', ''), 'ievalue': ev,
                                   'start': n1 + 1, 'end': n2, 'strand': strand})
        for c in cases:
            doms = sorted(found[c['family']], key=lambda r: r['ievalue']); kept = []
            for dm in doms:
                if all(dm['end'] < k2['start'] or dm['start'] > k2['end'] for k2 in kept): kept.append(dm)
            c['domains'] = sorted(kept, key=lambda r: r['start'])

    def near_protein(queries):
        """The tblastx comparison: nail of the six frames of each {name: internal region} against the six
        frames of the reference internal sequences, query positions back in nucleotides."""
        out = collections.defaultdict(list)
        if not queries or not ref_int: return out
        d = os.path.join(work, 'nail'); os.makedirs(d, exist_ok=True)
        names = list(ref_int); n_t = 0
        with open(f'{d}/query.fa', 'w') as fq:
            for n, s_ in queries.items():
                for k, aa in frames(s_):
                    if len(aa) >= 10: fq.write(f'>{n}|{k}\n{aa}\n')
        with open(f'{d}/target.fa', 'w') as ft:
            for i_, n in enumerate(names):
                for k, aa in frames(ref_int[n]):
                    if len(aa) >= 10: ft.write(f'>{i_}|{k}\n{aa}\n'); n_t += 1
        subprocess.run([a.nail, 'search', '-t', str(a.threads), '-E', '1e-5', '-Z', str(max(1, n_t)), '-X', '--tbl-format', 'blast',
                        '--mmseqs-path', a.mmseqs, '--tmp-dir', f'{d}/tmp', '--tbl-out', f'{d}/hits.tbl', f'{d}/query.fa', f'{d}/target.fa'],
                       check=True, capture_output=True)
        # query target pid length mismatch gapopen qstart qend tstart tend evalue score
        for l in open(f'{d}/hits.tbl'):
            if l.startswith('#') or not l.strip(): continue
            f = l.split(); qn, k = f[0].rsplit('|', 1); k = int(k); L = len(queries[qn])
            x, y = k % 3 + 3 * (int(f[6]) - 1) + 1, k % 3 + 3 * int(f[7])
            if k >= 3: x, y = L - y + 1, L - x + 1
            out[qn].append((names[int(f[1].rsplit('|', 1)[0])], float(f[2]), int(f[3]), x, y, float(f[11])))
        return out

    # ---- per family: the founding element -----------------------------------------------------------
    def prep(fam):
        rows = models[fam]
        pick = next((r for r in rows if r['kind'] == 'INT'), rows[0])
        e = element_of(pick['origin']) or element_of(rows[0]['origin'])
        if e is None or e not in cand: return None
        x = cand[e]; c = x['seq']
        es, ee = int(x['elem_start']), int(x['elem_end'])
        ls_, le = int(x['lltr_start']), int(x['lltr_end']); rs, re_ = int(x['rltr_start']), int(x['rltr_end'])
        elem = G.get(c, es, ee); lltr = G.get(c, ls_, le); rltr = G.get(c, rs, re_); intr = G.get(c, le + 1, rs - 1)
        d = os.path.join(work, fam); os.makedirs(d, exist_ok=True)
        for k, s_ in (('lltr', lltr), ('rltr', rltr)): open(f'{d}/{k}.fa', 'w').write(f'>{k}\n{s_}\n')
        # 5' vs 3' LTR: the best HSP, as rmblastn's pairwise text
        mx = subprocess.run(['bash', PICK_MATRIX, f'{d}/rltr.fa'], check=True, capture_output=True, text=True).stdout.strip()
        pair = subprocess.run([a.rmblastn, '--query', f'{d}/lltr.fa', '--db', f'{d}/rltr.fa', '--matrix', mx, '--gapopen', '24',
                               '--gapextend', '6', '--dust', 'no', '--num-threads', '1'], capture_output=True, text=True).stdout
        m = re.search(r'( Score = .*?)(?=\n\n Score = |\n\n\n|\Z)', pair, re.S)
        cl = clus.get(e); mem = members.get(cl, [e])
        mseqs = collections.Counter(re.match(r'(\S+):', y)[1] for y in mem)
        pbs = None
        if x.get('pbs_trna', '.') not in ('.', ''):
            pbs = {'trna': x['pbs_trna'], 'mm12': int(x['pbs_mm12']), 'edist18': x['pbs_edist18'], 'strand': x['pbs_strand'],
                   'offset': int(x['pbs_offset'])}
        return {
            'family': fam, 'id': e, 'seq': c, 'start': es, 'end': ee, 'len': ee - es + 1,
            'models': [{'name': r['new_name'], 'kind': r['kind'], 'len': int(r['len']), 'class': r['class']} for r in rows],
            'class': next((r['class'] for r in rows if r['class']), 'LTR'),
            'lltr': [ls_, le], 'rltr': [rs, re_], 'ltr_identity': float(x['identity']),
            'tsd_len': int(x['tsd_len']), 'tsd_seq': x['tsd_seq'], 'tsd_mm': x.get('tsd_mm', '0') not in ('0', '', '.'),
            'motif': x['motif'], 'motif_class': x['motif_class'],
            'ppt_offset': int(x['ppt_offset']), 'ppt_purines15': int(x['ppt_purines15']),
            'dom_nest': x.get('dom_nest', '.') not in ('.', '', '0'), 'evidence': x['evidence'],
            'n_runs': n_runs(elem), 'pbs': pbs,
            'flank5': G.get(c, es - FLANK, es - 1), 'ltr5_head': lltr[:30], 'ltr3_tail': rltr[-30:], 'flank3': G.get(c, ee + 1, ee + FLANK),
            'after5ltr': G.get(c, le + 1, le + 40), 'before3ltr': G.get(c, rs - 40, rs - 1), 'after3ltr_rc': rc(G.get(c, rs - 40, rs - 1)),
            'ltr_pair_aln': m.group(1).strip() if m else '',
            'cluster': {'elements': len(mem), 'seqs': len(mseqs), 'top_seq_share': round(max(mseqs.values()) / len(mem), 2),
                        'sd': cl in sd},
            'mask': fam_mask(fam) if a.mask else None,
            '_lltr': lltr, '_int': intr, '_elem': elem,
        }

    if a.mask and a.standard: fams = sorted(fams, key=lambda f: -fam_mask(f)['fp'])   # most FP first
    with ThreadPoolExecutor(a.threads) as ex: cases = [c for c in ex.map(prep, fams) if c]
    by_fam = {c['family']: c for c in cases}
    domains(cases)

    # ---- nearest families: one batch per search ---------------------------------------------------------
    lltrs = {c['family']: c['_lltr'] for c in cases}
    ints = {c['family']: c['_int'] for c in cases if c['_int'].strip('N')}
    own = rmb(lltrs, os.path.join(work, 'own_ltr.fa'), 'q_own')
    if a.reference:
        nl, ni, na = rmb(lltrs, os.path.join(work, 'ref.fa'), 'q_ltr'), rmb(ints, os.path.join(work, 'ref.fa'), 'q_int'), near_protein(ints)
    for c in cases:
        f = c['family']; lltr, intr = c.pop('_lltr'), c.pop('_int'); c.pop('_elem')
        c['near_ltr'] = best_family(nl[f], len(lltr)) if a.reference else []
        c['near_int'] = best_family(ni[f], len(intr)) if a.reference else []
        c['near_int_aa'] = best_family(na[f], len(intr)) if a.reference else []
        c['near_own'] = best_family([h for h in own[f] if not h[0].startswith(f + '_')], len(lltr))[:2]
    for c in cases: c['verdict'] = verdict(c, bool(a.reference))

    pbs_use = collections.Counter(); n_elem = 0
    for x in cand.values():
        if x.get('keep') != '1': continue
        n_elem += 1
        if x.get('pbs_trna', '.') not in ('.', ''): pbs_use[x['pbs_trna']] += 1
    meta = {'title': a.title or f'LTR family evidence: {os.path.basename(R)}', 'run': os.path.basename(R),
            'genome': os.path.basename(a.genome), 'families_total': len(models), 'shown': len(cases), 'selection': a.families,
            'reference': os.path.basename(a.reference) if a.reference else None, 'totals': totals,
            'pbs': {'elements': n_elem, 'with_pbs': sum(pbs_use.values()), 'top': pbs_use.most_common(8)}}
    t = open(os.path.join(HERE, 'evidence_template.html')).read()
    open(a.out, 'w').write(t.replace('/*DATA*/', json.dumps(cases)).replace('/*META*/', json.dumps(meta)))
    print(f'evidence_page: {len(cases)} families -> {a.out}; ' +
          ', '.join(f'{k} {v}' for k, v in collections.Counter(c['verdict'][1].split(' (')[0] for c in cases).most_common()))


def verdict(c, have_ref):
    """(level, label, reasons): a triage call from the evidence, most specific first."""
    doms = c['domains']; why_for, why_against = [], []
    te = [d for d in doms if not re.search(r'transposase|MULE|DDE|helitron', d['name'] + d['desc'], re.I)]
    strands = {d['strand'] for d in te}
    if c['dom_nest'] or len(strands) > 1:
        why_against.append('protein domains on both strands of the founding element: likely two elements, one nested in the other')
    if c['cluster']['sd']: why_against.append('its copies share both flanks: a segmental duplication')
    mk = c['mask']
    few = not mk or mk['ltr_full'] + mk['int_full'] < 10     # a local duplication leaves few copies genome-wide
    if c['cluster']['elements'] >= 2 and c['cluster']['seqs'] == 1 and c['cluster']['elements'] <= 3 and few:
        why_against.append(f"all {c['cluster']['elements']} full-length copies on one sequence, few elsewhere")
    if len(te) >= 2: why_for.append(f'{len(te)} retrotransposon protein domains')
    elif len(te) == 1: why_for.append('one retrotransposon protein domain')
    if c['tsd_len'] >= 4: why_for.append(f"{c['tsd_len']} bp TSD")
    if c['motif_class'] == 'TGCA': why_for.append('TG…CA ends')
    if c['pbs']: why_for.append(f"PBS for tRNA-{c['pbs']['trna']}")
    n = c['cluster']['elements']; s = c['cluster']['seqs']
    if n >= 3 and s >= 2: why_for.append(f'{n} full-length copies on {s} sequences')
    if mk and mk['ltr_full'] >= 10: why_for.append(f"{mk['ltr_full']} near-full-length LTR copies in the genome")
    rel = None
    if have_ref:
        for k, arr in (('LTR', c['near_ltr']), ('internal', c['near_int'])):
            h = arr[0] if arr else None
            if h and h['ident'] >= 80 and h['qcov'] >= 0.8: rel = (k, h); break
    reasons = why_for + why_against
    if why_against:
        return ['bad', 'Possible artifact', reasons]
    if rel:
        k, h = rel
        return ['warn', f"Known family ({h['target']})", [f"{k} {h['ident']}% identical to {h['target']} over {h['qcov']:.0%}"] + why_for]
    # lines of evidence: two or more domains, TSD, TG..CA, PBS, copies on several sequences, many genome copies.
    # Three suffice without any domain (non-autonomous elements keep the structure but lose the coding region)
    strong = sum(1 for w in why_for if not w.startswith('one '))
    if strong >= 3 or (len(te) >= 2 and strong >= 2):
        return ['ok', 'Likely novel family' if have_ref else 'Well supported', reasons]
    return ['warn', 'Weak evidence', reasons or ['no structural or coding evidence beyond the LTR pair']]


if __name__ == '__main__':
    main()
