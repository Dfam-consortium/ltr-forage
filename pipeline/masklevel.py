#!/usr/bin/env python3
"""Apply rmblastn's cross-subject --mask-level filter to hits searched with --mask-level 100 (no filter).

usage: masklevel.py filter <hits.tsv> <db.fa> <query.fa> <mask_level> > filtered.tsv
       masklevel.py changed <prev.fa> <cur.fa> <changed.fa>          records of cur.fa new or changed since prev.fa
       masklevel.py carry <prev.fa> <cur.fa> <raw_prev.tsv> > carried.tsv  hits on records left unchanged, renamed
       masklevel.py remap <cuts.tsv> <raw_prev.tsv> [join_tol [min_score [clipj [matrix [query.fa cut_db.fa]]]]] > remapped.tsv
                                                 hits moved onto the records excise_by_depth.py --cuts changed, no search

Records match by name up to '#' (excise_by_depth.py drops the #class suffix); "changed" means a different sequence.

hits.tsv is rmblastn --outfmt "6 score qseqid qstart qend sseqid sstart send sstrand" for every query of
query.fa against db.fa (it may be the concatenation of several searches against parts of db.fa).  The filter is
rmblast-lib engine.rs apply_mask_level in cross-subject mode: one query's HSPs from all subjects sorted by
score DESC, oid DESC (oid = the subject's index in db.fa), subject start ASC, subject end DESC, query offset ASC,
query end DESC; an HSP is dropped when one accepted HSP of at least its score covers >= mask_level % of its
query span (minus-strand query coordinates shifted by +1, as NCBI's interval tree).  Survivors are printed in
input order.  Every per-subject search is independent of the other subjects, so filtering the union of separate
searches gives the hits of one search against the whole db.
"""
import sys
from bisect import bisect_left, insort


def names_and_lengths(fa):
    order, lens, name = [], {}, None
    with open(fa) as fh:
        for line in fh:
            if line.startswith('>'):
                name = line[1:].split()[0]
                order.append(name); lens[name] = 0
            else:
                lens[name] += len(line.strip())
    return order, lens


def read_fa(fa):
    """{key: (name, SEQUENCE)} with key = the name up to '#'."""
    recs, name, buf = {}, None, []
    def put():
        if name is not None: recs[name.split('#')[0]] = (name, ''.join(buf).upper())
    with open(fa) as fh:
        for line in fh:
            if line.startswith('>'):
                put(); name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    put()
    return recs


def unchanged(prev_fa, cur_fa):
    """{previous name: current name} for records whose sequence did not change."""
    prev, cur = read_fa(prev_fa), read_fa(cur_fa)
    return {prev[k][0]: v[0] for k, v in cur.items() if k in prev and prev[k][1] == v[1]}


def changed(prev_fa, cur_fa, out_fa):
    prev, cur = read_fa(prev_fa), read_fa(cur_fa)
    n = bp = 0
    with open(out_fa, 'w') as fh:
        for k, (name, seq) in cur.items():
            if k not in prev or prev[k][1] != seq:
                fh.write(f'>{name}\n{seq}\n'); n += 1; bp += len(seq)
    print(f'{n} of {len(cur)} records changed ({bp:,} of {sum(len(v[1]) for v in cur.values()):,} bp)', file=sys.stderr)


def carry(prev_fa, cur_fa, raw_tsv):
    ren = unchanged(prev_fa, cur_fa)
    out = sys.stdout
    with open(raw_tsv) as fh:
        for line in fh:
            f = line.split('\t')
            if f[4] in ren:
                f[4] = ren[f[4]]; out.write('\t'.join(f))


def load_matrix(path):
    """(scores {(query base, subject base): score}, base frequencies, lambda) as rmblast-lib reads a matrix file
    (matrix.rs: lambda by bisection from the FREQS line)."""
    import math
    freqs, cols, M = {}, None, {}
    for line in open(path):
        t = line.split()
        if not t: continue
        if line.lstrip().startswith('#'):
            if 'FREQS' in t:
                i = t.index('FREQS'); freqs = {t[j]: float(t[j + 1]) for j in range(i + 1, len(t) - 1, 2)}
            continue
        if cols is None: cols = t; continue
        for c, v in zip(cols, t[1:]): M[(t[0], c)] = int(v)
    ev = lambda lam: sum(freqs[x] * freqs[y] * math.exp(lam * M[(x, y)]) for x in 'ACGT' for y in 'ACGT')
    lo, lam = 0.0, 0.5
    while ev(lam) < 1.0: lo = lam; lam *= 2
    hi = lam
    while hi - lo > 0.00001:
        lam = (lo + hi) / 2
        if ev(lam) >= 1.0: hi = lam
        else: lo = lam
    return M, freqs, lam


class Scorer:
    """rmblastn's score of an alignment given as its two gapped rows (outfmt qseq, sseq): matrix scores, gaps
    24 + 6 per base, then the complexity adjustment of complexity_adjust.rs on the query bases of aligned pairs.
    Reproduces the reported score of 99.5% of rice depth-search hits (2026-11-12)."""
    def __init__(self, path, go=24, ge=6):
        import re, math, collections
        self.M, self.freqs, self.lam = load_matrix(path)
        self.go, self.ge, self.gap, self.log, self.Counter = go, ge, re.compile('-+'), math.log, collections.Counter
        self.lf = {k: math.log(v) for k, v in self.freqs.items() if v > 0 and v != 1.0}

    def col(self, a, b):                       # one column, a gap column as an extension
        return -self.ge if a == '-' or b == '-' else self.M.get((a, b), -1)

    def score(self, q, s):
        M, raw, cnt = self.M, 0, {}
        for (a, b), c in self.Counter(zip(q, s)).items():
            if a == '-' or b == '-': continue
            raw += c * M.get((a, b), -1)
            if a in self.lf: cnt[a] = cnt.get(a, 0) + c
        for x in (q, s):
            if '-' in x:
                for m in self.gap.finditer(x): raw -= self.go + self.ge * (m.end() - m.start())
        if not cnt: return raw
        log = self.log; tf = ts = tc = 0.0
        for a, c in cnt.items():
            tf += c * log(c); ts += c * self.lf[a]; tc += c
        tf -= tc * log(tc); ts -= tf
        return max(0, int(raw + ts / self.lam + 0.999))


def remap(cuts_tsv, raw_tsv, tol=30, min_score=225, clipj=1, matrix=None, q_fa=None, db_fa=None, trim_win=64, ext_max=500, xdrop=112):
    """An approximation of searching the cut records again (review idea 3.3, 2026-11-11).  A hit on a dropped
    record goes; a hit on a cut record loses what lay inside a cut (pieces under min_score dropped) and its
    subject coordinates shift left by the bases cut before them.  A piece that ends within tol bp of a cut and
    one that starts within tol bp after it, same query and strand and collinear on the query, are a copy without
    the insertion that was cut: joined into one hit, as a search of the cut record would align it (tol 0: no
    joining).
    Without a matrix (8-column hits): a piece's score and query span are the hit's, prorated by the subject span
    kept; a joined hit has the summed score less a gap.
    With a matrix (10-column hits, ... qseq sseq): a piece is the hit's alignment columns over the kept subject
    bases, scored as rmblastn scores them; an end made by a cut is trimmed back to its best-scoring point within
    trim_win columns; a joined hit is the two alignments end to end (unaligned bases between them as N, a length
    difference as a gap) and scored whole.  The pieces keep their alignment rows for the next round.
    With the query batch and the cut library as well (q_fa, db_fa): an end made by a cut is instead extended
    without gaps across the junction while that scores more than trimming (X-drop xdrop, at most ext_max
    columns), as the search would run on into the sequence now beside it, and the bases between two joined
    alignments are aligned with the length difference as one gap where it scores best."""
    cuts = {}
    with open(cuts_tsv) as fh:
        for line in fh:
            n, v = line.rstrip('\n').split('\t')
            cuts[n] = None if v == 'drop' else [tuple(map(int, x.split('-'))) for x in v.split(',')]
    junc = {}
    for n, runs in cuts.items():
        if runs:
            js, rem = [], 0
            for i, j in runs:
                js.append(i - rem); rem += j - i
            junc[n] = js
    out = sys.stdout
    stats = dict(kept=0, dropped=0, clipped=0, pieces=0, joined=0)
    S = Scorer(matrix) if matrix else None
    if S: stats.update(rescored=0, ratio=0, trimmed=0, extended=0)
    Q = {k: v[1] for k, v in read_fa(q_fa).items()} if q_fa else None
    D = {k: v[1] for k, v in read_fa(db_fa).items()} if db_fa else None
    comp = str.maketrans('ACGTRYKMSWN', 'TGCAYRMKSWN')

    # a piece: [score, qa, qb, sa, sb, left_ok, right_ok] (half-open; *_ok: that subject end may join), with a
    # matrix + [q row, s row, left end made by a cut, right end made by a cut, score factor, rescore]; the rows
    # run along the subject (a minus-strand hit's rows reversed, so its query coordinates descend along them)
    def trim_end(p, minus, right, qn=None, key=None):
        """An end of p made by a cut: drop the run of up to trim_win columns with the lowest (negative) score, or
        extend it without gaps into the sequence beside it if that gains more."""
        q, s = p[7], p[8]; n = len(q); best, bk, r = 0, 0, 0
        for k in range(1, min(trim_win, n - 1) + 1):
            i = n - k if right else k - 1
            r += S.col(q[i], s[i])
            j = i - 1 if right else i + 1
            if q[j] == '-' or s[j] == '-': continue
            if q[i] == '-' or s[i] == '-': r -= S.go          # a whole gap goes
            if r < best: best, bk = r, k
        if Q is not None:
            qs, ds = Q[qn], D[key]; up = right != minus       # the query ascends away from this end
            xs = ds[p[4]:p[4] + ext_max] if right else ds[max(0, p[3] - ext_max):p[3]][::-1]
            xq = qs[p[2]:p[2] + ext_max] if up else qs[max(0, p[1] - ext_max):p[1]][::-1]
            if minus: xs = xs.translate(comp)
            M, r, eb, ek, k = S.M, 0, 0, 0, 0
            for x, y in zip(xq, xs):
                k += 1; r += M.get((x, y), -1)
                if r > eb: eb, ek = r, k
                elif eb - r > xdrop: break
            if ek and eb > -best:
                if right: p[7], p[8], p[4] = q + xq[:ek], s + xs[:ek], p[4] + ek
                else: p[7], p[8], p[3] = xq[:ek][::-1] + q, xs[:ek][::-1] + s, p[3] - ek
                if up: p[2] += ek
                else: p[1] -= ek
                p[11] = True; stats['extended'] += 1; return
        if not bk: return
        dq, ds = (q[n - bk:], s[n - bk:]) if right else (q[:bk], s[:bk])
        nq, ns = len(dq) - dq.count('-'), len(ds) - ds.count('-')
        p[7], p[8] = (q[:n - bk], s[:n - bk]) if right else (q[bk:], s[bk:])
        if right: p[4] -= ns
        else: p[3] += ns
        if right != minus: p[2] -= nq
        else: p[1] += nq
        p[11] = True; stats['trimmed'] += 1

    def join_rows(a, b, minus, qn=None, key=None):
        """b's rows appended to a's (a left of b on the subject); False if nothing of b is left."""
        qb_, sb_, bsa = b[7], b[8], b[3]
        over = (b[2] - a[1]) if minus else (a[2] - b[1])      # query bases the two share
        k = 0
        while over > 0 and k < len(qb_):
            if qb_[k] != '-': over -= 1
            if sb_[k] != '-': bsa += 1
            k += 1
        while k < len(qb_) and (qb_[k] == '-' or sb_[k] == '-'):
            if qb_[k] != '-': over -= 1                       # over < 0: a query gap opens between the two
            if sb_[k] != '-': bsa += 1
            k += 1
        if k >= len(qb_): return False
        gq, gs = -over, bsa - a[4]
        m = min(gq, gs)
        if Q is not None and (gq or gs):                      # the bases between the two, the gap where it scores best
            xq = Q[qn][a[1] - gq:a[1]][::-1] if minus else Q[qn][a[2]:a[2] + gq]
            xs = D[key][a[4]:bsa].translate(comp) if minus else D[key][a[4]:bsa]
            d = abs(gq - gs); M = S.M
            d0 = [M.get(z, -1) for z in zip(xq, xs)]
            d1 = [M.get(z, -1) for z in (zip(xq[d:], xs) if gq > gs else zip(xq, xs[d:]))]
            best, at, t = None, 0, sum(d1)
            for i in range(m + 1):
                if best is None or t > best: best, at = t, i
                if i < m: t += d0[i] - d1[i]
            if not d: at = m
            fq = xq[:at] + (xq[at:at + d] if gq > gs else '-' * d) + xq[at + (d if gq > gs else 0):]
            fs = xs[:at] + (xs[at:at + d] if gs > gq else '-' * d) + xs[at + (d if gs > gq else 0):]
        else:
            fq = 'N' * m + ('N' * (gq - m) if gq > gs else '-' * (gs - m))
            fs = 'N' * m + ('N' * (gs - m) if gs > gq else '-' * (gq - m))
        a[7], a[8] = a[7] + fq + qb_[k:], a[8] + fs + sb_[k:]
        return True

    def flush(q, groups):
        for (key, minus), ps in groups.items():
            if tol > 0 and len(ps) > 1:
                for J in junc[key]:
                    left = [p for p in ps if p[6] and J - tol <= p[4] <= J]
                    right = [p for p in ps if p[5] and J <= p[3] <= J + tol]
                    if not left or not right: continue
                    used = set()
                    for a in sorted(left, key=lambda p: -p[0]):
                        best = None
                        for b in right:
                            if id(b) in used or b is a: continue
                            gs = b[3] - a[4]
                            gq = (a[1] - b[2]) if minus else (b[1] - a[2])
                            d = abs(gq - gs)
                            if gq >= -tol and d <= tol and (best is None or d < best[0]): best = (d, b)
                        if best:
                            d, b = best
                            if S:
                                if not join_rows(a, b, minus, q, key): continue
                                a[10], a[11], a[12] = b[10], True, 1.0
                            else:
                                a[0] = a[0] + b[0] - (24 + 6 * d if d else 0)
                            used.add(id(b))
                            a[1], a[2] = min(a[1], b[1]), max(a[2], b[2])
                            a[4], a[6] = b[4], b[6]
                            ps.remove(b); stats['joined'] += 1
            for p in ps:
                if S:
                    if p[9]: trim_end(p, minus, False, q, key)
                    if p[10]: trim_end(p, minus, True, q, key)
                    if p[11]: p[0] = int(S.score(p[7], p[8]) * p[12] + 0.5); stats['rescored'] += 1
                sc, qa, qb, sa, sb = p[:5]
                if sc < min_score or qb <= qa or sb <= sa: continue
                a, b = (sb, sa + 1) if minus else (sa + 1, sb)
                rows = ('\t' + p[7][::-1] + '\t' + p[8][::-1] if minus else '\t' + p[7] + '\t' + p[8]) if S else ''
                out.write(f'{sc}\t{q}\t{qa + 1}\t{qb}\t{key}\t{a}\t{b}\t{"minus" if minus else "plus"}{rows}\n'); stats['pieces'] += 1

    cur_q, groups = None, {}
    with open(raw_tsv) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if f[1] != cur_q:                   # every query's hits stay together in the output
                if groups: flush(cur_q, groups)
                cur_q, groups = f[1], {}
            key = f[4].split('#')[0]
            if key not in cuts:
                if f[4] != key: f[4] = key; line = '\t'.join(f) + '\n'
                out.write(line); stats['kept'] += 1; continue
            runs = cuts[key]
            if runs is None:
                stats['dropped'] += 1; continue
            score, qa, qb, a, b = int(f[0]), int(f[2]) - 1, int(f[3]), int(f[5]), int(f[6])
            minus = f[7] == 'minus'
            lo, hi = min(a, b) - 1, max(a, b)
            span, qspan = hi - lo, qb - qa
            pos, rem, pieces = lo, 0, []       # rem: bases cut left of pos
            lclip = 0                           # how far the hit ran into the cut left of the current piece
            for i, j in runs:
                if j <= pos:
                    rem += j - i; continue
                if i >= hi: break
                if i > pos:
                    pieces.append((pos, i, rem, lclip, min(hi, j) - i))
                lclip = min(hi, j) - max(lo, i)
                rem += j - i; pos = j
            if pos < hi: pieces.append((pos, hi, rem, lclip, 0))
            whole = len(pieces) == 1 and pieces[0][:2] == (lo, hi)
            if not whole: stats['clipped'] += 1
            g = groups.setdefault((key, minus), [])
            if S:
                qr, sr = (f[8][::-1], f[9][::-1]) if minus else (f[8], f[9])
                if whole:
                    g.append([score, qa, qb, lo - pieces[0][2], hi - pieces[0][2], True, True, qr, sr, False, False, False, 1.0])
                    continue
                full = S.score(qr, sr)
                fac = 1.0 if full == score or full <= 0 else score / full   # a score rmblastn computed otherwise
                if fac != 1.0: stats['ratio'] += 1
                idx = None if '-' not in sr else [i for i, c in enumerate(sr) if c != '-']
                n = len(qr)
                for p0, p1, rm, lc, rc in pieces:
                    c0, c1 = (p0 - lo, p1 - lo) if idx is None else (idx[p0 - lo], idx[p1 - lo - 1] + 1)
                    while c0 < c1 and (qr[c0] == '-' or sr[c0] == '-'):      # no gap at a piece's end
                        if sr[c0] != '-': p0 += 1
                        c0 += 1
                    while c1 > c0 and (qr[c1 - 1] == '-' or sr[c1 - 1] == '-'):
                        if sr[c1 - 1] != '-': p1 -= 1
                        c1 -= 1
                    if c1 <= c0: continue
                    before = c0 - qr.count('-', 0, c0); inside = c1 - c0 - qr.count('-', c0, c1)
                    nqa, nqb = (qb - before - inside, qb - before) if minus else (qa + before, qa + before + inside)
                    g.append([0, nqa, nqb, p0 - rm, p1 - rm, clipj or lc <= tol, clipj or rc <= tol,
                              qr[c0:c1], sr[c0:c1], c0 > 0, c1 < n, False, fac])
                    g[-1][0] = int(S.score(g[-1][7], g[-1][8]) * fac + 0.5)   # set again if an end is trimmed or joined
                continue
            for p0, p1, rm, lc, rc in pieces:
                t0, t1 = (p0 - lo) / span, (p1 - lo) / span
                sc = score if (p0, p1) == (lo, hi) else int(score * (p1 - p0) / span)
                if minus: nqa, nqb = qb - round(qspan * t1), qb - round(qspan * t0)
                else: nqa, nqb = qa + round(qspan * t0), qa + round(qspan * t1)
                g.append([sc, nqa, nqb, p0 - rm, p1 - rm, clipj or lc <= tol, clipj or rc <= tol])
    if groups: flush(cur_q, groups)
    print('remap: ' + ', '.join(f'{v:,} {k}' for k, v in stats.items()), file=sys.stderr)


def main():
    mode = sys.argv[1]
    if mode == 'changed':
        return changed(*sys.argv[2:5])
    if mode == 'carry':
        return carry(*sys.argv[2:5])
    if mode == 'remap':
        return remap(sys.argv[2], sys.argv[3], *map(int, sys.argv[4:7]), *sys.argv[7:10])
    assert mode == 'filter', __doc__
    hits_f, db_f, q_f, level = sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
    oid = {n: i for i, n in enumerate(names_and_lengths(db_f)[0])}
    qlen = names_and_lengths(q_f)[1]
    by_q = {}
    with open(hits_f) as fh:
        for i, line in enumerate(fh):
            f = line.split('\t')
            score, q, qs1, qe, s, a, b, strand = int(f[0]), f[1], int(f[2]), int(f[3]), f[4], int(f[5]), int(f[6]), f[7].strip()
            qs = qs1 - 1
            ss, se = min(a, b) - 1, max(a, b)
            minus = strand == 'minus'
            L = qlen[q]
            qoff, qend = (L - qe, L - qs) if minus else (qs, qe)
            key = (-score, -oid[s], ss, -se, qoff, -qend)
            sq, eq = (qs + 1, qe + 1) if minus else (qs, qe)
            by_q.setdefault(q, []).append((key, sq, eq, score, i))
    keep = set()
    if level >= 100:
        keep = None
    else:
        for q, hs in by_q.items():
            hs.sort(key=lambda h: h[0])
            acc, starts, maxspan = [], [], 0
            for _, qs, qe, sc, i in hs:
                span = qe - qs
                if span == 0:
                    continue
                upper = bisect_left(starts, qe)
                lo = qs - maxspan
                masked = False
                j = upper
                while j > 0:
                    j -= 1
                    qsj, qej, scj = acc[j]
                    if qsj < lo:
                        break
                    if qej <= qs or scj < sc:
                        continue
                    ovlp = min(qej, qe) - max(qsj, qs)
                    if ovlp * 100 // span >= level:
                        masked = True; break
                if not masked:
                    p = bisect_left(starts, qs)
                    starts.insert(p, qs); acc.insert(p, (qs, qe, sc))
                    maxspan = max(maxspan, span)
                    keep.add(i)
    out = sys.stdout
    with open(hits_f) as fh:
        for i, line in enumerate(fh):
            if keep is None or i in keep:
                if line.count('\t') > 7: line = '\t'.join(line.split('\t')[:8]) + '\n'   # no alignment rows
                out.write(line)


if __name__ == '__main__':
    main()
