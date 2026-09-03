#!/usr/bin/env python3
"""Compare two LTRharvest-format prediction files, sequence by sequence.

Both inputs use the GenomeTools tabular layout (`s(ret) e(ret) l(ret)
s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim seq-nr`), which is what
`gt ltrharvest` and `ltr-forage --format gt` write.  Predictions are matched
within a sequence by best reciprocal overlap (at least half of the longer
element) and each matched pair is classified by how far its six coordinates
differ.

With --fasta, unmatched predictions are also checked for N inside either LTR,
which needs a samtools .fai index next to the FASTA.  With --out-prefix, the
unmatched predictions of each input are written to PREFIX.a_only.tsv and
PREFIX.b_only.tsv in the 8-column form rm_truth_recall.py accepts.

Usage:
    compare_predictions.py A.out B.out [--fasta genome.fa] [--out-prefix P]
"""
import argparse
import collections
import statistics


def parse(path):
    by_seq = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            c = line.split()
            if len(c) < 11 or not c[0].isdigit():
                continue
            e = dict(s=int(c[0]), e=int(c[1]), l1s=int(c[3]), l1e=int(c[4]),
                     l2s=int(c[6]), l2e=int(c[7]), sim=float(c[9]), seq=int(c[10]))
            by_seq[e['seq']].append(e)
    for v in by_seq.values():
        v.sort(key=lambda x: x['s'])
    return by_seq


class Fasta:
    """Random access into a FASTA through its samtools .fai index."""

    def __init__(self, path):
        self.fh = open(path, 'rb')
        self.idx = []
        with open(path + '.fai') as fai:
            for line in fai:
                name, ln, off, lb, lw = line.split()[:5]
                self.idx.append((name, int(ln), int(off), int(lb), int(lw)))

    def fetch(self, seqnum, start1, end1):
        _, ln, off, lb, lw = self.idx[seqnum]
        s0, e0 = start1 - 1, min(end1, ln)
        b0 = off + (s0 // lb) * lw + s0 % lb
        b1 = off + ((e0 - 1) // lb) * lw + (e0 - 1) % lb + 1
        self.fh.seek(b0)
        return self.fh.read(b1 - b0).replace(b'\n', b'').replace(b'\r', b'')


def match(a, b, threshold=0.5):
    pairs, a_only, b_only = [], [], []
    for seq in sorted(set(a) | set(b)):
        A, B = a.get(seq, []), b.get(seq, [])
        used, j0 = set(), 0
        for x in A:
            best, bi = -1.0, -1
            while j0 < len(B) and B[j0]['e'] < x['s'] - 20000:
                j0 += 1
            j = j0
            while j < len(B) and B[j]['s'] <= x['e']:
                y = B[j]
                if j not in used and y['e'] >= x['s']:
                    ov = min(x['e'], y['e']) - max(x['s'], y['s']) + 1
                    r = ov / max(x['e'] - x['s'] + 1, y['e'] - y['s'] + 1)
                    if r > best:
                        best, bi = r, j
                j += 1
            if bi >= 0 and best >= threshold:
                used.add(bi)
                pairs.append((x, B[bi]))
            else:
                a_only.append(x)
        b_only += [B[j] for j in range(len(B)) if j not in used]
    return pairs, a_only, b_only


def maxdiff(x, y):
    return max(abs(x[k] - y[k]) for k in ('s', 'e', 'l1e', 'l2s'))


def write_tsv(path, items):
    with open(path, 'w') as fh:
        for x in items:
            fh.write('\t'.join(str(x[k]) for k in ('seq', 's', 'e', 'l1s', 'l1e', 'l2s', 'l2e', 'sim')) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('a')
    ap.add_argument('b')
    ap.add_argument('--fasta', help='genome FASTA with a samtools .fai index, for N checks')
    ap.add_argument('--out-prefix', help='write unmatched predictions to PREFIX.a_only.tsv / PREFIX.b_only.tsv')
    args = ap.parse_args()

    a, b = parse(args.a), parse(args.b)
    na, nb = sum(map(len, a.values())), sum(map(len, b.values()))
    print(f'A={args.a}: {na}   B={args.b}: {nb}')
    pairs, a_only, b_only = match(a, b)

    same_outer = lambda x, y: x['s'] == y['s'] and x['e'] == y['e']
    exact = sum(1 for x, y in pairs if maxdiff(x, y) == 0)
    exact_sim = sum(1 for x, y in pairs if maxdiff(x, y) == 0 and abs(x['sim'] - y['sim']) < 0.005)
    outer = sum(1 for x, y in pairs if same_outer(x, y) and maxdiff(x, y) > 0)
    small = sum(1 for x, y in pairs if 0 < maxdiff(x, y) <= 20 and not same_outer(x, y))
    med = sum(1 for x, y in pairs if 20 < maxdiff(x, y) <= 200 and not same_outer(x, y))
    large = sum(1 for x, y in pairs if maxdiff(x, y) > 200 and not same_outer(x, y))
    print(f'matched pairs: {len(pairs)}  identical in all coordinates: {exact} (with identical similarity: {exact_sim})')
    print(f'  same outer boundaries, different inner boundaries: {outer}')
    print(f'  differ by <=20 bp: {small}   21-200 bp: {med}   >200 bp: {large}')
    print(f'A-only: {len(a_only)}   B-only: {len(b_only)}')

    if args.fasta:
        fa = Fasta(args.fasta)
        for lst, label in ((a_only, 'A-only'), (b_only, 'B-only')):
            with_n = 0
            for x in lst:
                l1 = fa.fetch(x['seq'], x['l1s'], x['l1e'])
                l2 = fa.fetch(x['seq'], x['l2s'], x['l2e'])
                if b'N' in l1.upper() or b'N' in l2.upper():
                    with_n += 1
            if lst:
                sims = [x['sim'] for x in lst]
                lens = [x['e'] - x['s'] + 1 for x in lst]
                print(f'{label}: {len(lst)} predictions, {with_n} with N inside an LTR; '
                      f'similarity median {statistics.median(sims):.2f}, element length median {statistics.median(lens):.0f}')

    if args.out_prefix:
        write_tsv(args.out_prefix + '.a_only.tsv', a_only)
        write_tsv(args.out_prefix + '.b_only.tsv', b_only)


if __name__ == '__main__':
    main()
