#!/usr/bin/env python3
"""Label each prediction by the RepeatMasker class under its two LTRs.

For every prediction, the dominant repeat class covering the left LTR and the
right LTR is found (Simple_repeat, Low_complexity and Satellite are pooled;
an LTR with under half its length annotated is "mostly-unannotated").  The
table shows how many predictions have both LTRs in each class, and how many
of the LTR-class ones have the same family name on both sides.  This is the
closest thing to a precision estimate available without curated truth.

Usage:
    rm_classify.py genome.fa.out genome.fa.fai predictions.out
"""
import bisect
import collections
import sys


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    rm_path, fai_path, pred_path = sys.argv[1:]
    names = [l.split()[0] for l in open(fai_path)]
    ann = collections.defaultdict(list)
    with open(rm_path) as fh:
        for line in fh:
            c = line.split()
            if len(c) < 15 or not c[0].isdigit():
                continue
            cls = c[10].split('/')[0]
            if cls in ('Simple_repeat', 'Low_complexity', 'Satellite'):
                cls = 'Simple/Sat'
            ann[c[4]].append((int(c[5]), int(c[6]), cls, c[9]))
    for v in ann.values():
        v.sort()
    starts = {c: [a[0] for a in v] for c, v in ann.items()}

    def dominant(chrom, s, e):
        v = ann.get(chrom, [])
        k = bisect.bisect_left(starts.get(chrom, []), s - 20000)
        cov, name = collections.Counter(), collections.defaultdict(collections.Counter)
        while k < len(v) and v[k][0] <= e:
            a = v[k]
            o = min(e, a[1]) - max(s, a[0]) + 1
            if o > 0:
                cov[a[2]] += o
                name[a[2]][a[3]] += o
            k += 1
        if not cov:
            return 'unannotated', ''
        cls, o = cov.most_common(1)[0]
        if o < 0.5 * (e - s + 1):
            return 'mostly-unannotated', ''
        return cls, name[cls].most_common(1)[0][0]   # the family covering most of it (was the last one overlapping)

    table, n, same_family = collections.Counter(), 0, 0
    with open(pred_path) as fh:
        for line in fh:
            c = line.split()
            if len(c) < 11 or not c[0].isdigit():
                continue
            chrom = names[int(c[10])]
            c1, n1 = dominant(chrom, int(c[3]), int(c[4]))
            c2, n2 = dominant(chrom, int(c[6]), int(c[7]))
            table[c1 if c1 == c2 else 'mixed(' + '/'.join(sorted({c1, c2})) + ')'] += 1
            n += 1
            if c1 == c2 == 'LTR' and n1 == n2:
                same_family += 1
    print(f'{pred_path}: {n} predictions; dominant RepeatMasker class under (left LTR, right LTR):')
    for k, v in table.most_common():
        print(f'  {k:<40} {v:>6}  {100 * v / n:5.1f}%')
    print(f'  (LTR class with the same family on both sides: {same_family})')


if __name__ == '__main__':
    main()
