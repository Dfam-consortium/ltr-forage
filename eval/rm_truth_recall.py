#!/usr/bin/env python3
"""Score LTR predictions against intact elements in a RepeatMasker annotation.

An intact element is an LTR-class annotation followed on the same strand by
an internal-region annotation (name ending in -int, _int, -I, _I or I) and
then the nearest annotation of the same LTR family, with both LTRs 100-1000 bp
long and their 5' ends 1000-15000 bp apart, the size window the tools were run
with.  A prediction hits an element when each predicted LTR overlaps the
corresponding annotated LTR by at least half of the shorter of the two.

Recall is reported by the larger of the two LTRs' divergence from consensus
(RepeatMasker's "perc div" column), because LTRharvest only finds pairs whose
LTRs are still similar to each other.

Prediction files are LTRharvest-format outputs (sequence numbers resolved
through the FASTA .fai, whose order matches the index) or the 8-column TSVs
written by compare_predictions.py --out-prefix.

Usage:
    rm_truth_recall.py genome.fa.out genome.fa.fai PRED [PRED ...]
"""
import bisect
import collections
import sys

INTERNAL_SUFFIXES = ('-int', '_int', '-I', '_I', 'I')
BINS = [(0, 5), (5, 10), (10, 15), (15, 20), (20, 100)]


def is_internal(name):
    return name.endswith(INTERNAL_SUFFIXES)


def load_truth(rm_path):
    ann = collections.defaultdict(list)
    with open(rm_path) as fh:
        for line in fh:
            c = line.split()
            if len(c) < 15 or not c[0].isdigit() or not c[10].startswith('LTR'):
                continue
            ann[c[4]].append((int(c[5]), int(c[6]), c[8], c[9], float(c[1])))
    truth = []
    for chrom, lst in ann.items():
        lst.sort()
        for i, (s1, e1, st1, nm1, dv1) in enumerate(lst):
            if is_internal(nm1) or not 100 <= e1 - s1 + 1 <= 1000:
                continue
            seen_int = False
            for s2, e2, st2, nm2, dv2 in lst[i + 1:]:
                if s2 - s1 > 15000:
                    break
                if nm2 == nm1 and st2 == st1:
                    if seen_int and s2 - s1 >= 1000 and 100 <= e2 - s2 + 1 <= 1000:
                        truth.append((chrom, s1, e1, s2, e2, max(dv1, dv2), nm1))
                    break
                if nm2 != nm1 and is_internal(nm2) and st2 == st1:
                    seen_int = True
    return truth


def bin_of(d):
    for lo, hi in BINS:
        if lo <= d < hi:
            return (lo, hi)
    return BINS[-1]


def overlap_frac(a1, a2, b1, b2):
    o = min(a2, b2) - max(a1, b1) + 1
    return o / min(a2 - a1 + 1, b2 - b1 + 1) if o > 0 else 0.0


def load_predictions(path, names):
    out = []
    with open(path) as fh:
        for line in fh:
            c = line.split()
            if len(c) >= 11 and c[0].isdigit():
                out.append((names[int(c[10])], int(c[3]), int(c[4]), int(c[6]), int(c[7])))
            elif len(c) == 8 and c[0].isdigit():
                out.append((names[int(c[0])], int(c[3]), int(c[4]), int(c[5]), int(c[6])))
    return out


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    rm_path, fai_path, preds = sys.argv[1], sys.argv[2], sys.argv[3:]
    names = [l.split()[0] for l in open(fai_path)]
    truth = load_truth(rm_path)
    print(f'intact elements in {rm_path}: {len(truth)}')
    total = collections.Counter(bin_of(t[5]) for t in truth)
    by_chrom = collections.defaultdict(list)
    for t in truth:
        by_chrom[t[0]].append(t)
    for v in by_chrom.values():
        v.sort(key=lambda t: t[1])
    starts = {c: [t[1] for t in v] for c, v in by_chrom.items()}

    hits = {}
    for path in preds:
        pred = load_predictions(path, names)
        hit = set()
        supported = 0
        for chrom, l1s, l1e, l2s, l2e in pred:
            v = by_chrom.get(chrom, [])
            k = bisect.bisect_left(starts.get(chrom, []), l1s - 2000)
            found = False
            while k < len(v) and v[k][1] <= l1e + 2000:
                t = v[k]
                if overlap_frac(l1s, l1e, t[1], t[2]) >= 0.5 and overlap_frac(l2s, l2e, t[3], t[4]) >= 0.5:
                    hit.add((chrom, t[1], t[3]))
                    found = True
                k += 1
            supported += found
        got = collections.Counter(bin_of(t[5]) for t in truth if (t[0], t[1], t[3]) in hit)
        print(f'\n{path}: {len(pred)} predictions, {supported} ({100 * supported / max(1, len(pred)):.1f}%) hit an intact element')
        print('  recall by max LTR divergence from consensus:')
        for b in BINS:
            print(f'    {b[0]:>2}-{b[1]:<3}%: {got[b]:>6} / {total[b]:<6} = {100 * got[b] / max(1, total[b]):5.1f}%')
        print(f'  elements hit: {len(hit)}')
        hits[path] = hit
    if len(preds) >= 2:
        a, b = hits[preds[0]], hits[preds[1]]
        print(f'\nhit only by {preds[0]}: {len(a - b)}; only by {preds[1]}: {len(b - a)}; by both: {len(a & b)}')


if __name__ == '__main__':
    main()
