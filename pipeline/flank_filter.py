#!/usr/bin/env python3
"""Flag element clusters whose copies are segmental duplications: copies of a family insert at unrelated
sites, so their flanks differ; copies of a duplicated block share flanking sequence on both sides.

usage: flank_filter.py write <clusters.json> <genome.fa> <flanks.fa> [--flank 1000]
       flank_filter.py bins  <clusters.json> <flanks.fa> <bin dir> [--bin-elements 1000]
       flank_filter.py score <clusters.json> <flanks.tsv> <sd_clusters.tsv> [--min-len 300] [--min-id 80]

write: the --flank bp left and right of every element in a cluster of >= 2 elements (needs <genome.fa>.fai;
       written if missing).  Search flanks.fa against itself (rmblastn) between the two calls.
bins:  split flanks.fa into files of whole clusters (about --bin-elements elements each) and one side,
       <bin dir>/b<N>.<L|R>.fa, each searched against itself.  score reads only same-cluster, same-side pairs,
       so the calls match one search of everything; it saves the all-vs-all of unrelated clusters (maize
       689 -> 48 min and 44 -> 1 GB of hits, soybean 33 -> 1.4 min, mouse 11 -> 2.1 min; identical calls on
       seven genomes, 2026-11-05).
score: a member pair shares a flank when its same-side flanks align over >= --min-len bp at >= --min-id %;
       a cluster is 'sd' when more than half of its member pairs share both flanks, or when every member
       shares both flanks with some other member (no independent insertion at all).  Writes one row per
       multi-element cluster: cluster, members, pairs sharing both flanks, pairs, call (sd|ok).
The second test catches a cluster made of several duplication groups, which dilute the pair vote (zebra
finch, stage B: a 6-copy tandem array plus two duplicated pairs in one cluster, 19 of 45 pairs; every other
cluster on zebra finch and rice has at most 88% of its members flank-paired).
On zebra finch (2026-10-19) 8 of 74 clusters were flagged, holding 9.7 of the 15.4 Mb masked outside the
curated LTR annotation (one 19-copy tandem array of a 13.5 kb block, flanks 99.6% identical); rice: 4 of
218, all two-copy pairs.  Single-element clusters cannot be tested this way.
Flanks are compared same side to same side, so a block duplicated in inverted orientation (its left flank
the reverse complement of the other copy's right flank) is not flagged; none has been seen.
"""
import argparse, collections, itertools, json, os


def fai(genome):
    path = genome + '.fai'
    if not os.path.exists(path):
        with open(genome, 'rb') as fh, open(path, 'w') as out:
            name = None; off = length = lb = lw = 0
            pos = 0
            for line in fh:
                if line.startswith(b'>'):
                    if name: out.write(f'{name}\t{length}\t{off}\t{lb}\t{lw}\n')
                    name = line[1:].split()[0].decode(); length = lb = lw = 0; off = pos + len(line)
                elif name:
                    if lb == 0: lb, lw = len(line.rstrip(b'\r\n')), len(line)
                    length += len(line.rstrip(b'\r\n'))
                pos += len(line)
            if name: out.write(f'{name}\t{length}\t{off}\t{lb}\t{lw}\n')
    return {l.split()[0]: list(map(int, l.split()[1:5])) for l in open(path)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('mode', choices=('write', 'bins', 'score'))
    ap.add_argument('clusters'); ap.add_argument('a'); ap.add_argument('b')
    ap.add_argument('--flank', type=int, default=1000)
    ap.add_argument('--bin-elements', type=int, default=1000)
    ap.add_argument('--min-len', type=int, default=300)
    ap.add_argument('--min-id', type=float, default=80)
    a = ap.parse_args()
    cl = json.load(open(a.clusters))
    mem = collections.defaultdict(list)
    for e, c in cl.items(): mem[c].append(e)
    multi = {c: v for c, v in mem.items() if len(v) > 1}

    if a.mode == 'write':
        idx = fai(a.a); fa = open(a.a, 'rb')
        def get(s, x, y):
            L, off, lb, lw = idx[s]; x = max(1, x); y = min(L, y)
            if y < x: return ''
            sa = off + (x - 1) // lb * lw + (x - 1) % lb; sb = off + (y - 1) // lb * lw + (y - 1) % lb
            fa.seek(sa); return fa.read(sb - sa + 1).decode().replace('\n', '').replace('\r', '')
        n = 0
        with open(a.b, 'w') as fh:
            for v in multi.values():
                for e in v:
                    s, r = e.rsplit(':', 1); x, y = map(int, r.split('..'))
                    for side, seq in (('L', get(s, x - a.flank, x - 1)), ('R', get(s, y + 1, y + a.flank))):
                        if len(seq) >= a.min_len: fh.write(f'>{e}|{side}\n{seq}\n'); n += 1
        print(f'{n} flanks from {sum(len(v) for v in multi.values())} elements in {len(multi)} multi-element clusters')
        return

    if a.mode == 'bins':
        os.makedirs(a.b, exist_ok=True)
        binof, nb, cur = {}, 0, 0   # largest clusters first; a cluster bigger than a bin gets one of its own
        for c, v in sorted(multi.items(), key=lambda x: (-len(x[1]), str(x[0]))):
            if cur and cur + len(v) > a.bin_elements: nb += 1; cur = 0
            binof[c] = nb; cur += len(v)
        fh = {}
        w = None
        for l in open(a.a):
            if l.startswith('>'):
                e, side = l[1:].split()[0].rsplit('|', 1)
                k = (binof[cl[e]], side)
                if k not in fh: fh[k] = open(os.path.join(a.b, f'b{k[0]}.{k[1]}.fa'), 'w')
                w = fh[k]
            w.write(l)
        for f in fh.values(): f.close()
        print(f'{len(fh)} bin files ({nb + 1 if multi else 0} bins of whole clusters, by side)')
        return

    share = collections.defaultdict(int)
    for l in open(a.a):
        q, s, pid, ln = l.split('\t')[:4]
        q = q.split('#')[0]; s = s.split('#')[0]
        (qe, qs), (se, ss) = q.rsplit('|', 1), s.rsplit('|', 1)
        if qe == se or qs != ss or float(pid) < a.min_id or cl.get(qe) != cl.get(se): continue
        k = (min(qe, se), max(qe, se), qs); share[k] = max(share[k], int(ln))
    nsd = esd = 0
    with open(a.b, 'w') as fh:
        fh.write('cluster\tmembers\tpairs_sharing_both_flanks\tpairs\tcall\n')
        for c, v in multi.items():
            n = both = 0; paired = set()
            for x, y in itertools.combinations(sorted(v), 2):
                n += 1
                if share.get((x, y, 'L'), 0) >= a.min_len and share.get((x, y, 'R'), 0) >= a.min_len:
                    both += 1; paired |= {x, y}
            call = 'sd' if both > n / 2 or len(paired) == len(v) else 'ok'
            if call == 'sd': nsd += 1; esd += len(v)
            fh.write(f'{c}\t{len(v)}\t{both}\t{n}\t{call}\n')
    print(f'{len(multi)} multi-element clusters; segmental-duplication-like: {nsd} holding {esd} elements')


if __name__ == '__main__':
    main()
