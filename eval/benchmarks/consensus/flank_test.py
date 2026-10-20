#!/usr/bin/env python3
"""2026-10-19: segmental-duplication (flank) test for element clusters.
usage: flank_test.py write <run_dir> <genome.fa> <out_dir> [flank=1000]   -> out_dir/flanks.fa (left and right flanks of
                                                                            every element in a multi-element cluster)
       flank_test.py score <run_dir> <out_dir> <flanks.tsv>              -> per cluster: member pairs whose left AND right
                                                                            flanks share >= 300 bp at >= 80% identity
flanks.tsv: rmblastn all-vs-all of flanks.fa (qseqid sseqid pident length ...). A cluster is SD-like when more than
half of its member pairs share both flanks (elements of a real family insert at unrelated sites)."""
import sys, os, json, collections, itertools
mode, run = sys.argv[1:3]
cl = json.load(open(os.path.join(run, 'elem', 'clusters_nested.json')))
mem = collections.defaultdict(list)
for e, c in cl.items(): mem[c].append(e)
multi = {c: v for c, v in mem.items() if len(v) > 1}
if mode == 'write':
    g, out = sys.argv[3:5]; fl = int(sys.argv[5]) if len(sys.argv) > 5 else 1000
    fai = {l.split()[0]: list(map(int, l.split()[1:4])) for l in open(g + '.fai')}
    fa = open(g, 'rb')
    def get(s, a, b):
        L, off, lb = fai[s][:3]; a = max(1, a); b = min(L, b)
        if b < a: return ''
        sa = off + (a - 1) // lb * (lb + 1) + (a - 1) % lb; sb = off + (b - 1) // lb * (lb + 1) + (b - 1) % lb
        fa.seek(sa); return fa.read(sb - sa + 1).decode().replace('\n', '')
    with open(os.path.join(out, 'flanks.fa'), 'w') as fh:
        for v in multi.values():
            for e in v:
                s, r = e.rsplit(':', 1); a, b = map(int, r.split('..'))
                for side, seq in (('L', get(s, a - fl, a - 1)), ('R', get(s, b + 1, b + fl))):
                    if len(seq) >= 300: fh.write(f'>{e}|{side}\n{seq}\n')
else:
    out, tsv = sys.argv[3:5]
    share = collections.defaultdict(int)  # (a, b, side) -> aligned bp
    for l in open(tsv):
        q, s, pid, ln = l.split('\t')[:4]
        q = q.split('#')[0]; s = s.split('#')[0]
        (qe, qs), (se, ss) = q.rsplit('|', 1), s.rsplit('|', 1)
        if qe == se or qs != ss or float(pid) < 80 or cl.get(qe) != cl.get(se): continue
        k = (min(qe, se), max(qe, se), qs); share[k] = max(share[k], int(ln))
    rows = []
    for c, v in multi.items():
        n = both = 0
        for a, b in itertools.combinations(sorted(v), 2):
            n += 1; both += share.get((a, b, 'L'), 0) >= 300 and share.get((a, b, 'R'), 0) >= 300
        rows.append((c, len(v), both, n))
    with open(os.path.join(out, 'flank_clusters.tsv'), 'w') as fh:
        fh.write('cluster\tmembers\tpairs_sharing_both_flanks\tpairs\n')
        for r in rows: fh.write('\t'.join(map(str, r)) + '\n')
    sd = [r for r in rows if r[2] > r[3] / 2]
    print(f'{len(rows)} multi-element clusters; SD-like (> half of pairs share both flanks): {len(sd)} holding {sum(r[1] for r in sd)} elements')
