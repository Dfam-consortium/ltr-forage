#!/usr/bin/env python3
"""Stage B helpers for ltr_library.py --stage-b: depth cleanup and re-clustering only where they are needed.

usage: stage_b.py split  <clusters_raw.json> <lib.fa> <reps.fa> <depth_lib.fa> <reps.json> [--min-cons 3]
       stage_b.py expand <clusters_raw.json> <reps.json> <clusters_red.json> <clean_nested.fa>
                         <clusters_nested.json> <exemplar.fa>
(env TSV: ltr-forage's cand.tsv, as step4.py)

The raw library is clustered first (step4.py on its all-vs-all).  A cluster of >= --min-cons elements is
represented from then on by its te-composer consensus (reps.fa, consensus_models.py --big-only), whose
insertion packing and column gate already drop the minority passengers depth cleanup exists to remove.
split:  the records depth cleanup and the second clustering work on: every record of the elements in
        smaller clusters, plus the consensus records of each big cluster (named after its exemplar
        element); reps.json maps each representative element to its raw cluster.
expand: step4.py on the reduced set gives clusters_red.json; every member of a big raw cluster takes
        the cluster of its representative.  Exemplar per cluster: the representative of the largest big
        raw cluster in it, else the member with the longest internal (exemplar_lib.py's rule).
"""
import argparse, collections, csv, json, os


def elements():
    recs = {}
    for r in csv.DictReader(open(os.environ['TSV']), delimiter='\t'):
        if r['keep'] != '1': continue
        s = r['seq']; e = f"{s}:{r['elem_start']}..{r['elem_end']}"
        ls, le, rs, re_ = int(r['lltr_start']), int(r['lltr_end']), int(r['rltr_start']), int(r['rltr_end'])
        recs[e] = [f'{s}:{ls}..{le}_LTR', f'{s}:{le+1}..{rs-1}_INT', f'{s}:{rs}..{re_}_LTR']
    return recs


def read_fa(path):
    seqs, names, n = {}, {}, None
    for l in open(path):
        if l.startswith('>'):
            full = l[1:].split()[0]; n = full.split('#')[0]; names[n] = full; seqs[n] = []
        else:
            seqs[n].append(l.strip())
    return {k: ''.join(v) for k, v in seqs.items()}, names


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('mode', choices=('split', 'expand'))
    ap.add_argument('files', nargs='+')
    ap.add_argument('--min-cons', type=int, default=3)
    a = ap.parse_args()
    erecs = elements()
    rec2elem = {r: e for e, v in erecs.items() for r in v}
    raw = json.load(open(a.files[0]))
    members = collections.defaultdict(list)
    for e, c in raw.items(): members[c].append(e)

    if a.mode == 'split':
        lib, reps_fa, out_fa, out_json = a.files[1:5]
        seqs, names = read_fa(lib)
        rseqs, rnames = read_fa(reps_fa)
        reps = {}
        for r in rseqs:
            e = rec2elem[r]; reps[e] = raw[e]
        big = {c for c in reps.values()}
        n = collections.Counter()
        with open(out_fa, 'w') as fh:
            for c, v in members.items():
                if c in big: continue
                for e in v:
                    for r in erecs[e]:
                        if r in seqs: fh.write(f'>{names[r]}\n{seqs[r]}\n'); n['small records'] += 1
            for r, s in rseqs.items():
                fh.write(f'>{rnames[r]}\n{s}\n'); n['representative records'] += 1
        json.dump(reps, open(out_json, 'w'))
        nb = sum(len(members[c]) for c in big)
        print(f'{len(members)} raw clusters; {len(big)} of >= {a.min_cons} elements ({nb} elements) represented by '
              f'their consensus; {len(raw) - nb} elements in smaller clusters; ' + ', '.join(f'{k} {v}' for k, v in n.items()))
        return

    reps_json, red_json, clean_fa, out_json, out_fa = a.files[1:6]
    reps = json.load(open(reps_json)); red = json.load(open(red_json))
    rep_of = {c: e for e, c in reps.items()}
    final = {}
    for e, c in raw.items():
        final[e] = red[rep_of[c]] if c in rep_of else red[e]
    json.dump(final, open(out_json, 'w'))
    seqs, names = read_fa(clean_fa)
    fmem = collections.defaultdict(list)
    for e, c in final.items(): fmem[c].append(e)
    n_out = merged = 0
    with open(out_fa, 'w') as fh:
        for c, v in fmem.items():
            rs = [e for e in v if e in reps]
            merged += len(rs) > 1
            if rs: best = max(rs, key=lambda e: len(members[reps[e]]))
            else: best = max(v, key=lambda e: len(seqs.get(erecs[e][1], '')))
            for r in erecs[best]:
                if r in seqs and len(seqs[r]) >= 50: fh.write(f'>{names[r]}\n{seqs[r]}\n'); n_out += 1
    nmulti = sum(1 for v in fmem.values() if len(v) > 1)
    print(f'{len(fmem)} clusters ({nmulti} multi-element); {sum(1 for v in fmem.values() if any(e in reps for e in v))} hold a '
          f'consensus representative ({merged} hold more than one); {n_out} exemplar records')


if __name__ == '__main__':
    main()
