#!/usr/bin/env python3
"""Drop raw singletons that step 21's tier would reject, before the depth cleanup (ltr_library.py --early-tier).
usage: early_tier.py <outdir> <tier: pbs|strict|all> [--restore]

A raw singleton (a pre-clustering cluster of one element) is dropped when its element fails the tier:
longer than 9 kb, or (pbs) LTR identity < 93% and no PBS, or (strict) identity < 90% or no PBS.
elem/clusters_raw.json, lib.fa and exemplar_raw.fa are rewritten from copies of the originals kept as
<file>.all, so a rerun starts from the full raw library; --restore puts the originals back."""
import sys, os, json, csv, shutil, collections
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from records import element_lookup
d, tier = sys.argv[1], sys.argv[2]
FILES = ['elem/clusters_raw.json', 'lib.fa', 'exemplar_raw.fa']
for f in FILES:
    p = os.path.join(d, f)
    if not os.path.exists(p + '.all'): shutil.copy(p, p + '.all')
    else: shutil.copy(p + '.all', p)
if '--restore' in sys.argv or tier == 'all':
    print('early tier: off'); sys.exit()

raw = json.load(open(os.path.join(d, 'elem/clusters_raw.json')))
size = collections.Counter(raw.values())
feat = {f"{x['seq']}:{x['elem_start']}..{x['elem_end']}": x for x in csv.DictReader(open(os.path.join(d, 'cand.tsv')), delimiter='\t')}
def admit(k):   # tier.py's rule
    x = feat[k]
    if int(x['elem_len']) > 9000: return False
    pbs = x['pbs_mm12'] != '.'
    return (float(x['identity']) >= 90 and pbs) if tier == 'strict' else (float(x['identity']) >= 93 or pbs)
drop = {k for k, c in raw.items() if size[c] == 1 and not admit(k)}

element = element_lookup(raw, os.path.join(d, 'cand.tsv'))   # the element a record (LTR or internal region) is of

def filter_fasta(p):
    n = kept = 0; keep = True
    with open(p + '.all') as fi, open(p, 'w') as fo:
        for l in fi:
            if l.startswith('>'):
                n += 1; keep = element(l[1:].split()[0]) not in drop; kept += keep
            if keep: fo.write(l)
    return n, kept
json.dump({k: c for k, c in raw.items() if k not in drop}, open(os.path.join(d, 'elem/clusters_raw.json'), 'w'))
n, kept = filter_fasta(os.path.join(d, 'lib.fa')); filter_fasta(os.path.join(d, 'exemplar_raw.fa'))
print(f'early tier ({tier}): {len(drop)} of {sum(1 for c in raw.values() if size[c] == 1)} raw singletons dropped, '
      f'{len(raw) - len(drop)} of {len(raw)} elements and {kept} of {n} records left')
