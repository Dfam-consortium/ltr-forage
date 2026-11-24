#!/usr/bin/env python3
"""Turn each library record to its element's coding strand (ltr_library.py step 23, after family.py).

usage: orient_models.py <library.fa> <families.tsv> <cand.tsv>

Every model is in the genome strand of the element it was founded on, its origin (a consensus model is
turned to that element's record when it is built).  The element's coding strand comes from ltr-forage's
table: its protein domains' strand, else its PBS's (the PBS follows the 5' LTR); none when the two disagree
or neither was found.  A record whose element is on the minus strand is reverse-complemented, in place.
families.tsv gets two columns: orientation, what the strand rests on (domains+PBS, domains, PBS, conflict,
unknown), and flipped (yes, no).  Records of unknown orientation stay in their element's genome strand.
Names and origin= tags do not change.
"""
import csv, os, sys
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from records import element_lookup

COMP = str.maketrans('ACGTRYKMSWBDHVNacgtrykmswbdhvn', 'TGCAYRMKSWVHDBNtgcayrmkswvhdbn')


def rc(s): return s[::-1].translate(COMP)


def strand_of(r):
    """(strand, what it rests on) for one ltr-forage table row."""
    d = r.get('dom_strand') if r.get('dom_strand') in ('+', '-') else None
    p = r.get('pbs_strand') if r.get('pbs_strand') in ('+', '-') else None
    if d and p and d != p: return None, 'conflict'
    if d: return d, 'domains+PBS' if p else 'domains'
    if p: return p, 'PBS'
    return None, 'unknown'


def main():
    lib_path, fam_path, cand_path = sys.argv[1:4]
    rows = {f"{r['seq']}:{r['elem_start']}..{r['elem_end']}": r
            for r in csv.DictReader(open(cand_path), delimiter='\t') if r['keep'] == '1'}
    element = element_lookup(rows, tsv=cand_path)
    fams = list(csv.DictReader(open(fam_path), delimiter='\t'))
    fields = list(fams[0].keys()) if fams else []
    for k in ('orientation', 'flipped'):
        if k not in fields: fields.append(k)
    turn = {}
    for f in fams:
        e = element(f['origin'])
        strand, why = strand_of(rows[e]) if e in rows else (None, 'unknown')
        f['orientation'], f['flipped'] = why, 'yes' if strand == '-' else 'no'
        turn[f['new_name']] = strand == '-'

    # library.fa: header lines kept, sequences of flipped records reverse-complemented
    recs, name, seq = [], None, []
    for l in open(lib_path):
        if l.startswith('>'):
            if name is not None: recs.append((name, ''.join(seq)))
            name, seq = l.rstrip('\n'), []
        else:
            seq.append(l.strip())
    if name is not None: recs.append((name, ''.join(seq)))
    with open(lib_path + '.tmp', 'w') as fo:
        for h, s in recs:
            key = h[1:].split()[0].split('#')[0]
            fo.write(f'{h}\n{rc(s) if turn.get(key) else s}\n')
    with open(fam_path + '.tmp', 'w') as fo:
        w = csv.DictWriter(fo, fieldnames=fields, delimiter='\t', lineterminator='\n'); w.writeheader()
        for f in fams: w.writerow(f)
    os.replace(lib_path + '.tmp', lib_path); os.replace(fam_path + '.tmp', fam_path)

    n = {k: sum(1 for f in fams if f['orientation'] == k) for k in ('domains+PBS', 'domains', 'PBS', 'conflict', 'unknown')}
    print(f"orientation: {sum(turn.values())} of {len(fams)} records reverse-complemented to their element's coding strand; "
          f"strand from domains and PBS {n['domains+PBS']}, domains {n['domains']}, PBS {n['PBS']}; "
          f"left as found: {n['unknown']} without either, {n['conflict']} where the two disagree")


if __name__ == '__main__':
    main()
