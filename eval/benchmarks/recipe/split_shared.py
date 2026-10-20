#!/usr/bin/env python3
"""Split an ltr-forage library into records whose element LTR_retriever also kept
(reciprocal overlap >= 0.8 with a pass-list element, or with an element named in
LTRlib.fa) and records only ltr-forage kept."""
import sys, re, collections
lib, passlist, lrlib, out_shared, out_only = sys.argv[1:6]
lr = collections.defaultdict(list)
for l in open(passlist):
    m = re.match(r'(\S+):(\d+)\.\.(\d+)', l)
    if m: lr[m[1]].append((int(m[2]), int(m[3])))
for l in open(lrlib):
    m = re.match(r'>(\S+):(\d+)\.\.(\d+)_', l)
    if m: lr[m[1]].append((int(m[2]), int(m[3])))
def shared(c, s, e):
    for a, b in lr.get(c, ()):
        ov = min(e, b) - max(s, a) + 1
        if ov > 0 and ov >= 0.8*(e-s+1) and ov >= 0.8*(b-a+1): return True
    return False
# element span per record: group LTR/INT records by element via TSV-free approach:
# records are named chr:s..e_LTR / _INT; the element is the union of an element's records,
# which appear consecutively (LTR, INT, LTR).
recs = []; name = None; seq = []
for l in open(lib):
    if l.startswith('>'):
        if name: recs.append((name, ''.join(seq)))
        name = l[1:].rstrip(); seq = []
    else: seq.append(l.rstrip())
if name: recs.append((name, ''.join(seq)))
# group records into elements: LTR, INT, LTR, each starting where the previous ended
elems = []
for n, s in recs:
    m = re.match(r'(\S+):(\d+)\.\.(\d+)_', n)
    if elems and '_INT' in n or (elems and '_LTR' in n and '_INT' in elems[-1][-1][0]):
        pm = re.match(r'(\S+):(\d+)\.\.(\d+)_', elems[-1][-1][0])
        if pm[1] == m[1] and int(m[2]) == int(pm[3]) + 1:
            elems[-1].append((n, s)); continue
    elems.append([(n, s)])
fs, fo = open(out_shared, 'w'), open(out_only, 'w')
ns = no = 0
for el in elems:
    coords = [re.match(r'(\S+):(\d+)\.\.(\d+)_', n) for n, _ in el]
    c = coords[0][1]; s = min(int(m[2]) for m in coords); e = max(int(m[3]) for m in coords)
    f = fs if shared(c, s, e) else fo
    if f is fs: ns += 1
    else: no += 1
    for n, sq in el: f.write(f'>{n}\n{sq}\n')
print(f'elements: shared {ns}, ltr-forage only {no}', file=sys.stderr)
