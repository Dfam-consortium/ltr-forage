"""Emit GT-format prediction lines from an LTR_retriever run: the .scn.adj rows for
loci whose defalse decision is pass or truncated (the ones that reach the library)."""
import sys, re
adj, defalse, want = sys.argv[1], sys.argv[2], set(sys.argv[3].split(','))
keep = {}
for blk in open(defalse).read().split('\n\n'):
    m = re.match(r'\S+:(\d+)\.\.(\d+)\t(\w+)\t', blk)
    if not m: continue
    s, e, dec = int(m.group(1)), int(m.group(2)), m.group(3)
    if dec in want: keep[(min(s,e), max(s,e))] = dec
n = 0
for line in open(adj):
    c = line.split()
    if len(c) < 11 or not c[0].isdigit(): continue
    s, e = int(c[0]), int(c[1])
    if (min(s,e), max(s,e)) in keep:
        print(' '.join(c[:11])); n += 1
print(f"# {n} of {len(keep)} kept loci written", file=sys.stderr)
