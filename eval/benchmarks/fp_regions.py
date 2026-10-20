#!/usr/bin/env python3
"""Where a library's false-positive masking comes from, region by region.

Usage: fp_regions.py <genome.fai> <standard.out> <full_annotation.out> <test.out> <library.fa> <clusters.json> <allvsall.tsv> [--top 30]

A masked base is a false positive when it lies outside the standard's LTR
mask.  Each RepeatMasker hit of the test library is projected onto its
library record, so every record gets a profile of false-positive hit depth
in 50 bp bins.  Runs of bins with at least five false-positive hits where
false positives are at least half the hits are the record's false-positive
regions.  For each region the script reports:

  - what the full annotation says the false-positive loci are (majority
    family and class), which is the evidence that they are not this LTR
    family's sequence;
  - how many genomic hits span the region and 100 bp of record on both
    sides (copies of the mosaic: the element was copied after the
    insertion) against hits confined to the region (copies of the insert
    elsewhere);
  - how many other library records carry the region's sequence, and the
    size of the record's element cluster (whether a partner existed to
    excise it by alignment).

The clusters file maps element names to cluster ids as step4.py writes it.
"""
import sys, json, collections, bisect, argparse
import numpy as np

ap = argparse.ArgumentParser()
for a in ('fai', 'standard', 'full', 'test', 'lib', 'clusters', 'ava'): ap.add_argument(a)
ap.add_argument('--top', type=int, default=30)
ap.add_argument('--bin', type=int, default=50)
ap.add_argument('--min-fp-hits', type=int, default=5)
ap.add_argument('--flank', type=int, default=100)
a = ap.parse_args()
BIN = a.bin

chrom_len = {}
for l in open(a.fai):
    c = l.split('\t'); chrom_len[c[0]] = int(c[1])
std = {c: np.zeros(n + 1, dtype=bool) for c, n in chrom_len.items()}
fam_id = {}; fam_name = []; fam_class = []
def fid(name, cls):
    if name not in fam_id:
        fam_id[name] = len(fam_name); fam_name.append(name); fam_class.append(cls)
    return fam_id[name]
fid('unannotated', 'none')
full = {c: np.zeros(n + 1, dtype=np.int32) for c, n in chrom_len.items()}
def rows(path):
    for line in open(path):
        c = line.split()
        if len(c) < 15 or not c[0].isdigit(): continue
        yield c
for c in rows(a.standard):
    if c[10].startswith('LTR') and c[4] in std:
        std[c[4]][int(c[5]):int(c[6]) + 1] = True
# the full annotation: later rows overwrite earlier ones where they overlap, which is fine at this resolution
for c in rows(a.full):
    if c[4] in full:
        full[c[4]][int(c[5]):int(c[6]) + 1] = fid(c[9], c[10])

seqlen = {}; name = None
for l in open(a.lib):
    if l.startswith('>'): name = l[1:].split()[0].split('#')[0]; seqlen[name] = 0
    else: seqlen[name] += len(l.strip())

def hit_fields(c):
    chrom, gs, ge = c[4], int(c[5]), int(c[6])
    rec = c[9].split('#')[0]
    if c[8] == '+': rs, re_ = int(c[11]), int(c[12])
    else: rs, re_ = int(c[13]), int(c[12])
    return chrom, gs, ge, c[8], rec, min(rs, re_), max(rs, re_)

# pass 1: per-record profiles
fp_depth = {}; all_depth = {}
total_masked = 0; total_fp = 0
for c in rows(a.test):
    chrom, gs, ge, strand, rec, rs, re_ = hit_fields(c)
    if chrom not in std or rec not in seqlen: continue
    L = seqlen[rec]; nb = (L - 1) // BIN + 1
    if rec not in fp_depth: fp_depth[rec] = np.zeros(nb, dtype=np.int32); all_depth[rec] = np.zeros(nb, dtype=np.int32)
    fp = ~std[chrom][gs:ge + 1]
    total_masked += ge - gs + 1; nfp = int(fp.sum()); total_fp += nfp
    b0, b1 = (rs - 1) // BIN, (re_ - 1) // BIN
    all_depth[rec][b0:b1 + 1] += 1
    if nfp == 0: continue
    # genome positions of false-positive bases, projected onto the record
    gpos = np.nonzero(fp)[0]
    frac = gpos / max(1, ge - gs)
    rpos = rs + frac * (re_ - rs) if strand == '+' else re_ - frac * (re_ - rs)
    bins = np.unique(((rpos - 1) // BIN).astype(int))
    bins = bins[(bins >= 0) & (bins < nb)]
    fp_depth[rec][bins] += 1

# regions
regions = {}
for rec, fpd in fp_depth.items():
    ad = all_depth[rec]
    hot = (fpd >= a.min_fp_hits) & (fpd * 2 >= ad)
    runs = []; i = 0; nb = len(hot)
    while i < nb:
        if hot[i]:
            j = i
            while j < nb and hot[j]: j += 1
            if (j - i) * BIN >= 2 * BIN: runs.append((i * BIN + 1, min(seqlen[rec], j * BIN)))
            i = j
        else: i += 1
    if runs: regions[rec] = runs

# pass 2: per-region evidence
reg_fp = collections.Counter(); reg_fam = collections.defaultdict(collections.Counter)
reg_span = collections.Counter(); reg_inside = collections.Counter(); reg_hits = collections.Counter()
for c in rows(a.test):
    chrom, gs, ge, strand, rec, rs, re_ = hit_fields(c)
    if rec not in regions or chrom not in std: continue
    fp = ~std[chrom][gs:ge + 1]
    ann = full[chrom][gs:ge + 1]
    for (ra, rb) in regions[rec]:
        if re_ < ra or rs > rb: continue
        key = (rec, ra, rb)
        reg_hits[key] += 1
        if rs <= ra - a.flank and re_ >= rb + a.flank: reg_span[key] += 1
        elif rs >= ra - 50 and re_ <= rb + 50: reg_inside[key] += 1
        # false-positive bases of this hit that fall inside the region (projected)
        if strand == '+':
            lo = gs + (ra - rs) * (ge - gs) / max(1, re_ - rs); hi = gs + (rb - rs) * (ge - gs) / max(1, re_ - rs)
        else:
            lo = gs + (re_ - rb) * (ge - gs) / max(1, re_ - rs); hi = gs + (re_ - ra) * (ge - gs) / max(1, re_ - rs)
        lo = int(max(gs, lo)); hi = int(min(ge, hi))
        if hi < lo: continue
        seg = slice(lo - gs, hi - gs + 1)
        f = fp[seg]; nfp = int(f.sum())
        if nfp == 0: continue
        reg_fp[key] += nfp
        ids, counts = np.unique(ann[seg][f], return_counts=True)
        for i, n in zip(ids, counts): reg_fam[key][int(i)] += int(n)

# sharing and cluster size
clus = json.load(open(a.clusters)); size = collections.Counter(clus.values())
by_seq = collections.defaultdict(list)
for e in clus:
    s, rest = e.split(':', 1); x, y = rest.split('..'); by_seq[s].append((int(x), int(y), e))
for s in by_seq: by_seq[s].sort()
def cluster_size(rec):
    s, rest = rec.split(':', 1); x, y = rest.split('_')[0].split('..'); x, y = int(x), int(y); lst = by_seq.get(s, [])
    i = bisect.bisect_right(lst, (x, float('inf'), '')) - 1
    while i >= 0 and lst[i][1] >= x:
        if lst[i][0] <= x and y <= lst[i][1]: return size[clus[lst[i][2]]]
        i -= 1
    return 0
shared = collections.defaultdict(set)
region_list = {rec: rs for rec, rs in regions.items()}
for l in open(a.ava):
    q, s, pid, ln, qs, qe, ss, se, ql, sl, bits = l.split('\t'); q = q.split('#')[0]; s = s.split('#')[0]
    if q == s or q not in region_list or float(pid) < 80: continue
    x, y = min(int(qs), int(qe)), max(int(qs), int(qe))
    for (ra, rb) in region_list[q]:
        if min(y, rb) - max(x, ra) + 1 >= 0.8 * (rb - ra + 1): shared[(q, ra, rb)].add(s)

# report
keys = sorted(reg_fp, key=lambda k: -reg_fp[k])
in_regions = sum(reg_fp.values())
print(f"masked {total_masked:,} bp, false-positive {total_fp:,} bp; {len(keys)} false-positive regions in {len(regions)} records hold {in_regions:,} bp ({100 * in_regions / max(1, total_fp):.1f}% of the false positives)")
print()
print("| record | region | length | FP bases | full annotation of the FP loci (majority; share) | hits spanning region + flanks | hits inside region only | other records carrying it | cluster size |")
print("|---|---|---|---|---|---|---|---|---|")
def majority(k):
    fam = reg_fam[k]
    if not fam: return ('none', 'none', 0.0)
    tot = sum(fam.values()); i, n = fam.most_common(1)[0]
    return (fam_name[i], fam_class[i], n / tot)
for k in keys[:a.top]:
    rec, ra, rb = k; fn, fc, sh = majority(k)
    print(f"| {rec} | {ra}-{rb} | {rb - ra + 1} | {reg_fp[k]:,} | {fn} ({fc}; {100 * sh:.0f}%) | {reg_span[k]} | {reg_inside[k]} | {len(shared[k])} | {cluster_size(rec)} |")
print()
# summaries over all regions, weighted by FP bases
def coarse(c):
    c = c.split('/')[0]
    return c
by_class = collections.Counter(); by_copy = collections.Counter(); by_share = collections.Counter(); by_pattern = collections.Counter()
for k in keys:
    rec, ra, rb = k; fn, fc, sh = majority(k); n = reg_fp[k]
    by_class[coarse(fc)] += n
    mosaic = reg_span[k] >= 3
    by_copy['mosaic copied (3+ hits span region and flanks)' if mosaic else ('few or no copies span it' )] += n
    sh_n = len(shared[k]); cs = cluster_size(rec)
    if sh_n >= 1: by_share['sequence present in other library records'] += n
    elif cs >= 2: by_share['unique to this record, cluster has partners'] += n
    else: by_share['unique to this record, singleton cluster (no partner to excise by alignment)'] += n
tot = max(1, in_regions)
print("False-positive bases in regions by the full annotation's class of the loci:")
for c, n in by_class.most_common(): print(f"  {100 * n / tot:5.1f}%  {c}")
print("By evidence that the element was copied after the insertion:")
for c, n in by_copy.most_common(): print(f"  {100 * n / tot:5.1f}%  {c}")
print("By whether alignment could have found it:")
for c, n in by_share.most_common(): print(f"  {100 * n / tot:5.1f}%  {c}")
