#!/usr/bin/env python3
"""Collect evidence for selected lotus elements into cases.json for the HTML page."""
import csv, json, os, re, subprocess, collections

W = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))   # LTRHarvest (holds runs/ and data/)
R = W + '/runs/lotus-ChineseLotus1.1'
G = W + '/data/genomes/lotus-ChineseLotus1.1/lotus-ChineseLotus1.1.fa'
JIANG = W + '/data/libraries/Jiang_lotus3.3_LTR'
HMM = '/usr/local/LTR_retriever-2.9.0/database/TEfam.hmm'
OUT = os.path.dirname(os.path.abspath(__file__))
PICK = ['NW_010729096.1:5829780..5835273', 'NW_010729082.1:5389356..5394739', 'NW_010729149.1:2264745..2269862',
        'NW_010729133.1:1402342..1407359', 'NW_010729194.1:69174..73461', 'NW_010729087.1:7825300..7834103',
        'NW_010729202.1:999007..1002955', 'NW_010729188.1:868007..871777']
FLANK = 30

fai = {l.split()[0]: list(map(int, l.split()[1:5])) for l in open(G + '.fai')}
gf = open(G, 'rb')
def get(c, s, e):
    L, off, lb, lB = fai[c]; s = max(1, s); e = min(L, e); s -= 1
    a = off + s // lb * lB + s % lb; b = off + (e - 1) // lb * lB + (e - 1) % lb
    gf.seek(a); return gf.read(b - a + 1).replace(b'\n', b'').decode().upper()

def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout

feat = {f"{x['seq']}:{x['elem_start']}..{x['elem_end']}": x for x in csv.DictReader(open(R + '/recipe/cand.tsv'), delimiter='\t')}

# HMM names/descriptions; press a private copy
hmm_desc = {}; name = None
for l in open(HMM):
    if l.startswith('NAME'): name = l.split()[1]
    elif l.startswith('ACC'): hmm_desc[name] = {'acc': l.split()[1]}
    elif l.startswith('DESC'): hmm_desc[name]['desc'] = l[6:].strip()
acc2name = {v['acc'].split('.')[0]: k for k, v in hmm_desc.items()}
hp = OUT + '/TEfam.hmm'
if not os.path.exists(hp + '.h3m'):
    sh(f'cp {HMM} {hp} && hmmpress -f {hp}')

CODON = {}
bases = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
i = 0
for a in bases:
    for b in bases:
        for c in bases: CODON[a + b + c] = aas[i]; i += 1
def rc(s): return s[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))
def orfs(seq):
    """6-frame ORFs (stop to stop) >= 60 aa; yields (id, aa, strand, frame, nt start offset)."""
    for strand, s in (('+', seq), ('-', rc(seq))):
        for fr in range(3):
            aa = ''.join(CODON.get(s[j:j + 3], 'X') for j in range(fr, len(s) - 2, 3))
            pos = 0
            for piece in aa.split('*'):
                if len(piece) >= 60: yield (f'{strand}{fr}_{pos}', piece, strand, fr, pos)
                pos += len(piece) + 1

# genome copies of LTR records (RepeatMasker run of the full recipe library)
cop = collections.Counter(); L = {}
for l in open(R + '/compare/recipe.fa'):
    if l.startswith('>'): n = l[1:].split()[0]; L[n] = 0
    else: L[n] += len(l.strip())
cov = collections.Counter(); hits_per = collections.Counter()
for l in open(R + '/compare/mask_recipe/lotus-ChineseLotus1.1.fa.out'):
    f = l.split()
    if len(f) < 14 or not f[0].isdigit(): continue
    b, e = (int(f[11]), int(f[12])) if f[8] == '+' else (int(f[13]), int(f[12]))
    cov[(f[9], f[-1] if f[-1] != '*' else f[-2])] += abs(e - b) + 1
    hits_per[f[9]] += 1
for (r, _), c in cov.items():
    if r in L and c >= 0.5 * L[r]: cop[r] += 1

# our own LTR records (for nearest other cluster)
own = OUT + '/own_ltr.fa'
if not os.path.exists(own + '.ndb'):
    sh(f"awk '/^>/{{p=/_LTR/}} p' {R}/compare/recipe.fa > {own} && makeblastdb -in {own} -dbtype nucl -out {own} > /dev/null")
jdb = OUT + '/jiang'
if not os.path.exists(jdb + '.ndb'):
    sh(f'makeblastdb -in {JIANG} -dbtype nucl -out {jdb} > /dev/null')
    sh(f"awk '/^>/{{p=/_INT/}} p' {JIANG} > {OUT}/jiang_int.fa && makeblastdb -in {OUT}/jiang_int.fa -dbtype nucl -out {OUT}/jiang_int > /dev/null")

def best_family(tsv, qlen):
    """Per curated family: summed bitscore, weighted identity, union query coverage, HSPs."""
    fam = collections.defaultdict(list)
    for l in tsv.splitlines():
        q, s, pid, ln, qs, qe, ss, se, ev, bits = l.split('\t')
        fam[s.split('#')[0]].append((int(qs), int(qe), int(ss), int(se), float(pid), int(ln), float(bits)))
    res = []
    for s, h in fam.items():
        iv = sorted((min(a, b), max(a, b)) for a, b, *_ in h); u = 0; end = 0
        for a, b in iv:
            if b > end: u += b - max(a, end + 1) + 1; end = b
        wid = sum(p * n for *_, p, n, _ in h) / sum(n for *_, p, n, _ in h)
        res.append({'target': s, 'bits': round(sum(x[6] for x in h), 1), 'ident': round(wid, 1), 'qcov': round(u / qlen, 3),
                    'hsps': [{'qs': min(a, b), 'qe': max(a, b), 'pid': p} for a, b, _, _, p, _, _ in h]})
    res.sort(key=lambda r: -r['bits'])
    return res[:3]

cases = []
for k in PICK:
    x = feat[k]; c = x['seq']
    es, ee = int(x['elem_start']), int(x['elem_end'])
    ls_, le = int(x['lltr_start']), int(x['lltr_end']); rs, re_ = int(x['rltr_start']), int(x['rltr_end'])
    elem = get(c, es, ee); lltr = get(c, ls_, le); rltr = get(c, rs, re_); intr = get(c, le + 1, rs - 1)
    tag = k.replace(':', '_').replace('..', '-')
    d = OUT + '/' + tag; os.makedirs(d, exist_ok=True)
    open(d + '/lltr.fa', 'w').write(f'>lltr\n{lltr}\n'); open(d + '/rltr.fa', 'w').write(f'>rltr\n{rltr}\n')
    open(d + '/int.fa', 'w').write(f'>int\n{intr}\n')
    # domains
    with open(d + '/orfs.fa', 'w') as o:
        meta = {}
        for oid, aa, strand, fr, pos in orfs(elem):
            o.write(f'>{oid}\n{aa}\n'); meta[oid] = (strand, fr, pos)
    sh(f'hmmscan --cpu 4 -E 1e-3 --domE 1e-3 --domtblout {d}/dom.tbl {hp} {d}/orfs.fa > /dev/null')
    doms = []
    for l in open(d + '/dom.tbl'):
        if l.startswith('#'): continue
        f = l.split()
        hname, oid, iev, a1, a2 = f[0], f[3], float(f[12]), int(f[17]), int(f[18])
        strand, fr, pos = meta[oid]
        n1, n2 = fr + 3 * (pos + a1 - 1), fr + 3 * (pos + a2)  # 0-based nt on the strand
        if strand == '-': n1, n2 = len(elem) - n2, len(elem) - n1
        doms.append({'name': hname, 'acc': hmm_desc[hname]['acc'], 'desc': hmm_desc[hname].get('desc', ''), 'ievalue': iev,
                     'start': n1 + 1, 'end': n2, 'strand': strand})
    doms.sort(key=lambda r: r['ievalue']); kept = []
    for dm in doms:  # drop weaker overlapping hits
        if all(dm['end'] < k2['start'] or dm['start'] > k2['end'] for k2 in kept): kept.append(dm)
    kept.sort(key=lambda r: r['start'])
    fmt = '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore'
    near_ltr = best_family(sh(f'blastn -task blastn -query {d}/lltr.fa -db {jdb} -evalue 1e-5 -max_target_seqs 50 -outfmt "{fmt}"'), len(lltr))
    near_int = best_family(sh(f'blastn -task blastn -query {d}/int.fa -db {jdb} -evalue 1e-5 -max_target_seqs 50 -outfmt "{fmt}"'), len(intr))
    near_int_aa = best_family(sh(f'tblastx -query {d}/int.fa -db {OUT}/jiang_int -evalue 1e-5 -max_target_seqs 20 -outfmt "{fmt}"'), len(intr))
    own_hits = sh(f'blastn -task blastn -query {d}/lltr.fa -db {own} -evalue 1e-5 -max_target_seqs 50 -outfmt "{fmt}"')
    own_hits = '\n'.join(l for l in own_hits.splitlines() if not l.split('\t')[1].startswith(c + ':') or
                         not (es <= int(re.search(r':(\d+)\.\.', l.split('\t')[1]).group(1)) <= ee))
    near_own = best_family(own_hits, len(lltr)) if own_hits else []
    ltr_pair = sh(f'blastn -task blastn -query {d}/lltr.fa -subject {d}/rltr.fa -outfmt 0 -max_hsps 1 -dust no')
    m = re.search(r'(Score = .*?)(?=\n\nLambda|\Z)', ltr_pair, re.S)
    recs = [n for n in L if n.startswith(c + ':') and es <= int(re.search(r':(\d+)\.\.', n).group(1)) <= ee]
    lrec = [n for n in recs if n.endswith('_LTR')]
    cases.append({
        'id': k, 'seq': c, 'start': es, 'end': ee, 'len': ee - es + 1,
        'lltr': [ls_, le], 'rltr': [rs, re_], 'ltr_identity': float(x['identity']),
        'tsd_len': int(x['tsd_len']), 'tsd_seq': x['tsd_seq'], 'motif': x['motif'], 'motif_class': x['motif_class'],
        'pbs_offset': int(x['pbs_offset']), 'ppt_offset': int(x['ppt_offset']), 'ppt_purines15': int(x['ppt_purines15']),
        'n_count': int(x['n_count']), 'int_n': intr.count('N'),
        'superfamily': x['dom_superfamily'], 'dom_order': x['dom_order'], 'dom_strand': x['dom_strand'], 'evidence': x['evidence'],
        'flank5': get(c, es - FLANK, es - 1), 'ltr5_head': lltr[:30], 'ltr5_tail': lltr[-30:],
        'ltr3_head': rltr[:30], 'ltr3_tail': rltr[-30:], 'flank3': get(c, ee + 1, ee + FLANK),
        'after5ltr': get(c, le + 1, le + 40), 'before3ltr': get(c, rs - 40, rs - 1),
        'domains': kept, 'near_ltr': near_ltr, 'near_int': near_int, 'near_int_aa': near_int_aa, 'near_own': near_own[:2],
        'ltr_pair_aln': m.group(1).strip() if m else '',
        'ltr_copies': {n: {'copies50': cop[n], 'hits': hits_per[n], 'len': L[n]} for n in lrec},
    })
    print(k, len(kept), 'domains;', 'LTR near', near_ltr[0]['target'] if near_ltr else '-', 'INT aa', near_int_aa[0]['ident'] if near_int_aa else '-')
json.dump(cases, open(OUT + '/cases.json', 'w'), indent=1)
