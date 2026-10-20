#!/usr/bin/env python3
"""Build an LTR retrotransposon family library from a genome assembly with ltr-forage.

usage: ltr_library.py <genome.fa> -o <outdir> [-t THREADS] [options]

Steps (logs in <outdir>/logs, timings in <outdir>/steps.log):
   0  index      sufr suffix array of the genome (skipped with --sufr)
   0  sample     depth sample: the whole genome up to --sample-above Mb, else a random
                 --sample-frac of its --window-kb windows (seeded); depth floors are scaled to it
 1-5  ltr-forage detect, re-extend, annotate (domains, TSD, PBS), decide, raw library
   5b precluster element clusters of the raw library (all-vs-all, as step 8)
   5c premodels  te-composer consensus for raw clusters of >= --min-cons elements; from here on
                 each big cluster is represented by its consensus, and steps 6-9 work on the elements of
                 smaller clusters plus those representatives (the members of a big cluster join the
                 cluster of its representative after step 9).  5b-5c are stage B, the default;
                 --stage-a (or --min-cons < 2) skips them and steps 6-9 work on every element
   6  depth      three rounds: rmblastn of the sample against the library, excise stretches far
                 above their record's own genome depth (nested passengers), drop high-depth LTRs.
                 A representative is cut only at 10x its own depth and never dropped: a big family's
                 consensus is high-copy itself, while a real insertion stands far above it
   7  nested     first element clustering (library all-vs-all), nested LTR excision by vote
 8-9  cluster    element clusters and one exemplar element per cluster
   9b flank      flag clusters whose copies share both flanks (segmental duplications; --no-flank):
                 most pairs share both flanks, or every copy shares both with some other copy
  10  contain    drop exemplars contained in another at >= 0.9
  11  excise     nested LTR excision from internals, two rounds
  12  pbs        primer binding site per element (ltr-forage's tRNA search)
  13  tier       keep multi-element clusters; singletons with LTR identity >= 93 or a PBS
                 (--tier strict: identity >= 90 and a PBS)
  13b models     te-composer consensus for clusters of >= --min-cons elements (oriented, columns held
                 by under half the median occupancy dropped; the exemplar where none can be built), the
                 exemplar for smaller ones; a stage B representative is kept as steps 6-11 left it;
                 flagged segmental-duplication clusters left out (--min-cons 0: exemplars only)
  14  families   80-80 families (LTRs and internals both), gap-aware, Gypsy/Copia class names

Outputs: library.fa (fam<N>_v<k>_LTR|INT#LTR/<class>), families.tsv, cand.tsv (every candidate
element, ltr-forage's table).  A rerun skips finished steps (<outdir>/.done); --from STEP
reruns from a step on.

Stage B against stage A (2026-10-20; sensitivity / precision of LTR masking against curated libraries):
rice 92.70 / 98.53 vs 92.51 / 98.29, zebra finch 74.76 / 83.03 vs 74.67 / 83.45 with steps 6-14 in 16 min
instead of ~109, gorilla 34.57 / 97.65 (HiRepEat 17.68 / 95.97).  Most of zebra finch's remaining FP is
outside the curated library and much of it looks like real families HiRepEat also finds: FP rules need
an artifact mechanism, not lower novelty.

External programs: ltr-forage, sufr, rmblastn (Rust) and te-composer (dfam-curator >= 0.6.6), found from
the options, the environment (LTRFORAGE, SUFR, RMBLASTN, TE_COMPOSER) or PATH; Python 3 with numpy."""
import argparse, os, random, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RECIPE = os.path.join(ROOT, 'eval', 'benchmarks', 'recipe')
BENCH = os.path.join(ROOT, 'eval', 'benchmarks')
STEPS = ['index', 'sample', 'ltr-forage', 'precluster', 'premodels', 'depth', 'nested', 'cluster', 'flank', 'contain', 'excise', 'pbs', 'tier',
         'models', 'families']
AVA = 'qseqid sseqid pident length qstart qend sstart send qlen slen bitscore'.split()


def find_tool(opt, env, name, fallbacks=(), check=None):
    """First executable from the option, the environment variable, PATH, then the fallbacks;
    with check, only one whose --help output contains that text (NCBI's rmblastn is often on
    PATH too, and only the Rust port takes the options used here)."""
    for p in (opt, os.environ.get(env), shutil.which(name), *fallbacks):
        if not (p and os.path.isfile(p) and os.access(p, os.X_OK)):
            continue
        if check and check not in subprocess.run([p, '--help'], capture_output=True, text=True).stdout:
            continue
        return os.path.abspath(p)
    sys.exit(f'ltr_library: cannot find {name}' + (f' ({check})' if check else '') + f'; give --{name} or set {env}')


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, ''.join(seq)


def genome_stats(path):
    total = ns = 0
    for _, s in read_fasta(path):
        total += len(s); ns += s.upper().count('N')
    return total, ns


def write_sample(genome, out, frac, window, seed):
    """A random frac of the genome's window-bp windows (named <seq>_<start>), seeded."""
    rng = random.Random(seed)
    kept = 0
    with open(out, 'w') as fh:
        for name, s in read_fasta(genome):
            for st in range(0, len(s), window):
                if rng.random() < frac:
                    w = s[st:st + window]
                    fh.write(f'>{name}_{st}\n')
                    for i in range(0, len(w), 80):
                        fh.write(w[i:i + 80] + '\n')
                    kept += len(w)
    return kept


class Run:
    def __init__(self, a):
        self.a = a
        self.out = os.path.abspath(a.outdir)
        os.makedirs(os.path.join(self.out, 'logs'), exist_ok=True)
        os.makedirs(os.path.join(self.out, 'elem'), exist_ok=True)
        os.makedirs(os.path.join(self.out, '.done'), exist_ok=True)
        self.env = dict(os.environ, TSV=os.path.join(self.out, 'cand.tsv'), JIANG='/dev/null',
                        RMBLASTN=a.rmblastn, WORD=str(a.word))
        if a.from_step:
            for s in STEPS[STEPS.index(a.from_step):]:
                p = self.done_path(s)
                if os.path.exists(p): os.remove(p)
        # stage A and B outputs from step 5b on don't mix: a rerun keeps the mode its finished steps were made in
        mode = 'B' if a.stage_b else 'A'
        mf = self.done_path('mode')
        made = open(mf).read().strip() if os.path.exists(mf) else ('B' if os.path.exists(self.path('reps.fa')) else 'A')
        later = [s for s in STEPS[STEPS.index('precluster'):] if os.path.exists(self.done_path(s))]
        if later and made != mode:
            sys.exit(f'ltr_library: {self.out} has steps from stage {made} ({later[0]} ...); rerun with '
                     f'{"--stage-a" if made == "A" else "no --stage-a"}, or --from precluster to redo them as stage {mode}')
        open(mf, 'w').write(mode + '\n')

    def done_path(self, step): return os.path.join(self.out, '.done', step)
    def path(self, *p): return os.path.join(self.out, *p)

    def log(self, msg):
        line = f'[{time.strftime("%H:%M:%S")}] {msg}'
        print(line, flush=True)
        with open(self.path('steps.log'), 'a') as fh: fh.write(line + '\n')

    def sh(self, cmd, log, stdout=None, threads=None):
        """Run one command in outdir; stderr (and stdout unless redirected) to logs/<log>."""
        env = dict(self.env, THREADS=str(threads or self.a.threads))
        with open(self.path('logs', log), 'a') as lf:
            lf.write('$ ' + ' '.join(map(str, cmd)) + '\n'); lf.flush()
            out = open(self.path(stdout), 'w') if stdout else lf
            try:
                r = subprocess.run([str(c) for c in cmd], cwd=self.out, env=env, stdout=out, stderr=lf)
            finally:
                if stdout: out.close()
        if r.returncode:
            sys.exit(f'ltr_library: {cmd[0]} failed (exit {r.returncode}); see {self.path("logs", log)}')

    def py(self, script, *args, log, stdout=None, where=RECIPE):
        self.sh([sys.executable, os.path.join(where, script), *args], log, stdout)

    def search(self, q, db, max_hsps, out, fields, log):
        self.sh(['bash', os.path.join(RECIPE, 'rmb_ava.sh'), q, db, max_hsps, out, *fields], log)

    def step(self, name, title, fn):
        if os.path.exists(self.done_path(name)):
            return
        self.log(title)
        fn()
        open(self.done_path(name), 'w').close()

    # ---- steps -----------------------------------------------------------
    def index(self):
        self.sh([self.a.sufr, '-t', self.a.threads, 'create', '--dna', '-o', self.path('genome.sufr'), self.genome], 'index.log')

    def sample(self):
        a = self.a
        total, ns = genome_stats(self.genome)
        mb = (total - ns) / 1e6
        if a.sample_file:
            frac = genome_stats(a.sample_file)[0] / total; src = os.path.abspath(a.sample_file)
            shutil.copy(src, self.path('depth_sample.fa'))
        elif mb <= a.sample_above:
            frac = 1.0; os.symlink(self.genome, self.path('depth_sample.fa'))
        else:
            kept = write_sample(self.genome, self.path('depth_sample.fa'), a.sample_frac, a.window_kb * 1000, a.seed)
            frac = kept / total
        # depth floors in whole-genome units, scaled with the non-N genome size (rice MSU7,
        # 373 Mb: 50 / 500) and capped at 50: every tuned genome used <= 50, and lotus (708 Mb,
        # uncapped 95 / 950) lost 2 points of precision (90.35 vs 92.36 at 50 / 500) because higher
        # floors trim less; excise_by_depth.py multiplies them by the sample's --scale
        mind = a.min_depth or min(50, max(5, round(50 * mb / 373)))
        absd = a.abs_depth or 10 * mind
        with open(self.path('sample.txt'), 'w') as fh:
            fh.write(f'genome_mb\t{mb:.1f}\nn_frac\t{ns / total:.4f}\nscale\t{frac:.4f}\nmin_depth\t{mind}\nabs_depth\t{absd}\n')
        self.log(f'  genome {mb:.1f} Mb non-N; depth sample {frac:.1%} of the genome; floors {mind} / {absd}')

    def params(self):
        return dict(l.rstrip('\n').split('\t') for l in open(self.path('sample.txt')))

    def ltr_forage(self):
        self.sh([self.a.ltr_forage, self.sufr, '--reextend', '--reextend-ltr-max', 7000, '--format', 'tsv',
                 '--domain-hmm', self.a.domain_hmm, '--library', 'lib.raw.fa', '-j', self.a.threads],
                'ltr-forage.log', stdout='cand.tsv')
        with open(self.path('lib.raw.fa')) as fi, open(self.path('lib.fa'), 'w') as fo:
            for l in fi: fo.write(l.split()[0] + '\n' if l.startswith('>') else l)

    def precluster(self):
        if not self.a.stage_b: return
        self.search('lib.fa', 'lib.fa', 5, 'raw_ava.tsv', AVA, 'precluster.log')
        self.py('step4.py', 'raw_ava.tsv', 'raw_ava.tsv', 'elem/clusters_raw.json', log='precluster.log', stdout='elem/raw.txt')
        self.py('exemplar_lib.py', 'elem/clusters_raw.json', 'lib.fa', 'exemplar_raw.fa', log='precluster.log')
        self.log('  ' + open(self.path('elem', 'raw.txt')).read().strip().split('\n')[-1])

    def premodels(self):
        if not self.a.stage_b: return
        a = self.a
        self.sh([sys.executable, os.path.join(HERE, 'consensus_models.py'), 'exemplar_raw.fa', 'lib.fa', 'elem/clusters_raw.json',
                 'reps.fa', '--big-only', '--min-cons', a.min_cons, '--jobs', max(1, a.threads // 2),
                 '--te-composer', a.te_composer, '--work', 'cons'], 'premodels.log', stdout='premodels.txt')
        self.sh([sys.executable, os.path.join(HERE, 'stage_b.py'), 'split', 'elem/clusters_raw.json', 'lib.fa', 'reps.fa',
                 'depth_lib.fa', 'elem/reps.json', '--min-cons', a.min_cons], 'premodels.log', stdout='split.txt')
        self.log('  ' + open(self.path('split.txt')).read().strip())

    def depth(self):
        p = self.params()
        mtx = subprocess.run(['bash', os.path.join(RECIPE, 'pick_matrix.sh'), 'depth_sample.fa'], cwd=self.out,
                             capture_output=True, text=True, check=True).stdout.strip()
        batches = 16
        self.env['BATCHES'] = str(batches)
        self.sh(['bash', os.path.join(BENCH, 'clean_library_rmb.sh'), 'depth_lib.fa' if self.a.stage_b else 'lib.fa', 'depth_sample.fa', 'clean', 3, p['scale'], mtx,
                 '--ratio', 3, '--min-depth', p['min_depth'], '--min-len', 50, '--abs-depth', p['abs_depth'], '--drop-frac', 0.7,
                 *(['--consensus', 'reps.fa'] if self.a.stage_b else [])],
                'depth.log', threads=max(1, self.a.threads // batches))
        if self.a.stage_b:
            # representatives stay in the search (so hits compete for the sample as in stage A) and take its
            # excisions, but one that depth cleanup dropped outright comes back as its consensus: a big family's
            # only LTR must survive (zebra finch TguERVK6: high-copy LTR dropped by --abs-depth, family lost its
            # LTR model; 5 of 110 representative LTRs dropped).  Restoring every representative instead put back
            # insertions depth had rightly cut (gorilla: a 5' Alu in one INT consensus, 294 Mb masked, precision 23%)
            reps = dict(read_fasta(self.path('reps.fa')))
            kept = list(read_fasta(self.path('clean', 'clean.fa')))
            have = {n.split('#')[0] for n, _ in kept}
            lost = [(n, s) for n, s in reps.items() if n.split('#')[0] not in have]
            with open(self.path('clean', 'clean.fa'), 'w') as fh:
                for n, s in kept + lost: fh.write(f'>{n}\n{s}\n')
            self.log(f'  {len(lost)} of {len(reps)} representatives dropped by depth cleanup restored as consensus')

    def nested(self):
        self.search('clean/clean.fa', 'clean/clean.fa', 5, 'clean_ava.tsv', AVA, 'nested.log')
        self.py('step4.py', 'clean_ava.tsv', 'clean_ava.tsv', 'elem/clusters_clean.json', log='nested.log', stdout='elem/clean.txt')
        self.py('nested_vote.py', 'clean_ava.tsv', 'elem/clusters_clean.json', 'clean/clean.fa', 'clean_nested.fa',
                log='nested.log', stdout='nested_cut.txt')

    def cluster(self):
        self.search('clean_nested.fa', 'clean_nested.fa', 5, 'nested_ava.tsv', AVA, 'cluster.log')
        if self.a.stage_b:
            self.py('step4.py', 'nested_ava.tsv', 'nested_ava.tsv', 'elem/clusters_red.json', log='cluster.log', stdout='elem/nested.txt')
            self.sh([sys.executable, os.path.join(HERE, 'stage_b.py'), 'expand', 'elem/clusters_raw.json', 'elem/reps.json',
                     'elem/clusters_red.json', 'clean_nested.fa', 'elem/clusters_nested.json', 'exemplar.fa'],
                    'cluster.log', stdout='exemplar.txt')
            self.log('  ' + open(self.path('exemplar.txt')).read().strip())
        else:
            self.py('step4.py', 'nested_ava.tsv', 'nested_ava.tsv', 'elem/clusters_nested.json', log='cluster.log', stdout='elem/nested.txt')
            self.py('exemplar_lib.py', 'elem/clusters_nested.json', 'clean_nested.fa', 'exemplar.fa', log='cluster.log', stdout='exemplar.txt')
        with open(self.path('exemplar.fa')) as fi, open(self.path('exemplar.ids.fa'), 'w') as fo:
            for l in fi: fo.write(l.split()[0] + '\n' if l.startswith('>') else l)

    def flank(self):
        if self.a.no_flank:
            open(self.path('sd_clusters.tsv'), 'w').write('cluster\tmembers\tpairs_sharing_both_flanks\tpairs\tcall\n')
            return
        here = os.path.join(HERE, 'flank_filter.py')
        self.sh([sys.executable, here, 'write', 'elem/clusters_nested.json', self.genome, 'flanks.fa'], 'flank.log')
        self.search('flanks.fa', 'flanks.fa', 5, 'flanks.tsv', AVA, 'flank.log')
        self.sh([sys.executable, here, 'score', 'elem/clusters_nested.json', 'flanks.tsv', 'sd_clusters.tsv'], 'flank.log',
                stdout='flank.txt')
        self.log('  ' + open(self.path('flank.txt')).read().strip())

    def contain(self):
        self.search('exemplar.ids.fa', 'exemplar.ids.fa', 20, 'ex_self.tsv', 'qseqid sseqid pident length qstart qend qlen'.split(), 'contain.log')
        self.py('containment.py', 'elem/clusters_nested.json', 'exemplar.ids.fa', 'ex_self.tsv', 'contained', 0.9,
                log='contain.log', stdout='containment.txt')

    def excise(self):
        shutil.copy(self.path('contained_90.fa'), self.path('r0.fa'))
        for i in (1, 2):
            p = f'r{i - 1}'
            self.search(f'{p}.fa', f'{p}.fa', 20, f'{p}.hits', 'qseqid sseqid pident length qstart qend sstart send qlen slen'.split(), 'excise.log')
            self.py('excise_contained.py', f'{p}.fa', 'elem/clusters_nested.json', f'{p}.hits', f'r{i}.fa', f'cuts{i}.tsv',
                    '--int-cov', 2, '--drop-multi', 2, log='excise.log', stdout=f'excise{i}.txt')
        shutil.copy(self.path('r2.fa'), self.path('final.fa'))

    def pbs(self):
        import csv
        with open(self.path('cand.tsv')) as fi, open(self.path('pbs.tsv'), 'w') as fo:
            fo.write('element\tpbs_mm12\tpbs_trna\tpbs_offset\n')
            for r in csv.DictReader(fi, delimiter='\t'):
                fo.write(f"{r['seq']}:{r['elem_start']}..{r['elem_end']}\t{r['pbs_mm12']}\t{r['pbs_trna']}\t{r['pbs_offset']}\n")

    def tier(self):
        opt = ['--min-id', 90, '--pbs', 'pbs.tsv', '--strict'] if self.a.tier == 'strict' else ['--min-id', 93, '--pbs', 'pbs.tsv']
        self.py('tier.py', 'final.fa', 'elem/clusters_nested.json', 'cand.tsv', 'core.fa', 'candidates.fa', *opt,
                log='tier.log', stdout='tier.txt')

    def models(self):
        a = self.a
        self.sh([sys.executable, os.path.join(HERE, 'consensus_models.py'), 'core.fa', 'lib.fa', 'elem/clusters_nested.json',
                 'models.fa', '--sd', 'sd_clusters.tsv', '--min-cons', a.min_cons if a.min_cons > 0 else 10 ** 9,
                 '--jobs', max(1, a.threads // 2), '--te-composer', a.te_composer, '--work', 'cons',
                 *(['--premade', 'reps.fa'] if a.stage_b else [])], 'models.log', stdout='models.txt')

    def families(self):
        self.search('models.fa', 'models.fa', 5, 'core_ava.tsv', AVA, 'families.log')
        self.py('family.py', 'models.fa', 'core_ava.tsv', 'library.fa', 'families.tsv', '--rule', 'and', '--cand', 'cand.tsv',
                '--gap-aware', log='families.log')

    def run(self):
        a = self.a
        self.genome = os.path.abspath(a.genome)
        if a.sufr_index:
            self.sufr = os.path.abspath(a.sufr_index); open(self.done_path('index'), 'w').close()
        else:
            self.sufr = self.path('genome.sufr')
        t0 = time.time()
        self.step('index', '0 index (sufr)', self.index)
        self.step('sample', '0 depth sample', self.sample)
        self.step('ltr-forage', '1-5 ltr-forage', self.ltr_forage)
        if a.stage_b:
            self.step('precluster', '5b precluster (raw library)', self.precluster)
            self.step('premodels', f'5c premodels (consensus for raw clusters >= {a.min_cons})', self.premodels)
        self.step('depth', '6 depth cleanup' + (' (small clusters + representatives)' if a.stage_b else ''), self.depth)
        self.step('nested', '7 nested vote', self.nested)
        self.step('cluster', '8-9 cluster + exemplar', self.cluster)
        self.step('flank', '9b flank (segmental duplications)', self.flank)
        self.step('contain', '10 containment', self.contain)
        self.step('excise', '11 exemplar nested excision', self.excise)
        self.step('pbs', '12 PBS', self.pbs)
        self.step('tier', f'13 tier ({a.tier})', self.tier)
        self.step('models', f'13b models (consensus for clusters >= {a.min_cons})', self.models)
        self.step('families', '14 families', self.families)
        nrec = sum(1 for l in open(self.path('library.fa')) if l.startswith('>'))
        nfam = len({l.split('\t')[0] for l in open(self.path('families.tsv'))} - {'family'})
        self.log(f'done: {nrec} records, {nfam} families ({time.time() - t0:.0f} s)')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__.split('\n', 2)[2])
    ap.add_argument('genome')
    ap.add_argument('-o', '--outdir', required=True)
    ap.add_argument('-t', '--threads', type=int, default=os.cpu_count())
    ap.add_argument('--sufr', dest='sufr_index', help='existing sufr index of the genome (skips step 0)')
    ap.add_argument('--tier', choices=('pbs', 'strict'), default='pbs')
    ap.add_argument('--from', dest='from_step', choices=STEPS, help='rerun from this step on')
    g = ap.add_argument_group('depth sample and floors')
    g.add_argument('--sample-above', type=float, default=200, help='sample genomes over this many non-N Mb (200)')
    g.add_argument('--sample-frac', type=float, default=0.1)
    g.add_argument('--window-kb', type=int, default=100)
    g.add_argument('--seed', type=int, default=1)
    g.add_argument('--sample-file', help='use this FASTA as the depth sample')
    g.add_argument('--min-depth', type=int, help='default 50 x genome non-N Mb / 373, capped at 50 (rice MSU7: 50)')
    g.add_argument('--abs-depth', type=int, help='default 10 x min depth')
    m = ap.add_argument_group('family models')
    m.add_argument('--min-cons', type=int, default=3, help='consensus for clusters of at least this many elements (3; 0: exemplars only)')
    m.add_argument('--stage-a', action='store_true',
                   help='depth-clean and re-cluster every element (the older flow) instead of clustering the raw library '
                        'first and representing clusters >= --min-cons by their consensus (stage B, the default)')
    m.add_argument('--stage-b', action='store_true', help=argparse.SUPPRESS)   # the default now; kept for old commands
    m.add_argument('--no-flank', action='store_true', help='skip the segmental-duplication flank filter')
    t = ap.add_argument_group('programs and data')
    t.add_argument('--ltr-forage', dest='ltr_forage')
    t.add_argument('--sufr-bin', dest='sufr')
    t.add_argument('--rmblastn')
    t.add_argument('--te-composer', dest='te_composer')
    t.add_argument('--domain-hmm', default=os.environ.get('HMM', '/usr/local/LTR_retriever-2.9.0/database/TEfam.hmm'))
    t.add_argument('--word', type=int, default=12, help='rmblastn seed word for the library all-vs-all (12)')
    a = ap.parse_args()
    a.ltr_forage = find_tool(a.ltr_forage, 'LTRFORAGE', 'ltr-forage', [os.path.join(ROOT, 'target', 'release', 'ltr-forage')])
    a.rmblastn = find_tool(a.rmblastn, 'RMBLASTN', 'rmblastn',
                           [os.path.join(ROOT, '..', '..', 'rmblast-port', 'rmblast', 'target', 'release', 'rmblastn')],
                           check='Rust port')
    if a.stage_a and a.stage_b: sys.exit('ltr_library: --stage-a and --stage-b are exclusive')
    a.stage_b = not a.stage_a and a.min_cons >= 2   # stage B needs consensus representatives
    if a.min_cons > 0:
        a.te_composer = find_tool(a.te_composer, 'TE_COMPOSER', 'te-composer',
                                  [os.path.join(ROOT, '..', '..', 'dfam-curator', 'target', 'release', 'te-composer')])
    if not a.sufr_index:
        a.sufr = find_tool(a.sufr, 'SUFR', 'sufr', [os.path.join(ROOT, '..', 'sufr', 'target', 'release', 'sufr')])
    if not os.path.isfile(a.domain_hmm):
        sys.exit(f'ltr_library: domain HMM file not found: {a.domain_hmm} (give --domain-hmm)')
    Run(a).run()


if __name__ == '__main__':
    main()
