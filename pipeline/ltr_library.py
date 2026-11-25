#!/usr/bin/env python3
"""Build an LTR retrotransposon family library from a genome assembly with ltr-forage.

usage: ltr_library.py <genome.fa> -o <outdir> [-t THREADS] [options]

Steps (they work in <outdir>/work, their commands' logs in <outdir>/work/logs; <outdir>/run.log has every step's
time and results, the screen a summary (-v: all of it); the numbers continue ltr-forage's log steps 1-11 and match
the boxes of the pipeline figure):
   0  sample     depth sample: the whole genome up to --sample-above Mb, else a random
                 --sample-frac of its --window-kb windows (seeded); depth floors are scaled to it
   1  index      sufr suffix array of the genome (skipped with --sufr)
1-11  ltr-forage detect, re-extend, annotate (domains, TSD, PBS), decide, raw library; tandem arrays
                 found by ltr-forage's FasTAN port (satellites, simple repeats, tandem copies; any unit up to
                 30 kb) are masked from the seed search (tandem.bed; --no-fastan), and --mask-bed adds
                 intervals of one's own.  TSD and TG..CA are evidence at the final boundaries (the TSD
                 search does not move the ends, as LTRharvest's does).  Calls whose re-extended ends agree
                 within 20 bp are merged, and overlapping calls are resolved by --overlaps
  12  precluster element clusters of the raw library (all-vs-all, as step 16).  A probe of 150 records of
                 each kind projects the all-vs-all hits first; above --greedy-above million (about 1.5-2.5 min
                 per million: mm10 45 M, 60-110 min; maize 2,365 M) the clustering is greedy instead
                 (greedy_cluster.py): batches of --greedy-batch elements all-vs-all, the rest searched
                 against at most --greedy-reps members per cluster
  13  premodels  te-composer consensus for raw clusters of >= --min-cons elements; from here on
                 each big cluster is represented by its consensus, and steps 14-16 work on the elements of
                 smaller clusters plus those representatives (the members of a big cluster join the
                 cluster of its representative after step 16).  12-13 are stage B, the default;
                 --stage-a (or --min-cons < 2) skips them and steps 14-16 work on every element
  14  depth      three rounds: rmblastn of the sample against the library, excise stretches far
                 above their record's own genome depth (nested passengers), drop high-depth LTRs.
                 Only the first round is a search: later rounds apply the last round's cuts to its
                 alignments (--no-depth-remap searches the cut records again instead).
                 A representative is cut only at 10x its own depth and never dropped: a big family's
                 consensus is high-copy itself, while a real insertion stands far above it
  15  nested     first element clustering (library all-vs-all), nested LTR excision by vote
  16  cluster    element clusters and one exemplar element per cluster
  17  flank      flag clusters whose copies share both flanks (segmental duplications; --no-flank):
                 most pairs share both flanks, or every copy shares both with some other copy
  18  contain    drop exemplars contained in another at >= 0.9
  19  excise     nested LTR excision from internals, two rounds
  20  pbs        primer binding site per element (ltr-forage's tRNA search)
  21  tier       keep multi-element clusters; singletons with LTR identity >= 93 or a PBS
                 (--tier strict: identity >= 90 and a PBS; --tier all: every singleton)
  22  models     te-composer consensus for clusters of >= --min-cons elements (oriented, columns held
                 by under half the median occupancy dropped; the exemplar where none can be built), the
                 exemplar for smaller ones; a stage B representative is kept as steps 14-19 left it;
                 flagged segmental-duplication clusters left out (--min-cons 0: exemplars only)
 22b  mdepth     one more depth round on the finished models, as step 14: a consensus model (stage B
                 representative, or one step 22 built) is cut at 10x its own depth, an exemplar at 5x (--model-depth-ratio; at 3x
                 conserved stretches of old families were cut: ath Gypsy-8_ALY, 3.1x); none is
                 dropped (without step 14's competition a singleton exemplar of a high-copy LTR family would
                 pass --abs-depth whole).  Step 14 runs with mask level 80,
                 so a non-LTR passenger shared by many records splits its sample copies among them and
                 each sees a low depth; a consensus model can also take back a passenger its founder
                 lost (mm10: L1 and B1/B2 pieces in five models masked 191 Mb outside annotated LTRs).
                 --no-model-depth skips it.  With --early-tier, an internal model with more than
                 --int-copies-max (2000) genome copies per element of its cluster is dropped: internals
                 are never solo, so it is mostly another repeat (mm10 --early-tier: an L1 exemplar
                 internal, 5,999 per element, 133 Mb outside the curated LTRs).  Without --early-tier the
                 check is off unless the option is given: it changed none of the six tuning genomes'
                 default runs.  The depth
                 search is then repeated on the cut models: an LTR model whose median depth rose
                 --depth-jump-fold (20) and whose element has neither a TSD nor TG..CA ends took over a
                 cut segment's copies and is dropped (mm10: two LTR models of L1 sequence, 101x and 60x,
                 44 Mb outside the curated LTRs; zebra finch fam0093, 170x with TSD and TG..CA, is real)
  23  families   80-80 families (LTRs and internals both), gap-aware, Gypsy/Copia class names; then each record
                 turned to its element's coding strand (orient_models.py: protein domains, else the PBS; families.tsv
                 orientation and flipped)
 23b  seeds      library.stk, a seed alignment for every record: its best --seed-copies non-overlapping
                 copies in the depth sample (the record's origin element when the sample has none),
                 aligned to the unchanged record by te-composer --fixed-cons, so #=GC RF is the record.
                 A te-composer without --fixed-cons (before 0.6.8) gives one-row alignments (the record alone)
  24  evidence   --evidence-page: evidence.html, one dossier per family from its founding element
                 (structure, domains, PBS, copies, verdict; pipeline/evidence_page.py, which also takes a
                 genome mask, a curated annotation and another tool's mask when run by hand)

Outputs in <outdir>: library.fa (fam<N>_v<k>_LTR|INT#LTR/<class>), library.stk, families.tsv, cand.tsv (every
candidate element, ltr-forage's table), evidence.html with --evidence-page, and run.log (the steps, their results, and
the programs, versions and paths used).  A successful run removes <outdir>/work unless --keep-work.  A rerun skips
finished steps (<outdir>/work/.done); --from STEP reruns from a step on (both need work/).

Stage B against stage A (2026-10-20; sensitivity / precision of LTR masking against curated libraries):
rice 92.70 / 98.53 vs 92.51 / 98.29, zebra finch 74.76 / 83.03 vs 74.67 / 83.45 with steps 14-23 (then numbered 6-14) in 16 min
instead of ~109, gorilla 34.57 / 97.65 (HiRepEat 17.68 / 95.97).  Most of zebra finch's remaining FP is
outside the curated library and much of it looks like real families HiRepEat also finds: FP rules need
an artifact mechanism, not lower novelty.

External programs: ltr-forage, sufr, rmblastn (Rust), te-composer (dfam-curator >= 0.6.6), and nail
(>= 0.7.1) with MMseqs2 for ltr-forage's domain search, found from the options, the environment (LTRFORAGE,
SUFR, RMBLASTN, TE_COMPOSER, NAIL, MMSEQS) or PATH; Python 3 with numpy."""
import argparse, concurrent.futures, glob, os, random, re, shutil, subprocess, sys, tempfile, time

HERE = os.path.dirname(os.path.realpath(__file__))   # through symlinks (an LTR_FORAGE_DIR of links)
ROOT = os.path.dirname(HERE)
STEPS = ['sample', 'index', 'ltr-forage', 'polya-family', 'precluster', 'premodels', 'depth', 'nested', 'cluster', 'flank', 'contain', 'excise', 'pbs', 'tier',
         'models', 'mdepth', 'families', 'seeds', 'evidence']
AVA = 'qseqid sseqid pident length qstart qend sstart send qlen slen bitscore'.split()
# the run's products, moved from <outdir>/work to <outdir> at the end (copied with --keep-work)
PRODUCTS = ('library.fa', 'library.stk', 'families.tsv', 'cand.tsv', 'evidence.html')
# the screen's headings: each step's group and, in the library building, its stage range
PREP, DETECT, BUILD = 'Sequence & Suffix Array Preparation (stages 0-1)', 'LTR Detection: ltr-forage (stages 1-11)', 'Library Building'
GROUPS = {'sample': (PREP, None), 'index': (PREP, None), 'ltr-forage': (DETECT, None), 'polya-family': (DETECT, None),
          **{s: (BUILD, 'stages 12-13: clustering and consensus building') for s in ('precluster', 'premodels')},
          **{s: (BUILD, 'stages 14-22: artifact cleanup') for s in
             ('depth', 'nested', 'cluster', 'flank', 'contain', 'excise', 'pbs', 'tier', 'models', 'mdepth')},
          **{s: (BUILD, 'stage    23: families and seed alignments') for s in ('families', 'seeds')},
          'evidence': (BUILD, 'stage    24: evidence page')}
# the options each step's output depends on (argparse names); a finished step is not skipped when one changed
STEP_OPTIONS = {
    'sample': ['sample_above', 'sample_frac', 'window_kb', 'seed', 'sample_file', 'min_depth', 'abs_depth'],
    'index': [],
    'ltr-forage': ['no_fastan', 'overlaps', 'forage_arg', 'mask_bed'],   # programs and profiles: in run.log
    'polya-family': ['polya_family', 'polya_pool', 'no_polya_rescue', 'seed'],
    'precluster': ['greedy_above', 'greedy_batch', 'greedy_reps', 'word'],
    'premodels': ['early_tier', 'tier', 'min_cons'],
    'depth': ['depth_rounds', 'no_depth_remap', 'depth_word'],
    'nested': ['word'], 'cluster': ['word'], 'flank': ['no_flank'], 'contain': ['word'], 'excise': ['word'], 'pbs': [],
    'tier': ['tier'],
    'models': ['min_cons'],
    'mdepth': ['no_model_depth', 'model_depth_ratio', 'int_copies_max', 'depth_jump_fold', 'depth_word', 'polya_family'],
    'families': ['word'],
    'seeds': ['seed_copies', 'mask_warn'],
    'evidence': ['evidence_reference'],
}


def tool_version(path, flag='--version'):
    """The first line a program prints for its version flag ('' if it prints none)."""
    try:
        r = subprocess.run([path, flag], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return ''
    return ((r.stdout or r.stderr).strip().split('\n') or [''])[0].strip()


def version_tuple(line):
    m = re.search(r'(\d+)\.(\d+)\.(\d+)', line)
    return tuple(map(int, m.groups())) if m else None


def find_tool(opt, env, name, check=None, min_version=None, flag='--version'):
    """First executable from the option, the environment variable, then PATH (no built-in install
    locations); with check, only one whose --help output contains that text (NCBI's rmblastn is often
    on PATH too, and only the Rust port takes the options used here); with min_version, only one at
    least that version."""
    seen = []
    for p in (opt, os.environ.get(env), shutil.which(name)):
        if not (p and os.path.isfile(p) and os.access(p, os.X_OK)):
            continue
        if check and check not in subprocess.run([p, '--help'], capture_output=True, text=True).stdout:
            continue
        if min_version:
            v = version_tuple(tool_version(p, flag))
            if v is None or v < min_version:
                seen.append(f'{p} ({".".join(map(str, v)) if v else "unknown version"})'); continue
        return os.path.abspath(p)
    need = f' >= {".".join(map(str, min_version))}' if min_version else ''
    sys.exit(f'ltr_library: cannot find {name}{need}' + (f' ({check})' if check else '') +
             (f'; too old: {", ".join(seen)}' if seen else '') + f'; give --{name.replace("_", "-")} or set {env}')


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


def genome_stats(path, check_names=False):
    """Total and N bases.  check_names: the genome's sequence names must be unique and free of '#', which
    separates a library record's name from its class (every helper keys records on name.split('#')[0])."""
    total = ns = 0; seen = set()
    for name, s in read_fasta(path):
        if check_names:
            if '#' in name:
                sys.exit(f"ltr_library: sequence name '{name}' in {path} contains '#', which library record names use "
                         f"as the class separator; rename the sequences (e.g. sed '/^>/s/#/_/g')")
            if name in seen:
                sys.exit(f"ltr_library: sequence name '{name}' occurs twice in {path}")
            seen.add(name)
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
        # the products and run.log go to the outdir; every step works in <outdir>/work, which a successful run
        # removes unless --keep-work (a stopped run keeps it, so a rerun resumes)
        self.top = os.path.abspath(a.outdir)
        self.out = os.path.join(self.top, 'work')
        if os.path.isdir(os.path.join(self.top, '.done')):
            sys.exit(f'ltr_library: {self.top} was made by an earlier ltr_library.py, with the intermediate files beside '
                     f'the library; use a new outdir')
        if not os.path.isdir(self.out) and os.path.exists(os.path.join(self.top, 'library.fa')):
            sys.exit(f'ltr_library: {self.top} holds a finished run whose intermediate files were removed; use a new outdir'
                     + (' (--keep-work keeps them, for a later --from)' if a.from_step else ''))
        self.group = self.sub = None; self.skipped = []
        os.makedirs(os.path.join(self.out, 'logs'), exist_ok=True)
        os.makedirs(os.path.join(self.out, 'elem'), exist_ok=True)
        os.makedirs(os.path.join(self.out, '.done'), exist_ok=True)
        # temporary files of every step go under the outdir, not the system /tmp, unless --tmpdir names another
        # place: sufr's index partitions need up to the index's size (C. elegans peak 0.58 GB, index 0.90 GB), sufr
        # does not check its writes, and a full /tmp gives a corrupt index without an error
        self.tmp = os.path.abspath(a.tmpdir) if a.tmpdir else os.path.join(self.out, 'tmp')
        os.makedirs(self.tmp, exist_ok=True)
        os.environ['TMPDIR'] = self.tmp; tempfile.tempdir = self.tmp
        self.env = dict(os.environ, TSV=os.path.join(self.out, 'cand.tsv'), RMBLASTN=a.rmblastn, WORD=str(a.word))
        redo = STEPS[STEPS.index(a.from_step):] if a.from_step else []
        # stage A and B outputs from step 12 on don't mix: a rerun keeps the mode its finished steps were made in.
        # Checked before --from removes any marker, so a refused rerun leaves the directory as it was
        mode = 'B' if a.stage_b else 'A'
        mf = self.done_path('mode')
        made = open(mf).read().strip() if os.path.exists(mf) else ('B' if os.path.exists(self.path('reps.fa')) else 'A')
        later = [s for s in STEPS[STEPS.index('precluster'):] if s not in redo and os.path.exists(self.done_path(s))]
        if later and made != mode:
            sys.exit(f'ltr_library: {self.out} has steps from stage {made} ({later[0]} ...); rerun with '
                     f'{"--stage-a" if made == "A" else "no --stage-a"}, or --from precluster to redo them as stage {mode}')
        # one genome per outdir: its path, size and modification time are kept with the markers
        self.genome = os.path.abspath(a.genome)
        st = os.stat(self.genome)
        ident = f'{self.genome}\t{st.st_size}\t{st.st_mtime_ns}\n'
        gf = self.done_path('genome')
        kept = [s for s in STEPS if s not in redo and os.path.exists(self.done_path(s))]
        if os.path.exists(gf) and open(gf).read() != ident and kept:
            was = open(gf).read().rstrip('\n').split('\t')
            sys.exit(f'ltr_library: {self.out} holds finished steps ({kept[0]} ...) for the genome {was[0]} '
                     f'({was[1]} bytes); now {self.genome} ({st.st_size} bytes'
                     f'{", modified since" if was[0] == self.genome and was[1] == str(st.st_size) else ""}). '
                     f'Use another outdir, or --from sample to redo every step')
        for s in redo:
            p = self.done_path(s)
            if os.path.exists(p): os.remove(p)
        open(gf, 'w').write(ident)
        open(mf, 'w').write(mode + '\n')

    def done_path(self, step): return os.path.join(self.out, '.done', step)
    def path(self, *p): return os.path.join(self.out, *p)

    def log(self, msg, level='detail', stamp=None):
        """One line to <outdir>/run.log, and to the screen for level 'user' and 'warn' (and 'detail' with --verbose).
        stamp: the line's own time (a program's log line), else now."""
        line = f'{stamp or time.strftime("[%H:%M:%S]")} {msg}'
        if level != 'detail' or self.a.verbose: print(line, flush=True)
        with open(os.path.join(self.top, 'run.log'), 'a') as fh: fh.write(line + '\n')

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

    def py(self, script, *args, log, stdout=None):
        self.sh([sys.executable, os.path.join(HERE, script), *args], log, stdout)

    def search(self, q, db, max_hsps, out, fields, log, clustering=False, by_kind=False):
        """rmb_ava.sh of q against db.  by_kind: two searches, LTR records against LTR records and internals
        against internals (q must be db), for readers that skip LTR-INT hits (step4.py, family.py).
        clustering (implies by_kind): also drop hits below 80% identity as they are written; step4.py reads
        only same-kind hits of >= 80% (mm10 raw library: 80.5 M hits, 9.7 GB, of which 45 M kept)."""
        env = dict(self.env)
        if clustering: self.env['MIN_PID'] = '80'; by_kind = True
        try:
            if not by_kind:
                self.sh(['bash', os.path.join(HERE, 'rmb_ava.sh'), q, db, max_hsps, out, *fields], log); return
            assert q == db
            base = os.path.splitext(out)[0]
            # the matrix of the whole library, as one search would pick (LTRs and internals differ in GC)
            if 'MATRIX' not in self.env:
                self.env['MATRIX'] = subprocess.run(['bash', os.path.join(HERE, 'pick_matrix.sh'), q], cwd=self.out,
                                                    capture_output=True, text=True, check=True).stdout.strip()
            with open(self.path(base + '.LTR.fa'), 'w') as fl, open(self.path(base + '.INT.fa'), 'w') as fi:
                for n, x in read_fasta(self.path(q)):
                    (fi if n.split('#')[0].endswith('_INT') else fl).write(f'>{n}\n{x}\n')
            for k in ('LTR', 'INT'):
                f = f'{base}.{k}.fa'
                self.sh(['bash', os.path.join(HERE, 'rmb_ava.sh'), f, f, max_hsps, f'{base}.{k}.tsv', *fields], log)
            with open(self.path(out), 'w') as fo:
                for k in ('LTR', 'INT'):
                    with open(self.path(f'{base}.{k}.tsv')) as fh: shutil.copyfileobj(fh, fo)
                    os.remove(self.path(f'{base}.{k}.tsv')); os.remove(self.path(f'{base}.{k}.fa'))
        finally: self.env = env

    def signature(self, name):
        """The options a step's output depends on, one 'option<TAB>value' line each (kept in its done-marker)."""
        a = self.a
        return ''.join(f'{k}\t{getattr(a, k)!r}\n' for k in STEP_OPTIONS[name]
                       if not (k == 'tier' and name == 'premodels' and not a.early_tier))

    def step(self, name, title, fn):
        p, sig = self.done_path(name), self.signature(name)
        if os.path.exists(p):
            # finished: skipped, unless it was made with other options (markers of older versions are empty)
            old = dict(l.rstrip('\n').split('\t', 1) for l in open(p) if '\t' in l)
            new = dict(l.rstrip('\n').split('\t', 1) for l in sig.splitlines())
            diff = [f'--{k.replace("_", "-")} {old[k]} -> {new[k]}' for k in new if k in old and old[k] != new[k]]
            if diff:
                sys.exit(f'ltr_library: step {name} in {self.out} was made with other options ({"; ".join(diff)}). '
                         f'Rerun with --from {name} to redo it and the steps after it, or with the options it was made with')
            if not (name == 'index' and self.a.sufr_index): self.skipped.append(name)
            return
        if self.skipped and self.group is None:
            self.log(f'Resuming: {len(self.skipped)} finished step(s) kept ({self.skipped[0]} to {self.skipped[-1]}); '
                     f'--from STEP redoes from a step', 'user')
        group, sub = GROUPS[name]
        if group != self.group: self.log(group, 'user'); self.group, self.sub = group, None
        if sub and sub != self.sub: self.log('  ' + sub, 'user'); self.sub = sub
        self.log(title)
        fn()
        open(p, 'w').write(sig)

    # ---- steps -----------------------------------------------------------
    def index(self):
        self.sh([self.a.sufr, '-t', self.a.threads, 'create', '--dna', '-o', self.path('genome.sufr'), self.genome], 'index.log')

    def sample(self):
        a = self.a
        total, ns = genome_stats(self.genome, check_names=True)
        if total == 0:
            sys.exit(f'ltr_library: no sequence in {self.genome}')
        mb = (total - ns) / 1e6
        # an earlier run's sample goes first: for a small genome it is a symlink to the genome, which a copy or
        # a new sample written to the same name would follow and overwrite
        ds = self.path('depth_sample.fa')
        if os.path.lexists(ds): os.remove(ds)
        if a.sample_file:
            frac = genome_stats(a.sample_file)[0] / total; src = os.path.abspath(a.sample_file); kind = 'file'
            shutil.copy(src, ds)
        elif mb <= a.sample_above:
            frac = 1.0; os.symlink(self.genome, ds); kind = 'genome'
        else:
            kept = write_sample(self.genome, ds, a.sample_frac, a.window_kb * 1000, a.seed); kind = 'windows'
            frac = kept / total
        # depth floors in whole-genome units, scaled with the non-N genome size (rice MSU7,
        # 373 Mb: 50 / 500) and capped at 50: every tuned genome used <= 50, and lotus (708 Mb,
        # uncapped 95 / 950) lost 2 points of precision (90.35 vs 92.36 at 50 / 500) because higher
        # floors trim less; excise_by_depth.py multiplies them by the sample's --scale
        mind = a.min_depth if a.min_depth is not None else min(50, max(5, round(50 * mb / 373)))
        absd = a.abs_depth if a.abs_depth is not None else 10 * mind
        with open(self.path('sample.txt'), 'w') as fh:
            fh.write(f'genome_mb\t{mb:.1f}\nn_frac\t{ns / total:.4f}\nscale\t{frac:.4f}\nmin_depth\t{mind}\nabs_depth\t{absd}\n'
                     f'kind\t{kind}\n')
        self.log(f'  genome {mb:.1f} Mb non-N; depth sample {frac:.1%} of the genome; floors {mind} / {absd}', 'user')

    def params(self):
        return dict(l.rstrip('\n').split('\t') for l in open(self.path('sample.txt')))

    def drop_full_copies(self, restore_lib=False):
        """Remove early_tier.py's copies of the full raw library (<file>.all), which it otherwise copies back over
        files a redone step 3 or 12 has just written.  restore_lib: lib.fa, step 12's input, first gets its full
        copy back (an --early-tier run left it without the dropped singletons)."""
        full = self.path('lib.fa.all')
        if restore_lib and os.path.exists(full): shutil.copy(full, self.path('lib.fa'))
        for f in ('lib.fa', 'exemplar_raw.fa', os.path.join('elem', 'clusters_raw.json')):
            if os.path.exists(self.path(f + '.all')): os.remove(self.path(f + '.all'))

    def ltr_forage(self):
        self.drop_full_copies()
        self.sh([self.a.ltr_forage, self.sufr, '-v', '--library', 'lib.raw.fa', '--threads', self.a.threads,
                 '--domain-hmm', self.a.domain_hmm, '--nail', self.a.nail, '--mmseqs', self.a.mmseqs, '--temp-dir', self.out,
                 '--overlaps', self.a.overlaps, *(['--no-fastan'] if self.a.no_fastan else ['--fastan-out', 'tandem.bed']),
                 *(['--mask-bed', os.path.abspath(self.a.mask_bed)] if self.a.mask_bed else []), *self.a.forage_arg],
                'ltr-forage.log', stdout='cand.tsv')
        with open(self.path('lib.raw.fa')) as fi, open(self.path('lib.fa'), 'w') as fo:
            for l in fi: fo.write(l.split()[0] + '\n' if l.startswith('>') else l)
        # ltr-forage's own log (this run's: after the last command line of its log) goes to run.log, indented
        lines = open(self.path('logs', 'ltr-forage.log')).read().split('\n')
        last = max((i for i, l in enumerate(lines) if l.startswith('$ ')), default=-1)
        for l in lines[last + 1:]:
            m = re.match(r'(\[\d\d:\d\d:\d\d\]) (.*)', l)
            if m: self.log('  ' + m.group(2), stamp=m.group(1))

    def polya_family(self):
        """Step 11b: veto kept elements whose LTRs are copies of the poly-A-terminated repeats ltr-forage vetoed."""
        self.drop_full_copies(restore_lib=True)
        self.py('polya_family.py', '.', self.genome, '--min-depth', self.a.polya_family, '--min-pool', self.a.polya_pool, '--seed', self.a.seed,
                *(['--no-rescue'] if self.a.no_polya_rescue else []), log='polya_family.log', stdout='polya_family.txt')
        self.log('  ' + open(self.path('polya_family.txt')).read().strip().split('\n')[-1])
        # the user's line: candidates (ltr-forage's step 9) and the elements left after the veto (one INT record each)
        cands = re.findall(r'\[9 decision\] kept \d+ of (\d+) candidate', open(self.path('logs', 'ltr-forage.log')).read())
        kept = sum(1 for n, _ in read_fasta(self.path('lib.fa')) if n.split('#')[0].endswith('_INT'))
        self.log(f'  {kept:,} LTR retrotransposon copies kept' + (f' of {int(cands[-1]):,} candidates' if cands else '')
                 + ' (cand.tsv)', 'user')

    def precluster(self):
        if not self.a.stage_b: return
        a = self.a
        self.drop_full_copies(restore_lib=True)
        if a.greedy_above >= 0:
            self.py('greedy_cluster.py', 'probe', 'lib.fa', '--work', 'greedy_probe', log='precluster.log', stdout='elem/probe.txt')
            proj = float(re.search(r'projected all-vs-all hits: ([0-9.]+) M', open(self.path('elem', 'probe.txt')).read()).group(1))
            self.log(f'  probe: {proj:.1f} M all-vs-all hits projected')
        if a.greedy_above >= 0 and proj > a.greedy_above:
            self.log(f'  greedy clustering (above {a.greedy_above:g} M): batches of {a.greedy_batch}, at most {a.greedy_reps} representatives per cluster', 'user')
            self.py('greedy_cluster.py', 'cluster', 'lib.fa', 'elem/clusters_raw.json', '--batch', a.greedy_batch,
                    '--reps', a.greedy_reps, '--work', 'greedy_work', log='precluster.log', stdout='elem/raw.txt')
        else:
            self.search('lib.fa', 'lib.fa', 5, 'raw_ava.tsv', AVA, 'precluster.log', clustering=True)
            self.py('step4.py', 'raw_ava.tsv', 'raw_ava.tsv', 'elem/clusters_raw.json', log='precluster.log', stdout='elem/raw.txt')
        self.py('exemplar_lib.py', 'elem/clusters_raw.json', 'lib.fa', 'exemplar_raw.fa', log='precluster.log')
        self.log('  ' + open(self.path('elem', 'raw.txt')).read().strip().split('\n')[-1])

    def premodels(self):
        if not self.a.stage_b: return
        a = self.a
        # --early-tier: raw singletons the tier would reject leave before the depth cleanup (without it the
        # full raw library is restored, so the option can be switched on a --from premodels rerun)
        self.py('early_tier.py', '.', a.tier if a.early_tier else 'all', log='premodels.log', stdout='early_tier.txt')
        if a.early_tier: self.log('  ' + open(self.path('early_tier.txt')).read().strip())
        self.sh([sys.executable, os.path.join(HERE, 'consensus_models.py'), 'exemplar_raw.fa', 'lib.fa', 'elem/clusters_raw.json',
                 'reps.fa', '--big-only', '--min-cons', a.min_cons, '--jobs', max(1, a.threads // 2),
                 '--te-composer', a.te_composer, '--work', 'cons'], 'premodels.log', stdout='premodels.txt')
        self.sh([sys.executable, os.path.join(HERE, 'stage_b.py'), 'split', 'elem/clusters_raw.json', 'lib.fa', 'reps.fa',
                 'depth_lib.fa', 'elem/reps.json', '--min-cons', a.min_cons], 'premodels.log', stdout='split.txt')
        self.log('  ' + open(self.path('split.txt')).read().strip())

    def depth(self):
        p = self.params()
        mtx = subprocess.run(['bash', os.path.join(HERE, 'pick_matrix.sh'), 'depth_sample.fa'], cwd=self.out,
                             capture_output=True, text=True, check=True).stdout.strip()
        # one single-thread batch per thread (was 16 batches of threads // 16: at -t 24, 8 cores idle); every sample
        # sequence is searched on its own, so the hits do not depend on the batching
        self.env['BATCHES'] = str(self.a.threads); self.env['DEPTH_WORD'] = str(self.a.depth_word)
        # rounds after the first apply the last round's cuts to its hits instead of searching again (--no-depth-remap)
        self.env['REMAP'] = '0' if self.a.no_depth_remap else '1'
        self.sh(['bash', os.path.join(HERE, 'clean_library_rmb.sh'), 'depth_lib.fa' if self.a.stage_b else 'lib.fa', 'depth_sample.fa', 'clean', self.a.depth_rounds, p['scale'], mtx,
                 '--ratio', 3, '--min-depth', p['min_depth'], '--min-len', 50, '--abs-depth', p['abs_depth'], '--drop-frac', 0.7,
                 *(['--consensus', 'reps.fa'] if self.a.stage_b else [])],
                'depth.log', threads=1)
        del self.env['REMAP']
        self.drop_raw_hits('clean')
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
        self.search('clean_nested.fa', 'clean_nested.fa', 5, 'nested_ava.tsv', AVA, 'cluster.log', clustering=True)
        if self.a.stage_b:
            self.py('step4.py', 'nested_ava.tsv', 'nested_ava.tsv', 'elem/clusters_red.json', '--reps', 'elem/reps.json', log='cluster.log', stdout='elem/nested.txt')
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
        # within-cluster searches only (flank_filter.py bins); the matrix is the one a single search would pick
        bd = self.path('flank_bins')
        if os.path.isdir(bd): shutil.rmtree(bd)
        self.sh([sys.executable, here, 'bins', 'elem/clusters_nested.json', 'flanks.fa', 'flank_bins'], 'flank.log')
        env = dict(self.env)
        try:
            self.env['MIN_PID'] = '80'
            self.env['MATRIX'] = subprocess.run(['bash', os.path.join(HERE, 'pick_matrix.sh'), 'flanks.fa'], cwd=self.out,
                                                capture_output=True, text=True, check=True).stdout.strip()
            with open(self.path('flanks.tsv'), 'w') as fo:
                for f in sorted(os.listdir(bd)):
                    if not f.endswith('.fa'): continue
                    b = os.path.join('flank_bins', f)
                    self.sh(['bash', os.path.join(HERE, 'rmb_ava.sh'), b, b, 5, b[:-3] + '.tsv', *AVA], 'flank.log')
                    with open(self.path(b[:-3] + '.tsv')) as fh: shutil.copyfileobj(fh, fo)
        finally: self.env = env
        shutil.rmtree(bd)
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
        opt = {'strict': ['--min-id', 90, '--pbs', 'pbs.tsv', '--strict'], 'pbs': ['--min-id', 93, '--pbs', 'pbs.tsv'],
               'all': ['--min-id', 0, '--max-len', 10 ** 9]}[self.a.tier]
        self.py('tier.py', 'final.fa', 'elem/clusters_nested.json', 'cand.tsv', 'core.fa', 'candidates.fa', *opt,
                log='tier.log', stdout='tier.txt')

    def models(self):
        a = self.a
        self.sh([sys.executable, os.path.join(HERE, 'consensus_models.py'), 'core.fa', 'lib.fa', 'elem/clusters_nested.json',
                 'models.nodepth.fa', '--sd', 'sd_clusters.tsv', '--min-cons', a.min_cons if a.min_cons > 0 else 10 ** 9,
                 '--jobs', max(1, a.threads // 2), '--te-composer', a.te_composer, '--work', 'cons',
                 *(['--premade', 'reps.fa'] if a.stage_b else [])], 'models.log', stdout='models.txt')

    def drop_raw_hits(self, *dirs):
        """The unfiltered hits of a finished depth step (gigabytes: every sample locus against every record it matches)."""
        for d in dirs:
            for f in glob.glob(self.path(d, 'batches', 'raw*.tsv')): os.remove(f)

    def mdepth(self):
        src = self.polya_cut()
        if self.a.no_model_depth:
            shutil.copy(self.path(src), self.path('models.fa')); return
        self.mdepth_passes(src)
        self.drop_raw_hits('mclean', 'mclean2', 'mclean3')

    def polya_cut(self):
        """The models without what step 11b's sampled poly-A-terminated repeats cover (polya_cut.py); the file the depth
        check reads."""
        if not self.a.polya_family or not os.path.exists(self.path('polya_family', 'polya.fa')): return 'models.nodepth.fa'
        self.py('polya_cut.py', 'models.nodepth.fa', 'polya_family/polya.fa', 'models.polyacut.fa', 'polya_family',
                '--min-depth', self.a.polya_family, '--cuts', 'polya_family/cuts.tsv', log='mdepth.log', stdout='polya_cut.txt')
        self.log('  ' + open(self.path('polya_cut.txt')).read().strip().split('\n')[-1])
        return 'models.polyacut.fa'

    def mdepth_passes(self, src='models.nodepth.fa'):
        p = self.params()
        mtx = subprocess.run(['bash', os.path.join(HERE, 'pick_matrix.sh'), 'depth_sample.fa'], cwd=self.out,
                             capture_output=True, text=True, check=True).stdout.strip()
        # consensus models: stage B representatives and records step 22 replaced by a consensus
        core = dict(read_fasta(self.path('core.fa')))
        reps = {n.split('#')[0] for n, _ in read_fasta(self.path('reps.fa'))} if self.a.stage_b else set()
        with open(self.path('models.cons.fa'), 'w') as fh:
            for n, s in read_fasta(self.path('models.nodepth.fa')):
                k = n.split('#')[0]
                if k in reps or core.get(n, '').upper() != s.upper(): fh.write(f'>{k}\n{s}\n')
        self.env['BATCHES'] = str(self.a.threads); self.env['DEPTH_WORD'] = str(self.a.depth_word)
        self.sh(['bash', os.path.join(HERE, 'clean_library_rmb.sh'), src, 'depth_sample.fa', 'mclean', 1, p['scale'],
                 mtx, '--ratio', self.a.model_depth_ratio, '--min-depth', p['min_depth'], '--min-len', 50, '--consensus', 'models.cons.fa'],
                'mdepth.log', threads=1)
        shutil.copy(self.path('mclean', 'clean.fa'), self.path('models.fa'))
        self.log('  ' + open(self.path('mclean', 'round1.log')).read().strip().split('\n')[-1])
        if self.a.int_copies_max:
            self.py('int_copies.py', 'models.fa', 'mclean/depth1.tsv', 'elem/clusters_nested.json', p['scale'],
                    self.a.int_copies_max, 'int_copies.tsv', log='mdepth.log', stdout='int_copies.txt')
            self.log('  ' + open(self.path('int_copies.txt')).read().strip())
        if self.a.depth_jump_fold:
            # the trimmed models differ from the first search's in a few records: only those are searched again
            self.env['PREV_RUN'] = 'mclean'
            self.sh(['bash', os.path.join(HERE, 'clean_library_rmb.sh'), 'models.fa', 'depth_sample.fa', 'mclean2', 1, p['scale'],
                     mtx, '--ratio', self.a.model_depth_ratio, '--min-depth', p['min_depth'], '--min-len', 50, '--consensus', 'models.cons.fa'],
                    'mdepth.log', threads=1)
            del self.env['PREV_RUN']
            self.py('depth_jump.py', 'models.fa', 'mclean/depth1.tsv', 'mclean2/depth1.tsv', 'cand.tsv', self.a.depth_jump_fold,
                    p['min_depth'], 'depth_jump.tsv', log='mdepth.log', stdout='depth_jump.txt')
            self.log('  ' + open(self.path('depth_jump.txt')).read().strip())
            if sum(1 for _ in open(self.path('depth_jump.tsv'))) > 1:
                # A dropped model's copies go to whichever model matches them next, as a cut segment's do.  A third depth
                # pass without the dropped models (no search: their hits leave the second pass's unfiltered hits and
                # the mask-level filter runs again); a model it cuts differently from the second pass, which still had
                # them, takes the third pass's cuts.  (mm10: an internal's L1 piece, depth 38 beside the two dropped
                # L1 models and 7,447 without them, otherwise masked 5.2 Mb outside the curated LTRs)
                self.env['PREV_RUN'] = 'mclean2'
                self.sh(['bash', os.path.join(HERE, 'clean_library_rmb.sh'), 'models.fa', 'depth_sample.fa', 'mclean3', 1, p['scale'],
                         mtx, '--ratio', self.a.model_depth_ratio, '--min-depth', p['min_depth'], '--min-len', 50, '--consensus', 'models.cons.fa'],
                        'mdepth.log', threads=1)
                del self.env['PREV_RUN']
                second = {n.split('#')[0]: s for n, s in read_fasta(self.path('mclean2', 'round1.fa'))}
                third = {n.split('#')[0]: s for n, s in read_fasta(self.path('mclean3', 'round1.fa'))}
                recs = list(read_fasta(self.path('models.fa'))); n = bp = gone = 0
                with open(self.path('models.fa'), 'w') as fh:
                    for name, s in recs:
                        k = name.split('#')[0]
                        if second.get(k) != third.get(k):
                            t = third.get(k)
                            if t is None: gone += 1; bp += len(s); continue
                            if t != s: n += 1; bp += len(s) - len(t); s = t
                        fh.write(f'>{name}\n{s}\n')
                self.log(f'  {n} model(s) trimmed where the dropped models\' copies went, {bp:,} bases removed, {gone} dropped')

    def families(self):
        self.search('models.fa', 'models.fa', 5, 'core_ava.tsv', AVA, 'families.log', by_kind=True)
        self.py('family.py', 'models.fa', 'core_ava.tsv', 'library.fa', 'families.tsv', '--rule', 'and', '--cand', 'cand.tsv',
                '--gap-aware', log='families.log')
        # each record to its element's coding strand, before the seed alignments are built on it
        self.py('orient_models.py', 'library.fa', 'families.tsv', 'cand.tsv', log='families.log', stdout='orient.txt')
        self.log('  ' + open(self.path('orient.txt')).read().strip())

    def seeds(self):
        a = self.a
        os.makedirs(self.path('seeds'), exist_ok=True)
        for f in glob.glob(self.path('seeds', '*.stk')): os.remove(f)   # a killed run's, numbered by its library's order
        self.search('library.fa', 'depth_sample.fa', 50, 'seeds/hits.tsv', ['qseqid', 'sseqid', 'sstart', 'send', 'bitscore'], 'seeds.log')
        hits = {}
        for l in open(self.path('seeds', 'hits.tsv')):
            q, sid, ss, se, bits = l.rstrip('\n').split('\t')
            hits.setdefault(q, []).append((float(bits), sid, min(int(ss), int(se)), max(int(ss), int(se))))
        recs, origin = [], {}
        with open(self.path('library.fa')) as fh:
            for l in fh:
                if l.startswith('>'):
                    f = l[1:].split(); recs.append(f[0])
                    o = [x[7:] for x in f[1:] if x.startswith('origin=')]
                    if o: origin[f[0]] = o[0]
        # each record's best copies, highest score first, none overlapping another already taken
        pick, want = {}, {}
        for r in recs:
            took = []
            for b, sid, lo, hi in sorted(hits.get(r, []), key=lambda h: -h[0]):
                if any(sid == t[0] and lo <= t[2] and t[1] <= hi for t in took): continue
                took.append((sid, lo, hi))
                if len(took) == a.seed_copies: break
            pick[r] = took
            for sid, lo, hi in took: want.setdefault(sid, []).append((lo, hi))
        # only a sample of the genome's windows names them <sequence>_<0-based start>; the whole genome and a
        # --sample-file keep their own names (sample.txt of older runs has no kind: a symlink is the genome)
        kind = self.params().get('kind', 'genome' if os.path.islink(self.path('depth_sample.fa')) else 'windows')
        def where(sid, lo, hi):
            m = re.match(r'^(.*)_(\d+)$', sid) if kind == 'windows' else None
            return f'{m.group(1)}:{int(m.group(2)) + lo}-{int(m.group(2)) + hi}_+' if m else f'{sid}:{lo}-{hi}_+'
        # a copy with 20 bp either side, named by the range taken (Smitten, 1-based)
        PAD = 20
        seqs = {}
        for name, x in read_fasta(self.path('depth_sample.fa')):
            for lo, hi in want.get(name, ()):
                a0, b0 = max(1, lo - PAD), min(len(x), hi + PAD)
                seqs[(name, lo, hi)] = (where(name, a0, b0), x[a0 - 1:b0])
        # records without a copy in the sample: the origin element, from the genome
        lone = {r: re.match(r'^(.*):(\d+)\.\.(\d+)(_[A-Z]+)?$', origin.get(r, '')) for r in recs if not pick[r]}
        lone = {r: m for r, m in lone.items() if m}
        if lone:
            # a record turned to its coding strand (step 23) gets its origin element turned the same way
            hdr = open(self.path('families.tsv')).readline().rstrip('\n').split('\t')
            flipped = {f[hdr.index('new_name')] for f in (l.rstrip('\n').split('\t') for l in list(open(self.path('families.tsv')))[1:])
                       if 'flipped' in hdr and f[hdr.index('flipped')] == 'yes'}
            comp = str.maketrans('ACGTRYKMSWBDHVNacgtrykmswbdhvn', 'TGCAYRMKSWVHDBNtgcayrmkswvhdbn')
            need = {m.group(1) for m in lone.values()}
            for name, x in read_fasta(self.genome):
                if name not in need: continue
                for r, m in lone.items():
                    if m.group(1) == name:
                        lo, hi = int(m.group(2)), int(m.group(3))
                        if r.split('#')[0] in flipped:
                            seqs[(name, lo, hi)] = (f'{name}:{lo}-{hi}_-', x[lo - 1:hi][::-1].translate(comp))
                        else:
                            seqs[(name, lo, hi)] = (f'{name}:{lo}-{hi}_+', x[lo - 1:hi])
                        pick[r] = [(name, lo, hi)]
        fixed = bool(a.te_composer) and '--fixed-cons' in subprocess.run([a.te_composer, '--help'], capture_output=True, text=True).stdout
        if not fixed:
            self.log('  warning: te-composer has no --fixed-cons (needs >= 0.6.8): one-row seed alignments', 'warn')
        lib = dict(read_fasta(self.path('library.fa')))
        failed = []
        def one(i, r):
            ident = r.split('#')[0]; out = self.path('seeds', f'{i}.stk')
            if fixed and pick[r]:
                cons, inst = self.path('seeds', f'{i}.cons.fa'), self.path('seeds', f'{i}.inst.fa')
                with open(cons, 'w') as fh: fh.write(f'>{ident}\n{lib[r]}\n')
                with open(inst, 'w') as fh:
                    for k in pick[r]:
                        fh.write('>%s\n%s\n' % seqs[k])
                with open(self.path('logs', 'seeds.log'), 'a') as lf:
                    rc = subprocess.run([a.te_composer, '--silent', '--fixed-cons', cons, '--name', ident, '--threads', '1', inst, out],
                                        stdout=lf, stderr=lf).returncode
                os.remove(cons); os.remove(inst)
                if rc == 0 and os.path.exists(out) and os.path.getsize(out): return len(pick[r])
                if rc: failed.append(ident)
            with open(out, 'w') as fh:   # the record alone
                fh.write(f'# STOCKHOLM 1.0\n#=GF ID {ident}\n#=GF SQ 1\n{ident} {lib[r]}\n#=GC RF {lib[r]}\n//\n')
            return 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=a.threads) as ex:
            rows = list(ex.map(lambda ir: one(*ir), enumerate(recs)))
        with open(self.path('library.stk'), 'w') as fo:
            for i in range(len(recs)):
                with open(self.path('seeds', f'{i}.stk')) as fh: shutil.copyfileobj(fh, fo)
                os.remove(self.path('seeds', f'{i}.stk'))
        self.log(f'  {len(recs)} seed alignments: {sum(1 for n in rows if n)} from genome copies ({sum(rows)} rows), '
                 f'{len(lone)} from the origin element alone, {sum(1 for n in rows if not n)} one-row')
        if failed:
            self.log(f'  warning: te-composer failed on {len(failed)} record(s), left as one-row alignments ({", ".join(failed[:5])}'
                     f'{" ..." if len(failed) > 5 else ""}); see work/logs/seeds.log (--keep-work keeps it)', 'warn')
        # report only: each record's estimated share of the genome masked, from the same hits (families.tsv
        # est_genome_pct); a record at or above --mask-warn percent is named in the log, nothing is changed
        p = self.params()
        self.py('mask_share.py', os.path.join('seeds', 'hits.tsv'), 'families.tsv', p['scale'], p['genome_mb'],
                '--warn', a.mask_warn, log='seeds.log', stdout='mask_share.txt')
        ms = open(self.path('mask_share.txt')).read().strip()
        self.log('  ' + ms, 'warn' if ms.startswith('WARNING') else 'detail')

    def evidence(self):
        a = self.a
        self.sh([sys.executable, os.path.join(HERE, 'evidence_page.py'), self.out, self.genome, '-o', 'evidence.html',
                 '--domain-hmm', a.domain_hmm, '--threads', a.threads, '--work', 'evidence',
                 '--name', os.path.basename(os.path.normpath(os.path.abspath(self.top))),
                 '--rmblastn', a.rmblastn, '--nail', a.nail, '--mmseqs', a.mmseqs,
                 *(['--reference', os.path.abspath(a.evidence_reference)] if a.evidence_reference else [])],
                'evidence.log', stdout='evidence.txt')
        self.log('  ' + open(self.path('evidence.txt')).read().strip())

    def run(self):
        a = self.a
        if a.sufr_index:
            self.sufr = os.path.abspath(a.sufr_index); open(self.done_path('index'), 'w').close()
        else:
            self.sufr = self.path('genome.sufr')
        t0 = time.time()
        with open(os.path.join(self.top, 'run.log'), 'a') as fh: fh.write(f'# {time.strftime("%Y-%m-%d %H:%M:%S")}\n')
        self.log('ltr_library.py ' + ' '.join(sys.argv[1:]), 'user')
        # what each run in this outdir used, for reproducing it (run.log has one block per invocation, the newest last)
        self.log('  command: ' + ' '.join([os.path.abspath(sys.argv[0])] + sys.argv[1:]))
        for k, v in (('genome', self.genome), ('sufr index', self.sufr), ('tmpdir', self.tmp), ('domain profiles', os.path.abspath(a.domain_hmm))):
            self.log(f'  {k}: {v}')
        for k, p, flag in (('ltr-forage', a.ltr_forage, '--version'), ('rmblastn', a.rmblastn, '-version'), ('nail', a.nail, '--version'),
                           ('mmseqs', a.mmseqs, 'version'), ('te-composer', a.te_composer, '--version'), ('sufr', a.sufr, '--version')):
            if p: self.log(f'  {k}: {p} ({tool_version(p, flag)})')
        self.step('sample', '0 depth sample', self.sample)
        self.step('index', '1 index (sufr)', self.index)
        self.step('ltr-forage', '1-11 ltr-forage', self.ltr_forage)
        self.step('polya-family', '11b poly-A family veto', self.polya_family)
        if os.path.getsize(self.path('lib.fa')) == 0:
            # nothing kept: an empty library is the result, not an error in a later step
            for f in ('library.fa', 'library.stk'): open(self.path(f), 'w').close()
            open(self.path('families.tsv'), 'w').write('family\tnew_name\torigin\tkind\tlen\tstatus\tclass\n')
            self.log('Results', 'user')
            self.log('  ltr-forage kept no element: the library is empty', 'user')
            self.finish(t0)
            return
        if a.stage_b:
            self.step('precluster', '12 pre-clustering (raw library)', self.precluster)
            self.step('premodels', f'13 consensus representatives (raw clusters >= {a.min_cons})', self.premodels)
        self.step('depth', '14 depth cleanup' + (' (small clusters + representatives)' if a.stage_b else ''), self.depth)
        self.step('nested', '15 nested-LTR vote', self.nested)
        self.step('cluster', '16 element clustering + exemplars', self.cluster)
        self.step('flank', '17 flank filter (segmental duplications)', self.flank)
        self.step('contain', '18 containment', self.contain)
        self.step('excise', '19 internal nested excision', self.excise)
        self.step('pbs', '20 PBS', self.pbs)
        self.step('tier', f'21 tiering ({a.tier})', self.tier)
        self.step('models', f'22 consensus models (clusters >= {a.min_cons})', self.models)
        self.step('mdepth', '22b depth check of the models', self.mdepth)
        self.step('families', '23 80-80 families', self.families)
        self.step('seeds', '23b seed alignments', self.seeds)
        if a.evidence_page:
            self.step('evidence', '24 evidence page', self.evidence)
        kinds = [n.split('#')[0].rsplit('_', 1)[-1] for n, _ in read_fasta(self.path('library.fa'))]
        nfam = len({l.split('\t')[0] for l in open(self.path('families.tsv'))} - {'family'})
        self.log('Results', 'user')
        self.log(f'  families: {nfam} ({kinds.count("LTR")} LTR and {kinds.count("INT")} internal models)', 'user')
        if a.evidence_page:
            ev = open(self.path('evidence.txt')).read().strip()
            self.log('  evidence: ' + (ev.split('; ', 1)[1] if '; ' in ev else ev), 'user')
        self.finish(t0)

    def finish(self, t0):
        """The products to the outdir: moved, and work/ removed, or copied with --keep-work."""
        a = self.a
        done = []
        for f in PRODUCTS:
            src, dst = self.path(f), os.path.join(self.top, f)
            if not os.path.exists(src) or (f == 'evidence.html' and not a.evidence_page): continue
            if os.path.lexists(dst): os.remove(dst)
            # a copy, not a link: a later --from rerun rewrites the work file in place
            if a.keep_work: shutil.copy2(src, dst)
            else: os.replace(src, dst)
            done.append(f)
        if a.keep_work:
            try: os.rmdir(self.path('tmp'))   # the default tmpdir, when the steps left it empty
            except OSError: pass
        else:
            shutil.rmtree(self.out, ignore_errors=True)
            if os.path.exists(self.out):
                self.log(f'  warning: could not remove all of {self.out}', 'warn')
        self.log(f'  output in {a.outdir}: {", ".join(done)}, run.log'
                 + ('; intermediate files in work/' if a.keep_work else ''), 'user')
        self.log(f'done: {(time.time() - t0) / 60:.1f} min', 'user')


USER_EPILOG = """\
Output in <outdir>: library.fa (the library: one LTR and one internal model per family,
fam<N>_v<k>_LTR|INT#LTR/<class>), library.stk (a seed alignment per model), families.tsv, cand.tsv
(ltr-forage's per-candidate table), evidence.html (--evidence-page) and run.log (every step, its results,
and the programs and versions used).

Programs: ltr-forage, sufr, rmblastn (Rust port), te-composer, nail and MMseqs2, from the options above, the
environment (LTRFORAGE, SUFR, RMBLASTN, TE_COMPOSER, NAIL, MMSEQS) or PATH; Python 3 with numpy.

--help-all also lists the development and tuning options and every step."""


def main():
    # -h lists the options a user sets; --help-all also the development and tuning options, which work either way
    full = '--help-all' in sys.argv
    if full: sys.argv = ['-h' if x == '--help-all' else x for x in sys.argv]
    dev = lambda h: h if full else argparse.SUPPRESS
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split('\n', 2)[2] if full else USER_EPILOG)
    ap.add_argument('genome', help='genome assembly, uncompressed FASTA')
    ap.add_argument('-o', '--outdir', required=True, help='output directory')
    ap.add_argument('-t', '--threads', type=int, default=os.cpu_count(), help=f'threads (default: all, {os.cpu_count()} here)')
    ap.add_argument('-v', '--verbose', action='store_true', help='print every step and its results (run.log always has them)')
    ap.add_argument('--keep-work', action='store_true', help='keep the intermediate files (<outdir>/work: the steps\' files, '
                    'their logs, the sufr index the run built); a successful run otherwise removes them, and with them --from')
    ap.add_argument('--tmpdir', help='directory for temporary files of all steps (default: <outdir>/work/tmp; '
                    'the sufr index step needs up to the size of the finished index there, ~9 bytes per genome base)')
    ap.add_argument('--sufr', dest='sufr_index', metavar='INDEX', help='existing sufr index of the genome (skips step 1; default: '
                    '<genome>.sufr beside the genome when it exists, with .fa/.fasta/.fna dropped)')
    ap.add_argument('--from', dest='from_step', choices=STEPS, metavar='STEP', help='rerun from this step on (needs the intermediate '
                    'files: a run with --keep-work, or one that stopped; --help-all lists the steps)')
    ap.add_argument('--help-all', action='store_true', help='also list the development and tuning options')
    lb = ap.add_argument_group('library')
    lb.add_argument('--tier', choices=('pbs', 'strict', 'all'), default='pbs', help='single-element families kept: LTR identity >= 93 '
                    'or a PBS (pbs, default), identity >= 90 and a PBS (strict), or every one (all)')
    lb.add_argument('--mask-bed', metavar='BED', help='regions to leave out of the search (sequence name, 0-based start, end)')
    lb.add_argument('--seed-copies', type=int, default=50, metavar='N', help='genome copies per seed alignment in library.stk (50)')
    lb.add_argument('--mask-warn', type=float, default=2.0, metavar='PCT',
                    help='warn about any record estimated to mask at least PCT percent of the genome: look at it before use '
                         '(a LINE or SINE called LTR, or a very abundant real family); report only, families.tsv has every '
                         'record\'s estimate (2.0; 0: no warning)')
    e = ap.add_argument_group('evidence page')
    e.add_argument('--evidence-page', action='store_true', help='also write evidence.html: one dossier per family')
    e.add_argument('--evidence-reference', metavar='LIB', help='curated library to measure each family against on the page')
    t = ap.add_argument_group('programs and data (when not found in the environment or on PATH)')
    t.add_argument('--ltr-forage', dest='ltr_forage', metavar='PATH', help='ltr-forage binary ($LTRFORAGE)')
    t.add_argument('--sufr-bin', dest='sufr', metavar='PATH', help='sufr binary ($SUFR), for building the index')
    t.add_argument('--rmblastn', metavar='PATH', help='rmblastn, the Rust port >= 3.0.9 ($RMBLASTN)')
    t.add_argument('--te-composer', dest='te_composer', metavar='PATH', help='te-composer >= 0.6.6 ($TE_COMPOSER)')
    t.add_argument('--nail', metavar='PATH', help='nail >= 0.7.1 ($NAIL)')
    t.add_argument('--mmseqs', metavar='PATH', help='MMseqs2 ($MMSEQS)')
    t.add_argument('--domain-hmm', default=os.environ.get('HMM', os.path.join(ROOT, 'data', 'TEdomains.hmm')), metavar='HMM',
                   help="protein profiles for ltr-forage's domain search (default: data/TEdomains.hmm, 122 Pfam profiles with their "
                        "class, see data/domains/README.md; data/TEfam.hmm is LTR_retriever's file)")
    d = ap.add_argument_group('detection (development)')
    d.add_argument('--no-fastan', action='store_true',
                   help=dev('search tandem arrays for seeds too (by default arrays found by the FasTAN port are masked: '
                            'satellite pairs otherwise dominate the run and flood re-extension; gorilla candidate detection '
                            '7.3 h -> 40 min at 8 threads, with more elements kept at higher curated support)'))
    d.add_argument('--overlaps', choices=('evidence', 'best'), default='evidence',
                   help=dev('overlapping calls: keep the one with the most evidence after annotation (evidence, default), '
                            'or LTRharvest\'s best seed similarity per overlapping group before re-extension (best; faster)'))
    d.add_argument('--forage-arg', action='append', default=[], metavar='ARG',
                   help=dev('pass ARG on to ltr-forage (repeatable; write --forage-arg=--option for options)'))
    g = ap.add_argument_group('depth sample and floors (development)')
    g.add_argument('--sample-above', type=float, default=200, help=dev('sample genomes over this many non-N Mb (200)'))
    g.add_argument('--sample-frac', type=float, default=0.1, help=dev('fraction of the genome sampled above --sample-above (0.1)'))
    g.add_argument('--window-kb', type=int, default=100, help=dev('sample window (100 kb)'))
    g.add_argument('--seed', type=int, default=1, help=dev('random seed of the depth sample and the poly-A family veto (1)'))
    g.add_argument('--sample-file', help=dev('use this FASTA as the depth sample'))
    g.add_argument('--min-depth', type=int, help=dev('default 50 x genome non-N Mb / 373, capped at 50 (rice MSU7: 50)'))
    g.add_argument('--abs-depth', type=int, help=dev('default 10 x min depth'))
    m = ap.add_argument_group('family models (development)')
    m.add_argument('--min-cons', type=int, default=3, help=dev('consensus for clusters of at least this many elements (3; 0: exemplars only)'))
    m.add_argument('--stage-a', action='store_true',
                   help=dev('depth-clean and re-cluster every element (the older flow) instead of clustering the raw library '
                        'first and representing clusters >= --min-cons by their consensus (stage B, the default)'))
    m.add_argument('--greedy-above', type=float, default=100,
                   help=dev('step 12 clusters greedily when the probe projects more than this many million all-vs-all hits '
                        '(100; 0: always greedy; -1: never, and no probe)'))
    m.add_argument('--greedy-batch', type=int, default=2000, help=dev('elements per all-vs-all batch of greedy step 12 (2000)'))
    m.add_argument('--greedy-reps', type=int, default=100, help=dev('representatives per cluster in greedy step 12 (100)'))
    m.add_argument('--no-flank', action='store_true', help=dev('skip the segmental-duplication flank filter'))
    m.add_argument('--no-model-depth', action='store_true', help=dev('skip step 22b, the depth check of the finished models'))
    m.add_argument('--depth-rounds', type=int, default=3, help=dev('rounds of the step 14 depth cleanup (3; 0 skips it)'))
    m.add_argument('--no-depth-remap', action='store_true',
                   help=dev('search the cut records again in step 14 rounds after the first.  By default those rounds search '
                        'nothing: the last round\'s cuts are applied to its hits (alignments clipped, rescored, joined '
                        'across a cut), which is close to the searched rounds but not identical'))
    m.add_argument('--depth-word', type=int, default=9, help=dev('rmblastn seed word size of the step 14 and 22b depth searches (9, as RepeatMasker)'))
    m.add_argument('--polya-family', type=int, default=10, metavar='N',
                   help=dev("veto a kept element when half of each LTR is covered by N of 2,000 sampled LTRs of poly-A-vetoed candidates (10; 0: off)"))
    m.add_argument('--polya-pool', type=int, default=2000, metavar='N',
                   help=dev("run the poly-A family veto (and the poly-A repeat cut) only when ltr-forage vetoed at least N candidates for poly-A "
                        "tails, not counting those with an LTR domain hit (2000: the sample is then always full).  Not measured: genomes "
                        "with 80 to 17,000 such candidates; the step's log line gives the count"))
    m.add_argument('--no-polya-rescue', action='store_true',
                   help=dev("do not keep again the elements that lost the overlap choice only to elements the poly-A family veto removes"))
    m.add_argument('--early-tier', action='store_true', help=dev("drop raw singletons the tier would reject before step 14"))
    m.add_argument('--int-copies-max', type=float,
                   help=dev('step 22b drops an internal model with more genome copies than this per element of its cluster '
                        '(2000 with --early-tier, otherwise off; 0: off)'))
    m.add_argument('--depth-jump-fold', type=float, default=20,
                   help=dev('step 22b drops an LTR model without TSD or TG..CA whose depth rose this many times after the cuts (20; 0: off)'))
    m.add_argument('--model-depth-ratio', type=float, default=5,
                   help=dev('step 22b cuts an exemplar model where depth is at least this many times its median (5)'))
    m.add_argument('--word', type=int, default=16, help=dev('rmblastn seed word for the library all-vs-all (16; 12 gave the same '
                   'scores within 0.05 points on rice, lotus and mouse, with step 12 up to 2x slower)'))
    a = ap.parse_args()
    if not os.path.isfile(a.genome):
        sys.exit(f'ltr_library: genome not found: {a.genome}')
    with open(a.genome, 'rb') as fh:
        if fh.read(2) == b'\x1f\x8b':
            sys.exit(f'ltr_library: {a.genome} is gzip-compressed; the genome must be an uncompressed FASTA file')
    if a.early_tier and (a.stage_a or a.min_cons < 2):
        sys.exit('ltr_library: --early-tier drops raw singletons of the step 12 pre-clustering, which --stage-a and '
                 '--min-cons below 2 skip')
    if a.int_copies_max is None: a.int_copies_max = 2000 if a.early_tier else 0
    if not a.sufr_index:
        idx = re.sub(r'\.(fa|fasta|fna)$', '', a.genome) + '.sufr'
        if os.path.isfile(idx):
            if os.path.getmtime(idx) >= os.path.getmtime(a.genome): a.sufr_index = idx
            else: print(f'ltr_library: not using {idx}: it is older than the genome (--sufr {idx} to use it anyway)', file=sys.stderr)
    a.ltr_forage = find_tool(a.ltr_forage, 'LTRFORAGE', 'ltr-forage')
    a.rmblastn = find_tool(a.rmblastn, 'RMBLASTN', 'rmblastn', check='Rust port', min_version=(3, 0, 9), flag='-version')
    a.nail = find_tool(a.nail, 'NAIL', 'nail', min_version=(0, 7, 1))
    a.mmseqs = find_tool(a.mmseqs, 'MMSEQS', 'mmseqs')
    v = subprocess.run([a.rmblastn, '-version'], capture_output=True, text=True).stdout.split()
    if v[1:2] == ['3.0.10'] and 'MT_MODE' not in os.environ:
        os.environ['MT_MODE'] = '0'   # 3.0.10's streaming --mt-mode 1 can deadlock (fixed in 3.0.11)
    a.stage_b = not a.stage_a and a.min_cons >= 2   # stage B needs consensus representatives
    # consensus models (step 22) and the seed alignments (step 23b)
    a.te_composer = find_tool(a.te_composer, 'TE_COMPOSER', 'te-composer', min_version=(0, 6, 6))
    if not a.sufr_index:
        a.sufr = find_tool(a.sufr, 'SUFR', 'sufr')
    if not os.path.isfile(a.domain_hmm):
        sys.exit(f'ltr_library: domain HMM file not found: {a.domain_hmm} (give --domain-hmm)')
    a.domain_hmm = os.path.abspath(a.domain_hmm)   # the steps run in the output directory
    Run(a).run()


if __name__ == '__main__':
    main()
