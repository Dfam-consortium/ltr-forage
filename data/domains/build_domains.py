#!/usr/bin/env python3
"""Build ltr-forage's annotated protein-domain profile file from Pfam.

usage: build_domains.py <Pfam-A.hmm[.gz]> [-p profiles.tsv] [-o ../TEdomains.hmm]

Takes the profiles listed in profiles.tsv (accession, class, role) out of a Pfam-A.hmm release and writes
them with their class and role in the profile file itself, which is how ltr-forage reads them
(src/domain.rs, ProfileTable):

    NAME  rve#LTR                      <Pfam name>#<RepeatMasker class>
    ACC   PF00665.33
    DESC  Integrase core domain [role=core]

The class is a RepeatMasker type or type/subtype.  A type of LTR is evidence for an LTR candidate, any
other type (DNA, LINE, RC, Other ...) evidence against it, and an LTR subtype of Gypsy or Copia counts
toward the superfamily call.  Role 'core' marks the LTR polyprotein domains (the coding-LTR veto and the
nesting test use them); 'accessory' LTR profiles and the not-LTR profiles carry no tag.  HMMER and nail
read the file as usual.

An optional sixth column, source, takes a profile from another profile file instead of the Pfam release
given on the command line: a path relative to ltr-forage's data/ directory (e.g. TEfam.hmm, for a 2015
build a later release made worse) or an absolute path.  Empty or 'pfam' means the release.

Profiles keep the version of the release they come from; an accession missing from the release is an
error unless --allow-missing is given (Pfam retires and merges families between releases).  Any HMMER3
profile file with ACC lines can stand in for Pfam-A.hmm, e.g. an older TEfam.hmm.
"""
import argparse, gzip, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROLES = {'core', 'accessory', '-'}


def read_table(path):
    rows = []
    for n, line in enumerate(open(path), 1):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 3:
            sys.exit(f'{path}:{n}: need accession, class and role')
        acc, cls, role = f[0], f[1], f[2]
        src = f[5].strip() if len(f) > 5 and f[5].strip() not in ('', 'pfam') else None
        ty = cls.split('/')[0]
        if not ty or ' ' in cls or '#' in cls:
            sys.exit(f'{path}:{n}: bad class {cls!r}')
        if role not in ROLES:
            sys.exit(f'{path}:{n}: role must be core, accessory or -, not {role!r}')
        if (ty == 'LTR') != (role in ('core', 'accessory')):
            sys.exit(f'{path}:{n}: LTR profiles are core or accessory, other classes -')
        if any(r[0] == acc for r in rows):
            sys.exit(f'{path}:{n}: {acc} listed twice')
        rows.append((acc, cls, role, src))
    return rows


def records(path):
    """HMMER3 profile records as lists of lines."""
    op = gzip.open if path.endswith('.gz') else open
    buf = []
    with op(path, 'rt') as fh:
        for line in fh:
            if not buf and not line.strip():
                continue            # TEfam.hmm starts with a blank line
            buf.append(line if line.endswith('\n') else line + '\n')   # TEfam.hmm's last line has none
            if line.startswith('//'):
                yield buf
                buf = []


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split('\n', 2)[2])
    ap.add_argument('pfam', help='Pfam-A.hmm or Pfam-A.hmm.gz (or another HMMER3 profile file)')
    ap.add_argument('-p', '--profiles', default=os.path.join(HERE, 'profiles.tsv'))
    ap.add_argument('-o', '--out', default=os.path.join(HERE, '..', 'TEdomains.hmm'))
    ap.add_argument('--allow-missing', action='store_true', help='leave out listed accessions the release lacks')
    a = ap.parse_args()

    table = read_table(a.profiles)
    sources = {}
    for acc, cls, role, src in table:
        path = a.pfam if src is None else os.path.join(HERE, '..', src) if not os.path.isabs(src) else src
        sources.setdefault(path, {})[acc] = (cls, role)
    found, origin = {}, {}
    for path, want in sources.items():
        for rec in records(path):
            acc = next((l.split()[1] for l in rec if l.startswith('ACC ')), None)
            if acc is None or acc.split('.')[0] not in want:
                continue
            cls, role = want[acc.split('.')[0]]
            out = []
            for l in rec:
                if l.startswith('NAME '):
                    name = l.split()[1]
                    if '#' in name:
                        name = name.split('#')[0]   # rebuilding from an annotated file
                    l = f'NAME  {name}#{cls}\n'
                elif l.startswith('DESC '):
                    desc = l[5:].strip().replace('[role=core]', '').strip()
                    l = f'DESC  {desc}' + (' [role=core]' if role == 'core' else '') + '\n'
                out.append(l)
            if not any(l.startswith('DESC ') for l in out):
                if role == 'core':
                    i = next(i for i, l in enumerate(out) if l.startswith('ACC '))
                    out.insert(i + 1, 'DESC  [role=core]\n')
            found[acc.split('.')[0]] = ''.join(out); origin[acc.split('.')[0]] = path

    missing = [acc for acc, *_ in table if acc not in found]
    if missing and not a.allow_missing:
        sys.exit(f'build_domains: not found in their source: {", ".join(missing)} (--allow-missing leaves them out)')
    with open(a.out, 'w') as fo:
        for acc, *_ in table:
            if acc in found:
                fo.write(found[acc])
    n_core = sum(1 for acc, _, r, _ in table if r == 'core' and acc in found)
    n_ltr = sum(1 for acc, c, _, _ in table if c.split('/')[0] == 'LTR' and acc in found)
    other = sorted(acc for acc in found if origin[acc] != a.pfam)
    if other:
        print(f'build_domains: from other files: {", ".join(other)}', file=sys.stderr)
    print(f'build_domains: {len(found)} profiles to {a.out} ({n_ltr} LTR, {n_core} of them core; '
          f'{len(found) - n_ltr} not LTR)' + (f'; missing: {", ".join(missing)}' if missing else ''), file=sys.stderr)


if __name__ == '__main__':
    main()
