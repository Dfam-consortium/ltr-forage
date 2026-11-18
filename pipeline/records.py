"""Library record names, shared by the pipeline helpers: a record is <seq>:<start>..<end>_LTR or _INT
(a '#class' may follow), an element <seq>:<start>..<end>."""
import bisect, collections, csv, os, re

NAME = re.compile(r'(.+):(\d+)\.\.(\d+)(?:_(LTR|INT))?$')


def base(name):
    """The record name without its '#class'."""
    return name.split('#')[0]


def is_int(name):
    return base(name).endswith('_INT')


def is_ltr(name):
    return base(name).endswith('_LTR')


def element_lookup(elements, tsv=None):
    """f(record name) -> the element of `elements` the record belongs to, or None.

    With the ltr-forage table (tsv, default the TSV environment variable): by the exact names of each kept
    element's records, the 5' LTR, the internal region and the 3' LTR.  A record not named there, or
    without a table: the smallest element that contains it, so the records of an element nested in another
    are its own and not the outer element's (the first element containing the record was taken before)."""
    elements = set(elements)
    exact = {}
    tsv = tsv or os.environ.get('TSV')
    if tsv and os.path.exists(tsv):
        for r in csv.DictReader(open(tsv), delimiter='\t'):
            if r['keep'] != '1': continue
            s = r['seq']; e = f"{s}:{r['elem_start']}..{r['elem_end']}"
            if e not in elements: continue
            ls, le, rs, re_ = int(r['lltr_start']), int(r['lltr_end']), int(r['rltr_start']), int(r['rltr_end'])
            for n in (f'{s}:{ls}..{le}_LTR', f'{s}:{le + 1}..{rs - 1}_INT', f'{s}:{rs}..{re_}_LTR'):
                exact.setdefault(n, e)
    by = collections.defaultdict(list); longest = collections.Counter()
    for e in elements:
        m = NAME.match(e); c, s, t = m[1], int(m[2]), int(m[3])
        by[c].append((s, t, e)); longest[c] = max(longest[c], t - s)
    for v in by.values(): v.sort()

    def lookup(name):
        n = base(name)
        if n in exact: return exact[n]
        m = NAME.match(n)
        if not m: return None
        c, s, t = m[1], int(m[2]), int(m[3]); v = by.get(c, ())
        i = bisect.bisect_right(v, (s, float('inf'), '')) - 1; best = None
        while i >= 0 and v[i][0] >= s - longest[c]:
            if v[i][1] >= t and (best is None or v[i][1] - v[i][0] < best[1] - best[0]): best = v[i]
            i -= 1
        return best[2] if best else None
    return lookup
