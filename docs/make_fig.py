# Builds the pipeline figure: python3 docs/make_fig.py docs/pipeline_figure.svg
# Step badges are the numbers ltr-forage (1-11) and pipeline/ltr_library.py (11b, 12-24, 22b, 23b) print in their logs.
import html, sys
out = []
E = html.escape
def t(x, y, s, size=12, weight=400, anchor='start', style=''):
    out.append(f'<text x="{x}" y="{y}" font-size="{size}" font-weight="{weight}" text-anchor="{anchor}"{style}>{E(s)}</text>')
def box(x, y, w, h, title, lines=(), filt=False, rx=4):
    out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{"#e4e4e4" if filt else "#ffffff"}" stroke="#222" stroke-width="1.1"/>')
    cx = x + w / 2
    if lines:
        t(cx, y + 19, title, 13.5, 600, 'middle')
        for i, l in enumerate(lines): t(cx, y + 36 + 14 * i, l, 11.5, 400, 'middle', ' fill="#333"')
    else:
        t(cx, y + h / 2 + 5, title, 13.5, 600, 'middle')
def pill(x, y, w, h, s):
    out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{h/2}" fill="#222" stroke="#222"/>')
    t(x + w / 2, y + h / 2 + 4.5, s, 12.5, 600, 'middle', ' fill="#fff"')
def arrow(pts, label=None, lx=0, ly=0, anchor='start'):
    d = 'M' + ' L'.join(f'{a},{b}' for a, b in pts)
    out.append(f'<path d="{d}" fill="none" stroke="#222" stroke-width="1.2" marker-end="url(#ah)"/>')
    if label: t(lx, ly, label, 11, 400, anchor, ' font-style="italic" fill="#444"')
def badge(x, y, s):   # step number, as ltr-forage and the driver print it in their logs
    r = 11 if len(s) < 3 else 13
    out.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="#222"/>')
    t(x, y + 4, s, 10.5 if len(s) < 3 else 9.5, 700, 'middle', ' fill="#fff"')
def hexa(x, y, w, h, lines):
    p = 12
    pts = [(x, y + h / 2), (x + p, y), (x + w - p, y), (x + w, y + h / 2), (x + w - p, y + h), (x + p, y + h)]
    out.append('<polygon points="' + ' '.join(f'{a},{b}' for a, b in pts) + '" fill="#ffffff" stroke="#222" stroke-width="1.1"/>')
    y0 = y + h / 2 - 7 * (len(lines) - 1) + 4
    for i, l in enumerate(lines): t(x + w / 2, y0 + 14 * i, l, 11, 400, 'middle')

W = 1000
# ---- panel a ----
t(30, 34, 'a', 18, 700); t(52, 34, 'Candidate detection (ltr-forage)', 15, 600)
AX, AW = 50, 400; ac = AX + AW / 2
pill(AX + 70, 52, AW - 140, 28, 'Genome assembly')
A = [  # y, h, title, lines, filter
    (104, 44, 'Suffix-array index', [], False),
    (172, 64, 'Tandem-array mask', ['FasTAN: 8-mer spacing + banded wave alignment;', 'arrays with units ≤ 30 kb hidden from seeding'], True),
    (260, 64, 'Seed repeat pairs', ['maximal repeats; LTR 100–1,000 bp,', 'LTR-start spacing 1–15 kb'], True),
    (348, 64, 'X-drop extension', ['gapped X-drop from the seed; length and spacing', 'limits; seed alignment identity ≥ 85%'], True),
    (436, 50, 'Duplicate removal', ['identical outer boundaries'], True),
    (510, 64, 'LTR re-extension', ['local alignment of the LTR pair sets boundaries and identity;', '2.5 kb window, 7 kb when an end nears the window edge'], False),
    (598, 50, 'Near-duplicate merge', ['final ends within 20 bp of a better call'], True),
]
prev = 80
for i, (y, h, ti, ls, f) in enumerate(A):
    arrow([(ac, prev), (ac, y)]); box(AX, y, AW, h, ti, ls, f); badge(AX - 17, y + h / 2, str(i + 1)); prev = y + h
# bracket: the steps that are LTRharvest's algorithm
bx, y0, y1 = AX + AW + 10, 260, 486
out.append(f'<path d="M{bx - 5},{y0} L{bx + 2},{y0} L{bx + 2},{y1} L{bx - 5},{y1}" fill="none" stroke="#222" stroke-width="1.1"/>')
out.append(f'<text x="{bx + 10}" y="{(y0 + y1) / 2}" font-size="11" text-anchor="middle" transform="rotate(90 {bx + 10} {(y0 + y1) / 2})" fill="#333">LTRharvest (Ellinghaus et al. 2008)</text>')
# evidence annotation box with four sub-boxes
ey, eh = 672, 104
arrow([(ac, prev), (ac, ey)])
out.append(f'<rect x="{AX}" y="{ey}" width="{AW}" height="{eh}" rx="4" fill="#ffffff" stroke="#222" stroke-width="1.1"/>')
t(ac, ey + 19, 'Evidence annotation (final boundaries)', 13.5, 600, 'middle'); badge(AX - 17, ey + eh / 2, '8')
sw, gap = 92, 6; sx0 = AX + (AW - 4 * sw - 3 * gap) / 2
for i, (n, a, b) in enumerate([('8a', 'LTR self-', 'alignment'), ('8b', 'Protein', 'domains'), ('8c', 'tRNA genes', '→ PBS'), ('8d', 'TSD and', 'TG..CA')]):
    sx = sx0 + i * (sw + gap)
    out.append(f'<rect x="{sx}" y="{ey + 32}" width="{sw}" height="{eh - 44}" rx="3" fill="#f4f4f4" stroke="#666" stroke-width="0.9"/>')
    t(sx + sw / 2, ey + 48, n, 10.5, 700, 'middle', ' fill="#444"')
    t(sx + sw / 2, ey + 64, a, 11.5, 400, 'middle'); t(sx + sw / 2, ey + 79, b, 11.5, 400, 'middle')
prev = ey + eh
dy = 800
arrow([(ac, prev), (ac, dy)]); box(AX, dy, AW, 50, 'Evidence decision', ['ordered rules (panel c)'], True); badge(AX - 17, dy + 25, '9')
oy = 874
arrow([(ac, dy + 50), (ac, oy)])
box(AX, oy, AW, 64, 'Overlap resolution', ['among overlapping kept calls, keep the one', 'with the most evidence'], True); badge(AX - 17, oy + 32, '10')
arrow([(ac, oy + 64), (ac, 962)])
pill(AX, 962, AW, 28, 'Kept elements: 5′ LTR · internal · 3′ LTR · strand'); badge(AX - 17, 976, '11')
AEND = 990

# ---- panel b ----
t(510, 34, 'b', 18, 700); t(532, 34, 'Family library construction', 15, 600)
BX, BW = 550, 400; bc = BX + BW / 2
pill(BX + 70, 52, BW - 140, 28, 'Element library')
# connector a -> b
arrow([(AX + AW, 976), (495, 976), (495, 66), (BX + 70, 66)])
B = [   # height, title, lines, filter, badge
    (64, 'Poly-A family veto', ['drop elements whose LTRs are copies of the poly-A-tailed', 'repeats ltr-forage vetoed (genomes with ≥ 2,000 of them)'], True),
    (92, 'Pre-clustering', ['LTR-vs-LTR and internal-vs-internal alignments; single-', 'linkage at ≥ 80% identity over 80% of both internals;', 'above 100 M projected hits, greedy: batches all-vs-all, the', 'rest against at most 100 members per cluster'], False),
    (64, 'Consensus representatives', ['clusters of ≥ 3 elements stand in as their consensus', '(at most 150 copies, a seeded random sample)'], False),
    (64, 'Depth cleanup (×3)', ['excise segments far above their record’s copy depth; one', 'search, later rounds move its alignments onto the cut records'], True),
    (50, 'Nested-LTR vote', ['excise LTRs absent from most cluster members'], True),
    (64, 'Element clustering', ['same single-linkage rule; singletons attached to the', 'best cluster; exemplar = member with longest internal'], False),
    (50, 'Flank filter', ['drop clusters whose copies share both flanks'], True),
    (50, 'Containment', ['drop exemplars ≥ 90% contained in another'], True),
    (50, 'Internal nested excision (×2)', ['remove nested LTRs from exemplar internals'], True),
    (50, 'Primer binding sites', ['PBS of each exemplar from the tRNA calls'], False),
    (50, 'Tiering', ['multi-copy clusters; singletons need identity ≥ 93% or a PBS'], True),
    (64, 'Consensus models', ['clusters of ≥ 3 (≤ 150 copies); LTR length outliers', '(> 2 × median) and sparsely occupied columns removed'], True),
    (92, 'Model depth check', ['cut segments ≥ 5× an exemplar’s own copy depth, ≥ 10× a', 'consensus model’s; drop LTRs whose depth rose ≥ 20× after', 'the cuts with neither a TSD nor TG..CA, then cut where', 'the dropped models’ copies went'], True),
    (64, '80–80 families', ['LTR and internal models grouped separately; each model', 'turned to its element’s coding strand (domains, else PBS)'], False),
    (50, 'Seed alignments', ['≤ 50 genome copies aligned to each fixed model'], False),
]
BADGES = ['11b'] + [str(n) for n in range(12, 23)] + ['22b', '23', '23b']
y = 104; prev = 80
for i, (h, ti, ls, f) in enumerate(B):
    arrow([(bc, prev), (bc, y)]); box(BX, y, BW, h, ti, ls, f); badge(BX - 17, y + h / 2, BADGES[i]); prev = y + h; y += h + 12
arrow([(bc, prev), (bc, y)])
pill(BX + 10, y, BW - 20, 28, 'Library: LTR + internal model per family (FASTA, Stockholm)')
# optional evidence page beside the library
ey2 = y + 28 + 14
out.append(f'<rect x="{BX + 60}" y="{ey2}" width="{BW - 120}" height="30" rx="4" fill="#ffffff" stroke="#222" stroke-width="1.1" stroke-dasharray="4 3"/>')
t(bc, ey2 + 19.5, 'Evidence page (optional): one dossier per family', 11.5, 400, 'middle')
badge(BX + 43, ey2 + 15, '24')
yb = ey2 + 30

# ---- panel c ----
cy = max(yb, AEND) + 50
t(30, cy, 'c', 18, 700); t(52, cy, 'Evidence decision (applied in order; the first rule that matches rejects)', 15, 600)
ny, nh, nw, step = cy + 26, 70, 124, 170
rules = [['Domains of a', 'non-LTR order?'], ['LTR ≥ 99%', 'self-repeat?'], ['Coding domain', 'inside an LTR?'],
         ['Both LTRs', 'poly-A tailed,', 'no TG..CA?'], ['No domain, and', 'no TG..CA with', 'TSD ≥ 4?']]
ry = ny + nh + 40
for i, r in enumerate(rules):
    x = 40 + i * step
    hexa(x, ny, nw, nh, r)
    nxt = x + step
    arrow([(x + nw, ny + nh / 2), (nxt - 2, ny + nh / 2)], 'no', x + nw + (step - nw) / 2, ny + nh / 2 - 6, 'middle')
    arrow([(x + nw / 2, ny + nh), (x + nw / 2, ry)], 'yes', x + nw / 2 + 5, ny + nh + 24)
kx = 40 + len(rules) * step
out.append(f'<rect x="{kx}" y="{ny + 12}" width="{960 - kx}" height="{nh - 24}" rx="4" fill="#222"/>')
t((kx + 960) / 2, ny + nh / 2 + 5, 'Keep', 13.5, 600, 'middle', ' fill="#fff"')
rw = (len(rules) - 1) * step + nw
out.append(f'<rect x="40" y="{ry}" width="{rw}" height="30" rx="4" fill="#e4e4e4" stroke="#222" stroke-width="1.1"/>')
t(40 + rw / 2, ry + 19.5, 'Reject (reason recorded in the candidate table)', 12.5, 600, 'middle')
# legend
ly = ry + 56
out.append(f'<rect x="40" y="{ly - 11}" width="22" height="14" rx="2" fill="#fff" stroke="#222"/>'); t(70, ly, 'builds or annotates', 11.5)
out.append(f'<rect x="210" y="{ly - 11}" width="22" height="14" rx="2" fill="#e4e4e4" stroke="#222"/>'); t(240, ly, 'removes candidates, elements or sequence', 11.5)
H = ly + 16

svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
       'font-family="Helvetica Neue, Helvetica, Arial, sans-serif" fill="#111">'
       '<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
       '<path d="M0,0 L10,5 L0,10 z" fill="#222"/></marker></defs>'
       f'<rect width="{W}" height="{H}" fill="#fff"/>' + ''.join(out) + '</svg>')
open(sys.argv[1], 'w').write(svg)
