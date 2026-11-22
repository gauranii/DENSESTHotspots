# -----------------------------------------------------------------------------
# File:    checks/figure_compare.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Compare two figures by the marks read by pdf_data.py.
#          - _fit: one affine map between two coordinate lists and the largest
#            leftover (pt).
#          - bars_vs_bars: histogram bars, series by series; heights must be
#            proportional.
#          - points_vs_points: point glyphs of a point-and-line plot, series
#            by series, after one affine map per axis.
#          - r_vs_python_bars, r_vs_python_points: a figure drawn by R against
#            the same figure drawn by matplotlib (python/zinull/plots.py):
#            the same marks one to one, in the same order and colours, at
#            the same places after one affine map per axis.
# -----------------------------------------------------------------------------
from collections import Counter

from pdf_data import marks

GLYPH = 15.0  # largest width or height of a point glyph (pt)


def _fit(u, v):
    n = len(u)
    mu, mv = sum(u) / n, sum(v) / n
    sxx = sum((a - mu) ** 2 for a in u)
    b = sum((a - mu) * (c - mv) for a, c in zip(u, v)) / sxx if sxx else 0.0
    a = mv - b * mu
    return a, b, max(abs(a + b * x - y) for x, y in zip(u, v))


# ----- Histograms -----

def bar_series(mk):
    """Bars of each barplot of a figure, {(fill, border): [(x, height)]} in
    drawing order. Bars are the filled rectangles on the common baseline; an
    overlaid barplot has another fill or border colour, which makes it a
    separate series."""
    rects = [r for r in mk['rects'] if r[4] is not None]
    if not rects:
        return []
    baseline = Counter(r[1] for r in rects).most_common(1)[0][0]
    series = {}
    for x, y, w, h, fill, stroke in rects:
        if y == baseline:
            series.setdefault((fill, stroke), []).append((x, h))
    return series


def counts_per_pt(mk):
    """Counts per pt of the vertical axis, from its tick labels (rotated
    numbers, centred on their ticks); None if the axis cannot be read."""
    ticks = []
    for x, y, size, string, font, rotated in mk['texts']:
        if rotated and string.isdigit():
            ticks.append((y + 0.278 * size * len(string), float(string)))  # Helvetica digit: 0.556 em
    if len(ticks) < 2:
        return None
    return abs(_fit([t[0] for t in ticks], [t[1] for t in ticks])[1])


def bars_vs_bars(pdf_a, pdf_b):
    """Bars of figure b against figure a (the published one). Returns
    (bars per series in a, bars per series in b, largest leftover of the bar
    positions after one affine map (pt), largest leftover of the bar heights
    after one scale, unit of that leftover: 'counts' or 'pt'). The leftovers
    are None when the numbers of bars differ."""
    ma = marks(pdf_a)
    a, b = bar_series(ma), bar_series(marks(pdf_b))
    if set(a) == set(b):  # same colours: pair the series by colour, else by drawing order
        b = {colour: b[colour] for colour in a}
    na, nb = [len(s) for s in a.values()], [len(s) for s in b.values()]
    if na != nb or not na:
        return na, nb, None, None, None
    a, b = [bar for s in a.values() for bar in s], [bar for s in b.values() for bar in s]
    position = _fit([bar[0] for bar in b], [bar[0] for bar in a])[2]
    ha, hb = [bar[1] for bar in a], [bar[1] for bar in b]
    sbb = sum(h * h for h in hb)
    scale = sum(x * y for x, y in zip(ha, hb)) / sbb if sbb else 0.0
    height = max(abs(x - scale * y) for x, y in zip(ha, hb))
    per_pt = counts_per_pt(ma)
    if per_pt is None:
        return na, nb, position, height, 'pt'
    return na, nb, position, height * per_pt, 'counts'


# ----- Point-and-line plots -----

def _glyphs(mk):
    """Point glyphs [(kind, x, y)] in drawing order. Circles are the
    ZapfDingbats character 'l' (or four Bezier curves), triangles and squares
    are small closed paths, stars and crosses are short strokes that share
    one midpoint."""
    found = []  # (order, kind, x, y)
    for i, (x, y, size, string, font, rotated) in enumerate(mk['texts']):
        if string == 'l' and font == '/F1':
            found.append((('t', i), 'circle', x + 0.396 * size, y + 0.347 * size))
    strokes = []
    for i, p in enumerate(mk['paths']):
        xs, ys = [q[0] for q in p['points']], [q[1] for q in p['points']]
        if max(xs) - min(xs) > GLYPH or max(ys) - min(ys) > GLYPH:
            continue
        centre = (max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2
        if p['curves'] == 4:
            found.append((('p', i), 'circle') + centre)
        elif p['closed'] and len(xs) == 3:
            found.append((('p', i), 'triangle', sum(xs) / 3, sum(ys) / 3))
        elif p['closed'] and len(xs) == 4:
            found.append((('p', i), 'square') + centre)
        elif len(xs) == 2 and not p['closed']:
            if strokes and max(abs(strokes[-1][1] - centre[0]), abs(strokes[-1][2] - centre[1])) < 0.02:
                strokes[-1][3] += 1
            else:
                strokes.append([('p', i)] + list(centre) + [1])
    found += [(order, 'star' if n > 2 else 'cross', x, y) for order, x, y, n in strokes if n > 1]
    return found


def point_series(pdf):
    """Points of each series of a point-and-line plot, {name: [(x, y)]} with
    the points ordered by x. A series is one glyph shape; its name is the
    legend label next to that shape, or the shape when there is no legend."""
    mk = marks(pdf)
    glyphs = _glyphs(mk)
    if not glyphs:
        return {}
    # the panel: the smallest rectangle that holds most glyphs (the rest are legend keys)
    def inside(r, g):
        return r[0] <= g[2] <= r[0] + r[2] and r[1] <= g[3] <= r[1] + r[3]
    boxes = [r[:4] for r in mk['rects']] + list(mk['clips'])
    boxes = [r for r in boxes if sum(inside(r, g) for g in glyphs) >= 0.6 * len(glyphs)]
    panel = min(boxes, key=lambda r: r[2] * r[3]) if boxes else None
    names = {}
    for order, kind, x, y in glyphs:
        if panel and not inside(panel, (order, kind, x, y)):
            labels = [t for t in mk['texts'] if t[0] > x and abs(t[1] + 0.35 * t[2] - y) < 6 and not t[5]]
            if labels:
                names.setdefault(kind, min(labels, key=lambda t: t[0])[3])
    series = {}
    for order, kind, x, y in glyphs:
        if panel is None or inside(panel, (order, kind, x, y)):
            series.setdefault(names.get(kind, kind), []).append((x, y))
    return {name: sorted(points) for name, points in series.items()}


def points_vs_points(pdf_a, pdf_b):
    """Points of figure b against figure a (the published one). Returns
    ([(series, points in a, points in b)], largest leftover (pt) after one
    affine map per axis from b to a, or None when the numbers of points
    differ). Series are paired by name, or by position when the names differ."""
    a, b = point_series(pdf_a), point_series(pdf_b)
    if set(a) == set(b):
        pairs = [(name, a[name], b[name]) for name in a]
    else:
        ka, kb = list(a), list(b)
        pairs = [(ka[i] if i < len(ka) else kb[i], a[ka[i]] if i < len(ka) else [],
                  b[kb[i]] if i < len(kb) else []) for i in range(max(len(ka), len(kb)))]
    counts = [(name, len(pa), len(pb)) for name, pa, pb in pairs]
    if not pairs or any(na != nb for _, na, nb in counts):
        return counts, None
    pa = [p for _, points, _ in pairs for p in points]
    pb = [p for _, _, points in pairs for p in points]
    leftover = max(_fit([p[0] for p in pb], [p[0] for p in pa])[2],
                   _fit([p[1] for p in pb], [p[1] for p in pa])[2])
    return counts, leftover


# ----- R against matplotlib -----

POSITION_TOLERANCE = 0.02  # pt of the R page; R prints its coordinates with two decimals
COLOUR_TOLERANCE = 0.002   # of a colour channel in [0, 1]; R prints three decimals
GLYPH_TO_VERTEX = 0.1      # pt between a point glyph of R and the line vertex under it


def _rgb(colour):
    return tuple(float(v) for v in colour.split())


def _same_colour(a, b):
    return max(abs(x - y) for x, y in zip(_rgb(a), _rgb(b))) <= COLOUR_TOLERANCE


def r_bars(pdf):
    """Bars of an R barplot in drawing order: [(left, width, height, fill,
    border)]; the filled rectangles on the common baseline."""
    rects = [r for r in marks(pdf)['rects'] if r[4] is not None]
    if not rects:
        return []
    baseline = Counter(r[1] for r in rects).most_common(1)[0][0]
    return [(x, w, h, fill, stroke) for x, y, w, h, fill, stroke in rects if y == baseline]


def python_bars(pdf):
    """Bars of a matplotlib bar plot in drawing order, as `r_bars()`: the
    closed four-cornered paths that are filled and stroked."""
    bars = []
    for p in marks(pdf)['paths']:
        if p['closed'] and len(p['points']) == 4 and p['fill'] and p['stroke'] and not p['curves']:
            xs, ys = [q[0] for q in p['points']], [q[1] for q in p['points']]
            bars.append((min(xs), max(xs) - min(xs), max(ys) - min(ys), p['fill'], p['stroke']))
    return bars


def r_vs_python_bars(r_pdf, py_pdf):
    """Bars of a histogram drawn by R and by matplotlib. Returns (bars in
    the R file, bars in the Python file, numbers compared, largest leftover
    in pt of the R page). The leftover is None when the bars do not
    correspond: another number, or a bar with another fill or border. Left
    and right ends are compared after one affine map, heights after one
    scale (heights must be proportional)."""
    r, py = r_bars(r_pdf), python_bars(py_pdf)
    if len(r) != len(py) or not r:
        return len(r), len(py), 0, None
    if not all(_same_colour(a[3], b[3]) and _same_colour(a[4], b[4]) for a, b in zip(r, py)):
        return len(r), len(py), 0, None
    ends = _fit([b[0] for b in py] + [b[0] + b[1] for b in py],
                [a[0] for a in r] + [a[0] + a[1] for a in r])[2]
    ha, hb = [a[2] for a in r], [b[2] for b in py]
    sbb = sum(h * h for h in hb)
    scale = sum(x * y for x, y in zip(ha, hb)) / sbb if sbb else 0.0
    heights = max(abs(x - scale * y) for x, y in zip(ha, hb))
    return len(r), len(py), 3 * len(r), max(ends, heights)


def _polylines(mk):
    """Stroked open paths of straight segments, not filled."""
    return [p['points'] for p in mk['paths']
            if p['stroke'] and not p['fill'] and not p['closed'] and not p['curves']]


def _distance(a, b):
    return max(max(abs(p[0] - q[0]), abs(p[1] - q[1])) for p, q in zip(a, b))


def r_point_lines(pdf):
    """Series of an R (ggplot2) point-and-line plot, in the order of their
    names: [(name, glyph centres ordered by x, vertices of the line through
    them)]; the line is None when no stroked path runs through the points."""
    lines = _polylines(marks(pdf))
    series = []
    for name, points in sorted(point_series(pdf).items()):
        through = [line for line in lines
                   if len(line) == len(points) and _distance(line, points) < GLYPH_TO_VERTEX]
        series.append((name, points, through[0] if through else None))
    return series


def python_point_lines(pdf):
    """Series of a matplotlib point-and-line plot in drawing order: [(marker
    name, marker centres, vertices of the line through them)]. A series is a
    marker form stamped at every vertex of one stroked line; the stamp in the
    legend lies on no vertex and is left out."""
    mk = marks(pdf)
    lines = _polylines(mk)
    stamps = {}
    for name, x, y in mk['stamps']:
        stamps.setdefault(name, []).append((x, y))
    series = []
    for name, points in stamps.items():
        through = [line for line in lines
                   if len(line) <= len(points) and _distance(line, points) < 1e-3]
        line = max(through, key=len) if through else None
        series.append((name, points[:len(line)] if line else points, line))
    return series


def r_vs_python_points(r_pdf, py_pdf):
    """Points and lines of a null-fit plot drawn by R and by matplotlib.
    Returns ([(series, points in R, points in Python)], numbers compared,
    largest leftover in pt of the R page after one affine map per axis). The
    leftover is None when the marks do not correspond: another number of
    series, of points in a series, or a series without its line. The series
    of R are taken in the order of their names, as matplotlib draws them."""
    r, py = r_point_lines(r_pdf), python_point_lines(py_pdf)
    counts = [(r[i][0] if i < len(r) else py[i][0], len(r[i][1]) if i < len(r) else 0,
               len(py[i][1]) if i < len(py) else 0) for i in range(max(len(r), len(py)))]
    if not r or len(r) != len(py) or any(a != b for _, a, b in counts):
        return counts, 0, None
    if any(line is None for _, _, line in r + py):
        return counts, 0, None
    leftover = 0.0
    for part in (1, 2):  # glyph centres, then line vertices
        a = [p for s in r for p in s[part]]
        b = [p for s in py for p in s[part]]
        for axis in (0, 1):
            leftover = max(leftover, _fit([p[axis] for p in b], [p[axis] for p in a])[2])
    return counts, 4 * sum(n for _, n, _ in counts), leftover
