# -----------------------------------------------------------------------------
# File:    checks/check_figures.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Compare every figure drawn by the scripts with the published file
#          by the marks in the two PDF files.
#          - Histograms (Figures 1 and 2, Web Figures S1 to S6): the bars of
#            each barplot; same number of bars in each series, same positions
#            and heights proportional (one scale per panel, largest leftover
#            below 0.05 counts).
#          - Point-and-line plots (Web Figures S7 to S12): the point glyphs of
#            each series; same number of points and largest leftover below
#            0.5 pt after one affine map per axis.
# Usage:   python3 checks/check_figures.py <published dir> <figures dir>
#          <checks output dir>
# Outputs: <checks output dir>/figures.csv; prints one line per figure (a
#          figure that was not drawn is skipped); exit status 1 if a figure
#          differs.
# -----------------------------------------------------------------------------
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_compare import bars_vs_bars, points_vs_points  # noqa: E402

published, drawn, out = sys.argv[1:4]
HEIGHT, POSITION, POINT = 0.05, 0.5, 0.5  # largest leftovers allowed: counts, pt, pt
rows = []


def label(name):
    """Fig1A_... -> 'Figure 1A', FigS10B_... -> 'Web Figure S10B'."""
    panel = name.split('_')[0]
    return panel.replace('FigS', 'Web Figure S') if panel.startswith('FigS') else panel.replace('Fig', 'Figure ')


def order(name):
    web, number, panel = re.match(r'Fig(S)?(\d+)(\w)', name).groups()
    return bool(web), int(number), panel


def report(name, kind, marks_a, marks_b, leftover, result):
    rows.append({'figure': label(name), 'file': name, 'kind': kind, 'published': marks_a,
                 'computed': marks_b, 'leftover': leftover, 'result': result})
    print(f"{label(name):16} {kind:10} {result:9} {marks_a} / {marks_b}  {leftover}")


# A drawn figure has the name of its published file.
for name in sorted((n for n in os.listdir(published) if n.endswith('.pdf')), key=order):
    a, b = os.path.join(published, name), os.path.join(drawn, name)
    kind = 'histogram' if '_histogram_' in name else 'points'
    if not os.path.exists(b):
        report(name, kind, '', '', '', 'not drawn')
    elif kind == 'histogram':
        na, nb, position, height, unit = bars_vs_bars(a, b)
        same = na == nb and position is not None and position < POSITION and height < HEIGHT
        leftover = '' if position is None else f'heights {height:.3f} {unit}, positions {position:.3f} pt'
        report(name, kind, ' + '.join(map(str, na)) + ' bars', ' + '.join(map(str, nb)) + ' bars',
               leftover, 'matches' if same else 'differs')
    else:
        counts, leftover = points_vs_points(a, b)
        same = leftover is not None and leftover < POINT
        report(name, kind, ', '.join(f'{s} {n}' for s, n, _ in counts) + ' points',
               ', '.join(f'{s} {n}' for s, _, n in counts) + ' points',
               '' if leftover is None else f'{leftover:.3f} pt', 'matches' if same else 'differs')

os.makedirs(out, exist_ok=True)
with open(os.path.join(out, 'figures.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
for kind in ('histogram', 'points'):
    of = [r for r in rows if r['kind'] == kind]
    print(f"{kind}: {sum(r['result'] == 'matches' for r in of)} of {len(of)} match, "
          f"{sum(r['result'] == 'differs' for r in of)} differ, "
          f"{sum(r['result'] == 'not drawn' for r in of)} not drawn")
sys.exit(1 if any(r['result'] == 'differs' for r in rows) else 0)
