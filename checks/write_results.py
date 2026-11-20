# -----------------------------------------------------------------------------
# File:    checks/write_results.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Turn the outputs of check_tables.py and check_figures.py into
#          checks/RESULTS.md.
# Usage:   python3 checks/write_results.py <checks output dir> > checks/RESULTS.md
# -----------------------------------------------------------------------------
import csv
import os
import sys

out = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))
tables = list(csv.DictReader(open(os.path.join(out, 'tables.csv'))))
figures = list(csv.DictReader(open(os.path.join(out, 'figures.csv'))))
paper = list(csv.DictReader(open(os.path.join(here, 'paper_values.csv'))))
ORDER = ['1', '2', '3'] + [f'S{n}' for n in range(2, 17)]
HEADS = {'two_stage': 'Two-stage', 'one_stage': 'One-stage', 'storey': "Storey's FDR"}


def name(table):
    return f'Web Table {table}' if table.startswith('S') else f'Table {table}'


def head(column):
    """two_stage_FDR -> 'Two-stage FDR', storey_R -> 'Storey R'; other names as they are."""
    for key, text in HEADS.items():
        if column.startswith(key):
            rest = column[len(key) + 1:]
            return text.replace("'s FDR", '') + ' ' + rest if rest else text
    return column


md = ['# Check results', '',
      'Written by `checks/run_checks.sh`: the outputs of the scripts in `scripts/`',
      '(run with their defaults) compared with the paper by code.', '',
      '## Tables 1 to 3 and Web Tables S2 to S16', '',
      'Every cell is compared with the number printed in the paper',
      '(`checks/paper_values.csv`), after rounding the computed value to the decimals',
      'of the printed cell. Each cell below reads paper / computed; a dash is a cell',
      'that is blank in the paper. Rows marked se are the standard errors, printed in',
      'parentheses in the paper.', '']
for n in ORDER:
    cells = [r for r in tables if r['table'] == n]
    if not cells:
        md += [f'### {name(n)}: not compared (results file not found)', '']
        continue
    columns = list(dict.fromkeys(r['column'] for r in cells))
    md += [f"### {name(n)}: {sum(r['equal'] == 'True' for r in cells)} of {len(cells)} cells equal", '',
           '| Row | ' + ' | '.join(head(c) for c in columns) + ' | Result |',
           '|---' * (len(columns) + 2) + '|']
    for row in dict.fromkeys(r['row'] for r in cells):
        of = {r['column']: r for r in cells if r['row'] == row}
        md.append(f'| {row} | ' +
                  ' | '.join(f"{of[c]['paper']} / {of[c]['computed']}" if c in of else '-' for c in columns) +
                  f" | {'equal' if all(r['equal'] == 'True' for r in of.values()) else 'differs'} |")
    md.append('')
md += ['## Figures 1 and 2 and Web Figures S1 to S12', '',
       'Every figure is compared with the published file (`archive/figures/`) by the',
       'marks in the two PDF files. Histograms: the number of bars of each barplot,',
       'and the largest leftover of the bar heights (after one scale per panel) and',
       'of the bar positions (after one affine map). Point-and-line plots: the number',
       'of points of each series, and the largest leftover of the point positions',
       'after one affine map per axis.', '']
for kind, title in [('histogram', 'Histograms'), ('points', 'Point-and-line plots')]:
    of = [r for r in figures if r['kind'] == kind]
    md += [f"{title}: {sum(r['result'] == 'matches' for r in of)} of {len(of)} match, "
           f"{sum(r['result'] == 'differs' for r in of)} differ, "
           f"{sum(r['result'] == 'not drawn' for r in of)} not drawn.", '']
md += ['| Figure | Published | Computed | Largest leftover | Result |', '|---|---|---|---|---|']
md += [f"| {r['figure']} | {r['published']} | {r['computed']} | {r['leftover']} | {r['result']} |" for r in figures]
print('\n'.join(md))
