# -----------------------------------------------------------------------------
# File:    checks/check_tables.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Compare every cell of Tables 1 to 3 and Web Tables S2 to S16
#          written by the table scripts with the number printed in the paper.
#          Each computed value is rounded to the decimals of the printed cell;
#          a printed NA is equal to a computed NA.
# Usage:   python3 checks/check_tables.py <results dir> <checks output dir>
# Inputs:  checks/paper_values.csv, <results dir>/Table<N>_*.csv and
#          <results dir>/WebTableS<N>_*.csv.
# Outputs: <checks output dir>/tables.csv; prints one line per table (a table
#          whose results file is missing is skipped); exit status 1 if a cell
#          differs.
# -----------------------------------------------------------------------------
import csv
import os
import sys

results, out = sys.argv[1:3]
here = os.path.dirname(os.path.abspath(__file__))
SETTINGS = ['geometric-pi0.80', 'binomial-pi0.35', 'geometric-pi0.85', 'binomial-pi0.80']
FILES = {'3': 'Table3_protein_rejections', 'S16': 'WebTableS16_protein_estimates'}
for tables, kind, c in [(['Table1', 'Table2', 'WebTableS6', 'WebTableS7'], 'simulation', 'C5'),
                        (['WebTableS2', 'WebTableS3', 'WebTableS4', 'WebTableS5'], 'bias', 'C5'),
                        (['WebTableS8', 'WebTableS9', 'WebTableS10', 'WebTableS11'], 'bias', 'C1'),
                        (['WebTableS12', 'WebTableS13', 'WebTableS14', 'WebTableS15'], 'simulation', 'C1')]:
    for table, setting in zip(tables, SETTINGS):
        FILES[table.replace('WebTable', '').replace('Table', '')] = f'{table}_{kind}_{setting}_{c}'
ORDER = ['1', '2', '3'] + [f'S{n}' for n in range(2, 17)]
MISSING = ('', 'NA', 'NaN')


def read(path):
    return list(csv.DictReader(open(path)))


def name(table):
    return f'Web Table {table}' if table.startswith('S') else f'Table {table}'


def plain(setting):
    """MZIP1 (the label of the C = 1 tables) and ZIP1 are the same setting."""
    return setting[1:] if setting.startswith('M') else setting


def index(table, rows):
    """Results rows by the labels of the paper cell they hold."""
    if table == '3':
        return {(r['domain'], r['cutoff'], r['procedure']): r for r in rows}
    if table == 'S16':
        return {(r['domain'], r['cutoff'], r['model']): r for r in rows}
    return {(plain(r['setting']), r['cutoff'], r['model'], r['statistic']): r for r in rows}


def computed_cell(table, rows, cell):
    """The computed value of a paper cell as text, 'missing' if not there."""
    if table == '3':
        row, column = rows.get((cell['setting'], cell['cutoff'], cell['column'])), cell['model']
    elif table == 'S16':
        row, column = rows.get((cell['setting'], cell['cutoff'], cell['model'])), cell['column']
    else:
        row = rows.get((plain(cell['setting']), cell['cutoff'], cell['model'], cell['statistic']))
        column = cell['column']
    if row is None or column not in row:
        return 'missing'
    value = row[column].strip()
    if value in MISSING:
        return 'NA'
    printed = cell['value']
    decimals = len(printed.split('.')[1]) if '.' in printed else 0
    if printed == 'NA':
        decimals = 4
    return f'{float(value):.{decimals}f}'


paper = read(os.path.join(here, 'paper_values.csv'))
rows = []
lines = []
for table in ORDER:
    path = os.path.join(results, FILES[table] + '.csv')
    if not os.path.exists(path):
        lines.append(f'{name(table)}: skipped, {FILES[table]}.csv not found')
        continue
    computed = index(table, read(path))
    cells = [c for c in paper if c['table'] == table]
    for c in cells:
        value = computed_cell(table, computed, c)
        rows.append({'table': table,
                     'row': ', '.join(c[k] for k in ('setting', 'cutoff', 'model', 'statistic')
                                      if not (k == 'statistic' and c[k] == 'value')),
                     'column': c['column'], 'paper': c['value'], 'computed': value,
                     'equal': value == c['value']})
    equal = sum(r['equal'] for r in rows if r['table'] == table)
    lines.append(f'{name(table)}: {equal} of {len(cells)} cells equal')

os.makedirs(out, exist_ok=True)
with open(os.path.join(out, 'tables.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['table', 'row', 'column', 'paper', 'computed', 'equal'])
    w.writeheader()
    w.writerows(rows)
print('\n'.join(lines))
sys.exit(0 if all(r['equal'] for r in rows) else 1)
