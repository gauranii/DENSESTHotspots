# -----------------------------------------------------------------------------
# File:    checks/extract_paper_values.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Pull every printed cell of Tables 1 to 3 (published article) and of
#          Web Tables S2 to S16 (LaTeX source of the web supplement) into
#          checks/paper_values.csv, one row per cell: table, setting (or
#          protein domain), cutoff, model, column, statistic, value.
#          Standard errors are stored without their parentheses; blank cells
#          are not stored. The LaTeX source of the article is read as a
#          cross-check of Tables 1 to 3; cells that differ from the published
#          text, and cells that could not be read, are listed on stderr.
# Usage:   python3 checks/extract_paper_values.py [article.txt [article.tex
#          [supplement.tex]]] > checks/paper_values.csv
#          article.txt is `pdftotext -layout` of the published PDF.
# -----------------------------------------------------------------------------
import csv
import os
import re
import sys

DRIVE = os.path.expanduser('~/NightShift/N26-drive/ms')
defaults = [os.path.join(DRIVE, 'pub.txt'),
            os.path.join(DRIVE, 'THESIS/PROPOSAL/BIOMETRICS/biometrics0710.tex'),
            os.path.join(DRIVE, 'BIOMETRICS/Web-based Supplementary Materials 0710.tex')]
pub_path, tex_path, supp_path = sys.argv[1:4] + defaults[len(sys.argv) - 1:]

MODELS = ['ZIGP', 'ZIP', 'GP', 'P']
SIMULATION = [f'{p}_{s}' for p in ('two_stage', 'one_stage', 'storey') for s in ('R', 'FDR', 'TPR')]
BIAS = ['eta', 'lambda', 'theta', 'pi0', 'C']
ESTIMATES = ['eta', 'lambda', 'theta', 'pi0', 'C', 'D']
PROCEDURES = ['one_stage', 'two_stage', 'storey']  # column groups of Table 3, left to right
NUMBER = r'-?\d+(?:\.\d+)?'
# Setting of each table: f1, pi0, C, settings in the rows (C = 1: MZIP1, MZIGP1, ...)
GROUPS = [('Geometric', '0.80', ['ZIP1', 'ZIGP1', 'ZIGP2']), ('Binomial', '0.35', ['ZIP2', 'ZIGP3', 'ZIGP4']),
          ('Geometric', '0.85', ['ZIP3', 'ZIGP5', 'ZIGP6']), ('Binomial', '0.80', ['ZIP4', 'ZIGP7', 'ZIGP8'])]
EXPECTED = {}
for tables, kind, c in [(['1', '2', 'S6', 'S7'], 'Numerical Comparison', '5'),
                        (['S2', 'S3', 'S4', 'S5'], 'Bias of Parameter Estimates', '5'),
                        (['S8', 'S9', 'S10', 'S11'], 'Bias of Parameter Estimates', '1'),
                        (['S12', 'S13', 'S14', 'S15'], 'Numerical Comparison', '1')]:
    for table, (f1, pi0, settings) in zip(tables, GROUPS):
        EXPECTED[table] = (kind, f1, pi0, c, [('M' if c == '1' else '') + s for s in settings])

cells = {}   # (table, setting, cutoff, model, column, statistic) -> printed string
notes = []   # cells that could not be read, caption mismatches


def store(table, setting, cutoff, model, column, statistic, value, into=cells):
    into[(table, setting, cutoff, model, column, statistic)] = value


# ----- LaTeX tables -----

def tex_rows(body):
    """Rows of a tabular as lists of cells, without rules and font commands."""
    body = re.sub(r'(?<!\\)%.*', '', body)
    body = re.sub(r'\\(?:textbf|texttt)\{([^{}]*)\}', r'\1', body)
    body = re.sub(r'\\hline|\\multicolumn\{\d+\}\{[^{}]*\}', '', body)
    body = re.sub(r'\(\s*(' + NUMBER + r')\s*\)', r'(\1)', body)
    return [[c.strip() for c in row.split('&')] for row in body.split('\\\\') if row.strip()]


def tex_grid(table, body, columns, into):
    """A simulation or bias table: value rows and standard error rows."""
    setting = cutoff = model = None
    means = {}
    for row in tex_rows(body):
        if len(row) != 3 + len(columns):
            if any(re.search(r'\d\.\d', c) for c in row):
                notes.append(f'Table {table}: row with {len(row)} cells not read: {" & ".join(row)}')
            continue
        if not any(row[3:]) or 'widehat' in ''.join(row):
            continue
        m = re.search(r'(M?ZIG?P)\$_\{?(\d+)\}?\$', row[0])
        setting = m.group(1) + m.group(2) if m else setting
        m = re.search(r'C\}?_\{?([12])', row[1])
        cutoff = 'C' + m.group(1) if m else cutoff
        statistic = 'mean' if row[2] else 'se'
        model = row[2] if row[2] else model
        where = f'Table {table}, {setting}, {cutoff}, {model}'
        if model not in MODELS or (row[2] and (row[0] or row[1]) and model != 'ZIGP'):
            notes.append(f'{where}: row labels not read: {" & ".join(row[:3])}')
            continue
        for column, cell in zip(columns, row[3:]):
            key = (setting, cutoff, model, column)
            if statistic == 'mean':
                means[key] = cell
                if cell == 'NA' or re.fullmatch(NUMBER, cell):
                    store(table, *key, 'mean', cell, into)
                else:
                    notes.append(f'{where}, {column}: value "{cell}" not read')
            elif re.fullmatch(r'\(' + NUMBER + r'\)', cell):
                store(table, *key, 'se', cell[1:-1], into)
            elif re.fullmatch(NUMBER, cell):
                store(table, *key, 'se', cell, into)
                notes.append(f'{where}, {column}: standard error {cell} printed without parentheses (stored)')
            elif cell == '' and means.get(key) == 'NA':
                pass  # no standard error under NA
            elif cell == '':
                notes.append(f'{where}, {column}: standard error blank under {means.get(key)} (not stored)')
            else:
                notes.append(f'{where}, {column}: standard error "{cell}" not read')


def tex_rejections(body, into):
    """Table 3: number of rejections, 12 counts per row."""
    domain = None
    for row in tex_rows(body):
        if len(row) == 14 and all(re.fullmatch(r'\d+', c) for c in row[2:]):
            domain = row[0] or domain
            cutoff = 'C' + re.search(r'C\}?_\{?([12])', row[1]).group(1)
            for i, cell in enumerate(row[2:]):
                store('3', domain, cutoff, MODELS[i % 4], PROCEDURES[i // 4], 'value', cell, into)


def tex_estimates(table, body):
    """Web Table S16: estimates with C1 (six columns), then with C2."""
    domain = None
    for row in tex_rows(body):
        if len(row) != 14 or row[1] not in MODELS:
            continue
        domain = row[0] or domain
        for i, cell in enumerate(row[2:]):
            if cell == 'NA' or re.fullmatch(NUMBER, cell):
                store(table, domain, f'C{i // 6 + 1}', row[1], ESTIMATES[i % 6], 'value', cell)
            else:
                notes.append(f'Table {table}, {domain}, {row[1]}, {ESTIMATES[i % 6]}: value "{cell}" not read')


def check_caption(table, caption, body_cells):
    """Compare the caption and the row labels of a table with EXPECTED."""
    kind, f1, pi0, c, settings = EXPECTED[table]
    m = re.search(r'shifted (\w+)\(.*?\\pi_0 = ([\d.]+)\$,? and \$C = (\d)\$', caption)
    read = (m.group(1), m.group(2), m.group(3)) if m else None
    if kind not in caption or read != (f1, pi0, c):
        notes.append(f'Table {table}: caption reads {read}, expected {kind}, {(f1, pi0, c)}')
    rows = sorted({k[1] for k in body_cells if k[0] == table})
    if rows != sorted(settings):
        notes.append(f'Table {table}: rows are {rows}, expected {settings}')


# ----- Web supplement: Web Tables S2 to S16 -----

supp = open(supp_path).read()
for m in re.finditer(r'\\begin\{tabular\}\{([^{}]*)\}(.*?)\\end\{tabular\}', supp, re.S):
    spec, body = m.groups()
    captions = list(re.finditer(r'Table (\d+):(.*)', supp[:m.start()]))
    number, caption = captions[-1].groups()
    table = 'S' + number
    if spec == 'lll|ccccc':
        tex_grid(table, body, BIAS, cells)
    elif spec == 'l|l|l|ccc|ccc|ccc':
        tex_grid(table, body, SIMULATION, cells)
    elif spec == 'l|l|rrrrrr|rrrrrr':
        tex_estimates(table, body)
        continue
    else:
        continue  # Web Table S1 (model list) and figure layouts
    check_caption(table, caption, cells)

# ----- Published article: Tables 1 to 3 -----

lines = open(pub_path).read().replace('\u02c6', '^').split('\n')
starts = {n: i for i, l in enumerate(lines) for n in '123' if re.fullmatch(r'\s*Table ' + n + r'\s*', l)}
for table in '12':
    setting = cutoff = model = None
    read = 0
    captions = ' '.join(lines[starts[table]:starts[table] + 6])
    for line in lines[starts[table]:]:
        tokens = line.split()
        if tokens and re.fullmatch(r'ZIG?P\d', tokens[0]):
            setting = tokens.pop(0)
        if tokens and re.fullmatch(r'C\^[12]', tokens[0]):
            cutoff = 'C' + tokens.pop(0)[-1]
        if len(tokens) == 10 and tokens[0] in MODELS and all(re.fullmatch(NUMBER, t) for t in tokens[1:]):
            model, statistic, values = tokens[0], 'mean', tokens[1:]
        elif len(tokens) == 9 and all(re.fullmatch(r'\(' + NUMBER + r'\)', t) for t in tokens):
            statistic, values = 'se', [t[1:-1] for t in tokens]
        else:
            continue
        for column, value in zip(SIMULATION, values):
            store(table, setting, cutoff, model, column, statistic, value)
        read += 1
        if read == 48:  # 3 settings x 2 cutoffs x 4 models x (value, standard error)
            break
    kind, f1, pi0, c, settings = EXPECTED[table]
    m = re.search(r'shifted (\w+)\(.*?π0 = ([\d.]+),? and C = (\d)', captions)
    if not m or m.groups() != (f1, pi0, c):
        notes.append(f'Table {table}: caption reads {m and m.groups()}, expected {(f1, pi0, c)}')
    rows = sorted({k[1] for k in cells if k[0] == table})
    if rows != sorted(settings) or read != 48:
        notes.append(f'Table {table}: {read} rows read for {rows}, expected 48 rows for {settings}')

domain = None
for line in lines[starts['3']:]:
    tokens = line.split()
    if tokens and re.fullmatch(r'(cd|pfam)\d+', tokens[0]):
        domain = tokens.pop(0)
    if len(tokens) == 13 and re.fullmatch(r'C\^[12]', tokens[0]) and all(t.isdigit() for t in tokens[1:]):
        for i, value in enumerate(tokens[1:]):
            store('3', domain, 'C' + tokens[0][-1], MODELS[i % 4], PROCEDURES[i // 4], 'value', value)
    if sum(k[0] == '3' for k in cells) == 144:  # 6 domains x 2 cutoffs x 12 counts
        break

# ----- Cross-check: LaTeX source of the article against the published text -----

tex = open(tex_path).read()
source = {}
for table, label in [('1', 'tab:zip1'), ('2', 'tab:zip2'), ('3', 'tab:prorej')]:
    body = tex[tex.index('\\label{' + label + '}'):]
    body = body[body.index('\\begin{tabular}'):body.index('\\end{tabular}')]
    if table == '3':
        tex_rejections(body, source)
    else:
        tex_grid(table, body, SIMULATION, source)
published = {k: v for k, v in cells.items() if k[0] in '123'}
for key in sorted(set(source) | set(published)):
    if source.get(key) != published.get(key):
        notes.append(f'Table {key[0]}, {", ".join(key[1:])}: LaTeX source {source.get(key)}, '
                     f'published {published.get(key)}')

# ----- Output -----

w = csv.writer(sys.stdout, lineterminator='\n')
w.writerow(['table', 'setting', 'cutoff', 'model', 'column', 'statistic', 'value'])
order = ['1', '2', '3'] + [f'S{n}' for n in range(2, 17)]
cells = {k: v for t in order for k, v in cells.items() if k[0] == t}
w.writerows(list(k) + [v] for k, v in cells.items())
counts = {}
for k in cells:
    counts[k[0]] = counts.get(k[0], 0) + 1
print('cells per table: ' + ', '.join(f'{t} {n}' for t, n in counts.items()), file=sys.stderr)
for note in notes:
    print(note, file=sys.stderr)
