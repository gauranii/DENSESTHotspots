# -----------------------------------------------------------------------------
# File:    checks/check_python.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Compare the outputs of the R scripts with those of the Python
#          scripts, output by output: every cell of every csv file under
#          results/ (tables, stored replications) and the marks of every
#          figure (bars, points and lines; figure_compare.py).
# Usage:   python3 checks/check_python.py archive|rerun <R output dir>
#          <Python output dir> <comparison csv> [runs]
#          Run from the repository root (scripts.nix and scripts/ are read
#          for the list of outputs).
# Inputs:  <dir>/results/, <dir>/figures/, <dir>/logs/<script>.done (time and
#          cores) of both implementations, written by checks/check_python.sh.
# Outputs: the comparison csv, one row per output of the paper (archive) or
#          per group of settings and true cutoff (rerun): files, cells and
#          marks compared, how many are not identical, the largest
#          difference, the verdict and the run times. Prints the same; exit
#          status 1 unless every output that was built is identical.
#          Tolerance: none for the csv files (the same header, the same rows,
#          every cell the same text; the files are also compared byte for
#          byte). The marks of a figure must correspond one to one, with the
#          same colours; their positions agree to 0.02 pt of the R page after
#          one straight-line map per axis (R prints its coordinates with two
#          decimals). A figure file "WebFigS..." and a file "FigS..." are the
#          same figure.
# -----------------------------------------------------------------------------
import csv
import glob
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_compare import POSITION_TOLERANCE, r_vs_python_bars, r_vs_python_points  # noqa: E402

mode, r_dir, py_dir, out_file = sys.argv[1:5]
runs = sys.argv[5] if len(sys.argv) > 5 else "1:2"


def scripts():
    """Command name -> script name, from scripts.nix."""
    with open("scripts.nix") as f:
        return dict(re.findall(r'^\s*(\w+) = "(\w+)";', f.read(), flags=re.M))


def script_text(script):
    with open(os.path.join("scripts", script + ".R")) as f:
        return f.read()


def natural(name):
    """Sort key: tables and figures of the paper, then the web ones, by number."""
    web, kind, number = re.fullmatch(r"(web)?(table|fig)S?(\d+)", name).groups()
    return bool(web), kind != "table", int(number)


def archive_outputs():
    """One output per script: (name, scripts, csv patterns, figure patterns,
    kind of figure, where the output comes from)."""
    outputs = []
    for name, script in sorted(scripts().items(), key=lambda item: natural(item[0])):
        label = (name[0].upper() + name[1:]).replace("fig", "Fig").replace("table", "Table")
        if "_simulation_histograms" in script:
            figure = label.replace("WebFig", "Fig")
            outputs.append((name, [script], [], [figure + "[A-Z]_*.pdf"], "bars",
                            "seeded simulated data"))
        elif "_protein_histograms" in script:
            outputs.append((name, [script], [], [label + "[A-Z]_*.pdf"], "bars",
                            "computed from data"))
        elif "_protein_null_fit" in script:
            figure = label.replace("WebFig", "Fig")
            outputs.append((name, [script], [], [figure + "[A-Z]_*.pdf"], "points",
                            "computed from data"))
        elif "simulation_table(" in script_text(script):
            outputs.append((name, [script], [label + "_*.csv"], [], None, "archived runs"))
        else:
            outputs.append((name, [script], [label + "_*.csv"], [], None, "computed from data"))
    return outputs


# Settings of the four groups of the simulation study (alternative, pi0), in
# the order of R/simulation.R.
GROUPS = {("geometric", "0.80"): ("ZIP1", "ZIGP1", "ZIGP2"),
          ("binomial", "0.35"): ("ZIP2", "ZIGP3", "ZIGP4"),
          ("geometric", "0.85"): ("ZIP3", "ZIGP5", "ZIGP6"),
          ("binomial", "0.80"): ("ZIP4", "ZIGP7", "ZIGP8")}


def rerun_outputs():
    """One output per group of settings and true cutoff: its two tables and
    the stored replications of its three settings."""
    groups = {}
    for name, script in scripts().items():
        call = re.search(r'simulation_table\("(\w+)", "\w+", f1 = "(\w+)", pi0 = ([\d.]+), '
                         r'C = (\d)', script_text(script))
        if call:
            table, f1, pi0, C = call.groups()
            order = (C != "5", list(GROUPS).index((f1, pi0)))
            groups.setdefault((order, f1, pi0, C), []).append((name, script, table))
    outputs = []
    for (_, f1, pi0, C), tables in sorted(groups.items()):
        tables.sort(key=lambda item: natural(item[0]))
        settings = GROUPS[(f1, pi0)]
        patterns = [table + "_*.csv" for _, _, table in tables]
        patterns += [f"simulation/{setting}_C{C}/run*.csv" for setting in settings]
        outputs.append(("+".join(name for name, _, _ in tables), [s for _, s, _ in tables],
                        patterns, [], None,
                        f"{f1}, pi0 = {pi0}, C = {C}: replications {runs} of "
                        + ", ".join(settings)))
    return outputs


def files_of(folder, patterns):
    """Files under a folder that match the patterns: {name with WebFigS
    written as FigS: relative path}."""
    found = {}
    for pattern in patterns:
        for prefix in ("", "Web"):
            for path in glob.glob(os.path.join(folder, prefix + pattern)):
                relative = os.path.relpath(path, folder)
                found[re.sub(r"^WebFigS", "FigS", relative)] = relative
    return found


def number(cell):
    """The number of a cell, NaN for a missing value, None for text."""
    if cell in ("NA", "NaN"):
        return math.nan
    try:
        return float(cell)
    except ValueError:
        return None


class Comparison:
    """Running totals of the comparison of one output."""

    def __init__(self):
        self.files = self.values = self.unequal = 0
        self.max_abs = 0.0
        self.problems = []

    def problem(self, text):
        self.problems.append(text)

    def csv_file(self, r_file, py_file, name):
        """Same header, same number of rows, every cell the same text, and
        the same bytes."""
        with open(r_file, newline="") as f:
            r_rows = list(csv.reader(f))
        with open(py_file, newline="") as f:
            py_rows = list(csv.reader(f))
        self.files += 1
        if r_rows[0] != py_rows[0] or len(r_rows) != len(py_rows):
            return self.problem(f"{name}: different columns or number of rows")
        for r_row, py_row in zip(r_rows[1:], py_rows[1:]):
            for a, b in zip(r_row, py_row):
                self.values += 1
                if a != b:
                    self.unequal += 1
                    x, y = number(a), number(b)
                    if x is not None and y is not None and not (math.isnan(x) or math.isnan(y)):
                        self.max_abs = max(self.max_abs, abs(x - y))
                    self.problem(f"{name}: cells differ")
        with open(r_file, "rb") as f, open(py_file, "rb") as g:
            if f.read() != g.read() and not self.problems:
                self.problem(f"{name}: same cells, other bytes (quotes or line ends)")

    def figure(self, r_file, py_file, name, kind):
        """The marks of a figure correspond and agree within the tolerance."""
        if kind == "bars":
            in_r, in_py, compared, difference = r_vs_python_bars(r_file, py_file)
        else:
            counts, compared, difference = r_vs_python_points(r_file, py_file)
            in_r, in_py = sum(n for _, n, _ in counts), sum(n for _, _, n in counts)
        self.files += 1
        self.values += compared
        if difference is None:
            self.unequal += max(in_r, in_py, 1)
            self.problem(f"{name}: {in_r} and {in_py} {kind} that do not correspond")
        else:
            self.max_abs = max(self.max_abs, difference)
            if difference > POSITION_TOLERANCE:
                self.unequal += compared
                self.problem(f"{name}: {kind} differ by {difference:.3g} pt")


def run_time(folder, script_names):
    """Seconds (summed over the scripts) and cores, from the .done files;
    None when a script was not run."""
    seconds, cores = 0.0, ""
    for script in script_names:
        try:
            with open(os.path.join(folder, "logs", script + ".done")) as f:
                time, cores = f.read().strip().split(",")
        except OSError:
            return None
        seconds += float(time)
    return f"{seconds:.1f}", cores


rows = []
for output, script_names, tables, figures, kind, source in (
        archive_outputs() if mode == "archive" else rerun_outputs()):
    r_time, py_time = run_time(r_dir, script_names), run_time(py_dir, script_names)
    total = Comparison()
    if r_time is None and py_time is None:  # the protein domain data are not there
        row = {"files": 0, "values": 0, "not_identical": 0, "max_abs_diff": "",
               "verdict": "skipped (no data)", "r_seconds": "", "python_seconds": "",
               "cores": "", "problems": ""}
    else:
        for folder, patterns in (("results", tables), ("figures", figures)):
            r_files = files_of(os.path.join(r_dir, folder), patterns)
            py_files = files_of(os.path.join(py_dir, folder), patterns)
            if sorted(r_files) != sorted(py_files):
                total.problem(f"the two implementations wrote different files in {folder}/")
            for file in sorted(set(r_files) & set(py_files)):
                paths = (os.path.join(r_dir, folder, r_files[file]),
                         os.path.join(py_dir, folder, py_files[file]))
                if folder == "results":
                    total.csv_file(*paths, file)
                else:
                    total.figure(*paths, file, kind)
        if total.files == 0:
            total.problem("no files to compare")
        row = {"files": total.files, "values": total.values, "not_identical": total.unequal,
               "max_abs_diff": f"{total.max_abs:.3g}" + (" pt" if figures else ""),
               "verdict": "DIFFERENT" if total.problems else "identical",
               "r_seconds": r_time[0] if r_time else "", "python_seconds": py_time[0] if py_time else "",
               "cores": (r_time or py_time)[1],
               "problems": "; ".join(sorted(set(total.problems))[:5])}
    rows.append({"output": output, "mode": mode, "source": source, **row})
    print(f"{output:24} {mode:8} files {row['files']:4d} values {row['values']:6d} "
          f"not identical {row['not_identical']:5d}  {row['verdict']}")
    for problem in sorted(set(total.problems))[:5]:
        print("        ", problem)

os.makedirs(os.path.dirname(os.path.abspath(out_file)), exist_ok=True)
with open(out_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
sys.exit(1 if any(row["verdict"] == "DIFFERENT" for row in rows) else 0)
