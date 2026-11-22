# -----------------------------------------------------------------------------
# File:    checks/write_python.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Write the comparison of the R and the Python implementation as
#          Markdown (checks/PYTHON.md) from the tables of
#          checks/check_python.py.
# Usage:   python3 checks/write_python.py <python_archive.csv>
#          <python_rerun.csv> > checks/PYTHON.md
#          A table that does not exist is left out.
# -----------------------------------------------------------------------------
import csv
import os
import re
import sys

INTRO = {
    "archive": """## With the defaults of the scripts

Both implementations run every script with its defaults, as `nix flake check`
runs them (`bash checks/check_python.sh archive <work dir>`): the 16
simulation tables from the archived runs (1000 replications per setting), the
histograms of simulated data from their seed, and Table 3, Web Table S16,
Figure 2 and Web Figures S7 to S12 from the protein domain counts
(`data/realdata.txt`; these are skipped when the file is not there).
""",
    "rerun": """## Recomputed by simulation

Both implementations run the 16 simulation tables with `--rerun` for the
first replications of all 24 combinations of a setting and a true cutoff,
with the seeds of the full runs (`bash checks/check_python.sh rerun <work
dir>`). Compared are the stored file of every replication (rejections, FDR
and TPR of the three procedures, estimates and their bias, with 17
significant digits) and the two tables of each group of settings. The times
are those of the two table scripts of a group together; the second one reads
the replications the first one stored.
""",
}


def name_of(output):
    """table1 -> Table 1, webfigS10 -> Web Figure S10; joined by a comma."""
    names = []
    for part in output.split("+"):
        web, kind, number = re.fullmatch(r"(web)?(table|fig)(S?\d+)", part).groups()
        names.append(("Web " if web else "") + ("Table " if kind == "table" else "Figure ") + number)
    return ", ".join(names)


def section(file):
    with open(file) as f:
        rows = list(csv.DictReader(f))
    mode = rows[0]["mode"]
    built = [row for row in rows if not row["verdict"].startswith("skipped")]
    identical = sum(row["verdict"] == "identical" for row in built)
    cores = next((row["cores"] for row in built if row["cores"]), "")
    summary = f"{identical} of {len(built)} outputs identical in every value ({cores} cores)."
    if len(built) < len(rows):
        summary += (f" {len(rows) - len(built)} outputs of the protein domain data were "
                    "skipped: the data file was not there.")
    lines = [INTRO[mode], summary, "",
             "| Output | From | Files | Values | Not identical | Largest difference | Verdict | "
             "R (s) | Python (s) |",
             "|---|---|---|---|---|---|---|---|---|"]
    for row in rows:
        lines.append(f"| {name_of(row['output'])} | {row['source']} | {row['files']} | "
                     f"{row['values']} | {row['not_identical']} | {row['max_abs_diff']} | "
                     f"{row['verdict']} | {row['r_seconds']} | {row['python_seconds']} |")
    problems = [f"- {name_of(row['output'])}: {row['problems']}" for row in rows if row["problems"]]
    if problems:
        lines += ["", "Differences:", ""] + problems
    return "\n".join(lines) + "\n"


print("""# R and Python

Written by `checks/check_python.sh`: the outputs of the R scripts (`scripts/`)
and of the Python scripts (`python/scripts/`) compared by code
(`checks/check_python.py`). Every csv file under `results/` is compared cell
by cell, without tolerance: the same header, the same rows, every cell the
same text, and the same bytes. Every figure is compared by the marks read from
the two PDF files (`checks/figure_compare.py`): the same number of bars, or of
points and line vertices in each series, in the same order and with the same
colours; bar ends, bar heights, point centres and line vertices agree to
0.02 pt of the R page after one straight-line map per axis (one scale for the
bar heights), since R and matplotlib lay their pages out differently and R
prints its coordinates with two decimals. "Values" counts the compared cells
and mark coordinates; "Largest difference" is in the units of the cells, or in
pt for a figure. Times are wall-clock seconds on an 8 core Apple M2 laptop
shared with other jobs; the Python times include loading its compiled (numba)
functions.

The tests in `python/tests/` (also run by `nix flake check`) compare every
function of `python/zinull/` with values computed by R
(`python/tests/make_r_values.R`), double for double.
""")
for file in sys.argv[1:3]:
    if os.path.exists(file):
        print(section(file))
