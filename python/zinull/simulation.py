# -----------------------------------------------------------------------------
# File:    python/zinull/simulation.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Tables 1 and 2, Web Tables S2 to S15.
# Purpose: Simulation study: the 12 null settings, data from the mixture of a
#          ZIGP null and a shifted alternative, one replication (null fits,
#          the three procedures, bias of the estimates), the driver that runs
#          replications in parallel and stores them, and the tables of means
#          and standard deviations over replications.
#          The Python counterpart of R/simulation.R, with the same results
#          and the same files.
# Inputs:  archive/simulation/<setting>_C<C>.csv.gz (all replications of a
#          setting; unless rerun).
# Outputs: results/<table>.csv; in a rerun also
#          results/simulation/<setting>_C<C>/runNNNN.csv (one per replication).
# Usage:   from zinull.setup import PAPER, PATHS
#          from zinull.simulation import simulation_table
#          simulation_table("Table1", "simulation", "geometric", 0.80, 5, opts, PAPER, PATHS)
# Source:  R/simulation.R.
# -----------------------------------------------------------------------------
import csv
import gzip
import os

import numpy as np
import pandas as pd
from numba import njit

from . import rbase, rng
from .cutoff import CUTOFFS, mask_absent_parameters, select_null_fits
from .em import MODELS
from .fdr import local_fdr_rates, storey_rates
from .setup import log_message, output_file, parse_index_set
from .zigp import _generate_zigp

# The 12 null settings: alternative f1 (geometric or binomial count added to
# C + 1 + shift), proportion of nulls pi0, and the ZIGP null (eta, lambda, theta).
SIMULATION_SETTINGS = pd.DataFrame({
    "setting": ["ZIP1", "ZIGP1", "ZIGP2", "ZIP2", "ZIGP3", "ZIGP4",
                "ZIP3", "ZIGP5", "ZIGP6", "ZIP4", "ZIGP7", "ZIGP8"],
    "f1": np.repeat(["geometric", "binomial", "geometric", "binomial"], 3),
    "shift": np.repeat([0.0, -40.0, 0.0, -40.0], 3),
    "pi0": np.repeat([0.80, 0.35, 0.85, 0.80], 3),
    "eta": np.repeat([0.8, 0.4, 0.4, 0.4], 3),
    "lambda": [1.5, 1.5, 3, 1.5, 1, 4, 1.5, 1, 3, 1.5, 1, 3],
    "theta": [0, 0.3, 0.3, 0, 0.3, 0.3, 0, 0.15, 0.15, 0, 0.2, 0.2],
})

# Alternative distributions: Binomial(size, prob) and Geometric(prob).
F1_BINOMIAL = {"size": 250, "prob": 0.20}
F1_GEOMETRIC = {"prob": 0.08}

# Columns of a replication file summarized in the two kinds of tables.
TABLE_COLUMNS = {
    "simulation": [procedure + rate for procedure in ("two_stage_", "one_stage_", "storey_")
                   for rate in ("R", "FDR", "TPR")],
    "bias": ["bias_" + name for name in ("eta", "lambda", "theta", "pi0", "C")],
}

# Columns of a replication file that hold text.
TEXT_COLUMNS = ("cutoff", "model")


@njit(cache=True, error_model="numpy")
def _labelled_counts(state, n, C, shift, binomial, size, prob_binomial, prob_geometric, pi0, mu,
                     phi, omega):
    """The draws of `simulate_labelled_counts()` from the stream `state`."""
    value = np.zeros(n)
    label = np.zeros(n)
    filled = 0
    while filled < n:
        nonnull = rng._unif(state) >= pi0
        if nonnull:
            count = (C + 1.0 + shift + rng._rbinom(state, size, prob_binomial) * binomial
                     + rng._rgeom(state, prob_geometric) * (1.0 - binomial)
                     + _generate_zigp(state, 1, mu, phi, omega)[0])
            if not count > C:
                continue
        else:
            count = _generate_zigp(state, 1, mu, phi, omega)[0]
        value[filled] = count
        label[filled] = 1.0 if nonnull else 0.0
        filled += 1
    return value, label


def simulate_labelled_counts(n, C, setting, seed):
    """Simulated counts from the mixture of the null and the alternative.

    A position is null with probability pi0. A non-null count is
    C + 1 + shift plus a draw of f1 plus a null draw, and is drawn again
    (with a new label) when it is not above C. Both the binomial and the
    geometric draw are always made, so the random numbers are used in a
    fixed order.

    Parameters: n (number of positions), C (true cutoff: non-null counts
    are above C), setting (row of `SIMULATION_SETTINGS`), seed (seed of the
    data set, set right before the draws).
    Returns: dict with `value` (counts) and `label` (1 = non-null), float
    arrays of length n.
    """
    rng.set_seed(seed)
    theta = float(setting["theta"])
    value, label = _labelled_counts(
        rng._state, int(n), float(C), float(setting["shift"]),
        1.0 if setting["f1"] == "binomial" else 0.0, float(F1_BINOMIAL["size"]),
        float(F1_BINOMIAL["prob"]), float(F1_GEOMETRIC["prob"]), float(setting["pi0"]),
        float(setting["lambda"]) / (1 - theta), 1 / (1 - theta), float(setting["eta"]))
    return {"value": value, "label": label}


def simulation_replication(setting, C, run, settings):
    """One replication of a setting.

    Parameters: setting (row of `SIMULATION_SETTINGS`), C (true cutoff),
    run (replication number; the data are drawn with seed
    `settings["seed"] + run`), settings (settings of the paper: `seed`, `n`,
    `fdr_level`).
    Returns: DataFrame with one row per criterion and model: `run`,
    `cutoff`, `model`, R, FDR and TPR of the three procedures, the estimates
    `eta`, `lambda`, `theta`, `pi0`, `C`, `D` and their bias `bias_*`
    (estimate minus truth; NA for a parameter the model does not have).
    """
    data = simulate_labelled_counts(settings["n"], C, setting, settings["seed"] + run)
    null = data["value"][data["label"] == 0]
    nonnull = data["value"][data["label"] == 1]
    fits = select_null_fits(np.concatenate((null, nonnull)))
    sorted_null, sorted_nonnull = np.sort(null), np.sort(nonnull)
    rates = []
    for _, fit in fits.iterrows():
        rates.append({**local_fdr_rates(sorted_null, sorted_nonnull, fit, settings["fdr_level"]),
                      **storey_rates(null, nonnull, fit, settings["fdr_level"])})
    truth = {**{name: float(setting[name]) for name in ("eta", "lambda", "theta", "pi0")},
             "C": float(C)}
    masked = mask_absent_parameters(fits)
    bias = {}
    for name, true_value in truth.items():
        estimate = masked[name].to_numpy(dtype=np.float64)
        difference = estimate - true_value
        difference[rbase.is_na(estimate)] = rbase.NA
        bias["bias_" + name] = difference
    return pd.concat([pd.DataFrame({"run": [int(run)] * len(fits)}), fits[["cutoff", "model"]],
                      pd.DataFrame(rates),
                      fits[["eta", "lambda", "theta", "pi0", "C", "D"]], pd.DataFrame(bias)],
                     axis=1)


# ----- Stored replications -----

def format_exact(x):
    """A double with 17 significant digits, as R's `sprintf("%.17g", x)`.

    Parameters: x (float).
    Returns: the text; NA for a missing value, NaN, Inf or -Inf otherwise.
    """
    x = float(x)
    if x != x:
        return "NA" if rbase.is_na(x) else "NaN"
    if x in (float("inf"), float("-inf")):
        return "Inf" if x > 0 else "-Inf"
    return "%.17g" % x


def write_exact_csv(table, file):
    """Write the numbers of a table with 17 significant digits.

    Parameters: table (DataFrame; float columns are written with
    `format_exact()`, the others as they are, nothing in quotes), file
    (output file).
    Returns: None.
    """
    columns = []
    for name in table.columns:
        values = table[name]
        if pd.api.types.is_float_dtype(values):
            columns.append([format_exact(value) for value in values.to_numpy()])
        else:
            columns.append([str(value) for value in values])
    with open(file, "w") as f:
        f.write(",".join(table.columns) + "\n")
        for row in zip(*columns):
            f.write(",".join(row) + "\n")


def parse_exact(text):
    """Numbers from text with 17 significant digits, exactly.

    Python reads a decimal number as the nearest double, which is the number
    written (R's `parse_exact()` has to correct `as.numeric()` for this).

    Parameters: text (list of texts written by `format_exact()`).
    Returns: float array equal to the numbers written; NA for "NA".
    """
    return np.array([rbase.NA if item == "NA" else float(item) for item in text],
                    dtype=np.float64)


def read_exact_csv(files):
    """Read replication files written by `write_exact_csv()`.

    Parameters: files (one file, or a list of files with the same columns;
    a file ending in .gz is decompressed).
    Returns: DataFrame of all rows, sorted by `run`, numbers read exactly.
    """
    if isinstance(files, (str, os.PathLike)):
        files = [files]
    header, rows = None, []
    for file in files:
        opener = gzip.open if str(file).endswith(".gz") else open
        with opener(file, "rt", newline="") as f:
            reader = csv.reader(f)
            header = next(reader)
            rows.extend(reader)
    columns = {}
    for i, name in enumerate(header):
        cells = [row[i] for row in rows]
        columns[name] = cells if name in TEXT_COLUMNS else parse_exact(cells)
    runs = pd.DataFrame(columns)
    order = np.argsort(runs["run"].to_numpy(), kind="stable")
    return runs.iloc[order].reset_index(drop=True)


def setting_tag(setting, C):
    """Folder name of a setting and a true cutoff, for example "ZIGP1_C5".

    Parameters: setting (name of the setting), C (true cutoff).
    Returns: the name.
    """
    return "%s_C%d" % (setting, C)


def run_replications(group, C, run_file, opts, settings):
    """Run the replications that are not yet stored.

    Each replication seeds itself, so the result does not depend on the
    number of parallel processes or on the order of the replications.

    Parameters: group (rows of `SIMULATION_SETTINGS` to run), C (true
    cutoff), run_file (function of a setting name and a run number giving
    the file of the replication), opts (options: `runs`, `cores`), settings
    (settings of the paper).
    Returns: None. Raises RuntimeError when a replication fails.
    """
    runs = parse_index_set(opts["runs"], count=True)
    jobs = [(row, run) for row in range(len(group)) for run in runs
            if not os.path.exists(run_file(group["setting"].iloc[row], run))]
    log_message("%s, C = %d: %d replications to run", " ".join(group["setting"]), C, len(jobs))

    def job(item):
        setting = group.iloc[item[0]]
        replication = simulation_replication(setting, C, item[1], settings)
        write_exact_csv(replication, run_file(setting["setting"], item[1]))

    status = rbase.parallel_map(job, jobs, cores=opts["cores"])
    failed = [result for result in status if isinstance(result, Exception)]
    if failed:
        raise RuntimeError(f"{len(failed)} replications failed: {failed[0]!r}")


def stored_replications(group, C, opts, settings, paths):
    """Replications of the settings of a table.

    All archived replications (archive/simulation/<setting>_C<C>.csv.gz) or,
    with `opts["rerun"]`, the requested ones, run here when not yet stored
    in results/simulation/<setting>_C<C>/.

    Parameters: group (rows of `SIMULATION_SETTINGS`), C (true cutoff), opts
    (options: `rerun`, `runs`, `cores`), settings (settings of the paper),
    paths (project paths: `archive`, `results`).
    Returns: dict (one element per setting) of DataFrames with the rows of
    `simulation_replication()`, sorted by `run`.
    """
    def run_file(setting, run):
        return os.path.join(paths["results"], "simulation", setting_tag(setting, C),
                            "run%04d.csv" % run)

    if opts["rerun"]:
        for setting in group["setting"]:
            os.makedirs(os.path.dirname(run_file(setting, 1)), exist_ok=True)
        run_replications(group, C, run_file, opts, settings)
    runs = {}
    for setting in group["setting"]:
        if opts["rerun"]:
            runs[setting] = read_exact_csv([run_file(setting, run) for run in
                                            parse_index_set(opts["runs"], count=True)])
            continue
        archived = os.path.join(paths["archive"], "simulation", setting_tag(setting, C) + ".csv.gz")
        if not os.path.exists(archived):
            raise SystemExit(f"Error: {archived} is missing; run the replications here "
                             "with --rerun")
        runs[setting] = read_exact_csv(archived)
    return runs


# ----- Tables -----

def summarize_runs(runs, setting, columns):
    """Mean and standard deviation over replications, per criterion and model.

    Missing values are left out (the FDR is NaN in a replication without
    rejections). An FDR cell with no value (no replication with rejections,
    or the standard deviation of a single one) is 0.

    Parameters: runs (replications of one setting, sorted by `run`), setting
    (name of the setting), columns (columns of `runs` to summarize).
    Returns: DataFrame with the columns `setting`, `cutoff`, `model`,
    `statistic` ("mean" or "se") and `columns`; NA when a value is missing
    in every replication.
    """
    statistics = {"mean": rbase.r_mean, "se": rbase.r_sd}
    cutoff, model = runs["cutoff"].to_numpy(), runs["model"].to_numpy()
    rows = []
    for cutoff_name in CUTOFFS:
        for model_name in MODELS:
            cell = (cutoff == cutoff_name) & (model == model_name)
            for statistic, function in statistics.items():
                row = {"setting": setting, "cutoff": cutoff_name, "model": model_name,
                       "statistic": statistic}
                for column in columns:
                    v = runs[column].to_numpy(dtype=np.float64)[cell]
                    value = rbase.NA if np.isnan(v).all() else function(v, na_rm=True)
                    if column.endswith("_FDR") and np.isnan(value):
                        value = 0.0
                    row[column] = value
                rows.append(row)
    return pd.DataFrame(rows)


def simulation_table(name, kind, f1, pi0, C, opts, settings, paths):
    """Build a simulation table (R, FDR, TPR) or a bias table of one group of
    three settings and write it to results/.

    Parameters: name (table name, for example "Table1" or "WebTableS8"),
    kind ("simulation": rates of the three procedures; "bias": bias of the
    estimates of eta, lambda, theta, pi0 and C), f1 (alternative of the
    group: "geometric" or "binomial"), pi0 (proportion of nulls of the
    group), C (true cutoff: 5, or 1 for the contaminated null), opts
    (options: `rerun`, `runs`, `cores`), settings (settings of the paper),
    paths (project paths: `archive`, `results`).
    Returns: the table (DataFrame).
    """
    group = SIMULATION_SETTINGS[(SIMULATION_SETTINGS["f1"] == f1)
                                & (SIMULATION_SETTINGS["pi0"] == pi0)]
    runs = stored_replications(group, C, opts, settings, paths)
    table = pd.concat([summarize_runs(runs[setting], setting, TABLE_COLUMNS[kind])
                       for setting in group["setting"]], ignore_index=True)
    table.columns = [column.removeprefix("bias_") for column in table.columns]
    file = "%s_%s_%s-pi%.2f_C%d.csv" % (name, kind, f1, pi0, C)
    rbase.write_csv(table, output_file(paths["results"], file))
    first = runs[group["setting"].iloc[0]]
    log_message("%s: %d replications per setting", file, len(np.unique(first["run"])))
    print(table.to_string(index=False))
    return table
