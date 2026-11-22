# -----------------------------------------------------------------------------
# File:    python/zinull/setup.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. All figures and tables.
# Purpose: Imported first by every script: finds the repository root and
#          defines the settings of the paper, the project paths and the
#          command line parser. The Python counterpart of R/setup.R.
# Inputs:  None.
# Outputs: None.
# Usage:   from zinull.setup import PAPER, PATHS, parse_args
# Source:  R/setup.R.
# -----------------------------------------------------------------------------
import os
import re
import sys
import time

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

# ----- Settings of the paper, passed to the functions as `settings` -----
PAPER = {
    "seed": 889191,              # replication k of a simulation uses seed + k
    "figure_seed": 889191 + 111,  # seed of the simulated data in the histograms
    "n": 1000,                   # positions per simulated data set
    "replications": 1000,        # replications per simulation setting
    "fdr_level": 0.05,           # FDR level of the simulation study
    "fdr_level_real": 0.10,      # FDR level of the protein domain data
}


# ----- Paths -----

def project_paths(root, output=None):
    """Folders and files of the project.

    Parameters: root (repository root), output (folder in which results/
    and figures/ are written; ZINULL_OUTPUT, or else the root, by default).
    Returns: dict of paths: `data` (ZINULL_DATA when set), `archive`,
    `results`, `figures`.
    """
    if output is None:
        output = os.environ.get("ZINULL_OUTPUT", root)
    return {"data": os.environ.get("ZINULL_DATA", os.path.join(root, "data", "realdata.txt")),
            "archive": os.path.join(root, "archive"),
            "results": os.path.join(output, "results"),
            "figures": os.path.join(output, "figures")}


PATHS = project_paths(REPO_ROOT)


def require_data(path):
    """Stop with a message when the protein domain data file is not there.

    Parameters: path (path of data/realdata.txt).
    Returns: `path`. Exits with status 1 and the message of R/setup.R when
    the file does not exist.
    """
    if not os.path.exists(path):
        raise SystemExit("Error: The protein domain counts (data/realdata.txt) are not included "
                         "in this repository. Place the file in data/ of the directory the "
                         "command runs from, or set ZINULL_DATA to its path.")
    return path


def output_file(folder, name):
    """Path of an output file, creating its folder when needed.

    Parameters: folder (output folder, for example `paths["results"]`), name
    (file name).
    Returns: the path folder/name.
    """
    os.makedirs(folder, exist_ok=True)
    return os.path.join(folder, name)


# ----- Command line -----

def show_help(args=None, script=None):
    """Answer `--help`: print the header of the script (file, purpose,
    inputs, outputs, usage) and exit.

    Parameters: args (command line arguments; those of the process by
    default), script (path of the script; the one that was started by
    default).
    Returns: None when `--help` is not among the arguments.
    """
    args = sys.argv[1:] if args is None else args
    if "--help" not in args and "-h" not in args:
        return
    script = sys.argv[0] if script is None else script
    with open(script) as f:
        for line in f:
            if not line.startswith("#"):
                break
            if not line.startswith("# ---"):
                print(line[2:].rstrip("\n"))
    raise SystemExit(0)


def parse_args(defaults, args=None):
    """Read `--name=value` arguments and bare `--name` flags.

    Every script accepts `--rerun`: recompute from scratch instead of
    reading the archived runs; and `--help` (see `show_help()`).

    Parameters: defaults (dict of options and their default values; a value
    given on the command line is converted to the type of its default, a
    flag sets a boolean default to True), args (command line arguments).
    Returns: `defaults` with the values given on the command line, and rerun.
    """
    args = sys.argv[1:] if args is None else args
    show_help(args)
    opts = {"rerun": False, **defaults}
    pairs = [re.fullmatch(r"--([a-z_]+)(=(.*))?", arg) for arg in args]
    if any(pair is None for pair in pairs):
        raise SystemExit("Error: arguments must look like --name=value or --name")
    for pair in pairs:
        name, is_flag = pair.group(1), pair.group(2) is None
        if name not in opts:
            raise SystemExit(f"Error: unknown argument --{name}")
        if is_flag != isinstance(opts[name], bool):
            problem = "needs a value" if is_flag else "takes no value"
            raise SystemExit(f"Error: --{name} {problem}")
        if is_flag:
            opts[name] = True
        elif isinstance(opts[name], int):
            opts[name] = int(float(pair.group(3)))
        else:
            opts[name] = type(opts[name])(pair.group(3))
    return opts


def parse_index_set(text, count=False):
    """Turn an index text into indices.

    Parameters: text ("3:7" for 3 to 7, "2,5,9", or a single number N),
    count (if True a single number N means 1 to N).
    Returns: list of integers.
    """
    if re.fullmatch(r"[0-9]+:[0-9]+", text):
        first, last = (int(end) for end in text.split(":"))
        step = 1 if last >= first else -1
        return list(range(first, last + step, step))
    index = [int(item) for item in text.split(",")]
    return list(range(1, index[0] + 1)) if count and len(index) == 1 else index


def log_message(msg, *values):
    """Print a time-stamped message.

    Parameters: msg (format of the message, with % fields), values (values
    for the format).
    Returns: None.
    """
    print(time.strftime("%Y-%m-%d %H:%M:%S"), "-", msg % values, flush=True)
