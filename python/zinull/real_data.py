# -----------------------------------------------------------------------------
# File:    python/zinull/real_data.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Table 3, Web Table S16.
# Purpose: Protein domain data: mutation counts per position of six domains,
#          the null fits of each domain, and the tables of the numbers of
#          rejected positions and of the estimates.
#          The Python counterpart of R/real_data.R, with the same results.
# Inputs:  data/realdata.txt (columns id and one per domain, NA padded).
# Outputs: None (the scripts write the tables).
# Usage:   from zinull.setup import PATHS
#          from zinull.real_data import protein_rejections
#          protein_rejections(PATHS["data"], 0.10)
# Source:  R/real_data.R.
# -----------------------------------------------------------------------------
import numpy as np
import pandas as pd

from .cutoff import CUTOFFS, mask_absent_parameters, select_null_fits
from .em import MODELS
from .fdr import one_stage_count, storey_count, two_stage_count
from .setup import require_data


def read_domains(path):
    """Mutation counts of the protein domains.

    Parameters: path (file with one column per domain after the `id`
    column, separated by white space, NA where a domain has no position).
    Returns: dict (one element per domain, in the order of the columns) of
    integer arrays of counts.
    """
    with open(require_data(path)) as f:
        header = f.readline().split()
        rows = [line.split() for line in f if line.strip()]
    return {domain: np.array([int(row[i]) for row in rows if row[i] != "NA"], dtype=np.int64)
            for i, domain in enumerate(header) if i > 0}


def fit_domain(x):
    """Null fits of one domain; the bound D is at most the largest count.

    Parameters: x (counts of the domain).
    Returns: DataFrame as returned by `select_null_fits()`.
    """
    return select_null_fits(x, cap=float(np.max(x)))


def protein_rejections(path, level):
    """Numbers of rejected positions of every domain (Table 3).

    Parameters: path (data file), level (FDR level).
    Returns: DataFrame with the columns `domain`, `cutoff`, `procedure`
    ("one_stage", "two_stage", "storey") and one column per null model.
    """
    procedures = {"one_stage": one_stage_count, "two_stage": two_stage_count,
                  "storey": storey_count}
    rows = []
    for domain, x in read_domains(path).items():
        fits = fit_domain(x)
        for cutoff in CUTOFFS:
            for procedure, count in procedures.items():
                row = {"domain": domain, "cutoff": cutoff, "procedure": procedure}
                for model in MODELS:
                    fit = fits[(fits["cutoff"] == cutoff) & (fits["model"] == model)].iloc[0]
                    row[model] = float(count(x, fit, level))
                rows.append(row)
    return pd.DataFrame(rows)


def protein_estimates(path):
    """Null estimates of every domain (Web Table S16).

    Parameters: path (data file).
    Returns: DataFrame with the columns `domain`, `cutoff`, `model`, `eta`,
    `lambda`, `theta`, `pi0`, `C`, `D`; NA for a parameter the model does
    not have.
    """
    tables = []
    for domain, x in read_domains(path).items():
        fits = mask_absent_parameters(fit_domain(x))
        fits.insert(0, "domain", domain)
        tables.append(fits)
    return pd.concat(tables, ignore_index=True)
