# -----------------------------------------------------------------------------
# File:    python/zinull/cutoff.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Tables 1 to 3, Web Tables S2 to S16.
# Purpose: Choice of the cutoff C below which the counts follow the null, by
#          the two criteria C1 and C2 (truncated null likelihood plus a
#          penalty for the counts above C), with the estimate of the null
#          proportion pi0 and the bound D of the two-stage procedure.
#          The Python counterpart of R/cutoff.R, with the same results.
# Inputs:  None.
# Outputs: None.
# Usage:   from zinull.cutoff import select_null_fits
#          select_null_fits(x)
# Source:  R/cutoff.R.
# Notes:   A table of null fits is a pandas DataFrame with the columns
#          cutoff, model, eta, lambda, theta, pi0, C, D; one fit is a row of
#          it (or any mapping with these names).
# -----------------------------------------------------------------------------
import math

import numpy as np
import pandas as pd
from numba import njit

from . import rbase
from .em import MODELS, count_up_to, em_estimates
from .rmath import _log, log
from .zigp import SMALL, _gp_pmf, _zigp_mass_below, gp_bound

# Cutoff criteria in the order of the tables.
CUTOFFS = ("C1", "C2")


@njit(cache=True, error_model="numpy")
def _truncated_neg_loglik(count, eta, lam, theta, C):
    zero = count[0] * (-_log(eta + (1.0 - eta) * math.exp(-lam)))
    positive = 0.0
    for x in range(1, C + 1):
        positive += count[x] * (-_log((1.0 - eta) * _gp_pmf(float(x), lam, theta)))
    return zero + positive


def truncated_neg_loglik(count, eta, lam, theta, C):
    """Negative log likelihood of the counts 0, ..., C under the null.

    Parameters: count (counts of 0, ..., max(x)), eta, lam, theta (null
    parameters), C (cutoff).
    Returns: the negative log likelihood.
    """
    return float(_truncated_neg_loglik(np.asarray(count, dtype=np.float64), float(eta),
                                       float(lam), float(theta), int(C)))


def cutoff_penalties(count):
    """Penalties of the two criteria for every candidate cutoff.

    C1 lets every count above the cutoff have its own probability; C2 uses
    only the numbers of counts below and above the cutoff. At the last
    cutoff the C1 penalty uses the frequency of max(x) - 1.

    Parameters: count (counts of 0, ..., K = max(x)).
    Returns: dict with `C1` and `C2`, each an array of length K (cutoffs
    1, ..., K).
    """
    count = np.asarray(count, dtype=np.float64)
    K = len(count) - 1
    M = rbase.r_sum(count)
    below = rbase.r_cumsum(count)[1:]
    share = rbase.r_cumsum(count + SMALL)[1:] / M
    information = count * log(count + SMALL)
    above = np.array([rbase.r_sum(information[C + 1:]) for C in range(1, K + 1)])
    above[K - 1] = information[K - 1]
    return {"C1": log(M) * (M - below) - above,
            "C2": below * (-log(share)) + (M - below) * (-log(1 + SMALL - share))}


def choose_cutoff(criterion):
    """Cutoff chosen by a criterion: its first local minimum.

    Parameters: criterion (values of the criterion for the cutoffs 1, 2, ...).
    Returns: the chosen cutoff (a whole number, 1 for the first value).
    """
    criterion = np.asarray(criterion, dtype=np.float64)
    with np.errstate(invalid="ignore"):
        turns = np.flatnonzero(np.diff(np.sign(np.diff(criterion))) > 0)
    candidates = [int(turn) + 3 for turn in turns]  # which(...) + 2, one-based
    smallest = rbase.which_min(criterion)
    if smallest is not None:
        candidates.append(smallest + 1)
    if not candidates:
        raise ValueError("the criterion has no value")
    return min(candidates)


def rejection_bound(model, C, lam, theta, n, cap=math.inf):
    """Bound D of the two-stage procedure: counts at or above D are rejected.

    Parameters: model (one of `MODELS`), C (chosen cutoff), lam, theta (null
    estimates at the cutoff), n (number of positions), cap (largest value
    allowed: max(x) for the protein domain data).
    Returns: the bound, a whole number (float).
    """
    if model in ("ZIGP", "GP"):
        null_bound = gp_bound(lam, theta, n)
    else:
        null_bound = rbase.r_max(lam, log(n))
    return rbase.r_min(cap, float(np.ceil(rbase.r_max(C + 1, null_bound))))


@njit(cache=True, error_model="numpy")
def _mass_and_loss(estimates, count):
    """P(X <= C) and the truncated negative log likelihood at the estimates
    of every cutoff C = 1, ..., K."""
    K = estimates.shape[1]
    mass = np.empty(K)
    loss = np.empty(K)
    for C in range(1, K + 1):
        eta, lam, theta = estimates[0, C - 1], estimates[1, C - 1], estimates[2, C - 1]
        mass[C - 1] = _zigp_mass_below(eta, lam, theta, C)
        loss[C - 1] = _truncated_neg_loglik(count, eta, lam, theta, C)
    return mass, loss


def select_model_fits(model, estimates, count, penalties, cap=math.inf):
    """Null fit of one model at the cutoff chosen by each criterion.

    Parameters: model (one of `MODELS`), estimates (3 x K array of EM
    estimates, from `em_estimates()`), count (counts of 0, ..., K),
    penalties (output of `cutoff_penalties()`), cap (as for
    `rejection_bound()`).
    Returns: DataFrame with one row per criterion and the columns `cutoff`,
    `model`, `eta`, `lambda`, `theta`, `pi0`, `C`, `D`.
    """
    count = np.asarray(count, dtype=np.float64)
    estimates = np.ascontiguousarray(estimates, dtype=np.float64)
    M = rbase.r_sum(count)
    mass, loss = _mass_and_loss(estimates, count)
    with np.errstate(all="ignore"):
        pi0 = rbase.pmin(1, rbase.r_cumsum(count)[1:] / (M * mass))
    rows = []
    for cutoff in CUTOFFS:
        C = choose_cutoff(penalties[cutoff] + loss)
        eta, lam, theta = (float(value) for value in estimates[:, C - 1])
        rows.append({"cutoff": cutoff, "model": model, "eta": eta, "lambda": lam, "theta": theta,
                     "pi0": float(pi0[C - 1]), "C": float(C),
                     "D": float(rejection_bound(model, C, lam, theta, M, cap))})
    return pd.DataFrame(rows)


def select_null_fits(x, cap=math.inf):
    """Null fits of the four models at the cutoffs chosen by the two criteria.

    Parameters: x (counts), cap (as for `rejection_bound()`).
    Returns: DataFrame with one row per criterion and model (criterion
    first) and the columns of `select_model_fits()`.
    """
    count = count_up_to(x, int(np.max(x)))
    penalties = cutoff_penalties(count)
    estimates = em_estimates(x)
    fits = pd.concat([select_model_fits(model, estimates[model], count, penalties, cap)
                      for model in MODELS], ignore_index=True)
    order = np.argsort([CUTOFFS.index(cutoff) for cutoff in fits["cutoff"]], kind="stable")
    return fits.iloc[order].reset_index(drop=True)


def mask_absent_parameters(fits):
    """Replace by NA the parameters a model does not have.

    Parameters: fits (DataFrame with the columns `model`, `eta`, `theta`).
    Returns: a copy of `fits` with `eta` NA for GP and P, `theta` NA for
    ZIP and P.
    """
    fits = fits.copy()
    for name, absent in (("eta", ("GP", "P")), ("theta", ("ZIP", "P"))):
        values = fits[name].to_numpy(dtype=np.float64, copy=True)
        values[fits["model"].isin(absent).to_numpy()] = rbase.NA
        fits[name] = values
    return fits
