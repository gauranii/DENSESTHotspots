# -----------------------------------------------------------------------------
# File:    python/zinull/fdr.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Tables 1 to 3, Web Tables S6, S7, S12 to S15.
# Purpose: The three multiple testing procedures: one-stage (local FDR),
#          two-stage (local FDR below the bound D, rejection at or above D)
#          and Storey's FDR (modified Benjamini-Hochberg with the estimated
#          pi0), for data with known labels (simulation) and without (real).
#          The Python counterpart of R/fdr.R, with the same results.
# Inputs:  None.
# Outputs: None.
# Usage:   from zinull.cutoff import select_null_fits
#          from zinull.fdr import one_stage_count
#          fit = select_null_fits(x).iloc[0]
#          one_stage_count(x, fit, level=0.10)
# Source:  R/fdr.R.
# Notes:   Positions returned by `storey_rejected()` are zero-based.
# -----------------------------------------------------------------------------
import math

import numpy as np

from . import rbase
from .em import count_up_to
from .zigp import gp_pmf, zigp_density


def _whole(x):
    """Counts as array indices."""
    return np.asarray(x).astype(np.int64)


def inverse_frequency(x, all):
    """Reciprocal of the empirical probability of each count.

    Parameters: x (counts of a subset of the data), all (counts of the whole
    data set).
    Returns: len(all) / (number of elements of `all` equal to x), per element.
    """
    all = _whole(all)
    return len(all) / np.bincount(all).astype(np.float64)[_whole(x)]


def local_fdr(x, weight, fit):
    """Local FDR: estimated null probability over empirical probability.

    Parameters: x (counts), weight (their `inverse_frequency()`), fit (null
    fit: mapping or row with `eta`, `lambda`, `theta`, `pi0`).
    Returns: array of the local FDR of each count.
    """
    x = np.asarray(x, dtype=np.float64)
    if len(x) == 0:
        return np.empty(0)
    return fit["pi0"] * zigp_density(x, fit["eta"], fit["lambda"], fit["theta"]) * weight


def count_below(values, level):
    """Number of values below a level, ignoring missing values.

    Parameters: values (numeric vector), level (threshold).
    Returns: the count (int).
    """
    values = np.asarray(values, dtype=np.float64)
    return int(np.count_nonzero(values[~np.isnan(values)] < level))


def count_between(x, lfdr, fit, level):
    """Number of local FDR values below the level among the counts above C
    and below D.

    Parameters: x (counts), lfdr (their local FDR), fit (null fit with the
    cutoff `C` and the bound `D`), level (FDR level).
    Returns: the difference of the two counts (int).
    """
    x = np.asarray(x, dtype=np.float64)
    lfdr = np.asarray(lfdr, dtype=np.float64)
    return count_below(lfdr[x < fit["D"]], level) - count_below(lfdr[x <= fit["C"]], level)


# ----- Data without labels -----

def one_stage_count(x, fit, level):
    """Rejections of the one-stage procedure: local FDR below the level.

    Parameters: x (counts), fit (null fit, a row of `select_null_fits()`),
    level (FDR level).
    Returns: number of rejected positions.
    """
    return count_below(local_fdr(x, inverse_frequency(x, x), fit), level)


def two_stage_count(x, fit, level):
    """Rejections of the two-stage procedure.

    Positions with a count at or above D are rejected; of those below D and
    above the cutoff C, the ones with a local FDR below the level.

    Parameters: x (counts), fit (null fit, a row of `select_null_fits()`),
    level (FDR level).
    Returns: number of rejected positions.
    """
    lfdr = local_fdr(x, inverse_frequency(x, x), fit)
    at_or_above = int(np.count_nonzero(np.asarray(x, dtype=np.float64) > fit["D"] - 1))
    return count_between(x, lfdr, fit, level) + at_or_above


def storey_count(x, fit, level):
    """Rejections of Storey's FDR procedure.

    The frequencies are read from the vector of the frequencies of
    0, ..., max(x), at the positions of the rejected counts among the
    distinct observed counts.

    Parameters: x (counts), fit (null fit, a row of `select_null_fits()`),
    level (FDR level).
    Returns: number of rejected positions.
    """
    K = int(np.max(x))
    count = count_up_to(x, K)
    observed = np.flatnonzero(count > 0)
    p_values = storey_p_values(K, fit)[observed]
    rejected = np.sort(storey_rejected(p_values, count[observed], fit["pi0"], level))
    return int(rbase.r_sum(count[rejected]))


# ----- Storey's FDR -----

def storey_p_values(K, fit):
    """P-values of the counts 0, ..., K under the null: P(X >= x).

    Parameters: K (largest count), fit (null fit with `eta`, `lambda`,
    `theta`).
    Returns: array of K + 1 p-values; the p-value of 0 is 1.
    """
    below = rbase.r_cumsum(gp_pmf(rbase.colon(0, K - 1), fit["lambda"], fit["theta"]))
    return np.concatenate(([1.0], rbase.pmax(0, (1 - fit["eta"]) * (1 - below))))


def storey_rejected(p_values, frequency, pi0, level):
    """Distinct counts rejected by the modified Benjamini-Hochberg procedure.

    The i-th smallest p-value is compared with the level times the share of
    positions whose count is at least the i-th smallest distinct count,
    over pi0.

    Parameters: p_values (p-values of the distinct counts, in increasing
    order of the counts), frequency (frequencies of the distinct counts),
    pi0 (estimated proportion of nulls), level (FDR level).
    Returns: int array of the positions in `p_values` (zero-based) of the
    rejected counts.
    """
    p_values = np.asarray(p_values, dtype=np.float64)
    frequency = np.asarray(frequency, dtype=np.float64)
    M = rbase.r_sum(frequency)
    at_least = np.concatenate(([M], M - rbase.r_cumsum(frequency)[:-1]))
    ordered = rbase.r_order(p_values)
    with np.errstate(all="ignore"):
        passed = np.flatnonzero(p_values[ordered] <= at_least * (level / np.float64(M * pi0)))
    if len(passed) == 0:
        return np.empty(0, dtype=np.int64)
    return ordered[:passed[-1] + 1]


# ----- Data with labels (simulation) -----

def discovery_rates(false, true, n_nonnull, R=None):
    """Number of rejections, FDR and TPR.

    The FDR of a replication without rejections is NaN, so that the mean
    over replications is taken over the replications with rejections.

    Parameters: false, true (numbers of false and true discoveries),
    n_nonnull (number of non-null positions), R (number of rejections; by
    default `true + false`).
    Returns: array (R, FDR, TPR).
    """
    if R is None:
        R = true + false
    tpr = true / n_nonnull if n_nonnull != 0 else (math.nan if true == 0 else math.inf)
    return np.array([R, math.nan if R == 0 else false / R, tpr], dtype=np.float64)


RATE_NAMES = ("R", "FDR", "TPR")


def local_fdr_rates(null, nonnull, fit, level):
    """One-stage and two-stage procedures on data with labels.

    The one-stage R is the number of positions whose local FDR is below the
    level; its false discoveries are counted above the cutoff C only, and
    its true discoveries (for the TPR) use, for the non-null counts at or
    above D, the empirical probabilities of the smallest non-null counts.

    Parameters: null, nonnull (sorted counts of the null and non-null
    positions), fit (null fit, a row of `select_null_fits()`), level (FDR
    level).
    Returns: dict: R, FDR, TPR of the two-stage, then of the one-stage
    procedure (`two_stage_R`, ..., `one_stage_TPR`).
    """
    null = np.asarray(null, dtype=np.float64)
    nonnull = np.asarray(nonnull, dtype=np.float64)
    all = np.concatenate((null, nonnull))
    lfdr_null = local_fdr(null, inverse_frequency(null, all), fit)
    weight = inverse_frequency(nonnull, all)
    lower = nonnull < fit["D"]
    upper = nonnull[nonnull > fit["D"] - 1]
    true_lower = count_below(local_fdr(nonnull[lower], weight[lower], fit), level)
    true_upper = count_below(local_fdr(upper, weight[:len(upper)], fit), level)

    false_two = (max(0, count_between(null, lfdr_null, fit, level))
                 + len(null) - int(np.count_nonzero(null < fit["D"])))
    true_two = true_lower + len(nonnull) - int(np.count_nonzero(lower))
    false_one = count_below(lfdr_null[null > fit["C"]], level)
    R_one = count_below(local_fdr(all, inverse_frequency(all, all), fit), level)
    rates = np.concatenate((discovery_rates(false_two, true_two, len(nonnull)),
                            discovery_rates(false_one, true_lower + true_upper, len(nonnull),
                                            R=R_one)))
    names = [stage + name for stage in ("two_stage_", "one_stage_") for name in RATE_NAMES]
    return dict(zip(names, (float(rate) for rate in rates)))


def storey_rates(null, nonnull, fit, level):
    """Storey's FDR procedure on data with labels.

    A single rejected count is not counted, and the FDR is NaN when nothing
    is rejected.

    Parameters: null, nonnull (counts of the null and non-null positions),
    fit (null fit, a row of `select_null_fits()`), level (FDR level).
    Returns: dict: `storey_R`, `storey_FDR`, `storey_TPR`.
    """
    null, nonnull = _whole(null), _whole(nonnull)
    all = np.concatenate((null, nonnull))
    distinct = np.unique(all)
    frequency = np.bincount(all)[distinct].astype(np.float64)
    p_values = storey_p_values(int(np.max(all)), fit)[distinct]
    rejected = storey_rejected(p_values, frequency, fit["pi0"], level)
    if len(rejected) < 2:
        rejected = np.empty(0, dtype=np.int64)
    size = int(np.max(all)) + 1
    false = int(np.bincount(null, minlength=size)[distinct][rejected].sum())
    true = int(np.bincount(nonnull, minlength=size)[distinct][rejected].sum())
    rates = discovery_rates(false, true, len(nonnull))
    return dict(zip(("storey_" + name for name in RATE_NAMES), (float(rate) for rate in rates)))
