# -----------------------------------------------------------------------------
# File:    python/zinull/em.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Tables 1 to 3, Web Tables S2 to S16.
# Purpose: EM algorithms for the null distribution fitted to the counts at or
#          below a cutoff C (the counts above C are treated as missing), for
#          four null models: ZIGP, ZIP, generalized Poisson (GP), Poisson (P).
#          The Python counterpart of R/em.R, with the same results.
# Inputs:  None.
# Outputs: None.
# Usage:   from zinull.em import em_estimates
#          em_estimates(x)["ZIGP"]
# Source:  R/em.R.
# Notes:   The estimates are arrays (eta, lambda, theta). The EM steps are
#          compiled (numba); they add every sum in order and round every
#          product and sum on its own, as R does.
# -----------------------------------------------------------------------------
import math

import numpy as np
from numba import njit

from . import rbase
from .rmath import sqrt
from .zigp import (SMALL, _gp_density, _gp_pmf, _gp_pmf_times_x, _r_max, _zigp_density,
                   _zigp_mass_below)

# Null models in the order of the tables.
MODELS = ("ZIGP", "ZIP", "GP", "P")

# The EM algorithm stops when no parameter moves by more than this.
EM_TOLERANCE = 10.0 ** -4

_ZIGP, _ZIP, _GP, _P = range(4)  # the models in the compiled functions


def _code(model):
    """Position of a model in `MODELS`."""
    if model not in MODELS:
        raise ValueError(f"model must be one of {MODELS}")
    return MODELS.index(model)


def count_up_to(x, C):
    """Counts of 0, 1, ..., C.

    Parameters: x (counts), C (largest value tabulated).
    Returns: float array of length C + 1.
    """
    return rbase.count_frequencies_from_zero(x, int(C))


def _mean_of(x):
    """`mean(x)` as R computes it for the type of x: whole numbers stored as
    integers are added and divided once, doubles get R's second pass."""
    x = np.asarray(x)
    return np.float64(rbase.r_mean_integer(x) if x.dtype.kind in "iub" else rbase.r_mean(x))


def _variance(x):
    """`var(x)` as a numpy double, so that a division by it follows IEEE 754."""
    return np.float64(rbase.r_var(x))


def _sqrt(x):
    """`sqrt(x)` as a numpy double: NaN below zero."""
    return np.float64(sqrt(x))


def em_start(x, model):
    """Starting values of the EM algorithm (moment estimates).

    Parameters: x (counts; an integer array is averaged as R averages
    integers, a float array as R averages doubles), model (one of `MODELS`).
    Returns: array (eta, lambda, theta).
    """
    _code(model)
    x = np.asarray(x)
    zero = np.float64(rbase.r_mean_integer(x == 0))
    nonzero = x[x != 0]
    with np.errstate(all="ignore"):  # a division by zero gives Inf or NaN, as in R
        if model == "ZIGP":
            lower = x[x <= rbase.quantile(x, 0.75)]
            ex = _mean_of(lower)
            phi = _sqrt(_variance(lower) / ex - zero / (1 - zero) * ex + SMALL)
            start = (rbase.r_max(0, zero), rbase.r_max(0, ex / ((1 - zero) * phi)),
                     rbase.r_min(rbase.r_max(0, 1 - 1 / phi), 1))
        elif model == "ZIP":
            start = (rbase.r_max(0, zero), rbase.r_max(0, _mean_of(x) / (1 - zero)), 0.0)
        elif model == "GP":
            lower = nonzero[nonzero <= rbase.quantile(nonzero, 0.5)]
            phi = _sqrt(_variance(lower) / (_mean_of(lower) + SMALL))
            start = (0.0, rbase.r_max(0, _mean_of(x) / phi), rbase.r_max(0, 1 - 1 / phi))
        else:
            start = (0.0, rbase.r_max(0, _mean_of(nonzero)), 0.0)
    return np.array(start, dtype=np.float64)


@njit(cache=True, error_model="numpy")
def _em_clamp(estimate, model):
    eta = _r_max(0.0, estimate[0]) if model == _ZIGP or model == _ZIP else 0.0
    if model == _ZIGP:
        theta = _r_max(SMALL, estimate[2])
    elif model == _GP:
        theta = _r_max(0.0, estimate[2])
    else:
        theta = 0.0
    return np.array([eta, _r_max(SMALL, estimate[1]), theta])


def em_clamp(estimate, model):
    """Keep the current estimates inside the parameter space.

    Parameters: estimate (eta, lambda, theta), model (one of `MODELS`).
    Returns: array (eta, lambda, theta) used in the next EM step.
    """
    return _em_clamp(np.asarray(estimate, dtype=np.float64), _code(model))


@njit(cache=True, error_model="numpy")
def _em_sums(hold, count, K, dispersion):
    """The sums of `em_sums()`: zero, gp, above, total_above, total, lambda,
    theta (the last two 0 without `dispersion`)."""
    eta, lam, theta = hold[0], hold[1], hold[2]
    C = len(count) - 1
    n = 0.0
    for v in range(C + 1):
        n += count[v]
    # the counts above C: C + 1, ..., K, or K + 1, K when C = K (R's (C + 1):K)
    n_missing = K - C if C < K else 2
    missing = np.empty(n_missing)
    for i in range(n_missing):
        missing[i] = C + 1.0 + i if C < K else K + 1.0 - i

    mass_below = _zigp_mass_below(eta, lam, theta, C)
    mixture_zero = _zigp_density(0.0, eta, lam, theta)
    gp = 0.0
    pmf_below = 0.0
    total_below = 0.0
    lagged_lambda = 0.0
    lagged_theta = 0.0
    for v in range(C + 1):
        x = float(v)
        from_gp = (1.0 - eta) * _gp_density(x, lam, theta) / _zigp_density(x, eta, lam, theta)
        gp += count[v] * from_gp
        pmf_below += _gp_pmf(x, lam, theta)
        total_below += x * count[v] * from_gp
        if dispersion:
            lagged = max(x - 1.0, 0.0) * count[v]
            lagged_lambda += lagged * (lam / (lam + theta * x))
            lagged_theta += lagged * (theta * x / (lam + theta * x))
    zero = count[0] * (eta / mixture_zero)
    above = n * ((1.0 - eta) * (1.0 - pmf_below)) / mass_below

    sum_x = 0.0        # sums over the counts above C
    sum_x_lambda = 0.0
    sum_lambda = 0.0
    sum_x_theta = 0.0
    sum_theta = 0.0
    for i in range(n_missing):
        x = missing[i]
        pmf_x = _gp_pmf_times_x(x, lam, theta)
        sum_x += pmf_x
        if dispersion:
            pmf = _gp_pmf(x, lam, theta)
            to_lambda = lam / (lam + theta * x)
            to_theta = theta * x / (lam + theta * x)
            sum_x_lambda += pmf_x * to_lambda
            sum_lambda += pmf * to_lambda
            sum_x_theta += pmf_x * to_theta
            sum_theta += pmf * to_theta
    total_above = n * ((1.0 - eta) * sum_x) / mass_below
    total = total_below + total_above
    sum_for_lambda = 0.0
    sum_for_theta = 0.0
    if dispersion:
        sum_for_lambda = (gp + above + lagged_lambda
                          + n * ((1.0 - eta) * sum_x_lambda) / mass_below
                          - n * ((1.0 - eta) * sum_lambda) / mass_below)
        sum_for_theta = (lagged_theta + n * ((1.0 - eta) * sum_x_theta) / mass_below
                         - n * ((1.0 - eta) * sum_theta) / mass_below)
    return zero, gp, above, total_above, total, sum_for_lambda, sum_for_theta


def em_sums(hold, count, K, dispersion):
    """Expected sums of the E-step for the counts at or below C and above C.

    The counts above C are replaced by their expectations under the current
    null, n P(X = v) / P(X <= C) for v = C + 1, ..., K. For C = K the range
    covers K + 1 and K, as R's (C + 1):K.

    Parameters: hold (current eta, lambda, theta), count (counts of
    0, ..., C), K (largest count of the data), dispersion (if True also the
    sums of the lambda and theta updates of the generalized Poisson models).
    Returns: dict: `zero` (expected structural zeros), `gp` and `above`
    (expected numbers from the Poisson part at or below and above C),
    `total` (their expected sum of counts), `total_above` (its part above
    C), and with `dispersion` the numerators `lambda` and `theta` of the
    updates of the generalized Poisson models.
    """
    values = _em_sums(np.asarray(hold, dtype=np.float64), np.asarray(count, dtype=np.float64),
                      int(K), bool(dispersion))
    names = ("zero", "gp", "above", "total_above", "total", "lambda", "theta")
    return dict(zip(names[:7 if dispersion else 5], (float(value) for value in values)))


@njit(cache=True, error_model="numpy")
def _em_step(hold, count, K, model):
    zero, gp, above, total_above, total, for_lambda, for_theta = _em_sums(
        hold, count, K, model == _ZIGP or model == _GP)
    eta = zero / (zero + gp + above + SMALL)
    if model == _ZIGP:
        return np.array([eta, (for_lambda + SMALL) / (gp + above + SMALL),
                         (for_theta + SMALL) / (total + SMALL)])
    if model == _ZIP:
        return np.array([eta, total / (gp + above + SMALL), 0.0])
    if model == _GP:
        return np.array([0.0, (for_lambda + SMALL) / (gp + above + SMALL),
                         (for_theta + SMALL) / (total + SMALL)])
    n = 0.0
    weighted = 0.0
    for v in range(len(count)):
        n += count[v]
        weighted += v * count[v]
    return np.array([0.0, (weighted + total_above) / (n + above + SMALL), 0.0])


def em_step(hold, count, K, model):
    """One EM step.

    Parameters: hold (current eta, lambda, theta, inside the parameter
    space), count (counts of 0, ..., C), K (largest count of the data),
    model (one of `MODELS`).
    Returns: array of the updated (eta, lambda, theta).
    """
    return _em_step(np.asarray(hold, dtype=np.float64), np.asarray(count, dtype=np.float64),
                    int(K), _code(model))


@njit(cache=True, error_model="numpy")
def _em_fit(count, K, model, start):
    """EM iterations until no parameter moves by more than the tolerance.
    Returns the estimates and whether a parameter became NaN (R stops with
    an error there)."""
    estimate = start
    while True:
        hold = _em_clamp(estimate, model)
        estimate = _em_step(hold, count, K, model)
        moved = False
        undefined = False
        for i in range(3):
            change = abs(hold[i] - estimate[i])
            if change > EM_TOLERANCE:
                moved = True
            elif math.isnan(change):
                undefined = True
        if not moved:
            return estimate, undefined


def em_fit(x, C, model, start):
    """EM estimates of one null model for one cutoff.

    Parameters: x (counts), C (cutoff: the counts at or below C are used),
    model (one of `MODELS`), start (starting values, from `em_start()`).
    Returns: array (eta, lambda, theta); 0 for a parameter the model does
    not have. Raises ValueError when a parameter becomes NaN, where R stops
    with an error.
    """
    estimate, undefined = _em_fit(count_up_to(x, C), int(np.max(x)), _code(model),
                                  np.asarray(start, dtype=np.float64))
    if undefined:
        raise ValueError(f"the EM algorithm of {model} gave a missing value at C = {C}")
    return estimate


@njit(cache=True, error_model="numpy")
def _em_fits(count, model, start):
    """`_em_fit()` for every cutoff C = 1, ..., K; `count` holds the counts
    of 0, ..., K. Returns the 3 x K estimates and the first cutoff at which
    a parameter became NaN (0 when none)."""
    K = len(count) - 1
    estimates = np.empty((3, K))
    for C in range(1, K + 1):
        estimate, undefined = _em_fit(count[:C + 1], K, model, start)
        if undefined:
            return estimates, C
        estimates[:, C - 1] = estimate
    return estimates, 0


def em_estimates(x):
    """EM estimates of the four null models for every candidate cutoff.

    Parameters: x (counts).
    Returns: dict (one element per model) of 3 x max(x) arrays with rows
    eta, lambda, theta and one column per cutoff C = 1, ..., max(x).
    """
    count = count_up_to(x, int(np.max(x)))
    estimates = {}
    for model in MODELS:
        estimates[model], failed = _em_fits(count, _code(model), em_start(x, model))
        if failed:
            raise ValueError(f"the EM algorithm of {model} gave a missing value at C = {failed}")
    return estimates
