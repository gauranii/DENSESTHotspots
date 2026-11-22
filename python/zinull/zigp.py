# -----------------------------------------------------------------------------
# File:    python/zinull/zigp.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. All simulated data.
# Purpose: Zero-inflated generalized Poisson (ZIGP) distribution: random
#          numbers (inversion method), probability functions, the mass of
#          the null below a cutoff and the bound D of the two-stage procedure.
#          The Python counterpart of R/zigp.R, with the same results.
# Inputs:  None.
# Outputs: None.
# Usage:   from zinull.zigp import generate_zigp
#          generate_zigp(10, mu=1.5, phi=1, omega=0.8, seed=1)
# Source:  R/zigp.R.
# Notes:   `lam` stands for the parameter lambda (a reserved word in Python).
#          The functions whose names start with an underscore are compiled
#          (numba) and work on one number; every sum is added in order, and
#          every product and sum is rounded on its own, as in R.
# -----------------------------------------------------------------------------
import math

import numpy as np
from numba import njit

from . import rng
from .rmath import _exp, _lgammafn, _log, _r_pow

# Constant added to probabilities and denominators to keep them positive.
SMALL = 10.0 ** -10


# ----- Random numbers -----

@njit(cache=True, error_model="numpy")
def _zigp_probability(i, mu, phi, omega):
    """P(X = i) of the ZIGP distribution for a count i >= 1."""
    return ((1.0 - omega) * mu * _r_pow(mu + (phi - 1.0) * i, i - 1.0) / _exp(_lgammafn(i + 1.0))
            * _r_pow(phi, -i) * _exp(-1.0 / phi * (mu + (phi - 1.0) * i)))


@njit(cache=True, error_model="numpy")
def _generate_zigp(state, n, mu, phi, omega):
    """`n` ZIGP counts by inversion from the stream `state`."""
    x = np.zeros(n)
    u = np.empty(n)
    upper = -np.inf
    for j in range(n):
        u[j] = rng._unif(state)
        if u[j] > upper:
            upper = u[j]
    s = np.zeros(1000)  # cumulative probabilities of 0, 1, 2, ...
    p = omega + (1.0 - omega) * _exp(-mu / phi)  # P(X = 0)
    s[0] = p
    if upper > 0.0:
        recursive = False
        log_p_previous = 0.0
        i = 1.0
        while s[int(i) - 1] < upper:
            if not recursive:
                p = _zigp_probability(i, mu, phi, omega)
            if p == np.inf:  # switch to the recursion on the log scale
                recursive = True
                log_p_previous = _log(_zigp_probability(i - 1.0, mu, phi, omega))
            if recursive:
                log_p = _log((mu + (i - 1.0) * (phi - 1.0)) / (phi * i)
                             * _r_pow(1.0 + (phi - 1.0) / (mu + (i - 1.0) * (phi - 1.0)), i - 1.0)
                             * _exp(1.0 / phi - 1.0)) + log_p_previous
                log_p_previous = log_p
                p = _exp(log_p)
            if int(i) + 1 > len(s):
                s = np.concatenate((s, np.zeros(1000)))
            s[int(i)] = s[int(i) - 1] + p
            i += 1.0
    for j in range(n):
        i = 0
        while u[j] > s[i]:
            i += 1
        x[j] = i
    return x


def generate_zigp(n, mu, phi, omega, seed=None):
    """Random numbers from the ZIGP distribution by inversion.

    The cumulative probabilities are built up to the largest uniform drawn,
    and each uniform is mapped to the first count whose cumulative
    probability reaches it. Gives the same numbers as `generate_zigp()` of
    R/zigp.R for the same seed.

    Parameters: n (number of draws), mu (mean parameter of the generalized
    Poisson part, lambda / (1 - theta)), phi (dispersion parameter,
    1 / (1 - theta)), omega (zero-inflation probability eta, in [0, 1]),
    seed (when not None, the seed set right before the draws; otherwise the
    draws continue the stream).
    Returns: float array of `n` counts.
    """
    if omega < 0 or omega > 1:
        raise ValueError("omega has to be in [0, 1]")
    if seed is not None:
        rng.set_seed(seed)
    return _generate_zigp(rng._state, int(n), float(mu), float(phi), float(omega))


# ----- Probability functions -----

@njit(cache=True, error_model="numpy")
def _gp_pmf(x, lam, theta):
    return _exp(_log(lam) + (x - 1.0) * _log(lam + theta * x) - _lgammafn(x + 1.0) - lam
                - theta * x)


@njit(cache=True, error_model="numpy")
def _gp_pmf_times_x(x, lam, theta):
    return _exp(_log(x * lam) + (x - 1.0) * _log(lam + theta * x) - _lgammafn(x + 1.0) - lam
                - theta * x)


@njit(cache=True, error_model="numpy")
def _gp_density(x, lam, theta):
    return _gp_pmf(x, lam, theta) + SMALL


@njit(cache=True, error_model="numpy")
def _zigp_density(x, eta, lam, theta):
    zero = 1.0 if x == 0.0 else 0.0
    return eta * zero + (1.0 - eta) * _gp_density(x, lam, theta) + SMALL


@njit(cache=True, error_model="numpy")
def _zigp_mass_below(eta, lam, theta, C):
    total = 0.0
    for x in range(C + 1):
        total += _gp_pmf(float(x), lam, theta)
    return eta + (1.0 - eta) * total


@njit(cache=True, error_model="numpy")
def _over_counts(kind, x, eta, lam, theta):
    """One of the probability functions at every count of an array."""
    out = np.empty(len(x))
    for i in range(len(x)):
        if kind == 0:
            out[i] = _gp_pmf(x[i], lam, theta)
        elif kind == 1:
            out[i] = _gp_pmf_times_x(x[i], lam, theta)
        elif kind == 2:
            out[i] = _gp_density(x[i], lam, theta)
        else:
            out[i] = _zigp_density(x[i], eta, lam, theta)
    return out


def _vectorized(kind, x, eta, lam, theta):
    """A probability function at one count (float) or at each of an array."""
    values = _over_counts(kind, np.atleast_1d(np.asarray(x, dtype=np.float64)), float(eta),
                          float(lam), float(theta))
    return float(values[0]) if np.ndim(x) == 0 else values


def gp_pmf(x, lam, theta):
    """Probability function of the generalized Poisson distribution.

    Parameters: x (count or array of counts), lam, theta (parameters;
    theta = 0 gives the Poisson distribution).
    Returns: P(X = x), a float or one value per element of `x`.
    """
    return _vectorized(0, x, 0.0, lam, theta)


def gp_pmf_times_x(x, lam, theta):
    """x P(X = x) for the generalized Poisson distribution, on the log scale.

    Parameters: x (positive count or array of counts), lam, theta.
    Returns: x times `gp_pmf(x, lam, theta)`.
    """
    return _vectorized(1, x, 0.0, lam, theta)


def gp_density(x, lam, theta):
    """Generalized Poisson probabilities kept positive.

    Parameters: x (count or array of counts), lam, theta.
    Returns: `gp_pmf()` plus `SMALL`.
    """
    return _vectorized(2, x, 0.0, lam, theta)


def zigp_density(x, eta, lam, theta):
    """Probability function of the ZIGP distribution, kept positive.

    Parameters: x (count or array of counts), eta (zero-inflation
    probability), lam, theta (parameters of the generalized Poisson part).
    Returns: eta 1(x = 0) + (1 - eta) `gp_density()`, plus `SMALL`.
    """
    return _vectorized(3, x, eta, lam, theta)


def zigp_mass_below(eta, lam, theta, C):
    """Probability that a ZIGP count is at most C.

    Parameters: eta, lam, theta (parameters of the ZIGP distribution), C
    (cutoff, a whole number >= 0).
    Returns: P(X <= C).
    """
    return float(_zigp_mass_below(float(eta), float(lam), float(theta), int(C)))


@njit(cache=True, error_model="numpy")
def _r_max(a, b):
    """`max(a, b)` of R: NaN when one of the two is."""
    if math.isnan(a) or math.isnan(b):
        return np.nan
    return a if a > b else b


@njit(cache=True, error_model="numpy")
def _gp_bound(lam, theta, n):
    return _r_max(lam / (_exp(theta - 1.0) - theta), _log(n) / (theta - 1.0 - _log(theta)))


def gp_bound(lam, theta, n):
    """Bound D_n above which a generalized Poisson null gives no count among n.

    Parameters: lam, theta (parameters of the generalized Poisson
    distribution), n (number of positions).
    Returns: the bound (not rounded).
    """
    return float(_gp_bound(float(lam), float(theta), float(n)))
