# -----------------------------------------------------------------------------
# File:    python/zinull/rbase.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. All tables and figures.
# Purpose: What base R, stats and parallel supply to the code in R/, with R's
#          order of operations: sums, cumulative sums, means and variances
#          that add in R's order, quantile(), min() and max() with R's
#          treatment of missing values, tabulate(), order(), which.min(),
#          the sequence a:b, write.csv() with numbers written as R writes
#          them, and mclapply().
# Inputs:  None.
#          R adds its sums and means in the C type long double, which is a
#          double on arm64 macOS; the sums here are added in doubles, so they
#          round as R's do there.
# Outputs: None.
# Usage:   from zinull.rbase import r_sum, r_mean, r_var, quantile, write_csv
# Source:  R 4.6.1: src/main/summary.c (sum, mean, min, max), cum.c (cumsum),
#          src/library/stats/src/cov.c (var), stats::quantile (type 7),
#          src/main/format.c (formatReal, scientific), utils::write.csv,
#          parallel::mclapply.
# -----------------------------------------------------------------------------
import math
import multiprocessing
import struct

import numpy as np
import pandas as pd
from numba import njit

from .rmath import muladd

# ----- Missing values -----

# R's missing value of a double: a NaN whose low word is 1954. Kept apart
# from the NaN of an undefined result because sprintf() of R writes the two differently
# ("NA" and "NaN").
NA = struct.unpack("<d", struct.pack("<Q", 0x7FF80000000007A2))[0]


def is_na(x):
    """`is.na(x) & !is.nan(x)`: which values are R's missing value NA
    (and not another NaN).

    Parameters: x (float or array). Returns: bool, or bool array.
    """
    x = np.asarray(x, dtype=np.float64)
    bits = np.ascontiguousarray(x).view(np.uint64)
    na = np.isnan(x) & ((bits & np.uint64(0xFFFFFFFF)) == np.uint64(1954))
    return bool(na) if na.ndim == 0 else na


# ----- Sums, means, variances -----

@njit(cache=True)
def _sum(x):
    total = 0.0
    for v in x:
        total += v
    return total


def r_sum(x):
    """`sum(x)`: the elements added in order. Returns: float."""
    return float(_sum(np.asarray(x, dtype=np.float64).ravel()))


@njit(cache=True)
def _mean(x):
    n = len(x)
    s = 0.0
    for v in x:
        s += v
    s /= n
    if np.isfinite(s):
        t = 0.0
        for v in x:
            t += v - s
        s += t / n
    return s


@njit(cache=True)
def _cumsum(x):
    out = np.empty(len(x))
    total = 0.0
    for i in range(len(x)):
        total += x[i]
        out[i] = total
    return out


def r_cumsum(x):
    """`cumsum(x)`: the running sums, added in order. Returns: float array."""
    return _cumsum(np.asarray(x, dtype=np.float64))


def r_mean(x, na_rm=False):
    """`mean(x, na.rm = na_rm)` of doubles: the mean with R's second pass.

    Parameters: x (vector), na_rm (drop missing values). Returns: float.
    """
    x = np.asarray(x, dtype=np.float64)
    if na_rm:
        x = x[~np.isnan(x)]
    if len(x) == 0:
        return float("nan")
    return float(_mean(x))


def r_mean_integer(x, na_rm=False):
    """`mean(x, na.rm = na_rm)` of integers or logicals: the sum divided by
    the length, without the second pass R makes for doubles.

    Parameters: x (vector of whole numbers; NaN is a missing value), na_rm.
    Returns: float.
    """
    x = np.asarray(x, dtype=np.float64)
    if na_rm:
        x = x[~np.isnan(x)]
    if len(x) == 0:
        return float("nan")
    return float(_sum(x) / len(x))


@njit(cache=True)
def _var(x, complete):
    """Variance of a vector without missing values, as R's cov.c adds it.

    With `complete` False (cov_na_1, `var(x)`): R's compiled loop adds the
    squared deviations eight at a time, each square rounded before it is
    added; the ones left over at the end (fewer than eight) are added with
    one rounding each (muladd) on arm64. With `complete` True (cov_complete1,
    `var(x, na.rm = TRUE)`): every square is added with one rounding."""
    n = len(x)
    mean = _mean(x)
    blocks = 0 if complete else n - n % 8
    total = 0.0
    for k in range(n):
        d = x[k] - mean
        if k < blocks:
            total += d * d
        else:
            total = muladd(d, d, total)
    return total / (n - 1)


def r_var(x, na_rm=False):
    """`var(x, na.rm = na_rm)`: the variance of a vector.

    Parameters: x (vector), na_rm (drop missing values).
    Returns: float; NA for fewer than two values, NaN when a value is
    missing and `na_rm` is False.
    """
    x = np.asarray(x, dtype=np.float64)
    if na_rm:
        x = x[~np.isnan(x)]
    elif np.isnan(x).any():
        return NA
    if len(x) <= 1:
        return NA
    return float(_var(x, na_rm))


def r_sd(x, na_rm=False):
    """`sd(x, na.rm = na_rm)`. Returns: float."""
    variance = r_var(x, na_rm)
    return variance if math.isnan(variance) else float(np.sqrt(variance))


def quantile(x, prob):
    """`quantile(x, probs = prob)[[1]]`, type 7. Returns: float."""
    x = np.sort(np.asarray(x, dtype=np.float64))
    index = 1 + max(len(x) - 1, 0) * prob
    lo, hi = int(np.floor(index)), int(np.ceil(index))
    q = x[lo - 1]
    if index > lo and x[hi - 1] != q:
        h = index - lo
        q = (1 - h) * q + h * x[hi - 1]
    return float(q)


# ----- Smallest and largest, missing values as in R -----

def r_min(*values):
    """`min(...)` of numbers: NaN when one is missing. Returns: float."""
    values = [float(v) for v in values]
    return float("nan") if any(math.isnan(v) for v in values) else min(values)


def r_max(*values):
    """`max(...)` of numbers: NaN when one is missing. Returns: float."""
    values = [float(v) for v in values]
    return float("nan") if any(math.isnan(v) for v in values) else max(values)


def pmin(bound, x):
    """`pmin(bound, x)` for a number and a vector: the smaller of the two
    at every element, NaN where x is missing. Returns: float array."""
    x = np.asarray(x, dtype=np.float64)
    return np.where(np.isnan(x), np.nan, np.where(x > bound, bound, x))


def pmax(bound, x):
    """`pmax(bound, x)`: the larger of the two at every element, NaN where
    one is missing. Parameters: bound (number or array), x (array).
    Returns: float array."""
    x = np.asarray(x, dtype=np.float64)
    bound = np.broadcast_to(np.asarray(bound, dtype=np.float64), x.shape)
    return np.where(np.isnan(x) | np.isnan(bound), np.nan, np.where(x < bound, bound, x))


def tabulate(x, nbins):
    """`tabulate(x, nbins)`: how often each of 1, ..., nbins occurs in x;
    other values are ignored.

    Parameters: x (whole numbers), nbins. Returns: float array.
    """
    x = np.asarray(x, dtype=np.float64)
    x = x[~np.isnan(x)].astype(np.int64)
    x = x[(x >= 1) & (x <= nbins)]
    return np.bincount(x, minlength=int(nbins) + 1)[1:].astype(np.float64)


def count_frequencies_from_zero(x, largest):
    """`c(sum(x == 0), tabulate(x, largest))`: how often each of 0, ...,
    largest occurs in x.

    Parameters: x (whole numbers), largest. Returns: float array.
    """
    x = np.asarray(x, dtype=np.float64)
    return np.concatenate(([float(np.sum(x == 0))], tabulate(x, largest)))


def r_order(x):
    """`order(x)`, zero-based: the stable order of a vector, missing values
    last. Returns: int array."""
    return np.argsort(np.asarray(x, dtype=np.float64), kind="stable")


def which_min(x):
    """`which.min(x)`, zero-based: the position of the first smallest value,
    missing values left out; None when every value is missing or there is
    none.

    Parameters: x (vector). Returns: int or None.
    """
    x = np.asarray(x, dtype=np.float64)
    if len(x) == 0 or np.isnan(x).all():
        return None
    return int(np.nanargmin(x))


def colon(first, last):
    """`first:last`: the whole numbers from first to last, counting down
    when last is smaller (so 3:2 is 3, 2).

    Parameters: first, last (integers). Returns: float array.
    """
    step = 1 if last >= first else -1
    return np.arange(first, last + step, step, dtype=np.float64)


# ----- Tables in files -----

_EXACT_POWER = 22  # largest power of ten that is an exact double


def _scaled_integer(r, kp, digits):
    """r / 10^kp rounded to the nearest integer (ties to even), and kp, one
    lower when the quotient has fewer than `digits` digits: a step of R's
    scientific()."""
    if abs(kp) <= _EXACT_POWER:
        scaled = r / 10.0 ** kp if kp >= 0 else r * 10.0 ** -kp
    elif kp <= -308:
        scaled = (r * 1e+303) / math.pow(10.0, kp + 303)
    else:
        scaled = r / math.pow(10.0, kp)
    if scaled < 10.0 ** (digits - 1):
        scaled *= 10.0
        kp -= 1
    return float(round(scaled)), kp


def format_number(x, digits=15):
    """A number as R writes it with `digits` significant digits
    (`as.character()`, `write.csv()`).

    R finds the significant digits a number needs by scaling it to an
    integer of `digits` digits and dropping the zeros at its end, then
    writes it in fixed or scientific notation, whichever is shorter.

    Parameters: x (float), digits (significant digits, at most 15).
    Returns: the text; NA for a missing value or a NaN, as write.csv().
    """
    x = float(x)
    if math.isnan(x):
        return "NA"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    if x == 0:
        return "0"
    negative, r = int(x < 0), abs(x)
    kp = math.floor(math.log10(r)) - digits + 1
    alpha, kp = _scaled_integer(r, kp, digits)
    significant = digits
    for _ in range(digits):
        alpha /= 10.0
        if alpha != math.floor(alpha):
            break
        significant -= 1
    if significant == 0:
        significant = 1
        kp += 1
    power = kp + digits - 1
    left = power + 1
    cut = min(max(digits - power, 0), _EXACT_POWER)  # the decimal place rounding cuts off
    if 0 < power <= _EXACT_POWER and r < 10.0 ** power - 0.5 / 10.0 ** cut:
        left -= 1
    decimals = max(0, significant - power - 1)
    fixed_width = negative + max(left, 1) + decimals + (decimals > 0)
    mantissa = significant - 1
    scientific_width = negative + (mantissa > 0) + mantissa + 4 + (2 if abs(power) >= 100 else 1)
    return f"{x:.{decimals}f}" if fixed_width <= scientific_width else f"{x:.{mantissa}e}"


def write_csv(table, file):
    """`write.csv(table, file, row.names = FALSE)`: numbers with 15
    significant digits as R writes them, text in quotes, NA for a missing
    value.

    Parameters: table (DataFrame), file (path). Returns: None.
    """
    def quoted(text):
        return '"' + str(text).replace('"', '""') + '"'

    columns = []
    for name in table.columns:
        values = table[name]
        if pd.api.types.is_bool_dtype(values):
            columns.append(["TRUE" if value else "FALSE" for value in values])
        elif pd.api.types.is_integer_dtype(values):
            columns.append([str(value) for value in values])
        elif pd.api.types.is_float_dtype(values):
            columns.append([format_number(float(value)) for value in values])
        else:
            columns.append(["NA" if value is None or (isinstance(value, float) and math.isnan(value))
                            else quoted(value) for value in values])
    with open(file, "w") as f:
        f.write(",".join(quoted(name) for name in table.columns) + "\n")
        for row in zip(*columns):
            f.write(",".join(row) + "\n")


# ----- Parallel processes -----

_job = None


def _run_job(item):
    try:
        return _job(item)
    except Exception as error:  # as mclapply: the error is the result
        return error


def parallel_map(function, items, cores=1):
    """`parallel::mclapply(items, function, mc.cores = cores,
    mc.preschedule = FALSE)`.

    The processes are forked, so each starts with the memory and the random
    number stream of this one. A job that fails gives its exception as its
    result, as mclapply gives the error.

    Parameters: function (of one item), items (list), cores (number of
    parallel processes).
    Returns: list of the results, in the order of `items`.
    """
    global _job
    items = list(items)
    cores = int(cores)
    _job = function
    if cores <= 1 or len(items) <= 1:
        return [_run_job(item) for item in items]
    with multiprocessing.get_context("fork").Pool(min(cores, len(items))) as pool:
        return pool.map(_run_job, items, chunksize=1)
