# -----------------------------------------------------------------------------
# File:    python/zinull/figure_data.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. Figures 1 and 2, Web Figures S1 to S12.
# Purpose: Data of the figures: the settings of the 24 simulated data sets,
#          the generators of a null / non-null mixture of counts, the bar
#          heights of a histogram of counts, and the observed and fitted
#          null probabilities of a protein domain (Web Figures S7 to S12).
#          The Python counterpart of R/figure_data.R, with the same results.
# Inputs:  data/realdata.txt (protein domain counts, read by protein_counts()).
# Outputs: None.
# Usage:   from zinull.figure_data import HISTOGRAM_SETTINGS, simulate_histogram_data
#          data = simulate_histogram_data(HISTOGRAM_SETTINGS.iloc[0], 1000, 889302)
# Source:  R/figure_data.R.
# -----------------------------------------------------------------------------
import io

import numpy as np
import pandas as pd

from . import rng
from .em import count_up_to
from .real_data import read_domains
from .zigp import generate_zigp, zigp_density

CUTOFF = 5           # C: non-null counts are built as C + 1 + D + signal + null count
BINOMIAL_SIZE = 250  # number of trials of the binomial signal

# Settings of the simulated data sets, one row per panel.
#   original  file name of the panel in the original code
#   signal    distribution of the signal of a non-null count, with parameter p
#   D         shift of the non-null counts
#   min_kept  a non-null count below this value is drawn again
#   order     the counts drawn first; the other ones are drawn on top
#   pi0       proportion of null counts
#   eta, lambda, theta  zero inflation, rate and dispersion of the null counts
HISTOGRAM_SETTINGS = pd.read_csv(io.StringIO("""
  figure    panel model  original signal    D   min_kept order    pi0  p    eta lambda theta
  Fig1      A     ZIP1   data3    geometric 0   0        non_null 0.80 0.08 0.8 1.5    0
  Fig1      B     ZIGP1  data1    geometric 0   0        non_null 0.80 0.08 0.8 1.5    0.30
  Fig1      C     ZIGP2  data2    geometric 0   0        non_null 0.80 0.08 0.8 3      0.30
  Fig1      D     ZIP2   bar6     binomial  -40 6        non_null 0.35 0.20 0.4 1.5    0
  Fig1      E     ZIGP3  bar4     binomial  -40 6        non_null 0.35 0.20 0.4 1      0.30
  Fig1      F     ZIGP4  bar5     binomial  -40 6        non_null 0.35 0.20 0.4 4      0.30
  FigS1     A     ZIGP5  data5    geometric 0   0        non_null 0.85 0.08 0.4 1      0.15
  FigS1     B     ZIGP6  data6    geometric 0   0        non_null 0.85 0.08 0.4 3      0.15
  FigS1     C     ZIP3   data7    geometric 0   0        non_null 0.85 0.08 0.4 1.5    0
  FigS2     A     ZIGP7  bar1     binomial  -40 6        non_null 0.80 0.20 0.4 1      0.20
  FigS2     B     ZIGP8  bar2     binomial  -40 6        non_null 0.80 0.20 0.4 3      0.20
  FigS2     C     ZIP4   bar3     binomial  -40 6        non_null 0.80 0.20 0.4 1.5    0
  FigS3     A     MZIGP1 mixed1   geometric -5  0        null     0.80 0.08 0.8 1.5    0.30
  FigS3     B     MZIGP2 mixed2   geometric -5  0        null     0.80 0.08 0.8 3      0.30
  FigS3     C     MZIP1  mixed3   geometric -5  0        null     0.80 0.08 0.8 1.5    0
  FigS4     A     MZIGP3 mixed4   binomial  -40 0        null     0.35 0.20 0.4 1      0.30
  FigS4     B     MZIGP4 mixed5   binomial  -40 0        null     0.35 0.20 0.4 4      0.30
  FigS4     C     MZIP2  mixed6   binomial  -40 0        null     0.35 0.20 0.4 1.5    0
  FigS5     A     MZIGP5 mixed7   geometric -5  0        null     0.85 0.08 0.4 1      0.15
  FigS5     B     MZIGP6 mixed8   geometric -5  0        null     0.85 0.08 0.4 3      0.15
  FigS5     C     MZIP3  mixed9   geometric -5  0        null     0.85 0.08 0.4 1.5    0
  FigS6     A     MZIGP7 mixed10  binomial  -40 0        null     0.80 0.20 0.4 1      0.20
  FigS6     B     MZIGP8 mixed11  binomial  -40 0        null     0.80 0.20 0.4 3      0.20
  FigS6     C     MZIP4  mixed12  binomial  -40 0        null     0.80 0.20 0.4 1.5    0
"""), sep=r"\s+", keep_default_na=False)

# Protein domains of Figure 2, in the order of the panels.
PROTEIN_DOMAINS = ("cd00031", "cd00180", "cd00204", "cd00882", "cd00054", "pfam00001")


def draw_null_count(eta, lam, theta):
    """One null count: a ZIGP draw in the parametrization of the paper.

    Parameters: eta (zero-inflation probability), lam (rate), theta
    (dispersion; 0 gives the zero-inflated Poisson).
    Returns: one count (float).
    """
    return float(generate_zigp(1, mu=lam / (1 - theta), phi=1 / (1 - theta), omega=eta)[0])


def simulate_mixture(n, draw_signal, offset, min_kept, pi0, eta, lam, theta, seed=None):
    """Simulate a mixture of null and non-null counts.

    Each position is non-null when a uniform draw is at least `pi0`. A null
    count is a ZIGP draw. A non-null count is `offset` plus a signal plus a
    ZIGP draw (drawn in this order); it is dropped, and the position drawn
    again, when it is below `min_kept`.

    Parameters: n (number of counts), draw_signal (function without
    arguments that draws one signal), offset (constant added to every
    non-null count), min_kept (smallest non-null count that is kept), pi0
    (proportion of null counts), eta, lam, theta (parameters of the null
    counts), seed (seed of the random numbers, set right before the draws;
    None continues the stream).
    Returns: DataFrame with `n` rows: `count` and `non_null` (0 or 1).
    """
    if seed is not None:
        rng.set_seed(seed)
    count, non_null = [], []
    while len(count) < n:
        is_non_null = bool(rng.runif(1)[0] >= pi0)
        signal = offset + draw_signal() if is_non_null else 0
        value = signal + draw_null_count(eta, lam, theta)
        if not is_non_null or value >= min_kept:
            count.append(float(value))
            non_null.append(float(is_non_null))
    return pd.DataFrame({"count": count, "non_null": non_null})


def simulate_geometric_mixture(p, seed=None, **mixture):
    """Mixture whose non-null counts carry a geometric signal.

    Parameters: p (success probability of the geometric signal), seed (seed
    of the random numbers; None continues the stream), mixture (the other
    arguments of `simulate_mixture()`).
    Returns: DataFrame with columns `count` and `non_null`.
    """
    return simulate_mixture(draw_signal=lambda: float(rng.rgeom(1, p)[0]), seed=seed, **mixture)


def simulate_binomial_mixture(p, size=BINOMIAL_SIZE, seed=None, **mixture):
    """Mixture whose non-null counts carry a binomial signal.

    Parameters: p (success probability of the binomial signal), size (number
    of trials of the binomial signal), seed (seed of the random numbers;
    None continues the stream), mixture (the other arguments of
    `simulate_mixture()`).
    Returns: DataFrame with columns `count` and `non_null`.
    """
    return simulate_mixture(draw_signal=lambda: float(rng.rbinom(1, size, p)[0]), seed=seed,
                            **mixture)


def simulate_histogram_data(setting, n, seed, cutoff=CUTOFF):
    """Simulate the data set of one panel.

    Parameters: setting (one row of `HISTOGRAM_SETTINGS`), n (number of
    counts), seed (seed of the random numbers), cutoff (the constant C of
    the non-null counts).
    Returns: DataFrame with columns `count` and `non_null`.
    """
    simulate = {"geometric": simulate_geometric_mixture,
                "binomial": simulate_binomial_mixture}[setting["signal"]]
    return simulate(p=float(setting["p"]), n=n, offset=cutoff + 1 + int(setting["D"]),
                    min_kept=int(setting["min_kept"]), pi0=float(setting["pi0"]),
                    eta=float(setting["eta"]), lam=float(setting["lambda"]),
                    theta=float(setting["theta"]), seed=seed)


def bar_heights(x):
    """Bar heights of a histogram of counts: the frequency of 0, 1, ..., max(x).

    Parameters: x (non-negative counts).
    Returns: pandas Series of the frequencies, indexed by the counts.
    """
    largest = int(np.max(x))
    return pd.Series(count_up_to(x, largest), index=range(largest + 1))


def protein_counts(domain, data_file):
    """Counts of one protein domain.

    Parameters: domain (name of the domain, a column of the data file),
    data_file (path of data/realdata.txt).
    Returns: integer array of the counts at the positions of the domain.
    """
    return read_domains(data_file)[domain]


# Protein domains of Web Figures S7 to S12, one figure per domain; panel A
# uses the cutoff criterion C1, panel B the criterion C2.
NULL_FIT_DOMAINS = {"FigS7": "cd00031", "FigS8": "cd00054", "FigS9": "cd00180",
                    "FigS10": "cd00204", "FigS11": "cd00882", "FigS12": "pfam00001"}


def null_fit_series(x, fits, criterion):
    """Observed and fitted null probabilities of a protein domain.

    For one cutoff criterion: the observed share of positions with count
    0, 1, ..., K ("fHat"), and pi0 times the fitted ZIGP and ZIP probability
    of each count up to the cutoff chosen for that model ("zigp", "zip"),
    where K is the larger of the two cutoffs.

    Parameters: x (counts of the domain), fits (null fits of the domain,
    from `fit_domain()`), criterion (one of `CUTOFFS`).
    Returns: DataFrame with the columns `fEstimation` (series), `C` (count)
    and `values` (probability).
    """
    def fit(model):
        return fits[(fits["cutoff"] == criterion) & (fits["model"] == model)].iloc[0]

    def fitted(f):
        counts = np.arange(int(f["C"]) + 1, dtype=np.float64)
        return f["pi0"] * zigp_density(counts, f["eta"], f["lambda"], f["theta"])

    zigp, zip_ = fit("ZIGP"), fit("ZIP")
    K = int(max(zigp["C"], zip_["C"]))
    series = {"fHat": count_up_to(x, K) / len(x), "zigp": fitted(zigp), "zip": fitted(zip_)}
    return pd.DataFrame({
        "fEstimation": [name for name, values in series.items() for _ in values],
        "C": [count for values in series.values() for count in range(len(values))],
        "values": np.concatenate(list(series.values())),
    })
