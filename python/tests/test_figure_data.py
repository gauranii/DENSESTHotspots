# -----------------------------------------------------------------------------
# File:    python/tests/test_figure_data.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Tests of figure_data.py: the simulated data sets of the
#          histograms, the bar heights and the series of the null-fit plots,
#          against known values and against R, double for double.
# Usage:   python -m pytest python/tests/test_figure_data.py
# -----------------------------------------------------------------------------
import numpy as np

from conftest import numbers, same
from inputs import COUNTS
from zinull import figure_data, rng
from zinull.real_data import fit_domain, read_domains

MIXTURE = {"n": 40, "offset": 6, "min_kept": 0, "pi0": 0.8, "eta": 0.8, "lam": 1.5, "theta": 0.3}


def test_settings():
    panels = figure_data.HISTOGRAM_SETTINGS
    assert len(panels) == 24 and set(panels["order"]) == {"non_null", "null"}
    assert list(panels["figure"].unique()) == ["Fig1"] + ["FigS%d" % i for i in range(1, 7)]
    assert figure_data.CUTOFF == 5 and figure_data.BINOMIAL_SIZE == 250
    assert len(figure_data.PROTEIN_DOMAINS) == 6
    assert list(figure_data.NULL_FIT_DOMAINS) == ["FigS%d" % i for i in range(7, 13)]
    assert set(figure_data.NULL_FIT_DOMAINS.values()) == set(figure_data.PROTEIN_DOMAINS)


def test_draw_null_count(r):
    rng.set_seed(31)
    draws = [figure_data.draw_null_count(0.4, 1.5, 0.3), figure_data.draw_null_count(0.4, 1.5, 0.3),
             figure_data.draw_null_count(0, 3, 0)]
    same(draws, r("figure_draw_null_count"))


def test_simulate_mixture(r):
    mixture = figure_data.simulate_mixture(40, lambda: 2, offset=6, min_kept=9, pi0=0.5, eta=0.4,
                                           lam=1.5, theta=0.3, seed=5)
    same(numbers(mixture), r("figure_mixture"))
    assert (mixture["count"][mixture["non_null"] == 1] >= 9).all()
    # without a seed the draws continue the stream: set.seed() gives the same data
    rng.set_seed(5)
    again = figure_data.simulate_mixture(40, lambda: 2, offset=6, min_kept=9, pi0=0.5, eta=0.4,
                                         lam=1.5, theta=0.3)
    same(numbers(again), numbers(mixture))


def test_simulate_geometric_mixture_and_simulate_binomial_mixture(r):
    geometric = figure_data.simulate_geometric_mixture(0.08, seed=6, **MIXTURE)
    same(numbers(geometric), r("figure_geometric"))
    binomial = figure_data.simulate_binomial_mixture(0.2, n=40, offset=-34, min_kept=6, pi0=0.35,
                                                     eta=0.4, lam=1, theta=0.3, seed=6)
    same(numbers(binomial), r("figure_binomial"))
    rng.set_seed(6)  # the seed argument is set.seed() right before the draws
    same(numbers(figure_data.simulate_geometric_mixture(0.08, **MIXTURE)), numbers(geometric))


def test_simulate_histogram_data(r):
    for row in (1, 4, 13):  # a geometric, a binomial and a mixed panel
        data = figure_data.simulate_histogram_data(figure_data.HISTOGRAM_SETTINGS.iloc[row - 1],
                                                   50, 889302)
        assert len(data) == 50
        same(numbers(data), r("figure_histogram_data_%d" % row))


def test_bar_heights(r):
    heights = figure_data.bar_heights(COUNTS)
    assert list(heights.index) == list(range(21)) and heights.sum() == len(COUNTS)
    same(heights, r("figure_bar_heights"))
    assert list(figure_data.bar_heights([0, 0, 3])) == [2, 0, 0, 1]


def test_protein_counts(r, data_file):
    counts = figure_data.protein_counts("domB", data_file)
    assert len(counts) == 120 and counts.dtype.kind == "i"
    same(counts, r("real_protein_counts"))


def test_null_fit_series(r, data_file):
    x = read_domains(data_file)["domA"]
    series = figure_data.null_fit_series(x, fit_domain(x), "C1")
    assert list(series.columns) == ["fEstimation", "C", "values"]
    assert list(series["fEstimation"]) == r.lines("figure_null_fit_series_names")
    same(numbers(series), r("figure_null_fit_series"))
    observed = series[series["fEstimation"] == "fHat"]
    assert observed["values"].iloc[0] == np.mean(x == 0)
