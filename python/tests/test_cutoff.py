# -----------------------------------------------------------------------------
# File:    python/tests/test_cutoff.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Tests of cutoff.py: the truncated likelihood, the penalties, the
#          chosen cutoff, the bound D and the selected null fits, against
#          known values and against R, double for double.
# Usage:   python -m pytest python/tests/test_cutoff.py
# -----------------------------------------------------------------------------
import math

import numpy as np
import pandas as pd

from conftest import numbers, same
from inputs import COUNTS, WHOLE
from zinull import cutoff, rbase
from zinull.em import count_up_to, em_estimates

COUNT = count_up_to(COUNTS, 20)


def test_truncated_neg_loglik(r):
    values = [cutoff.truncated_neg_loglik(COUNT, 0.4, 1.5, 0.3, 4),
              cutoff.truncated_neg_loglik(COUNT, 0, 1.5, 0, 20)]
    same(values, r("cutoff_neg_loglik"))
    # one count of zero under a Poisson null: -log(exp(-lambda)) = lambda
    assert cutoff.truncated_neg_loglik([1.0, 0.0], 0.0, 2.0, 0.0, 1) == 2.0


def test_cutoff_penalties(r):
    penalties = cutoff.cutoff_penalties(COUNT)
    assert len(penalties["C1"]) == len(penalties["C2"]) == 20
    same(penalties["C1"], r("cutoff_penalties_C1"))
    same(penalties["C2"], r("cutoff_penalties_C2"))


def test_choose_cutoff(r):
    """The cutoff after the first turn upwards (R's which(diff(sign(diff)))
    > 0) + 2), or the overall minimum when that comes first."""
    chosen = [cutoff.choose_cutoff(criterion) for criterion in
              ([5, 4, 3, 4, 1], [1, 2, 3], [3, 2, 1], [7], [2, 2, 1, 5])]
    assert chosen == [4, 1, 3, 1, 3]
    same(chosen, r("cutoff_choose"))


def test_rejection_bound(r):
    bounds = [cutoff.rejection_bound("ZIGP", 5, 1.5, 0.3, 1000),
              cutoff.rejection_bound("GP", 5, 1.5, 0.3, 1000, cap=9),
              cutoff.rejection_bound("ZIP", 5, 1.5, 0, 1000),
              cutoff.rejection_bound("P", 9, 1.5, 0, 1000)]
    same(bounds, r("cutoff_rejection_bound"))
    assert bounds[1] == 9  # the cap
    assert bounds[2] == math.ceil(math.log(1000)) and bounds[3] == 10  # max(lambda, log n), C + 1


def test_select_model_fits(r):
    fits = cutoff.select_model_fits("ZIGP", em_estimates(COUNTS)["ZIGP"], COUNT,
                                    cutoff.cutoff_penalties(COUNT))
    assert list(fits.columns) == ["cutoff", "model", "eta", "lambda", "theta", "pi0", "C", "D"]
    assert list(fits["cutoff"]) == ["C1", "C2"]
    same(numbers(fits), r("cutoff_model_fits"))


def test_select_null_fits(r):
    fits = cutoff.select_null_fits(COUNTS)
    assert [f"{c} {m}" for c, m in zip(fits["cutoff"], fits["model"])] \
        == r.lines("cutoff_null_fits_names")
    same(numbers(fits), r("cutoff_null_fits"))
    capped = cutoff.select_null_fits(WHOLE, cap=12)
    assert (capped["D"] <= 12).all()
    same(numbers(capped), r("cutoff_null_fits_capped"))


def test_mask_absent_parameters(r):
    fits = cutoff.select_null_fits(COUNTS)
    masked = cutoff.mask_absent_parameters(fits)
    same(numbers(masked[["eta", "theta"]]), r("cutoff_masked"))
    assert rbase.is_na(masked["eta"].to_numpy())[masked["model"].isin(["GP", "P"])].all()
    assert rbase.is_na(masked["theta"].to_numpy())[masked["model"].isin(["ZIP", "P"])].all()
    assert not np.isnan(fits["eta"]).any()  # the argument is not changed
    small = pd.DataFrame({"model": ["ZIGP", "P"], "eta": [0.1, 0.0], "theta": [0.2, 0.0]})
    assert np.isnan(cutoff.mask_absent_parameters(small)[["eta", "theta"]].to_numpy()).tolist() \
        == [[False, False], [True, True]]
