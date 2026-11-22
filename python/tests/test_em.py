# -----------------------------------------------------------------------------
# File:    python/tests/test_em.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Tests of em.py: starting values, the sums of the E-step, one EM
#          step, the fit at one cutoff and at every cutoff, against R, double
#          for double, for the four null models.
# Usage:   python -m pytest python/tests/test_em.py
# -----------------------------------------------------------------------------
import numpy as np
import pytest

from conftest import same
from inputs import COUNTS, HOLD, SETTING, WHOLE
from zinull import em
from zinull.simulation import simulate_labelled_counts


def test_constants(r):
    assert em.MODELS == ("ZIGP", "ZIP", "GP", "P")
    same([em.EM_TOLERANCE], r("em_tolerance"))


def test_count_up_to(r):
    same(em.count_up_to(COUNTS, 7), r("em_count_up_to"))
    np.testing.assert_array_equal(em.count_up_to([0, 0, 1, 3, 9], 3), [2, 1, 0, 1])


@pytest.mark.parametrize("model", em.MODELS)
def test_em_start(model, r):
    """Doubles and integers are averaged differently in R; both are kept."""
    same(np.concatenate((em.em_start(COUNTS, model), em.em_start(WHOLE, model))),
         r("em_start_" + model))


def test_em_start_of_a_model_without_a_parameter_is_zero():
    assert em.em_start(COUNTS, "P")[0] == 0 and em.em_start(COUNTS, "P")[2] == 0
    assert em.em_start(COUNTS, "ZIP")[2] == 0 and em.em_start(COUNTS, "GP")[0] == 0
    with pytest.raises(ValueError):
        em.em_start(COUNTS, "NB")


@pytest.mark.parametrize("model", em.MODELS)
def test_em_clamp(model, r):
    same(np.concatenate((em.em_clamp([-0.2, -1, -0.5], model), em.em_clamp(HOLD, model))),
         r("em_clamp_" + model))


def test_em_sums(r):
    sums = em.em_sums(HOLD, em.count_up_to(COUNTS, 3), 20, dispersion=True)
    assert list(sums) == ["zero", "gp", "above", "total_above", "total", "lambda", "theta"]
    same(list(sums.values()), r("em_sums"))
    # at the last cutoff the range of the counts above C is K + 1, K, as in R
    last = em.em_sums(HOLD, em.count_up_to(COUNTS, 20), 20, dispersion=True)
    same(list(last.values()), r("em_sums_last"))
    plain = em.em_sums([0.4, 1.5, 0], em.count_up_to(COUNTS, 3), 20, dispersion=False)
    assert list(plain) == ["zero", "gp", "above", "total_above", "total"]
    same(list(plain.values()), r("em_sums_plain"))


@pytest.mark.parametrize("model", em.MODELS)
def test_em_step(model, r):
    hold = em.em_clamp(HOLD, model)
    steps = [em.em_step(hold, em.count_up_to(COUNTS, C), 20, model) for C in (3, 20)]
    same(np.concatenate(steps), r("em_step_" + model))


@pytest.mark.parametrize("model", em.MODELS)
def test_em_fit_and_em_estimates(model, r):
    same(em.em_fit(COUNTS, 4, model, em.em_start(COUNTS, model)), r("em_fit_" + model))
    estimates = em.em_estimates(COUNTS)[model]
    assert estimates.shape == (3, 20)
    same(estimates.ravel(order="F"), r("em_estimates_" + model))
    same(estimates[:, 3], r("em_fit_" + model))


def test_em_stops_where_r_stops(r):
    """On this data set a parameter becomes NaN and R stops with an error."""
    assert r("em_stops")[0] == 1
    failing = simulate_labelled_counts(150, 5, SETTING, 11)["value"]
    with pytest.raises(ValueError):
        em.em_estimates(failing)
