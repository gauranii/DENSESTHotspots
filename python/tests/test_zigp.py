# -----------------------------------------------------------------------------
# File:    python/tests/test_zigp.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Tests of zigp.py: the ZIGP draws and probability functions against
#          known values and against R, double for double.
# Usage:   python -m pytest python/tests/test_zigp.py
# -----------------------------------------------------------------------------
import math

import numpy as np
import pytest

from conftest import same
from zinull import rng, zigp

COUNTS = np.arange(9)


def test_small_is_the_constant_of_r(r):
    same([zigp.SMALL], r("zigp_small"))


def test_generate_zigp_draws_as_r(r):
    first = zigp.generate_zigp(12, mu=1.5 / 0.7, phi=1 / 0.7, omega=0.4, seed=21)
    same(first, r("zigp_generate"))
    # without a seed the draws continue the stream
    same(zigp.generate_zigp(3, mu=3, phi=1, omega=0), r("zigp_generate_next"))
    # the seed argument is set.seed() right before the draws
    rng.set_seed(21)
    same(zigp.generate_zigp(12, mu=1.5 / 0.7, phi=1 / 0.7, omega=0.4), first)


def test_generate_zigp_known_cases():
    assert (zigp.generate_zigp(50, mu=2.0, phi=1.0, omega=1.0, seed=3) == 0).all()
    draws = zigp.generate_zigp(4000, mu=2.0, phi=1.0, omega=0.0, seed=3)  # Poisson(2)
    assert abs(draws.mean() - 2.0) < 0.1 and (draws == np.floor(draws)).all()
    with pytest.raises(ValueError):
        zigp.generate_zigp(1, mu=1.0, phi=1.0, omega=1.5)


def test_gp_pmf(r):
    same(np.concatenate((zigp.gp_pmf(COUNTS, 1.5, 0.3), zigp.gp_pmf(COUNTS, 2, 0))),
         r("zigp_gp_pmf"))
    # theta = 0 is the Poisson distribution
    assert zigp.gp_pmf(3, 2.0, 0.0) == pytest.approx(math.exp(-2) * 8 / 6, rel=1e-14)
    assert zigp.gp_pmf(np.arange(60), 1.5, 0.3).sum() == pytest.approx(1.0, abs=1e-9)


def test_gp_pmf_times_x(r):
    same(zigp.gp_pmf_times_x(COUNTS[1:], 1.5, 0.3), r("zigp_gp_pmf_times_x"))
    assert zigp.gp_pmf_times_x(4, 1.5, 0.3) == pytest.approx(4 * zigp.gp_pmf(4, 1.5, 0.3),
                                                             rel=1e-14)


def test_densities_are_kept_positive(r):
    same(zigp.gp_density(COUNTS, 1.5, 0.3), r("zigp_gp_density"))
    same(zigp.zigp_density(COUNTS, 0.4, 1.5, 0.3), r("zigp_zigp_density"))
    assert zigp.gp_density(2, 1.5, 0.3) == zigp.gp_pmf(2, 1.5, 0.3) + zigp.SMALL
    assert zigp.zigp_density(0, 1.0, 1.5, 0.3) == 1.0 + zigp.SMALL


def test_zigp_mass_below(r):
    same([zigp.zigp_mass_below(0.4, 1.5, 0.3, 5), zigp.zigp_mass_below(0, 2, 0, 1)],
         r("zigp_mass_below"))
    assert zigp.zigp_mass_below(0, 2, 0, 1) == pytest.approx(3 * math.exp(-2), rel=1e-14)


def test_gp_bound(r):
    same([zigp.gp_bound(1.5, 0.3, 1000), zigp.gp_bound(0.2, 0.6, 150)], r("zigp_gp_bound"))
