# -----------------------------------------------------------------------------
# File:    python/tests/test_rng.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Tests of rng.py: the same stream as R after set.seed(), draw for
#          draw, for every sampler.
# Usage:   python -m pytest python/tests/test_rng.py
# -----------------------------------------------------------------------------
import numpy as np

from conftest import same
from zinull import rng


def test_set_seed_restarts_the_stream():
    rng.set_seed(1)
    first = rng.runif(3)
    rng.set_seed(1)
    same(rng.runif(3), first)
    # the first uniform numbers of R after set.seed(1)
    np.testing.assert_allclose(first, [0.2655087, 0.3721239, 0.5728534], atol=5e-8)


def test_one_stream_of_every_sampler_as_r(r):
    """The draws of make_r_values.R, in its order, from the same seed."""
    rng.set_seed(889192)
    same(rng.runif(5), r("rng_runif"))
    same(rng.rnorm(4, 1, 2), r("rng_rnorm"))
    same(rng.rexp(4), r("rng_rexp"))
    binomial = np.concatenate((rng.rbinom(6, 250, 0.2), rng.rbinom(6, 10, 0.3),
                               rng.rbinom(4, 40, 0.9)))
    same(binomial, r("rng_rbinom"))
    same(np.concatenate((rng.rpois(6, 3.5), rng.rpois(6, 40))), r("rng_rpois"))
    same(rng.rgeom(8, 0.08), r("rng_rgeom"))


def test_draws_are_in_range():
    rng.set_seed(5)
    u = rng.runif(1000)
    assert ((u > 0) & (u < 1)).all()
    assert (rng.rexp(100) > 0).all()
    b = rng.rbinom(200, 250, 0.2)
    assert ((b >= 0) & (b <= 250) & (b == np.floor(b))).all()
    assert (rng.rpois(100, 2.0) >= 0).all() and (rng.rgeom(100, 0.3) >= 0).all()
    assert abs(rng.rnorm(4000).mean()) < 0.1


def test_parameters_out_of_range_give_nan():
    assert np.isnan(rng.rbinom(1, 10, 1.5)[0])
    assert np.isnan(rng.rgeom(1, 0.0)[0])
    assert np.isnan(rng.rpois(1, -1.0)[0])
