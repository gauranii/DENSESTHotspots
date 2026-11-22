# -----------------------------------------------------------------------------
# File:    python/tests/test_rmath.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Tests of rmath.py: known values of every function, and the same
#          doubles as R.
# Usage:   python -m pytest python/tests/test_rmath.py
# -----------------------------------------------------------------------------
import math

import numpy as np
from numba import njit

from conftest import same
from inputs import DECIMALS
from zinull import rmath


def test_muladd_adds_a_product_with_one_rounding():
    @njit
    def fused(a, b, c):
        return rmath.muladd(a, b, c)

    # 0.1 * 10 is 1 after rounding; the exact product is 1 + 5.55e-17.
    expected = 5.551115123125783e-17 if rmath.FUSED else 0.0
    assert fused(0.1, 10.0, -1.0) == expected
    assert fused(2.0, 3.0, 4.0) == 10.0


def test_lgamma_and_lfactorial(r):
    assert rmath.lgamma(1.0) == 0.0
    assert rmath.lgamma(0.0) == math.inf
    same(rmath.lgamma(np.array([0.5, 1, 3.5, 11, 150, 1e6])), r("rmath_lgamma"))
    same([rmath.lgamma(3.5)], r("rmath_lgamma")[2:3])
    same(rmath.lfactorial(np.arange(9)), r("rmath_lfactorial"))
    assert rmath.lfactorial(1) == 0.0


def test_log_exp_and_sqrt_answer_every_double(r):
    assert rmath.log(0.0) == -math.inf
    assert math.isnan(rmath.log(-1.0))
    assert rmath.exp(1000.0) == math.inf
    assert rmath.exp(0.0) == 1.0
    assert math.isnan(rmath.sqrt(-1.0)) and rmath.sqrt(4.0) == 2.0
    same(rmath.log(DECIMALS), r("rmath_log"))
    same(rmath.exp(-DECIMALS), r("rmath_exp"))
    same([rmath.sqrt(value) for value in DECIMALS], r("rmath_sqrt"))


def test_powers_as_r(r):
    powers = [rmath.r_pow(10, -10), rmath.r_pow(10, -4), rmath.r_pow(2.5, 3.7),
              rmath.r_pow(1.7, 2), rmath.r_pow(0.3, -6), rmath.r_pow(7, 0), rmath.r_pow(0, 2),
              rmath.r_pow(1, 55), rmath.r_pow(3, 5)]
    same(powers, r("rmath_pow"))
    assert rmath.r_pow(0.0, -1.0) == math.inf
    assert rmath.r_pow(1.7, 2) == 1.7 * 1.7
    assert rmath.pow_di(2.0, 10) == 1024.0 and rmath.pow_di(2.0, -2) == 0.25
