# -----------------------------------------------------------------------------
# File:    python/tests/test_fdr.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Tests of fdr.py: the local FDR, the three procedures on data with
#          and without labels, against known values and against R, double
#          for double.
# Usage:   python -m pytest python/tests/test_fdr.py
# -----------------------------------------------------------------------------
import math

import numpy as np

from conftest import same
from inputs import COUNTS, FIT, NONNULL, NULL
from zinull import fdr
from zinull.cutoff import select_null_fits

PART = np.array([0, 1, 1, 5, 20.0])


def test_inverse_frequency(r):
    same(fdr.inverse_frequency(PART, COUNTS), r("fdr_inverse_frequency"))
    same(fdr.inverse_frequency([0, 2], [0, 0, 0, 2]), [4 / 3, 4])


def test_local_fdr(r):
    same(fdr.local_fdr(PART, fdr.inverse_frequency(PART, COUNTS), FIT), r("fdr_local_fdr"))
    assert len(fdr.local_fdr([], [], FIT)) == 0


def test_count_below_and_count_between(r):
    lfdr = fdr.local_fdr(COUNTS, fdr.inverse_frequency(COUNTS, COUNTS), FIT)
    assert fdr.count_below([0.01, np.nan, 0.2, 0.04, np.nan], 0.05) == 2
    same([2, fdr.count_below(lfdr, 0.05)], r("fdr_count_below"))
    same([fdr.count_between(COUNTS, lfdr, FIT, 0.9)], r("fdr_count_between"))


def test_rejections_without_labels(r):
    counts = [fdr.one_stage_count(COUNTS, FIT, 0.05), fdr.two_stage_count(COUNTS, FIT, 0.05),
              fdr.storey_count(COUNTS, FIT, 0.05)]
    same(counts, r("fdr_counts"))
    fitted = [procedure(COUNTS, fit, 0.1) for _, fit in select_null_fits(COUNTS).iterrows()
              for procedure in (fdr.one_stage_count, fdr.two_stage_count, fdr.storey_count)]
    same(fitted, r("fdr_counts_fitted"))


def test_storey_p_values(r):
    p_values = fdr.storey_p_values(8, FIT)
    assert p_values[0] == 1 and (np.diff(p_values) <= 0).all()
    same(p_values, r("fdr_storey_p_values"))


def test_storey_rejected(r):
    rejected = fdr.storey_rejected([1, 0.5, 0.04, 0.001, 0.0002], [50, 20, 10, 5, 15], 0.8, 0.05)
    same(rejected + 1, r("fdr_storey_rejected"))  # R counts positions from 1
    assert len(fdr.storey_rejected([1, 0.9], [5, 5], 1, 0.05)) == r("fdr_storey_rejected_none")[0]


def test_discovery_rates(r):
    rates = np.concatenate((fdr.discovery_rates(2, 8, 10), fdr.discovery_rates(0, 0, 10),
                            fdr.discovery_rates(1, 4, 10, R=7)))
    same(rates, r("fdr_discovery_rates"))
    same(fdr.discovery_rates(2, 8, 10), [10, 0.2, 0.8])
    assert math.isnan(fdr.discovery_rates(0, 0, 10)[1])  # no rejections: no FDR


def test_rates_with_labels(r):
    local = fdr.local_fdr_rates(np.sort(NULL), np.sort(NONNULL), FIT, 0.05)
    assert list(local) == ["two_stage_R", "two_stage_FDR", "two_stage_TPR", "one_stage_R",
                           "one_stage_FDR", "one_stage_TPR"]
    same(list(local.values()), r("fdr_local_fdr_rates"))
    storey = fdr.storey_rates(NULL, NONNULL, FIT, 0.05)
    assert list(storey) == ["storey_R", "storey_FDR", "storey_TPR"]
    same(list(storey.values()), r("fdr_storey_rates"))
    fitted = []
    for _, fit in select_null_fits(COUNTS).iterrows():
        fitted += list(fdr.local_fdr_rates(np.sort(NULL), np.sort(NONNULL), fit, 0.05).values())
        fitted += list(fdr.storey_rates(NULL, NONNULL, fit, 0.05).values())
    same(fitted, r("fdr_rates_fitted"))
