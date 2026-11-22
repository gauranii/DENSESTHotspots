# -----------------------------------------------------------------------------
# File:    python/tests/test_rbase.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Tests of rbase.py: known values of every function, and the same
#          doubles and the same text as R.
# Usage:   python -m pytest python/tests/test_rbase.py
# -----------------------------------------------------------------------------
import math

import numpy as np
import pandas as pd
import pytest

from conftest import same
from inputs import COUNTS, DECIMALS, WHOLE
from zinull import rbase


def test_missing_values():
    assert math.isnan(rbase.NA) and rbase.is_na(rbase.NA)
    assert not rbase.is_na(float("nan")) and not rbase.is_na(1.0)
    np.testing.assert_array_equal(rbase.is_na(np.array([1.0, rbase.NA, np.nan])),
                                  [False, True, False])


def test_sums_and_means(r):
    assert rbase.r_sum([1, 2, 3.5]) == 6.5
    same([rbase.r_sum(DECIMALS), rbase.r_sum(DECIMALS / 7)], r("rbase_sum"))
    same(rbase.r_cumsum(DECIMALS / 7), r("rbase_cumsum"))
    np.testing.assert_array_equal(rbase.r_cumsum([1, 2, 3]), [1, 3, 6])
    means = [rbase.r_mean(DECIMALS), rbase.r_mean(DECIMALS / 7), rbase.r_mean(COUNTS),
             rbase.r_mean_integer(WHOLE)]
    same(means, r("rbase_mean"))
    assert rbase.r_mean([1, np.nan, 3], na_rm=True) == 2.0
    assert rbase.r_mean_integer([1, 2, 4]) == 7 / 3
    assert math.isnan(rbase.r_mean([], na_rm=True))


def test_variances(r):
    """var(x) and var(x, na.rm = TRUE) add in different orders in R."""
    variances = [rbase.r_var(DECIMALS), rbase.r_var(COUNTS), rbase.r_var(DECIMALS[:9]),
                 rbase.r_var(np.append(DECIMALS, np.nan), na_rm=True),
                 rbase.r_var(COUNTS / 7, na_rm=True)]
    same(variances, r("rbase_var"))
    sds = [rbase.r_sd(DECIMALS), rbase.r_sd(np.append(DECIMALS, np.nan), na_rm=True),
           rbase.r_sd(COUNTS / 7, na_rm=True)]
    same(sds, r("rbase_sd"))
    assert rbase.r_var([1, 2, 3, 4]) == 5 / 3
    assert rbase.is_na(rbase.r_var([1.0])) and rbase.is_na(rbase.r_sd([np.nan], na_rm=True))
    assert math.isnan(rbase.r_var([1.0, np.nan]))


def test_quantile(r):
    assert rbase.quantile([1, 2, 3, 4], 0.5) == 2.5
    assert rbase.quantile([3, 1, 2], 1.0) == 3.0
    quantiles = [rbase.quantile(DECIMALS, 0.75), rbase.quantile(COUNTS, 0.5),
                 rbase.quantile(DECIMALS, 0.31)]
    same(quantiles, r("rbase_quantile"))


def test_smallest_and_largest_with_missing_values():
    assert rbase.r_min(2, 1, 3) == 1 and rbase.r_max(2, 1, 3) == 3
    assert math.isnan(rbase.r_min(1, math.nan)) and math.isnan(rbase.r_max(math.nan, 1))
    same(rbase.pmin(1, [0.5, 2, np.nan]), [0.5, 1, np.nan])
    same(rbase.pmax(0, [-1, 2, np.nan]), [0, 2, np.nan])
    assert rbase.which_min([3, 1, np.nan, 1]) == 1
    assert rbase.which_min([np.nan]) is None and rbase.which_min([]) is None


def test_tabulate_order_and_sequences():
    np.testing.assert_array_equal(rbase.tabulate([1, 3, 3, 7, 0], 4), [1, 0, 2, 0])
    np.testing.assert_array_equal(rbase.count_frequencies_from_zero([0, 0, 2, 5], 3), [2, 0, 1, 0])
    np.testing.assert_array_equal(rbase.r_order([2, 1, np.nan, 1]), [1, 3, 0, 2])
    np.testing.assert_array_equal(rbase.colon(2, 5), [2, 3, 4, 5])
    np.testing.assert_array_equal(rbase.colon(3, 2), [3, 2])  # R counts down
    np.testing.assert_array_equal(rbase.colon(4, 4), [4])


NUMBERS = [0.1, 1 / 3, 200.913, 1e-5, 123456789.123456789, 1e15, 1e-15, 0.00418141098991527,
           -2.5, 100, 1e5, 1e22, 0.1 + 0.2, 5e-324]


def test_format_number_writes_as_r(r):
    assert [rbase.format_number(value) for value in NUMBERS] == r.lines("rbase_format_number")
    assert rbase.format_number(math.nan) == "NA" and rbase.format_number(rbase.NA) == "NA"
    assert rbase.format_number(math.inf) == "Inf" and rbase.format_number(0.0) == "0"


def test_write_csv_writes_the_file_of_r(r, tmp_path):
    table = pd.DataFrame({"name": ["a", "b"], "value": [1 / 3, rbase.NA], "other": [np.nan, 2.0]})
    rbase.write_csv(table, tmp_path / "table.csv")
    assert (tmp_path / "table.csv").read_text().splitlines() == r.lines("rbase_write_csv")


def _square(item):
    if item == 3:
        raise ValueError("three")
    return item * item


@pytest.mark.parametrize("cores", [1, 2])
def test_parallel_map_keeps_the_order_and_returns_errors(cores):
    results = rbase.parallel_map(_square, [1, 2, 3, 4], cores=cores)
    assert results[:2] == [1, 4] and results[3] == 16
    assert isinstance(results[2], ValueError)
