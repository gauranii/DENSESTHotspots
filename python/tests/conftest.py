# -----------------------------------------------------------------------------
# File:    python/tests/conftest.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Shared set-up of the unit tests: the import path, the values of R
#          stored by python/tests/make_r_values.R, and the comparison of
#          doubles without tolerance.
# Usage:   python -m pytest python/tests
# -----------------------------------------------------------------------------
import json
import os
import sys

import numpy as np
import pytest

TESTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(TESTS, ".."))
sys.path.insert(0, TESTS)

from zinull import rbase  # noqa: E402

SPECIAL = {"NA": rbase.NA, "NaN": float("nan"), "Inf": float("inf"), "-Inf": float("-inf")}


class RValues:
    """The values R computed, by name."""

    def __init__(self, file):
        with open(file) as f:
            self.text = json.load(f)

    def __call__(self, name):
        """The doubles stored under a name (NA and NaN are both NaN)."""
        return np.array([SPECIAL[item] if item in SPECIAL else float(item)
                         for item in self.text[name]], dtype=np.float64)

    def lines(self, name):
        """The text stored under a name (lines of a file, names)."""
        return list(self.text[name])


@pytest.fixture(scope="session")
def r():
    return RValues(os.path.join(TESTS, "r_values.json"))


def same(actual, expected):
    """Assert that two sets of doubles are equal, double for double (a
    missing value equals a missing value). No tolerance."""
    actual = np.asarray(actual, dtype=np.float64).ravel()
    expected = np.asarray(expected, dtype=np.float64).ravel()
    assert actual.shape == expected.shape
    equal = (actual == expected) | (np.isnan(actual) & np.isnan(expected))
    assert equal.all(), (f"first difference at element {int(np.argmin(equal))}: "
                         f"{actual[np.argmin(equal)]!r} and {expected[np.argmin(equal)]!r}")


def numbers(table):
    """The numeric columns of a DataFrame, one after the other (as
    make_r_values.R stores a data frame)."""
    columns = [table[name].to_numpy(dtype=np.float64) for name in table.columns
               if table[name].dtype.kind in "fiu"]
    return np.concatenate(columns)


@pytest.fixture
def data_file():
    """The counts of two simulated "domains", in the layout of data/realdata.txt."""
    return os.path.join(TESTS, "data", "small_domains.txt")
