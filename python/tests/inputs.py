# -----------------------------------------------------------------------------
# File:    python/tests/inputs.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Small inputs shared by the tests, the same as at the top of
#          python/tests/make_r_values.R.
# Usage:   from inputs import COUNTS, FIT
# -----------------------------------------------------------------------------
import numpy as np

from zinull.simulation import SIMULATION_SETTINGS

# Counts: a null part near zero and a well-separated part.
COUNTS = np.array([0] * 30 + [1] * 14 + [2] * 8 + [3] * 4 + [4] * 2 + [5]
                  + [14, 15, 15, 16, 16, 16, 17, 17, 18, 20], dtype=np.float64)
WHOLE = COUNTS.astype(np.int64)  # the same counts stored as integers
NULL, NONNULL = COUNTS[:59], COUNTS[59:]
DECIMALS = np.array([0.1, 2.5, 3.7, 1e-3, 12.25, 7.3, 0.6, 5.5, 9.1, 4.4, 8.8, 1.9])
HOLD = np.array([0.4, 1.5, 0.3])  # eta, lambda, theta
FIT = {"eta": 0.4, "lambda": 1.2, "theta": 0.1, "pi0": 0.8, "C": 5, "D": 9}
SMALL = {"seed": 889191, "n": 150, "fdr_level": 0.05}  # a short simulation
SETTING = SIMULATION_SETTINGS.iloc[1]    # ZIGP1
BINOMIAL = SIMULATION_SETTINGS.iloc[4]   # ZIGP3
