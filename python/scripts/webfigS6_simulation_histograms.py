# -----------------------------------------------------------------------------
# File:    python/scripts/webfigS6_simulation_histograms.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. Web Figure S6.
# Purpose: Histograms of three simulated data sets of 1000 counts in which the
#          non-null counts (black, binomial signal) overlap the null counts
#          (gray): ZIGP (A, B) and ZIP (C) null counts, 80 % null counts.
# Inputs:  None (simulated data, seed 889191 + 111).
# Outputs: figures/FigS6{A,B,C}_simulation_histogram_{MZIGP7,MZIGP8,MZIP4}-binomial.pdf
# Usage:   python python/scripts/webfigS6_simulation_histograms.py
#          Options and outputs as in the R script of the same name;
#          --help prints this header.
# Source:  scripts/webfigS6_simulation_histograms.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from zinull.setup import PAPER, PATHS, show_help  # noqa: E402

show_help()

from zinull.plots import save_simulation_histograms  # noqa: E402

save_simulation_histograms("FigS6", PAPER, PATHS)
