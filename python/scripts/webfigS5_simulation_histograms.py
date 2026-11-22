# -----------------------------------------------------------------------------
# File:    python/scripts/webfigS5_simulation_histograms.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. Web Figure S5.
# Purpose: Histograms of three simulated data sets of 1000 counts in which the
#          non-null counts (black, geometric signal) overlap the null counts
#          (gray): ZIGP (A, B) and ZIP (C) null counts, 85 % null counts.
# Inputs:  None (simulated data, seed 889191 + 111).
# Outputs: figures/FigS5{A,B,C}_simulation_histogram_{MZIGP5,MZIGP6,MZIP3}-geometric.pdf
# Usage:   python python/scripts/webfigS5_simulation_histograms.py
#          Options and outputs as scripts/webfigS5_simulation_histograms.R, which this file follows;
#          --help prints this header.
# Source:  scripts/webfigS5_simulation_histograms.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from zinull.setup import PAPER, PATHS, show_help  # noqa: E402

show_help()

from zinull.plots import save_simulation_histograms  # noqa: E402

save_simulation_histograms("FigS5", PAPER, PATHS)
