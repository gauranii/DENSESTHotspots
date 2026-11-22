# -----------------------------------------------------------------------------
# File:    python/scripts/fig1_simulation_histograms.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. Figure 1.
# Purpose: Histograms of six simulated data sets of 1000 counts, null counts over
#          non-null counts (black): geometric signal with ZIP (A) and ZIGP
#          (B, C) null counts, binomial signal with ZIP (D) and ZIGP (E, F)
#          null counts.
# Inputs:  None (simulated data, seed 889191 + 111).
# Outputs: figures/Fig1{A,B,C}_simulation_histogram_{ZIP1,ZIGP1,ZIGP2}-geometric.pdf,
#          figures/Fig1{D,E,F}_simulation_histogram_{ZIP2,ZIGP3,ZIGP4}-binomial.pdf
# Usage:   python python/scripts/fig1_simulation_histograms.py
#          Options and outputs as in the R script of the same name;
#          --help prints this header.
# Source:  scripts/fig1_simulation_histograms.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from zinull.setup import PAPER, PATHS, show_help  # noqa: E402

show_help()

from zinull.plots import save_simulation_histograms  # noqa: E402

save_simulation_histograms("Fig1", PAPER, PATHS)
