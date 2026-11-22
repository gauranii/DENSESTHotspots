# -----------------------------------------------------------------------------
# File:    python/scripts/webfigS9_protein_null_fit.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. Web Figure S9.
# Purpose: Protein domain cd00180: observed share of positions with each count
#          (fHat) and pi0 times the fitted ZIGP and ZIP null probabilities up
#          to the chosen cutoff, for the cutoff criteria C1 (A) and C2 (B).
# Inputs:  data/realdata.txt
# Outputs: figures/FigS9{A,B}_protein_null-fit_cd00180_{C1,C2}.pdf
# Usage:   python python/scripts/webfigS9_protein_null_fit.py
#          Options and outputs as scripts/webfigS9_protein_null_fit.R, which this file follows;
#          --help prints this header.
# Source:  scripts/webfigS9_protein_null_fit.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from zinull.setup import PATHS, show_help  # noqa: E402

show_help()

from zinull.plots import save_null_fit_plots  # noqa: E402

save_null_fit_plots("FigS9", PATHS)
