# -----------------------------------------------------------------------------
# File:    python/scripts/fig2_protein_histograms.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. Figure 2.
# Purpose: Histograms of the counts at the positions of six protein domains.
# Inputs:  data/realdata.txt
# Outputs: figures/Fig2{A,...,F}_protein_histogram_<domain>.pdf, domains cd00031,
#          cd00180, cd00204, cd00882, cd00054, pfam00001
# Usage:   python python/scripts/fig2_protein_histograms.py
#          Options and outputs as in the R script of the same name;
#          --help prints this header.
# Source:  scripts/fig2_protein_histograms.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from zinull.setup import PATHS, show_help  # noqa: E402

show_help()

from zinull.figure_data import PROTEIN_DOMAINS  # noqa: E402
from zinull.plots import save_protein_histograms  # noqa: E402

save_protein_histograms("Fig2", PROTEIN_DOMAINS, PATHS)
