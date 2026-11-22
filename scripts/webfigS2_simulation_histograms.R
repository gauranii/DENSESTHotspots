# -----------------------------------------------------------------------------
# File:    scripts/webfigS2_simulation_histograms.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. Web Figure S2.
# Purpose: Histograms of three simulated data sets of 1000 counts with a
#          binomial signal, null counts over non-null counts (black): ZIGP
#          (A, B) and ZIP (C) null counts.
# Inputs:  None (simulated data, seed 889191 + 111).
# Outputs: figures/FigS2{A,B,C}_simulation_histogram_{ZIGP7,ZIGP8,ZIP4}-binomial.pdf
# Usage:   Rscript scripts/webfigS2_simulation_histograms.R
# Source:  PLOTS.R of 24 January 2017 (bar1, bar2, bar3).
# -----------------------------------------------------------------------------

source("R/setup.R")

save_simulation_histograms("FigS2", PAPER, PATHS)
