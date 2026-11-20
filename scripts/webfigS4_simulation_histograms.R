# -----------------------------------------------------------------------------
# File:    scripts/webfigS4_simulation_histograms.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. Web Figure S4.
# Purpose: Histograms of three simulated data sets of 1000 counts in which the
#          non-null counts (black, binomial signal) overlap the null counts
#          (gray): ZIGP (A, B) and ZIP (C) null counts, 35 % null counts.
# Inputs:  None (simulated data, seed 889191 + 111).
# Outputs: figures/WebFigS4{A,B,C}_simulation_histogram_{MZIGP3,MZIGP4,MZIP2}-binomial.pdf
# Usage:   Rscript scripts/webfigS4_simulation_histograms.R
# Source:  Last block of the R history of 24 February 2017 (mixed4, mixed5, mixed6).
# -----------------------------------------------------------------------------

source("R/setup.R")

save_simulation_histograms("WebFigS4", PAPER, PATHS)
