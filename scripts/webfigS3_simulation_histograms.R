# -----------------------------------------------------------------------------
# File:    scripts/webfigS3_simulation_histograms.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. Web Figure S3.
# Purpose: Histograms of three simulated data sets of 1000 counts in which the
#          non-null counts (black, geometric signal) overlap the null counts
#          (gray): ZIGP (A, B) and ZIP (C) null counts, 80 % null counts.
# Inputs:  None (simulated data, seed 889191 + 111).
# Outputs: figures/FigS3{A,B,C}_simulation_histogram_{MZIGP1,MZIGP2,MZIP1}-geometric.pdf
# Usage:   Rscript scripts/webfigS3_simulation_histograms.R
# Source:  Last block of the R history of 24 February 2017 (mixed1, mixed2, mixed3).
# -----------------------------------------------------------------------------

source("R/setup.R")

save_simulation_histograms("FigS3", PAPER, PATHS)
