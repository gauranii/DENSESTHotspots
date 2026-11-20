# -----------------------------------------------------------------------------
# File:    scripts/webfigS1_simulation_histograms.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. Web Figure S1.
# Purpose: Histograms of three simulated data sets of 1000 counts with a
#          geometric signal, null counts over non-null counts (black): ZIGP
#          (A, B) and ZIP (C) null counts.
# Inputs:  None (simulated data, seed 889191 + 111).
# Outputs: figures/WebFigS1{A,B,C}_simulation_histogram_{ZIGP5,ZIGP6,ZIP3}-geometric.pdf
# Usage:   Rscript scripts/webfigS1_simulation_histograms.R
# Source:  PLOTS.R of 30 June 2016 (data5, data6, data7).
# -----------------------------------------------------------------------------

source("R/setup.R")

save_simulation_histograms("WebFigS1", PAPER, PATHS)
