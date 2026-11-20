# -----------------------------------------------------------------------------
# File:    scripts/fig1_simulation_histograms.R
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
# Usage:   Rscript scripts/fig1_simulation_histograms.R
# Source:  PLOTS.R of 30 June 2016 (data3, data1, data2), PLOTS.R of 24 January 2017
#          (bar6, bar4, bar5).
# -----------------------------------------------------------------------------

source("R/setup.R")

save_simulation_histograms("Fig1", PAPER, PATHS)
