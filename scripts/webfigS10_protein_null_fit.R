# -----------------------------------------------------------------------------
# File:    scripts/webfigS10_protein_null_fit.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. Web Figure S10.
# Purpose: Protein domain cd00204: observed share of positions with each count
#          (fHat) and pi0 times the fitted ZIGP and ZIP null probabilities up
#          to the chosen cutoff, for the cutoff criteria C1 (A) and C2 (B).
# Inputs:  data/realdata.txt
# Outputs: figures/WebFigS10{A,B}_protein_null-fit_cd00204_{C1,C2}.pdf
# Usage:   Rscript scripts/webfigS10_protein_null_fit.R
# Source:  FINALREAL 0418.R of 18 April 2017 (plot); null fits as in
#          FINAL R CODES.R (RealDataFDR).
# -----------------------------------------------------------------------------

source("R/setup.R")

save_null_fit_plots("WebFigS10", PATHS)
