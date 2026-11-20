# -----------------------------------------------------------------------------
# File:    scripts/fig2_protein_histograms.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. Figure 2.
# Purpose: Histograms of the counts at the positions of six protein domains.
# Inputs:  data/realdata.txt
# Outputs: figures/Fig2{A,...,F}_protein_histogram_<domain>.pdf, domains cd00031,
#          cd00180, cd00204, cd00882, cd00054, pfam00001
# Usage:   Rscript scripts/fig2_protein_histograms.R
# Source:  realdata8.Rhistory (lines 393-404).
# -----------------------------------------------------------------------------

source("R/setup.R")

save_protein_histograms("Fig2", PROTEIN_DOMAINS, PATHS)
