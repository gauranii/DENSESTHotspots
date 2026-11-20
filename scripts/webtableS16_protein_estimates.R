# -----------------------------------------------------------------------------
# File:    scripts/webtableS16_protein_estimates.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Web Table S16.
# Purpose: Estimates of eta, lambda, theta, pi0, the cutoff C and the bound D in six
#          protein domains, for four null models and two cutoff criteria.
# Inputs:  data/realdata.txt.
# Outputs: results/WebTableS16_protein_estimates.csv.
# Usage:   Rscript scripts/webtableS16_protein_estimates.R
# Source:  FINAL R CODES.R (RealDataFDR for k = 1, ..., 6; the Estimates matrices).
# -----------------------------------------------------------------------------

source("R/setup.R")

table <- protein_estimates(PATHS$data)
write.csv(table, output_file(PATHS$results, "WebTableS16_protein_estimates.csv"), row.names = FALSE)
print(table, row.names = FALSE)
