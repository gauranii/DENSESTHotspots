# -----------------------------------------------------------------------------
# File:    scripts/webtableS3_bias_binomial.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Web Table S3.
# Purpose: Bias of the estimates of eta, lambda, theta, pi0 and the cutoff C
#          (mean and standard deviation over replications)
#          for the settings ZIP2, ZIGP3, ZIGP4: binomial alternative, pi0 = 0.35,
#          true cutoff C = 5, four null models, two cutoff criteria.
# Inputs:  archive/simulation/<setting>_C5.csv.gz (all replications of the
#          three settings); none with --rerun (simulated data).
# Outputs: results/WebTableS3_bias_binomial-pi0.35_C5.csv; with --rerun also
#          results/simulation/<setting>_C5/runNNNN.csv.
# Usage:   Rscript scripts/webtableS3_bias_binomial.R [--rerun] [--runs=1000] [--cores=4]
#          --rerun: run the replications here instead of reading the archive.
#          --runs: N (runs 1 to N), a:b or a comma list; finished runs are kept.
# Source:  FINAL R CODES.R (FDRCalculation; calls with C <- 5).
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(list(runs = as.character(PAPER$replications), cores = 4))

simulation_table("WebTableS3", "bias", f1 = "binomial", pi0 = 0.35, C = 5,
                 opts, PAPER, PATHS)
