# -----------------------------------------------------------------------------
# File:    scripts/webtableS14_simulation_geometric.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Web Table S14.
# Purpose: Number of rejections, FDR and TPR of the two-stage, one-stage and
#          Storey procedures (mean and standard deviation over replications)
#          for the settings ZIP3, ZIGP5, ZIGP6: geometric alternative, pi0 = 0.85,
#          true cutoff C = 1 (contaminated null), four null models, two cutoff criteria.
# Inputs:  archive/simulation/<setting>_C1.csv.gz (all replications of the
#          three settings); none with --rerun (simulated data).
# Outputs: results/WebTableS14_simulation_geometric-pi0.85_C1.csv; with --rerun also
#          results/simulation/<setting>_C1/runNNNN.csv.
# Usage:   Rscript scripts/webtableS14_simulation_geometric.R [--rerun] [--runs=1000] [--cores=4]
#          --rerun: run the replications here instead of reading the archive.
#          --runs: N (runs 1 to N), a:b or a comma list; finished runs are kept.
# Source:  FINAL R CODES.R (FDRCalculation; calls with C <- 1, "CONTAMINATION").
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(list(runs = as.character(PAPER$replications), cores = 4))

simulation_table("WebTableS14", "simulation", f1 = "geometric", pi0 = 0.85, C = 1,
                 opts, PAPER, PATHS)
