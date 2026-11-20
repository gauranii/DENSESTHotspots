# -----------------------------------------------------------------------------
# File:    scripts/table1_simulation_geometric.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Table 1.
# Purpose: Number of rejections, FDR and TPR of the two-stage, one-stage and
#          Storey procedures (mean and standard deviation over replications)
#          for the settings ZIP1, ZIGP1, ZIGP2: geometric alternative, pi0 = 0.80,
#          true cutoff C = 5, four null models, two cutoff criteria.
# Inputs:  archive/simulation/<setting>_C5.csv.gz (all replications of the
#          three settings); none with --rerun (simulated data).
# Outputs: results/Table1_simulation_geometric-pi0.80_C5.csv; with --rerun also
#          results/simulation/<setting>_C5/runNNNN.csv.
# Usage:   Rscript scripts/table1_simulation_geometric.R [--rerun] [--runs=1000] [--cores=4]
#          --rerun: run the replications here instead of reading the archive.
#          --runs: N (runs 1 to N), a:b or a comma list; finished runs are kept.
# Source:  FINAL R CODES.R (FDRCalculation; calls with C <- 5).
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(list(runs = as.character(PAPER$replications), cores = 4))

simulation_table("Table1", "simulation", f1 = "geometric", pi0 = 0.80, C = 5,
                 opts, PAPER, PATHS)
