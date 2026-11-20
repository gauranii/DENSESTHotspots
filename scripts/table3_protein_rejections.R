# -----------------------------------------------------------------------------
# File:    scripts/table3_protein_rejections.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Table 3.
# Purpose: Numbers of positions rejected at FDR level 0.10 in six protein domains by
#          the one-stage, two-stage and Storey procedures, for four null models and
#          two cutoff criteria.
# Inputs:  data/realdata.txt.
# Outputs: results/Table3_protein_rejections.csv.
# Usage:   Rscript scripts/table3_protein_rejections.R
# Source:  FINAL R CODES.R (RealDataFDR for k = 1, ..., 6; the Rejection matrices).
# -----------------------------------------------------------------------------

source("R/setup.R")

table <- protein_rejections(PATHS$data, PAPER$fdr_level_real)
write.csv(table, output_file(PATHS$results, "Table3_protein_rejections.csv"), row.names = FALSE)
print(table, row.names = FALSE)
