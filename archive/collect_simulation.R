# -----------------------------------------------------------------------------
# File:    archive/collect_simulation.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Tables 1 and 2, Web Tables S2 to S15.
# Purpose: Combine the finished replications of every simulation setting into
#          the one file per setting kept in archive/. Numbers are copied as
#          text (17 significant digits), so no digit changes.
# Inputs:  results/simulation/<setting>_C<C>/runNNNN.csv.
# Outputs: archive/simulation/<setting>_C<C>.csv.gz.
# Usage:   Rscript archive/collect_simulation.R
# -----------------------------------------------------------------------------

source("R/setup.R")

for (folder in list.dirs(file.path(PATHS$results, "simulation"), recursive = FALSE)) {
  files <- sort(list.files(folder, pattern = "^run[0-9]+\\.csv$", full.names = TRUE))
  runs <- do.call(rbind, lapply(files, read.csv, colClasses = "character"))
  archived <- output_file(file.path(PATHS$archive, "simulation"),
                          paste0(basename(folder), ".csv.gz"))
  write.csv(runs, gzfile(archived), row.names = FALSE, quote = FALSE)
  log_message("%s: %d replications", basename(archived), length(files))
}
