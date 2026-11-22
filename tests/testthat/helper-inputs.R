# -----------------------------------------------------------------------------
# File:    tests/testthat/helper-inputs.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Inputs shared by the unit tests (testthat loads this file before
#          the tests): temporary folders, made-up counts with a known null,
#          a made-up file of protein domain counts, made-up replications.
# Inputs:  None.
# Outputs: None.
# Usage:   loaded by testthat; see tests/testthat.R
# -----------------------------------------------------------------------------

if (!exists("REPO_ROOT")) {
  stop("run the tests with tests/testthat.R, which sources R/setup.R first")
}

#' A new empty temporary folder.
#'
#' @return Its path.
new_folder <- function() {
  folder <- tempfile("zinull-test-")
  dir.create(folder)
  folder
}

#' Run code without what it prints.
#'
#' @param code Expression that prints tables or time-stamped messages.
#' @return The value of `code`.
quietly <- function(code) {
  capture.output(value <- code)
  value
}

#' A data frame without its row names, to compare the rows of two tables.
#'
#' @param table Data frame.
#' @return `table` with the row names 1, 2, ...
plain <- function(table) {
  rownames(table) <- NULL
  table
}

#' TRUE if a file starts as a PDF file does and is not empty of marks.
#'
#' @param file Path of the file.
#' @return TRUE or FALSE.
is_pdf <- function(file) {
  file.exists(file) && file.size(file) > 1000 &&
    identical(readBin(file, "raw", 4), charToRaw("%PDF"))
}

#' Made-up counts: 300 draws of a zero-inflated Poisson null (eta = 0.4,
#' lambda = 1.5) and five large counts.
#'
#' @param seed Seed of the null draws.
#' @return Vector of 305 counts.
made_up_counts <- function(seed) {
  c(generate_zigp(300, mu = 1.5, phi = 1, omega = 0.4, seed = seed), 15, 16, 18, 20, 25)
}

# Two made-up protein domains and the null fits of the first one.
DOMAIN_COUNTS <- list(made1 = made_up_counts(11), made2 = made_up_counts(13)[1:250])
MADE_FITS <- select_null_fits(DOMAIN_COUNTS$made1)

#' Write the made-up domains in the layout of data/realdata.txt: a column
#' `id` and one column per domain, tab separated, the shorter ones NA padded.
#'
#' @param folder Folder of the file.
#' @return The path of the file.
write_domain_file <- function(folder) {
  rows <- max(lengths(DOMAIN_COUNTS))
  table <- data.frame(id = seq_len(rows),
                      lapply(DOMAIN_COUNTS, function(x) c(x, rep(NA, rows - length(x)))))
  file <- file.path(folder, "realdata.txt")
  write.table(table, file, sep = "\t", quote = FALSE, row.names = FALSE)
  file
}

# A null fit with round numbers and 100 counts, for the cases worked by hand.
# Local FDR of each count, 0.9 f0(x) 100 / frequency:
#   0: 0.80   1: 0.72   2: 1.51   6: 0.19   8: 0.0038   9: 0.0013   12: 2e-6
HAND_FIT <- list(eta = 0.4, lambda = 1.5, theta = 0, pi0 = 0.9, C = 2, D = 10)
HAND_NULL <- c(rep(0, 60), rep(1, 25), rep(2, 9), 6, 8)
HAND_NONNULL <- c(8, 9, 12, 12)
HAND_COUNTS <- c(HAND_NULL, HAND_NONNULL)

#' Made-up replications of one setting in the layout of
#' `simulation_replication()`: every column is a simple function of the run,
#' the criterion and the model.
#'
#' @param runs Run numbers.
#' @return Data frame with one row per run, criterion and model.
made_up_replications <- function(runs) {
  cells <- expand.grid(model = MODELS, cutoff = CUTOFFS, run = runs, stringsAsFactors = FALSE)
  base <- cells$run + 10 * match(cells$model, MODELS) + 100 * match(cells$cutoff, CUTOFFS)
  columns <- c(TABLE_COLUMNS$simulation, "eta", "lambda", "theta", "pi0", "C", "D",
               TABLE_COLUMNS$bias)
  values <- vapply(seq_along(columns), function(j) base + j / 7, numeric(nrow(cells)))
  colnames(values) <- columns
  data.frame(cells[c("run", "cutoff", "model")], values)
}
