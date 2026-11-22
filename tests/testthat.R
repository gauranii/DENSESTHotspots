# -----------------------------------------------------------------------------
# File:    tests/testthat.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. All figures and tables.
# Purpose: Run the unit tests of the functions in R/: every file
#          tests/testthat/test-*.R, on small synthetic inputs.
# Inputs:  tests/testthat/test-*.R; the tests do not read data/ (one reads
#          the first archived replication) and write only to temporary
#          folders.
# Outputs: prints one line per test file; exit status 1 if a test fails.
# Usage:   Rscript tests/testthat.R      # from anywhere, inside `nix develop`
#          nix run .#test                # the same, with the pinned R
# -----------------------------------------------------------------------------

library(testthat)

runner <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
tests_dir <- file.path(dirname(normalizePath(runner)), "testthat")

# R/setup.R finds the repository from the file it is sourced from, so it is
# sourced here, once, and the test files use its functions.
source(file.path(tests_dir, "..", "..", "R", "setup.R"))

test_dir(tests_dir, reporter = "summary", stop_on_failure = TRUE, stop_on_warning = FALSE)
