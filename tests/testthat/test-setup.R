# -----------------------------------------------------------------------------
# File:    tests/testthat/test-setup.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Unit tests of R/setup.R: the settings of the paper, the project
#          paths, the message for the missing data file, the output folder,
#          the command line parser and the index sets.
# Inputs:  None.
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

test_that("PAPER holds the settings of the paper", {
  expect_equal(PAPER$seed, 889191)
  expect_equal(PAPER$figure_seed, 889302)
  expect_equal(c(PAPER$n, PAPER$replications), c(1000, 1000))
  expect_equal(c(PAPER$fdr_level, PAPER$fdr_level_real), c(0.05, 0.10))
  expect_true(file.exists(file.path(REPO_ROOT, "R", "setup.R")))
})

test_that("project_paths puts archive under the root and the outputs under output", {
  paths <- withr::with_envvar(c(ZINULL_DATA = NA), project_paths("/repo", output = "/out"))
  expect_equal(paths, list(data = "/repo/data/realdata.txt", archive = "/repo/archive",
                           results = "/out/results", figures = "/out/figures"))
})

test_that("project_paths reads ZINULL_DATA and ZINULL_OUTPUT", {
  paths <- withr::with_envvar(c(ZINULL_DATA = "/elsewhere/counts.txt", ZINULL_OUTPUT = "/work"),
                              project_paths("/repo"))
  expect_equal(paths$data, "/elsewhere/counts.txt")
  expect_equal(paths$results, "/work/results")
  expect_equal(paths$archive, "/repo/archive")
  paths <- withr::with_envvar(c(ZINULL_DATA = NA, ZINULL_OUTPUT = NA), project_paths("/repo"))
  expect_equal(paths$figures, "/repo/figures")
})

test_that("require_data returns the path of a file that is there and stops otherwise", {
  file <- tempfile()
  writeLines("id", file)
  expect_equal(require_data(file), file)
  expect_invisible(require_data(file))
  expect_error(require_data(file.path(new_folder(), "realdata.txt")),
               "not included in this repository")
})

test_that("output_file creates the folder and returns the path in it", {
  folder <- file.path(new_folder(), "results", "deep")
  expect_false(dir.exists(folder))
  expect_equal(output_file(folder, "Table.csv"), file.path(folder, "Table.csv"))
  expect_true(dir.exists(folder))
  expect_equal(output_file(folder, "Other.csv"), file.path(folder, "Other.csv"))  # folder kept
})

test_that("parse_args keeps the defaults and adds rerun", {
  expect_equal(parse_args(list(runs = "1000", cores = 4), args = character(0)),
               list(rerun = FALSE, runs = "1000", cores = 4))
})

test_that("parse_args converts a value to the type of its default", {
  opts <- parse_args(list(runs = "1000", cores = 4, level = 0.05, draft = FALSE),
                     args = c("--runs=1:10", "--cores=2", "--rerun", "--level=0.1", "--draft"))
  expect_equal(opts, list(rerun = TRUE, runs = "1:10", cores = 2, level = 0.1, draft = TRUE))
  expect_type(opts$cores, "double")
  expect_type(opts$runs, "character")
})

test_that("parse_args stops at arguments it does not know or cannot read", {
  defaults <- list(runs = "1000", cores = 4)
  expect_error(parse_args(defaults, args = "runs=3"), "must look like")
  expect_error(parse_args(defaults, args = "-cores=3"), "must look like")
  expect_error(parse_args(defaults, args = "--seed=3"), "unknown argument --seed")
  expect_error(parse_args(defaults, args = "--cores"), "--cores needs a value")
  expect_error(parse_args(defaults, args = "--rerun=yes"), "--rerun takes no value")
})

test_that("parse_index_set reads a range, a list and a single number", {
  expect_equal(parse_index_set("3:7"), 3:7)
  expect_equal(parse_index_set("2,5,9"), c(2L, 5L, 9L))
  expect_equal(parse_index_set("4"), 4L)
  expect_equal(parse_index_set("4", count = TRUE), 1:4)
  expect_equal(parse_index_set("2,5", count = TRUE), c(2L, 5L))
  expect_equal(parse_index_set("7:7", count = TRUE), 7L)
})

test_that("log_message prints the time and the formatted message", {
  printed <- capture.output(value <- log_message("wrote %d files to %s", 3, "results/"))
  expect_null(value)
  expect_match(printed, "^\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2} - wrote 3 files to results/ $")
})
