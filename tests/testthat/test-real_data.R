# -----------------------------------------------------------------------------
# File:    tests/testthat/test-real_data.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Unit tests of R/real_data.R on two made-up protein domains: the
#          counts read from a file in the layout of data/realdata.txt, the
#          null fits of a domain, and the tables of the numbers of rejected
#          positions and of the estimates against their parts.
# Inputs:  None (a made-up file of domain counts in a temporary folder).
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

DATA_FILE <- write_domain_file(new_folder())
REJECTIONS <- protein_rejections(DATA_FILE, 0.10)

test_that("read_domains gives the counts of each domain without the padding", {
  domains <- read_domains(DATA_FILE)
  expect_named(domains, c("made1", "made2"))
  expect_equal(lengths(domains), c(made1 = 305, made2 = 250))
  expect_equal(domains$made1, DOMAIN_COUNTS$made1)
  expect_equal(domains$made2, DOMAIN_COUNTS$made2)
})

test_that("read_domains stops with a message when the data file is not there", {
  expect_error(read_domains(file.path(new_folder(), "realdata.txt")),
               "not included in this repository")
})

test_that("fit_domain is the null fits with D at most the largest count", {
  x <- DOMAIN_COUNTS$made2
  fits <- fit_domain(x)
  expect_equal(fits, select_null_fits(x, cap = max(x)))
  expect_true(all(fits$D <= max(x)))
  expect_equal(dim(fits), c(8, 8))
})

test_that("protein_rejections has one row per domain, criterion and procedure", {
  expect_named(REJECTIONS, c("domain", "cutoff", "procedure", MODELS))
  expect_equal(nrow(REJECTIONS), 12)
  expect_equal(REJECTIONS$domain, rep(c("made1", "made2"), each = 6))
  expect_equal(REJECTIONS$cutoff, rep(rep(CUTOFFS, each = 3), 2))
  expect_equal(REJECTIONS$procedure, rep(c("one_stage", "two_stage", "storey"), 4))
  counts <- as.matrix(REJECTIONS[MODELS])
  expect_true(all(counts >= 0 & counts <= 305 & counts == round(counts)))
})

test_that("protein_rejections holds the counts of the three procedures", {
  x <- DOMAIN_COUNTS$made1
  fits <- fit_domain(x)
  expect_equal(fits, MADE_FITS)   # no bound above the largest count, 25
  cell <- function(criterion, procedure, model) {
    REJECTIONS[REJECTIONS$domain == "made1" & REJECTIONS$cutoff == criterion &
                 REJECTIONS$procedure == procedure, model]
  }
  for (i in seq_len(nrow(fits))) {
    fit <- fits[i, ]
    expect_equal(cell(fit$cutoff, "one_stage", fit$model), one_stage_count(x, fit, 0.10))
    expect_equal(cell(fit$cutoff, "two_stage", fit$model), two_stage_count(x, fit, 0.10))
    expect_equal(cell(fit$cutoff, "storey", fit$model), storey_count(x, fit, 0.10))
  }
  # the five large counts (15 to 25) are found with the zero-inflated Poisson null
  expect_equal(cell("C1", "one_stage", "ZIP"), 5)
  expect_gte(cell("C1", "two_stage", "ZIP"), 5)
})

test_that("protein_estimates has the null fits of every domain, absent parameters NA", {
  estimates <- protein_estimates(DATA_FILE)
  expect_named(estimates, c("domain", "cutoff", "model", "eta", "lambda", "theta", "pi0",
                            "C", "D"))
  expect_equal(estimates$domain, rep(c("made1", "made2"), each = 8))
  expect_equal(estimates$model, rep(MODELS, 4))
  expect_equal(plain(estimates[1:8, -1]), plain(mask_absent_parameters(MADE_FITS)))
  expect_equal(is.na(estimates$eta), estimates$model %in% c("GP", "P"))
  expect_equal(is.na(estimates$theta), estimates$model %in% c("ZIP", "P"))
  expect_true(all(estimates$D[9:16] <= max(DOMAIN_COUNTS$made2)))
})
