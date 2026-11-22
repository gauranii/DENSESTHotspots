# -----------------------------------------------------------------------------
# File:    tests/testthat/test-figure_data.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Unit tests of R/figure_data.R: the settings of the panels, the
#          simulated mixtures against the order of their random numbers and
#          under a seed, the bar heights, the counts of a domain, and the
#          observed and fitted probabilities of the null-fit plots.
# Inputs:  None (a made-up file of domain counts in a temporary folder).
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

test_that("the settings of the figures are those of the paper", {
  expect_equal(c(CUTOFF, BINOMIAL_SIZE), c(5, 250))
  expect_equal(dim(HISTOGRAM_SETTINGS), c(24, 13))
  expect_equal(unique(HISTOGRAM_SETTINGS$figure), c("Fig1", paste0("FigS", 1:6)))
  expect_equal(as.vector(table(HISTOGRAM_SETTINGS$figure)[c("Fig1", "FigS1")]), c(6, 3))
  expect_true(all(HISTOGRAM_SETTINGS$signal %in% c("geometric", "binomial")))
  expect_true(all(HISTOGRAM_SETTINGS$order %in% c("non_null", "null")))
  expect_equal(anyDuplicated(HISTOGRAM_SETTINGS$original), 0)
  expect_equal(PROTEIN_DOMAINS,
               c("cd00031", "cd00180", "cd00204", "cd00882", "cd00054", "pfam00001"))
  expect_named(NULL_FIT_DOMAINS, paste0("FigS", 7:12))
  expect_setequal(NULL_FIT_DOMAINS, PROTEIN_DOMAINS)
})

test_that("draw_null_count is one ZIGP draw in the parametrization of the paper", {
  set.seed(31)
  draws <- vapply(1:200, function(i) draw_null_count(0.4, 1.5, 0.3), numeric(1))
  expect_identical(draws, generate_zigp(200, mu = 1.5 / 0.7, phi = 1 / 0.7, omega = 0.4,
                                        seed = 31))
  expect_equal(draw_null_count(1, 1.5, 0.3), 0)
  set.seed(32)
  poisson <- vapply(1:200, function(i) draw_null_count(0, 2, 0), numeric(1))
  set.seed(32)
  expect_equal(poisson, qpois(runif(200), 2))
})

test_that("simulate_mixture labels a position non-null when its uniform is at least pi0", {
  # eta = 1: every null count is 0, so a non-null count is offset + signal;
  # each position uses two uniforms (the label, the null count)
  data <- simulate_mixture(50, draw_signal = function() 100, offset = 7, min_kept = 0,
                           pi0 = 0.7, eta = 1, lambda = 1.5, theta = 0, seed = 41)
  set.seed(41)
  label <- as.numeric(runif(100)[seq(1, 99, by = 2)] >= 0.7)
  expect_named(data, c("count", "non_null"))
  expect_equal(data$non_null, label)
  expect_equal(data$count, 107 * label)
  expect_true(any(label == 1))
})

test_that("simulate_mixture draws a non-null count below min_kept again", {
  kept <- simulate_mixture(30, function() 3, offset = 0, min_kept = 3, pi0 = 0.5,
                           eta = 1, lambda = 1, theta = 0, seed = 42)
  expect_equal(sort(unique(kept$count)), c(0, 3))
  dropped <- simulate_mixture(30, function() 3, offset = 0, min_kept = 4, pi0 = 0.5,
                              eta = 1, lambda = 1, theta = 0, seed = 42)
  expect_equal(nrow(dropped), 30)
  expect_true(all(dropped$non_null == 0) && all(dropped$count == 0))
})

test_that("simulate_mixture depends on the seed only", {
  draw <- function(seed) {
    simulate_mixture(40, function() rgeom(1, 0.08), offset = 6, min_kept = 0, pi0 = 0.8,
                     eta = 0.4, lambda = 1.5, theta = 0.3, seed = seed)
  }
  first <- draw(5)
  runif(2)
  expect_identical(draw(5), first)
  expect_false(identical(draw(6), first))
})

test_that("simulate_geometric_mixture adds a geometric signal to the non-null counts", {
  data <- simulate_geometric_mixture(p = 0.08, n = 60, offset = 6, min_kept = 0, pi0 = 0.8,
                                     eta = 0.4, lambda = 1.5, theta = 0.3, seed = 43)
  same <- simulate_mixture(60, function() rgeom(1, prob = 0.08), offset = 6, min_kept = 0,
                           pi0 = 0.8, eta = 0.4, lambda = 1.5, theta = 0.3, seed = 43)
  expect_identical(data, same)
  expect_true(all(data$count[data$non_null == 1] >= 6))
  # p = 1: the signal is 0, a non-null count is the offset plus a null count (0 here)
  flat <- simulate_geometric_mixture(p = 1, n = 40, offset = 6, min_kept = 0, pi0 = 0.5,
                                     eta = 1, lambda = 1.5, theta = 0, seed = 44)
  expect_equal(flat$count, 6 * flat$non_null)
  expect_error(simulate_geometric_mixture(p = 0.08, n = 5, offset = 6, min_kept = 0, pi0 = 0.8,
                                          eta = 0.4, lambda = 1.5, theta = 0))   # no seed
})

test_that("simulate_binomial_mixture adds a binomial signal to the non-null counts", {
  data <- simulate_binomial_mixture(p = 0.2, n = 60, offset = -34, min_kept = 6, pi0 = 0.35,
                                    eta = 0.4, lambda = 1.5, theta = 0, seed = 45)
  same <- simulate_mixture(60, function() rbinom(1, size = 250, prob = 0.2), offset = -34,
                           min_kept = 6, pi0 = 0.35, eta = 0.4, lambda = 1.5, theta = 0,
                           seed = 45)
  expect_identical(data, same)
  expect_true(all(data$count[data$non_null == 1] >= 6))
  # p = 1 and 7 trials: the signal is 7
  flat <- simulate_binomial_mixture(p = 1, size = 7, n = 40, offset = 2, min_kept = 0,
                                    pi0 = 0.5, eta = 1, lambda = 1.5, theta = 0, seed = 46)
  expect_equal(flat$count, 9 * flat$non_null)
})

test_that("simulate_histogram_data simulates the data of a panel", {
  geometric <- HISTOGRAM_SETTINGS[HISTOGRAM_SETTINGS$original == "mixed1", ]   # D = -5
  data <- simulate_histogram_data(geometric, n = 80, seed = 47)
  expect_identical(data, simulate_geometric_mixture(p = 0.08, n = 80, offset = 5 + 1 - 5,
                                                    min_kept = 0, pi0 = 0.80, eta = 0.8,
                                                    lambda = 1.5, theta = 0.30, seed = 47))
  binomial <- HISTOGRAM_SETTINGS[HISTOGRAM_SETTINGS$original == "bar6", ]      # D = -40
  data <- simulate_histogram_data(binomial, n = 80, seed = 48)
  expect_identical(data, simulate_binomial_mixture(p = 0.20, n = 80, offset = 5 + 1 - 40,
                                                   min_kept = 6, pi0 = 0.35, eta = 0.4,
                                                   lambda = 1.5, theta = 0, seed = 48))
  expect_equal(dim(data), c(80, 2))
  expect_equal(mean(data$non_null), 0.65, tolerance = 0.2)
  other <- simulate_histogram_data(binomial, n = 80, seed = 48, cutoff = 9)
  expect_true(mean(other$count) > mean(data$count))   # the non-null counts move up by 4
})

test_that("bar_heights is the frequency of 0, 1, ..., max(x)", {
  expect_equal(bar_heights(c(0, 0, 2, 3, 3, 3)), c("0" = 2, "1" = 0, "2" = 1, "3" = 3))
  expect_equal(bar_heights(c(4, 1)), c("0" = 0, "1" = 1, "2" = 0, "3" = 0, "4" = 1))
  x <- DOMAIN_COUNTS$made1
  expect_equal(sum(bar_heights(x)), length(x))
  expect_equal(unname(bar_heights(x)), count_up_to(x, max(x)))
})

test_that("protein_counts reads the counts of one domain", {
  file <- write_domain_file(new_folder())
  expect_equal(protein_counts("made2", file), DOMAIN_COUNTS$made2)   # padding left out
  expect_equal(protein_counts("made1", file), read_domains(file)$made1)
  expect_error(protein_counts("made1", file.path(new_folder(), "realdata.txt")),
               "not included in this repository")
})

test_that("null_fit_series matches a case worked by hand", {
  x <- c(0, 0, 0, 0, 1, 1, 2, 3, 3, 9)
  fits <- data.frame(cutoff = c("C1", "C1", "C2", "C2"), model = c("ZIGP", "ZIP", "ZIGP", "ZIP"),
                     eta = c(0.4, 0.3, 0.5, 0.5), lambda = c(1.5, 1, 1, 1),
                     theta = c(0.1, 0, 0.2, 0), pi0 = c(0.9, 0.8, 1, 1), C = c(1, 3, 2, 2))
  series <- null_fit_series(x, fits, "C1")
  expect_named(series, c("fEstimation", "C", "values"))
  # fHat up to the larger cutoff (3), zigp up to 1, zip up to 3
  expect_equal(series$fEstimation, rep(c("fHat", "zigp", "zip"), c(4, 2, 4)))
  expect_equal(series$C, c(0:3, 0:1, 0:3))
  expect_equal(series$values[1:4], c(0.4, 0.2, 0.1, 0.2))
  expect_equal(series$values[5:6], 0.9 * c(0.4 + 0.6 * exp(-1.5), 0.6 * 1.5 * exp(-1.6)),
               tolerance = 1e-8)
  expect_equal(series$values[7:10], 0.8 * (0.3 * (0:3 == 0) + 0.7 * dpois(0:3, 1)),
               tolerance = 1e-8)
  expect_equal(nrow(null_fit_series(x, fits, "C2")), 9)
})

test_that("null_fit_series follows the null fits of a domain", {
  x <- DOMAIN_COUNTS$made1
  series <- null_fit_series(x, MADE_FITS, "C2")
  zigp <- MADE_FITS[MADE_FITS$cutoff == "C2" & MADE_FITS$model == "ZIGP", ]
  zip <- MADE_FITS[MADE_FITS$cutoff == "C2" & MADE_FITS$model == "ZIP", ]
  expect_equal(as.vector(table(series$fEstimation)[c("fHat", "zigp", "zip")]),
               c(max(zigp$C, zip$C), zigp$C, zip$C) + 1)
  observed <- series$values[series$fEstimation == "fHat"]
  expect_equal(observed[1], mean(x == 0))
  fitted <- series$values[series$fEstimation == "zip"]
  expect_equal(fitted, zip$pi0 * zigp_density(0:zip$C, zip$eta, zip$lambda, 0))
  expect_equal(fitted[1], observed[1], tolerance = 0.1)   # the null fits the zeros
})
