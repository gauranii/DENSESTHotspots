# -----------------------------------------------------------------------------
# File:    tests/testthat/test-zigp.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Unit tests of R/zigp.R: the ZIGP random numbers against the
#          inverse of the distribution function and under a seed, the
#          generalized Poisson and ZIGP probabilities against the Poisson,
#          their sum and values worked by hand, the mass below a cutoff and
#          the bound D.
# Inputs:  None.
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

test_that("generate_zigp with omega = 0 and phi = 1 is the Poisson quantile of the uniforms", {
  set.seed(5)
  draws <- generate_zigp(300, mu = 3, phi = 1, omega = 0)
  set.seed(5)
  expect_equal(draws, qpois(runif(300), 3))
})

test_that("generate_zigp inverts the ZIGP distribution function", {
  # mu = lambda / (1 - theta), phi = 1 / (1 - theta), omega = eta
  eta <- 0.4; lambda <- 1.5; theta <- 0.3
  draws <- generate_zigp(500, mu = lambda / (1 - theta), phi = 1 / (1 - theta), omega = eta,
                         seed = 21)
  set.seed(21)
  u <- runif(500)
  cdf <- cumsum(eta * (0:200 == 0) + (1 - eta) * gp_pmf(0:200, lambda, theta))
  expect_equal(draws, vapply(u, function(v) which(cdf >= v)[1] - 1, numeric(1)))
  expect_true(all(draws == round(draws)) && all(draws >= 0))
})

test_that("generate_zigp sets the seed it is given and leaves the generator alone otherwise", {
  first <- generate_zigp(40, mu = 2, phi = 1.3, omega = 0.3, seed = 9)
  set.seed(9)
  expect_identical(generate_zigp(40, mu = 2, phi = 1.3, omega = 0.3), first)
  expect_identical(generate_zigp(40, mu = 2, phi = 1.3, omega = 0.3, seed = 9), first)
  # without a seed the draws go on from the current state: one uniform per count
  set.seed(9)
  parts <- c(generate_zigp(15, mu = 2, phi = 1.3, omega = 0.3),
             generate_zigp(25, mu = 2, phi = 1.3, omega = 0.3))
  expect_identical(parts, first)
  next_uniform <- runif(1)
  set.seed(9)
  expect_identical(runif(41)[41], next_uniform)
})

test_that("generate_zigp gives zeros when omega = 1 and checks omega", {
  expect_equal(generate_zigp(25, mu = 2, phi = 1.2, omega = 1, seed = 3), rep(0, 25))
  expect_error(generate_zigp(5, mu = 2, phi = 1.2, omega = 1.1), "omega has to be in")
  expect_error(generate_zigp(5, mu = 2, phi = 1.2, omega = -0.1), "omega has to be in")
})

test_that("generate_zigp has the ZIGP share of zeros and mean", {
  draws <- generate_zigp(20000, mu = 2, phi = 1.25, omega = 0.3, seed = 77)
  # lambda = mu / phi = 1.6; P(0) = 0.3 + 0.7 exp(-1.6); mean = 0.7 * 2
  expect_equal(mean(draws == 0), 0.3 + 0.7 * exp(-1.6), tolerance = 0.02)
  expect_equal(mean(draws), 1.4, tolerance = 0.03)
})

test_that("gp_pmf is the Poisson when theta = 0 and sums to one", {
  expect_equal(gp_pmf(0:8, 2.5, 0), dpois(0:8, 2.5))
  p <- gp_pmf(0:400, 2, 0.3)
  expect_true(all(p >= 0))
  expect_equal(sum(p), 1, tolerance = 1e-10)
  expect_equal(sum(0:400 * p), 2 / (1 - 0.3), tolerance = 1e-10)   # mean lambda / (1 - theta)
})

test_that("gp_pmf matches values worked by hand", {
  # P(0) = exp(-lambda), P(1) = lambda exp(-lambda - theta),
  # P(2) = lambda (lambda + 2 theta) exp(-lambda - 2 theta) / 2
  expect_equal(gp_pmf(c(0, 1, 2), 2, 0.5), c(exp(-2), 2 * exp(-2.5), 2 * 3 * exp(-3) / 2))
  expect_equal(gp_pmf(2, 2, 0.5), 0.14936120510359183)
})

test_that("gp_pmf_times_x is x times gp_pmf", {
  x <- 1:30
  expect_equal(gp_pmf_times_x(x, 1.5, 0.2), x * gp_pmf(x, 1.5, 0.2))
  expect_equal(gp_pmf_times_x(3, 2, 0), 3 * dpois(3, 2))
  expect_equal(sum(gp_pmf_times_x(1:400, 1.5, 0.2)), 1.5 / 0.8, tolerance = 1e-10)
})

test_that("gp_density adds SMALL to gp_pmf", {
  expect_equal(SMALL, 1e-10)
  expect_identical(gp_density(0:5, 1.5, 0.1), gp_pmf(0:5, 1.5, 0.1) + SMALL)
  expect_gt(gp_density(500, 1.5, 0), 0)   # positive where the Poisson probability is 0
  expect_equal(gp_pmf(500, 1.5, 0), 0)
})

test_that("zigp_density is eta at zero plus (1 - eta) times the GP part", {
  # worked by hand: P(0) = eta + (1 - eta) exp(-lambda), P(1) = (1 - eta) lambda exp(-lambda - theta)
  expect_equal(zigp_density(c(0, 1), 0.25, 2, 0.5),
               c(0.25 + 0.75 * exp(-2), 0.75 * 2 * exp(-2.5)), tolerance = 1e-9)
  x <- c(3, 0, 7, 0, 1)
  expect_equal(zigp_density(x, 0.4, 1.5, 0), 0.4 * (x == 0) + 0.6 * dpois(x, 1.5),
               tolerance = 1e-9)
  expect_identical(zigp_density(x, 0.4, 1.5, 0.2),
                   0.4 * (x == 0) + 0.6 * gp_density(x, 1.5, 0.2) + SMALL)
  expect_equal(zigp_density(0:6, 0, 2.5, 0), dpois(0:6, 2.5), tolerance = 1e-9)
  expect_equal(sum(zigp_density(0:300, 0.3, 2, 0.2)), 1, tolerance = 1e-7)
})

test_that("zigp_mass_below is the ZIGP distribution function", {
  expect_equal(zigp_mass_below(0.4, 1.5, 0, 3), 0.4 + 0.6 * ppois(3, 1.5))
  expect_equal(zigp_mass_below(0, 2, 0, 0), exp(-2))
  expect_equal(zigp_mass_below(0.25, 2, 0.5, 1), 0.25 + 0.75 * (exp(-2) + 2 * exp(-2.5)))
  expect_equal(zigp_mass_below(0.3, 2, 0.3, 400), 1, tolerance = 1e-10)
  expect_equal(zigp_mass_below(1, 2, 0.3, 2), 1)
})

test_that("gp_bound is the larger of its two terms", {
  # theta = 0.3, n = 1000: lambda / (exp(-0.7) - 0.3) and log(1000) / (-0.7 - log(0.3))
  first <- function(lambda) lambda / (exp(-0.7) - 0.3)
  second <- log(1000) / (-0.7 - log(0.3))
  expect_equal(second, 13.706603, tolerance = 1e-7)
  expect_equal(first(1.5), 7.6302754, tolerance = 1e-6)
  expect_equal(gp_bound(1.5, 0.3, 1000), second)       # the second term is larger
  expect_equal(gp_bound(4, 0.3, 1000), first(4))       # 20.3: the first term is larger
  expect_gt(gp_bound(1.5, 0.3, 10^6), gp_bound(1.5, 0.3, 1000))
})
