# -----------------------------------------------------------------------------
# File:    tests/testthat/test-em.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Unit tests of R/em.R: the tabulated counts, the starting values
#          and the clamp on cases worked by hand, the E-step sums against
#          their definitions, the EM step at a fixed point (expected counts
#          of a known null) and the EM fit on counts drawn from a known null.
# Inputs:  None.
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

#' Expected counts of 0, ..., C among N draws of a ZIGP null.
expected_counts <- function(N, eta, lambda, theta, C) {
  N * (eta * (0:C == 0) + (1 - eta) * gp_pmf(0:C, lambda, theta))
}

test_that("MODELS and EM_TOLERANCE are those of the paper", {
  expect_equal(MODELS, c("ZIGP", "ZIP", "GP", "P"))
  expect_equal(EM_TOLERANCE, 1e-4)
})

test_that("count_up_to tabulates 0, 1, ..., C", {
  x <- c(0, 0, 1, 3, 3, 5)
  expect_equal(count_up_to(x, 4), c(2, 1, 0, 2, 0))     # the 5 is above C
  expect_equal(count_up_to(x, 5), c(2, 1, 0, 2, 0, 1))
  expect_equal(count_up_to(c(2, 2), 2), c(0, 0, 2))
  expect_equal(sum(count_up_to(x, max(x))), length(x))
})

test_that("em_start gives the moment estimates of the ZIP and Poisson models", {
  x <- c(0, 0, 1, 3)   # half zeros, mean 1, mean of the non-zero counts 2
  expect_equal(em_start(x, "ZIP"), c(0.5, 2, 0))
  expect_equal(em_start(x, "P"), c(0, 2, 0))
})

test_that("em_start follows its definition for the generalized Poisson models", {
  x <- c(0, 0, 0, 1, 1, 2, 2, 3, 5, 9)
  # ZIGP: counts up to the third quartile (2.75): 0 0 0 1 1 2 2
  lower <- c(0, 0, 0, 1, 1, 2, 2)
  phi <- sqrt(var(lower) / mean(lower) - 0.3 / 0.7 * mean(lower) + 1e-10)
  expect_equal(em_start(x, "ZIGP"), c(0.3, mean(lower) / (0.7 * phi), 0))  # phi < 1: theta 0
  expect_lt(phi, 1)
  # GP: non-zero counts up to their median (2): 1 1 2 2
  phi <- sqrt(var(c(1, 1, 2, 2)) / (1.5 + 1e-10))
  expect_equal(em_start(x, "GP"), c(0, mean(x) / phi, 0))
  # overdispersed counts give a positive theta = 1 - 1 / phi
  wide <- c(0, 1, 1, 1, 6, 6, 6, 20)
  start <- em_start(wide, "GP")
  expect_equal(start[3], 1 - 1 / sqrt(var(c(1, 1, 1, 6, 6, 6)) / (3.5 + 1e-10)))
  expect_true(start[3] > 0 && start[3] < 1)
})

test_that("em_clamp keeps the estimates inside the parameter space of the model", {
  outside <- c(-0.1, 0, -0.2)
  expect_equal(em_clamp(outside, "ZIGP"), c(0, SMALL, SMALL))
  expect_equal(em_clamp(outside, "ZIP"), c(0, SMALL, 0))
  expect_equal(em_clamp(outside, "GP"), c(0, SMALL, 0))
  expect_equal(em_clamp(outside, "P"), c(0, SMALL, 0))
  inside <- c(0.3, 2, 0.2)
  expect_equal(em_clamp(inside, "ZIGP"), inside)
  expect_equal(em_clamp(inside, "ZIP"), c(0.3, 2, 0))   # no theta
  expect_equal(em_clamp(inside, "GP"), c(0, 2, 0.2))    # no eta
  expect_equal(em_clamp(inside, "P"), c(0, 2, 0))
})

test_that("em_sums matches its definition for a zero-inflated Poisson null", {
  hold <- c(0.4, 1.5, 0)
  count <- c(50, 20, 10)   # counts of 0, 1, 2: C = 2, n = 80
  sums <- em_sums(hold, count, K = 30, dispersion = FALSE)
  expect_named(sums, c("zero", "gp", "above", "total_above", "total"))
  p0 <- 0.4 + 0.6 * exp(-1.5)
  mass <- 0.4 + 0.6 * ppois(2, 1.5)
  tol <- 1e-8   # SMALL is added to the probabilities
  expect_equal(sums$zero, 50 * 0.4 / p0, tolerance = tol)
  expect_equal(sums$gp, 50 * 0.6 * exp(-1.5) / p0 + 30, tolerance = tol)
  expect_equal(sums$zero + sums$gp, 80, tolerance = tol)
  expect_equal(sums$above, 80 * 0.6 * (1 - ppois(2, 1.5)) / mass, tolerance = tol)
  expect_equal(sums$total_above, 80 * 0.6 * sum(3:30 * dpois(3:30, 1.5)) / mass, tolerance = tol)
  expect_equal(sums$total, 20 + 2 * 10 + sums$total_above, tolerance = tol)
})

test_that("em_sums adds the sums of the lambda and theta updates when asked", {
  hold <- c(0.4, 1.5, 0.2)
  count <- c(50, 20, 10, 4)
  plain <- em_sums(hold, count, K = 40, dispersion = FALSE)
  full <- em_sums(hold, count, K = 40, dispersion = TRUE)
  expect_named(full, c("zero", "gp", "above", "total_above", "total", "lambda", "theta"))
  expect_equal(full[names(plain)], plain)
  # lambda / (lambda + theta v) + theta v / (lambda + theta v) = 1, so the two
  # sums add up to: gp + above + sum of (v - 1) count_v for v >= 1 +
  # expected (sum of counts - number) above C
  mass <- zigp_mass_below(0.4, 1.5, 0.2, 3)
  above <- 84 * 0.6 * sum(gp_pmf(4:40, 1.5, 0.2)) / mass
  expect_equal(full$lambda + full$theta,
               full$gp + full$above + sum(c(0, 0, 1, 2) * count) + full$total_above - above)
  expect_gt(full$theta, 0)
})

test_that("em_step leaves the true ZIP and Poisson parameters where they are", {
  # with the expected counts of the null as data, the truth is a fixed point
  zip <- c(0.4, 1.5, 0)
  expect_equal(em_step(zip, expected_counts(1000, 0.4, 1.5, 0, C = 6), K = 60, "ZIP"), zip,
               tolerance = 1e-7)
  poisson <- c(0, 2, 0)
  expect_equal(em_step(poisson, expected_counts(1000, 0, 2, 0, C = 6), K = 60, "P"), poisson,
               tolerance = 1e-7)
})

test_that("em_step moves towards the truth from elsewhere", {
  count <- expected_counts(1000, 0.4, 1.5, 0, C = 6)
  step <- em_step(c(0.2, 1, 0), count, K = 60, "ZIP")
  expect_true(step[1] > 0.2 && step[1] < 0.4)
  expect_true(step[2] > 1 && step[2] < 1.5)
  expect_equal(step[3], 0)
})

test_that("em_step updates are the ratios of the E-step sums", {
  hold <- c(0.3, 1.2, 0.15)
  count <- c(40, 25, 12, 6, 2)
  s <- em_sums(hold, count, K = 25, dispersion = TRUE)
  eta <- s$zero / (s$zero + s$gp + s$above + SMALL)
  zigp <- em_step(hold, count, K = 25, "ZIGP")
  expect_equal(zigp, c(eta, (s$lambda + SMALL) / (s$gp + s$above + SMALL),
                       (s$theta + SMALL) / (s$total + SMALL)))
  expect_true(all(zigp > 0) && all(zigp[c(1, 3)] < 1))
  gp <- em_step(c(0, 1.2, 0.15), count, K = 25, "GP")
  expect_equal(gp[1], 0)
  expect_true(gp[2] > 0 && gp[3] > 0 && gp[3] < 1)
  expect_equal(em_step(c(0.3, 1.2, 0), count, K = 25, "ZIP")[3], 0)
  expect_equal(em_step(c(0, 1.2, 0), count, K = 25, "P")[c(1, 3)], c(0, 0))
})

test_that("em_fit recovers a known zero-inflated Poisson null", {
  x <- generate_zigp(4000, mu = 1.5, phi = 1, omega = 0.4, seed = 1)
  fit <- em_fit(x, C = 5, "ZIP", em_start(x, "ZIP"))
  expect_equal(fit[1], 0.4, tolerance = 0.1)
  expect_equal(fit[2], 1.5, tolerance = 0.1)
  expect_equal(fit[3], 0)
  # the fit is a fixed point of the EM step up to the tolerance
  count <- count_up_to(x, 5)
  expect_lt(max(abs(em_step(em_clamp(fit, "ZIP"), count, max(x), "ZIP") - fit)), EM_TOLERANCE)
})

test_that("em_fit of the Poisson model is close to the sample mean", {
  set.seed(4)
  x <- rpois(3000, 2)
  fit <- em_fit(x, C = 5, "P", em_start(x, "P"))
  expect_equal(fit[c(1, 3)], c(0, 0))
  expect_equal(fit[2], mean(x), tolerance = 0.01)
  expect_equal(fit[2], 2, tolerance = 0.05)
})

test_that("em_estimates has one column per cutoff and one matrix per model", {
  x <- c(generate_zigp(150, mu = 1.5, phi = 1, omega = 0.4, seed = 6), 9)
  estimates <- em_estimates(x)
  expect_named(estimates, MODELS)
  for (model in MODELS) {
    expect_equal(dim(estimates[[model]]), c(3, max(x)))
    expect_equal(rownames(estimates[[model]]), c("eta", "lambda", "theta"))
    expect_true(all(is.finite(estimates[[model]])))
  }
  expect_equal(unname(estimates$ZIP[, 3]), em_fit(x, 3, "ZIP", em_start(x, "ZIP")))
  expect_equal(unname(estimates$ZIGP[, 4]), em_fit(x, 4, "ZIGP", em_start(x, "ZIGP")))
  expect_true(all(estimates$ZIP["theta", ] == 0))
  expect_true(all(estimates$GP["eta", ] == 0))
  expect_true(all(estimates$P[c("eta", "theta"), ] == 0))
  expect_true(all(estimates$ZIGP["eta", ] >= 0 & estimates$ZIGP["eta", ] <= 1))
})
