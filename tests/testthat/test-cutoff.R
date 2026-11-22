# -----------------------------------------------------------------------------
# File:    tests/testthat/test-cutoff.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Unit tests of R/cutoff.R: the truncated likelihood, the two
#          penalties, the chosen cutoff and the bound D on cases worked by
#          hand, and the selected null fits on made-up counts with a known
#          null (300 zero-inflated Poisson draws and five large counts).
# Inputs:  None.
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

X <- DOMAIN_COUNTS$made1

test_that("truncated_neg_loglik matches a case worked by hand", {
  count <- c(5, 3, 2, 7)   # the 7 counts of 3 are above C = 2
  by_hand <- -(5 * log(0.2 + 0.8 * exp(-1)) + 3 * log(0.8 * dpois(1, 1)) +
                 2 * log(0.8 * dpois(2, 1)))
  expect_equal(truncated_neg_loglik(count, eta = 0.2, lambda = 1, theta = 0, C = 2), by_hand)
  expect_equal(by_hand, 11.025039, tolerance = 1e-7)
  # with theta the positive counts use the generalized Poisson probability
  expect_equal(truncated_neg_loglik(count, 0, 1, 0.3, C = 1),
               5 * 1 - 3 * log(gp_pmf(1, 1, 0.3)))
})

test_that("cutoff_penalties matches a case worked by hand", {
  # counts of 0, 1, 2 among M = 10; cutoffs 1 and 2
  penalties <- cutoff_penalties(c(6, 3, 1))
  expect_named(penalties, c("C1", "C2"))
  # C1: log(M) times the number above the cutoff, minus the sum of n log n
  # above the cutoff (at the last cutoff: of the frequency of max(x) - 1)
  expect_equal(penalties$C1, c(log(10) * 1 - 1 * log(1), 0 - 3 * log(3)), tolerance = 1e-8)
  # C2: binomial likelihood of 9 below and 1 above, then of 10 below and 0 above
  expect_equal(penalties$C2, c(-9 * log(0.9) - log(0.1), 0), tolerance = 1e-8)
})

test_that("cutoff_penalties has one value per cutoff", {
  count <- count_up_to(X, max(X))
  penalties <- cutoff_penalties(count)
  expect_equal(lengths(penalties), c(C1 = max(X), C2 = max(X)))
  expect_true(all(is.finite(unlist(penalties))))
  expect_equal(penalties$C1[3], log(305) * sum(X > 3) - sum((count * log(count + SMALL))[5:26]))
})

test_that("choose_cutoff takes the cutoff after the first turn upwards, or the minimum", {
  expect_equal(choose_cutoff(c(5, 3, 4, 1)), 3)      # turns upwards from the second to the third
  expect_equal(choose_cutoff(c(9, 7, 5, 6, 2, 8)), 4)
  expect_equal(choose_cutoff(c(4, 3, 2, 1)), 4)      # no turn: the minimum
  expect_equal(choose_cutoff(c(1, 2, 3)), 1)
  expect_equal(choose_cutoff(c(3, 1, 2)), 2)         # the minimum comes before the turn
})

test_that("rejection_bound is the larger of C + 1 and the bound of the null, rounded up", {
  # ZIP and P: max(lambda, log n) = log(1000) = 6.9
  expect_equal(rejection_bound("ZIP", C = 3, lambda = 1.5, theta = 0, n = 1000), 7)
  expect_equal(rejection_bound("P", C = 3, lambda = 9.2, theta = 0, n = 1000), 10)
  expect_equal(rejection_bound("ZIP", C = 9, lambda = 1.5, theta = 0, n = 1000), 10)
  expect_equal(rejection_bound("ZIP", C = 3, lambda = 1.5, theta = 0, n = 1000, cap = 5), 5)
  # ZIGP and GP: gp_bound = 13.7
  expect_equal(rejection_bound("ZIGP", C = 3, lambda = 1.5, theta = 0.3, n = 1000), 14)
  expect_equal(rejection_bound("GP", C = 3, lambda = 1.5, theta = 0.3, n = 1000),
               ceiling(gp_bound(1.5, 0.3, 1000)))
  expect_equal(rejection_bound("GP", C = 20, lambda = 1.5, theta = 0.3, n = 1000), 21)
})

test_that("select_null_fits has one row per criterion and model, criterion first", {
  expect_equal(CUTOFFS, c("C1", "C2"))
  expect_named(MADE_FITS, c("cutoff", "model", "eta", "lambda", "theta", "pi0", "C", "D"))
  expect_equal(MADE_FITS$cutoff, rep(CUTOFFS, each = 4))
  expect_equal(MADE_FITS$model, rep(MODELS, 2))
  expect_true(all(MADE_FITS$pi0 > 0 & MADE_FITS$pi0 <= 1))
  expect_true(all(MADE_FITS$D >= MADE_FITS$C + 1))
  expect_true(all(MADE_FITS$eta[MADE_FITS$model %in% c("GP", "P")] == 0))
  expect_true(all(MADE_FITS$theta[MADE_FITS$model %in% c("ZIP", "P")] == 0))
})

test_that("select_null_fits bounds D by the cap", {
  capped <- select_null_fits(X, cap = 4)
  expect_equal(capped$D, pmin(4, MADE_FITS$D))
  expect_true(any(MADE_FITS$D > 4))
  expect_equal(capped[names(capped) != "D"], MADE_FITS[names(MADE_FITS) != "D"])
})

test_that("select_null_fits finds the known null of the made-up counts", {
  # null: eta 0.4, lambda 1.5, 300 of 305 counts; the non-null counts start at 15
  for (criterion in CUTOFFS) {
    zip <- MADE_FITS[MADE_FITS$cutoff == criterion & MADE_FITS$model == "ZIP", ]
    expect_equal(zip$eta, 0.4, tolerance = 0.1)
    expect_equal(zip$lambda, 1.5, tolerance = 0.25)
    expect_equal(zip$pi0, 300 / 305, tolerance = 0.03)
    expect_lt(zip$C, 15)
    expect_lt(zip$D, 15)   # every non-null count is rejected by the bound
  }
})

test_that("select_null_fits takes the estimates and pi0 at the chosen cutoff", {
  estimates <- em_estimates(X)
  count <- count_up_to(X, max(X))
  for (i in seq_len(nrow(MADE_FITS))) {
    fit <- MADE_FITS[i, ]
    expect_equal(unname(unlist(fit[c("eta", "lambda", "theta")])),
                 unname(estimates[[fit$model]][, fit$C]))
    expect_equal(fit$pi0, min(1, sum(X <= fit$C) /
                                (305 * zigp_mass_below(fit$eta, fit$lambda, fit$theta, fit$C))))
    expect_equal(fit$D, rejection_bound(fit$model, fit$C, fit$lambda, fit$theta, 305))
  }
  # select_model_fits gives the rows of one model
  zip <- select_model_fits("ZIP", estimates$ZIP, count, cutoff_penalties(count), cap = Inf)
  rows <- MADE_FITS[MADE_FITS$model == "ZIP", ]
  rownames(rows) <- NULL
  expect_equal(zip, rows)
  expect_equal(zip$cutoff, CUTOFFS)
})

test_that("select_model_fits chooses the cutoff from the penalty plus the likelihood", {
  estimates <- em_estimates(X)$P
  count <- count_up_to(X, max(X))
  penalties <- cutoff_penalties(count)
  fits <- select_model_fits("P", estimates, count, penalties, cap = Inf)
  loss <- vapply(seq_len(max(X)), function(C) {
    truncated_neg_loglik(count, estimates[1, C], estimates[2, C], estimates[3, C], C)
  }, numeric(1))
  expect_equal(fits$C, c(choose_cutoff(penalties$C1 + loss), choose_cutoff(penalties$C2 + loss)))
  # a cap bounds D
  capped <- select_model_fits("P", estimates, count, penalties, cap = 3)
  expect_equal(capped$D, pmin(3, fits$D))
  expect_equal(capped[names(capped) != "D"], fits[names(fits) != "D"])
})

test_that("mask_absent_parameters sets the parameters a model does not have to NA", {
  masked <- mask_absent_parameters(MADE_FITS)
  expect_equal(is.na(masked$eta), masked$model %in% c("GP", "P"))
  expect_equal(is.na(masked$theta), masked$model %in% c("ZIP", "P"))
  kept <- c("cutoff", "model", "lambda", "pi0", "C", "D")
  expect_equal(masked[kept], MADE_FITS[kept])
  expect_equal(masked$eta[masked$model == "ZIGP"], MADE_FITS$eta[MADE_FITS$model == "ZIGP"])
})
