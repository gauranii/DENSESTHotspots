# -----------------------------------------------------------------------------
# File:    tests/testthat/test-fdr.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Unit tests of R/fdr.R: the local FDR, the rejections of the
#          one-stage, two-stage and Storey procedures and their R, FDR and
#          TPR on 100 counts with a round null fit, worked by hand (the
#          local FDR of each count is listed in helper-inputs.R).
# Inputs:  None.
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

with_fit <- function(...) modifyList(HAND_FIT, list(...))

test_that("inverse_frequency is the number of counts over the frequency of each count", {
  expect_equal(inverse_frequency(c(0, 2), all = c(0, 0, 1, 2)), c(2, 4))
  expect_equal(inverse_frequency(c(1, 1, 0), all = c(0, 0, 1, 2)), c(4, 4, 2))
  expect_equal(inverse_frequency(c(12, 8, 0), HAND_COUNTS), c(50, 50, 100 / 60))
})

test_that("local_fdr is pi0 f0(x) over the empirical probability", {
  weight <- inverse_frequency(HAND_COUNTS, HAND_COUNTS)
  lfdr <- local_fdr(HAND_COUNTS, weight, HAND_FIT)
  expect_length(lfdr, 100)
  expect_equal(lfdr[1], 0.9 * (0.4 + 0.6 * exp(-1.5)) / 0.60, tolerance = 1e-8)   # count 0
  expect_equal(lfdr[61], 0.9 * 0.6 * dpois(1, 1.5) / 0.25, tolerance = 1e-8)      # count 1
  expect_equal(lfdr[95], 0.9 * 0.6 * dpois(6, 1.5) / 0.01, tolerance = 1e-6)      # count 6
  expect_equal(round(lfdr[c(1, 61, 86, 95, 96, 98, 99)], 4),
               c(0.8008, 0.7229, 1.5061, 0.1906, 0.0038, 0.0013, 0))
  expect_equal(local_fdr(3, 2, list(eta = 0, lambda = 2, theta = 0, pi0 = 0.5)), dpois(3, 2),
               tolerance = 1e-8)
})

test_that("count_below counts the values below the level and skips missing ones", {
  expect_equal(count_below(c(0.01, NA, 0.2, 0.04, NaN), 0.05), 2)
  expect_equal(count_below(c(0.05, 0.06), 0.05), 0)   # strictly below
  expect_equal(count_below(numeric(0), 0.05), 0)
})

test_that("count_between counts the small local FDR values above C and below D", {
  x <- c(0, 1, 5, 6, 9)
  lfdr <- c(0.01, 0.5, 0.01, 0.02, 0.001)
  expect_equal(count_between(x, lfdr, list(C = 1, D = 8), 0.05), 2)   # 5 and 6; not 0, not 9
  expect_equal(count_between(x, lfdr, list(C = 5, D = 8), 0.05), 1)
  expect_equal(count_between(x, lfdr, list(C = 1, D = 10), 0.05), 3)
})

test_that("one_stage_count rejects the counts with a local FDR below the level", {
  expect_equal(one_stage_count(HAND_COUNTS, HAND_FIT, 0.05), 5)    # 8, 8, 9, 12, 12
  expect_equal(one_stage_count(HAND_COUNTS, HAND_FIT, 0.20), 6)    # and the 6
  expect_equal(one_stage_count(HAND_COUNTS, HAND_FIT, 0.003), 3)   # 9, 12, 12
  expect_equal(one_stage_count(HAND_COUNTS, HAND_FIT, 0.75), 31)   # and the 25 ones
})

test_that("two_stage_count rejects at or above D, and by the local FDR between C and D", {
  expect_equal(two_stage_count(HAND_COUNTS, HAND_FIT, 0.05), 5)                # 8, 8, 9 + 12, 12
  expect_equal(two_stage_count(HAND_COUNTS, with_fit(D = 6), 0.05), 6)         # all from 6 on
  expect_equal(two_stage_count(HAND_COUNTS, with_fit(D = 9), 0.003), 3)        # 9, 12, 12 by D
  expect_equal(two_stage_count(HAND_COUNTS, with_fit(D = 13), 0.003), 3)       # by the local FDR
  # a count at or below C is never rejected: the ones pass 0.75, C = 2 keeps them
  expect_equal(two_stage_count(HAND_COUNTS, HAND_FIT, 0.75), 6)
  expect_equal(two_stage_count(HAND_COUNTS, with_fit(C = 0), 0.75), 31)
})

test_that("storey_p_values are P(X >= x) under the null", {
  p <- storey_p_values(4, HAND_FIT)
  expect_equal(p, c(1, 0.6 * (1 - ppois(0:3, 1.5))))
  expect_true(all(diff(p) < 0))
  expect_equal(storey_p_values(2, list(eta = 0, lambda = 2, theta = 0)),
               c(1, 1 - exp(-2), 1 - 3 * exp(-2)))
  gp <- storey_p_values(3, list(eta = 0.2, lambda = 1, theta = 0.3))
  expect_equal(gp[3], 0.8 * (1 - gp_pmf(0, 1, 0.3) - gp_pmf(1, 1, 0.3)))
})

test_that("storey_rejected matches a case worked by hand", {
  # four distinct counts with frequencies 70, 20, 6, 4: at least 100, 30, 10, 4
  # thresholds at level 0.05 and pi0 0.8: 0.0625, 0.01875, 0.00625, 0.0025
  # sorted p-values 0.001, 0.01, 0.3, 1: the first two pass
  p <- c(1, 0.3, 0.01, 0.001)
  frequency <- c(70, 20, 6, 4)
  expect_equal(storey_rejected(p, frequency, pi0 = 0.8, level = 0.05), c(4, 3))
  # the smallest p-value belongs to the third count; the next one (0.07) fails
  expect_equal(storey_rejected(c(1, 0.3, 0.018, 0.07), frequency, 0.8, 0.05), 3)
  expect_equal(storey_rejected(c(1, 0.3, 0.07, 0.08), frequency, 0.8, 0.05), integer(0))
  expect_equal(storey_rejected(p, frequency, pi0 = 0.8, level = 1e-4), integer(0))
  # level 1: thresholds 1.25, 0.375, 0.125, 0.05; with pi0 = 0.01 all pass
  expect_equal(storey_rejected(p, frequency, pi0 = 0.8, level = 1), c(4, 3))
  expect_equal(storey_rejected(p, frequency, pi0 = 0.01, level = 1), c(4, 3, 2, 1))
})

test_that("storey_count adds the frequencies of the rejected counts", {
  # every count from 0 to 4 is observed: counts 3 and 4 are rejected
  x <- c(rep(0, 60), rep(1, 25), rep(2, 10), 3, 4, 4, 4, 4)
  fit <- with_fit(lambda = 0.5)
  p <- storey_p_values(4, fit)
  expect_equal(storey_rejected(p, c(60, 25, 10, 1, 4), 0.9, 0.05), c(5, 4))
  expect_equal(storey_count(x, fit, 0.05), 5)
  expect_equal(storey_count(x, fit, 1e-6), 0)
})

test_that("storey_count reads the frequencies at the positions among the distinct counts", {
  # distinct counts 0 1 2 6 8 9 12: the last four are rejected (positions 4 to
  # 7); the frequencies read are those of the counts 3 to 6: 0, 0, 0, 1
  count <- count_up_to(HAND_COUNTS, 12)
  observed <- which(count > 0)
  expect_equal(observed - 1, c(0, 1, 2, 6, 8, 9, 12))
  rejected <- storey_rejected(storey_p_values(12, HAND_FIT)[observed], count[observed], 0.9, 0.05)
  expect_equal(sort(rejected), 4:7)
  expect_equal(storey_count(HAND_COUNTS, HAND_FIT, 0.05), 1)
})

test_that("discovery_rates gives R, FDR and TPR", {
  expect_equal(discovery_rates(false = 2, true = 8, n_nonnull = 20), c(10, 0.2, 0.4))
  expect_equal(discovery_rates(0, 0, 20), c(0, NaN, 0))            # no rejection: no FDR
  expect_equal(discovery_rates(1, 4, 4, R = 6), c(6, 1 / 6, 1))    # R given
})

test_that("local_fdr_rates matches a case worked by hand", {
  # null 6 and 8, non-null 8, 9, 12, 12; local FDR below 0.05 for 8, 9, 12
  rates <- local_fdr_rates(sort(HAND_NULL), sort(HAND_NONNULL), HAND_FIT, 0.05)
  expect_named(rates, c("two_stage_R", "two_stage_FDR", "two_stage_TPR",
                        "one_stage_R", "one_stage_FDR", "one_stage_TPR"))
  # D = 10: the null 8 is the one false discovery of both procedures
  expect_equal(unname(rates), c(5, 0.2, 1, 5, 0.2, 1))
  # D = 6: the null 6 is rejected by the bound too
  rates <- local_fdr_rates(sort(HAND_NULL), sort(HAND_NONNULL), with_fit(D = 6), 0.05)
  expect_equal(unname(rates[1:3]), c(6, 2 / 6, 1))
  expect_equal(unname(rates[4:6]), c(5, 0.2, 1))
  # level 0.003: only 9, 12, 12 pass; D = 13 leaves all to the local FDR
  rates <- local_fdr_rates(sort(HAND_NULL), sort(HAND_NONNULL), with_fit(D = 13), 0.003)
  expect_equal(unname(rates), c(3, 0, 0.75, 3, 0, 0.75))
})

test_that("local_fdr_rates agrees with the counts of the procedures without labels", {
  for (D in c(6, 10, 13)) {
    fit <- with_fit(D = D)
    rates <- local_fdr_rates(sort(HAND_NULL), sort(HAND_NONNULL), fit, 0.05)
    expect_equal(rates[["two_stage_R"]], two_stage_count(HAND_COUNTS, fit, 0.05))
    expect_equal(rates[["one_stage_R"]], one_stage_count(HAND_COUNTS, fit, 0.05))
  }
})

test_that("storey_rates matches a case worked by hand", {
  # rejected: the counts 6, 8, 9, 12; false: the null 6 and 8
  rates <- storey_rates(HAND_NULL, HAND_NONNULL, HAND_FIT, 0.05)
  expect_named(rates, c("storey_R", "storey_FDR", "storey_TPR"))
  expect_equal(unname(rates), c(6, 2 / 6, 1))
  # nothing rejected: R 0, no FDR
  expect_equal(unname(storey_rates(HAND_NULL, HAND_NONNULL, HAND_FIT, 1e-12)), c(0, NaN, 0))
})

test_that("storey_rates does not count a single rejected count", {
  # one distinct count passes (12) at this level
  p <- storey_p_values(12, HAND_FIT)[c(0, 1, 2, 6, 8, 9, 12) + 1]
  level <- 1e-7
  expect_equal(storey_rejected(p, c(60, 25, 9, 1, 2, 1, 2), 0.9, level), 7)
  expect_equal(unname(storey_rates(HAND_NULL, HAND_NONNULL, HAND_FIT, level)), c(0, NaN, 0))
})
