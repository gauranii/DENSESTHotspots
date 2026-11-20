# -----------------------------------------------------------------------------
# File:    R/fdr.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Tables 1 to 3, Web Tables S6, S7, S12 to S15.
# Purpose: The three multiple testing procedures: one-stage (local FDR),
#          two-stage (local FDR below the bound D, rejection at or above D)
#          and Storey's FDR (modified Benjamini-Hochberg with the estimated
#          pi0), for data with known labels (simulation) and without (real).
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R"); fit <- select_null_fits(x)[1, ]
#          one_stage_count(x, fit, level = 0.10)
# Source:  FINAL R CODES.R (lfdr, f12, f13, f14, ModifiedBH, ModifiedBHFDR,
#          RealMBHFDR; the P, U, W, RV, FD, TD parts of FDRCalculation and the
#          LFDR, U1, U2 parts of RealDataFDR).
# -----------------------------------------------------------------------------

#' Reciprocal of the empirical probability of each count.
#'
#' @param x Counts of a subset of the data.
#' @param all Counts of the whole data set.
#' @return length(all) / (number of elements of `all` equal to x), per element.
inverse_frequency <- function(x, all) {
  length(all) / tabulate(all + 1)[x + 1]
}

#' Local FDR: estimated null probability over empirical probability.
#'
#' @param x Counts.
#' @param weight Their `inverse_frequency()`.
#' @param fit Null fit: list or one-row data frame with `eta`, `lambda`,
#'   `theta`, `pi0`.
#' @return The local FDR of each count.
local_fdr <- function(x, weight, fit) {
  fit$pi0 * zigp_density(x, fit$eta, fit$lambda, fit$theta) * weight
}

#' Number of values below a level, ignoring missing values.
#'
#' @param values Numeric vector.
#' @param level Threshold.
#' @return The count.
count_below <- function(values, level) {
  sum(values[!is.na(values)] < level)
}

#' Number of local FDR values below the level among the counts above C and
#' below D.
#'
#' @param x Counts.
#' @param lfdr Their local FDR.
#' @param fit Null fit with the cutoff `C` and the bound `D`.
#' @param level FDR level.
#' @return The difference of the two counts.
count_between <- function(x, lfdr, fit, level) {
  count_below(lfdr[which(x < fit$D)], level) - count_below(lfdr[which(x <= fit$C)], level)
}

# ----- Data without labels -----

#' Rejections of the one-stage procedure: local FDR below the level.
#'
#' @param x Counts.
#' @param fit Null fit, a row of `select_null_fits()`.
#' @param level FDR level.
#' @return Number of rejected positions.
one_stage_count <- function(x, fit, level) {
  count_below(local_fdr(x, inverse_frequency(x, x), fit), level)
}

#' Rejections of the two-stage procedure.
#'
#' Positions with a count at or above D are rejected; of those below D and
#' above the cutoff C, the ones with a local FDR below the level.
#'
#' @param x Counts.
#' @param fit Null fit, a row of `select_null_fits()`.
#' @param level FDR level.
#' @return Number of rejected positions.
two_stage_count <- function(x, fit, level) {
  lfdr <- local_fdr(x, inverse_frequency(x, x), fit)
  count_between(x, lfdr, fit, level) + length(which(x > fit$D - 1))
}

#' Rejections of Storey's FDR procedure.
#'
#' The frequencies are read from the vector of the frequencies of
#' 0, ..., max(x), at the positions of the rejected counts among the distinct
#' observed counts.
#'
#' @param x Counts.
#' @param fit Null fit, a row of `select_null_fits()`.
#' @param level FDR level.
#' @return Number of rejected positions.
storey_count <- function(x, fit, level) {
  count <- count_up_to(x, max(x))
  observed <- which(count > 0)
  p_values <- storey_p_values(max(x), fit)[observed]
  sum(count[sort(storey_rejected(p_values, count[observed], fit$pi0, level))])
}

# ----- Storey's FDR -----

#' P-values of the counts 0, ..., K under the null: P(X >= x).
#'
#' @param K Largest count.
#' @param fit Null fit with `eta`, `lambda`, `theta`.
#' @return Vector of K + 1 p-values; the p-value of 0 is 1.
storey_p_values <- function(K, fit) {
  below <- cumsum(gp_pmf(0:(K - 1), fit$lambda, fit$theta))
  c(1, pmax(0, (1 - fit$eta) * (1 - below)))
}

#' Distinct counts rejected by the modified Benjamini-Hochberg procedure.
#'
#' The i-th smallest p-value is compared with the level
#' times the share of positions whose count is at least the i-th smallest
#' distinct count, over pi0.
#'
#' @param p_values P-values of the distinct counts, in increasing order of
#'   the counts.
#' @param frequency Frequencies of the distinct counts.
#' @param pi0 Estimated proportion of nulls.
#' @param level FDR level.
#' @return Positions in `p_values` of the rejected counts.
storey_rejected <- function(p_values, frequency, pi0, level) {
  M <- sum(frequency)
  at_least <- c(M, M - cumsum(frequency)[-length(frequency)])
  ordered <- order(p_values)
  passed <- which(p_values[ordered] <= at_least * (level / (M * pi0)))
  if (length(passed) == 0) integer(0) else ordered[seq_len(max(passed))]
}

# ----- Data with labels (simulation) -----

#' Number of rejections, FDR and TPR.
#'
#' @param false,true Numbers of false and true discoveries.
#' @param n_nonnull Number of non-null positions.
#' @param empty FDR reported when nothing is rejected.
#' @return c(R, FDR, TPR).
discovery_rates <- function(false, true, n_nonnull, empty = 0) {
  R <- true + false
  c(R, if (R == 0) empty else false / R, true / n_nonnull)
}

#' One-stage and two-stage procedures on data with labels.
#'
#' The one-stage false discoveries are counted above
#' the cutoff C only, and the local FDR of the non-null counts at or above D
#' uses the empirical probabilities of the smallest non-null counts.
#'
#' @param null,nonnull Sorted counts of the null and non-null positions.
#' @param fit Null fit, a row of `select_null_fits()`.
#' @param level FDR level.
#' @return Named vector: R, FDR, TPR of the two-stage, then of the one-stage
#'   procedure.
local_fdr_rates <- function(null, nonnull, fit, level) {
  all <- c(null, nonnull)
  lfdr_null <- local_fdr(null, inverse_frequency(null, all), fit)
  weight <- inverse_frequency(nonnull, all)
  lower <- which(nonnull < fit$D)
  upper <- nonnull[which(nonnull > fit$D - 1)]
  true_lower <- count_below(local_fdr(nonnull[lower], weight[lower], fit), level)
  true_upper <- count_below(local_fdr(upper, weight[seq_along(upper)], fit), level)

  false_two <- max(0, count_between(null, lfdr_null, fit, level)) +
    length(null) - length(which(null < fit$D))
  true_two <- true_lower + length(nonnull) - length(lower)
  false_one <- count_below(lfdr_null[which(null > fit$C)], level)
  rates <- c(discovery_rates(false_two, true_two, length(nonnull)),
             discovery_rates(false_one, true_lower + true_upper, length(nonnull)))
  setNames(rates, paste0(rep(c("two_stage_", "one_stage_"), each = 3), c("R", "FDR", "TPR")))
}

#' Storey's FDR procedure on data with labels.
#'
#' A single rejected count is not counted, and the FDR
#' is NaN when nothing is rejected.
#'
#' @param null,nonnull Counts of the null and non-null positions.
#' @param fit Null fit, a row of `select_null_fits()`.
#' @param level FDR level.
#' @return Named vector: R, FDR, TPR.
storey_rates <- function(null, nonnull, fit, level) {
  all <- c(null, nonnull)
  distinct <- sort(unique(all))
  frequency <- tabulate(match(all, distinct))
  p_values <- storey_p_values(max(all), fit)[distinct + 1]
  rejected <- storey_rejected(p_values, frequency, fit$pi0, level)
  if (length(rejected) < 2) rejected <- integer(0)
  false <- sum(tabulate(match(null, distinct), length(distinct))[rejected])
  true <- sum(tabulate(match(nonnull, distinct), length(distinct))[rejected])
  setNames(discovery_rates(false, true, length(nonnull), empty = NaN),
           paste0("storey_", c("R", "FDR", "TPR")))
}
