# -----------------------------------------------------------------------------
# File:    R/cutoff.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Tables 1 to 3, Web Tables S2 to S16.
# Purpose: Choice of the cutoff C below which the counts follow the null, by
#          the two criteria C1 and C2 (truncated null likelihood plus a
#          penalty for the counts above C), with the estimate of the null
#          proportion pi0 and the bound D of the two-stage procedure.
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R"); select_null_fits(x)
# Source:  FINAL R CODES.R (negLogLh; the pi, psi, cut, CO and AD parts of
#          FDRCalculation and RealDataFDR).
# -----------------------------------------------------------------------------

# Cutoff criteria in the order of the tables.
CUTOFFS <- c("C1", "C2")

#' Negative log likelihood of the counts 0, ..., C under the null.
#'
#' @param count Counts of 0, ..., max(x).
#' @param eta,lambda,theta Null parameters.
#' @param C Cutoff.
#' @return The negative log likelihood.
truncated_neg_loglik <- function(count, eta, lambda, theta, C) {
  zero <- count[1] * (-log(eta + (1 - eta) * exp(-lambda)))
  positive <- count[2:(C + 1)] * (-log((1 - eta) * gp_pmf(1:C, lambda, theta)))
  zero + sum(positive)
}

#' Penalties of the two criteria for every candidate cutoff.
#'
#' C1 lets every count above the cutoff have its own probability; C2 uses
#' only the numbers of counts below and above the cutoff. At the last cutoff
#' the C1 penalty uses the frequency of max(x) - 1.
#'
#' @param count Counts of 0, ..., K = max(x).
#' @return List with `C1` and `C2`, each of length K (cutoffs 1, ..., K).
cutoff_penalties <- function(count) {
  K <- length(count) - 1
  M <- sum(count)
  below <- cumsum(count)[-1]
  share <- cumsum(count + SMALL)[-1] / M
  information <- count * log(count + SMALL)
  above <- vapply(seq_len(K), function(C) sum(information[(C + 2):(K + 1)]), numeric(1))
  above[K] <- sum(information[K])
  list(C1 = log(M) * (M - below) - above,
       C2 = below * (-log(share)) + (M - below) * (-log(1 + SMALL - share)))
}

#' Cutoff chosen by a criterion: its first local minimum.
#'
#' @param criterion Values of the criterion for the cutoffs 1, 2, ...
#' @return The chosen cutoff.
choose_cutoff <- function(criterion) {
  min(which(diff(sign(diff(criterion))) > 0) + 2, which.min(criterion))
}

#' Bound D of the two-stage procedure: counts at or above D are rejected.
#'
#' @param model One of `MODELS`.
#' @param C Chosen cutoff.
#' @param lambda,theta Null estimates at the cutoff.
#' @param n Number of positions.
#' @param cap Largest value allowed (max(x) for the protein domain data).
#' @return The bound, a whole number.
rejection_bound <- function(model, C, lambda, theta, n, cap = Inf) {
  null_bound <- if (model %in% c("ZIGP", "GP")) gp_bound(lambda, theta, n) else max(lambda, log(n))
  min(cap, ceiling(max(C + 1, null_bound)))
}

#' Null fit of one model at the cutoff chosen by each criterion.
#'
#' @param model One of `MODELS`.
#' @param estimates 3 x K matrix of EM estimates, from `em_estimates()`.
#' @param count Counts of 0, ..., K.
#' @param penalties Output of `cutoff_penalties()`.
#' @param cap As for `rejection_bound()`.
#' @return Data frame with one row per criterion and the columns `cutoff`,
#'   `model`, `eta`, `lambda`, `theta`, `pi0`, `C`, `D`.
select_model_fits <- function(model, estimates, count, penalties, cap) {
  M <- sum(count)
  cutoffs <- seq_len(ncol(estimates))
  at <- function(f) vapply(cutoffs, function(C) {
    f(estimates[1, C], estimates[2, C], estimates[3, C], C)
  }, numeric(1))
  pi0 <- pmin(1, cumsum(count)[-1] / (M * at(zigp_mass_below)))
  loss <- at(function(eta, lambda, theta, C) truncated_neg_loglik(count, eta, lambda, theta, C))
  C <- vapply(penalties[CUTOFFS], function(penalty) choose_cutoff(penalty + loss), numeric(1))
  D <- vapply(C, function(chosen) {
    rejection_bound(model, chosen, estimates[2, chosen], estimates[3, chosen], M, cap)
  }, numeric(1))
  data.frame(cutoff = CUTOFFS, model = model, eta = estimates[1, C], lambda = estimates[2, C],
             theta = estimates[3, C], pi0 = pi0[C], C = C, D = D, row.names = NULL)
}

#' Null fits of the four models at the cutoffs chosen by the two criteria.
#'
#' @param x Counts.
#' @param cap As for `rejection_bound()`.
#' @return Data frame with one row per criterion and model (criterion first)
#'   and the columns of `select_model_fits()`.
select_null_fits <- function(x, cap = Inf) {
  count <- count_up_to(x, max(x))
  penalties <- cutoff_penalties(count)
  estimates <- em_estimates(x)
  fits <- do.call(rbind, lapply(MODELS, function(model) {
    select_model_fits(model, estimates[[model]], count, penalties, cap)
  }))
  fits[order(match(fits$cutoff, CUTOFFS)), ]
}

#' Replace by NA the parameters a model does not have.
#'
#' @param fits Data frame with the columns `model`, `eta`, `theta`.
#' @return `fits` with `eta` NA for GP and P, `theta` NA for ZIP and P.
mask_absent_parameters <- function(fits) {
  fits$eta[fits$model %in% c("GP", "P")] <- NA
  fits$theta[fits$model %in% c("ZIP", "P")] <- NA
  fits
}
