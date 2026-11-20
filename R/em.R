# -----------------------------------------------------------------------------
# File:    R/em.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Tables 1 to 3, Web Tables S2 to S16.
# Purpose: EM algorithms for the null distribution fitted to the counts at or
#          below a cutoff C (the counts above C are treated as missing), for
#          four null models: ZIGP, ZIP, generalized Poisson (GP), Poisson (P).
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R"); em_estimates(x)$ZIGP
# Source:  FINAL R CODES.R (zigpEMAlgorithm, zipEMAlgorithm, gpEMAlgorithm,
#          pEMAlgorithm, tau0 to tau3, f, f1 to f11).
# -----------------------------------------------------------------------------

# Null models in the order of the tables.
MODELS <- c("ZIGP", "ZIP", "GP", "P")

# The EM algorithm stops when no parameter moves by more than this.
EM_TOLERANCE <- 10^(-4)

#' Counts of 0, 1, ..., C.
#'
#' @param x Counts.
#' @param C Largest value tabulated.
#' @return Integer vector of length C + 1.
count_up_to <- function(x, C) {
  c(length(x[which(x == 0)]), tabulate(x, nbins = C))
}

#' Starting values of the EM algorithm (moment estimates).
#'
#' @param x Counts.
#' @param model One of `MODELS`.
#' @return c(eta, lambda, theta).
em_start <- function(x, model) {
  zero <- mean(x == 0)
  nonzero <- x[which(x != 0)]
  switch(model,
         ZIGP = {
           lower <- x[which(x <= quantile(x, 0.75))]
           ex <- mean(lower)
           phi <- sqrt(var(lower) / ex - zero / (1 - zero) * ex + SMALL)
           c(max(0, zero), max(0, ex / ((1 - zero) * phi)), min(max(0, 1 - 1 / phi), 1))
         },
         ZIP = c(max(0, zero), max(0, mean(x) / (1 - zero)), 0),
         GP = {
           lower <- nonzero[which(nonzero <= quantile(nonzero, 0.5))]
           phi <- sqrt(var(lower) / (mean(lower) + SMALL))
           c(0, max(0, mean(x) / phi), max(0, 1 - 1 / phi))
         },
         P = c(0, max(0, mean(nonzero)), 0))
}

#' Keep the current estimates inside the parameter space.
#'
#' @param estimate c(eta, lambda, theta).
#' @param model One of `MODELS`.
#' @return c(eta, lambda, theta) used in the next EM step.
em_clamp <- function(estimate, model) {
  eta <- if (model %in% c("ZIGP", "ZIP")) max(0, estimate[1]) else 0
  theta <- switch(model, ZIGP = max(SMALL, estimate[3]), GP = max(0, estimate[3]), 0)
  c(eta, max(SMALL, estimate[2]), theta)
}

#' Expected sums of the E-step for the counts at or below C and above C.
#'
#' The counts above C are replaced by their expectations under the current
#' null, n P(X = v) / P(X <= C) for v = C + 1, ..., K. For C = K the range
#' (C + 1):K covers K + 1 and K.
#'
#' @param hold Current c(eta, lambda, theta).
#' @param count Counts of 0, ..., C.
#' @param K Largest count of the data.
#' @param dispersion If TRUE also the sums of the lambda and theta updates of
#'   the generalized Poisson models.
#' @return List: `zero` (expected structural zeros), `gp` and `above`
#'   (expected numbers from the Poisson part at or below and above C),
#'   `total` (their expected sum of counts), `total_above` (its part above
#'   C), and with `dispersion` the numerators `lambda` and `theta` of the
#'   updates of the generalized Poisson models.
em_sums <- function(hold, count, K, dispersion) {
  eta <- hold[1]
  lambda <- hold[2]
  theta <- hold[3]
  C <- length(count) - 1
  n <- sum(count)
  observed <- 0:C
  missing <- (C + 1):K
  mixture <- zigp_density(observed, eta, lambda, theta)
  from_gp <- (1 - eta) * gp_density(observed, lambda, theta) / mixture
  mass_below <- zigp_mass_below(eta, lambda, theta, C)
  expected <- function(terms) n * ((1 - eta) * sum(terms)) / mass_below

  sums <- list(zero = count[1] * (eta / mixture[1]),
               gp = sum(count * from_gp),
               above = n * ((1 - eta) * (1 - sum(gp_pmf(observed, lambda, theta)))) / mass_below,
               total_above = expected(gp_pmf_times_x(missing, lambda, theta)))
  sums$total <- sum(observed * count * from_gp) + sums$total_above
  if (!dispersion) return(sums)

  pmf <- gp_pmf(missing, lambda, theta)
  pmf_x <- gp_pmf_times_x(missing, lambda, theta)
  lagged <- c(0, seq(0, C - 1, 1)) * count
  to_lambda <- function(v) lambda / (lambda + theta * v)
  to_theta <- function(v) (theta * v) / (lambda + theta * v)
  sums$lambda <- sums$gp + sums$above + sum(lagged * to_lambda(observed)) +
    expected(pmf_x * to_lambda(missing)) - expected(pmf * to_lambda(missing))
  sums$theta <- sum(lagged * to_theta(observed)) +
    expected(pmf_x * to_theta(missing)) - expected(pmf * to_theta(missing))
  sums
}

#' One EM step.
#'
#' @param hold Current c(eta, lambda, theta), inside the parameter space.
#' @param count Counts of 0, ..., C.
#' @param K Largest count of the data.
#' @param model One of `MODELS`.
#' @return Updated c(eta, lambda, theta).
em_step <- function(hold, count, K, model) {
  s <- em_sums(hold, count, K, dispersion = model %in% c("ZIGP", "GP"))
  eta <- s$zero / (s$zero + s$gp + s$above + SMALL)
  switch(model,
         ZIGP = c(eta, (s$lambda + SMALL) / (s$gp + s$above + SMALL),
                  (s$theta + SMALL) / (s$total + SMALL)),
         ZIP = c(eta, s$total / (s$gp + s$above + SMALL), 0),
         GP = c(0, (s$lambda + SMALL) / (s$gp + s$above + SMALL),
                (s$theta + SMALL) / (s$total + SMALL)),
         P = c(0, (sum((0:(length(count) - 1)) * count) + s$total_above) /
                 (sum(count) + s$above + SMALL), 0))
}

#' EM estimates of one null model for one cutoff.
#'
#' @param x Counts.
#' @param C Cutoff: the counts at or below C are used.
#' @param model One of `MODELS`.
#' @param start Starting values, from `em_start()`.
#' @return c(eta, lambda, theta); 0 for a parameter the model does not have.
em_fit <- function(x, C, model, start) {
  count <- count_up_to(x, C)
  K <- max(x)
  estimate <- start
  repeat {
    hold <- em_clamp(estimate, model)
    estimate <- em_step(hold, count, K, model)
    if (!any(abs(hold - estimate) > EM_TOLERANCE)) return(estimate)
  }
}

#' EM estimates of the four null models for every candidate cutoff.
#'
#' @param x Counts.
#' @return Named list (one element per model) of 3 x max(x) matrices with rows
#'   eta, lambda, theta and one column per cutoff C = 1, ..., max(x).
em_estimates <- function(x) {
  lapply(setNames(MODELS, MODELS), function(model) {
    start <- em_start(x, model)
    estimates <- vapply(seq_len(max(x)), function(C) em_fit(x, C, model, start), numeric(3))
    rownames(estimates) <- c("eta", "lambda", "theta")
    estimates
  })
}
