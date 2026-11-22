# -----------------------------------------------------------------------------
# File:    R/zigp.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. All simulated data.
# Purpose: Zero-inflated generalized Poisson (ZIGP) distribution: random
#          numbers (inversion method), probability functions, the mass of
#          the null below a cutoff and the bound D of the two-stage procedure.
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R"); generate_zigp(10, mu = 1.5, phi = 1, omega = 0.8)
# Source:  FINAL R CODES.R (generateZigp, the same algorithm as rzigp() of the
#          R package ZIGP 3.8 (V. Erhardt, 2010); gpd, gpdMixture, f0, Dn).
# -----------------------------------------------------------------------------

#' Random numbers from the ZIGP distribution by inversion.
#'
#' The cumulative probabilities are built up to the largest uniform drawn,
#' and each uniform is mapped to the first count whose cumulative probability
#' reaches it. Gives the same numbers as `ZIGP::rzigp()` for the same seed.
#'
#' @param n Number of draws.
#' @param mu Mean parameter of the generalized Poisson part (lambda / (1 - theta)).
#' @param phi Dispersion parameter (1 / (1 - theta)).
#' @param omega Zero-inflation probability (eta), in [0, 1].
#' @param seed Seed set before the uniforms are drawn; NULL (the default)
#'   goes on from the current state of the random number generator.
#' @return Numeric vector of `n` counts.
generate_zigp <- function(n, mu, phi, omega, seed = NULL) {
  if (omega < 0 || omega > 1) stop("omega has to be in [0, 1]")
  if (!is.null(seed)) set.seed(seed)
  x <- double(n)
  u <- runif(n, 0, 1)
  upper <- max(u)
  s <- double(1000)
  p <- omega + (1 - omega) * exp(-mu / phi)  # P(X = 0)
  s[1] <- p
  if (upper > 0) {
    recursive <- FALSE
    i <- 1
    while (s[i] < upper) {
      if (!recursive) {  # P(X = i)
        p <- (1 - omega) * mu * (mu + (phi - 1) * i)^(i - 1) / exp(lgamma(i + 1)) *
          phi^(-i) * exp(-1 / phi * (mu + (phi - 1) * i))
      }
      if (p == Inf) {  # switch to the recursion on the log scale
        recursive <- TRUE
        log_p_previous <- log((1 - omega) * mu * (mu + (phi - 1) * (i - 1))^(i - 2) /
                                exp(lgamma(i - 1 + 1)) * phi^(-(i - 1)) *
                                exp(-1 / phi * (mu + (phi - 1) * (i - 1))))
      }
      if (recursive) {
        log_p <- log((mu + (i - 1) * (phi - 1)) / (phi * i) *
                       (1 + (phi - 1) / (mu + (i - 1) * (phi - 1)))^(i - 1) *
                       exp(1 / phi - 1)) + log_p_previous
        log_p_previous <- log_p
        p <- exp(log_p)
      }
      if (ceiling(i / 1000) == floor(i / 1000)) s <- c(s, double(1000))
      s[i + 1] <- s[i] + p
      i <- i + 1
    }
  }
  for (j in seq_along(u)) {
    i <- 1
    while (u[j] > s[i]) i <- i + 1
    x[j] <- i - 1
  }
  x
}

# Constant added to probabilities and denominators to keep them positive.
SMALL <- 10^(-10)

#' Probability function of the generalized Poisson distribution.
#'
#' @param x Counts.
#' @param lambda,theta Parameters; `theta = 0` gives the Poisson distribution.
#' @return P(X = x), one value per element of `x`.
gp_pmf <- function(x, lambda, theta) {
  exp(log(lambda) + (x - 1) * log(lambda + theta * x) - lfactorial(x) - lambda - theta * x)
}

#' x P(X = x) for the generalized Poisson distribution, on the log scale.
#'
#' @param x Counts (positive).
#' @param lambda,theta Parameters.
#' @return x times `gp_pmf(x, lambda, theta)`.
gp_pmf_times_x <- function(x, lambda, theta) {
  exp(log(x * lambda) + (x - 1) * log(lambda + theta * x) - lfactorial(x) - lambda - theta * x)
}

#' Generalized Poisson probabilities kept positive.
#'
#' @param x Counts.
#' @param lambda,theta Parameters.
#' @return `gp_pmf()` plus `SMALL`.
gp_density <- function(x, lambda, theta) {
  gp_pmf(x, lambda, theta) + SMALL
}

#' Probability function of the ZIGP distribution, kept positive.
#'
#' @param x Counts.
#' @param eta Zero-inflation probability.
#' @param lambda,theta Parameters of the generalized Poisson part.
#' @return eta 1(x = 0) + (1 - eta) `gp_density()`, plus `SMALL`.
zigp_density <- function(x, eta, lambda, theta) {
  eta * (1 * (x == 0)) + (1 - eta) * gp_density(x, lambda, theta) + SMALL
}

#' Probability that a ZIGP count is at most C.
#'
#' @param eta,lambda,theta Parameters of the ZIGP distribution.
#' @param C Cutoff.
#' @return P(X <= C).
zigp_mass_below <- function(eta, lambda, theta, C) {
  eta + (1 - eta) * sum(gp_pmf(0:C, lambda, theta))
}

#' Bound D_n above which a generalized Poisson null gives no count among n.
#'
#' @param lambda,theta Parameters of the generalized Poisson distribution.
#' @param n Number of positions.
#' @return The bound (not rounded).
gp_bound <- function(lambda, theta, n) {
  max(lambda / (exp(theta - 1) - theta), log(n) / (theta - 1 - log(theta)))
}
