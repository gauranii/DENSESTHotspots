# -----------------------------------------------------------------------------
# File:    R/zigp.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. All simulated data.
# Purpose: Zero-inflated generalized Poisson (ZIGP) distribution: random
#          numbers (inversion method).
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R"); generate_zigp(10, mu = 1.5, phi = 1, omega = 0.8)
# Source:  FINAL R CODES.R (generateZigp), the same algorithm as rzigp() of the
#          R package ZIGP 3.8 (V. Erhardt, 2010).
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
#' @return Numeric vector of `n` counts.
generate_zigp <- function(n, mu, phi, omega) {
  if (omega < 0 || omega > 1) stop("omega has to be in [0, 1]")
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
