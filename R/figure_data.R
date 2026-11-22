# -----------------------------------------------------------------------------
# File:    R/figure_data.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. Figures 1 and 2, Web Figures S1 to S12.
# Purpose: Data of the figures: the settings of the 24 simulated data sets,
#          the generators of a null / non-null mixture of counts, the bar
#          heights of a histogram of counts, and the observed and fitted
#          null probabilities of a protein domain (Web Figures S7 to S12).
# Inputs:  data/realdata.txt (protein domain counts, read by protein_counts()).
# Outputs: None.
# Usage:   source("R/setup.R")
#          data <- simulate_histogram_data(HISTOGRAM_SETTINGS[1, ], 1000, 889302)
# Source:  PLOTS.R of 30 June 2016 (mixture1; data<k>.pdf), PLOTS.R of
#          24 January 2017 (mixture2; bar<k>.pdf), last block of the R history
#          of 24 February 2017 (mixture3, mixture4; mixed<k>.pdf),
#          realdata8.Rhistory (protein histograms); FINALREAL 0418.R of
#          18 April 2017 (null-fit plots; the code of the published version
#          of July 2017 is not kept, see null_fit_series()).
# -----------------------------------------------------------------------------

CUTOFF        <- 5    # C: non-null counts are built as C + 1 + D + signal + null count
BINOMIAL_SIZE <- 250  # number of trials of the binomial signal

# Settings of the simulated data sets, one row per panel.
#   original  file name of the panel in the original code
#   signal    distribution of the signal of a non-null count, with parameter p
#   D         shift of the non-null counts
#   min_kept  a non-null count below this value is drawn again
#   order     the counts drawn first; the other ones are drawn on top
#   pi0       proportion of null counts
#   eta, lambda, theta  zero inflation, rate and dispersion of the null counts
HISTOGRAM_SETTINGS <- read.table(header = TRUE, stringsAsFactors = FALSE, text = "
  figure    panel model  original signal    D   min_kept order    pi0  p    eta lambda theta
  Fig1      A     ZIP1   data3    geometric 0   0        non_null 0.80 0.08 0.8 1.5    0
  Fig1      B     ZIGP1  data1    geometric 0   0        non_null 0.80 0.08 0.8 1.5    0.30
  Fig1      C     ZIGP2  data2    geometric 0   0        non_null 0.80 0.08 0.8 3      0.30
  Fig1      D     ZIP2   bar6     binomial  -40 6        non_null 0.35 0.20 0.4 1.5    0
  Fig1      E     ZIGP3  bar4     binomial  -40 6        non_null 0.35 0.20 0.4 1      0.30
  Fig1      F     ZIGP4  bar5     binomial  -40 6        non_null 0.35 0.20 0.4 4      0.30
  WebFigS1  A     ZIGP5  data5    geometric 0   0        non_null 0.85 0.08 0.4 1      0.15
  WebFigS1  B     ZIGP6  data6    geometric 0   0        non_null 0.85 0.08 0.4 3      0.15
  WebFigS1  C     ZIP3   data7    geometric 0   0        non_null 0.85 0.08 0.4 1.5    0
  WebFigS2  A     ZIGP7  bar1     binomial  -40 6        non_null 0.80 0.20 0.4 1      0.20
  WebFigS2  B     ZIGP8  bar2     binomial  -40 6        non_null 0.80 0.20 0.4 3      0.20
  WebFigS2  C     ZIP4   bar3     binomial  -40 6        non_null 0.80 0.20 0.4 1.5    0
  WebFigS3  A     MZIGP1 mixed1   geometric -5  0        null     0.80 0.08 0.8 1.5    0.30
  WebFigS3  B     MZIGP2 mixed2   geometric -5  0        null     0.80 0.08 0.8 3      0.30
  WebFigS3  C     MZIP1  mixed3   geometric -5  0        null     0.80 0.08 0.8 1.5    0
  WebFigS4  A     MZIGP3 mixed4   binomial  -40 0        null     0.35 0.20 0.4 1      0.30
  WebFigS4  B     MZIGP4 mixed5   binomial  -40 0        null     0.35 0.20 0.4 4      0.30
  WebFigS4  C     MZIP2  mixed6   binomial  -40 0        null     0.35 0.20 0.4 1.5    0
  WebFigS5  A     MZIGP5 mixed7   geometric -5  0        null     0.85 0.08 0.4 1      0.15
  WebFigS5  B     MZIGP6 mixed8   geometric -5  0        null     0.85 0.08 0.4 3      0.15
  WebFigS5  C     MZIP3  mixed9   geometric -5  0        null     0.85 0.08 0.4 1.5    0
  WebFigS6  A     MZIGP7 mixed10  binomial  -40 0        null     0.80 0.20 0.4 1      0.20
  WebFigS6  B     MZIGP8 mixed11  binomial  -40 0        null     0.80 0.20 0.4 3      0.20
  WebFigS6  C     MZIP4  mixed12  binomial  -40 0        null     0.80 0.20 0.4 1.5    0
")

# Protein domains of Figure 2, in the order of the panels.
PROTEIN_DOMAINS <- c("cd00031", "cd00180", "cd00204", "cd00882", "cd00054", "pfam00001")

#' One null count: a ZIGP draw in the parametrization of the paper.
#'
#' @param eta Zero-inflation probability.
#' @param lambda Rate.
#' @param theta Dispersion (0 gives the zero-inflated Poisson).
#' @return One count.
draw_null_count <- function(eta, lambda, theta) {
  generate_zigp(1, mu = lambda / (1 - theta), phi = 1 / (1 - theta), omega = eta)
}

#' Simulate a mixture of null and non-null counts.
#'
#' Each position is non-null when a uniform draw is at least `pi0`. A null
#' count is a ZIGP draw. A non-null count is `offset` plus a signal plus a
#' ZIGP draw (drawn in this order); it is dropped, and the position drawn
#' again, when it is below `min_kept`.
#'
#' @param n Number of counts.
#' @param draw_signal Function without arguments that draws one signal.
#' @param offset Constant added to every non-null count.
#' @param min_kept Smallest non-null count that is kept.
#' @param pi0 Proportion of null counts.
#' @param eta,lambda,theta Parameters of the null counts.
#' @param seed Seed of the random numbers.
#' @return Data frame with `n` rows: `count` and `non_null` (0 or 1).
simulate_mixture <- function(n, draw_signal, offset, min_kept, pi0, eta, lambda, theta, seed) {
  set.seed(seed)
  count <- numeric(0)
  non_null <- numeric(0)
  while (length(count) < n) {
    is_non_null <- runif(1) >= pi0
    signal <- if (is_non_null) offset + draw_signal() else 0
    value <- signal + draw_null_count(eta, lambda, theta)
    if (!is_non_null || value >= min_kept) {
      count <- c(count, value)
      non_null <- c(non_null, as.numeric(is_non_null))
    }
  }
  data.frame(count = count, non_null = non_null)
}

#' Mixture whose non-null counts carry a geometric signal.
#'
#' @param p Success probability of the geometric signal.
#' @param ... Other arguments of `simulate_mixture()`.
#' @param seed Seed of the random numbers, set by `simulate_mixture()` before
#'   the first draw.
#' @return Data frame with columns `count` and `non_null`.
simulate_geometric_mixture <- function(p, ..., seed) {
  simulate_mixture(draw_signal = function() rgeom(1, prob = p), ..., seed = seed)
}

#' Mixture whose non-null counts carry a binomial signal.
#'
#' @param p Success probability of the binomial signal.
#' @param size Number of trials of the binomial signal.
#' @param ... Other arguments of `simulate_mixture()`.
#' @param seed Seed of the random numbers, set by `simulate_mixture()` before
#'   the first draw.
#' @return Data frame with columns `count` and `non_null`.
simulate_binomial_mixture <- function(p, size = BINOMIAL_SIZE, ..., seed) {
  simulate_mixture(draw_signal = function() rbinom(1, size = size, prob = p), ..., seed = seed)
}

#' Simulate the data set of one panel.
#'
#' @param setting One row of `HISTOGRAM_SETTINGS`.
#' @param n Number of counts.
#' @param seed Seed of the random numbers.
#' @param cutoff The constant C of the non-null counts.
#' @return Data frame with columns `count` and `non_null`.
simulate_histogram_data <- function(setting, n, seed, cutoff = CUTOFF) {
  simulate <- switch(setting$signal, geometric = simulate_geometric_mixture,
                     binomial = simulate_binomial_mixture)
  simulate(p = setting$p, n = n, offset = cutoff + 1 + setting$D, min_kept = setting$min_kept,
           pi0 = setting$pi0, eta = setting$eta, lambda = setting$lambda,
           theta = setting$theta, seed = seed)
}

#' Bar heights of a histogram of counts: the frequency of 0, 1, ..., max(x).
#'
#' @param x Vector of non-negative counts.
#' @return Frequencies, named by the counts.
bar_heights <- function(x) {
  heights <- c(sum(x == 0), tabulate(x, nbins = max(x)))
  names(heights) <- 0:max(x)
  heights
}

#' Counts of one protein domain.
#'
#' @param domain Name of the domain (a column of the data file).
#' @param data_file Path of data/realdata.txt.
#' @return Vector of the counts at the positions of the domain.
protein_counts <- function(domain, data_file) {
  counts <- read.delim(require_data(data_file))[[domain]]
  counts[!is.na(counts)]
}

# Protein domains of Web Figures S7 to S12, one figure per domain; panel A
# uses the cutoff criterion C1, panel B the criterion C2.
NULL_FIT_DOMAINS <- c(WebFigS7 = "cd00031", WebFigS8 = "cd00054", WebFigS9 = "cd00180",
                      WebFigS10 = "cd00204", WebFigS11 = "cd00882", WebFigS12 = "pfam00001")

#' Observed and fitted null probabilities of a protein domain.
#'
#' For one cutoff criterion: the observed share of positions with count
#' 0, 1, ..., K ("fHat"), and pi0 times the fitted ZIGP and ZIP probability
#' of each count up to the cutoff chosen for that model ("zigp", "zip"),
#' where K is the larger of the two cutoffs. The April 2017 draft of the
#' plot (FINALREAL 0418.R) drew fHat only up to the ZIP cutoff; the
#' published panels draw it up to K, as here.
#'
#' @param x Counts of the domain.
#' @param fits Null fits of the domain, from `fit_domain()`.
#' @param criterion One of `CUTOFFS`.
#' @return Data frame with the columns `fEstimation` (series), `C` (count)
#'   and `values` (probability).
null_fit_series <- function(x, fits, criterion) {
  fit <- function(model) fits[fits$cutoff == criterion & fits$model == model, ]
  zigp <- fit("ZIGP")
  zip <- fit("ZIP")
  K <- max(zigp$C, zip$C)
  fitted <- function(f) f$pi0 * zigp_density(0:f$C, f$eta, f$lambda, f$theta)
  series <- list(fHat = count_up_to(x, K) / length(x), zigp = fitted(zigp), zip = fitted(zip))
  data.frame(fEstimation = rep(names(series), lengths(series)),
             C = unlist(lapply(series, function(v) seq_along(v) - 1)),
             values = unlist(series), row.names = NULL)
}
