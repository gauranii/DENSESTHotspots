# -----------------------------------------------------------------------------
# File:    R/plots.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. Figures 1 and 2, Web Figures S1 to S6.
# Purpose: Histograms of counts: null and non-null counts of a simulated data
#          set drawn over each other, and the counts of a protein domain.
# Inputs:  None.
# Outputs: figures/<name>.pdf (written by save_histogram()).
# Usage:   source("R/setup.R"); save_simulation_histograms("Fig1", PAPER, PATHS)
# Source:  PLOTS.R of 30 June 2016 and of 24 January 2017, last block of the
#          R history of 24 February 2017, realdata8.Rhistory.
# -----------------------------------------------------------------------------

#' Histogram with the non-null counts (black) under the null counts.
#'
#' @param null,non_null Bar heights of the null and the non-null counts.
#' @return NULL, invisibly.
draw_non_null_first <- function(null, non_null) {
  barplot(non_null, cex.names = 0.8, cex = 0.8, col = "black", space = 0,
          ylim = c(0, max(null) + 20), border = NULL)
  barplot(null, cex.names = 0.8, cex = 0.8, add = TRUE, space = 0, border = "gray")
  invisible(NULL)
}

#' Histogram with the null counts (gray) under the non-null counts (black).
#'
#' @param null,non_null Bar heights of the null and the non-null counts.
#' @return NULL, invisibly.
draw_null_first <- function(null, non_null) {
  largest_count <- max(length(null), length(non_null)) - 1
  barplot(null, cex.names = 0.7, cex = 0.7, col = "gray", space = 0,
          xlim = c(0, largest_count), ylim = c(0, max(null) + 20), border = NULL)
  barplot(non_null, cex.names = 0.7, cex = 0.7, add = TRUE, space = 0, col = "black",
          border = "black")
  invisible(NULL)
}

#' Histogram of the counts of a protein domain.
#'
#' @param heights Bar heights, named by the counts.
#' @return NULL, invisibly.
draw_protein_histogram <- function(heights) {
  barplot(heights, names.arg = names(heights), cex.names = 0.7, cex = 0.7)
  invisible(NULL)
}

#' Save a histogram as pdf (7 x 7 in, the default of `pdf()`).
#'
#' @param name File name.
#' @param figures_dir Folder of the figures.
#' @param draw Function without arguments that draws the histogram.
#' @return NULL, invisibly.
save_histogram <- function(name, figures_dir, draw) {
  pdf(output_file(figures_dir, name))
  on.exit(dev.off())
  draw()
  log_message("wrote figures/%s", name)
}

#' Simulate and save the histograms of one figure.
#'
#' The files are named <figure><panel>_simulation_histogram_<model>-<signal>.pdf.
#'
#' @param figure Figure of the paper (for example "WebFigS3").
#' @param settings Settings of the paper (`n`, `figure_seed`).
#' @param paths Project paths (`figures`).
#' @param panels Table of the panels, as `HISTOGRAM_SETTINGS`.
#' @return NULL, invisibly.
save_simulation_histograms <- function(figure, settings, paths, panels = HISTOGRAM_SETTINGS) {
  panels <- panels[panels$figure == figure, ]
  for (i in seq_len(nrow(panels))) {
    panel <- panels[i, ]
    data <- simulate_histogram_data(panel, settings$n, settings$figure_seed)
    null <- bar_heights(data$count[data$non_null == 0])
    non_null <- bar_heights(data$count[data$non_null == 1])
    draw <- switch(panel$order, non_null = draw_non_null_first, null = draw_null_first)
    name <- sprintf("%s%s_simulation_histogram_%s-%s.pdf", figure, panel$panel, panel$model,
                    panel$signal)
    save_histogram(name, paths$figures, function() draw(null, non_null))
  }
}

#' Save the histograms of the protein domains.
#'
#' The files are named <figure><panel>_protein_histogram_<domain>.pdf.
#'
#' @param figure Figure of the paper ("Fig2").
#' @param domains Names of the domains, in the order of the panels.
#' @param paths Project paths (`data`, `figures`).
#' @return NULL, invisibly.
save_protein_histograms <- function(figure, domains, paths) {
  for (i in seq_along(domains)) {
    heights <- bar_heights(protein_counts(domains[i], paths$data))
    name <- sprintf("%s%s_protein_histogram_%s.pdf", figure, LETTERS[i], domains[i])
    save_histogram(name, paths$figures, function() draw_protein_histogram(heights))
  }
}
