# -----------------------------------------------------------------------------
# File:    tests/testthat/test-plots.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Unit tests of R/plots.R: each drawing function draws on a pdf
#          device in a temporary folder, the saving functions write the
#          files of a figure under their names, from simulated data and
#          from a made-up file of domain counts.
# Inputs:  None (a made-up file of domain counts in a temporary folder).
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

#' Draw on a pdf device in a temporary folder and return the file.
drawn <- function(draw) {
  file <- file.path(new_folder(), "drawn.pdf")
  pdf(file)
  tryCatch(draw(), finally = dev.off())
  file
}

NULL_HEIGHTS <- bar_heights(c(0, 0, 0, 1, 1, 2))
NON_NULL_HEIGHTS <- bar_heights(c(6, 7, 7, 9))
PATHS_MADE <- function() {
  folder <- new_folder()
  list(data = write_domain_file(folder), figures = file.path(folder, "figures"))
}

test_that("draw_non_null_first draws the two histograms", {
  file <- drawn(function() expect_null(draw_non_null_first(NULL_HEIGHTS, NON_NULL_HEIGHTS)))
  expect_true(is_pdf(file))
  empty <- drawn(function() plot.new())
  expect_gt(file.size(file), file.size(empty))
})

test_that("draw_null_first draws the two histograms", {
  file <- drawn(function() expect_null(draw_null_first(NULL_HEIGHTS, NON_NULL_HEIGHTS)))
  expect_true(is_pdf(file))
  empty <- drawn(function() plot.new())
  expect_gt(file.size(file), file.size(empty))
})

test_that("draw_protein_histogram draws one bar per count", {
  file <- drawn(function() expect_null(draw_protein_histogram(bar_heights(c(0, 0, 1, 4)))))
  expect_true(is_pdf(file))
})

test_that("draw_null_fit returns the plot of the three series", {
  series <- null_fit_series(DOMAIN_COUNTS$made1, MADE_FITS, "C1")
  plot <- draw_null_fit(series)
  expect_s3_class(plot, "ggplot")
  expect_equal(plot$labels$x, "C")
  expect_equal(plot$labels$y, "estimated probabilities")
  expect_equal(plot$data, series)
  expect_length(plot$layers, 2)   # points and lines
  expect_true(is_pdf(drawn(function() print(plot))))
})

test_that("save_figure writes a pdf file in the folder of the figures", {
  folder <- file.path(new_folder(), "figures")
  printed <- capture.output(save_figure("Fig9A_made.pdf", folder, function() plot(1:3)))
  expect_match(printed, "wrote figures/Fig9A_made.pdf")
  expect_true(is_pdf(file.path(folder, "Fig9A_made.pdf")))
  expect_null(dev.list())   # the device is closed
  # also when the drawing fails
  expect_error(quietly(save_figure("Fig9B_made.pdf", folder, function() stop("no data"))),
               "no data")
  expect_null(dev.list())
})

test_that("save_simulation_histograms writes one file per panel of the figure", {
  paths <- list(figures = file.path(new_folder(), "figures"))
  settings <- list(n = 120, figure_seed = 5)
  quietly(save_simulation_histograms("Fig1", settings, paths))
  expect_equal(sort(list.files(paths$figures)),
               c("Fig1A_simulation_histogram_ZIP1-geometric.pdf",
                 "Fig1B_simulation_histogram_ZIGP1-geometric.pdf",
                 "Fig1C_simulation_histogram_ZIGP2-geometric.pdf",
                 "Fig1D_simulation_histogram_ZIP2-binomial.pdf",
                 "Fig1E_simulation_histogram_ZIGP3-binomial.pdf",
                 "Fig1F_simulation_histogram_ZIGP4-binomial.pdf"))
  expect_true(all(vapply(list.files(paths$figures, full.names = TRUE), is_pdf, logical(1))))
})

test_that("save_simulation_histograms takes the panels it is given", {
  paths <- list(figures = file.path(new_folder(), "figures"))
  panels <- HISTOGRAM_SETTINGS[HISTOGRAM_SETTINGS$original %in% c("mixed1", "mixed4"), ]
  panels$figure <- "FigS99"
  panels$panel <- c("A", "B")
  quietly(save_simulation_histograms("FigS99", list(n = 120, figure_seed = 5), paths, panels))
  expect_equal(sort(list.files(paths$figures)),
               c("FigS99A_simulation_histogram_MZIGP1-geometric.pdf",
                 "FigS99B_simulation_histogram_MZIGP3-binomial.pdf"))
  # a figure without panels writes nothing
  quietly(save_simulation_histograms("Fig7", list(n = 120, figure_seed = 5), paths))
  expect_length(list.files(paths$figures), 2)
})

test_that("save_protein_histograms writes one file per domain, lettered in order", {
  paths <- PATHS_MADE()
  quietly(save_protein_histograms("Fig2", c("made2", "made1"), paths))
  expect_equal(sort(list.files(paths$figures)),
               c("Fig2A_protein_histogram_made2.pdf", "Fig2B_protein_histogram_made1.pdf"))
  expect_true(all(vapply(list.files(paths$figures, full.names = TRUE), is_pdf, logical(1))))
  paths$data <- file.path(new_folder(), "realdata.txt")
  expect_error(save_protein_histograms("Fig2", "made1", paths), "not included")
})

test_that("save_null_fit_plots writes one file per cutoff criterion", {
  paths <- PATHS_MADE()
  value <- quietly(save_null_fit_plots("FigS7", paths, domains = c(FigS7 = "made1")))
  expect_null(value)
  expect_equal(sort(list.files(paths$figures)),
               c("FigS7A_protein_null-fit_made1_C1.pdf", "FigS7B_protein_null-fit_made1_C2.pdf"))
  expect_true(all(vapply(list.files(paths$figures, full.names = TRUE), is_pdf, logical(1))))
})
