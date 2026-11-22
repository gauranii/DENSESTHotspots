# -----------------------------------------------------------------------------
# File:    python/tests/make_r_values.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Compute with the R code (R/) and with base R the values the Python
#          unit tests compare with, on the small inputs of
#          python/tests/inputs.py, and store them: every double with 17
#          significant digits, which is the double itself.
# Inputs:  None (the inputs are written out below; the same inputs are in
#          python/tests/inputs.py).
# Outputs: python/tests/r_values.json (name -> list of values, as text);
#          python/tests/data/small_domains.txt (counts of two simulated
#          "domains" in the layout of data/realdata.txt).
# Usage:   Rscript python/tests/make_r_values.R   # from the repository root
#          Run again only when R/ changes; the tests read the stored file.
# -----------------------------------------------------------------------------

source("R/setup.R")

values <- list()
put <- function(name, x) {
  if (is.data.frame(x)) x <- unlist(lapply(x[vapply(x, is.numeric, logical(1))], as.double))
  if (!is.character(x)) {
    x <- as.double(x)
    x <- ifelse(is.na(x) & !is.nan(x), "NA", sprintf("%.17g", x))
  }
  stopifnot(!name %in% names(values))
  values[[name]] <<- unname(x)
}

# ----- Inputs (as python/tests/inputs.py) -----
counts <- c(rep(0, 30), rep(1, 14), rep(2, 8), rep(3, 4), rep(4, 2), 5,
            14, 15, 15, 16, 16, 16, 17, 17, 18, 20)
whole <- as.integer(counts)
decimals <- c(0.1, 2.5, 3.7, 1e-3, 12.25, 7.3, 0.6, 5.5, 9.1, 4.4, 8.8, 1.9)
hold <- c(0.4, 1.5, 0.3)
fit <- list(eta = 0.4, lambda = 1.2, theta = 0.1, pi0 = 0.8, C = 5, D = 9)
small <- list(seed = 889191, n = 150, fdr_level = 0.05)
setting <- SIMULATION_SETTINGS[2, ]   # ZIGP1
binomial <- SIMULATION_SETTINGS[5, ]  # ZIGP3

# ----- rng.py -----
set.seed(889192)
put("rng_runif", runif(5))
put("rng_rnorm", rnorm(4, 1, 2))
put("rng_rexp", rexp(4))
put("rng_rbinom", c(rbinom(6, 250, 0.2), rbinom(6, 10, 0.3), rbinom(4, 40, 0.9)))
put("rng_rpois", c(rpois(6, 3.5), rpois(6, 40)))
put("rng_rgeom", rgeom(8, 0.08))

# ----- rmath.py -----
put("rmath_lgamma", lgamma(c(0.5, 1, 3.5, 11, 150, 1e6)))
put("rmath_lfactorial", lfactorial(0:8))
put("rmath_log", log(decimals))
put("rmath_exp", exp(-decimals))
put("rmath_sqrt", sqrt(decimals))
put("rmath_pow", c(10^(-10), 10^(-4), 2.5^3.7, 1.7^2, 0.3^(-6), 7^0, 0^2, 1^55, 3^5))

# ----- rbase.py -----
put("rbase_sum", c(sum(decimals), sum(decimals / 7)))
put("rbase_cumsum", cumsum(decimals / 7))
put("rbase_mean", c(mean(decimals), mean(decimals / 7), mean(counts), mean(whole)))
put("rbase_var", c(var(decimals), var(counts), var(decimals[1:9]), var(c(decimals, NA), na.rm = TRUE),
                   var(counts / 7, na.rm = TRUE)))
put("rbase_sd", c(sd(decimals), sd(c(decimals, NaN), na.rm = TRUE), sd(counts / 7, na.rm = TRUE)))
put("rbase_quantile", c(quantile(decimals, 0.75), quantile(counts, 0.5), quantile(decimals, 0.31)))
numbers <- c(0.1, 1 / 3, 200.913, 1e-5, 123456789.123456789, 1e15, 1e-15, 0.00418141098991527,
             -2.5, 100, 1e5, 1e22, 0.1 + 0.2, 5e-324)
put("rbase_format_number", as.character(numbers))
file <- tempfile()
write.csv(data.frame(name = c("a", "b"), value = c(1 / 3, NA), other = c(NaN, 2)), file,
          row.names = FALSE)
put("rbase_write_csv", readLines(file))

# ----- zigp.py -----
set.seed(21)
put("zigp_generate", generate_zigp(12, mu = 1.5 / 0.7, phi = 1 / 0.7, omega = 0.4))
put("zigp_generate_next", generate_zigp(3, mu = 3, phi = 1, omega = 0))
put("zigp_gp_pmf", c(gp_pmf(0:8, 1.5, 0.3), gp_pmf(0:8, 2, 0)))
put("zigp_gp_pmf_times_x", gp_pmf_times_x(1:8, 1.5, 0.3))
put("zigp_gp_density", gp_density(0:8, 1.5, 0.3))
put("zigp_zigp_density", zigp_density(0:8, 0.4, 1.5, 0.3))
put("zigp_mass_below", c(zigp_mass_below(0.4, 1.5, 0.3, 5), zigp_mass_below(0, 2, 0, 1)))
put("zigp_gp_bound", c(gp_bound(1.5, 0.3, 1000), gp_bound(0.2, 0.6, 150)))
put("zigp_small", SMALL)

# ----- em.py -----
put("em_count_up_to", count_up_to(counts, 7))
for (model in MODELS) {
  put(paste0("em_start_", model), c(em_start(counts, model), em_start(whole, model)))
  put(paste0("em_clamp_", model), c(em_clamp(c(-0.2, -1, -0.5), model), em_clamp(hold, model)))
  put(paste0("em_step_", model), c(em_step(em_clamp(hold, model), count_up_to(counts, 3), 20, model),
                                   em_step(em_clamp(hold, model), count_up_to(counts, 20), 20, model)))
  put(paste0("em_fit_", model), em_fit(counts, 4, model, em_start(counts, model)))
  put(paste0("em_estimates_", model), em_estimates(counts)[[model]])
}
put("em_sums", unlist(em_sums(hold, count_up_to(counts, 3), 20, dispersion = TRUE)))
put("em_sums_last", unlist(em_sums(hold, count_up_to(counts, 20), 20, dispersion = TRUE)))
put("em_sums_plain", unlist(em_sums(c(0.4, 1.5, 0), count_up_to(counts, 3), 20, dispersion = FALSE)))
put("em_tolerance", EM_TOLERANCE)
# A data set on which the EM algorithm of R stops (a parameter becomes NaN).
failing <- simulate_labelled_counts(150, 5, setting, 11)$value
put("em_stops", inherits(try(em_estimates(failing), silent = TRUE), "try-error"))

# ----- cutoff.py -----
count <- count_up_to(counts, max(counts))
put("cutoff_neg_loglik", c(truncated_neg_loglik(count, 0.4, 1.5, 0.3, 4),
                           truncated_neg_loglik(count, 0, 1.5, 0, 20)))
penalties <- cutoff_penalties(count)
put("cutoff_penalties_C1", penalties$C1)
put("cutoff_penalties_C2", penalties$C2)
put("cutoff_choose", c(choose_cutoff(c(5, 4, 3, 4, 1)), choose_cutoff(c(1, 2, 3)),
                       choose_cutoff(c(3, 2, 1)), choose_cutoff(7), choose_cutoff(c(2, 2, 1, 5))))
put("cutoff_rejection_bound", c(rejection_bound("ZIGP", 5, 1.5, 0.3, 1000),
                                rejection_bound("GP", 5, 1.5, 0.3, 1000, cap = 9),
                                rejection_bound("ZIP", 5, 1.5, 0, 1000),
                                rejection_bound("P", 9, 1.5, 0, 1000)))
estimates <- em_estimates(counts)
put("cutoff_model_fits", select_model_fits("ZIGP", estimates$ZIGP, count, penalties, Inf))
fits <- select_null_fits(counts)
put("cutoff_null_fits", fits)
put("cutoff_null_fits_names", paste(fits$cutoff, fits$model))
put("cutoff_null_fits_capped", select_null_fits(whole, cap = 12))
put("cutoff_masked", mask_absent_parameters(fits)[c("eta", "theta")])

# ----- fdr.py -----
part <- c(0, 1, 1, 5, 20)
put("fdr_inverse_frequency", inverse_frequency(part, counts))
put("fdr_local_fdr", local_fdr(part, inverse_frequency(part, counts), fit))
lfdr <- local_fdr(counts, inverse_frequency(counts, counts), fit)
put("fdr_count_below", c(count_below(c(0.01, NA, 0.2, 0.04, NaN), 0.05), count_below(lfdr, 0.05)))
put("fdr_count_between", count_between(counts, lfdr, fit, 0.9))
put("fdr_counts", c(one_stage_count(counts, fit, 0.05), two_stage_count(counts, fit, 0.05),
                    storey_count(counts, fit, 0.05)))
put("fdr_counts_fitted", unlist(lapply(seq_len(nrow(fits)), function(i) {
  c(one_stage_count(counts, fits[i, ], 0.1), two_stage_count(counts, fits[i, ], 0.1),
    storey_count(counts, fits[i, ], 0.1))
})))
put("fdr_storey_p_values", storey_p_values(8, fit))
put("fdr_storey_rejected", storey_rejected(c(1, 0.5, 0.04, 0.001, 0.0002), c(50, 20, 10, 5, 15), 0.8, 0.05))
put("fdr_storey_rejected_none", length(storey_rejected(c(1, 0.9), c(5, 5), 1, 0.05)))
put("fdr_discovery_rates", c(discovery_rates(2, 8, 10), discovery_rates(0, 0, 10),
                             discovery_rates(1, 4, 10, R = 7)))
null <- counts[1:59]
nonnull <- counts[60:69]
put("fdr_local_fdr_rates", local_fdr_rates(sort(null), sort(nonnull), fit, 0.05))
put("fdr_storey_rates", storey_rates(null, nonnull, fit, 0.05))
put("fdr_rates_fitted", unlist(lapply(seq_len(nrow(fits)), function(i) {
  c(local_fdr_rates(sort(null), sort(nonnull), fits[i, ], 0.05),
    storey_rates(null, nonnull, fits[i, ], 0.05))
})))

# ----- simulation.py -----
data <- simulate_labelled_counts(60, 5, setting, 7)
put("simulation_counts", c(data$value, data$label))
data <- simulate_labelled_counts(60, 1, binomial, 8)
put("simulation_counts_binomial", c(data$value, data$label))
replication <- simulation_replication(setting, 5, 3, small)
put("simulation_replication", replication)
write_exact_csv(replication, file)
put("simulation_replication_file", readLines(file))
put("simulation_parse_exact", suppressWarnings(parse_exact(c("0.10000000000000001", "0.30000000000000004", "NA", "12"))))
runs <- do.call(rbind, lapply(1:4, function(run) simulation_replication(binomial, 1, run, small)))
folder <- tempfile()
dir.create(folder)
for (run in 1:4) write_exact_csv(runs[runs$run == run, ], file.path(folder, sprintf("run%04d.csv", run)))
read <- suppressWarnings(read_exact_csv(file.path(folder, sprintf("run%04d.csv", c(3, 1, 4, 2)))))
put("simulation_read_exact", read)
for (kind in names(TABLE_COLUMNS)) {
  table <- summarize_runs(read, "ZIGP3", TABLE_COLUMNS[[kind]])
  put(paste0("simulation_summary_", kind), table)
  write.csv(table, file, row.names = FALSE)
  put(paste0("simulation_summary_file_", kind), readLines(file))
}
put("simulation_setting_tag", setting_tag("ZIGP1", 5))
out <- tempfile()
opts <- list(rerun = TRUE, runs = "1:2", cores = 1)
paths <- list(archive = file.path(REPO_ROOT, "archive"), results = file.path(out, "results"))
invisible(capture.output(table <- suppressWarnings(simulation_table("Table2", "simulation", "binomial", 0.35, 1,
                                                    opts, small, paths))))
put("simulation_table", table)
put("simulation_table_file", readLines(file.path(out, "results",
                                                 "Table2_simulation_binomial-pi0.35_C1.csv")))
put("simulation_run_file", readLines(file.path(out, "results", "simulation", "ZIGP4_C1",
                                               "run0002.csv")))

# ----- figure_data.py -----
set.seed(31)
put("figure_draw_null_count", c(draw_null_count(0.4, 1.5, 0.3), draw_null_count(0.4, 1.5, 0.3),
                                draw_null_count(0, 3, 0)))
mixture <- simulate_mixture(40, function() 2, offset = 6, min_kept = 9, pi0 = 0.5, eta = 0.4,
                            lambda = 1.5, theta = 0.3, seed = 5)
put("figure_mixture", mixture)
put("figure_geometric", simulate_geometric_mixture(0.08, n = 40, offset = 6, min_kept = 0, pi0 = 0.8,
                                                   eta = 0.8, lambda = 1.5, theta = 0.3, seed = 6))
put("figure_binomial", simulate_binomial_mixture(0.2, n = 40, offset = -34, min_kept = 6, pi0 = 0.35,
                                                 eta = 0.4, lambda = 1, theta = 0.3, seed = 6))
for (row in c(1, 4, 13)) {
  put(paste0("figure_histogram_data_", row),
      simulate_histogram_data(HISTOGRAM_SETTINGS[row, ], 50, 889302))
}
put("figure_bar_heights", bar_heights(counts))

# ----- real_data.py, figure_data.py: two simulated "domains" -----
domA <- simulate_labelled_counts(150, 5, setting, 13)$value
domB <- simulate_labelled_counts(120, 5, SIMULATION_SETTINGS[7, ], 12)$value
domains <- data.frame(id = 1:150, domA = as.integer(domA), domB = c(as.integer(domB), rep(NA, 30)))
data_file <- file.path(REPO_ROOT, "python", "tests", "data", "small_domains.txt")
write.table(domains, data_file, sep = "\t", quote = FALSE, row.names = FALSE)
read <- read_domains(data_file)
put("real_domains", c(length(read$domA), length(read$domB), sum(read$domA), sum(read$domB)))
put("real_fit_domain", fit_domain(read$domA))
put("real_rejections", protein_rejections(data_file, 0.10))
put("real_estimates", protein_estimates(data_file))
put("real_protein_counts", protein_counts("domB", data_file))
series <- null_fit_series(read$domA, fit_domain(read$domA), "C1")
put("figure_null_fit_series", series)
put("figure_null_fit_series_names", series$fEstimation)

# ----- Write -----
quote_text <- function(x) paste0('"', gsub('"', '\\\\"', gsub("\\\\", "\\\\\\\\", x)), '"')
lines <- vapply(names(values), function(name) {
  paste0("  ", quote_text(name), ": [", paste(quote_text(values[[name]]), collapse = ", "), "]")
}, character(1))
writeLines(c("{", paste(lines, collapse = ",\n"), "}"),
           file.path(REPO_ROOT, "python", "tests", "r_values.json"))
cat(length(values), "sets of values written\n")
