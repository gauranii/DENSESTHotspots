# -----------------------------------------------------------------------------
# File:    R/simulation.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Tables 1 and 2, Web Tables S2 to S15.
# Purpose: Simulation study: the 12 null settings, data from the mixture of a
#          ZIGP null and a shifted alternative, one replication (null fits,
#          the three procedures, bias of the estimates), the driver that runs
#          replications in parallel and stores them, and the tables of means
#          and standard deviations over replications.
# Inputs:  archive/simulation/<setting>_C<C>.csv.gz (all replications of a
#          setting; unless rerun).
# Outputs: results/<table>.csv; in a rerun also
#          results/simulation/<setting>_C<C>/runNNNN.csv (one per replication).
# Usage:   source("R/setup.R")
#          simulation_table("Table1", "simulation", "geometric", 0.80, 5, opts, PAPER, PATHS)
# Source:  FINAL R CODES.R (mixture, zigpSimulation, FDRCalculation and its
#          calls in "OVERALL PARAMETER SPECIFICATION" and after).
# -----------------------------------------------------------------------------

# The 12 null settings: alternative f1 (geometric or binomial count added to
# C + 1 + shift), proportion of nulls pi0, and the ZIGP null (eta, lambda, theta).
SIMULATION_SETTINGS <- data.frame(
  setting = c("ZIP1", "ZIGP1", "ZIGP2", "ZIP2", "ZIGP3", "ZIGP4",
              "ZIP3", "ZIGP5", "ZIGP6", "ZIP4", "ZIGP7", "ZIGP8"),
  f1      = rep(c("geometric", "binomial", "geometric", "binomial"), each = 3),
  shift   = rep(c(0, -40, 0, -40), each = 3),
  pi0     = rep(c(0.80, 0.35, 0.85, 0.80), each = 3),
  eta     = rep(c(0.8, 0.4, 0.4, 0.4), each = 3),
  lambda  = c(1.5, 1.5, 3, 1.5, 1, 4, 1.5, 1, 3, 1.5, 1, 3),
  theta   = c(0, 0.3, 0.3, 0, 0.3, 0.3, 0, 0.15, 0.15, 0, 0.2, 0.2)
)

# Alternative distributions: Binomial(size, prob) and Geometric(prob).
F1_BINOMIAL  <- list(size = 250, prob = 0.20)
F1_GEOMETRIC <- list(prob = 0.08)

# Columns of a replication file summarized in the two kinds of tables.
TABLE_COLUMNS <- list(
  simulation = paste0(rep(c("two_stage_", "one_stage_", "storey_"), each = 3),
                      c("R", "FDR", "TPR")),
  bias       = paste0("bias_", c("eta", "lambda", "theta", "pi0", "C"))
)

#' Simulated counts from the mixture of the null and the alternative.
#'
#' A position is null with probability pi0. A non-null count is C + 1 + shift
#' plus a draw of f1 plus a null draw, and is drawn again (with a new label)
#' when it is not above C. Both the binomial and the geometric draw are
#' always made, so the random numbers are used in a fixed order.
#'
#' @param n Number of positions.
#' @param C True cutoff: non-null counts are above C.
#' @param setting Row of `SIMULATION_SETTINGS`.
#' @param seed Seed of the data set.
#' @return List with `value` (counts) and `label` (1 = non-null), length n.
simulate_labelled_counts <- function(n, C, setting, seed) {
  set.seed(seed)
  binomial <- as.numeric(setting$f1 == "binomial")
  draw_null <- function() {
    generate_zigp(1, mu = setting$lambda / (1 - setting$theta), phi = 1 / (1 - setting$theta),
                  omega = setting$eta)
  }
  value <- label <- numeric(n)
  filled <- 0
  while (filled < n) {
    nonnull <- runif(1) >= setting$pi0
    count <- if (nonnull) {
      C + 1 + setting$shift + rbinom(1, F1_BINOMIAL$size, F1_BINOMIAL$prob) * binomial +
        rgeom(1, F1_GEOMETRIC$prob) * (1 - binomial) + draw_null()
    } else {
      draw_null()
    }
    if (nonnull && !(count > C)) next
    filled <- filled + 1
    value[filled] <- count
    label[filled] <- as.numeric(nonnull)
  }
  list(value = value, label = label)
}

#' One replication of a setting.
#'
#' @param setting Row of `SIMULATION_SETTINGS`.
#' @param C True cutoff.
#' @param run Replication number; the data are drawn with seed
#'   `settings$seed + run`.
#' @param settings Settings of the paper (`seed`, `n`, `fdr_level`).
#' @return Data frame with one row per criterion and model: `run`, `cutoff`,
#'   `model`, R, FDR and TPR of the three procedures, the estimates `eta`,
#'   `lambda`, `theta`, `pi0`, `C`, `D` and their bias `bias_*` (estimate
#'   minus truth; NA for a parameter the model does not have).
simulation_replication <- function(setting, C, run, settings) {
  data <- simulate_labelled_counts(settings$n, C, setting, settings$seed + run)
  null <- data$value[data$label == 0]
  nonnull <- data$value[data$label == 1]
  fits <- select_null_fits(c(null, nonnull))
  rates <- t(vapply(seq_len(nrow(fits)), function(i) {
    c(local_fdr_rates(sort(null), sort(nonnull), fits[i, ], settings$fdr_level),
      storey_rates(null, nonnull, fits[i, ], settings$fdr_level))
  }, numeric(9)))
  truth <- c(unlist(setting[c("eta", "lambda", "theta", "pi0")]), C = C)
  bias <- sweep(as.matrix(mask_absent_parameters(fits)[names(truth)]), 2, truth)
  colnames(bias) <- paste0("bias_", names(truth))
  data.frame(run = run, fits[c("cutoff", "model")], rates,
             fits[c("eta", "lambda", "theta", "pi0", "C", "D")], bias, row.names = NULL)
}

# ----- Stored replications -----

#' Write the numbers of a table with 17 significant digits.
#'
#' @param table Data frame.
#' @param file Output file.
#' @return NULL, invisibly.
write_exact_csv <- function(table, file) {
  numeric_columns <- vapply(table, is.double, logical(1))
  table[numeric_columns] <- lapply(table[numeric_columns], sprintf, fmt = "%.17g")
  write.csv(table, file, row.names = FALSE, quote = FALSE)
}

#' Numbers from text with 17 significant digits, exactly.
#'
#' `as.numeric()` can be off by a few units in the last place; the neighbour
#' that prints as the text is returned.
#'
#' @param text Character vector written with `sprintf("%.17g", x)`.
#' @return Numeric vector equal to the numbers written.
parse_exact <- function(text) {
  x <- as.numeric(text)
  off <- which(is.finite(x) & sprintf("%.17g", x) != text)
  half_ulp <- 2^(floor(log2(abs(x[off]))) - 53)
  exact <- rep(NA_real_, length(off))
  for (step in -8:8) {
    candidate <- x[off] + step * half_ulp
    hit <- sprintf("%.17g", candidate) == text[off]
    exact[hit] <- candidate[hit]
  }
  if (anyNA(exact)) stop("could not read exactly: ", text[off][is.na(exact)][1])
  x[off] <- exact
  x
}

#' Read replication files written by `write_exact_csv()`.
#'
#' @param files One file, or several files with the same columns.
#' @return Data frame of all rows, sorted by `run`, numbers read exactly.
read_exact_csv <- function(files) {
  runs <- do.call(rbind, lapply(files, read.csv, colClasses = "character"))
  numeric_columns <- !names(runs) %in% c("cutoff", "model")
  runs[numeric_columns] <- lapply(runs[numeric_columns], parse_exact)
  runs[order(runs$run), ]
}

#' Folder name of a setting and a true cutoff, for example "ZIGP1_C5".
#'
#' @param setting Name of the setting.
#' @param C True cutoff.
#' @return The name.
setting_tag <- function(setting, C) {
  sprintf("%s_C%d", setting, C)
}

#' Run the replications that are not yet stored.
#'
#' Each replication seeds itself, so the result does not depend on the number
#' of parallel processes or on the order of the replications.
#'
#' @param group Rows of `SIMULATION_SETTINGS` to run.
#' @param C True cutoff.
#' @param run_file Function of a setting name and run numbers giving the
#'   files of the replications.
#' @param opts Options: `runs`, `cores`.
#' @param settings Settings of the paper.
#' @return NULL, invisibly.
run_replications <- function(group, C, run_file, opts, settings) {
  jobs <- expand.grid(run = parse_index_set(opts$runs, count = TRUE), row = seq_len(nrow(group)))
  jobs <- jobs[!file.exists(run_file(group$setting[jobs$row], jobs$run)), ]
  log_message("%s, C = %d: %d replications to run", paste(group$setting, collapse = " "), C,
              nrow(jobs))
  status <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
    setting <- group[jobs$row[i], ]
    replication <- simulation_replication(setting, C, jobs$run[i], settings)
    write_exact_csv(replication, run_file(setting$setting, jobs$run[i]))
  }, mc.cores = opts$cores, mc.preschedule = FALSE)
  failed <- vapply(status, inherits, logical(1), what = "try-error")
  if (any(failed)) stop(sum(failed), " replications failed: ", status[failed][[1]])
  invisible(NULL)
}

#' Replications of the settings of a table.
#'
#' All archived replications (archive/simulation/<setting>_C<C>.csv.gz) or,
#' with `opts$rerun`, the requested ones, run here when not yet stored in
#' results/simulation/<setting>_C<C>/.
#'
#' @param group Rows of `SIMULATION_SETTINGS`.
#' @param C True cutoff.
#' @param opts Options: `rerun`, `runs`, `cores`.
#' @param settings Settings of the paper.
#' @param paths Project paths (`archive`, `results`).
#' @return Named list (one element per setting) of data frames with the rows
#'   of `simulation_replication()`, sorted by `run`.
stored_replications <- function(group, C, opts, settings, paths) {
  run_file <- function(setting, run) {
    file.path(paths$results, "simulation", setting_tag(setting, C), sprintf("run%04d.csv", run))
  }
  if (opts$rerun) {
    for (folder in dirname(run_file(group$setting, 1))) dir.create(folder, FALSE, TRUE)
    run_replications(group, C, run_file, opts, settings)
  }
  lapply(setNames(group$setting, group$setting), function(setting) {
    if (opts$rerun) {
      return(read_exact_csv(run_file(setting, parse_index_set(opts$runs, count = TRUE))))
    }
    archived <- file.path(paths$archive, "simulation", paste0(setting_tag(setting, C), ".csv.gz"))
    if (!file.exists(archived)) {
      stop(archived, " is missing; run the replications here with --rerun")
    }
    read_exact_csv(archived)
  })
}

# ----- Tables -----

#' Mean and standard deviation over replications, per criterion and model.
#'
#' Missing values are left out (the FDR of Storey's
#' procedure is NaN in a replication without rejections).
#'
#' @param runs Replications of one setting, sorted by `run`.
#' @param setting Name of the setting.
#' @param columns Columns of `runs` to summarize.
#' @return Data frame with the columns `setting`, `cutoff`, `model`,
#'   `statistic` ("mean" or "se") and `columns`; NA when a value is missing
#'   in every replication.
summarize_runs <- function(runs, setting, columns) {
  statistics <- list(mean = mean, se = sd)
  cells <- expand.grid(statistic = names(statistics), model = MODELS, cutoff = CUTOFFS,
                       stringsAsFactors = FALSE)
  values <- t(vapply(seq_len(nrow(cells)), function(i) {
    rows <- runs[runs$cutoff == cells$cutoff[i] & runs$model == cells$model[i], columns]
    vapply(rows, function(v) {
      if (all(is.na(v))) NA_real_ else statistics[[cells$statistic[i]]](v, na.rm = TRUE)
    }, numeric(1))
  }, numeric(length(columns))))
  data.frame(setting = setting, cells[c("cutoff", "model", "statistic")], values)
}

#' Build a simulation table (R, FDR, TPR) or a bias table of one group of
#' three settings and write it to results/.
#'
#' @param name Table name, for example "Table1" or "WebTableS8".
#' @param kind "simulation" (rates of the three procedures) or "bias" (bias
#'   of the estimates of eta, lambda, theta, pi0 and C).
#' @param f1 Alternative of the group: "geometric" or "binomial".
#' @param pi0 Proportion of nulls of the group.
#' @param C True cutoff: 5, or 1 for the contaminated null.
#' @param opts Options: `rerun`, `runs`, `cores`.
#' @param settings Settings of the paper.
#' @param paths Project paths (`archive`, `results`).
#' @return The table, invisibly.
simulation_table <- function(name, kind, f1, pi0, C, opts, settings, paths) {
  group <- SIMULATION_SETTINGS[SIMULATION_SETTINGS$f1 == f1 & SIMULATION_SETTINGS$pi0 == pi0, ]
  runs <- stored_replications(group, C, opts, settings, paths)
  table <- do.call(rbind, lapply(group$setting, function(setting) {
    summarize_runs(runs[[setting]], setting, TABLE_COLUMNS[[kind]])
  }))
  names(table) <- sub("^bias_", "", names(table))
  file <- sprintf("%s_%s_%s-pi%.2f_C%d.csv", name, kind, f1, pi0, C)
  write.csv(table, output_file(paths$results, file), row.names = FALSE)
  log_message("%s: %d replications per setting", file, length(unique(runs[[1]]$run)))
  print(table, row.names = FALSE)
  invisible(table)
}
