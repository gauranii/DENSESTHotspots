# -----------------------------------------------------------------------------
# File:    tests/testthat/test-simulation.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Unit tests of R/simulation.R: the simulated counts against the
#          order of their random numbers and under a seed, one replication
#          against its parts, the exact reading and writing of numbers, the
#          stored replications (a rerun of two small replications and a
#          made-up archive in a temporary folder) and the tables of means
#          and standard deviations.
# Inputs:  archive/simulation/ZIP1_C5.csv.gz (the first replication).
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

ZIP1 <- SIMULATION_SETTINGS[SIMULATION_SETTINGS$setting == "ZIP1", ]   # geometric signal
ZIP2 <- SIMULATION_SETTINGS[SIMULATION_SETTINGS$setting == "ZIP2", ]   # binomial signal
SMALL_RUN <- list(seed = 500, n = 150, fdr_level = 0.05)
REPLICATION <- simulation_replication(ZIP1, C = 5, run = 2, SMALL_RUN)

#' Write made-up replications as archive/simulation/<setting>_C<C>.csv.gz.
write_archive <- function(folder, setting, C, runs) {
  dir.create(file.path(folder, "simulation"), showWarnings = FALSE)
  file <- file.path(folder, "simulation", paste0(setting_tag(setting, C), ".csv.gz"))
  connection <- gzfile(file, "w")
  write_exact_csv(runs, connection)
  close(connection)
  file
}

test_that("the settings of the simulation study are those of the paper", {
  expect_equal(dim(SIMULATION_SETTINGS), c(12, 7))
  expect_equal(sum(grepl("^ZIP", SIMULATION_SETTINGS$setting)), 4)
  expect_true(all(SIMULATION_SETTINGS$theta[grepl("^ZIP", SIMULATION_SETTINGS$setting)] == 0))
  expect_equal(unique(SIMULATION_SETTINGS[c("f1", "pi0")])$pi0, c(0.80, 0.35, 0.85, 0.80))
  expect_equal(F1_BINOMIAL, list(size = 250, prob = 0.20))
  expect_equal(F1_GEOMETRIC, list(prob = 0.08))
  expect_equal(lengths(TABLE_COLUMNS), c(simulation = 9, bias = 5))
})

test_that("simulate_labelled_counts gives n labelled counts, the non-null ones above C", {
  data <- simulate_labelled_counts(400, C = 5, ZIP1, seed = 3)
  expect_named(data, c("value", "label"))
  expect_equal(lengths(data), c(value = 400, label = 400))
  expect_true(all(data$label %in% c(0, 1)))
  expect_true(all(data$value[data$label == 1] > 5))
  expect_true(all(data$value >= 0 & data$value == round(data$value)))
  expect_equal(mean(data$label == 0), 0.80, tolerance = 0.08)   # pi0
  # null counts: zero-inflated Poisson with eta 0.8
  expect_equal(mean(data$value[data$label == 0] == 0), 0.8 + 0.2 * exp(-1.5), tolerance = 0.06)
})

test_that("simulate_labelled_counts depends on the seed only", {
  first <- simulate_labelled_counts(60, C = 5, ZIP2, seed = 8)
  runif(3)
  expect_identical(simulate_labelled_counts(60, C = 5, ZIP2, seed = 8), first)
  expect_false(identical(simulate_labelled_counts(60, C = 5, ZIP2, seed = 9), first))
})

test_that("simulate_labelled_counts uses its random numbers in a fixed order", {
  # per position: the uniform of the label; for a non-null count the binomial
  # and the geometric draw (both, whatever the signal); the uniform of the
  # null count, inverted through the ZIP distribution function (eta 0.4,
  # lambda 1.5). A non-null count that is not above C is drawn again.
  set.seed(12)
  value <- label <- numeric(0)
  while (length(value) < 40) {
    nonnull <- runif(1) >= 0.35
    signal <- 0
    if (nonnull) {
      signal <- 5 + 1 - 40 + rbinom(1, 250, 0.20)
      rgeom(1, 0.08)
    }
    count <- signal + which(0.4 + 0.6 * ppois(0:100, 1.5) >= runif(1))[1] - 1
    if (nonnull && count <= 5) next
    value <- c(value, count)
    label <- c(label, as.numeric(nonnull))
  }
  expect_equal(simulate_labelled_counts(40, C = 5, ZIP2, seed = 12),
               list(value = value, label = label))
  expect_true(any(label == 1) && any(label == 0))
})

test_that("simulation_replication has one row per criterion and model", {
  expect_equal(dim(REPLICATION), c(8, 23))
  expect_named(REPLICATION, c("run", "cutoff", "model", TABLE_COLUMNS$simulation,
                              "eta", "lambda", "theta", "pi0", "C", "D", TABLE_COLUMNS$bias))
  expect_true(all(REPLICATION$run == 2))
  expect_equal(REPLICATION$cutoff, rep(CUTOFFS, each = 4))
  expect_equal(REPLICATION$model, rep(MODELS, 2))
})

test_that("simulation_replication is the null fits and the rates of its data", {
  data <- simulate_labelled_counts(150, 5, ZIP1, seed = 502)   # seed + run
  null <- data$value[data$label == 0]
  nonnull <- data$value[data$label == 1]
  fits <- select_null_fits(c(null, nonnull))
  estimates <- c("eta", "lambda", "theta", "pi0", "C", "D")
  expect_equal(plain(REPLICATION[estimates]), plain(fits[estimates]))
  for (i in c(2, 5)) {
    expect_equal(unlist(REPLICATION[i, TABLE_COLUMNS$simulation]),
                 c(local_fdr_rates(sort(null), sort(nonnull), fits[i, ], 0.05),
                   storey_rates(null, nonnull, fits[i, ], 0.05)))
  }
  # bias: estimate minus truth (eta 0.8, lambda 1.5, theta 0, pi0 0.80, C 5)
  expect_equal(REPLICATION$bias_lambda, REPLICATION$lambda - 1.5)
  expect_equal(REPLICATION$bias_pi0, REPLICATION$pi0 - 0.80)
  expect_equal(REPLICATION$bias_C, REPLICATION$C - 5)
  zigp <- REPLICATION$model == "ZIGP"
  expect_equal(REPLICATION$bias_eta[zigp], REPLICATION$eta[zigp] - 0.8)
  expect_equal(is.na(REPLICATION$bias_eta), REPLICATION$model %in% c("GP", "P"))
  expect_equal(is.na(REPLICATION$bias_theta), REPLICATION$model %in% c("ZIP", "P"))
})

test_that("simulation_replication gives the first archived replication of ZIP1", {
  archived <- read_exact_csv(file.path(REPO_ROOT, "archive", "simulation", "ZIP1_C5.csv.gz"))
  first <- simulation_replication(ZIP1, C = 5, run = 1, PAPER)
  # the EM algorithm stops at a tolerance of 1e-4
  expect_equal(plain(first), plain(archived[archived$run == 1, ]), tolerance = 1e-4)
})

test_that("parse_exact reads back what sprintf('%.17g') wrote", {
  x <- c(0.1, 1 / 3, pi, 123456789.123456789, -2.5e10, 0, 7, 1e-30 * pi, 1e30 * pi,
         0.039682539682539680, 0.59024390243902436)
  expect_identical(parse_exact(sprintf("%.17g", x)), x)
  set.seed(15)
  y <- c(runif(500), rnorm(500) * 10^sample(-20:20, 500, replace = TRUE))
  expect_identical(parse_exact(sprintf("%.17g", y)), y)
  # a missing value as read.csv gives it, and the values that are not finite
  expect_identical(parse_exact(c(NA, "NaN", "Inf", "2")), c(NA, NaN, Inf, 2))
})

test_that("parse_exact stops at text that is not a number with 17 digits", {
  expect_error(parse_exact("0.1"), "could not read exactly: 0.1")
})

test_that("write_exact_csv and read_exact_csv keep every number", {
  file <- file.path(new_folder(), "run0002.csv")
  write_exact_csv(REPLICATION, file)
  header <- readLines(file, n = 1)
  expect_equal(header, paste(names(REPLICATION), collapse = ","))
  back <- read_exact_csv(file)
  expect_identical(back$cutoff, REPLICATION$cutoff)
  for (column in setdiff(names(REPLICATION), c("cutoff", "model"))) {
    expect_identical(as.numeric(back[[column]]), as.numeric(REPLICATION[[column]]))
  }
  # numbers are written with 17 significant digits
  table <- data.frame(run = 1L, cutoff = "C1", model = "P", value = 1 / 3)
  write_exact_csv(table, file)
  expect_equal(readLines(file), c("run,cutoff,model,value", "1,C1,P,0.33333333333333331"))
})

test_that("read_exact_csv binds several files and sorts by run", {
  folder <- new_folder()
  runs <- made_up_replications(1:3)
  files <- file.path(folder, sprintf("run%04d.csv", 3:1))
  for (i in 1:3) write_exact_csv(runs[runs$run == 4 - i, ], files[i])
  back <- read_exact_csv(files)
  expect_equal(back$run, rep(1:3, each = 8))
  expect_equal(plain(back), plain(runs))
})

test_that("setting_tag names the folder of a setting and a cutoff", {
  expect_equal(setting_tag("ZIGP1", 5), "ZIGP1_C5")
  expect_equal(setting_tag(c("ZIP3", "ZIGP8"), 1), c("ZIP3_C1", "ZIGP8_C1"))
})

test_that("run_replications writes the replications that are not yet stored", {
  folder <- new_folder()
  run_file <- function(setting, run) file.path(folder, sprintf("%s_run%04d.csv", setting, run))
  opts <- list(runs = "2", cores = 1)
  printed <- capture.output(value <- run_replications(ZIP1, 5, run_file, opts, SMALL_RUN))
  expect_null(value)
  expect_match(printed, "ZIP1, C = 5: 2 replications to run")
  expect_equal(sort(list.files(folder)), c("ZIP1_run0001.csv", "ZIP1_run0002.csv"))
  stored <- read_exact_csv(run_file("ZIP1", 2))
  expect_equal(plain(stored), plain(REPLICATION))
  # stored runs are kept: nothing left to run, the files are not written again
  writeLines("kept", run_file("ZIP1", 1))
  printed <- capture.output(run_replications(ZIP1, 5, run_file, opts, SMALL_RUN))
  expect_match(printed, "0 replications to run")
  expect_equal(readLines(run_file("ZIP1", 1)), "kept")
})

test_that("stored_replications reruns into results/simulation/<setting>_C<C>/", {
  folder <- new_folder()
  paths <- list(results = file.path(folder, "results"), archive = file.path(folder, "archive"))
  opts <- list(rerun = TRUE, runs = "2:2", cores = 1)
  runs <- quietly(stored_replications(ZIP1, 5, opts, SMALL_RUN, paths))
  expect_named(runs, "ZIP1")
  expect_true(file.exists(file.path(paths$results, "simulation", "ZIP1_C5", "run0002.csv")))
  expect_equal(plain(runs$ZIP1), plain(REPLICATION))
})

test_that("stored_replications reads the archive unless rerun", {
  folder <- new_folder()
  paths <- list(results = file.path(folder, "results"), archive = folder)
  made <- made_up_replications(c(3, 1, 2))
  write_archive(folder, "ZIP1", 5, made)
  opts <- list(rerun = FALSE, runs = "1000", cores = 1)
  runs <- stored_replications(ZIP1, 5, opts, SMALL_RUN, paths)
  expect_named(runs, "ZIP1")
  expect_equal(runs$ZIP1$run, rep(1:3, each = 8))   # sorted by run
  expect_equal(plain(runs$ZIP1), plain(made[order(made$run), ]))
  expect_false(dir.exists(paths$results))           # nothing is run or written
  expect_error(stored_replications(ZIP2, 5, opts, SMALL_RUN, paths),
               "ZIP2_C5.csv.gz is missing; run the replications here with --rerun")
})

test_that("summarize_runs gives the mean and the standard deviation per criterion and model", {
  runs <- made_up_replications(1:4)
  table <- summarize_runs(runs, "ZIP1", TABLE_COLUMNS$simulation)
  expect_named(table, c("setting", "cutoff", "model", "statistic", TABLE_COLUMNS$simulation))
  expect_equal(nrow(table), 16)
  expect_equal(table$cutoff, rep(CUTOFFS, each = 8))
  expect_equal(table$model, rep(rep(MODELS, each = 2), 2))
  expect_equal(table$statistic, rep(c("mean", "se"), 8))
  expect_true(all(table$setting == "ZIP1"))
  # made-up values: run + 10 model + 100 criterion + column / 7, runs 1 to 4
  expect_equal(table$two_stage_R[1:2], c(2.5 + 10 + 100 + 1 / 7, sd(1:4)))
  row <- table$cutoff == "C2" & table$model == "GP"
  expect_equal(table$storey_TPR[row], c(2.5 + 30 + 200 + 9 / 7, sd(1:4)))
})

test_that("summarize_runs leaves out missing values and writes 0 for an FDR without value", {
  runs <- made_up_replications(1:3)
  first <- runs$cutoff == "C1" & runs$model == "ZIGP"
  runs$two_stage_FDR[first] <- c(0.2, NaN, 0.4)          # one replication without rejections
  runs$one_stage_FDR[first] <- c(NaN, NaN, 0.3)          # a single value: no standard deviation
  runs$storey_FDR[first] <- NaN                          # no replication with rejections
  runs$storey_TPR[first] <- NA
  table <- summarize_runs(runs, "ZIP1", TABLE_COLUMNS$simulation)[1:2, ]
  expect_equal(table$two_stage_FDR, c(0.3, sd(c(0.2, 0.4))))
  expect_equal(table$one_stage_FDR, c(0.3, 0))
  expect_equal(table$storey_FDR, c(0, 0))
  expect_equal(table$storey_TPR, c(NA_real_, NA_real_))
})

test_that("simulation_table writes the table of the three settings of a group", {
  folder <- new_folder()
  paths <- list(results = file.path(folder, "results"), archive = folder)
  opts <- list(rerun = FALSE, runs = "1000", cores = 1)
  group <- c("ZIP1", "ZIGP1", "ZIGP2")   # geometric signal, pi0 = 0.80
  for (i in 1:3) {
    made <- made_up_replications(1:4)
    made[-(1:3)] <- made[-(1:3)] + 1000 * i
    write_archive(folder, group[i], 5, made)
  }
  printed <- capture.output(
    table <- simulation_table("Table1", "simulation", "geometric", 0.80, 5, opts, SMALL_RUN, paths)
  )
  expect_true(any(grepl("Table1_simulation_geometric-pi0.80_C5.csv: 4 replications per setting",
                        printed)))
  file <- file.path(paths$results, "Table1_simulation_geometric-pi0.80_C5.csv")
  expect_true(file.exists(file))
  expect_equal(dim(table), c(48, 13))
  expect_equal(table$setting, rep(group, each = 16))
  expect_equal(table$two_stage_R[c(1, 2, 17, 33)],
               c(1112.5 + 1 / 7, sd(1:4), 2112.5 + 1 / 7, 3112.5 + 1 / 7))
  written <- read.csv(file)
  expect_equal(plain(written), plain(table))
})

test_that("simulation_table of the bias names its columns by the parameters", {
  folder <- new_folder()
  paths <- list(results = file.path(folder, "results"), archive = folder)
  opts <- list(rerun = FALSE, runs = "1000", cores = 1)
  for (setting in c("ZIP4", "ZIGP7", "ZIGP8")) {   # binomial signal, pi0 = 0.80
    write_archive(folder, setting, 1, made_up_replications(1:2))
  }
  table <- quietly(simulation_table("WebTableS9", "bias", "binomial", 0.80, 1, opts, SMALL_RUN,
                                    paths))
  expect_named(table, c("setting", "cutoff", "model", "statistic",
                        "eta", "lambda", "theta", "pi0", "C"))
  expect_true(file.exists(file.path(paths$results, "WebTableS9_bias_binomial-pi0.80_C1.csv")))
  # bias_eta is column 16 of the made-up values
  expect_equal(table$eta[1:2], c(1.5 + 10 + 100 + 16 / 7, sd(1:2)))
})
