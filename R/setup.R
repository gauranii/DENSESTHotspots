# -----------------------------------------------------------------------------
# File:    R/setup.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. All figures and tables.
# Purpose: Loader sourced first by every script: finds the repository root,
#          sources the library files in R/, and defines the settings of the
#          paper, the project paths and the command line parser.
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R")   # from a script started in the repository root
# Source:  FINAL R CODES.R (2017), section "OVERALL PARAMETER SPECIFICATION".
# -----------------------------------------------------------------------------

REPO_ROOT <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))

for (file in c("zigp.R", "em.R", "cutoff.R", "fdr.R", "simulation.R",
               "real_data.R", "figure_data.R", "plots.R")) {
  source(file.path(REPO_ROOT, "R", file))
}

# ----- Settings of the paper, passed to the functions as `settings` -----
PAPER <- list(
  seed          = 889191,        # replication k of a simulation uses seed + k
  figure_seed   = 889191 + 111,  # seed of the simulated data in the histograms
  n             = 1000,          # positions per simulated data set
  replications  = 1000,          # replications per simulation setting
  fdr_level     = 0.05,          # FDR level of the simulation study
  fdr_level_real = 0.10          # FDR level of the protein domain data
)

# ----- Paths -----

#' Folders and files of the project.
#'
#' @param root Repository root.
#' @param output Folder in which results/ and figures/ are written.
#' @return Named list of paths: `data` (ZINULL_DATA when set), `archive`,
#'   `results`, `figures`.
project_paths <- function(root, output = Sys.getenv("ZINULL_OUTPUT", root)) {
  list(data    = Sys.getenv("ZINULL_DATA", file.path(root, "data", "realdata.txt")),
       archive = file.path(root, "archive"),
       results = file.path(output, "results"),
       figures = file.path(output, "figures"))
}

PATHS <- project_paths(REPO_ROOT)

#' Stop with a message when the protein domain data file is not there.
#'
#' @param path Path of data/realdata.txt.
#' @return `path`, invisibly.
require_data <- function(path) {
  if (!file.exists(path))
    stop("The protein domain counts (data/realdata.txt) are not included in ",
         "this repository. Place the file in data/ of the directory the ",
         "command runs from, or set ZINULL_DATA to its path.",
         call. = FALSE)
  invisible(path)
}

#' Path of an output file, creating its folder when needed.
#'
#' @param folder Output folder (for example `paths$results`).
#' @param name File name.
#' @return The path `folder/name`.
output_file <- function(folder, name) {
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  file.path(folder, name)
}

# ----- Command line -----

#' Read `--name=value` arguments and bare `--name` flags.
#'
#' Every script accepts `--rerun`: recompute from scratch instead of reading
#' the archived runs.
#'
#' @param defaults Named list of options and their default values. A value
#'   given on the command line is converted to the type of its default; a
#'   flag sets a logical default to TRUE.
#' @param args Command line arguments.
#' @return `defaults` with the values given on the command line, and `rerun`.
parse_args <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  defaults <- c(list(rerun = FALSE), defaults)
  pairs <- regmatches(args, regexec("^--([a-z_]+)(=(.*))?$", args))
  if (any(lengths(pairs) != 4)) stop("arguments must look like --name=value or --name")
  for (pair in pairs) {
    name <- pair[2]
    if (!name %in% names(defaults)) stop("unknown argument --", name)
    is_flag <- pair[3] == ""
    if (is_flag != is.logical(defaults[[name]])) {
      stop("--", name, if (is_flag) " needs a value" else " takes no value")
    }
    defaults[[name]] <- if (is_flag) TRUE else methods::as(pair[4], class(defaults[[name]]))
  }
  defaults
}

#' Turn an index text into indices.
#'
#' @param text "3:7" (3 to 7), "2,5,9", or a single number N.
#' @param count If TRUE a single number N means 1 to N.
#' @return Integer vector.
parse_index_set <- function(text, count = FALSE) {
  if (grepl("^[0-9]+:[0-9]+$", text)) {
    ends <- as.integer(strsplit(text, ":")[[1]])
    return(ends[1]:ends[2])
  }
  index <- as.integer(strsplit(text, ",")[[1]])
  if (count && length(index) == 1) seq_len(index) else index
}

#' Print a time-stamped message.
#'
#' @param msg Format of the message, as for `sprintf()`.
#' @param ... Values for the format.
#' @return NULL, invisibly.
log_message <- function(msg, ...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", sprintf(msg, ...), "\n")
}
