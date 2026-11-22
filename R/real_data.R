# -----------------------------------------------------------------------------
# File:    R/real_data.R
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Table 3, Web Table S16.
# Purpose: Protein domain data: mutation counts per position of six domains,
#          the null fits of each domain, and the tables of the numbers of
#          rejected positions and of the estimates.
# Inputs:  data/realdata.txt (columns id and one per domain, NA padded).
# Outputs: None (the scripts write the tables).
# Usage:   source("R/setup.R"); protein_rejections(PATHS$data, 0.10)
# Source:  FINAL R CODES.R (RealDataFDR and its calls for k = 1, ..., 6).
# -----------------------------------------------------------------------------

#' Mutation counts of the protein domains.
#'
#' @param path File with one column per domain after the `id` column.
#' @return Named list (one element per domain) of count vectors.
read_domains <- function(path) {
  protein <- read.delim(require_data(path), sep = "", header = TRUE)
  lapply(protein[-1], function(counts) counts[!is.na(counts)])
}

#' Null fits of one domain; the bound D is at most the largest count.
#'
#' @param x Counts of the domain.
#' @return Data frame as returned by `select_null_fits()`.
fit_domain <- function(x) {
  select_null_fits(x, cap = max(x))
}

#' Numbers of rejected positions of every domain (Table 3).
#'
#' @param path Data file.
#' @param level FDR level.
#' @return Data frame with the columns `domain`, `cutoff`, `procedure`
#'   ("one_stage", "two_stage", "storey") and one column per null model.
protein_rejections <- function(path, level) {
  procedures <- list(one_stage = one_stage_count, two_stage = two_stage_count,
                     storey = storey_count)
  domains <- read_domains(path)
  do.call(rbind, lapply(names(domains), function(domain) {
    x <- domains[[domain]]
    fits <- fit_domain(x)
    cells <- expand.grid(procedure = names(procedures), cutoff = CUTOFFS,
                         stringsAsFactors = FALSE)
    counts <- t(vapply(seq_len(nrow(cells)), function(i) {
      vapply(MODELS, function(model) {
        fit <- fits[fits$cutoff == cells$cutoff[i] & fits$model == model, ]
        procedures[[cells$procedure[i]]](x, fit, level)
      }, numeric(1))
    }, numeric(length(MODELS))))
    data.frame(domain = domain, cells[c("cutoff", "procedure")], counts)
  }))
}

#' Null estimates of every domain (Web Table S16).
#'
#' @param path Data file.
#' @return Data frame with the columns `domain`, `cutoff`, `model`, `eta`,
#'   `lambda`, `theta`, `pi0`, `C`, `D`; NA for a parameter the model does
#'   not have.
protein_estimates <- function(path) {
  domains <- read_domains(path)
  do.call(rbind, lapply(names(domains), function(domain) {
    data.frame(domain = domain, mask_absent_parameters(fit_domain(domains[[domain]])),
               row.names = NULL)
  }))
}
