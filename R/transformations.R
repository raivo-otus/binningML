# R/transformations.R
# Wrappers around mia::transformAssay for all four transformations.
#
# Shared signature: transform_*(tse, assay_type = "counts", name = <method>, ...)
# Each returns the TSE with a new assay added.
#
# Binning note: we tune `n_bins` (not the break points themselves), so it is
# fine to call transform_binning() inside CV folds — applied separately to
# train and test partitions.


# ── Standard transformations ──────────────────────────────────────────────────

#' TSS (relative abundance) transformation
#'
#' @param tse        A TreeSummarizedExperiment.
#' @param assay_type Character. Source assay name. Default "counts".
#' @param name       Character. Name for the new assay. Default "relabundance".
#' @return TSE with a new "relabundance" assay (or `name`).
transform_tss <- function(tse, assay_type = "counts", name = "relabundance") {
  mia::transformAssay(
    tse,
    method     = "relabundance",
    assay.type = assay_type,
    name       = name,
    MARGIN     = "samples"
  )
}


#' rCLR (robust centred log-ratio) transformation
#'
#' Zero-aware; no pseudocount needed. Requires vegan >= 2.7-1.
#'
#' @param tse        A TreeSummarizedExperiment.
#' @param assay_type Character. Source assay name. Default "counts".
#' @param name       Character. Name for the new assay. Default "rclr".
#' @return TSE with a new "rclr" assay (or `name`).
transform_rclr <- function(tse, assay_type = "counts", name = "rclr") {
  mia::transformAssay(
    tse,
    method     = "rclr",
    assay.type = assay_type,
    name       = name,
    MARGIN     = "samples"
  )
}


#' Presence/absence transformation
#'
#' @param tse        A TreeSummarizedExperiment.
#' @param assay_type Character. Source assay name. Default "counts".
#' @param name       Character. Name for the new assay. Default "pa".
#' @return TSE with a new "pa" assay (or `name`).
transform_pa <- function(tse, assay_type = "counts", name = "pa") {
  mia::transformAssay(
    tse,
    method     = "pa",
    assay.type = assay_type,
    name       = name,
    MARGIN     = "samples"
  )
}


#' Quantile binning transformation
#'
#' Zeros are assigned to bin 0; non-zero values are split into `n_bins`
#' equal-count quantile groups per feature. Safe to call inside CV folds
#' since we only tune `n_bins`, not the break points themselves.
#'
#' @param tse        A TreeSummarizedExperiment.
#' @param n_bins     Integer. Number of quantile bins for non-zero values.
#' @param assay_type Character. Source assay name. Default "counts".
#' @param name       Character. Name for the new assay. Default "binning".
#' @return TSE with a new "binning" assay (or `name`).
transform_binning <- function(tse, n_bins, assay_type = "counts", name = "binning") {
  stopifnot(is.numeric(n_bins), n_bins >= 2)
  mia::transformAssay(
    tse,
    method     = "binning",
    assay.type = assay_type,
    name       = name,
    nbins      = as.integer(n_bins),
    MARGIN     = "samples"
  )
}


# ── Batch helper ──────────────────────────────────────────────────────────────

#' Add all leakage-free transformed assays to every TSE in a list
#'
#' Applies TSS, rCLR, and presence/absence to each TSE. Binning is
#' deliberately excluded — it requires tuning inside the CV train/test
#' split (see compute_bin_breaks() / apply_bin_breaks()).
#'
#' @param tse_list Named list of TreeSummarizedExperiment objects.
#' @return Named list of TSEs with added assays: "relabundance", "rclr", "pa".
add_all_transformations <- function(tse_list) {
  lapply(tse_list, function(tse) {
    tse <- mia::agglomerateByRank(tse, rank = "genus")
    tse <- transform_tss(tse)
    tse <- transform_rclr(tse)
    tse <- transform_pa(tse)
    tse
  })
}