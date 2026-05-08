# R/transformations.R
# Wrappers around mia::transformAssay for all four transformations.
#
# Shared signature: transform_*(tse, assay_type = "counts", name = <method>, ...)
# Each returns the TSE with a new assay added.
#
# All four transformations operate sample-wise (MARGIN = "samples"): each
# sample's transformed values depend only on that sample's own counts, with no
# cross-sample statistics. They are therefore safe to compute on a full TSE
# before train/test splitting. 

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
#' Sample-wise: for each sample, non-zero counts are split into `n_bins`
#' equal-count quantile groups; zeros map to bin 0. Because the breaks for
#' each sample are computed from that sample's own values, the transformation
#' itself does not leak information across samples — but `n_bins` is a
#' hyperparameter and must still be tuned via inner CV on the training
#' partition only (see modeling code).
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

#' Add sample-wise transformed assays to every TSE in a list
#'
#' Applies TSS, rCLR, and presence/absence to each TSE. Binning is excluded
#' because `n_bins` is tuned per (dataset, model) inside the modeling path;
#' the binning transformation itself is sample-wise and safe to apply at any
#' stage, but pre-computing it for a single fixed `n_bins` would defeat the
#' tuning grid.
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
