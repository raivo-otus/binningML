# R/recipe_steps.R
# Custom recipes step: step_bin_samplewise()
#
# Wraps mia::transformAssay(method = "binning", MARGIN = "samples") as a
# tidymodels recipe step so that n_bins can be tuned via tune::tune_grid().
# Binning is sample-wise (each sample's non-zero feature values are ranked
# and assigned to quantile bins independently of other samples), so prep()
# has no cross-sample state to learn — it is a near-no-op that only resolves
# the selected columns.


# ── dials parameter ───────────────────────────────────────────────────────────

#' Number-of-bins dials parameter for step_bin_samplewise
#'
#' @param range Integer vector of length 2. Lower and upper bound.
#' @return A `quant_param` object.
n_bins_param <- function(range = c(1L, 50L)) {
  dials::new_quant_param(
    type      = "integer",
    range     = log2(range),
    inclusive = c(TRUE, TRUE),
    trans     = scales::log2_trans(),
    label     = c(n_bins = "Number of Quantile Bins"),
    finalize  = NULL
  )
}


# ── step constructor ──────────────────────────────────────────────────────────

#' Sample-wise quantile binning recipe step
#'
#' For each sample, non-zero predictor values are split into `n_bins`
#' equal-count quantile groups; zeros map to bin 0. Because the breaks are
#' derived per-sample from that sample's own values, no information is shared
#' across samples and the step is leakage-free by construction. `n_bins` is
#' exposed as a tunable parameter for `tune::tune_grid()`.
#'
#' @param recipe   A `recipe` object.
#' @param ...      Column selector (tidyselect); usually `all_predictors()`.
#' @param n_bins   Integer or `tune::tune()`. Number of quantile bins.
#' @param role     Not used; retained for recipes API compatibility.
#' @param trained  Logical. Set to `TRUE` by `prep()`.
#' @param columns  Character. Resolved column names, set by `prep()`.
#' @param skip     Logical. Whether to skip during `bake(new_data = testing())`.
#' @param id       Character. Unique step identifier.
#' @return The recipe with this step appended.
step_bin_samplewise <- function(recipe,
                                ...,
                                n_bins  = tune::tune(),
                                role    = NA,
                                trained = FALSE,
                                columns = NULL,
                                skip    = FALSE,
                                id      = recipes::rand_id("bin_samplewise")) {
  recipes::add_step(
    recipe,
    new_step_bin_samplewise(
      terms   = rlang::enquos(...),
      n_bins  = n_bins,
      role    = role,
      trained = trained,
      columns = columns,
      skip    = skip,
      id      = id
    )
  )
}

new_step_bin_samplewise <- function(terms, n_bins, role, trained,
                                    columns, skip, id) {
  recipes::step(
    terms    = terms,
    n_bins   = n_bins,
    role     = role,
    trained  = trained,
    columns  = columns,
    skip     = skip,
    id       = id,
    subclass = "bin_samplewise"
  )
}


# ── prep ─────────────────────────────────────────────────────────────────────

#' @export
prep.step_bin_samplewise <- function(x, training, info = NULL, ...) {
  col_names <- recipes::recipes_eval_select(x$terms, training, info)
  recipes::check_type(training[, col_names], quant = TRUE)
  new_step_bin_samplewise(
    terms   = x$terms,
    n_bins  = x$n_bins,
    role    = x$role,
    trained = TRUE,
    columns = col_names,
    skip    = x$skip,
    id      = x$id
  )
}


# ── bake ─────────────────────────────────────────────────────────────────────

#' @export
bake.step_bin_samplewise <- function(object, new_data, ...) {
  cols <- object$columns
  if (length(cols) == 0L) return(new_data)

  # Transpose to features x samples for mia (SE assay convention)
  pred_mat <- t(as.matrix(new_data[, cols, drop = FALSE]))

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = pred_mat)
  )
  se <- mia::transformAssay(
    se,
    method     = "binning",
    assay.type = "counts",
    name       = "binning",
    nbins      = as.integer(object$n_bins),
    MARGIN     = "samples"
  )

  # Transpose back to samples x features
  binned <- t(SummarizedExperiment::assay(se, "binning"))

  for (j in seq_along(cols)) {
    new_data[[cols[j]]] <- binned[, j]
  }
  tibble::as_tibble(new_data)
}


# ── print / tidy / tunable ───────────────────────────────────────────────────

#' @export
print.step_bin_samplewise <- function(x, width = max(20, getOption("width") - 30),
                                      ...) {
  title <- "Sample-wise quantile binning"
  recipes::print_step(x$columns, x$terms, x$trained, title, width)
  invisible(x)
}

#' @export
tidy.step_bin_samplewise <- function(x, ...) {
  tibble::tibble(n_bins = x$n_bins, id = x$id)
}

#' @export
tunable.step_bin_samplewise <- function(x, ...) {
  tibble::tibble(
    name         = "n_bins",
    call_info    = list(list(pkg = NULL, fun = "n_bins_param")),
    source       = "recipe",
    component    = "step_bin_samplewise",
    component_id = x$id
  )
}
