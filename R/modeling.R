#' Modeling functions for microbiome ML benchmarking
#'
#' Thin wrappers around mikropml::run_ml. For each (dataset, transformation)
#' pair we hand mikropml a feature matrix + outcome and let it do the
#' train/test split, inner CV, RF fit, and metric computation.
#'
#' Binning is the only transformation with a tunable preprocessing parameter
#' (n_bins). We tune it by running mikropml once per candidate n_bins on the
#' same seed (so the same train/test split is reused), picking the n_bins
#' with the best train-CV performance, and reporting the corresponding
#' held-out test metrics.


#' Build a mikropml-ready data frame from a TSE assay
#'
#' @param tse        A TreeSummarizedExperiment.
#' @param assay_name Assay to pull as the feature matrix.
#' @param label_col  colData column holding the outcome.
#' @return data.frame with outcome in column 1, features in remaining columns.
tse_to_ml_df <- function(tse, assay_name, label_col) {
  x <- t(SummarizedExperiment::assay(tse, assay_name))
  y <- SummarizedExperiment::colData(tse)[[label_col]]
  data.frame(outcome = as.factor(y), as.data.frame(x), check.names = FALSE)
}


#' Run mikropml on a single (TSE, assay) combination
#'
#' @param tse        A TSE that already contains `assay_name`.
#' @param assay_name Assay to use as the feature matrix.
#' @param label_col  colData column with class labels.
#' @param kfold      Inner CV folds for hyperparameter tuning.
#' @param seed       Random seed (controls the train/test split).
#' @return The mikropml::run_ml() result list.
run_ml_on_assay <- function(tse,
                            assay_name,
                            method = "rf",
                            label_col = "group",
                            kfold = 5L,
                            seed = 42L) {
  df <- tse_to_ml_df(tse, assay_name = assay_name, label_col = label_col)
  preprocessed <- mikropml::preprocess_data(df, outcome_colname = "outcome")
  mikropml::run_ml(
    dataset         = preprocessed$dat_transformed,
    method          = method,
    outcome_colname = "outcome",
    kfold           = kfold,
    cv_times        = 1L,
    seed            = seed
  )
}


#' Tune n_bins and evaluate the binning transformation
#'
#' For each candidate n_bins we re-bin the full TSE and call run_ml with the
#' same seed (so the train/test split is identical across candidates). We
#' pick the n_bins with the best inner-CV performance and return its held-out
#' test metrics together with the full tuning trace.
#'
#' @param tse         A TSE with a "counts" assay.
#' @param n_bins_grid Integer vector of candidate bin counts.
#' @param label_col   colData column with class labels.
#' @param kfold       Inner CV folds.
#' @param seed        Random seed.
#' @return List with `best_n_bins`, `performance` (best model's row), and
#'   `tuning` (one row per candidate n_bins).
evaluate_binning <- function(tse,
                             n_bins_grid,
                             method = "rf",
                             label_col = "group",
                             kfold = 5L,
                             seed = 42L) {
  runs <- lapply(n_bins_grid, function(nb) {
    binned <- transform_binning(tse, n_bins = nb)
    run_ml_on_assay(binned, assay_name = "binning", method = method,
                    label_col = label_col, kfold = kfold, seed = seed)
  })

  tuning <- purrr::map2_dfr(n_bins_grid, runs, function(nb, r) {
    dplyr::mutate(r$performance, n_bins = nb, .before = 1)
  })

  best_idx <- which.max(tuning$cv_metric_AUC)
  list(
    best_n_bins = n_bins_grid[best_idx],
    performance = tuning[best_idx, ],
    tuning      = tuning
  )
}


#' Evaluate one transformation on one dataset
#'
#' Dispatches to `evaluate_binning()` for "binning" and to a single
#' `run_ml_on_assay()` call for the pre-computed assays.
#'
#' @param tse            A TSE from `tse_list_transformed`.
#' @param transformation One of "relabundance", "rclr", "pa", "binning".
#' @param n_bins_grid    Used only when transformation == "binning".
#' @param label_col      colData column with class labels.
#' @param kfold          Inner CV folds.
#' @param seed           Random seed.
#' @return A one-row tibble of test-set performance with `transformation`
#'   and `selected_n_bins` (NA for non-binning) columns prepended.
evaluate_transformation <- function(tse,
                                    transformation,
                                    method = "rf",
                                    n_bins_grid = c(3L, 4L, 5L, 7L, 10L,
                                                    15L, 20L),
                                    label_col = "group",
                                    kfold = 5L,
                                    seed = 42L) {
  if (transformation == "binning") {
    res <- evaluate_binning(tse, n_bins_grid = n_bins_grid, method = method,
                            label_col = label_col, kfold = kfold, seed = seed)
    dplyr::mutate(res$performance,
                  method          = method,
                  transformation  = "binning",
                  selected_n_bins = res$best_n_bins,
                  .before = 1)
  } else {
    r <- run_ml_on_assay(tse, assay_name = transformation, method = method,
                         label_col = label_col, kfold = kfold, seed = seed)
    dplyr::mutate(r$performance,
                  method          = method,
                  transformation  = transformation,
                  selected_n_bins = NA_integer_,
                  .before = 1)
  }
}
