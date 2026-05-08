# R/modeling.R
# Tidymodels-based modeling for microbiome ML benchmarking.
#
# For each (dataset, transformation) pair, we build a workflow with a recipe
# (preprocessing) and a parsnip model spec, tune hyperparameters via inner
# k-fold CV on the training partition, and evaluate once on the held-out test
# partition via tune::last_fit().
#
# Binning uses step_bin_samplewise() from R/recipe_steps.R; n_bins is tuned
# jointly with the model hyperparameters in the inner CV so the selected
# n_bins is never informed by test-set performance.


#' Build a samples x features tibble from a TSE assay
#'
#' @param tse        A TreeSummarizedExperiment.
#' @param assay_name Assay to extract as the feature matrix.
#' @param label_col  colData column holding the outcome.
#' @return Tibble with `outcome` as first column (factor), then one column per
#'   feature. Column names are taxon names from the assay rownames.
tse_to_tibble <- function(tse, assay_name, label_col) {
  x <- t(SummarizedExperiment::assay(tse, assay_name))
  y <- SummarizedExperiment::colData(tse)[[label_col]]
  dplyr::bind_cols(
    tibble::tibble(outcome = as.factor(y)),
    tibble::as_tibble(as.data.frame(x, check.names = FALSE))
  )
}


#' Build a parsnip model spec for a given classifier name
#'
#' @param method Character. `"rf"` (random forest via ranger) or `"glmnet"`.
#' @return A parsnip model specification with tunable hyperparameters.
model_spec <- function(method) {
  switch(method,
    rf = parsnip::rand_forest(mtry = tune::tune(), min_n = tune::tune(),
                              trees = 500L) |>
      parsnip::set_engine("ranger") |>
      parsnip::set_mode("classification"),
    glmnet = parsnip::logistic_reg(penalty = tune::tune(),
                                   mixture = tune::tune()) |>
      parsnip::set_engine("glmnet") |>
      parsnip::set_mode("classification"),
    svm_rbf = parsnip::svm_rbf(cost = tune::tune(), rbf_sigma = tune::tune()) |>
      parsnip::set_engine("kernlab") |>
      parsnip::set_mode("classification"),
    stop("Unsupported method: ", method)
  )
}


#' Yardstick metric set used throughout the modeling path
#'
#' @return A `metric_set` with roc_auc, accuracy, and kap.
metric_set_clf <- function() {
  yardstick::metric_set(yardstick::roc_auc, yardstick::accuracy, yardstick::kap)
}


#' Build a tuning grid for a given classifier, optionally crossed with n_bins
#'
#' For RF, `mtry` is bounded by the number of predictor columns in `train_x`
#' (computed before any nzv filtering, so values are conservative). For glmnet,
#' penalty and mixture are varied on a regular log/linear grid.
#'
#' @param method       `"rf"` or `"glmnet"`.
#' @param train_x      Data frame of predictor columns from the training set.
#' @param n_bins_levels Integer or `NULL`. If provided, dials generates this
#'   many log2-spaced n_bins candidates via `grid_regular(n_bins_param())`.
#' @return A data frame of hyperparameter combinations.
make_tune_grid <- function(method, train_x, n_bins_levels = NULL, n_bins_range = NULL) {
  p <- ncol(train_x)

  model_grid <- switch(method,
    rf = {
      mtry_vals <- unique(as.integer(c(
        max(1L, floor(sqrt(p))),
        max(1L, floor(p / 3L)),
        max(1L, floor(p / 2L))
      )))
      tidyr::expand_grid(mtry = mtry_vals, min_n = c(1L, 5L, 10L))
    },
    glmnet = dials::grid_regular(
      dials::penalty(range = c(-4, 0)),
      dials::mixture(),
      levels = c(penalty = 20L, mixture = 5L)
    ),
    svm_rbf = dials::grid_regular(
      dials::cost(),
      dials::rbf_sigma(),
      levels = c(cost = 5L, rbf_sigma = 5L)
    ),
    stop("Unsupported method: ", method)
  )

  if (!is.null(n_bins_levels)) {
    bins_param <- if (!is.null(n_bins_range)) n_bins_param(range = n_bins_range) else n_bins_param()
    tidyr::crossing(dials::grid_regular(bins_param, levels = n_bins_levels),
                    model_grid)
  } else {
    model_grid
  }
}


#' Run tidymodels on a single (TSE, assay) combination
#'
#' Outer split: stratified 80/20. Inner CV: k-fold stratified on the training
#' partition. Preprocessing: zero-variance filtering + normalisation (fit on each inner
#' training fold only). Classifier hyperparameters are tuned by inner-CV AUC;
#' the best set is refitted on the full training partition and evaluated once
#' on the held-out test partition.
#'
#' @param tse        A TSE that already contains `assay_name`.
#' @param assay_name Assay to use as the feature matrix.
#' @param method     `"rf"` or `"glmnet"`.
#' @param label_col  colData column with class labels.
#' @param kfold      Number of inner CV folds.
#' @param seed       Random seed (controls the outer split and inner folds).
#' @return One-row tibble with columns `Accuracy`, `Kappa`, `AUC`.
run_tidymodels_on_assay <- function(tse,
                                    assay_name,
                                    method    = "rf",
                                    label_col = "group",
                                    kfold     = 5L,
                                    seed      = 42L) {
  df <- tse_to_tibble(tse, assay_name = assay_name, label_col = label_col)

  set.seed(seed)
  split <- rsample::initial_split(df, prop = 0.8, strata = "outcome")
  train <- rsample::training(split)
  folds <- rsample::vfold_cv(train, v = kfold, strata = "outcome")

  rec <- recipes::recipe(outcome ~ ., data = train) |>
    recipes::step_zv(recipes::all_predictors())
  if (method %in% c("glmnet", "svm_rbf")) {
    rec <- rec |> recipes::step_normalize(recipes::all_predictors())
  }

  wf <- workflows::workflow(rec, model_spec(method))

  train_x <- recipes::bake(
    recipes::prep(rec, training = train),
    new_data = NULL, recipes::all_predictors()
  )
  grid <- make_tune_grid(method, train_x)

  tuned <- tune::tune_grid(
    wf,
    resamples = folds,
    grid      = grid,
    metrics   = metric_set_clf(),
    control   = tune::control_grid(save_pred = FALSE, verbose = FALSE)
  )

  best     <- tune::select_best(tuned, metric = "roc_auc")
  final_wf <- tune::finalize_workflow(wf, best)
  final    <- tune::last_fit(final_wf, split = split, metrics = metric_set_clf())

  tune::collect_metrics(final) |>
    dplyr::select(.metric, .estimate) |>
    tidyr::pivot_wider(names_from = .metric, values_from = .estimate) |>
    dplyr::rename(Accuracy = accuracy, Kappa = kap, AUC = roc_auc)
}


#' Tune n_bins and evaluate the binning transformation via tidymodels
#'
#' Uses `step_bin_samplewise(n_bins = tune())` in the recipe so that
#' `tune::tune_grid()` selects `n_bins` from inner-CV AUC on the training
#' partition only — the test partition is never seen during selection.
#' The binning transformation itself is sample-wise (no cross-sample state),
#' so the mia call inside `bake()` is leakage-free regardless of when it runs.
#'
#' @param tse          A TSE with a "counts" assay (raw, genus-agglomerated).
#' @param n_bins_levels Integer. Number of log2-spaced n_bins candidates
#'   generated by `dials::grid_regular(n_bins_param())`.
#' @param method       `"rf"` or `"glmnet"`.
#' @param label_col    colData column with class labels.
#' @param kfold        Number of inner CV folds.
#' @param seed         Random seed.
#' @return List with `best_n_bins` (integer), `performance` (one-row tibble
#'   with `Accuracy`, `Kappa`, `AUC`), and `tuning` (inner-CV roc_auc per
#'   combination of n_bins and model hyperparameters).
evaluate_binning_tidymodels <- function(tse,
                                        n_bins_levels,
                                        n_bins_range  = NULL,
                                        method    = "rf",
                                        label_col = "group",
                                        kfold     = 5L,
                                        seed      = 42L) {
  df <- tse_to_tibble(tse, assay_name = "counts", label_col = label_col)

  set.seed(seed)
  split <- rsample::initial_split(df, prop = 0.8, strata = "outcome")
  train <- rsample::training(split)
  folds <- rsample::vfold_cv(train, v = kfold, strata = "outcome")

  rec <- recipes::recipe(outcome ~ ., data = train) |>
    step_bin_samplewise(recipes::all_predictors(), n_bins = tune::tune()) |>
    recipes::step_zv(recipes::all_predictors())
  if (method %in% c("glmnet", "svm_rbf")) {
    rec <- rec |> recipes::step_normalize(recipes::all_predictors())
  }

  wf <- workflows::workflow(rec, model_spec(method))

  rec_sizing <- recipes::recipe(outcome ~ ., data = train) |>
    step_bin_samplewise(recipes::all_predictors(), n_bins = 5L) |>
    recipes::step_zv(recipes::all_predictors())
  train_x <- recipes::bake(
    recipes::prep(rec_sizing, training = train),
    new_data = NULL, recipes::all_predictors()
  )
  grid <- make_tune_grid(method, train_x, n_bins_levels = n_bins_levels, n_bins_range = n_bins_range)

  tuned <- tune::tune_grid(
    wf,
    resamples = folds,
    grid      = grid,
    metrics   = metric_set_clf(),
    control   = tune::control_grid(save_pred = FALSE, verbose = FALSE)
  )

  best     <- tune::select_best(tuned, metric = "roc_auc")
  final_wf <- tune::finalize_workflow(wf, best)
  final    <- tune::last_fit(final_wf, split = split, metrics = metric_set_clf())

  perf <- tune::collect_metrics(final) |>
    dplyr::select(.metric, .estimate) |>
    tidyr::pivot_wider(names_from = .metric, values_from = .estimate) |>
    dplyr::rename(Accuracy = accuracy, Kappa = kap, AUC = roc_auc)

  list(
    best_n_bins = best$n_bins,
    performance = perf,
    tuning      = tune::collect_metrics(tuned) |>
      dplyr::filter(.metric == "roc_auc")
  )
}


#' Evaluate one transformation on one dataset
#'
#' Dispatches to `evaluate_binning_tidymodels()` for "binning" and to
#' `run_tidymodels_on_assay()` for the pre-computed assays.
#'
#' @param tse            A TSE from `tse_list_transformed`.
#' @param transformation One of "relabundance", "rclr", "pa", "binning".
#' @param method         `"rf"` or `"glmnet"`.
#' @param n_bins_levels  Used only when transformation == "binning". Number of
#'   log2-spaced n_bins candidates for `dials::grid_regular(n_bins_param())`.
#' @param label_col      colData column with class labels.
#' @param kfold          Inner CV folds.
#' @param seed           Random seed.
#' @return A one-row tibble with columns `method`, `transformation`,
#'   `selected_n_bins`, `Accuracy`, `Kappa`, `AUC`.
evaluate_transformation <- function(tse,
                                    transformation,
                                    method        = "rf",
                                    n_bins_levels = 6L,
                                    n_bins_range  = NULL,
                                    label_col     = "group",
                                    kfold         = 5L,
                                    seed          = 42L) {
  if (transformation == "binning") {
    res <- evaluate_binning_tidymodels(
      tse, n_bins_levels = n_bins_levels, n_bins_range = n_bins_range,
      method = method, label_col = label_col, kfold = kfold, seed = seed
    )
    dplyr::mutate(res$performance,
                  method          = method,
                  transformation  = "binning",
                  selected_n_bins = res$best_n_bins,
                  .before         = 1)
  } else {
    perf <- run_tidymodels_on_assay(
      tse, assay_name = transformation, method = method,
      label_col = label_col, kfold = kfold, seed = seed
    )
    dplyr::mutate(perf,
                  method          = method,
                  transformation  = transformation,
                  selected_n_bins = NA_integer_,
                  .before         = 1)
  }
}
