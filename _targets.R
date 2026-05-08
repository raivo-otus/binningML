library(targets)
library(tarchetypes)
library(crew)

# Source all function scripts from R/
tar_source("R/")

# Pipeline-wide options
tar_option_set(
  packages = c(
    "dials",
    "dplyr",
    "ggsci",
    "ggplot2",
    "mia",
    "parsnip",
    "patchwork",
    "purrr",
    "kernlab",
    "ranger",
    "recipes",
    "rlang",
    "qs2",
    "rsample",
    "tibble",
    "tidyr",
    "tune",
    "workflows",
    "yardstick"
  ),
  format = "qs",  # Faster then rds for object cache
  controller = crew_controller_local(workers = 6)
)

# Global parameters — tweak here, everything downstream updates
params <- list(
  min_class_n  = 50L,        # Minimum samples per class for a dataset to be included
  n_bins_levels = 9L,        # Log2-spaced n_bins candidates
  n_bins_range = c(1L, 50L), # Min and max n_bins explored during tuning
  cv_folds     = 5L,         # Folds for both outer and inner CV loops
  random_seed  = 42L
)

list(
  # ── Parameters ────────────────────────────────────────────────────────────
  tar_target(min_class_n, params$min_class_n),
  tar_target(n_bins_levels, params$n_bins_levels),
  tar_target(n_bins_range, params$n_bins_range),
  tar_target(cv_folds, params$cv_folds),
  tar_target(random_seed, params$random_seed),

  # ── Data –––––––––––––––––––––––––─────────────────────────────────────────
  tar_target(data_file, "data/tses_85_280126.Rdata", format = "file"),
  tar_target(tse_list_raw, load_tse_list(data_file)),
  tar_target(tse_list, filter_tse_list(tse_list_raw, min_class_n)),

  # ── Transformations ───────────────────────────────────────────────────────
  # One list of TSEs; each TSE gains assays: relabundance, rclr, pa.
  tar_target(tse_list_transformed, add_all_transformations(tse_list)),

  # ── Modeling ──────────────────────────────────────────────────────────────
  # Dynamic branching over the (dataset, transformation) grid. Each branch
  # calls evaluate_transformation(); the binning branch tunes n_bins via
  # inner CV (tune::tune_grid) so n_bins selection never sees the test set.
  tar_target(transformation_names,
             c("relabundance", "rclr", "pa", "binning")),
  tar_target(classifier_names, c("rf", "glmnet", "svm_rbf")),
  tar_target(
    model_grid,
    tidyr::expand_grid(
      dataset_name   = names(tse_list_transformed),
      transformation = transformation_names,
      classifier     = classifier_names
    )
  ),
  tar_target(
    model_results,
    evaluate_transformation(
      tse            = tse_list_transformed[[model_grid$dataset_name]],
      transformation = model_grid$transformation,
      method         = model_grid$classifier,
      n_bins_levels  = n_bins_levels,
      n_bins_range   = n_bins_range,
      kfold          = cv_folds,
      seed           = random_seed
    ) |>
      dplyr::mutate(dataset = model_grid$dataset_name, .before = 1),
    pattern = map(model_grid)
  ),


  # ── Reporting ─────────────────────────────────────────────────────────────
  tar_quarto(report, path = "analysis/report.qmd", quiet = FALSE)

)
