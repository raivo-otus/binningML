library(targets)
library(tarchetypes)
library(crew)

# Source all function scripts from R/
tar_source("R/")

# Pipeline-wide options
tar_option_set(
  packages = c(
    "dplyr",
    "ggsci",
    "ggplot2",
    "mia",
    "mikropml",
    "ranger",
    "patchwork",
    "purrr",
    "tidyr"
  ),
  format = "rds",
  # controller = crew_controller_local(workers = 4) # disable parallel for now
)

# Global parameters — tweak here, everything downstream updates
params <- list(
  min_class_n  = 50L, # Minimum samples per class for a dataset to be included
  n_bins_grid  = c(1, 2, 3, 5, 8, 13, 21, 34, 55), # Bin count search grid
  cv_folds     = 5L, # Folds for both outer and inner CV loops
  random_seed  = 42L
)

list(
  # ── Parameters ────────────────────────────────────────────────────────────
  tar_target(min_class_n, params$min_class_n),
  tar_target(n_bins_grid, params$n_bins_grid),
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
  # runs mikropml::run_ml; the binning branch additionally tunes n_bins.
  tar_target(transformation_names,
             c("relabundance", "rclr", "pa", "binning")),
  tar_target(classifier_names, c("rf", "rpart2", "svmRadial")),
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
      n_bins_grid    = n_bins_grid,
      kfold          = cv_folds,
      seed           = random_seed
    ) |>
      dplyr::mutate(dataset = model_grid$dataset_name, .before = 1),
    pattern = map(model_grid)
  ),


  # ── Reporting ─────────────────────────────────────────────────────────────
  tar_quarto(
    report,
    path = "analysis/report.qmd",
    quiet = FALSE
  )
)
