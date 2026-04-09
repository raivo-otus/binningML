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
  controller = crew_controller_local(workers = 4)
)

# Global parameters — tweak here, everything downstream updates
params <- list(
  min_class_n  = 30L, # Minimum samples per class for a dataset to be included
  n_bins_grid  = c(3, 4, 5, 7, 10, 15, 20), # Bin count search grid
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
  tar_target(data_file, "data/tses_85_280126.rds", format = "file"),
  tar_target(tse_list_raw, load_tse_list(data_file)),
  tar_target(tse_list, filter_tse_list(tse_list_raw, min_class_n)),

  # ── Transformations (Step 3 — dynamic branching placeholder) ──────────────
  # tar_target(transformation_names, c("binning", "rclr", "relabundance", "pa")),

  # ── Modeling (Steps 4–5 — wired in later) ─────────────────────────────────

  # ── Reporting ─────────────────────────────────────────────────────────────
  tar_quarto(
    report,
    path = "analysis/report.qmd",
    quiet = FALSE
  )
)
