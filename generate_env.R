library(rix)

rix(
  r_ver = "latest-upstream",
  r_pkgs = c(
    "crew",
    "devtools",
    "dials",
    "dplyr",
    "ggExtra",
    "ggplot2",
    "ggsci",
    "glmnet",
    "kernlab",
    "mia",
    "miaViz",
    "patchwork",
    "parsnip",
    "purrr",
    "qs2",
    "quarto",
    "ranger",
    "recipes",
    "rix",
    "rlang",
    "rsample",
    "scater",
    "stringr",
    "tarchetypes",
    "targets",
    "tibble",
    "tidyr",
    "tune",
    "usethis",
    "workflows",
    "yardstick"
  ),
  system_pkgs = c(
    "quarto"
  ),
  git_pkgs = NULL,
  ide = "code",
  project_path = ".",
  overwrite = TRUE
)
