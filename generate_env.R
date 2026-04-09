library(rix)

rix(
  r_ver = "latest-upstream",
  r_pkgs = c(
    "devtools",
    "dplyr",
    "ggExtra",
    "ggplot2",
    "ggsci",
    "mia",
    "miaViz",
    "targets",
    "tarchetypes",
    "mikropml",
    "ranger",
    "patchwork",
    "quarto",
    "rix",
    "scater",
    "stringr",
    "tidyr",
    "purrr",
    "crew",
    "usethis"
  ),
  system_pkgs = c(
    "quarto"
  ),
  git_pkgs = NULL,
  ide = "code",
  project_path = ".",
  overwrite = TRUE
)
