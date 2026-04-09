#' Load the full curated list of TreeSummarizedExperiment objects
#'
#' @return A named list of TSE objects with colData containing case/control labels.
load_tse_list <- function(file) {
  e <- new.env()
  load(file, envir = e)
  e$tses_raw
}


#' Filter TSE list to datasets with sufficient per-class sample size
#'
#' Drops any dataset where at least one class has fewer than `min_n` samples.
#' Expects a column named `group` in colData with exactly two levels
#' (case / control).
#'
#' @param tse_list Named list of TSE objects.
#' @param min_n    Integer. Minimum number of samples required per class.
#'                 Datasets failing this threshold are dropped entirely.
#' @param label_col Character. Name of the colData column holding class labels.
#'                  Defaults to "group".
#'
#' @return A filtered named list of TSE objects. Emits a message summarising
#'         which datasets were kept and which were dropped.
filter_tse_list <- function(tse_list, min_n, label_col = "group") {
  stopifnot(is.list(tse_list), length(tse_list) > 0)
  if (is.null(names(tse_list))) {
    names(tse_list) <- paste0("dataset_", seq_along(tse_list))
  }

  has_enough_samples <- function(tse, name) {
    labels <- SummarizedExperiment::colData(tse)[[label_col]]

    if (is.null(labels)) {
      warning(sprintf("Dataset '%s': column '%s' not found — dropping.", name, label_col))
      return(FALSE)
    }

    counts <- table(labels)

    if (length(counts) != 2L) {
      warning(sprintf(
        "Dataset '%s': expected 2 classes, found %d — dropping.", name, length(counts)
      ))
      return(FALSE)
    }

    all(counts >= min_n)
  }

  keep <- vapply(
    seq_along(tse_list),
    function(i) isTRUE(has_enough_samples(tse_list[[i]], names(tse_list)[[i]])),
    logical(1)
  )
  names(keep) <- names(tse_list)

  dropped <- names(keep)[!keep]
  retained <- names(keep)[keep]

  if (length(dropped) > 0) {
    message(sprintf(
      "Dropped %d dataset(s) (min per-class n < %d): %s",
      length(dropped), min_n, paste(dropped, collapse = ", ")
    ))
  }
  message(sprintf("Retained %d dataset(s): %s", length(retained), paste(retained, collapse = ", ")))

  tse_list[keep]
}


#' Summarise dataset metadata for reporting
#'
#' Extracts a tidy one-row-per-dataset tibble with basic metadata.
#' Useful as a pipeline target for the Quarto report overview table.
#'
#' @param tse_list Named, filtered list of TSE objects.
#' @param label_col Character. colData column holding class labels.
#'
#' @return A tibble with columns: dataset, n_samples, n_features, n_case, n_control.
summarise_datasets <- function(tse_list, label_col = "group") {
  purrr::imap(tse_list, function(tse, name) {
    labels <- SummarizedExperiment::colData(tse)[[label_col]]
    counts <- table(labels)
    tibble::tibble(
      dataset    = name,
      n_samples  = ncol(tse),
      n_features = nrow(tse),
      n_case     = counts[[2]],
      n_control  = counts[[1]]
    )
  }) |>
    purrr::list_rbind()
}
