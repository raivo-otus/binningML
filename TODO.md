# Migration: mikropml → caret + fix data leakage

## Problem

Two sources of data leakage in the current pipeline:

### 1. Binning applied to the full dataset before train/test split

In `evaluate_binning()` ([R/modeling.R:74-78](R/modeling.R#L74)), `transform_binning(tse, n_bins)` is called on the entire TSE — all samples including the held-out test set. The quantile bin breakpoints are therefore computed from the full sample distribution. When `run_ml()` later splits the pre-binned data, the test set was already binned using its own quantile information. This inflates performance estimates for binning relative to the other transformations.

### 2. `mikropml::preprocess_data()` fitted on the full dataset

In `run_ml_on_assay()` ([R/modeling.R:42](R/modeling.R#L42)), `preprocess_data(df)` (near-zero-variance filtering, centering, scaling) is called on the full df before `run_ml()` does the train/test split. Feature selection and scaling statistics therefore incorporate test-set information for all four transformations, not just binning.

---

## Proposed Fix

Replace mikropml with caret and implement a proper two-level CV structure:

```
Outer split (fixed seed, stratified 80/20)
  └── Training set
        ├── Inner k-fold CV (n_bins tuning)
        │     For each fold:
        │       fit_bin_breaks()  ← on inner-train samples only
        │       apply to inner-val
        │       fit classifier (caret::train)
        │       record val AUC per n_bins candidate
        │     best_n_bins = argmax(mean val AUC across folds)
        │
        └── Final refit on full training set with best_n_bins
              fit_bin_breaks()  ← on full training set only
              apply to test set using those breaks
              fit caret model
              evaluate on test set → report Accuracy, Kappa, AUC
```

For the non-binning transformations (relabundance, rclr, pa), which are already
pre-computed on the full dataset but are sample-wise transforms with no
cross-sample statistics, the simplified fix is:

```
Outer split (same fixed seed)
  Training set → preProcess(method = c("nzv", "center", "scale"))  ← fit here
  Apply preProcess to test set
  caret::train() with inner CV for classifier hyperparams
  Evaluate on test set
```

---

## Tasks

### Phase 1 — Low-level binning utilities (R/transformations.R)

- [ ] Add `fit_bin_breaks(count_mat, n_bins)`:
  - Input: `count_mat` is features × samples (raw counts, training samples only)
  - For each feature: compute `n_bins` quantile breakpoints from non-zero values only
  - Zeros are special (bin 0); quantiles are fit on `x[x > 0]`
  - Return a named list of break vectors (one per feature; `NULL` if feature is all-zero)

- [ ] Add `apply_bin_breaks(count_mat, breaks)`:
  - Input: `count_mat` features × samples, `breaks` from `fit_bin_breaks()`
  - Zeros → 0L; non-zeros → `cut()` using the supplied breaks (no refitting)
  - Return integer matrix features × samples

- [ ] Remove or deprecate `transform_binning()` — it will no longer be called from
  the modeling path (but keep it if mia-based binning is still used elsewhere)

### Phase 2 — Replace mikropml modeling with caret (R/modeling.R)

- [ ] Remove `mikropml` from `tar_option_set(packages = ...)` in `_targets.R`; add
  `caret`, `pROC`, and the backend packages (`ranger`, `glmnet` already present)

- [ ] Add helper `make_train_control(kfold)` returning a `caret::trainControl` with:
  - `method = "cv"`, `number = kfold`
  - `classProbs = TRUE`, `summaryFunction = twoClassSummary`
  - `savePredictions = "final"`

- [ ] Add helper `caret_method_name(method)` mapping `"rf"` → `"ranger"` and
  `"glmnet"` → `"glmnet"` (caret method strings)

- [ ] Replace `run_ml_on_assay()` with `run_caret_on_assay(tse, assay_name, method, label_col, kfold, seed)`:
  1. Extract feature matrix X (samples × features) and outcome y from the TSE assay
  2. Stratified outer split via `caret::createDataPartition(y, p = 0.8)` with `set.seed(seed)`
  3. `preProcess(train_X, method = c("nzv", "center", "scale"))` — fit on train only
  4. Apply preProcess to both train and test X
  5. `caret::train(x = train_X, y = y_train, method = ..., trControl = make_train_control(kfold))`
  6. `predict()` on test_X; compute Accuracy, Kappa, AUC via `caret::confusionMatrix()` + `pROC::roc()`
  7. Return a one-row tibble matching the current `r$performance` column layout

- [ ] Replace `evaluate_binning()` with `evaluate_binning_caret(tse, n_bins_grid, method, label_col, kfold, seed)`:
  1. Extract raw counts matrix and outcome y from TSE
  2. Fixed outer split identical to above (`createDataPartition`, same seed)
  3. n_bins tuning — inner CV loop on training samples only:
     - `createFolds(y_train, k = kfold)` to create inner fold indices
     - For each `nb` in `n_bins_grid`:
       - For each inner fold:
         - `fit_bin_breaks(train_mat[, -val_idx], nb)` ← fit on inner-train
         - `apply_bin_breaks()` on inner-train and inner-val
         - `preProcess()` on inner-train binned matrix, apply to inner-val
         - `caret::train()` on inner-train; predict inner-val; record AUC
       - Compute mean AUC across folds for this `nb`
     - `best_n_bins = n_bins_grid[which.max(mean_aucs)]`
  4. Final refit:
     - `fit_bin_breaks(train_mat, best_n_bins)` ← full training set
     - `apply_bin_breaks()` on training mat and test mat
     - `preProcess()` on training binned mat, apply to test
     - `caret::train()` on training set; evaluate on test set
  5. Return the same tibble structure as `run_caret_on_assay()`, plus `tuning` trace
     (mean_auc per n_bins candidate) and `best_n_bins`

- [ ] Update `evaluate_transformation()` to call the new caret functions and
  preserve the existing output column contract (`method`, `transformation`,
  `selected_n_bins`, `Accuracy`, `Kappa`, plus `AUC` if not already present)

### Phase 3 — Pipeline cleanup (_targets.R)

- [ ] Update `tar_option_set(packages = ...)`: remove `mikropml`, add `caret`, `pROC`
- [ ] Confirm `model_grid` branching and `map(model_grid)` pattern still work
  (no structural changes needed — only `evaluate_transformation()` internals change)
- [ ] Invalidate the targets cache (`targets::tar_invalidate(everything())`) and do
  a full pipeline rerun to regenerate results

### Phase 4 — Report update (analysis/report.qmd)

- [ ] Add AUC as a reported metric alongside Accuracy and Kappa if not already present
- [ ] Re-render report after pipeline rerun

---

## Notes

- TSS, rCLR, and P/A are per-sample transforms (no cross-sample statistics), so
  pre-computing them in `add_all_transformations()` is fine and requires no change.
- The `fit_bin_breaks` / `apply_bin_breaks` split mirrors the caret `preProcess` /
  `predict.preProcess` pattern — fit an object on training data, apply it elsewhere.
- Keep the same outer-split seed (42) and 80/20 ratio so results are comparable
  to the mikropml baseline during the transition.
- caret's `"ranger"` method tunes `mtry`, `splitrule`, and `min.node.size` in the
  inner CV automatically; no changes needed to the classifier tuning grid.
