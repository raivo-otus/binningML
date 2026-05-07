# Migration: mikropml → tidymodels + fix data leakage

## Problem

Two sources of data leakage in the current pipeline:

### 1. Binning applied to the full dataset before train/test split

In `evaluate_binning()` ([R/modeling.R:74-78](R/modeling.R#L74)), `transform_binning(tse, n_bins)` is called on the entire TSE — all samples including the held-out test set. The quantile bin breakpoints are therefore computed from the full sample distribution. When `run_ml()` later splits the pre-binned data, the test set was already binned using its own quantile information. This inflates performance estimates for binning relative to the other transformations.

### 2. `mikropml::preprocess_data()` fitted on the full dataset

In `run_ml_on_assay()` ([R/modeling.R:42](R/modeling.R#L42)), `preprocess_data(df)` (near-zero-variance filtering, centering, scaling) is called on the full df before `run_ml()` does the train/test split. Feature selection and scaling statistics therefore incorporate test-set information for all four transformations, not just binning.

---

## Proposed Fix

Replace mikropml with tidymodels (rsample / recipes / parsnip / workflows /
tune / yardstick) and implement a proper two-level CV structure. Preprocessing
is encoded as a `recipe` so prep/bake happen automatically per resample —
preventing leakage by construction.

For binning, we add a custom recipe step `step_quantile_bin()` that wraps the
Phase 1 utilities: `prep()` calls `fit_bin_breaks()` on the resample's analysis
set; `bake()` calls `apply_bin_breaks()`. `n_bins` is exposed as a tunable
parameter so `tune::tune_grid()` handles inner CV automatically.

```
rsample::initial_split(strata = outcome, prop = 0.8)   ← fixed seed
  └── training()
        ├── rsample::vfold_cv(strata, v = kfold)       ← inner CV
        │     workflow = recipe + model spec
        │       recipe: step_quantile_bin(n_bins = tune())
        │               step_nzv() → step_normalize()
        │       spec:   rand_forest / logistic_reg
        │     tune::tune_grid(workflow, resamples, grid = n_bins_grid,
        │                     metrics = metric_set(roc_auc, accuracy, kap))
        │     ← prep/bake fitted per fold on analysis() only
        │
        └── tune::select_best(metric = "roc_auc")
              tune::finalize_workflow(best_params)
              tune::last_fit(split)                    ← refit on training(),
                                                        evaluate on testing()
              yardstick metrics: accuracy, kap, roc_auc
```

For the non-binning transformations (relabundance, rclr, pa) — already
pre-computed on the full TSE but sample-wise with no cross-sample statistics —
the recipe drops `step_quantile_bin()` and keeps only `step_nzv()` +
`step_normalize()`:

```
rsample::initial_split(strata, prop = 0.8)             ← same seed
  workflow = recipe(step_nzv → step_normalize) + spec
  tune::tune_grid() over classifier hyperparameters    ← rsample::vfold_cv
  tune::last_fit() on the outer split
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

### Phase 2 — Replace mikropml modeling with tidymodels (R/modeling.R)

- [ ] Remove `mikropml` from `tar_option_set(packages = ...)` in `_targets.R`; add
  `tidymodels` (or the constituent packages: `rsample`, `recipes`, `parsnip`,
  `workflows`, `tune`, `yardstick`, `dials`). Backends `ranger` and `glmnet`
  are already present.

- [ ] Add a custom recipe step `step_quantile_bin()` (R/recipe_steps.R):
  - Constructor stores `n_bins` (tunable via `dials::new_quant_param()`) and
    selected columns
  - `prep.step_quantile_bin()`: pull analysis-set predictors as a features ×
    samples matrix, call `fit_bin_breaks(mat, n_bins)`, store the breaks list
    on the step
  - `bake.step_quantile_bin()`: apply `apply_bin_breaks()` using the stored
    breaks; replace the columns in the baked tibble
  - Register a tunable parameter `n_bins()` so `tune::tune_grid()` recognises it
  - This wraps the Phase 1 utilities — no new binning math here.

- [ ] Add helper `model_spec(method)` returning a parsnip spec:
  - `"rf"` → `parsnip::rand_forest(mtry = tune(), min_n = tune(), trees = 500)`
    `|> set_engine("ranger") |> set_mode("classification")`
  - `"glmnet"` → `parsnip::logistic_reg(penalty = tune(), mixture = tune())`
    `|> set_engine("glmnet") |> set_mode("classification")`

- [ ] Add helper `metric_set_clf()` returning
  `yardstick::metric_set(roc_auc, accuracy, kap)` (replaces the caret
  AUC/confusionMatrix combo).

- [ ] Replace `run_ml_on_assay()` with
  `run_tidymodels_on_assay(tse, assay_name, method, label_col, kfold, seed)`:
  1. Build a tibble with outcome + features from the TSE assay.
  2. `set.seed(seed)`; `split <- rsample::initial_split(df, prop = 0.8, strata = outcome)`.
  3. Recipe: `recipe(outcome ~ ., data = training(split)) |> step_nzv(all_predictors()) |> step_normalize(all_predictors())`.
  4. `wf <- workflow() |> add_recipe(rec) |> add_model(model_spec(method))`.
  5. Inner CV: `folds <- rsample::vfold_cv(training(split), v = kfold, strata = outcome)`.
  6. `tuned <- tune::tune_grid(wf, resamples = folds, grid = ..., metrics = metric_set_clf())`.
  7. `best <- tune::select_best(tuned, metric = "roc_auc")`;
     `final_wf <- tune::finalize_workflow(wf, best)`.
  8. `final_fit <- tune::last_fit(final_wf, split = split, metrics = metric_set_clf())`.
  9. Return `tune::collect_metrics(final_fit) |> tidyr::pivot_wider(...)` shaped
     to match the existing column contract (`Accuracy`, `Kappa`, `AUC`).

- [ ] Replace `evaluate_binning()` with
  `evaluate_binning_tidymodels(tse, n_bins_grid, method, label_col, kfold, seed)`:
  1. Build a tibble with outcome + raw counts from the TSE.
  2. Same outer split (`initial_split`, same seed).
  3. Recipe prepends `step_quantile_bin(all_predictors(), n_bins = tune())`
     before `step_nzv()` and `step_normalize()` — every other step is identical
     to `run_tidymodels_on_assay()`.
  4. Inner CV: `vfold_cv(training(split), v = kfold, strata = outcome)`.
  5. Tuning grid combines `n_bins_grid` with the model's hyperparameters
     (e.g. via `dials::grid_regular()` or an explicit `tidyr::crossing()`).
     `tune_grid()` re-preps the recipe per fold automatically — no leakage.
  6. `select_best(metric = "roc_auc")` → `finalize_workflow()` → `last_fit()`.
  7. Return the same shaped row as `run_tidymodels_on_assay()`, plus a `tuning`
     trace (`tune::collect_metrics(tuned, summarize = TRUE)` filtered to
     `roc_auc`) and `best_n_bins` extracted from `best$n_bins`.

- [ ] Update `evaluate_transformation()` to call the new tidymodels functions
  and preserve the existing output column contract (`method`,
  `transformation`, `selected_n_bins`, `Accuracy`, `Kappa`, plus `AUC` if
  not already present).

### Phase 3 — Pipeline cleanup (_targets.R)

- [ ] Update `tar_option_set(packages = ...)`: remove `mikropml`, add
  `tidymodels` (or the explicit set: `rsample`, `recipes`, `parsnip`,
  `workflows`, `tune`, `yardstick`, `dials`)
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
- The `fit_bin_breaks` / `apply_bin_breaks` split mirrors the recipes
  `prep()` / `bake()` pattern — fit on training data, apply elsewhere — and
  is what `step_quantile_bin()` wraps internally.
- Keep the same outer-split seed (42) and 80/20 ratio so results are comparable
  to the mikropml baseline during the transition.
- For `rand_forest` with the `ranger` engine, parsnip exposes `mtry`, `min_n`,
  and `trees` as tuning parameters; `splitrule` is fixed to `"gini"` for
  classification and is not tuned.
- `tune::last_fit()` is the tidymodels primitive for "refit on the training
  split, evaluate once on the test split" — it takes the `rsample::initial_split`
  object directly, so the outer split is reused without re-derivation.
