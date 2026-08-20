# SHAPTURE pipeline

## Overview
This pipeline is designed to identify features that distinguish two biological conditions within each analysis subset, or one target subset versus the pooled background, then rank interpretable features using SHAP importance, patient-level support, model reproducibility, expression signal, and value-SHAP relationship. It supports both **single-cell** and **metacell** analysis units and can use either **XGBoost** or **LightGBM** as the classifier.

At a high level, the workflow:

1. starts from a prepared `dataset_all` table,
2. optionally builds metacells,
3. analyzes each selected subset separately, or each one-vs-others target separately,
4. repeatedly trains balanced classifiers for a two-condition comparison,
5. computes SHAP values and model metrics across iterations,
6. aggregates the results into robust feature rankings,
7. generates plots and summary tables for interpretation.

---

## Expected input structure
The main script assumes that `dataset_all` already exists in memory before the modeling section starts. The helper scripts show two ways to create it.

### Required columns in `dataset_all`
- `condition`: the two-group outcome to compare, for example `low` vs `high`
- `sample_id`: sample, donor, or patient identifier
- optional `split_by_col`: subset / cluster / cell-type / stratification label used for subset-specific analysis
- feature columns: genes, proteins, or other quantitative variables used as predictors

### Input preparation options
#### 1. From a SingleCellExperiment object
`getting_data_from_sce.R`:
- loads an SCE object
- selects a clustering resolution
- selects a condition column and sample ID column
- extracts marker features from `type_markers(sce)` and `state_markers(sce)`
- optionally removes unwanted markers
- creates `dataset_all`

#### 2. From a Seurat object
`getting_data_from_Seurat.R`:
- loads a Seurat object
- joins assay layers and normalizes the selected assay
- optionally subsets selected condition levels
- finds variable features
- extracts normalized expression values from the `data` layer
- creates `dataset_all`

---

## Main analysis settings
The main script defines a set of user-controlled parameters.

### Comparison setup
- `split_by_col`: metadata column used for subset-specific analysis, or the target-defining column in `one_vs_others`
- `subset_levels`: optional `split_by_col` levels to keep before task construction; `NULL` keeps all levels
- `split_levels`: optional task levels to analyze; in `one_vs_others`, these are the target levels that get result folders
- `comparison_mode`: `condition_vs_condition` or `one_vs_others`
- `what_0`, `what_1`: the two condition labels to compare in `condition_vs_condition`
- `what_genes`: feature set used by the model
- `feature_min_pct`: task-specific detection filter from original single-cell task data; default is `0.25`, and `0` disables it
- `feature_min_log2FC`: task-specific expression log2FC filter from original single-cell task data; default is `0.25`, and `0` disables it
- `shap_expr_cor_cutoff`: minimum absolute value-to-SHAP correlation required for `Model_expression_pattern`; default is `0.25`
- `default_prior_scores_weights`: prior weights for the five score components
- `if_adaptive_weights`: whether to adapt those weights to component spread in each task; default is `FALSE`
- `score_scaling_method`: `"winsorized"` or `"minmax"`; default is `"winsorized"`
- `score_scaling_probs`: lower/upper quantiles used by winsorized component scaling and adaptive-weight spread estimates; default is `c(0.05, 0.95)`

### Repeated modeling
- `train_ratio`, `test_ratio`: training and testing fractions
- `iteration_no`: number of repeated model fits
- `patient_split_strategy`: `shared_patients` (default) or `unseen_patients`
- `top_feature_label`: number of features to label in some plots
- `top_feature_show`: number of top features to display and save for post-analysis plotting

### Model choice
- `model_type`: either `xgboost` or `lightgbm`

### Analysis unit choice
- `analysis_unit_type`: either `single_cell` or `metacell`

### Memory mode
- `memory_mode`:
  - `full`
  - `safe`

`safe` mode writes iteration-level SHAP payloads to temporary disk chunks, rebuilds pooled SHAP summaries one chunk at a time, and deletes the temporary files automatically at the end. It also avoids sending full per-iteration model-fit objects back through parallel workers. It is designed to reduce peak RAM use without changing the final result tables or figures.

Both `full` and `safe` modes use the same interpretation rule: if one cell/metacell is evaluated in multiple held-out test sets, its SHAP values and model-input feature values are averaged before feature-level interpretation summaries are calculated. When a gene was not included in a feature-bagged iteration, that missing cell-gene SHAP value remains missing and is not treated as zero.

### Post-analysis plotting
- Analysis always saves result tables and `docs/visualization_inputs.rds`; figures are made afterward with `pipeline/run_shapture_visualization.R`.
- `pipeline/run_shapture_visualization.R` is an editable figure recipe, not a black-box rendering function. Each figure has a visible plot-function call and a visible `pdf()` or `png()` call where width, height, and dot-label size can be changed directly.
- `save_plot_objects`: `FALSE` by default inside the visualization script. Set it to `TRUE` there only when you want the reusable ggplot object list saved to `docs/plots.rds`; leaving it off makes large runs faster and uses less disk.

### Global visualization across split levels
- `pipeline/run_shapture_global_visualization.R` reads multiple completed result folders without rerunning model fitting or SHAP calculations.
- `top_feature_show` selects the top-ranked features separately from each folder. The union of those per-folder lists is used in the global figures, so a strong feature specific to one split level is not lost because it is absent elsewhere.
- The cross-level landscape displays final scores, model-expression pattern symbols, and average-linkage hierarchical clustering of both features and split levels.
- The direction-aware Jaccard similarity is the number of shared selected features with matching normalized direction divided by the union of selected feature names. Opposite directions do not count as a match.
- Condition-vs-condition folders must use the same condition labels. One-vs-others patterns are normalized to target versus others before comparison.
- Width, height, colors, borders, labels, and rasterization behavior remain visible and editable in the global visualization script.

### Patient balancing choice
- `patient_balance_method`:
  - `hard_median_downsample`
  - `weighted_sampling`
  - `model_weights`

### Patient train/test strategy
- `shared_patients`: preserves the original cell/metacell-level split; a test
  unit is never a training unit in the same iteration, but the same patient may
  contribute different units to both sets
- `unseen_patients`: assigns eligible patients to globally disjoint train and
  test groups using rounded 60%/40% counts separately for each condition

The unseen-patient strategy requires at least five eligible patients per
condition. Eligibility requires at least 20 cells per patient-condition in
single-cell mode or 5 metacells per patient-condition in metacell mode. Each
fold must retain at least 50 training cells or 20 training metacells per
condition. If any requirement cannot be met, the complete task falls back to
`shared_patients`, emits a warning, and saves the reason.

### Metacell settings
Used only when `analysis_unit_type == "metacell"`:
- minimum cells needed to attempt metacell construction
- minimum and maximum cells per metacell
- number of variable features
- number of PCs
- aggregation function such as `mean` or `sum`
- parallelization settings

### Temporary environment overrides
The script also supports environment-variable overrides for quick test runs:
- `TEMP_RUN_TAG`
- `TEMP_ITERATION_NO`
- `TEMP_ONLY_SPLIT`
- `TEMP_ONLY_CLUSTER` (legacy alias)
- `TEMP_ANALYSIS_UNIT_TYPE`
- `TEMP_PATIENT_BALANCE_METHOD`

---

## Step-by-step pipeline logic

## Step 1. Start from a prepared analysis table
The workflow begins with `dataset_all`, a table where rows are cells or later metacells, and columns contain:
- outcome label (`condition`)
- sample identity (`sample_id`)
- optional subset label (`split_by_col`)
- numeric predictors (`what_genes`)

The script records the number of observations per condition and subset before entering the main loop.

---

## Step 2. Optionally convert cells into metacells
If `analysis_unit_type` is set to `metacell`, the script calls `build_grouped_metacells(...)`.

Based on the arguments shown, metacells are built:
- within each sample,
- within each sample, subset, and condition,
- while preserving condition labels,
- using variable features and PCA,
- then aggregating cells into metacell-level expression profiles.

Outputs:
- `dataset_all_analysis`: the table that will actually be modeled
- `metacell_summary`: summary information about metacell construction

If `analysis_unit_type == "single_cell"`, then `dataset_all_analysis <- dataset_all`.

If `memory_mode == "safe"`, this step does not change the modeling inputs. It only changes how large iteration-level SHAP outputs are stored and reconstructed later.

---

## Step 3. Loop over analysis tasks
Before tasks are built, `subset_levels` can optionally reduce `dataset_all` to only selected levels of `split_by_col`. `split_levels` then chooses which tasks are run from the retained data.

Each analysis task is processed separately:

```r
for (temp_split in split_levels) {
  ...
}
```

This means SHAPTURE is not training one global model across all cell types. Instead, it asks a subset-specific question:

**Within this analyzed subset or target comparison, which features best distinguish condition 0 from condition 1?**

In `one_vs_others`, `subset_levels` controls the universe used to construct the pooled “others” class. For example, `subset_levels = c(0, 1, 2)` and `split_levels = c(0, 1)` runs `0` vs `1+2` and `1` vs `0+2`; level `2` is available as background but does not get its own folder unless it is also included in `split_levels`.

Before model fitting, `feature_min_pct` and `feature_min_log2FC` reduce the feature set inside each task. Both filters are calculated from the original single-cell task data, even when the model itself is run on metacells. A feature is kept when it is detected as `> 0` in more than `feature_min_pct` of cells in at least one task group and has absolute expression log2FC at least `feature_min_log2FC`. The feature-level summary is written to `docs/feature_prefilter_summary.csv`, and each successful task writes a compact task-level `Integrated_adaptive_scoring_results/run_summary.csv`.

---

## Step 4. Balance the data across samples and conditions
A central concern in the pipeline is unequal representation across patients or samples. The workflow provides three strategies.

### A. `hard_median_downsample`
Each condition is downsampled separately before modeling.

- in `condition_vs_condition`, each sample is capped at the median sample size within its condition
- in `one_vs_others`, the target class is capped by `sample_id`, while the pooled non-target class is capped by `sample_id × split_by_col`

### B. `weighted_sampling`
All rows remain available, but train and test sampling probabilities are adjusted within each condition so overrepresented groups are sampled less often.

- in `condition_vs_condition`, the balancing group is `sample_id`
- in `one_vs_others`, the target class uses `sample_id`, while the pooled non-target class uses `sample_id × split_by_col`

### C. `model_weights`
All rows remain in the dataset, and the training rows receive model weights that are inversely related to group-level representation.

- in `condition_vs_condition`, the balancing group is `sample_id`
- in `one_vs_others`, the target class uses `sample_id`, while the pooled non-target class uses `sample_id × split_by_col`

This step is meant to reduce the chance that one patient with many cells dominates the classifier, and in `one_vs_others` mode it also reduces domination by a very large non-target split level.

---

## Step 4B. Configure plot colors and post-analysis visualization
The plotting layer now supports two user-facing color controls:

- `group_colors`
  - an unnamed two-color vector
  - first color = left / `what_0` / `others`
  - second color = right / `what_1` / target
- `continuous_colors`
  - any valid R color vector with at least two colors
  - used for SHAP, support, expression, and other continuous gradients

Examples that work:
- `c("navy", "firebrick")`
- `c("navy", "white", "firebrick")`
- `viridisLite::viridis(100)`
- `RColorBrewer::brewer.pal(11, "RdBu")`

During analysis, each completed result folder saves `docs/visualization_inputs.rds`. This is the reusable plot-input bundle used by `pipeline/run_shapture_visualization.R`; it lets you change top-feature counts, colors, sizes, and label sizes without rerunning model fitting or SHAP calculation. For unit-level plots such as beeswarm and density plots, decreasing top-N works from saved inputs; increasing beyond the originally saved top features requires rerunning analysis with a larger `top_feature_show`.

When `save_plot_objects = TRUE` in the visualization step, the result folder also saves `docs/plots.rds`, which contains the generated plot objects. This makes it possible to reopen a plot later, print it at a different size, or modify themes, fonts, labels, and colors without rerunning the analysis. The default is `FALSE` because saving large plot-object lists can slow down large runs.

Saved figure dimensions and dot-label sizes are controlled directly inside `pipeline/run_shapture_visualization.R`. Edit variables such as `density_plot_width`, `density_plot_height`, or `shap_vs_expression_dot_label_size` immediately before the corresponding `pdf()`/`png()` call. The comparison name is kept in the output folder, so integrated figure filenames are intentionally shorter, for example `ROC_plot.pdf`, `feature_importance_top.N.pdf`, `score_decomposition_top.N.pdf`, `evidence_heatmap_top.N.pdf`, `shap_vs_expression.pdf`, and `valueSHAP_by_feature.pdf`.

---

## Step 5. Build train/test sets for each iteration
The function `prepare_train_test_data(...)` in `helper_function.R` handles the split.

What it does:
- uses a coverage-aware test split plan built for the whole task
- separates metadata from feature columns
- enforces exactly two condition levels
- samples a balanced number of rows per condition for training
- samples a balanced number per condition for testing, prioritizing cells/metacells with the fewest previous test appearances
- optionally uses weighted sampling probabilities
- optionally scales the data

The coverage summary is saved as `docs/sampling_coverage_summary.csv`. When the total number of requested held-out rows is large enough, the goal is for every available cell/metacell in each condition to appear in the test set at least once across iterations.

When `patient_split_strategy = "unseen_patients"`, patient separation happens
before row selection. Sixty percent of eligible patients are assigned to
training and forty percent to testing after rounding; all eligible patients are
used in one of the two groups. Patient IDs are separated globally, including in
one-vs-others tasks, and patients with fewer previous test appearances are
prioritized in later iterations. Test rows are distributed across held-out
patients so a patient with many more cells does not dominate simply because of
its size.

The requested and applied strategy is saved to
`docs/patient_split_summary.csv`. Patient eligibility and rotation are saved to
`docs/patient_sampling_coverage_summary.csv`.

After feature filtering, SHAPTURE also builds a coverage-aware gene-subset plan. Each iteration uses 20% of the candidate genes, and genes with the fewest previous inclusions are prioritized for the next iteration. The gene coverage summary is saved as `docs/feature_bagging_summary.csv`.

For unseen-patient tasks, feature selection is coordinated with the current
held-out patients. Genes with less prior SHAP coverage in those patients are
prioritized while overall gene inclusion remains balanced. The achieved
patient-condition coverage for every gene is saved to
`docs/patient_feature_coverage_summary.csv`. Missing patient-gene observations
remain missing and are never replaced with zero. Its denominator includes all
eligible patient-condition units, including eligible patients that were never
held out because too few iterations were requested.

### Preprocessing options
- `none`: pass raw values through
- `scale`: perform train-only preprocessing:
  1. z-score using training mean and SD
  2. rescale to 0–1 using training min and max
  3. apply the same transform to test data
  4. clip out-of-range test values to 0–1

This avoids leakage from the test set into the training transform.

---

## Step 6. Train the classifier
For each iteration and cluster, the pipeline trains a binary classifier using either:
- **XGBoost**, or
- **LightGBM**

The exact fitting code is in the main script beyond the excerpt you shared, but the imported packages and structure show that the model is repeatedly fit on balanced train/test splits.

The target is the binary `condition` variable.

---

## Step 7. Evaluate predictive performance
The script imports `pROC` and computes classification metrics across repeated fits. In the broader pipeline, this typically includes:
- AUC
- ROC-related summaries
- confusion-matrix-derived metrics such as accuracy, sensitivity, specificity, precision, recall, and F1

Because the model is repeated many times, the pipeline can later summarize both average performance and stability.

AUC and the other predictive metrics are calculated from the same test units
used for SHAP. Under `shared_patients`, these are unseen cells or metacells from
patients who may also contribute training units. Under `unseen_patients`, they
come only from patients completely absent from training.

Important interpretation note:
- `stability` is displayed as **Model reproducibility** and should be read as resampling consistency, not as independent biological replication, because different iterations reuse overlapping underlying data.

---

## Step 8. Compute SHAP values
After model fitting, the pipeline uses `shapviz` to calculate SHAP values for each feature.

These SHAP values are then used to summarize:
- how important a feature is for classification,
- the direction and magnitude of feature contribution,
- how stable that contribution is across iterations,
- how feature values relate to SHAP values.

This is the interpretability core of the pipeline.

If the same cell/metacell is held out in more than one iteration, SHAPTURE averages that unit's SHAP values and scaled model-input feature values before calculating SHAP-expression correlations, patient-level SHAP summaries, beeswarm inputs, condition-level waterfall summaries, and SHAP importance used for ranking. This avoids treating repeated test appearances of the same biological unit as independent observations in interpretation outputs.

Because each model sees only a coverage-aware subset of genes, SHAP values are averaged only over observed cell-gene pairs. If a gene was absent from a model in an iteration, that missing SHAP value is ignored rather than filled with zero.

---

## Step 9. Aggregate feature evidence across iterations
The pipeline does not rely on a single model run. Instead, it combines evidence across repeated iterations.

Model performance and model reproducibility remain iteration-aware. Feature interpretation summaries are unit- and gene-aware first, so a cell/metacell contributes once to biological interpretation for each evaluated gene even if it was evaluated in several held-out iterations.

The final ranking incorporates multiple scaled components:
- `importance_scaled`
- `patient_level_support_scaled`
- `model_reproducibility_scaled`
- `expression_score`
- `shap_expr_cor_abs_scaled`

These represent a multi-criteria feature ranking strategy rather than pure SHAP magnitude alone. By default, each component is winsorized at `score_scaling_probs` and then rescaled to 0-1 before scoring. If `if_adaptive_weights = FALSE`, the configured prior weights are used directly. If `if_adaptive_weights = TRUE`, SHAPTURE estimates component spread after the same winsorization step, shrinks those data-driven weights toward the prior weights, and normalizes the final weights to sum to 1.

### Meaning of the components
- **importance**: average feature importance or SHAP-based importance
- **patient-level support**: a patient-count-aware component; it uses descriptive `bootstrap_support` in low-patient settings and permutation-adjusted p-value evidence when patient count is large enough
- **model reproducibility**: reproducibility across repeated subsampling/model fits
- **expression**: support from how strongly, broadly, and consistently within each condition the feature is expressed. It combines scaled absolute expression log2FC, scaled absolute detection-rate difference, and scaled inverse average within-condition expression variance.
- **value_shap**: relationship between feature value and SHAP contribution

### Model Expression Pattern
`Model_expression_pattern` is saved in `combined_results.csv` and used as an interpretation label in plots.

It is driven by:
- the sign of `expression_log2FC`
- the sign of `shap_expr_cor`

So it is best interpreted as:
- **Which condition do higher feature values push the model toward, when the expression contrast and value-to-SHAP relationship are both strong and internally consistent?**

The current cutoffs are:
- `abs(expression_log2FC) >= feature_min_log2FC`
- `abs(shap_expr_cor) >= shap_expr_cor_cutoff`

Its values are:
- `Higher feature values push toward <condition_1 label>`
- `Higher feature values push toward <condition_0 label>`
- `Insufficient directional evidence`

`Insufficient directional evidence` is assigned when either value is weak (`abs(expression_log2FC) < feature_min_log2FC` or `abs(shap_expr_cor) < shap_expr_cor_cutoff`) or when `expression_log2FC` and `shap_expr_cor` have different signs.

After SHAP values are computed, features with `Insufficient directional evidence` are removed before final scoring and ranking. This means final ranks are assigned only among retained, directionally interpretable features.

The scoring function can:
- check that all required summary columns exist,
- apply prior weights,
- winsorize extreme component values before scaling when `score_scaling_method = "winsorized"`,
- optionally adapt weights to component spread when `if_adaptive_weights = TRUE`,
- combine multiple evidence streams into a final score.

This makes the ranking more robust and less dependent on one noisy metric.

### Patient-count-aware interpretation
The current pipeline tracks the smallest patient count across the two conditions:

- if `min_patients < 3`, the run is treated as descriptive only
- in that setting, `bootstrap_support` is used as the support component
- if `min_patients >= 3`, permutation-adjusted p-values can be used as the patient-level support component
- `pval_min_achievable` may be `NA`, which means an exact patient-label permutation floor is not informative for the current patient counts

This state is written to:

- `docs/stats_interpretation.csv`
- `docs/patient_level_support_method_used.csv`
- `docs/shapture_config_used.csv`

---

## Step 10. Produce final ranked marker tables
For each analyzed subset or one-vs-others target, the pipeline combines per-iteration information into a final ranked result table.

Although the exact write-out code is not visible in the snippet, the script structure strongly indicates outputs such as:
- ranked marker tables per analysis task
- performance summaries
- fitted values / SHAP summaries
- interpretation summaries such as `stats_interpretation.csv`
- patient-level support mode records such as `patient_level_support_method_used.csv`
- a full resolved parameter record in `shapture_config_used.csv`
- publication-style plots

The `feature_importance` plot title intentionally displays the positive/right level first, as `condition_1/target vs condition_0/others`, so the title matches the model direction being prioritized.

---

## Step 11. Generate interpretation plots after analysis
The package imports plotting tools commonly used in this script:
- `ggplot2`
- `ggrepel`
- `cowplot`
- `reshape2`
- `paletteer`

After the analysis finishes, set `result_dir` in `pipeline/run_shapture_visualization.R` to one completed result folder and run the script. The visualization script reads `docs/combined_results.csv` and `docs/visualization_inputs.rds`, then calls each plot function explicitly and writes or overwrites the figure files.

Expected outputs include:
- SHAP feature-importance plots
- beeswarm-like SHAP plots
- density plots
- ROC plots
- score decomposition plots for top features and all retained features
- support matrix / evidence heatmap plots
- five-component heatmaps for top features and all retained features
- value-SHAP by-feature plots
- SHAP effect versus expression log2FC plots

Beeswarm point colors use the preprocessed/scaled feature values that were passed to the model. This keeps `full` and `safe` memory modes consistent and makes the beeswarm reflect what the model saw.

The support matrix summarizes final score, model reproducibility, and patient-level support for the top-ranked retained features. Confidence-interval support is no longer part of the final score or kept figures.

---

## Conceptual output of the pipeline
For each analyzed subset or one-vs-others target, SHAPTURE aims to answer:

1. **Can this cluster distinguish the two conditions?**
2. **Which features drive that distinction?**
3. **Are those features reproducible across repeated runs?**
4. **Are they expressed strongly enough to be biologically useful?**
5. **Do feature values and SHAP contributions move together in a meaningful way?**

So the output is not just a classifier. It is a **subset-specific, stability-aware, SHAP-based marker discovery workflow**.

---

## Practical interpretation
This pipeline is especially useful when:
- you want subset-level condition markers,
- sample sizes differ across donors,
- you want to reduce patient dominance,
- you care about robustness across repeated subsampling,
- you want interpretable machine-learning markers rather than only differential expression.

It is well suited to:
- scRNA-seq
- CITE-seq
- CyTOF
- other single-cell or metacell structured datasets

as long as the data can be represented as a feature matrix with sample IDs, condition labels, and any subset labels you want to analyze.

---

## Minimal required objects before running
At minimum, the user should define:
- `dataset_all`
- `what_genes`
- optional `split_by_col`
- `what_0`
- `what_1`
- `sample_id` column inside `dataset_all`
- a valid working directory and helper script path

---

## Suggested repository structure
The package repository is organized as:

```text
SHAPTURE/
├── R/                         # package functions
├── man/                       # generated function documentation
├── pipeline/
│   ├── helper_function.R
│   ├── run_shapture_analysis.R
│   ├── run_shapture_visualization.R
│   ├── run_shapture_visualization_one.vs.others.R
│   ├── getting_data_from_sce.R
│   ├── getting_data_from_Seurat.R
│   └── shapture_config_full_reference.R
├── tests/testthat/            # automated package tests
├── DESCRIPTION
├── NAMESPACE
├── README.md
└── pipeline.md
```

---

## One-sentence summary
SHAPTURE is a repeated, subset-specific, SHAP-based machine-learning pipeline that identifies robust condition-associated markers while accounting for sample imbalance and optionally aggregating cells into metacells.
