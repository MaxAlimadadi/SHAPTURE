# SHAPTURE pipeline

## Overview
This pipeline is designed to identify features that distinguish two biological conditions within each analysis subset, or one target subset versus the pooled background, then rank interpretable features using SHAP importance, Cross-patient support, model reproducibility, expression signal, and value-SHAP relationship. It supports both **single-cell** and **metacell** analysis units and can use either **XGBoost** or **LightGBM** as the classifier.

At a high level, the workflow:

1. starts from a prepared `dataset_all` table,
2. optionally builds metacells,
3. analyzes each selected subset separately, or each one-vs-others target separately,
4. repeatedly trains balanced classifiers for a two-condition comparison,
5. computes SHAP values and model metrics across iterations,
6. aggregates the results into robust feature rankings,
7. generates plots and summary tables for interpretation.

---

## Package workflow
After installing the dependencies listed in the README, install or refresh the
current GitHub version with:

```r
remotes::install_github("MaxAlimadadi/SHAPTURE", ref = "main", force = TRUE)
```

Use the editable templates from the same `main` branch. The installer may ask
whether to update dependency packages.

Load it with `library(SHAPTURE)`. The editable data-preparation,
analysis, and visualization scripts load the installed package; its internal
helper functions do not need to be sourced separately.

## Expected input structure
The main script assumes that `dataset_all` already exists in memory before the modeling section starts. The data-preparation templates show two ways to create it.

### Required columns in `dataset_all`
- `condition`: the two-group outcome to compare, for example `low` vs `high`
- `sample_id`: sample, donor, or patient identifier
- optional `split_level_column`: subset / cluster / cell-type / stratification label used for subset-specific analysis
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
- returns the prepared-data split name in `dataset_res$split_level_column`

#### 2. From a Seurat object
`getting_data_from_Seurat.R`:
- loads a Seurat object
- joins assay layers and normalizes the selected assay
- optionally subsets selected condition levels
- finds variable features
- extracts normalized expression values from the `data` layer
- creates `dataset_all`
- returns the prepared-data split name in `dataset_res$split_level_column`

For both adapters, `condition_col`, `sample_col`, and `split_level_column` refer to
columns in the source object's metadata and may have arbitrary names. SHAPTURE
standardizes the prepared outcome and sample roles to `condition` and
`sample_id`. If the source split name conflicts with either standardized name,
the adapter stores it under a deterministic internal alias. The analysis
config must therefore use `dataset_res$split_level_column`, as shown in the template
scripts. Nonconflicting split names are preserved exactly.

---

## Main analysis settings
The main script defines a set of user-controlled parameters.

### Comparison setup
- `split_level_column`: metadata column used for subset-specific analysis, or the target-defining column in `one_vs_others`
- `split_levels_to_analyze`: optional task levels to analyze; in `one_vs_others`, these are the target levels that get result folders
- `comparison_mode`: `condition_vs_condition` or `one_vs_others`
- `what_0`, `what_1`: the two condition labels to compare in `condition_vs_condition`
- `what_genes`: feature set used by the model
- `feature_min_pct`: task-specific detection filter from original single-cell task data; default is `0.25`, and `0` disables it
- `feature_min_log2FC`: task-specific expression log2FC filter from original single-cell task data; default is `0.25`, and `0` disables it
- `shap_expr_cor_cutoff`: minimum absolute value-to-SHAP correlation required for `Model_expression_pattern`; default is `0.25`
- `default_prior_scores_weights`: prior weights for the five score components
- `if_adaptive_weights`: whether to adapt those weights to component spread in each task; default is `FALSE`
- `score_scaling_method`: default is `"median_anchor"`; legacy `"winsorized"` and `"minmax"` values remain accepted for script compatibility, while evidence scoring uses the fixed median-anchor/native-bounded strategy
- `score_scaling_probs`: lower/upper quantiles used only for robust adaptive-weight spread estimates; default is `c(0.05, 0.95)`

### Repeated modeling
- `train_ratio`, `test_ratio`: training and testing fractions
- `iteration_no`: number of repeated model fits; default is `50`
- `gene_fraction_per_iteration`: fraction of filtered candidate features used
  in each coverage-aware model fit; default is `0.30`
- `cross_patient_bootstraps`: number of patient-block bootstrap replicates used
  for Cross-patient Support; default is `500L`
- `n_cores`: requested maximum workers for repeated model fitting and
  Cross-patient bootstrapping; the bootstrap stage can automatically use fewer
  workers when CPU or RAM is limited
- `analysis_seed`: reproducibility root for all stochastic stages; default is
  `123L`
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

### Parallel execution safeguards
Model-fitting PSOCK workers are task-specific and are stopped before
Cross-patient bootstrapping starts. For the bootstrap stage, SHAPTURE treats
`n_cores` as a maximum rather than a promise: it checks portable CPU limits,
installed and currently available RAM, the main R-process memory use, and a
task-specific worker-memory estimate. It reserves at least 25% of installed RAM
or 4 GiB, whichever is larger.

When at least two workers appear safe, SHAPTURE runs one real worker health
check and recalculates the safe worker count from the measured worker memory.
PSOCK workers are then recycled after at most five bootstrap replicates to
limit memory growth. The same controller is used on macOS, Windows, and Linux.
If cluster setup, the health check, a worker batch, or result validation fails,
SHAPTURE discards partial parallel output and replays every deterministic
bootstrap plan serially. No score definition, seed, resampling unit, or
bootstrap formula changes between the two execution paths.

Both `full` and `safe` modes use the same interpretation rule: if one cell/metacell is evaluated in multiple held-out test sets, its SHAP values and model-input feature values are averaged before feature-level interpretation summaries are calculated. When a gene was not included in a feature-bagged iteration, that missing cell-gene SHAP value remains missing and is not treated as zero.

### Reproducible random stages
`analysis_seed` is converted into stable task-specific and stage-specific child
seeds. Cell/metacell sampling, coverage-aware feature selection, model fitting,
hard downsampling, metacell construction, and Cross-patient bootstrapping do
not rely on the ambient R random state. Task seeds are derived from the task
identity rather than loop position, and iteration seeds are derived directly
from their iteration number. Therefore, adding a split level or increasing the
requested iterations does not alter the existing seed plan.

Each successful task saves the complete plan in
`docs/reproducibility_seeds.csv`. Exact reproduction still requires the same
input data and row order, package versions, and analysis parameters.

### Post-analysis plotting
- Analysis always saves result tables and `docs/visualization_inputs.rds`; figures are made afterward with `pipeline/run_shapture_visualization.R`.
- `pipeline/run_shapture_visualization.R` is one editable, loop-based recipe for one or several completed result folders, not a black-box rendering function. `top_feature_show`, `top_feature_label`, `top_feature_show_extended`, `width_per_feature`, and `height_per_feature` are defined once before the loop.
- Set `analysis_root` to a completed result folder or its parent. If that directory contains `docs/combined_results.csv` and `docs/visualization_inputs.rds`, only it is plotted. Otherwise, the script applies the same settings to every immediate child containing both files. Discovery is independent of comparison mode and folder name and does not search deeper levels.
- Condition-level waterfall plots use `top_feature_show` for their displayed row count and label every row with its average feature value where available. Any collapsed "other features" row retains its summary label. `top_feature_label` applies to selective labels in force and scatter plots only; it never hides or merges waterfall rows. The waterfall helper still accepts the older `top_feature_label` argument for compatibility but ignores it.
- Feature-dependent dimensions are calculated from the actual number displayed. Plot-specific base sizes and bounds remain in `calculate_shapture_figure_size()`, while every visible `pdf()` or `png()` call can still be changed to a fixed width or height directly. Density plots use the actual saved density-feature count to form a near-square grid with consistent panel space and additional horizontal room for the legend.
- The original beeswarm remains unchanged. The visualization recipe also creates an equal-count binned beeswarm and one dependency PDF per feature actually displayed in the beeswarm. Both new plots reuse the same feature-specific groups.
- `save_plot_objects`: `FALSE` by default inside the visualization script. Set it to `TRUE` there only when you want the reusable ggplot object list saved to `docs/plots.rds`; leaving it off makes large runs faster and uses less disk.

### Global visualization across split levels
- `pipeline/run_shapture_global_visualization.R` reads multiple completed result folders without rerunning model fitting or SHAP calculations.
- `top_feature_show` selects the top-ranked features separately from each folder. The union of those per-folder lists is used in the global figures, so a strong feature specific to one split level is not lost because it is absent elsewhere.
- The cross-level landscape displays final scores, model-expression pattern symbols, and average-linkage hierarchical clustering of both features and split levels.
- The direction-aware Jaccard similarity is the number of shared selected features with matching normalized direction divided by the union of selected feature names. Opposite directions do not count as a match.
- Condition-vs-condition folders must use the same condition labels. One-vs-others patterns are normalized to target versus others before comparison.
- Global PDF dimensions are calculated from the actual matrices after loading: the landscape uses its selected-feature and split-level counts, and the Jaccard matrix uses its split-level count on both axes. Visible `width` and `height` values default to `NULL`; enter an exact size in inches to override either dimension. Colors, borders, labels, and rasterization behavior remain editable in the script.

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
- optional subset label (`split_level_column`)
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
Filter the source object or `dataset_all` before running SHAPTURE to restrict
which levels are retained. `split_level_column` names the grouping column;
`split_levels_to_analyze` chooses tasks from the supplied data, with `NULL`
analyzing all available levels. It does not filter input rows.

Each analysis task is processed separately:

```r
for (temp_split in split_levels_to_analyze) {
  ...
}
```

This means SHAPTURE is not training one global model across all cell types. Instead, it asks a subset-specific question:

**Within this analyzed subset or target comparison, which features best distinguish condition 0 from condition 1?**

In `one_vs_others`, all non-target levels present in the supplied data form the
pooled “others” class. For example, an input containing levels `0`, `1` and `2`
with `split_levels_to_analyze = c(0, 1)` runs `0` vs `1+2` and `1` vs `0+2`.
Level `2` remains background without receiving its own result folder. To remove
it from the background, filter it out of the input before calling the pipeline.

Before model fitting, `feature_min_pct` and `feature_min_log2FC` reduce the feature set inside each task. Both filters are calculated from the original single-cell task data, even when the model itself is run on metacells. A feature is kept when it is detected as `> 0` in more than `feature_min_pct` of cells in at least one task group and has absolute expression log2FC at least `feature_min_log2FC`. The feature-level summary is written to `docs/feature_prefilter_summary.csv`, and each successful task writes a compact task-level `Integrated_adaptive_scoring_results/run_summary.csv`.

---

## Step 4. Balance the data across samples and conditions
A central concern in the pipeline is unequal representation across patients or samples. The workflow provides three strategies.

### A. `hard_median_downsample`
Each condition is downsampled separately before modeling.

- in `condition_vs_condition`, each sample is capped at the median sample size within its condition
- in `one_vs_others`, the target class is capped by `sample_id`, while the pooled non-target class is capped by `sample_id × split_level_column`

### B. `weighted_sampling`
All rows remain available, but train and test sampling probabilities are adjusted within each condition so overrepresented groups are sampled less often.

- in `condition_vs_condition`, the balancing group is `sample_id`
- in `one_vs_others`, the target class uses `sample_id`, while the pooled non-target class uses `sample_id × split_level_column`

### C. `model_weights`
All rows remain in the dataset, and the training rows receive model weights that are inversely related to group-level representation.

- in `condition_vs_condition`, the balancing group is `sample_id`
- in `one_vs_others`, the target class uses `sample_id`, while the pooled non-target class uses `sample_id × split_level_column`

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

The unified visualization recipe calculates dimensions automatically. Edit `width_per_feature` and `height_per_feature` near the top to change how most feature-dependent figures grow. Density plots instead calculate `ceiling(sqrt(n_features))` columns from the actual saved features and derive the required rows and device dimensions automatically, keeping facet dimensions approximately consistent. To use an exact size for one figure, replace `calculated_size$width` or `calculated_size$height` with a number in that figure's visible `pdf()`/`png()` call. Dot-label sizes remain visible beside the corresponding plot call. The comparison name is kept in the output folder, so integrated figure filenames are intentionally shorter, for example `ROC_plot.pdf`, `feature_importance_top.N.pdf`, `beeswarm_plot_top.N.features.png`, `binned_beeswarm_plot_top.N.features.png`, `score_decomposition_top.N.pdf`, `evidence_heatmap_top.N.pdf`, `shap_vs_expression.pdf`, and `valueSHAP_by_feature.pdf`. Here, `N` is the number actually shown. Dependency plots are written separately to `Integrated_adaptive_scoring_results/shap_value_dependency_plots/`.

---

## Step 5. Build train/test sets for each iteration
The package function `prepare_train_test_data(...)`, defined in
`R/shapture_processing.R`, handles the split.

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

After feature filtering, SHAPTURE also builds a coverage-aware gene-subset plan.
Each iteration uses the fraction specified by `gene_fraction_per_iteration`
(default `0.20`), and genes with the fewest previous inclusions are prioritized
for the next iteration. The gene coverage summary is saved as
`docs/feature_bagging_summary.csv`.

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
- **Model reproducibility** is a 0-1 resampling-consistency score, not
  independent biological replication, because different iterations reuse
  overlapping underlying data. It equally combines recurrence-adjusted
  raw mean absolute SHAP-importance retention, signed delta-SHAP retention,
  and Value-SHAP direction-and-threshold retention.

The three Model reproducibility measurements are calculated as follows:

- **SHAP-importance retention** is the fraction of valid feature-containing
  iterations with nonzero importance multiplied by the 10th-percentile-to-median
  ratio among positive raw mean absolute SHAP importances. Intentionally
  excluded feature iterations are not included in the denominator.
- **Delta-SHAP retention** is the fraction of valid iterations retaining the
  aggregate delta-SHAP direction multiplied by the 10th-percentile-to-median
  ratio of relative delta-SHAP magnitude among direction-matching iterations.
- **Value-SHAP retention** is the fraction of valid feature-containing
  iterations retaining the aggregate correlation direction and meeting
  `shap_expr_cor_cutoff`.

The three bounded results are averaged with equal weight. The final Model
reproducibility score is already on a 0-1 scale and is therefore not rescaled
across features before scoring.

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
- `cross_patient_support`
- `model_reproducibility_scaled`
- `expression_score`
- `shap_expr_cor_abs_scaled`

These represent a multi-criteria feature ranking strategy rather than pure SHAP magnitude alone. SHAPTURE does not stretch naturally bounded evidence across the observed feature range. Instead, unbounded SHAP importance and absolute expression log2FC use `x / (x + k)`, where `k` is the positive median of that measurement across the fixed candidate-feature set for the task. Absolute detection-rate difference, within-condition expression consistency, absolute Value-SHAP correlation, Cross-patient support, and model reproducibility keep their natural 0-1 values. If `if_adaptive_weights = FALSE`, the configured prior weights are used directly. If `if_adaptive_weights = TRUE`, SHAPTURE estimates component spread after winsorization at `score_scaling_probs`, shrinks those data-driven weights toward the prior weights, and normalizes the final weights to sum to 1.

### Meaning of the components
- **importance**: mean absolute SHAP importance transformed as `I / (I + k_I)`, where `k_I` is the median positive importance among fixed candidate features
- **Cross-patient support**: a 0-1 patient-block bootstrap retention score that combines support for SHAP importance, Expression signal, and the Value-SHAP relationship
- **model reproducibility**: the equal-weight mean of recurrence-adjusted
  raw mean absolute SHAP-importance retention, signed delta-SHAP retention,
  and Value-SHAP direction-and-threshold retention across repeated model fits.
  It is already bounded from 0 to 1 and is not rescaled across features.
- **expression**: support from how strongly, broadly, and consistently within each condition the feature is expressed. It averages median-anchored absolute expression log2FC, natural absolute detection-rate difference, and natural within-condition expression consistency.
- **value_shap**: natural absolute Spearman correlation between feature value and SHAP contribution

Within-condition expression consistency is calculated from the original task-level cells for both single-cell and metacell runs. For each feature and condition, zeros remain unchanged. If at least 20 cells are positive, only positive values above their condition-specific 99th percentile are capped. The condition-specific pseudocount is `max(0.01 * median(positive values), 1e-6)`, and robust relative variability is calculated as `SD / (abs(mean) + pseudocount)`. SHAPTURE converts this to `1 / (1 + CV)` and averages the two condition scores equally. A condition in which every cell is zero receives a consistency score of 1 for homogeneous absence; a feature absent in both conditions is removed. `within_condition_consistency_raw` stores the bounded score. The compatibility column `within_condition_consistency_scaled` now contains the same natural value rather than an across-feature rescaling.

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
- fit fixed task-level median anchors before direction and Model-expression pattern filtering,
- preserve naturally bounded evidence without across-feature rescaling,
- optionally adapt weights to component spread when `if_adaptive_weights = TRUE`,
- combine multiple evidence streams into a final score.

This makes the ranking more robust and less dependent on one noisy metric.

### Cross-patient support
SHAPTURE resamples patients as complete blocks while keeping the final feature
list fixed. `cross_patient_bootstraps` controls the number of replicates and
defaults to 500. For condition-vs-condition tasks, patients are resampled
separately within each condition and the original number of patient blocks is
drawn in each arm. For one-vs-others tasks, patients are resampled globally so
all target and background rows from a selected patient remain together.

Each bootstrap recalculates the same three main evidence streams used in ranking:

- raw mean absolute SHAP importance
- Expression signal from oriented expression log2FC, oriented detection-rate difference, and within-condition expression consistency
- Value-SHAP correlation

SHAP-importance support uses a smoother retention calculation. For raw
bootstrap importance values \(I_b\), values above \(\epsilon=10^{-8}\) are
treated as positive. SHAPTURE calculates the positive-value frequency
\(f_+\), the ratio between the positive-value fifth percentile and median
\(R_+\), and then combines them as \(S_{importance}=f_+R_+\). This penalizes
importance that is intermittent or loses magnitude without forcing support to
zero only because some bootstraps have zero importance. If no positive values
are observed, support is zero.

Within every patient bootstrap, oriented absolute expression log2FC is
transformed with the same fixed full-task log2FC anchor used by the main score.
Oriented absolute detection difference and within-condition consistency keep
their natural 0-1 values. Their average forms the bootstrap Expression signal,
which uses the fifth-percentile-to-median retention ratio across all valid
bootstrap values. Value-SHAP support remains the
fraction of valid bootstraps that retain the full-data correlation direction
and meet `shap_expr_cor_cutoff`. The equal-weight mean of the three support
measurements is saved as `cross_patient_support` on its natural 0-1 scale; it
is not rescaled across features.

Cross-patient support is available with at least two patients in each condition for condition-vs-condition, or at least two patients in total for one-vs-others. If the minimum is not met, the value remains `NA`, its weight is set to zero, and the other four weights are renormalized. The state and diagnostics are written to:

- `docs/stats_interpretation.csv`
- `docs/cross_patient_support_method_used.csv`
- `docs/cross_patient_support_by_feature.csv`
- `docs/cross_patient_bootstrap_runtime.csv`
- `docs/model_reproducibility_by_feature.csv`
- `docs/model_reproducibility_by_iteration.csv`
- `docs/score_anchor_summary.csv`
- `docs/shapture_config_used.csv`

---

## Step 10. Produce final ranked marker tables
For each analyzed subset or one-vs-others target, the pipeline combines per-iteration information into a final ranked result table.

Each completed analysis task writes:
- ranked marker tables per analysis task
- performance summaries
- fitted values / SHAP summaries
- interpretation summaries such as `stats_interpretation.csv`
- Cross-patient support method and per-feature bootstrap diagnostics
- Cross-patient bootstrap backend, resource guard, timing, and fallback audit
- a full resolved parameter record in `shapture_config_used.csv`
- publication-style plots

The `feature_importance` plot displays the top features in final-score order. Feature names use the same saved `Model_expression_pattern` colors as the beeswarm: `group_colors[2]` for condition 1/target and `group_colors[1]` for condition 0/others. Its legend lists condition 1 first, with the full pattern labels on two lines. The title shows the top-feature count; the comparison name remains in the result folder name.

---

## Step 11. Generate interpretation plots after analysis
The package imports plotting tools commonly used in this script:
- `ggplot2`
- `ggrepel`
- `cowplot`
- `reshape2`
- `paletteer`

After the analysis finishes, set `analysis_root` in `pipeline/run_shapture_visualization.R` to one completed result folder or the immediate parent of several, then run the script. A completed input folder is plotted directly; otherwise, each completed immediate child is plotted using the same settings. The visualization script reads `docs/combined_results.csv` and `docs/visualization_inputs.rds`, then calls each plot function explicitly and writes or overwrites the figure files. This unified recipe is available on the development branch after v0.1.0; install the matching package version when using its templates.

Expected outputs include:
- SHAP feature-importance plots
- the original unit-level beeswarm and an equal-count binned beeswarm
- one SHAP-value dependency PDF per feature shown in the beeswarm
- density plots
- ROC plots
- score decomposition plots for top features and all retained features
- support matrix / evidence heatmap plots
- five-component heatmaps when Cross-patient support is available, or four-component heatmaps when it is unavailable
- value-SHAP by-feature plots
- SHAP effect versus expression log2FC plots

Beeswarm point colors use the preprocessed/scaled feature values that were passed to the model. This keeps `full` and `safe` memory modes consistent and makes the beeswarm reflect what the model saw.

The binned beeswarm and dependency plot are display summaries. For each
feature, finite unit-level SHAP and model-seen feature values are mapped to
comparable 0-1 ranges for grouping only. A deterministic recursive partition
then splits whichever local dimension has the larger robust 5th-to-95th
percentile spread. Group capacities are balanced so they differ by no more
than one unit. The target is 500 groups with at least 10 units per group; when
fewer than 500 valid units exist, the original units are displayed instead.
Both plots reuse these exact groups and show group means, while all analysis
statistics remain unbinned. Each dependency smooth fits
`SHAP ~ s(feature value)` to the original unbinned units, not to the hexagons.
The group audit is saved as `docs/equal_count_shap_group_summary.csv`.

The support matrix summarizes final score, model reproducibility, and Cross-patient support for the top-ranked retained features. If Cross-patient support is unavailable, its column is omitted rather than shown as zero or a neutral midpoint.

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
- optional `split_level_column`
- `what_0`
- `what_1`
- `sample_id` column inside `dataset_all`
- an installed SHAPTURE package and a valid output directory

---

## Suggested repository structure
The package repository is organized as:

```text
SHAPTURE/
├── R/                         # package functions
├── man/                       # generated function documentation
├── pipeline/
│   ├── run_shapture_analysis.R
│   ├── run_shapture_visualization.R
│   ├── run_shapture_global_visualization.R
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
