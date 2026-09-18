# SHAPTURE

SHAPTURE is an R package for sample-aware, SHAP-based feature discovery in
single-cell and metacell data. It is designed for studies where the main
question is not only whether a model can separate two groups, but which
features drive that separation in a way that is reproducible, interpretable,
and connected to expression biology.

The package supports two common analysis patterns:

- condition-vs-condition comparisons, either globally or within selected
  subsets
- one-vs-others comparisons, where each selected level is compared with the
  remaining levels

SHAPTURE can start from a prepared table, a Seurat object, or a
SingleCellExperiment / CATALYST object.

## Why SHAPTURE

Single-cell datasets often contain uneven numbers of cells per sample, repeated
measurements from the same patient or donor, and biological effects that differ
across subsets. A simple differential-expression test can miss model-level
patterns, while a single machine-learning model can overstate unstable signals.

SHAPTURE tries to sit between those worlds. It repeatedly trains binary
classifiers, computes SHAP values on held-out cells or metacells, summarizes
the model behavior at the feature level, and ranks features using five evidence
components:

- SHAP importance
- expression signal
- value-SHAP relationship
- Cross-patient support
- model reproducibility

The expression signal combines absolute expression log2FC, absolute
detection-rate difference, and within-condition expression consistency. The
consistency term is calculated separately in each condition from the original
cell values using a zero-preserving robust coefficient of variation, converted
to a 0-1 score, and averaged with equal condition weight. This calculation uses
the original cells even when metacells are used for model fitting.

Cross-patient support asks whether the three main evidence streams persist when
patients are resampled as complete blocks. By default, SHAPTURE runs 500
patient-block bootstraps, recalculates SHAP importance, Expression signal, and
the Value-SHAP relationship for the fixed retained feature set, and combines
their 0-1 retention measures with equal weight. SHAP-importance support is the
frequency of positive bootstrap importance multiplied by the fifth-percentile-
to-median ratio among positive values. Expression signal keeps its
fifth-percentile-to-median retention measure, while Value-SHAP support remains
the frequency of bootstraps that retain the required direction and correlation
threshold. The number of replicates is controlled by
`cross_patient_bootstraps`. If there are too few patients, this component is
left unavailable and the remaining four score weights are renormalized.

Cross-patient bootstrap calculations can run in parallel without changing the
patient resampling plans or score calculations. `n_cores` is the requested
maximum for both repeated model fitting and this bootstrap stage. Before
bootstrapping, SHAPTURE checks the CPU and currently available RAM, keeps at
least 25% of installed memory or 4 GiB in reserve, and may lower the worker
count or use serial execution. A one-worker health check measures real memory
use, workers are recycled after at most five bootstrap replicates, and any
parallel setup, worker, or output-validation failure triggers a complete serial
replay. This avoids mixing partial parallel and serial results.

Model reproducibility is calculated across repeated model fits from three
equally weighted measurements: recurrence-adjusted raw mean absolute
SHAP-importance retention, signed delta-SHAP retention, and Value-SHAP
direction-and-threshold retention. Its natural 0-1 score enters the final
ranking without across-feature rescaling.

SHAPTURE keeps naturally bounded evidence on its original 0-1 scale. Absolute
detection-rate difference, within-condition expression consistency, absolute
Value-SHAP correlation, Cross-patient support, and model reproducibility are
therefore not stretched across features. The two unbounded measurements use a
fixed median-anchor transformation, `x / (x + k)`: `k` is the median positive
SHAP importance or median nonzero absolute expression log2FC among the fixed
candidate features for that task. A value equal to its task anchor receives
0.5. Anchors are fitted after detection/log2FC prefiltering and before the
direction and Model-expression pattern filters, then reused inside the
Cross-patient calculation where applicable.

Cross-patient SHAP-importance retention is calculated from raw bootstrap
importance, without the nonlinear main-score transformation. The fixed
log2FC anchor is reused when constructing each bootstrap Expression signal.
Anchors are task-specific references: the resulting scores are not
probabilities, and changing the candidate-feature set can change the anchors.

The goal is a ranked marker table that is useful for biological follow-up, not
just a model accuracy report.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("MaxAlimadadi/SHAPTURE")
```

Some workflows use Bioconductor packages such as `SingleCellExperiment`,
`SummarizedExperiment`, and `CATALYST`. Install those with Bioconductor if your
analysis starts from SCE / CATALYST objects.

```r
install.packages("BiocManager")
BiocManager::install(c(
  "SingleCellExperiment",
  "SummarizedExperiment",
  "CATALYST",
  "ComplexHeatmap"
))
```

## Input data

The modeling functions expect a table with one row per cell or metacell. At
minimum, the table should contain:

- `condition`: the two biological groups to compare
- `sample_id`: patient, donor, or sample identifier
- optional `split_by_col`: subset, cluster, cell type, sex, batch, or another
  stratifying label
- feature columns: genes, proteins, markers, or other numeric predictors

Example structure:

| condition | sample_id | split_by_col | GeneA | GeneB | GeneC |
|---|---|---|---:|---:|---:|
| control | S1 | subset_1 | 0.14 | 1.22 | 0.00 |
| disease | S2 | subset_1 | 0.51 | 0.92 | 0.07 |

Template scripts for preparing data and running analyses are in `pipeline/`.
They are meant to be copied or edited for a project, not treated as fixed
package internals.

## Quick start

After preparing `dataset_all` and `what_genes`, define a SHAPTURE config and
run the pipeline.

```r
library(SHAPTURE)

shapture_config <- list(
  split_by_col = "your_subset_column",
  subset_levels = NULL,
  split_levels = NULL,
  comparison_mode = "condition_vs_condition",
  what_0 = "control",
  what_1 = "disease",
  iteration_no = 50,
  gene_fraction_per_iteration = 0.30,
  cross_patient_bootstraps = 500L,
  analysis_seed = 123L,
  top_feature_label = 30,
  top_feature_show = 50,
  feature_min_pct = 0.25,
  feature_min_log2FC = 0.25,
  shap_expr_cor_cutoff = 0.25,
  if_adaptive_weights = FALSE,
  score_scaling_method = "median_anchor",
  score_scaling_probs = c(0.05, 0.95),
  analysis_unit_type = "single_cell",
  memory_mode = "safe",
  model_type = "xgboost",
  preprocessing_method = "scale",
  patient_balance_method = "model_weights",
  patient_split_strategy = "shared_patients",
  shap_expr_cor_direction_filter = "both",
  group_colors = c("darkblue", "darkred"),
  continuous_colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
  n_cores = 4,
  project_name = NULL
)

pipeline_config <- resolve_shapture_script_config(
  envir = list2env(list(shapture_config = shapture_config)),
  dataset_all = dataset_all,
  what_genes = what_genes,
  main_dir = getwd()
)

do.call(run_shapture_pipeline, pipeline_config)
```

For a complete parameter reference, see:

```r
pipeline/shapture_config_full_reference.R
```

### Reproducibility

`analysis_seed = 123L` is the default reproducibility root. SHAPTURE derives
separate deterministic seeds for every task, iteration, cell/metacell split,
feature subset, model fit, metacell group, balancing step, and Cross-patient
bootstrap. Existing task seeds do not depend on task order, so adding another
split level does not change previously selected tasks. Increasing
`iteration_no` also preserves the seed plan for the earlier iterations.

Keep the same input data and row order, package versions, settings, and
`analysis_seed` when an analysis must be reproduced exactly. Change
`analysis_seed` deliberately when running an independent stochastic repeat.
Every completed task records its complete seed plan in
`docs/reproducibility_seeds.csv`.

### Parallel execution

Set `n_cores` to the maximum number of workers you want SHAPTURE to use. Model
workers are stopped before Cross-patient bootstrapping begins, so the two
parallel stages do not coexist in memory. The bootstrap controller then selects
the smaller of the requested, CPU-available, and memory-safe worker counts.
This controller uses PSOCK workers, which are portable across macOS, Windows,
and Linux. If safe parallel execution cannot be confirmed, the task continues
serially with the same deterministic bootstrap plans.

The applied backend, resource estimates, worker count, timing, and fallback
reason are saved in `docs/cross_patient_bootstrap_runtime.csv`. These execution
details are also summarized in `Integrated_adaptive_scoring_results/run_summary.csv`.

Every task also records the fixed SHAP-importance and absolute-log2FC anchors
in `docs/score_anchor_summary.csv`. The legacy `score_scaling_method` values
`"winsorized"` and `"minmax"` remain accepted so older scripts do not fail,
but evidence scoring follows the median-anchor/native-bounded strategy.
`score_scaling_probs` is used only when optional adaptive weighting estimates
robust component spread.

### Patient train/test strategies

`patient_split_strategy = "shared_patients"` is the default and preserves the
original SHAPTURE behavior: test cells or metacells are held out, but their
patients may also contribute different units to training.

`patient_split_strategy = "unseen_patients"` keeps patient IDs globally
separate between training and testing. Within each task, eligible patients are
split using rounded 60% training and 40% testing counts. The method requires at
least five eligible patients per condition, at least 20 cells per
patient-condition in single-cell mode or 5 metacells in metacell mode, and at
least 50 training cells or 20 training metacells per condition. Patient
selection rotates across iterations so held-out appearances remain balanced.

If a task cannot satisfy these requirements, that whole task uses
`shared_patients` and records the reason. A one-vs-others split is always global:
no cells from a held-out patient can appear in training, regardless of whether
they belong to the target or pooled background.

## Analysis modes

### Condition vs condition

Use this when you want to compare two biological groups across all cells or
within each selected level of `split_by_col`.

### One vs others

Use this when each selected level should become the positive class once. For
example, if `split_levels = c("A", "B", "C")`, SHAPTURE creates separate
result folders for:

- `A` vs all other retained levels
- `B` vs all other retained levels
- `C` vs all other retained levels

`subset_levels` controls which levels are retained before the task is built.
`split_levels` controls which retained levels are analyzed as targets.

## Feature filtering and model-expression pattern

Before modeling each task, SHAPTURE can filter features using:

- detection rate in the original single-cell task data
- absolute expression log2FC in the original single-cell task data

After SHAP values are calculated, SHAPTURE keeps features with interpretable
model-expression patterns. A retained feature should have enough expression
change and enough value-to-SHAP correlation to support one of these statements:

- higher feature values push the model toward condition 1 or target
- higher feature values push the model toward condition 0 or reference

Features without enough directional evidence are removed before final ranking.

## Outputs

Each analysis task creates a result folder with:

- `docs/combined_results.csv`: ranked feature table
- `docs/shapture_config_used.csv`: resolved parameters used for that run
- `docs/feature_prefilter_summary.csv`: feature filtering summary
- `docs/sampling_coverage_summary.csv`: held-out cell or metacell coverage
- `docs/feature_bagging_summary.csv`: gene-subsetting coverage
- `docs/reproducibility_seeds.csv`: analysis, task, stage, and per-iteration
  seeds used to reproduce stochastic operations
- `docs/patient_split_summary.csv`: requested and applied patient split strategy
- `docs/patient_sampling_coverage_summary.csv`: patient eligibility and rotation
- `docs/patient_feature_coverage_summary.csv`: held-out patient coverage by gene
- `docs/cross_patient_bootstrap_runtime.csv`: applied serial or parallel
  bootstrap backend, CPU/RAM safeguards, timing, and fallback status
- `docs/equal_count_shap_group_summary.csv`: feature-specific equal-count
  groups used by the binned beeswarm and dependency plots
- `Integrated_adaptive_scoring_results/run_summary.csv`: compact task settings,
  patient split outcome, coverage, feature counts, model performance, and
  Cross-patient bootstrap runtime summary
- `Integrated_adaptive_scoring_results/`: summary figures
- `Integrated_adaptive_scoring_results/shap_value_dependency_plots/`: one
  SHAP-value dependency PDF for every feature shown in the beeswarm
- `SHAP_only/`: SHAP-only condition summaries

The pipeline also saves `docs/visualization_inputs.rds`, which lets you redraw
figures later without rerunning the model.

## Redrawing figures after analysis

After a run finishes, edit:

```r
pipeline/run_shapture_visualization.R
```

Set `result_dir` to a completed result folder, then rerun the script. This
lets you change top-feature counts, figure sizes, colors, and label sizes
without retraining models.

The editable settings near the top of the script include
`top_feature_show`, `top_feature_label`, and `top_feature_show_extended`.
Figures whose layout depends on feature count are sized automatically from the
actual number displayed. `width_per_feature` and `height_per_feature` control
the global growth rate; plot-specific base sizes and minimum/maximum safeguards
are handled by `calculate_shapture_figure_size()`. To force any figure to an
exact size, replace the calculated width or height directly in its visible
`pdf()` or `png()` call.

Density plots use a separate balanced layout: their actual available features
are arranged with approximately equal numbers of rows and columns, and the PDF
allocates consistent physical space to each facet plus horizontal room for the
legend. The feature-specific axis ranges remain free. For example, 10, 20, 30,
and 50 features are arranged as 4 by 3, 5 by 4, 6 by 5, and 8 by 7 grids.

To redraw several completed result folders with the same settings, edit:

```r
pipeline/run_shapture_visualization_loop_multiple_folder.R
```

Set `analysis_root` to their immediate parent directory. The script detects
subfolders containing both `docs/combined_results.csv` and
`docs/visualization_inputs.rds`, so it works for condition-vs-condition and
one-vs-others results without depending on folder names.

Feature names in the beeswarm and feature-importance plots use the same
`Model_expression_pattern` colors: `group_colors[2]` for condition 1/target
and `group_colors[1]` for condition 0/others. Both legends list condition 1
first. These colors describe the saved class assignment, not the score magnitude.

The original beeswarm remains an unbinned unit-level plot. A second binned
beeswarm and the feature-specific dependency PDFs use the same deterministic
groups for each feature. The script aims for 500 groups, keeps at least 10
valid cells or metacells in each group, and displays original units when fewer
than 500 are available. Groups are formed jointly from SHAP and model-seen
feature values after mapping both dimensions to 0-1 for grouping only. Each
filled geometric hexagon is placed at its group's mean SHAP and mean feature
value. Set `binned_beeswarm_hexagon_size` in either editable visualization
script to make the binned-beeswarm symbols smaller or larger; the default is
`0.98`. This binning is only for display: scores, correlations, rankings, and
dependency smooths continue to use the original unbinned values. Output names
use the number of
features actually displayed rather than the requested number when fewer saved
features are available.

Dependency plots reserve their explanatory sentence above the coordinate panel
and add a 3% display margin around both axes so boundary hexagons remain fully
visible. This margin changes only the displayed limits.

## Global visualization across split levels

When an analysis produces several split-level result folders, edit:

```r
pipeline/run_shapture_global_visualization.R
```

Set `results_root` to the folder containing the completed result folders. The
script reads each `docs/combined_results.csv` and does not rerun any models.
`top_feature_show` selects the top-ranked features independently from every
folder; their union is used in two global figures:

- a hierarchically clustered cross-level final-score landscape
- a direction-aware Jaccard similarity matrix

The direction-aware Jaccard numerator contains features shared by two folders
with the same normalized model-expression pattern. Its denominator is the
union of their selected feature names. Therefore, a feature found in both
folders but pointing in opposite directions does not count as a directional
match. Condition-vs-condition folders must use the same condition labels. For
one-vs-others results, directions are normalized to target versus others before
folders are compared.

Both figures use average-linkage hierarchical clustering by default. Their PDF
dimensions are calculated after loading the results, using the actual plotted
feature and split-level counts. The cross-level landscape grows with both
matrix dimensions, while the square Jaccard view grows with the split-level
count. The visible `width` and `height` controls default to `NULL`; replace
either one with an exact size in inches to override only that dimension.
Palette, tile borders, feature-name display, direction symbols, and similarity
value labels are also directly editable in the visualization script.
The default palette is:

```r
rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
```

For very large feature-by-level matrices, the heatmap body is rasterized
automatically while labels and dendrograms remain vector elements. Set
`use_raster = FALSE` in the editable script to force an entirely vector PDF.

## Interpretation notes

SHAPTURE is a discovery workflow. It is built to prioritize features that are
important to the model, directionally connected to expression, supported across
patients when possible, and reproducible across repeated resampling. It does
not replace experimental validation or fully sample-level inference.

Cross-patient support requires at least two patients in each condition for a
condition-vs-condition task, or at least two patients in total for a
one-vs-others task. This is a resampling-based support measure, not a p-value or
independent validation.

## Development status

SHAPTURE is under active development. The current focus is making the workflow
cleaner as an R package while preserving the existing analysis behavior.

## License

MIT. See `LICENSE.md`.
