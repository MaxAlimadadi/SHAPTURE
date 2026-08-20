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
- patient-level support
- model reproducibility

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
  iteration_no = 30,
  top_feature_label = 30,
  top_feature_show = 50,
  feature_min_pct = 0.25,
  feature_min_log2FC = 0.25,
  shap_expr_cor_cutoff = 0.25,
  if_adaptive_weights = FALSE,
  score_scaling_method = "winsorized",
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
- `docs/patient_split_summary.csv`: requested and applied patient split strategy
- `docs/patient_sampling_coverage_summary.csv`: patient eligibility and rotation
- `docs/patient_feature_coverage_summary.csv`: held-out patient coverage by gene
- `Integrated_adaptive_scoring_results/run_summary.csv`: compact task settings,
  patient split outcome, coverage, feature counts, and model performance
- `Integrated_adaptive_scoring_results/`: summary figures
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

Both figures use average-linkage hierarchical clustering by default. Width,
height, palette, tile borders, feature-name display, direction symbols, and
similarity value labels are directly editable in the visualization script.
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

For very small patient counts, permutation p-values have a hard resolution
limit. SHAPTURE records this in the output files and treats those runs more
cautiously.

## Development status

SHAPTURE is under active development. The current focus is making the workflow
cleaner as an R package while preserving the existing analysis behavior.

## License

MIT. See `LICENSE.md`.
