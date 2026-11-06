#####
# Prepare SHAPTURE input from a Seurat object
#
# Edit this file for your dataset, run it first, then run
# pipeline/run_shapture_analysis.R.
#####

library(Seurat)
library(dplyr)
library(RColorBrewer)
library(SHAPTURE)

#####
# 1. Load your object
#####
seurat_path <- "path/to/your_seurat_object.rds"
my_obj <- readRDS(seurat_path)

#####
# 2. Define metadata columns and condition levels
#####
condition_col <- "condition"      # binary biological condition column in meta.data
sample_col <- "sample_id"         # patient / donor / sample identifier column
split_level_column <- "subset_label"    # use NULL for one global analysis

# Put the reference/left condition first and the target/right condition second.
condition_levels <- c("condition_0", "condition_1")

#####
# 3. Choose assay, layer, and features
#####
assay <- "RNA"
layer <- "data"

DefaultAssay(my_obj) <- assay
# Keeps legacy assays unchanged; joins layers only for assay classes that support it.
my_obj <- join_seurat_layers_if_supported(my_obj, assay = assay)

# Use your preferred feature set here. Examples:
# what_genes <- VariableFeatures(my_obj)
# what_genes <- rownames(my_obj)
# what_genes <- c("GeneA", "GeneB", "GeneC")
what_genes <- VariableFeatures(my_obj)
if (length(what_genes) == 0) {
  what_genes <- rownames(my_obj)
}

#####
# 4. Build SHAPTURE input table
#####
data_prep_config <- list(
  condition_col = condition_col,
  sample_col = sample_col,
  split_level_column = split_level_column,
  assay = assay,
  layer = layer,
  genes = what_genes,
  condition_levels = condition_levels
)

dataset_res <- do.call(
  build_dataset_from_seurat,
  c(list(seurat_obj = my_obj), data_prep_config)
)

dataset_all <- dataset_res$dataset_all
what_genes <- dataset_res$what_genes
# Always use the prepared-data split name returned by the adapter. It is the
# original name unless that name would collide with an internal metadata field.
analysis_split_level_column <- dataset_res$split_level_column

#####
# 5. Define analysis settings
#
# Keep one unified pipeline schema underneath, but expose two compact
# user-facing presets here. Pick one at the bottom and assign it to
# shapture_config.
#####

# Result-folder prefix for this run; change it to keep analyses separate.
# NULL uses automatic comparison names. Reusing a name in the same output
# directory reuses the corresponding result folders.
project_name <- "analysis_name"

shapture_config_condition_vs_condition <- list(
  # subset column; use NULL for one global condition-vs-condition run
  split_level_column = analysis_split_level_column,
  
  # optional task levels to analyze; NULL uses all levels present in the input
  # filter the input data first if you want to remove levels entirely
  split_levels_to_analyze = NULL,
  
  comparison_mode = "condition_vs_condition",
  what_0 = data_prep_config$condition_levels[1],
  what_1 = data_prep_config$condition_levels[2],
  iteration_no = 50,

  # fraction of filtered candidate features included in each model iteration
  gene_fraction_per_iteration = 0.30,

  # number of patient-block bootstrap replicates for Cross-patient Support
  cross_patient_bootstraps = 500L,

  # reproducibility root; keep the same value to reproduce stochastic stages
  analysis_seed = 123L,
  
  # number of features to label in scatter-type plots
  top_feature_label = 30,
  
  # number of top features to show/save in visualization outputs
  top_feature_show = 50,
  
  # keep features detected (>0) in more than this fraction of cells in at least one task group; 0 disables this filter
  feature_min_pct = 0.25,
  
  # keep features with absolute expression log2FC at least this large in original single-cell task data; 0 disables this filter
  feature_min_log2FC = 0.25,
  
  # require at least this absolute value-SHAP correlation for Model_expression_pattern classification
  shap_expr_cor_cutoff = 0.25,
  
  # prior weights used by adaptive scoring before data-driven adjustment
  default_prior_scores_weights = c(
    SHAP_importance_prior_weight = 0.30,
    expression_signal_prior_weight = 0.30,
    value_SHAP_relationship_prior_weight = 0.20,
    cross_patient_support_prior_weight = 0.10,
    model_reproducibility_prior_weight = 0.10
  ),
  
  # FALSE uses the prior weights directly; TRUE adapts them to component spread in this task
  if_adaptive_weights = FALSE,
  
  # median-anchor scoring for unbounded evidence; legacy values remain accepted for script compatibility
  score_scaling_method = "median_anchor",
  
  # lower/upper quantiles used only for adaptive-weight spread estimates
  score_scaling_probs = c(0.05, 0.95),
  
  # "single_cell" or "metacell"
  analysis_unit_type = "single_cell",
  
  # "full" keeps pooled SHAP payloads in RAM; "safe" uses temporary disk chunks
  memory_mode = "safe",
  
  # "xgboost" or "lightgbm"
  model_type = "xgboost",
  
  # "none" or "scale"
  preprocessing_method = "scale",
  
  # "hard_median_downsample", "weighted_sampling", or "model_weights"
  patient_balance_method = "model_weights",

  # "shared_patients" (default) or "unseen_patients"
  # unseen_patients uses globally separate 60%/40% patient train/test folds
  # and automatically falls back per task when patient/unit counts are too low
  patient_split_strategy = "shared_patients",
  
  # "both", "positive", or "negative"
  shap_expr_cor_direction_filter = "both",
  
  # first = left/what_0, second = right/what_1
  group_colors = c("darkblue", "darkred"),
  
  # any valid R color vector with at least 2 colors
  continuous_colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
  
  # requested maximum for model and Cross-patient bootstrap workers;
  # bootstrap workers are automatically reduced when CPU or RAM is limited
  n_cores = 4,
  # Use the run name defined above as the result-folder prefix.
  project_name = project_name
)

shapture_config_one_vs_others <- list(
  # target-defining column for one-vs-others mode
  split_level_column = analysis_split_level_column,
  
  # optional targets to analyze; NULL uses all levels present in the input
  # every other level in the input remains available to the pooled "others" class
  split_levels_to_analyze = NULL,
  
  comparison_mode = "one_vs_others",
  iteration_no = 50,

  # fraction of filtered candidate features included in each model iteration
  gene_fraction_per_iteration = 0.30,

  # number of patient-block bootstrap replicates for Cross-patient Support
  cross_patient_bootstraps = 500L,

  # reproducibility root; keep the same value to reproduce stochastic stages
  analysis_seed = 123L,
  
  # number of features to label in scatter-type plots
  top_feature_label = 30,
  
  # number of top features to show/save in visualization outputs
  top_feature_show = 50,
  
  # keep features detected (>0) in more than this fraction of cells in at least one task group; 0 disables this filter
  feature_min_pct = 0.25,
  
  # keep features with absolute expression log2FC at least this large in original single-cell task data; 0 disables this filter
  feature_min_log2FC = 0.25,
  
  # require at least this absolute value-SHAP correlation for Model_expression_pattern classification
  shap_expr_cor_cutoff = 0.25,
  
  # prior weights used by adaptive scoring before data-driven adjustment
  default_prior_scores_weights = c(
    SHAP_importance_prior_weight = 0.30,
    expression_signal_prior_weight = 0.30,
    value_SHAP_relationship_prior_weight = 0.20,
    cross_patient_support_prior_weight = 0.10,
    model_reproducibility_prior_weight = 0.10
  ),
  
  # FALSE uses the prior weights directly; TRUE adapts them to component spread in this task
  if_adaptive_weights = FALSE,
  
  # median-anchor scoring for unbounded evidence; legacy values remain accepted for script compatibility
  score_scaling_method = "median_anchor",
  
  # lower/upper quantiles used only for adaptive-weight spread estimates
  score_scaling_probs = c(0.05, 0.95),
  
  # "single_cell" or "metacell"
  analysis_unit_type = "single_cell",
  
  # "full" keeps pooled SHAP payloads in RAM; "safe" uses temporary disk chunks
  memory_mode = "safe",
  
  # "xgboost" or "lightgbm"
  model_type = "xgboost",
  
  # "none" or "scale"
  preprocessing_method = "scale",
  
  # "hard_median_downsample", "weighted_sampling", or "model_weights"
  patient_balance_method = "model_weights",

  # "shared_patients" (default) or "unseen_patients"
  # unseen_patients keeps every test patient completely absent from training
  patient_split_strategy = "shared_patients",
  
  # "both", "positive", or "negative"
  shap_expr_cor_direction_filter = "both",
  
  # first = others, second = target
  group_colors = c("darkblue", "darkred"),
  
  # any valid R color vector with at least 2 colors
  continuous_colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
  
  # requested maximum for model and Cross-patient bootstrap workers;
  # bootstrap workers are automatically reduced when CPU or RAM is limited
  n_cores = 4,
  # Use the run name defined above as the result-folder prefix.
  project_name = project_name
)

# Pick the active config here.
shapture_config <- shapture_config_condition_vs_condition

if (!is.null(shapture_config$split_level_column)) {
  dataset_all[[shapture_config$split_level_column]] <- factor(dataset_all[[shapture_config$split_level_column]])
}

rm(my_obj)
