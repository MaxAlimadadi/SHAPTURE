#####
# Prepare SHAPTURE input from a SingleCellExperiment / CATALYST object
#
# Edit this file for your dataset, run it first, then run
# pipeline/run_shapture_analysis.R.
#####

library(CATALYST)
library(dplyr)
library(RColorBrewer)
library(SHAPTURE)

#####
# 1. Load your object
#####
sce_path <- "path/to/your_sce_object.rds"
sce <- readRDS(sce_path)

#####
# 2. Define metadata columns and condition levels
#####
condition_col <- "condition"       # binary biological condition column in colData(sce)
sample_col <- "sample_id"          # patient / donor / sample identifier column
split_by_col <- "subset_label"     # use NULL for one global analysis

# If split_by_col should be populated from CATALYST cluster_ids(), set the
# clustering resolution key here. Otherwise use NULL and read split_by_col
# directly from colData(sce).
cluster_resolution <- NULL

# Put the reference/left condition first and the target/right condition second.
condition_levels <- c("condition_0", "condition_1")

#####
# 3. Choose features
#####
# NULL uses CATALYST type_markers(sce) and state_markers(sce).
# You can also provide a custom vector, for example:
# feature_cols <- c("MarkerA", "MarkerB", "MarkerC")
feature_cols <- NULL

#####
# 4. Build SHAPTURE input table
#####
data_prep_config <- list(
  split_by_col = split_by_col,
  cluster_resolution = cluster_resolution,
  condition_col = condition_col,
  sample_col = sample_col,
  feature_cols = feature_cols
)

dataset_res <- do.call(
  build_dataset_from_sce,
  c(list(sce = sce), data_prep_config)
)

dataset_all <- dataset_res$dataset_all
what_genes <- dataset_res$what_genes

dataset_all <- dataset_all[
  dataset_all$condition %in% condition_levels,
  ,
  drop = FALSE
]
dataset_all$condition <- factor(dataset_all$condition, levels = condition_levels)

#####
# 5. Define analysis settings
#
# Keep one unified pipeline schema underneath, but expose two compact
# user-facing presets here. Pick one at the bottom and assign it to
# shapture_config.
#####

shapture_config_condition_vs_condition <- list(
  # subset column; use NULL for one global condition-vs-condition run
  split_by_col = data_prep_config$split_by_col,
  
  # optional split_by_col levels to keep before analysis; NULL keeps all levels
  subset_levels = NULL,
  
  # optional task levels to analyze; NULL uses all levels after subset_levels
  split_levels = NULL,
  
  comparison_mode = "condition_vs_condition",
  what_0 = condition_levels[1],
  what_1 = condition_levels[2],
  iteration_no = 30,
  
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
  
  # "winsorized" caps extreme component values before 0-1 scaling; "minmax" uses direct min-max scaling
  score_scaling_method = "winsorized",
  
  # lower/upper quantiles for winsorized component scaling and adaptive-weight spread estimates
  score_scaling_probs = c(0.05, 0.95),
  
  # "single_cell" or "metacell"
  analysis_unit_type = "metacell",
  
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
  
  n_cores = 4,
  project_name = NULL
)

shapture_config_one_vs_others <- list(
  # target-defining column for one-vs-others mode
  split_by_col = data_prep_config$split_by_col,
  
  # optional split_by_col levels to keep before analysis; NULL keeps all levels.
  # In one_vs_others, this defines the universe used to build the pooled "others" class.
  subset_levels = NULL,
  
  # optional target levels to analyze; NULL uses all levels after subset_levels
  split_levels = NULL,
  
  comparison_mode = "one_vs_others",
  iteration_no = 30,
  
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
  
  # "winsorized" caps extreme component values before 0-1 scaling; "minmax" uses direct min-max scaling
  score_scaling_method = "winsorized",
  
  # lower/upper quantiles for winsorized component scaling and adaptive-weight spread estimates
  score_scaling_probs = c(0.05, 0.95),
  
  # "single_cell" or "metacell"
  analysis_unit_type = "metacell",
  
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
  
  n_cores = 4,
  project_name = NULL
)

# Pick the active config here.
shapture_config <- shapture_config_condition_vs_condition

if (!is.null(shapture_config$split_by_col)) {
  dataset_all[[shapture_config$split_by_col]] <- factor(dataset_all[[shapture_config$split_by_col]])
}

rm(sce)
