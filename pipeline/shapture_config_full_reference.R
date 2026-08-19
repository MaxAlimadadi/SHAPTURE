#####
# SHAPTURE full config reference
#####

# This file is a user-facing reference for the full unified `shapture_config`
# schema. The pipeline still uses one shared config structure internally, but
# the getter scripts expose shorter mode-specific presets for day-to-day use.
#
# Typical workflow:
# 1. prepare `dataset_all` and `what_genes` using one of the getter scripts
# 2. copy or adapt this list if you need the full configuration surface
# 3. assign your chosen list to `shapture_config`

shapture_config_full <- list(
  #
  # ---- Analysis structure ----
  #
  # metadata column used to split the analysis; use NULL for one global run
  # in one_vs_others mode, this should be the target-defining column
  split_by_col = NULL,
  
  # values of split_by_col to keep before analysis; NULL keeps all levels
  # in one_vs_others mode, this defines the universe used to build "others"
  subset_levels = NULL,
  
  # values of split_by_col to analyze; NULL means use all levels after subset_levels
  # in one_vs_others mode, these are the target levels that each get a folder
  split_levels = NULL,
  
  # "condition_vs_condition" or "one_vs_others"
  comparison_mode = "condition_vs_condition",
  
  # optional explicit target-defining column for one-vs-others mode; usually
  # leave NULL and let the pipeline use split_by_col automatically
  one_vs_others_col = NULL,
  
  # optional explicit target levels for one-vs-others mode; usually leave NULL
  # and let the pipeline use split_levels automatically
  one_vs_others_levels = NULL,
  
  # negative-class label used in one-vs-others mode
  one_vs_others_negative_label = "others",
  
  # labels used in condition_vs_condition mode
  what_0 = "low",
  what_1 = "high",
  
  # optional output prefix; NULL builds one automatically
  project_name = NULL,
  
  #
  # ---- Resampling ----
  #
  # proportion of sampled rows per condition used for training
  train_ratio = 0.7,
  
  # proportion of sampled rows per condition used for testing
  test_ratio = 0.3,
  
  # number of repeated train/test/modeling iterations
  iteration_no = 30,
  
  #
  # ---- Ranking and display ----
  #
  # number of features to label in scatter-type plots
  top_feature_label = 30,
  
  # number of top features to show/save in visualization outputs; often length(what_genes)
  top_feature_show = length(what_genes),
  
  # whether to adapt score weights to component spread in each task; FALSE uses
  # default_prior_scores_weights directly
  if_adaptive_weights = FALSE,

  # "winsorized" caps extreme component values before 0-1 scaling; "minmax"
  # uses direct min-max scaling
  score_scaling_method = "winsorized",

  # lower/upper quantiles used by winsorized component scaling and by
  # adaptive-weight spread estimates when adaptive weights are enabled
  score_scaling_probs = c(0.05, 0.95),
  
  # keep features detected (>0) in more than this fraction of cells in at least
  # one task group; 0 disables this filter
  feature_min_pct = 0.25,

  # keep features with absolute expression log2FC at least this large in the
  # original single-cell task data; 0 disables this filter
  feature_min_log2FC = 0.25,

  # require at least this absolute value-SHAP correlation for
  # Model_expression_pattern classification
  shap_expr_cor_cutoff = 0.25,

  # prior weights used by adaptive scoring before data-driven adjustment
  default_prior_scores_weights = c(
    SHAP_importance_prior_weight = 0.25,
    expression_signal_prior_weight = 0.25,
    value_SHAP_relationship_prior_weight = 0.30,
    patient_level_support_prior_weight = 0.10,
    model_reproducibility_prior_weight = 0.10
  ),
  
  #
  # ---- Modeling ----
  #
  # "none" keeps input values as-is; "scale" scales features using training data only
  preprocessing_method = "scale",
  
  # "xgboost" or "lightgbm"
  model_type = "lightgbm",
  
  # "single_cell" or "metacell"
  analysis_unit_type = "metacell",
  
  # "full" keeps pooled SHAP payloads in memory; "safe" uses temporary
  # iteration chunks on disk, rebuilds pooled SHAP summaries one chunk at a
  # time, and deletes the temporary files automatically at the end
  memory_mode = "full",
  
  #
  # ---- Sample-aware balancing ----
  #
  # "hard_median_downsample", "weighted_sampling", or "model_weights"
  # in one_vs_others mode, the target level is balanced by sample_id and the
  # pooled non-target class is balanced by sample_id + split_by_col
  patient_balance_method = "model_weights",
  
  # "both", "positive", or "negative"
  # "positive" keeps shap_expr_cor > 0; "negative" keeps shap_expr_cor < 0
  shap_expr_cor_direction_filter = "both",
  
  #
  # ---- Plot colors ----
  #
  # two unnamed colors: first = left/condition_0/others, second = right/condition_1/target
  group_colors = c("darkblue", "darkred"),
  
  # continuous palette for SHAP/expression/support gradients; any valid R color
  # vector works, for example:
  # c("navy", "firebrick")
  # c("navy", "white", "firebrick")
  # viridisLite::viridis(100)
  # RColorBrewer::brewer.pal(11, "RdBu")
  continuous_colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
  
  #
  # ---- Parallel execution ----
  #
  # number of parallel model-fitting workers; use 1 for sequential execution
  n_cores = 14,
  
  #
  # ---- Metacell settings (used only when analysis_unit_type = "metacell") ----
  #
  # minimum cells needed in one sample/subset block before metacells are built
  metacell_min_cells_to_build = 40,
  
  # minimum cells per metacell
  metacell_min_cells_per_metacell = 5,
  
  # maximum cells per metacell
  metacell_max_cells_per_metacell = 25,
  
  # number of variable features used internally for metacell neighborhood building
  metacell_n_variable_features = 2000,
  
  # number of PCs used for metacell neighborhood building
  metacell_n_pcs = 20,
  
  # "mean" or "sum"
  metacell_aggregation_fun = "mean",
  
  # parallelize metacell construction
  parallel_metacells = TRUE,
  
  # number of workers for metacell construction
  metacell_n_cores = 10,
  
  # "serial", "mclapply", or "psock"
  metacell_parallel_backend = "psock"
)

# Example:
# shapture_config <- shapture_config_full
