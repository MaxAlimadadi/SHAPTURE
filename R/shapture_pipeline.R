#' Null-coalescing helper
#'
#' @param x Primary value.
#' @param y Fallback value.
#'
#' @return `x` when not `NULL`, otherwise `y`.
#' @noRd
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


#' Load the runtime packages used by SHAPTURE
#'
#' @return Invisibly returns `TRUE`.
#' @keywords internal
load_shapture_runtime <- function() {
  required_packages <- c(
    "shapviz", "ggplot2", "xgboost", "dplyr", "scales", "pROC",
    "doParallel", "ggbeeswarm", "ggrepel", "ggtext", "tidyr", "cowplot", "reshape2",
    "foreach", "lightgbm"
  )
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0) {
    stop(
      "Missing required packages for SHAPTURE runtime: ",
      paste(missing_packages, collapse = ", ")
    )
  }
  
  invisible(TRUE)
}


#' Build default SHAPTURE script configuration
#'
#' @param what_genes Feature vector used to determine default display counts.
#'
#' @return Named list of default configuration values.
#' @keywords internal
default_shapture_script_config <- function(what_genes) {
  list(
    split_by_col = NULL,
    subset_levels = NULL,
    split_levels = NULL,
    comparison_mode = "condition_vs_condition",
    one_vs_others_col = NULL,
    one_vs_others_levels = NULL,
    one_vs_others_negative_label = "others",
    what_0 = "low",
    what_1 = "high",
    project_name = NULL,
    train_ratio = 0.7,
    test_ratio = 0.3,
    iteration_no = 50,
    gene_fraction_per_iteration = 0.30,
    cross_patient_bootstraps = 500L,
    analysis_seed = 123L,
    top_feature_label = 30,
    top_feature_show = length(what_genes),
    if_adaptive_weights = FALSE,
    feature_min_pct = 0.25,
    feature_min_log2FC = 0.25,
    shap_expr_cor_cutoff = 0.25,
    score_scaling_method = "median_anchor",
    score_scaling_probs = c(0.05, 0.95),
    default_prior_scores_weights = c(
      SHAP_importance_prior_weight = 0.25,
      expression_signal_prior_weight = 0.25,
      value_SHAP_relationship_prior_weight = 0.30,
      cross_patient_support_prior_weight = 0.10,
      model_reproducibility_prior_weight = 0.10
    ),
    preprocessing_method = "scale",
    model_type = "lightgbm",
    analysis_unit_type = "metacell",
    memory_mode = "full",
    patient_balance_method = "model_weights",
    patient_split_strategy = "shared_patients",
    shap_expr_cor_direction_filter = "both",
    group_colors = c("darkblue", "darkred"),
    continuous_colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
    n_cores = 14,
    metacell_min_cells_to_build = 40,
    metacell_min_cells_per_metacell = 5,
    metacell_max_cells_per_metacell = 25,
    metacell_n_variable_features = 2000,
    metacell_n_pcs = 20,
    metacell_aggregation_fun = "mean",
    parallel_metacells = TRUE,
    metacell_n_cores = 10,
    metacell_parallel_backend = "psock"
  )
}


#' Validate a user-supplied color vector
#'
#' @param colors Character vector of R colors.
#' @param arg_name Argument name used in error messages.
#' @param min_length Minimum allowed length.
#' @param exact_length Optional exact required length.
#'
#' @return Unnamed character vector of validated colors.
#' @keywords internal
validate_color_vector <- function(colors, arg_name, min_length = 1, exact_length = NULL) {
  if (!is.character(colors)) {
    stop(arg_name, " must be a character vector of valid R colors.")
  }
  if (!is.null(exact_length) && length(colors) != exact_length) {
    stop(arg_name, " must contain exactly ", exact_length, " colors.")
  }
  if (length(colors) < min_length) {
    stop(arg_name, " must contain at least ", min_length, " colors.")
  }
  tryCatch(
    grDevices::col2rgb(colors),
    error = function(e) {
      stop(arg_name, " contains invalid R color values.")
    }
  )
  unname(colors)
}


#' Build a two-group color map from unnamed left/right colors
#'
#' @param what_0 Left or negative-class label.
#' @param what_1 Right or positive-class label.
#' @param group_colors Unnamed two-color vector.
#'
#' @return Named color vector for the two labels.
#' @keywords internal
build_group_color_map <- function(what_0, what_1, group_colors) {
  stats::setNames(unname(group_colors[seq_len(2)]), c(as.character(what_0), as.character(what_1)))
}


#' Normalize public comparison-mode aliases
#'
#' @param comparison_mode Comparison mode string.
#'
#' @return Canonical comparison mode string.
#' @keywords internal
normalize_comparison_mode <- function(comparison_mode) {
  if (identical(comparison_mode, "one_vs_rest")) {
    return("one_vs_others")
  }
  comparison_mode
}


#' Clean shapviz label placeholders that trigger ggplot warnings
#'
#' @param p A ggplot-compatible shapviz plot.
#'
#' @return Plot with sanitized axis labels.
#' @keywords internal
sanitize_shapviz_plot <- function(p) {
  if (inherits(p, "ggplot")) {
    if (!is.null(p$labels$y) && inherits(p$labels$y, "element_blank")) {
      p <- p + ggplot2::labs(y = NULL)
    }
    if (!is.null(p$labels$x) && inherits(p$labels$x, "element_blank")) {
      p <- p + ggplot2::labs(x = NULL)
    }
  }
  p
}


#' Serialize config values for human-readable export
#'
#' @param x Any config value.
#'
#' @return Character scalar.
#' @keywords internal
serialize_config_value <- function(x) {
  if (is.null(x)) {
    return("NULL")
  }
  if (length(x) == 0) {
    return("")
  }
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (inherits(x, "Date") || inherits(x, "POSIXt")) {
    return(paste(as.character(x), collapse = ", "))
  }
  if (is.atomic(x)) {
    if (!is.null(names(x)) && any(nzchar(names(x)))) {
      pieces <- paste0(names(x), "=", as.character(x))
      return(paste(pieces, collapse = "; "))
    }
    return(paste(as.character(x), collapse = ", "))
  }
  if (is.data.frame(x)) {
    return(paste0("data.frame[", nrow(x), " rows x ", ncol(x), " columns]"))
  }
  if (is.list(x)) {
    if (!is.null(names(x)) && any(nzchar(names(x)))) {
      pieces <- mapply(
        function(nm, val) paste0(nm, "=", serialize_config_value(val)),
        names(x),
        x,
        SIMPLIFY = TRUE,
        USE.NAMES = FALSE
      )
      return(paste(pieces, collapse = "; "))
    }
    return(paste(vapply(x, serialize_config_value, character(1)), collapse = "; "))
  }
  paste(capture.output(str(x, give.attr = FALSE)), collapse = " ")
}


#' Write resolved SHAPTURE config to a CSV file
#'
#' @param config_list Named list of resolved config values.
#' @param output_path CSV file path.
#'
#' @return Invisibly returns the output path.
#' @keywords internal
write_shapture_config_csv <- function(config_list, output_path) {
  if ("default_prior_scores_weights" %in% names(config_list)) {
    weight_name_map <- c(
      importance = "SHAP_importance_prior_weight",
      expression = "expression_signal_prior_weight",
      value_shap = "value_SHAP_relationship_prior_weight",
      cross_patient_support = "cross_patient_support_prior_weight",
      model_reproducibility = "model_reproducibility_prior_weight"
    )
    score_weights <- config_list$default_prior_scores_weights
    if (is.numeric(score_weights) &&
        !is.null(names(score_weights)) &&
        all(names(score_weights) %in% names(weight_name_map))) {
      score_weights <- score_weights[c(
        "importance",
        "expression",
        "value_shap",
        "cross_patient_support",
        "model_reproducibility"
      )]
      names(score_weights) <- unname(weight_name_map[names(score_weights)])
      config_list$default_prior_scores_weights <- score_weights
    }
  }

  config_table <- data.frame(
    parameter_name = names(config_list),
    parameter_used = vapply(config_list, serialize_config_value, character(1)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  utils::write.csv(config_table, output_path, row.names = FALSE)
  invisible(output_path)
}


#' Resolve script-level SHAPTURE configuration
#'
#' Reads user-defined objects from a script environment, fills defaults, applies
#' temporary environment overrides, validates options, and returns a clean
#' argument list for `run_shapture_pipeline()`.
#'
#' @param envir Environment to read user-specified objects from.
#' @param dataset_all Input data frame.
#' @param what_genes Feature vector.
#' @param main_dir Project directory.
#'
#' @return Named list of resolved pipeline arguments.
#' @export
resolve_shapture_script_config <- function(
    envir,
    dataset_all,
    what_genes,
    main_dir
) {
  config <- default_shapture_script_config(what_genes = what_genes)
  
  if (exists("shapture_config", envir = envir, inherits = FALSE)) {
    user_config <- get("shapture_config", envir = envir, inherits = FALSE)
    if (!is.list(user_config)) {
      stop("shapture_config must be a list.")
    }
    if (!is.null(user_config$one_vs_rest_col) && is.null(user_config$one_vs_others_col)) {
      user_config$one_vs_others_col <- user_config$one_vs_rest_col
    }
    if (!is.null(user_config$one_vs_rest_levels) && is.null(user_config$one_vs_others_levels)) {
      user_config$one_vs_others_levels <- user_config$one_vs_rest_levels
    }
    if (!is.null(user_config$one_vs_rest_negative_label) &&
        is.null(user_config$one_vs_others_negative_label)) {
      user_config$one_vs_others_negative_label <- user_config$one_vs_rest_negative_label
    }
    recognized_names <- intersect(names(user_config), names(config))
    config[recognized_names] <- user_config[recognized_names]
  }
  
  for (nm in names(config)) {
    if (exists(nm, envir = envir, inherits = FALSE)) {
      config[[nm]] <- get(nm, envir = envir, inherits = FALSE)
    }
  }
  if (exists("one_vs_rest_col", envir = envir, inherits = FALSE) &&
      is.null(config$one_vs_others_col)) {
    config$one_vs_others_col <- get("one_vs_rest_col", envir = envir, inherits = FALSE)
  }
  if (exists("one_vs_rest_levels", envir = envir, inherits = FALSE) &&
      is.null(config$one_vs_others_levels)) {
    config$one_vs_others_levels <- get("one_vs_rest_levels", envir = envir, inherits = FALSE)
  }
  if (exists("one_vs_rest_negative_label", envir = envir, inherits = FALSE) &&
      identical(config$one_vs_others_negative_label, "others")) {
    config$one_vs_others_negative_label <- get("one_vs_rest_negative_label", envir = envir, inherits = FALSE)
  }
  
  if (exists("condition_levels", envir = envir, inherits = FALSE)) {
    condition_levels <- get("condition_levels", envir = envir, inherits = FALSE)
  } else {
    condition_levels <- c(config$what_0, config$what_1)
  }
  
  temp_run_tag <- Sys.getenv("TEMP_RUN_TAG", unset = "")
  temp_iteration_override <- Sys.getenv("TEMP_ITERATION_NO", unset = "")
  temp_split_override <- Sys.getenv("TEMP_ONLY_SPLIT", unset = "")
  temp_cluster_override <- Sys.getenv("TEMP_ONLY_CLUSTER", unset = "")
  temp_analysis_override <- Sys.getenv("TEMP_ANALYSIS_UNIT_TYPE", unset = "")
  temp_balance_override <- Sys.getenv("TEMP_PATIENT_BALANCE_METHOD", unset = "")
  temp_shap_expr_cor_direction_override <- Sys.getenv("TEMP_SHAP_EXPR_COR_DIRECTION_FILTER", unset = "")
  temp_one_vs_others_override <- Sys.getenv("TEMP_ONE_VS_OTHERS_LEVEL", unset = "")
  if (!nzchar(temp_one_vs_others_override)) {
    temp_one_vs_others_override <- Sys.getenv("TEMP_ONE_VS_REST_LEVEL", unset = "")
  }
  
  if (nzchar(temp_iteration_override)) {
    config$iteration_no <- as.integer(temp_iteration_override)
  }
  if (nzchar(temp_analysis_override)) {
    config$analysis_unit_type <- temp_analysis_override
  }
  if (nzchar(temp_balance_override)) {
    config$patient_balance_method <- temp_balance_override
  }
  if (nzchar(temp_shap_expr_cor_direction_override)) {
    config$shap_expr_cor_direction_filter <- temp_shap_expr_cor_direction_override
  }

  # Object adapters record a mapping only when a source split name collides
  # with SHAPTURE's prepared condition/sample columns. This also keeps older
  # copied scripts that pass the source name into shapture_config working.
  config$split_by_col <- resolve_mapped_prepared_column(
    column_name = config$split_by_col,
    dataset_all = dataset_all
  )
  config$one_vs_others_col <- resolve_mapped_prepared_column(
    column_name = config$one_vs_others_col,
    dataset_all = dataset_all
  )
  
  config$comparison_mode <- normalize_comparison_mode(config$comparison_mode)
  
  if (!is.null(config$subset_levels) && is.null(config$split_by_col)) {
    stop("subset_levels requires split_by_col.")
  }
  
  if (!is.null(config$subset_levels)) {
    if (!config$split_by_col %in% colnames(dataset_all)) {
      stop("split_by_col is not present in dataset_all: ", config$split_by_col)
    }
    dataset_all <- dataset_all[
      as.character(dataset_all[[config$split_by_col]]) %in% as.character(config$subset_levels),
      ,
      drop = FALSE
    ]
    if (nrow(dataset_all) == 0) {
      stop("No rows remain after applying subset_levels.")
    }
  }
  
  if (identical(config$comparison_mode, "condition_vs_condition") && !is.null(config$split_by_col)) {
    available_split_levels <- as.character(config$split_levels %||% unique(dataset_all[[config$split_by_col]]))
    if (nzchar(temp_split_override)) {
      config$split_levels <- intersect(available_split_levels, temp_split_override)
    } else if (nzchar(temp_cluster_override)) {
      config$split_levels <- intersect(available_split_levels, temp_cluster_override)
    }
  }

  if (identical(config$comparison_mode, "one_vs_others")) {
    config$one_vs_others_col <- config$one_vs_others_col %||% config$split_by_col
    config$one_vs_others_levels <- config$one_vs_others_levels %||% config$split_levels
  }
  
  if (identical(config$comparison_mode, "one_vs_others") &&
      !is.null(config$one_vs_others_col) &&
      nzchar(temp_one_vs_others_override)) {
    config$one_vs_others_levels <- intersect(
      as.character(config$one_vs_others_levels %||% unique(dataset_all[[config$one_vs_others_col]])),
      temp_one_vs_others_override
    )
  }
  
  valid_balance_methods <- c("hard_median_downsample", "weighted_sampling", "model_weights")
  valid_patient_split_strategies <- c("shared_patients", "unseen_patients")
  valid_shap_expr_cor_direction_filters <- c("both", "positive", "negative")
  valid_memory_modes <- c("full", "safe")
  valid_score_scaling_methods <- c("median_anchor", "minmax", "winsorized")
  
  if (!config$patient_balance_method %in% valid_balance_methods) {
    stop("patient_balance_method must be one of: ", paste(valid_balance_methods, collapse = ", "))
  }
  if (!config$patient_split_strategy %in% valid_patient_split_strategies) {
    stop(
      "patient_split_strategy must be one of: ",
      paste(valid_patient_split_strategies, collapse = ", ")
    )
  }
  
  if (!config$shap_expr_cor_direction_filter %in% valid_shap_expr_cor_direction_filters) {
    stop("shap_expr_cor_direction_filter must be one of: ", paste(valid_shap_expr_cor_direction_filters, collapse = ", "))
  }
  
  if (!config$memory_mode %in% valid_memory_modes) {
    stop("memory_mode must be one of: ", paste(valid_memory_modes, collapse = ", "))
  }
  if (!config$score_scaling_method %in% valid_score_scaling_methods) {
    stop("score_scaling_method must be one of: ", paste(valid_score_scaling_methods, collapse = ", "))
  }
  
  if (!is.numeric(config$feature_min_pct) ||
      length(config$feature_min_pct) != 1 ||
      is.na(config$feature_min_pct) ||
      config$feature_min_pct < 0 ||
      config$feature_min_pct >= 1) {
    stop("feature_min_pct must be one numeric value >= 0 and < 1.")
  }
  if (!is.numeric(config$feature_min_log2FC) ||
      length(config$feature_min_log2FC) != 1 ||
      is.na(config$feature_min_log2FC) ||
      config$feature_min_log2FC < 0) {
    stop("feature_min_log2FC must be one numeric value >= 0.")
  }
  if (!is.numeric(config$shap_expr_cor_cutoff) ||
      length(config$shap_expr_cor_cutoff) != 1 ||
      is.na(config$shap_expr_cor_cutoff) ||
      config$shap_expr_cor_cutoff < 0 ||
      config$shap_expr_cor_cutoff > 1) {
    stop("shap_expr_cor_cutoff must be one numeric value between 0 and 1.")
  }
  if (!is.numeric(config$score_scaling_probs) ||
      length(config$score_scaling_probs) != 2 ||
      any(!is.finite(config$score_scaling_probs)) ||
      any(config$score_scaling_probs < 0 | config$score_scaling_probs > 1) ||
      config$score_scaling_probs[1] >= config$score_scaling_probs[2]) {
    stop("score_scaling_probs must be two increasing numeric probabilities between 0 and 1.")
  }
  config$analysis_seed <- validate_analysis_seed(config$analysis_seed)
  config$default_prior_scores_weights <- normalize_score_weights(config$default_prior_scores_weights)
  
  config$group_colors <- validate_color_vector(
    colors = config$group_colors,
    arg_name = "group_colors",
    exact_length = 2
  )
  config$continuous_colors <- validate_color_vector(
    colors = config$continuous_colors,
    arg_name = "continuous_colors",
    min_length = 2
  )
  
  if (identical(config$comparison_mode, "condition_vs_condition")) {
    dataset_all$condition <- factor(dataset_all$condition, levels = condition_levels)
  }
  
  if (is.null(config$project_name)) {
    config$project_name <- if (identical(config$comparison_mode, "condition_vs_condition")) {
      paste0(config$what_0, "_", config$what_1)
    } else {
      "one_vs_others"
    }
  }
  
  if (nzchar(temp_run_tag)) {
    config$project_name <- paste0(config$project_name, "_", temp_run_tag)
  }
  
  config$dataset_all <- dataset_all
  config$what_genes <- what_genes
  config$main_dir <- main_dir
  
  config
}


#' Average in-memory SHAP summaries by original cell/metacell unit
#'
#' @param shap_values Numeric SHAP matrix.
#' @param feature_values Numeric matrix/data frame of model-input feature values
#'   aligned to `shap_values`.
#' @param metadata Data frame aligned to the SHAP rows.
#' @param feature_names Feature names to keep.
#' @param unit_col Stable cell/metacell identifier column in `metadata`.
#'
#' @return List with averaged SHAP values, feature values, metadata, and
#'   per-unit evaluation counts.
#' @keywords internal
average_shap_by_unit_in_memory <- function(
    shap_values,
    feature_values,
    metadata,
    feature_names,
    unit_col = ".shapture_unit_id"
) {
  feature_names <- unique(as.character(feature_names))
  shap_values <- as.matrix(shap_values[, feature_names, drop = FALSE])
  feature_values <- as.matrix(feature_values[, feature_names, drop = FALSE])
  metadata <- as.data.frame(metadata, check.names = FALSE)
  
  if (!unit_col %in% colnames(metadata)) {
    stop("Unit-averaged SHAP summaries require metadata column: ", unit_col)
  }
  
  unit_ids <- as.character(metadata[[unit_col]])
  valid_rows <- !is.na(unit_ids) & nzchar(unit_ids)
  if (!all(valid_rows)) {
    shap_values <- shap_values[valid_rows, , drop = FALSE]
    feature_values <- feature_values[valid_rows, , drop = FALSE]
    metadata <- metadata[valid_rows, , drop = FALSE]
    unit_ids <- unit_ids[valid_rows]
  }
  if (length(unit_ids) == 0) {
    stop("No valid unit IDs found for unit-averaged SHAP summaries.")
  }
  
  unit_factor <- factor(unit_ids, levels = unique(unit_ids))
  shap_observed <- !is.na(shap_values)
  x_observed <- !is.na(feature_values)
  shap_values[!shap_observed] <- 0
  feature_values[!x_observed] <- 0
  shap_sum <- rowsum(shap_values, group = unit_factor, reorder = FALSE)
  x_sum <- rowsum(feature_values, group = unit_factor, reorder = FALSE)
  shap_n <- rowsum(shap_observed * 1, group = unit_factor, reorder = FALSE)
  x_n <- rowsum(x_observed * 1, group = unit_factor, reorder = FALSE)
  unit_metadata <- metadata[match(levels(unit_factor), unit_ids), , drop = FALSE]
  
  S <- shap_sum
  X <- x_sum
  S[shap_n > 0] <- shap_sum[shap_n > 0] / shap_n[shap_n > 0]
  X[x_n > 0] <- x_sum[x_n > 0] / x_n[x_n > 0]
  S[shap_n == 0] <- NA_real_
  X[x_n == 0] <- NA_real_
  row.names(S) <- levels(unit_factor)
  row.names(X) <- levels(unit_factor)
  row.names(unit_metadata) <- levels(unit_factor)
  
  list(
    S = S,
    X = as.data.frame(X, check.names = FALSE),
    metadata = unit_metadata,
    unit_counts = rowMeans(shap_n, na.rm = TRUE)
  )
}


#' Row-bind named numeric vectors with absent features as missing
#'
#' @param x_list List of named numeric vectors.
#' @param feature_names Feature names defining the output columns.
#'
#' @return Numeric matrix with one row per list element.
#' @keywords internal
rbind_named_vectors <- function(x_list, feature_names) {
  out <- matrix(
    NA_real_,
    nrow = length(x_list),
    ncol = length(feature_names),
    dimnames = list(seq_along(x_list), feature_names)
  )
  for (i in seq_along(x_list)) {
    x <- x_list[[i]]
    if (is.null(names(x))) {
      stop("Named numeric vectors are required for feature-aware row binding.")
    }
    common <- intersect(names(x), feature_names)
    out[i, common] <- x[common]
  }
  out
}


#' Row-bind matrices with absent feature columns as missing
#'
#' @param x_list List of matrices/data frames.
#' @param feature_names Feature names defining the output columns.
#'
#' @return Numeric matrix with all requested feature columns.
#' @keywords internal
rbind_feature_matrices <- function(x_list, feature_names) {
  n_rows <- vapply(x_list, nrow, integer(1))
  out <- matrix(
    NA_real_,
    nrow = sum(n_rows),
    ncol = length(feature_names),
    dimnames = list(NULL, feature_names)
  )
  row_start <- 1L
  for (i in seq_along(x_list)) {
    x <- as.matrix(x_list[[i]])
    row_end <- row_start + nrow(x) - 1L
    common <- intersect(colnames(x), feature_names)
    out[row_start:row_end, common] <- x[, common, drop = FALSE]
    row_start <- row_end + 1L
  }
  out
}


#' Average safe-mode SHAP chunks by original cell/metacell unit
#'
#' @param chunk_files Vector of chunk-file paths.
#' @param feature_names Feature names to keep.
#' @param unit_col Stable cell/metacell identifier column.
#'
#' @return List with averaged SHAP matrix, averaged feature matrix, metadata,
#'   baseline, and per-unit evaluation counts.
#' @keywords internal
average_shap_by_unit_from_chunks <- function(
    chunk_files,
    feature_names,
    unit_col = ".shapture_unit_id"
) {
  feature_names <- unique(as.character(feature_names))
  metadata_env <- new.env(parent = emptyenv())
  unit_order <- character(0)
  baselines <- numeric(length(chunk_files))
  
  for (i in seq_along(chunk_files)) {
    chunk <- readRDS(chunk_files[i])
    metadata <- as.data.frame(chunk$metadata, check.names = FALSE)
    if (!unit_col %in% colnames(metadata)) {
      stop("Unit-averaged SHAP summaries require metadata column in chunks: ", unit_col)
    }
    baselines[i] <- chunk$baseline
    
    unit_ids <- as.character(metadata[[unit_col]])
    valid_rows <- !is.na(unit_ids) & nzchar(unit_ids)
    if (!all(valid_rows)) {
      metadata <- metadata[valid_rows, , drop = FALSE]
      unit_ids <- unit_ids[valid_rows]
    }
    
    for (unit_id in setdiff(unique(unit_ids), unit_order)) {
        assign(unit_id, metadata[match(unit_id, unit_ids), , drop = FALSE], envir = metadata_env)
        unit_order <- c(unit_order, unit_id)
    }
    
    rm(chunk, metadata, unit_ids)
    gc(verbose = FALSE)
  }
  
  if (length(unit_order) == 0) {
    stop("No valid unit IDs found in safe-mode SHAP chunks.")
  }
  
  S_sum <- matrix(0, nrow = length(unit_order), ncol = length(feature_names), dimnames = list(unit_order, feature_names))
  X_sum <- matrix(0, nrow = length(unit_order), ncol = length(feature_names), dimnames = list(unit_order, feature_names))
  N <- matrix(0, nrow = length(unit_order), ncol = length(feature_names), dimnames = list(unit_order, feature_names))
  
  for (i in seq_along(chunk_files)) {
    chunk <- readRDS(chunk_files[i])
    metadata <- as.data.frame(chunk$metadata, check.names = FALSE)
    unit_ids <- as.character(metadata[[unit_col]])
    valid_rows <- !is.na(unit_ids) & nzchar(unit_ids)
    unit_ids <- unit_ids[valid_rows]
    present_features <- intersect(feature_names, colnames(chunk$shap_values))
    present_features <- intersect(present_features, colnames(chunk$feature_values))
    if (length(unit_ids) == 0 || length(present_features) == 0) {
      next
    }
    
    shap_matrix <- as.matrix(chunk$shap_values[valid_rows, present_features, drop = FALSE])
    x_matrix <- as.matrix(chunk$feature_values[valid_rows, present_features, drop = FALSE])
    group_factor <- factor(unit_ids, levels = unique(unit_ids))
    shap_sum_chunk <- rowsum(shap_matrix, group = group_factor, reorder = FALSE)
    x_sum_chunk <- rowsum(x_matrix, group = group_factor, reorder = FALSE)
    counts_by_unit <- as.numeric(table(group_factor))
    row_idx <- match(rownames(shap_sum_chunk), unit_order)
    col_idx <- match(present_features, feature_names)
    
    S_sum[row_idx, col_idx] <- S_sum[row_idx, col_idx, drop = FALSE] + shap_sum_chunk[, present_features, drop = FALSE]
    X_sum[row_idx, col_idx] <- X_sum[row_idx, col_idx, drop = FALSE] + x_sum_chunk[, present_features, drop = FALSE]
    N[row_idx, col_idx] <- N[row_idx, col_idx, drop = FALSE] +
      matrix(counts_by_unit, nrow = length(row_idx), ncol = length(col_idx))
    
    rm(chunk, metadata, unit_ids, shap_matrix, x_matrix, shap_sum_chunk, x_sum_chunk)
    gc(verbose = FALSE)
  }
  
  S <- S_sum
  X <- X_sum
  observed <- N > 0
  S[observed] <- S_sum[observed] / N[observed]
  X[observed] <- X_sum[observed] / N[observed]
  S[!observed] <- NA_real_
  X[!observed] <- NA_real_
  metadata <- do.call(rbind, lapply(unit_order, function(unit_id) {
    get(unit_id, envir = metadata_env, inherits = FALSE)
  }))
  row.names(metadata) <- unit_order
  
  list(
    S = S,
    X = as.data.frame(X, check.names = FALSE),
    metadata = metadata,
    baseline = mean(baselines, na.rm = TRUE),
    unit_counts = rowMeans(N, na.rm = TRUE)
  )
}


#' Accumulate patient-level SHAP means from iteration chunk files
#'
#' @param chunk_files Vector of chunk-file paths.
#' @param feature_names Feature names to keep.
#' @param sample_col Sample identifier column in chunk metadata.
#' @param condition_col Condition column in chunk metadata.
#'
#' @return List with `patient_means` and `patient_units`.
#' @keywords internal
build_patient_means_from_chunks <- function(
    chunk_files,
    feature_names,
    sample_col = "sample_id",
    condition_col = "condition"
) {
  unit_average <- average_shap_by_unit_from_chunks(
    chunk_files = chunk_files,
    feature_names = feature_names
  )
  metadata <- unit_average$metadata
  shap_matrix <- unit_average$S
  
  valid_rows <- !is.na(metadata[[sample_col]]) & !is.na(metadata[[condition_col]])
  metadata <- metadata[valid_rows, , drop = FALSE]
  shap_matrix <- shap_matrix[valid_rows, , drop = FALSE]
  
  patient_keys <- paste(metadata[[sample_col]], metadata[[condition_col]], sep = "___")
  patient_unit_keys <- unique(patient_keys)
  if (length(patient_unit_keys) == 0) {
    stop("No patient-level SHAP data were found after unit averaging.")
  }
  
  patient_means <- t(vapply(
    patient_unit_keys,
    function(key) {
      idx <- patient_keys == key
      colMeans(shap_matrix[idx, , drop = FALSE], na.rm = TRUE)
    },
    numeric(length(feature_names))
  ))
  rownames(patient_means) <- patient_unit_keys
  colnames(patient_means) <- feature_names
  
  patient_unit_split <- strsplit(patient_unit_keys, "___", fixed = TRUE)
  patient_units <- data.frame(
    sample_id = vapply(patient_unit_split, `[`, character(1), 1),
    condition = vapply(patient_unit_split, `[`, character(1), 2),
    stringsAsFactors = FALSE
  )
  
  list(patient_means = patient_means, patient_units = patient_units)
}


#' Compute SHAP-value correlations from iteration chunk files
#'
#' @param chunk_files Vector of chunk-file paths.
#' @param feature_names Feature names to process. Chunk feature values should
#'   contain the preprocessed model-input values used for prediction.
#' @param block_size Number of features to process per block.
#'
#' @return Named numeric vector of Spearman correlations.
#' @keywords internal
compute_shap_expr_cor_from_chunks <- function(chunk_files, feature_names, block_size = 32L) {
  unit_average <- average_shap_by_unit_from_chunks(
    chunk_files = chunk_files,
    feature_names = feature_names
  )
  raw_mat <- as.matrix(unit_average$X[, feature_names, drop = FALSE])
  shap_mat <- as.matrix(unit_average$S[, feature_names, drop = FALSE])
  block_size <- max(1L, as.integer(block_size))
  shap_expr_cor <- numeric(length(feature_names))
  names(shap_expr_cor) <- feature_names
  
  block_starts <- seq.int(1L, length(feature_names), by = block_size)
  for (block_start in block_starts) {
    block_end <- min(length(feature_names), block_start + block_size - 1L)
    feature_block <- feature_names[block_start:block_end]
    raw_block <- raw_mat[, feature_block, drop = FALSE]
    shap_block <- shap_mat[, feature_block, drop = FALSE]
    
    block_cor <- vapply(seq_along(feature_block), function(j) {
      suppressWarnings(
        stats::cor(
          raw_block[, j],
          shap_block[, j],
          method = "spearman",
          use = "pairwise.complete.obs"
        )
      )
    }, numeric(1))
    block_cor[!is.finite(block_cor)] <- 0
    shap_expr_cor[feature_block] <- block_cor
    
    rm(raw_block, shap_block, block_cor, feature_block)
    gc(verbose = FALSE)
  }
  
  shap_expr_cor
}


#' Collect pooled SHAP data for selected features from iteration chunks
#'
#' @param chunk_files Vector of chunk-file paths.
#' @param feature_names Feature names to collect.
#'
#' @return List with pooled `S`, `X`, `metadata`, and `baseline`.
#' @keywords internal
collect_shap_plot_data_from_chunks <- function(chunk_files, feature_names) {
  unit_average <- average_shap_by_unit_from_chunks(
    chunk_files = chunk_files,
    feature_names = feature_names
  )
  list(
    S = unit_average$S,
    X = unit_average$X,
    metadata = unit_average$metadata,
    baseline = unit_average$baseline,
    unit_counts = unit_average$unit_counts
  )
}


#' Build condition-level mean SHAP and feature summaries from chunk files
#'
#' @param chunk_files Vector of chunk-file paths.
#' @param feature_names Feature names to aggregate.
#' @param condition_col Condition column in chunk metadata.
#'
#' @return Named list keyed by condition value.
#' @keywords internal
build_condition_shap_aggregates_from_chunks <- function(
    chunk_files,
    feature_names,
    condition_col = "condition"
) {
  unit_average <- average_shap_by_unit_from_chunks(
    chunk_files = chunk_files,
    feature_names = feature_names
  )
  metadata <- unit_average$metadata
  shap_values <- as.matrix(unit_average$S[, feature_names, drop = FALSE])
  feature_values <- as.matrix(unit_average$X[, feature_names, drop = FALSE])
  
  condition_keys <- sort(unique(as.character(metadata[[condition_col]])))
  aggregates <- setNames(vector("list", length(condition_keys)), condition_keys)
  for (key in condition_keys) {
    idx <- as.character(metadata[[condition_col]]) == key
    aggregates[[key]] <- list(
      S = colMeans(shap_values[idx, , drop = FALSE], na.rm = TRUE),
      X = colMeans(feature_values[idx, , drop = FALSE], na.rm = TRUE)
    )
  }
  attr(aggregates, "baseline") <- unit_average$baseline
  aggregates
}


#' Prepare split metadata for SHAPTURE analyses
#'
#' @param dataset_all Input data frame.
#' @param split_by_col Optional subset column. In `one_vs_others` mode, this is
#'   also the column whose levels define the positive class.
#'
#' @return List with a possibly augmented dataset and the effective split column.
#' @keywords internal
prepare_split_context <- function(dataset_all, split_by_col) {
  internal_all_data_col <- ".shapture_all_data"
  
  if (is.null(split_by_col)) {
    dataset_all[[internal_all_data_col]] <- "all_data"
    return(list(
      dataset = dataset_all,
      effective_split_by_col = internal_all_data_col,
      effective_split_levels = "all_data"
    ))
  }
  
  if (!split_by_col %in% colnames(dataset_all)) {
    stop("split_by_col is not present in dataset_all: ", split_by_col)
  }
  
  split_levels <- as.character(unique(dataset_all[[split_by_col]]))
  split_levels <- split_levels[!is.na(split_levels) & nzchar(split_levels)]
  
  list(
    dataset = dataset_all,
    effective_split_by_col = split_by_col,
    effective_split_levels = split_levels
  )
}


#' Build the analysis task table for SHAPTURE
#'
#' @param dataset_all Input data frame.
#' @param comparison_mode One of `"condition_vs_condition"` or
#'   `"one_vs_others"`.
#' @param split_by_col Optional subset column.
#' @param split_levels Optional selected split levels.
#' @param one_vs_others_col Optional one-vs-others grouping column.
#' @param one_vs_others_levels Optional target levels for one-vs-others.
#' @param what_0 Negative class label.
#' @param what_1 Positive class label.
#' @param one_vs_others_negative_label Label used for the pooled negative class
#'   in one-vs-others mode.
#'
#' @return A list containing the task table, effective split column, and
#'   prepared input data.
#' @keywords internal
build_shapture_tasks <- function(
    dataset_all,
    comparison_mode,
    split_by_col,
    split_levels,
    one_vs_others_col,
    one_vs_others_levels,
    what_0,
    what_1,
    one_vs_others_negative_label
) {
  comparison_mode <- normalize_comparison_mode(comparison_mode)
  comparison_mode <- match.arg(comparison_mode, c("condition_vs_condition", "one_vs_others"))
  
  if (comparison_mode == "condition_vs_condition") {
    split_context <- prepare_split_context(dataset_all = dataset_all, split_by_col = split_by_col)
    effective_split_levels <- split_context$effective_split_levels
    
    if (!is.null(split_levels)) {
      effective_split_levels <- intersect(as.character(effective_split_levels), as.character(split_levels))
    }
    
    if (length(effective_split_levels) == 0) {
      stop("No split levels remain to analyze.")
    }
    
    task_table <- data.frame(
      temp_split = as.character(effective_split_levels),
      task_what_0 = what_0,
      task_what_1 = what_1,
      stringsAsFactors = FALSE
    )
    
    return(list(
      dataset = split_context$dataset,
      effective_split_by_col = split_context$effective_split_by_col,
      task_metacell_split_col = split_context$effective_split_by_col,
      task_table = task_table
    ))
  }
  
  one_vs_others_col <- one_vs_others_col %||% split_by_col
  one_vs_others_levels <- one_vs_others_levels %||% split_levels

  if (is.null(one_vs_others_col) || !one_vs_others_col %in% colnames(dataset_all)) {
    stop(
      "For one_vs_others mode, split_by_col must be set to the column that defines the target groups."
    )
  }
  
  target_levels <- as.character(unique(dataset_all[[one_vs_others_col]]))
  target_levels <- target_levels[!is.na(target_levels) & nzchar(target_levels)]
  
  if (!is.null(one_vs_others_levels)) {
    target_levels <- intersect(target_levels, as.character(one_vs_others_levels))
  }
  
  if (length(target_levels) == 0) {
    stop("No one-vs-others target levels remain to analyze.")
  }
  
  task_table <- data.frame(
    temp_split = target_levels,
    task_what_0 = one_vs_others_negative_label,
    task_what_1 = target_levels,
    stringsAsFactors = FALSE
  )
  
  list(
    dataset = dataset_all,
    effective_split_by_col = NULL,
    task_metacell_split_col = one_vs_others_col,
    task_table = task_table
  )
}


#' Materialize one SHAPTURE analysis task
#'
#' @param dataset Analysis data set.
#' @param comparison_mode Analysis comparison mode.
#' @param temp_split Current split/target level.
#' @param split_by_col Effective split column for condition-vs-condition mode.
#' @param one_vs_others_col Column used in one-vs-others mode.
#' @param task_what_0 Negative class label for this task.
#' @param task_what_1 Positive class label for this task.
#'
#' @return Task-specific data frame with a binary condition column.
#' @keywords internal
materialize_task_dataset <- function(
    dataset,
    comparison_mode,
    temp_split,
    split_by_col,
    one_vs_others_col,
    task_what_0,
    task_what_1
) {
  comparison_mode <- normalize_comparison_mode(comparison_mode)
  comparison_mode <- match.arg(comparison_mode, c("condition_vs_condition", "one_vs_others"))
  
  if (comparison_mode == "condition_vs_condition") {
    if (is.null(split_by_col) || identical(split_by_col, ".shapture_all_data")) {
      dataset_task <- dataset
    } else {
      dataset_task <- dataset[dataset[[split_by_col]] == temp_split, , drop = FALSE]
    }
    dataset_task$condition <- factor(dataset_task$condition, levels = c(task_what_0, task_what_1))
    return(dataset_task)
  }
  
  dataset_task <- dataset
  dataset_task$condition <- ifelse(dataset_task[[one_vs_others_col]] == temp_split, task_what_1, task_what_0)
  dataset_task$condition <- factor(dataset_task$condition, levels = c(task_what_0, task_what_1))
  dataset_task
}


#' Filter task features by detection percentage and expression log2FC
#'
#' @param dataset Task-level data frame.
#' @param feature_cols Candidate feature columns.
#' @param condition_col Binary condition column.
#' @param condition_levels Two condition labels in negative/positive order.
#' @param feature_min_pct Minimum fraction detected in at least one group.
#' @param feature_min_log2FC Minimum absolute expression log2FC between the two
#'   task groups. Set to `0` to disable this filter.
#' @param expression_threshold Values greater than this are treated as detected.
#' @param chunk_size Number of features processed per chunk.
#'
#' @return A list with kept feature names and a summary data frame.
#' @keywords internal
filter_features_by_min_pct <- function(
    dataset,
    feature_cols,
    condition_col = "condition",
    condition_levels,
    feature_min_pct = 0.25,
    feature_min_log2FC = 0.25,
    expression_threshold = 0,
    chunk_size = 1000
) {
  feature_cols <- intersect(feature_cols, colnames(dataset))
  if (length(feature_cols) == 0) {
    return(list(
      features = character(0),
      summary = data.frame()
    ))
  }
  
  condition_values <- as.character(dataset[[condition_col]])
  idx_0 <- condition_values == as.character(condition_levels[1])
  idx_1 <- condition_values == as.character(condition_levels[2])
  
  pct_0 <- stats::setNames(rep(NA_real_, length(feature_cols)), feature_cols)
  pct_1 <- stats::setNames(rep(NA_real_, length(feature_cols)), feature_cols)
  mean_0 <- stats::setNames(rep(NA_real_, length(feature_cols)), feature_cols)
  mean_1 <- stats::setNames(rep(NA_real_, length(feature_cols)), feature_cols)
  chunk_starts <- seq.int(1L, length(feature_cols), by = chunk_size)
  
  for (chunk_start in chunk_starts) {
    chunk_end <- min(chunk_start + chunk_size - 1L, length(feature_cols))
    chunk_features <- feature_cols[chunk_start:chunk_end]
    
    if (any(idx_0)) {
      mat_0 <- as.matrix(dataset[idx_0, chunk_features, drop = FALSE])
      pct_0[chunk_features] <- colMeans(
        mat_0 > expression_threshold,
        na.rm = TRUE
      )
      mean_0[chunk_features] <- colMeans(mat_0, na.rm = TRUE)
    }
    if (any(idx_1)) {
      mat_1 <- as.matrix(dataset[idx_1, chunk_features, drop = FALSE])
      pct_1[chunk_features] <- colMeans(
        mat_1 > expression_threshold,
        na.rm = TRUE
      )
      mean_1[chunk_features] <- colMeans(mat_1, na.rm = TRUE)
    }
  }
  
  max_pct <- pmax(
    ifelse(is.na(pct_0), 0, pct_0),
    ifelse(is.na(pct_1), 0, pct_1)
  )
  
  if (!is.finite(feature_min_pct) || feature_min_pct <= 0) {
    kept_by_min_pct <- rep(TRUE, length(feature_cols))
  } else {
    kept_by_min_pct <- max_pct > feature_min_pct
  }
  names(kept_by_min_pct) <- feature_cols
  
  log_pseudocount <- compute_dynamic_pseudocount(
    x = c(mean_0, mean_1),
    fraction = 0.01,
    fallback = 1e-6
  )
  expression_log2FC <- log2((mean_1 + log_pseudocount) / (mean_0 + log_pseudocount))
  expression_log2FC_abs <- abs(expression_log2FC)
  
  if (!is.finite(feature_min_log2FC) || feature_min_log2FC <= 0) {
    kept_by_min_log2FC <- rep(TRUE, length(feature_cols))
  } else {
    kept_by_min_log2FC <- expression_log2FC_abs >= feature_min_log2FC
  }
  names(kept_by_min_log2FC) <- feature_cols
  
  kept <- kept_by_min_pct & kept_by_min_log2FC
  names(kept) <- feature_cols
  
  summary <- data.frame(
    feature = feature_cols,
    pct_condition0 = unname(pct_0[feature_cols]),
    pct_condition1 = unname(pct_1[feature_cols]),
    max_pct = unname(max_pct[feature_cols]),
    feature_min_pct = feature_min_pct,
    kept_by_min_pct = unname(kept_by_min_pct[feature_cols]),
    expression_log2FC_singlecell = unname(expression_log2FC[feature_cols]),
    expression_log2FC_abs_singlecell = unname(expression_log2FC_abs[feature_cols]),
    feature_min_log2FC = feature_min_log2FC,
    kept_by_min_log2FC = unname(kept_by_min_log2FC[feature_cols]),
    kept_by_prefilter = unname(kept[feature_cols]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  list(
    features = feature_cols[kept],
    summary = summary
  )
}


#' Run one SHAPTURE analysis task
#'
#' @param dataset_all Original input data.
#' @param dataset_all_analysis Analysis data after optional metacell generation.
#' @param what_genes Feature list.
#' @param temp_split Current split or one-vs-others target level.
#' @param task_what_0 Negative class label for this task.
#' @param task_what_1 Positive class label for this task.
#' @param main_dir Output root directory.
#' @param project_name Project prefix.
#' @param comparison_mode Comparison mode.
#' @param effective_split_by_col Effective split column for condition-vs-condition mode.
#' @param metacell_split_col Column used during metacell construction.
#' @param one_vs_others_col Column defining the positive class in
#'   one-vs-others mode.
#' @param metacell_summary Optional metacell summary data frame.
#' @param continuous_colors Continuous color palette.
#' @param group_colors Two unnamed colors used for left/right group displays.
#' @param train_ratio Training fraction.
#' @param test_ratio Test fraction.
#' @param iteration_no Number of repeated iterations.
#' @param gene_fraction_per_iteration Fraction of task features included in
#'   each coverage-aware model iteration.
#' @param cross_patient_bootstraps Number of patient-block bootstrap replicates
#'   used for Cross-patient Support.
#' @param n_cores Requested maximum workers for model fitting and Cross-patient
#'   bootstrapping. Bootstrap safeguards may reduce this value or run serially.
#' @param analysis_seed Positive integer reproducibility root for all stochastic
#'   task stages.
#' @param metacell_seed Seed used for grouped metacell construction, or `NA`
#'   when single cells are analyzed.
#' @param top_feature_label Number of labeled features.
#' @param top_feature_show Number of top features to show and save for
#'   post-analysis plotting.
#' @param if_adaptive_weights Whether to use adaptive weights.
#' @param feature_min_pct Minimum fraction of cells with feature value greater
#'   than zero in at least one task group. Set to `0` to disable this filter.
#' @param feature_min_log2FC Minimum absolute expression log2FC in original
#'   single-cell task data. Set to `0` to disable this filter.
#' @param shap_expr_cor_cutoff Minimum absolute value-SHAP correlation needed
#'   to assign a model-expression pattern label.
#' @param score_scaling_method Backward-compatible configuration field. The
#'   evidence scores use the fixed median-anchor and native bounded strategy;
#'   legacy `"winsorized"` and `"minmax"` values remain accepted so older
#'   scripts continue to run.
#' @param score_scaling_probs Quantiles used only when winsorizing component
#'   values for optional adaptive-weight estimation.
#' @param default_prior_scores_weights Named prior weights for final score
#'   components.
#' @param feature_min_log2FC Minimum absolute expression log2FC in original
#'   single-cell task data. Set to `0` to disable this filter.
#' @param preprocessing_method Preprocessing method.
#' @param model_type Model backend.
#' @param analysis_unit_type Either `"single_cell"` or `"metacell"`.
#' @param memory_mode Either `"full"` or `"safe"`.
#' @param patient_balance_method Patient balancing method.
#' @param patient_split_strategy Either `"shared_patients"` for the existing
#'   cell/metacell-level split or `"unseen_patients"` for globally disjoint
#'   patient-level train/test folds with automatic task-level fallback.
#' @param shap_expr_cor_direction_filter One of `"both"`, `"positive"`,
#'   `"negative"`. Filters retained features by the sign of `shap_expr_cor`.
#' @param default_prior_scores_weights Named prior weights for final score
#'   components.
#' @param resolved_config_for_docs Named list of resolved run configuration
#'   values to save into each task's docs folder.
#'
#' @return A list summarizing the task output.
#' @keywords internal
run_shapture_task <- function(
    dataset_all,
    dataset_all_analysis,
    what_genes,
    temp_split,
    task_what_0,
    task_what_1,
    main_dir,
    project_name,
    comparison_mode,
    effective_split_by_col,
    metacell_split_col,
    one_vs_others_col,
    metacell_summary,
    continuous_colors,
    group_colors,
    train_ratio,
    test_ratio,
    iteration_no,
    gene_fraction_per_iteration,
    cross_patient_bootstraps,
    n_cores,
    analysis_seed,
    metacell_seed,
    top_feature_label,
    top_feature_show,
    if_adaptive_weights,
    feature_min_pct,
    feature_min_log2FC,
    shap_expr_cor_cutoff,
    score_scaling_method,
    score_scaling_probs,
    preprocessing_method,
    model_type,
    analysis_unit_type,
    memory_mode,
    patient_balance_method,
    patient_split_strategy,
    shap_expr_cor_direction_filter,
    default_prior_scores_weights,
    resolved_config_for_docs
) {
  memory_mode <- match.arg(memory_mode, c("full", "safe"))
  patient_split_strategy <- match.arg(
    patient_split_strategy,
    c("shared_patients", "unseen_patients")
  )
  analysis_seed <- validate_analysis_seed(analysis_seed)
  seed_plan <- build_shapture_task_seed_plan(
    analysis_seed = analysis_seed,
    comparison_mode = comparison_mode,
    split_by_col = effective_split_by_col,
    temp_split = temp_split,
    task_what_0 = task_what_0,
    task_what_1 = task_what_1,
    iteration_no = iteration_no
  )
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  
  setwd(main_dir)
  
  dataset_task <- materialize_task_dataset(
    dataset = dataset_all_analysis,
    comparison_mode = comparison_mode,
    temp_split = temp_split,
    split_by_col = effective_split_by_col,
    one_vs_others_col = one_vs_others_col,
    task_what_0 = task_what_0,
    task_what_1 = task_what_1
  )
  
  dataset_task_singlecell <- materialize_task_dataset(
    dataset = dataset_all,
    comparison_mode = comparison_mode,
    temp_split = temp_split,
    split_by_col = effective_split_by_col,
    one_vs_others_col = one_vs_others_col,
    task_what_0 = task_what_0,
    task_what_1 = task_what_1
  )
  n_rows_task_model_initial <- nrow(dataset_task)
  n_rows_task_singlecell <- nrow(dataset_task_singlecell)
  dataset_task$.shapture_unit_id <- make.unique(as.character(row.names(dataset_task)))
  dataset_task_singlecell$.shapture_unit_id <- make.unique(as.character(row.names(dataset_task_singlecell)))
  
  sub_folder <- paste0(project_name, "_", temp_split)
  dir.create(sub_folder, showWarnings = FALSE)
  setwd(sub_folder)
  dir.create("docs", showWarnings = FALSE)
  dir.create("SHAP_only", showWarnings = FALSE)
  dir.create("Integrated_adaptive_scoring_results", showWarnings = FALSE)
  
  task_feature_cols <- intersect(what_genes, colnames(dataset_task))
  n_candidate_features <- length(task_feature_cols)
  feature_filter_res <- filter_features_by_min_pct(
    dataset = dataset_task_singlecell,
    feature_cols = task_feature_cols,
    condition_col = "condition",
    condition_levels = c(task_what_0, task_what_1),
    feature_min_pct = feature_min_pct,
    feature_min_log2FC = feature_min_log2FC
  )
  what_genes <- feature_filter_res$features
  n_features_after_prefilter <- length(what_genes)
  feature_prefilter_summary <- feature_filter_res$summary
  if (nrow(feature_prefilter_summary) > 0) {
    write.csv(feature_prefilter_summary, "docs/feature_prefilter_summary.csv", row.names = FALSE)
  }
  
  if (length(what_genes) == 0) {
    cat(paste0("no features passed feature_min_pct/feature_min_log2FC for ", temp_split, " skipping the analysis", "\n"))
    return(list(
      temp_split = temp_split,
      analyzed = FALSE,
      each_condition_sample_no = NA_real_
    ))
  }
  
  metadata_cols_for_split <- ".shapture_unit_id"
  if ("sample_id" %in% colnames(dataset_task)) {
    metadata_cols_for_split <- c(metadata_cols_for_split, "sample_id")
  }
  sampling_prob_col <- NULL
  balance_col <- NULL

  if ("sample_id" %in% colnames(dataset_task)) {
    dataset_task <- build_balance_groups(
      dataset = dataset_task,
      comparison_mode = comparison_mode,
      sample_col = "sample_id",
      condition_col = "condition",
      split_by_col = one_vs_others_col,
      positive_label = task_what_1,
      balance_col = "balance_group"
    )
    balance_col <- "balance_group"
  }
  
  if (patient_balance_method == "hard_median_downsample") {
    if (!is.null(balance_col)) {
      dataset_task <- downsample_by_condition_and_group(
        dataset = dataset_task,
        condition_col = "condition",
        balance_col = balance_col,
        seed = seed_plan$balance_seed
      )
    }
  } else if (patient_balance_method == "weighted_sampling") {
    if (!is.null(balance_col)) {
      dataset_task <- add_inverse_group_probabilities(
        dataset = dataset_task,
        group_col = balance_col,
        condition_col = "condition",
        prob_col = "sampling_prob"
      )
      sampling_prob_col <- "sampling_prob"
    }
  } else if (patient_balance_method == "model_weights") {
    metadata_cols_for_split <- unique(c(".shapture_unit_id", "sample_id", balance_col))
  }
  metadata_cols_for_split <- intersect(metadata_cols_for_split, colnames(dataset_task))
  
  dataset_task$condition <- ifelse(dataset_task$condition == task_what_0, 0, 1)
  
  keep_cols <- c("condition", metadata_cols_for_split, sampling_prob_col, what_genes)
  keep_cols <- unique(keep_cols[!is.na(keep_cols) & nzchar(keep_cols)])
  dataset_task <- dataset_task[, keep_cols, drop = FALSE]
  
  feature_cols_in_task <- setdiff(colnames(dataset_task), c("condition", metadata_cols_for_split, sampling_prob_col))
  nonzero_feature_cols <- feature_cols_in_task[colSums(dataset_task[, feature_cols_in_task, drop = FALSE]) > 0]
  if (length(nonzero_feature_cols) == 0) {
    cat(paste0("no non-zero features for ", temp_split, " skipping the analysis", "\n"))
    return(list(
      temp_split = temp_split,
      analyzed = FALSE,
      each_condition_sample_no = NA_real_
    ))
  }
  dataset_task <- dataset_task[, c("condition", metadata_cols_for_split, sampling_prob_col, nonzero_feature_cols), drop = FALSE]
  
  min_n_cells <- min(table(dataset_task$condition))
  
  if (analysis_unit_type == "metacell") {
    min_required_per_condition <- 50
    each.condition_total_sample_No <- min(floor(0.5 * min_n_cells), 5000)
  } else {
    min_required_per_condition <- 100
    each.condition_total_sample_No <- min(floor(0.5 * min_n_cells), 5000)
  }
  
  if (min_n_cells < min_required_per_condition) {
    cat(paste0("very low number of rows for ", temp_split, " skipping the analysis", "\n"))
    return(list(
      temp_split = temp_split,
      analyzed = FALSE,
      each_condition_sample_no = each.condition_total_sample_No
    ))
  }
  
  top_feature_show_local <- min(top_feature_show, length(what_genes))
  seed_manifest <- seed_plan$iteration_seeds
  random_seeds <- seed_manifest$cell_sampling_seed
  feature_selection_seeds <- seed_manifest$feature_selection_seed
  model_seeds <- seed_manifest$model_seed
  feature_bagging_fraction <- gene_fraction_per_iteration
  patient_split_strategy_requested <- patient_split_strategy
  patient_split_strategy_used <- "shared_patients"
  patient_split_fallback_reason <- NA_character_
  patient_coverage_summary <- NULL
  patient_feature_coverage_summary <- NULL
  eligible_patient_condition_keys <- NULL
  min_unseen_units_per_patient <- if (analysis_unit_type == "single_cell") 20L else 5L
  min_unseen_train_units_per_condition <- if (analysis_unit_type == "single_cell") 50L else 20L

  if (patient_split_strategy_requested == "unseen_patients") {
    unseen_plan_res <- build_unseen_patient_split_plan(
      dataset = dataset_task,
      random_seeds = random_seeds,
      train_ratio = train_ratio,
      test_ratio = test_ratio,
      each_condition_total_sample_no = each.condition_total_sample_No,
      condition_col = "condition",
      sample_col = "sample_id",
      unit_col = ".shapture_unit_id",
      sampling_prob_col = sampling_prob_col,
      min_patients_per_condition = 5L,
      patient_train_fraction = 0.60,
      min_units_per_patient = min_unseen_units_per_patient,
      min_train_units_per_condition = min_unseen_train_units_per_condition
    )
    patient_coverage_summary <- unseen_plan_res$patient_coverage_summary
    if (isTRUE(unseen_plan_res$success)) {
      split_plan_res <- unseen_plan_res
      patient_split_strategy_used <- "unseen_patients"
      eligible_rows <- unseen_plan_res$patient_coverage_summary$eligible
      eligible_patient_condition_keys <- paste(
        unseen_plan_res$patient_coverage_summary$condition[eligible_rows],
        unseen_plan_res$patient_coverage_summary$patient_id[eligible_rows],
        sep = "::"
      )
    } else {
      patient_split_fallback_reason <- unseen_plan_res$reason
      warning(
        "Task '", temp_split, "' requested unseen_patients but will use ",
        "shared_patients. Reason: ", patient_split_fallback_reason,
        call. = FALSE
      )
    }
  }

  if (patient_split_strategy_used == "shared_patients") {
    split_plan_res <- build_coverage_aware_split_plan(
      dataset = dataset_task,
      random_seeds = random_seeds,
      train_ratio = train_ratio,
      test_ratio = test_ratio,
      each_condition_total_sample_no = each.condition_total_sample_No,
      condition_col = "condition",
      unit_col = ".shapture_unit_id",
      sampling_prob_col = sampling_prob_col
    )
  }
  split_plans <- split_plan_res$plans
  sampling_coverage_summary <- split_plan_res$coverage_summary
  if (patient_split_strategy_used == "unseen_patients") {
    feature_plan_res <- build_patient_coverage_aware_feature_plan(
      feature_names = nonzero_feature_cols,
      random_seeds = feature_selection_seeds,
      test_patient_condition_keys = split_plan_res$test_patient_condition_keys,
      eligible_patient_condition_keys = eligible_patient_condition_keys,
      fraction = feature_bagging_fraction
    )
    patient_feature_coverage_summary <- feature_plan_res$patient_feature_coverage
    feature_plan_res$patient_feature_counts <- NULL
  } else {
    feature_plan_res <- build_coverage_aware_feature_plan(
      feature_names = nonzero_feature_cols,
      random_seeds = feature_selection_seeds,
      fraction = feature_bagging_fraction
    )
  }
  feature_plans <- feature_plan_res$plans
  feature_bagging_summary <- feature_plan_res$summary
  if (!is.null(patient_coverage_summary)) {
    patient_coverage_summary$condition_label <- ifelse(
      as.character(patient_coverage_summary$condition) == "0",
      task_what_0,
      task_what_1
    )
  }
  
  task_config_for_docs <- c(
    resolved_config_for_docs,
    list(
      temp_split = temp_split,
      task_what_0 = task_what_0,
      task_what_1 = task_what_1,
      task_seed = seed_plan$task_seed,
      metacell_seed = metacell_seed,
      cross_patient_bootstrap_seed = seed_plan$cross_patient_bootstrap_seed,
      patient_split_strategy_requested = patient_split_strategy_requested,
      patient_split_strategy_used = patient_split_strategy_used,
      patient_split_fallback_reason = patient_split_fallback_reason
    )
  )
  write_shapture_config_csv(task_config_for_docs, "docs/shapture_config_used.csv")
  temp_chunk_dir <- file.path(getwd(), "docs", ".tmp_iteration_chunks")
  if (memory_mode == "safe") {
    dir.create(temp_chunk_dir, showWarnings = FALSE, recursive = TRUE)
    temp_chunk_dir <- normalizePath(temp_chunk_dir, winslash = "/", mustWork = TRUE)
    on.exit(unlink(temp_chunk_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }
  
  cat(paste0("Analyzing ", temp_split, "\n"))
  
  write.csv(data.frame(included.genes = what_genes), "docs/included.genes.csv", row.names = FALSE)
  write.csv(data.frame(random_seeds = random_seeds), "docs/random_seeds.csv", row.names = FALSE)
  reproducibility_seeds <- data.frame(
    analysis_seed = analysis_seed,
    task_seed = seed_plan$task_seed,
    balance_seed = seed_plan$balance_seed,
    metacell_seed = metacell_seed,
    cross_patient_bootstrap_seed = seed_plan$cross_patient_bootstrap_seed,
    seed_manifest,
    check.names = FALSE
  )
  write.csv(
    reproducibility_seeds,
    "docs/reproducibility_seeds.csv",
    row.names = FALSE
  )
  write.csv(sampling_coverage_summary, "docs/sampling_coverage_summary.csv", row.names = FALSE)
  write.csv(feature_bagging_summary, "docs/feature_bagging_summary.csv", row.names = FALSE)
  patient_split_summary <- data.frame(
    patient_split_strategy_requested = patient_split_strategy_requested,
    patient_split_strategy_used = patient_split_strategy_used,
    fallback_reason = patient_split_fallback_reason,
    patient_train_fraction = if (patient_split_strategy_requested == "unseen_patients") 0.60 else NA_real_,
    patient_test_fraction = if (patient_split_strategy_requested == "unseen_patients") 0.40 else NA_real_,
    min_patients_per_condition = if (patient_split_strategy_requested == "unseen_patients") 5L else NA_integer_,
    min_units_per_patient = if (patient_split_strategy_requested == "unseen_patients") min_unseen_units_per_patient else NA_integer_,
    min_train_units_per_condition = if (patient_split_strategy_requested == "unseen_patients") min_unseen_train_units_per_condition else NA_integer_,
    train_patients_condition0 = if (patient_split_strategy_used == "unseen_patients") split_plan_res$n_train_patients[[1]] else NA_integer_,
    train_patients_condition1 = if (patient_split_strategy_used == "unseen_patients") split_plan_res$n_train_patients[[2]] else NA_integer_,
    test_patients_condition0 = if (patient_split_strategy_used == "unseen_patients") split_plan_res$n_test_patients[[1]] else NA_integer_,
    test_patients_condition1 = if (patient_split_strategy_used == "unseen_patients") split_plan_res$n_test_patients[[2]] else NA_integer_,
    train_rows_per_condition = if (patient_split_strategy_used == "unseen_patients") split_plan_res$train_rows_per_condition else NA_integer_,
    test_rows_per_condition = if (patient_split_strategy_used == "unseen_patients") split_plan_res$test_rows_per_condition else NA_integer_,
    stringsAsFactors = FALSE
  )
  write.csv(patient_split_summary, "docs/patient_split_summary.csv", row.names = FALSE)
  if (!is.null(patient_coverage_summary)) {
    write.csv(
      patient_coverage_summary,
      "docs/patient_sampling_coverage_summary.csv",
      row.names = FALSE
    )
  }
  if (!is.null(patient_feature_coverage_summary)) {
    write.csv(
      patient_feature_coverage_summary,
      "docs/patient_feature_coverage_summary.csv",
      row.names = FALSE
    )
  }
  if (!is.null(metacell_summary)) {
    write.csv(
      metacell_summary[metacell_summary$split_level == temp_split, , drop = FALSE],
      "docs/metacell_summary.csv",
      row.names = FALSE
    )
  }
  
  worker_exports <- c(
    "prepare_train_test_data",
    "compute_model_weights_by_group",
    "compute_patient_balanced_shap_effect"
  )
  worker_export_env <- environment()
  missing_worker_exports <- worker_exports[!vapply(
    worker_exports,
    function(x) exists(x, envir = worker_export_env, inherits = TRUE),
    logical(1)
  )]
  if (length(missing_worker_exports) > 0) {
    stop(
      "Missing foreach exports in run_shapture_task(): ",
      paste(missing_worker_exports, collapse = ", ")
    )
  }
  
  file.create("docs/model performance.txt")

  model_cluster <- NULL
  model_workers_stopped <- TRUE
  on.exit({
    safe_stop_shapture_cluster(model_cluster)
    foreach::registerDoSEQ()
  }, add = TRUE)
  if (n_cores > 1L) {
    model_cluster <- parallel::makePSOCKcluster(n_cores)
    doParallel::registerDoParallel(model_cluster)
  } else {
    foreach::registerDoSEQ()
  }
  
  iter_results <- foreach::foreach(
    temp_iter = seq_len(iteration_no),
    .packages = c("xgboost", "lightgbm", "shapviz", "ggplot2", "pROC", "caret", "dplyr"),
    .export = worker_exports
  ) %dopar% {
    worker_stage <- "initializing worker iteration"
    tryCatch({
    temp_random_seed <- random_seeds[temp_iter]
    model_seed <- model_seeds[temp_iter]
    
    worker_stage <- "preparing train/test data"
    prep_res <- prepare_train_test_data(
      dataset = dataset_task,
      temp_random_seed = temp_random_seed,
      train_ratio = train_ratio,
      test_ratio = test_ratio,
      each_condition_total_sample_no = each.condition_total_sample_No,
      condition_col = "condition",
      preprocess_method = preprocessing_method,
      metadata_cols = metadata_cols_for_split,
      sampling_prob_col = sampling_prob_col,
      precomputed_split = split_plans[[temp_iter]]
    )
    
    worker_stage <- "extracting train/test matrices"
    train_data <- prep_res$train_data
    test_data <- prep_res$test_data
    feature <- intersect(feature_plans[[temp_iter]], prep_res$feature)
    if (length(feature) > 0) {
      train_metadata_cols <- setdiff(colnames(train_data), prep_res$feature)
      test_metadata_cols <- setdiff(colnames(test_data), prep_res$feature)
      train_data <- train_data[, c(train_metadata_cols, feature), drop = FALSE]
      test_data <- test_data[, c(test_metadata_cols, feature), drop = FALSE]
      train_raw_metadata_cols <- setdiff(colnames(prep_res$train_raw), prep_res$feature)
      test_raw_metadata_cols <- setdiff(colnames(prep_res$test_raw), prep_res$feature)
      prep_res$train_raw <- prep_res$train_raw[, c(train_raw_metadata_cols, feature), drop = FALSE]
      prep_res$test_raw <- prep_res$test_raw[, c(test_raw_metadata_cols, feature), drop = FALSE]
    }
    
    if (length(feature) == 0) {
      stop(paste0("No features available for split ", temp_split, " in iteration ", temp_iter, "."))
    }
    
    if (nrow(train_data) == 0 || nrow(test_data) == 0) {
      stop(paste0(
        "Empty train/test split for ", temp_split,
        " in iteration ", temp_iter,
        ". train rows = ", nrow(train_data),
        ", test rows = ", nrow(test_data), "."
      ))
    }
    
    x_train_mat <- data.matrix(train_data[, feature, drop = FALSE])
    x_test_mat <- data.matrix(test_data[, feature, drop = FALSE])
    
    if (nrow(x_train_mat) == 0 || nrow(x_test_mat) == 0 ||
        ncol(x_train_mat) == 0 || ncol(x_test_mat) == 0) {
      stop(paste0(
        "Empty matrix passed to modeling for ", temp_split,
        " in iteration ", temp_iter,
        ". x_train = ", nrow(x_train_mat), "x", ncol(x_train_mat),
        ", x_test = ", nrow(x_test_mat), "x", ncol(x_test_mat), "."
      ))
    }
    
    train_weights <- NULL
    if (patient_balance_method == "model_weights") {
      worker_stage <- "computing model weights"
      if ("balance_group" %in% colnames(train_data)) {
        train_weights <- compute_model_weights_by_group(train_data$balance_group)
      } else if ("sample_id" %in% colnames(train_data)) {
        train_weights <- compute_model_weights_by_group(train_data$sample_id)
      } else {
        stop("sample_id or balance_group metadata is required for model_weights balancing.")
      }
    }
    
    if (model_type == "xgboost") {
      worker_stage <- "training xgboost model"
      set.seed(model_seed)
      dtrain <- xgboost::xgb.DMatrix(data = x_train_mat, label = train_data$condition, weight = train_weights)
      n_cells <- nrow(train_data)
      n_features <- length(feature)
      xgb_params <- list(
        objective = "binary:logistic",
        eval_metric = "logloss",
        eta = ifelse(n_cells < 1000, 0.15, 0.1),
        max_depth = ifelse(n_cells < 500, 3, ifelse(n_cells < 5000, 5, 6)),
        min_child_weight = ifelse(n_cells < 500, 2, 1),
        subsample = ifelse(n_cells < 500, 0.9, 0.8),
        colsample_bytree = ifelse(n_features < 100, 1.0, ifelse(n_features < 2000, 0.5, 0.2)),
        tree_method = ifelse(n_cells > 2000 || n_features > 2000, "hist", "auto"),
        nthread = 1,
        seed = model_seed
      )
      xgb_nrounds <- ifelse(n_cells < 500, 20, ifelse(n_cells < 5000, 50, 100))
      fit <- xgboost::xgb.train(params = xgb_params, data = dtrain, nrounds = xgb_nrounds)
      worker_stage <- "predicting xgboost outputs"
      dtest <- xgboost::xgb.DMatrix(data = x_test_mat, label = test_data$condition)
      pred <- predict(fit, dtest)
      shap_contrib <- predict(fit, dtest, predcontrib = TRUE)
    } else if (model_type == "lightgbm") {
      worker_stage <- "training lightgbm model"
      set.seed(model_seed)
      dtrain <- lightgbm::lgb.Dataset(data = x_train_mat, label = train_data$condition, weight = train_weights)
      n_cells <- nrow(train_data)
      n_features <- length(feature)
      lgb_params <- list(
        objective = "binary",
        metric = "binary_logloss",
        learning_rate = ifelse(n_cells < 1000, 0.15, 0.1),
        max_depth = ifelse(n_cells < 500, 3, ifelse(n_cells < 5000, 5, 6)),
        num_leaves = ifelse(n_cells < 500, 7, ifelse(n_cells < 5000, 31, 63)),
        min_data_in_leaf = ifelse(n_cells < 500, 10, ifelse(n_cells < 5000, 20, 50)),
        feature_fraction = ifelse(n_features < 100, 1.0, ifelse(n_features < 2000, 0.5, 0.2)),
        bagging_fraction = ifelse(n_cells < 500, 0.9, 0.8),
        bagging_freq = 1,
        num_threads = 1,
        verbosity = -1,
        seed = model_seed,
        bagging_seed = model_seed,
        feature_fraction_seed = model_seed,
        data_random_seed = model_seed
      )
      lgb_nrounds <- ifelse(n_cells < 500, 20, ifelse(n_cells < 5000, 50, 100))
      fit <- lightgbm::lgb.train(params = lgb_params, data = dtrain, nrounds = lgb_nrounds)
      worker_stage <- "predicting lightgbm outputs"
      pred <- predict(fit, x_test_mat)
      shap_contrib <- predict(fit, x_test_mat, type = "contrib")
    } else {
      stop("model_type must be either 'xgboost' or 'lightgbm'")
    }
    
    if (is.null(dim(shap_contrib))) {
      shap_contrib <- matrix(shap_contrib, nrow = nrow(x_test_mat), byrow = TRUE)
    } else {
      shap_contrib <- as.matrix(shap_contrib)
    }
    
    if (nrow(shap_contrib) == 0 || ncol(shap_contrib) < 2) {
      stop(paste0("Invalid SHAP contribution matrix for ", temp_split, " in iteration ", temp_iter, "."))
    }
    
    if (is.null(colnames(shap_contrib)) || length(colnames(shap_contrib)) != ncol(shap_contrib)) {
      colnames(shap_contrib) <- c(feature, "BIAS")
    }
    
    baseline_value <- mean(shap_contrib[, ncol(shap_contrib)], na.rm = TRUE)
    shap_values <- shap_contrib[, -ncol(shap_contrib), drop = FALSE]
    colnames(shap_values) <- feature
    shap_values[!is.finite(shap_values)] <- 0
    
    prediction <- ifelse(pred > 0.5, task_what_1, task_what_0)
    
    worker_stage <- "summarizing model performance"
    tbl_1 <- table(
      Predicted = prediction,
      Actual = ifelse(test_data$condition == 0, task_what_0, task_what_1)
    )
    
    confusionMatrix_table <- cbind(iteration = temp_iter, as.data.frame(tbl_1))
    ROC_data_table <- data.frame(iteration = temp_iter, condition = test_data$condition, pred = pred)
    
    roc_iter <- pROC::roc(
      response = test_data$condition,
      predictor = pred,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )
    
    auc_iter <- as.numeric(pROC::auc(roc_iter))
    
    con_ACC <- caret::confusionMatrix(
      data = factor(prediction),
      reference = factor(ifelse(test_data$condition == 0, task_what_0, task_what_1)),
      positive = task_what_1
    )
    
    model_performance <- t(data.frame(round(
      c(
        con_ACC$overall["Accuracy"],
        con_ACC$byClass[c("Sensitivity", "Specificity", "Precision", "Recall", "F1")]
      ),
      digits = 2
    )))
    row.names(model_performance) <- temp_iter
    
    iteration_importance <- colMeans(abs(shap_values), na.rm = TRUE)
    total_iteration_importance <- sum(iteration_importance, na.rm = TRUE)
    importance_data <- data.frame(
      iteration = temp_iter,
      feature = colnames(shap_values),
      `average|SHAP value|` = iteration_importance,
      check.names = FALSE
    )
    
    mean_shap_0 <- colMeans(shap_values[test_data$condition == 0, , drop = FALSE])
    mean_shap_1 <- colMeans(shap_values[test_data$condition == 1, , drop = FALSE])
    delta_shap <- mean_shap_1 - mean_shap_0
    patient_delta <- compute_patient_balanced_shap_effect(
      shap_matrix = shap_values,
      metadata = test_data[, c("sample_id", "condition"), drop = FALSE],
      sample_col = "sample_id",
      condition_col = "condition"
    )
    patient_delta_values <- stats::setNames(
      patient_delta$effect_size,
      patient_delta$feature
    )
    value_shap_cor_iteration <- vapply(feature, function(feature_name) {
      feature_values <- as.numeric(test_data[[feature_name]])
      feature_shap <- as.numeric(shap_values[, feature_name])
      valid <- is.finite(feature_values) & is.finite(feature_shap)
      if (sum(valid) < 3L ||
          length(unique(feature_values[valid])) < 2L ||
          length(unique(feature_shap[valid])) < 2L) {
        return(NA_real_)
      }
      suppressWarnings(stats::cor(
        feature_values[valid],
        feature_shap[valid],
        method = "spearman"
      ))
    }, numeric(1))
    model_reproducibility_iteration <- data.frame(
      iteration = temp_iter,
      feature = feature,
      shap_importance = iteration_importance[feature],
      delta_shap = patient_delta_values[feature],
      relative_delta_shap = abs(patient_delta_values[feature]) /
        max(total_iteration_importance, .Machine$double.eps),
      value_shap_cor = value_shap_cor_iteration[feature],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    if (memory_mode == "safe") {
      worker_stage <- "writing safe-mode iteration chunk"
      chunk_file <- file.path(temp_chunk_dir, sprintf("iter_%03d_chunk.rds", temp_iter))
      chunk_metadata <- test_data[, intersect(c(".shapture_unit_id", "sample_id", "condition"), colnames(test_data)), drop = FALSE]
      saveRDS(
        list(
          shap_values = shap_values,
          feature_values = as.data.frame(test_data[, feature, drop = FALSE], check.names = FALSE),
          metadata = chunk_metadata,
          baseline = baseline_value
        ),
        file = chunk_file
      )
      
      rm(train_data, test_data, x_train_mat, x_test_mat, shap_contrib, shap_values, dtrain)
      rm(mean_shap_0, mean_shap_1)
      gc(verbose = FALSE)
      
      return(list(
        iteration = temp_iter,
        fit = NULL,
        pred = pred,
        confusionMatrix_table = confusionMatrix_table,
        ROC_data_table = ROC_data_table,
        auc_table = data.frame(iteration = temp_iter, auc = auc_iter),
        con_ACC = con_ACC,
        model_performance = model_performance,
        importance_data = importance_data,
        delta_shap = delta_shap,
        model_reproducibility_iteration = model_reproducibility_iteration,
        chunk_file = chunk_file
      ))
    }
    
    worker_stage <- "returning full-mode iteration result"
    list(
      iteration = temp_iter,
      fit = fit,
      pred = pred,
      confusionMatrix_table = confusionMatrix_table,
      ROC_data_table = ROC_data_table,
      auc_table = data.frame(iteration = temp_iter, auc = auc_iter),
      con_ACC = con_ACC,
      model_performance = model_performance,
      importance_data = importance_data,
      shp_S = shap_values,
      shp_X = as.data.frame(test_data[, feature, drop = FALSE], check.names = FALSE),
      shp_X_unscaled = as.data.frame(prep_res$test_raw[, feature, drop = FALSE], check.names = FALSE),
      test_data = test_data,
      baseline = baseline_value,
      delta_shap = delta_shap,
      model_reproducibility_iteration = model_reproducibility_iteration
    )
    }, error = function(e) {
      call_txt <- paste(deparse(conditionCall(e)), collapse = " ")
      if (!nzchar(call_txt)) {
        call_txt <- "unknown call"
      }
      stop(
        paste0(
          "split ", temp_split,
          ", iteration ", temp_iter,
          " failed during ", worker_stage,
          ". Worker cwd: ", getwd(),
          ". temp_chunk_dir: ", temp_chunk_dir,
          ". temp_chunk_dir_exists: ", dir.exists(temp_chunk_dir),
          ". temp_chunk_dir_writable: ", ifelse(dir.exists(temp_chunk_dir), file.access(temp_chunk_dir, 2) == 0, FALSE),
          ". Error call: ", call_txt,
          ". Error message: ", conditionMessage(e)
        ),
        call. = FALSE
      )
    })
  }

  model_workers_stopped <- safe_stop_shapture_cluster(model_cluster)
  if (!isTRUE(model_workers_stopped)) {
    stop("Model workers could not be stopped safely before Cross-patient bootstrapping.")
  }
  model_cluster <- NULL
  foreach::registerDoSEQ()
  gc(verbose = FALSE)
  
  fit_list <- lapply(iter_results, `[[`, "fit")
  pred_list <- lapply(iter_results, `[[`, "pred")
  delta_list <- lapply(iter_results, `[[`, "delta_shap")
  confusionMatrix_all <- do.call(rbind, lapply(iter_results, `[[`, "confusionMatrix_table"))
  ROC_all <- do.call(rbind, lapply(iter_results, `[[`, "ROC_data_table"))
  AUC_all <- do.call(rbind, lapply(iter_results, `[[`, "auc_table"))
  model_performance_all <- do.call(rbind, lapply(iter_results, `[[`, "model_performance"))
  importance_data_all <- do.call(rbind, lapply(iter_results, `[[`, "importance_data"))
  model_reproducibility_iteration_all <- do.call(
    rbind,
    lapply(iter_results, `[[`, "model_reproducibility_iteration")
  )
  if (memory_mode == "full") {
    shp_S <- rbind_feature_matrices(lapply(iter_results, `[[`, "shp_S"), nonzero_feature_cols)
    shp_X <- as.data.frame(
      rbind_feature_matrices(lapply(iter_results, `[[`, "shp_X"), nonzero_feature_cols),
      check.names = FALSE
    )
    shp_X_unscaled <- as.data.frame(
      rbind_feature_matrices(lapply(iter_results, `[[`, "shp_X_unscaled"), nonzero_feature_cols),
      check.names = FALSE
    )
    shp_testdata <- do.call(rbind, lapply(iter_results, function(res) {
      res$test_data[, intersect(c(".shapture_unit_id", "sample_id", "condition"), colnames(res$test_data)), drop = FALSE]
    }))
    shp_baseline <- unlist(lapply(iter_results, `[[`, "baseline"))
    chunk_files <- character(0)
    if (".shapture_unit_id" %in% colnames(shp_testdata)) {
      unit_avg_full <- average_shap_by_unit_in_memory(
        shap_values = shp_S,
        feature_values = shp_X,
        metadata = shp_testdata,
        feature_names = colnames(shp_S),
        unit_col = ".shapture_unit_id"
      )
      shp_S <- unit_avg_full$S
      shp_X <- unit_avg_full$X
      shp_testdata <- unit_avg_full$metadata
      rm(unit_avg_full)
      gc(verbose = FALSE)
    }
  } else {
    shp_S <- NULL
    shp_X <- NULL
    shp_X_unscaled <- NULL
    shp_testdata <- NULL
    shp_baseline <- NULL
    chunk_files <- vapply(iter_results, `[[`, character(1), "chunk_file")
  }
  
  for (res in iter_results) {
    sink("docs/model performance.txt", append = TRUE)
    cat("iteration ", res$iteration, "\n")
    print(res$con_ACC)
    cat(".....................", "\n")
    sink()
  }
  
  if (patient_split_strategy_used == "unseen_patients" &&
      !is.null(patient_coverage_summary)) {
    held_out_patient_counts <- patient_coverage_summary[
      patient_coverage_summary$eligible &
        patient_coverage_summary$test_iterations > 0,
      ,
      drop = FALSE
    ]
    patient_count_condition0 <- length(unique(
      held_out_patient_counts$patient_id[
        as.character(held_out_patient_counts$condition) == "0"
      ]
    ))
    patient_count_condition1 <- length(unique(
      held_out_patient_counts$patient_id[
        as.character(held_out_patient_counts$condition) == "1"
      ]
    ))
    min_patients <- min(patient_count_condition0, patient_count_condition1)
  } else if ("sample_id" %in% colnames(dataset_task)) {
    n_patients_per_condition <- dataset_task %>%
      dplyr::group_by(condition) %>%
      dplyr::summarise(n_patients = dplyr::n_distinct(sample_id), .groups = "drop")
    
    patient_count_condition0 <- n_patients_per_condition$n_patients[n_patients_per_condition$condition == 0]
    patient_count_condition1 <- n_patients_per_condition$n_patients[n_patients_per_condition$condition == 1]
    patient_count_condition0 <- ifelse(length(patient_count_condition0) == 0, NA_integer_, patient_count_condition0)
    patient_count_condition1 <- ifelse(length(patient_count_condition1) == 0, NA_integer_, patient_count_condition1)
    min_patients <- min(c(patient_count_condition0, patient_count_condition1), na.rm = TRUE)
  } else {
    patient_count_condition0 <- NA_integer_
    patient_count_condition1 <- NA_integer_
    min_patients <- NA_integer_
  }
  
  delta_mat <- rbind_named_vectors(delta_list, nonzero_feature_cols)
  stability <- apply(delta_mat, 2, function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) {
      return(NA_real_)
    }
    mean(sign(x) == sign(mean(x)))
  })
  
  if (memory_mode == "full") {
    unit_shap_matrix <- as.matrix(shp_S)
    unit_feature_matrix <- as.matrix(shp_X)
    unit_shap_metadata <- as.data.frame(shp_testdata, check.names = FALSE)
  } else {
    unit_average_all <- average_shap_by_unit_from_chunks(
      chunk_files = chunk_files,
      feature_names = nonzero_feature_cols
    )
    unit_shap_matrix <- as.matrix(unit_average_all$S)
    unit_feature_matrix <- as.matrix(unit_average_all$X)
    unit_shap_metadata <- as.data.frame(unit_average_all$metadata, check.names = FALSE)
    rm(unit_average_all)
  }

  patient_shap_effect <- compute_patient_balanced_shap_effect(
    shap_matrix = unit_shap_matrix,
    metadata = unit_shap_metadata,
    sample_col = "sample_id",
    condition_col = "condition"
  )
  
  stats_results <- data.frame(
    feature = patient_shap_effect$feature,
    effect_size = patient_shap_effect$effect_size,
    stability = stability[patient_shap_effect$feature],
    n_patients_condition0 = patient_count_condition0,
    n_patients_condition1 = patient_count_condition1,
    min_patients = min_patients
  )
  
  expression_matrix <- as.matrix(dataset_task[, nonzero_feature_cols, drop = FALSE])
  mean_expr_0 <- colMeans(expression_matrix[dataset_task$condition == 0, , drop = FALSE])
  mean_expr_1 <- colMeans(expression_matrix[dataset_task$condition == 1, , drop = FALSE])
  log_pseudocount <- compute_dynamic_pseudocount(x = expression_matrix, fraction = 0.01, fallback = 1e-6)
  expression_log2FC <- log2((mean_expr_1 + log_pseudocount) / (mean_expr_0 + log_pseudocount))
  missing_logfc_features <- setdiff(stats_results$feature, names(expression_log2FC))
  if (length(missing_logfc_features) > 0) {
    stop(
      "Missing expression_log2FC values for modeled features in ", temp_split, ": ",
      paste(head(missing_logfc_features, 10), collapse = ", ")
    )
  }
  stats_results$expression_log2FC <- expression_log2FC[stats_results$feature]
  stats_results$expression_log2FC_abs <- abs(stats_results$expression_log2FC)
  rm(expression_matrix)
  gc(verbose = FALSE)
  
  single_cell_features <- intersect(nonzero_feature_cols, colnames(dataset_task_singlecell))
  dataset_task_singlecell$condition <- ifelse(dataset_task_singlecell$condition == task_what_0, 0, 1)
  binary_mat <- as.matrix(dataset_task_singlecell[, single_cell_features, drop = FALSE]) > 0
  det_rate_0 <- colMeans(binary_mat[dataset_task_singlecell$condition == 0, , drop = FALSE])
  det_rate_1 <- colMeans(binary_mat[dataset_task_singlecell$condition == 1, , drop = FALSE])
  det_diff <- det_rate_1 - det_rate_0
  missing_det_diff_features <- setdiff(stats_results$feature, names(det_diff))
  if (length(missing_det_diff_features) > 0) {
    warning(
      "Detection-rate differences are unavailable for ",
      length(missing_det_diff_features),
      " modeled feature(s) in ",
      temp_split,
      " after single-cell feature filtering. These values will be set to NA."
    )
  }
  stats_results$det_diff <- det_diff[stats_results$feature]
  stats_results$det_diff_abs <- abs(stats_results$det_diff)
  rm(binary_mat)
  gc(verbose = FALSE)

  # Keep expression consistency on original cells, independent of modeling unit.
  within_condition_consistency <- compute_within_condition_expression_consistency(
    expression_data = dataset_task_singlecell[, single_cell_features, drop = FALSE],
    condition = dataset_task_singlecell$condition,
    condition_levels = c(0, 1),
    upper_quantile = 0.99,
    min_positive_cells = 20L,
    pseudocount_fraction = 0.01,
    pseudocount_fallback = 1e-6
  )
  both_zero_features <- within_condition_consistency$feature[
    within_condition_consistency$both_conditions_all_zero
  ]
  if (length(both_zero_features) > 0) {
    warning(
      "Removing ",
      length(both_zero_features),
      " feature(s) with homogeneous absence in both task conditions: ",
      paste(head(both_zero_features, 10), collapse = ", ")
    )
    stats_results <- stats_results[
      !stats_results$feature %in% both_zero_features,
      ,
      drop = FALSE
    ]
  }
  if (nrow(stats_results) == 0) {
    stop("No modeled features remain after removing features absent in both conditions.")
  }

  consistency_match <- match(
    stats_results$feature,
    within_condition_consistency$feature
  )
  if (anyNA(consistency_match)) {
    missing_consistency_features <- stats_results$feature[is.na(consistency_match)]
    stop(
      "Missing within-condition expression consistency values for modeled features in ",
      temp_split,
      ": ",
      paste(head(missing_consistency_features, 10), collapse = ", ")
    )
  }
  consistency_columns <- setdiff(
    colnames(within_condition_consistency),
    c("feature", "both_conditions_all_zero")
  )
  stats_results[, consistency_columns] <- within_condition_consistency[
    consistency_match,
    consistency_columns,
    drop = FALSE
  ]
  rm(within_condition_consistency)
  gc(verbose = FALSE)
  
  common_features <- intersect(colnames(unit_shap_matrix), colnames(unit_feature_matrix))
  shap_expr_cor <- vapply(common_features, function(f) {
    suppressWarnings(
      stats::cor(
        unit_feature_matrix[, f],
        unit_shap_matrix[, f],
        method = "spearman",
        use = "pairwise.complete.obs"
      )
    )
  }, numeric(1))
  shap_expr_cor[!is.finite(shap_expr_cor)] <- 0
  shap_expr_cor_abs <- abs(shap_expr_cor)
  missing_shap_expr_cor_features <- setdiff(stats_results$feature, names(shap_expr_cor))
  if (length(missing_shap_expr_cor_features) > 0) {
    stop(
      "Missing SHAP-expression correlations for modeled features in ", temp_split, ": ",
      paste(head(missing_shap_expr_cor_features, 10), collapse = ", ")
    )
  }
  stats_results$shap_expr_cor <- shap_expr_cor[stats_results$feature]
  stats_results$shap_expr_cor_abs <- shap_expr_cor_abs[stats_results$feature]

  model_reproducibility_result <- compute_model_reproducibility(
    iteration_data = model_reproducibility_iteration_all,
    full_stats = stats_results,
    shap_expr_cor_cutoff = shap_expr_cor_cutoff,
    lower_probability = 0.10
  )
  reproducibility_match <- match(
    stats_results$feature,
    model_reproducibility_result$feature
  )
  if (anyNA(reproducibility_match)) {
    stop("Model reproducibility did not cover every modeled feature in ", temp_split, ".")
  }
  stats_results$legacy_delta_direction_stability <- stats_results$stability
  stats_results$model_reproducibility <- model_reproducibility_result$model_reproducibility[
    reproducibility_match
  ]
  stats_results$stability <- stats_results$model_reproducibility
  reproducibility_columns <- setdiff(
    colnames(model_reproducibility_result),
    c("feature", "model_reproducibility")
  )
  stats_results[, reproducibility_columns] <- model_reproducibility_result[
    reproducibility_match,
    reproducibility_columns,
    drop = FALSE
  ]
  write.csv(
    model_reproducibility_result,
    "docs/model_reproducibility_by_feature.csv",
    row.names = FALSE
  )
  write.csv(
    model_reproducibility_iteration_all,
    "docs/model_reproducibility_by_iteration.csv",
    row.names = FALSE
  )

  average_importance_for_ranking <- data.frame(
    feature = colnames(unit_shap_matrix),
    `average|SHAP value|` = colMeans(abs(unit_shap_matrix), na.rm = TRUE),
    check.names = FALSE
  )
  score_anchor_features <- intersect(
    stats_results$feature,
    average_importance_for_ranking$feature
  )
  if (length(score_anchor_features) == 0L) {
    stop("No fixed candidate features are available for score anchors in ", temp_split, ".")
  }
  score_anchor_importance_match <- match(
    score_anchor_features,
    average_importance_for_ranking$feature
  )
  score_anchor_stats_match <- match(score_anchor_features, stats_results$feature)
  score_anchor_importance_values <- average_importance_for_ranking$`average|SHAP value|`[
    score_anchor_importance_match
  ]
  score_anchor_log2fc_values <- abs(
    stats_results$expression_log2FC[score_anchor_stats_match]
  )
  shap_importance_median_anchor <- compute_positive_median_anchor(
    score_anchor_importance_values,
    component_name = "SHAP importance"
  )
  expression_log2FC_abs_median_anchor <- compute_positive_median_anchor(
    score_anchor_log2fc_values,
    component_name = "absolute expression log2FC"
  )
  score_anchor_summary <- data.frame(
    component = c("SHAP importance", "Absolute expression log2FC"),
    transformation = "x / (x + positive median anchor)",
    median_anchor = c(
      shap_importance_median_anchor,
      expression_log2FC_abs_median_anchor
    ),
    candidate_features_used = length(score_anchor_features),
    positive_values_used = c(
      sum(is.finite(score_anchor_importance_values) & score_anchor_importance_values > 1e-8),
      sum(is.finite(score_anchor_log2fc_values) & score_anchor_log2fc_values > 1e-8)
    ),
    anchor_fitted_before_direction_and_pattern_filters = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  n_features_before_direction_filter <- nrow(stats_results)
  if (shap_expr_cor_direction_filter == "positive") {
    stats_results <- stats_results[stats_results$shap_expr_cor > 0, , drop = FALSE]
  } else if (shap_expr_cor_direction_filter == "negative") {
    stats_results <- stats_results[stats_results$shap_expr_cor < 0, , drop = FALSE]
  }
  
  n_features_after_direction_filter <- nrow(stats_results)
  if (n_features_after_direction_filter == 0) {
    stop(
      "No features remain after applying shap_expr_cor_direction_filter = '",
      shap_expr_cor_direction_filter,
      "' in ",
      temp_split,
      "."
    )
  }
  
  support_pattern <- compute_model_expression_pattern(
    expression_log2FC = stats_results$expression_log2FC,
    shap_expr_cor = stats_results$shap_expr_cor,
    what_0 = task_what_0,
    what_1 = task_what_1,
    expression_log2FC_cutoff = feature_min_log2FC,
    shap_expr_cor_cutoff = shap_expr_cor_cutoff
  )
  insufficient_pattern <- unname(
    model_expression_pattern_labels(task_what_0, task_what_1)["insufficient"]
  )
  support_feature_names <- stats_results$feature[support_pattern != insufficient_pattern]
  if (length(support_feature_names) == 0) {
    stop("No features remain for Cross-patient support after model-expression pattern filtering.")
  }

  cross_patient_result <- compute_cross_patient_support(
    shap_matrix = unit_shap_matrix,
    shap_feature_values = unit_feature_matrix,
    shap_metadata = unit_shap_metadata,
    model_expression_data = dataset_task[, c("sample_id", "condition", support_feature_names), drop = FALSE],
    single_cell_expression_data = dataset_task_singlecell[, c("sample_id", "condition", support_feature_names), drop = FALSE],
    full_stats = stats_results,
    average_importance = average_importance_for_ranking,
    feature_names = support_feature_names,
    comparison_mode = comparison_mode,
    feature_min_log2FC = feature_min_log2FC,
    shap_expr_cor_cutoff = shap_expr_cor_cutoff,
    expression_log2FC_abs_median_anchor = expression_log2FC_abs_median_anchor,
    n_bootstraps = cross_patient_bootstraps,
    seed = seed_plan$cross_patient_bootstrap_seed,
    n_cores = n_cores,
    worker_log_file = file.path(getwd(), "docs", "cross_patient_bootstrap_workers.log")
  )
  cross_patient_result$runtime$model_workers_stopped_before_bootstrap <-
    model_workers_stopped
  support_match <- match(stats_results$feature, cross_patient_result$support$feature)
  stats_results$cross_patient_support <- cross_patient_result$support$cross_patient_support[support_match]
  stats_results$cross_patient_support_available <- cross_patient_result$available
  stats_results$cross_patient_support_unavailable_reason <- cross_patient_result$reason
  
  mean_auc <- mean(AUC_all$auc, na.rm = TRUE)
  sd_auc <- stats::sd(AUC_all$auc, na.rm = TRUE)
  auc_ci <- stats::quantile(AUC_all$auc, probs = c(0.025, 0.975), na.rm = TRUE)
  
  roc_plot_all_iteration <- pROC::roc(
    response = ROC_all$condition,
    predictor = ROC_all$pred,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  combined_results <- compute_combined_shap_ranking(
    average_importance = average_importance_for_ranking,
    stats_results = stats_results,
    adaptive_weights = if_adaptive_weights,
    min_adaptive_weight = 0.03,
    default_prior_scores_weights = default_prior_scores_weights,
    what_0 = task_what_0,
    what_1 = task_what_1,
    expression_log2FC_cutoff = feature_min_log2FC,
    shap_expr_cor_cutoff = shap_expr_cor_cutoff,
    shap_importance_median_anchor = shap_importance_median_anchor,
    expression_log2FC_abs_median_anchor = expression_log2FC_abs_median_anchor,
    score_scaling_method = score_scaling_method,
    score_scaling_probs = score_scaling_probs
  )
  
  combined_results$rank <- seq_len(nrow(combined_results))
  top_feature_label_actual <- min(top_feature_label, nrow(combined_results))
  top_feature_show_actual <- min(top_feature_show_local, nrow(combined_results))
  retained_feature_names <- intersect(as.character(combined_results$feature), nonzero_feature_cols)
  if (length(retained_feature_names) == 0) {
    stop("No retained ranked features are available for SHAP-only summaries in ", temp_split, ".")
  }
  min_patients_value <- unique(combined_results$min_patients)[1]
  cross_patient_support_available <- isTRUE(cross_patient_result$available)
  cross_patient_support_note <- if (cross_patient_support_available) {
    paste0(
      "Cross-patient support used ", cross_patient_result$n_bootstraps_valid_task,
      " valid patient-block bootstraps out of ", cross_patient_result$n_bootstraps_requested, "."
    )
  } else {
    "Cross-patient support unavailable: insufficient patients. Final scores use four renormalized components."
  }
  
  stats_interpretation <- data.frame(
    split_level = temp_split,
    comparison_mode = comparison_mode,
    split_by_col = ifelse(is.null(effective_split_by_col), NA_character_, effective_split_by_col),
    one_vs_others_col = ifelse(is.null(one_vs_others_col), NA_character_, one_vs_others_col),
    n_patients_condition0 = patient_count_condition0,
    n_patients_condition1 = patient_count_condition1,
    min_patients = min_patients,
    cross_patient_eligible_patients_condition0 = cross_patient_result$n_patients_condition0,
    cross_patient_eligible_patients_condition1 = cross_patient_result$n_patients_condition1,
    cross_patient_eligible_patients_total = cross_patient_result$n_patients_total,
    patient_split_strategy_requested = patient_split_strategy_requested,
    patient_split_strategy_used = patient_split_strategy_used,
    patient_split_fallback_reason = patient_split_fallback_reason,
    cross_patient_support_available = cross_patient_support_available,
    cross_patient_support_unavailable_reason = cross_patient_result$reason,
    cross_patient_support_method = "patient_block_bootstrap",
    cross_patient_support_bootstraps_requested = cross_patient_result$n_bootstraps_requested,
    cross_patient_support_bootstraps_valid = cross_patient_result$n_bootstraps_valid_task,
    cross_patient_support_note = cross_patient_support_note,
    shap_expr_cor_direction_filter = shap_expr_cor_direction_filter,
    n_features_before_direction_filter = n_features_before_direction_filter,
    n_features_after_direction_filter = n_features_after_direction_filter,
    n_features_before_model_expression_pattern_filter = unique(combined_results$n_features_before_model_expression_pattern_filter)[1],
    n_features_after_model_expression_pattern_filter = unique(combined_results$n_features_after_model_expression_pattern_filter)[1],
    interpretation = cross_patient_support_note,
    stringsAsFactors = FALSE
  )
  
  write.csv(stats_interpretation, "docs/stats_interpretation.csv", row.names = FALSE)
  write.csv(score_anchor_summary, "docs/score_anchor_summary.csv", row.names = FALSE)
  write.csv(
    data.frame(
      cross_patient_support_method = "patient_block_bootstrap",
      cross_patient_support_available = cross_patient_support_available,
      cross_patient_support_unavailable_reason = cross_patient_result$reason,
      n_bootstraps_requested = cross_patient_result$n_bootstraps_requested,
      n_bootstraps_valid = cross_patient_result$n_bootstraps_valid_task,
      cross_patient_bootstrap_seed = seed_plan$cross_patient_bootstrap_seed,
      patient_split_strategy_requested = patient_split_strategy_requested,
      patient_split_strategy_used = patient_split_strategy_used,
      shap_expr_cor_direction_filter = shap_expr_cor_direction_filter,
      scoring_strategy = "median_anchor_unbounded_natural_bounded",
      shap_importance_median_anchor = shap_importance_median_anchor,
      expression_log2FC_abs_median_anchor = expression_log2FC_abs_median_anchor,
      cross_patient_support_note = cross_patient_support_note,
      stringsAsFactors = FALSE
    ),
    "docs/cross_patient_support_method_used.csv",
    row.names = FALSE
  )
  write.csv(
    cross_patient_result$diagnostics,
    "docs/cross_patient_support_by_feature.csv",
    row.names = FALSE
  )
  write.csv(
    cross_patient_result$runtime,
    "docs/cross_patient_bootstrap_runtime.csv",
    row.names = FALSE
  )
  
  weight_table <- data.frame(
    component = c("SHAP_importance", "cross_patient_support", "model_reproducibility", "expression_signal", "value_SHAP_relationship"),
    weight = c(
      combined_results$weight_importance[1],
      combined_results$weight_cross_patient_support[1],
      combined_results$weight_model_reproducibility[1],
      combined_results$weight_expression[1],
      combined_results$weight_value_shap[1]
    )
  )
  write.csv(weight_table, paste0("docs/adaptive_weights_used_", if_adaptive_weights, ".csv"), row.names = FALSE)
  
  write.csv(AUC_all, "docs/AUC_by_iteration.csv", row.names = FALSE)
  write.csv(
    data.frame(
      pooled_auc = as.numeric(roc_plot_all_iteration$auc),
      mean_iteration_auc = mean_auc,
      sd_iteration_auc = sd_auc,
      ci_lower_iteration_auc = auc_ci[1],
      ci_upper_iteration_auc = auc_ci[2]
    ),
    "docs/AUC_summary.csv",
    row.names = FALSE
  )
  
  average_model_performance <- round(c("AUC" = mean_auc, colMeans(model_performance_all)), 2)
  write.csv(data.frame(average_model_performance), "docs/average_model_performance.csv")

  rows_by_condition <- table(dataset_task$condition)
  get_condition_row_count <- function(level) {
    val <- rows_by_condition[as.character(level)]
    if (length(val) == 0 || is.na(val)) {
      return(0L)
    }
    as.integer(val)
  }
  get_coverage_value <- function(level, column) {
    val <- sampling_coverage_summary[
      as.character(sampling_coverage_summary$condition) == as.character(level),
      column,
      drop = TRUE
    ]
    if (length(val) == 0 || is.na(val)) {
      return(NA_real_)
    }
    val[[1]]
  }
  run_summary <- data.frame(
    split_level = temp_split,
    comparison_mode = comparison_mode,
    split_by_col = ifelse(is.null(effective_split_by_col), NA_character_, effective_split_by_col),
    condition_0 = task_what_0,
    condition_1 = task_what_1,
    analysis_seed = analysis_seed,
    task_seed = seed_plan$task_seed,
    metacell_seed = metacell_seed,
    cross_patient_bootstrap_seed = seed_plan$cross_patient_bootstrap_seed,
    analysis_unit_type = analysis_unit_type,
    model_type = model_type,
    preprocessing_method = preprocessing_method,
    patient_balance_method = patient_balance_method,
    patient_split_strategy_requested = patient_split_strategy_requested,
    patient_split_strategy_used = patient_split_strategy_used,
    patient_split_fallback_reason = patient_split_fallback_reason,
    unseen_patient_train_fraction = if (patient_split_strategy_used == "unseen_patients") 0.60 else NA_real_,
    unseen_patient_test_fraction = if (patient_split_strategy_used == "unseen_patients") 0.40 else NA_real_,
    unseen_min_patients_per_condition = if (patient_split_strategy_requested == "unseen_patients") 5L else NA_integer_,
    unseen_min_units_per_patient = if (patient_split_strategy_requested == "unseen_patients") min_unseen_units_per_patient else NA_integer_,
    unseen_min_train_units_per_condition = if (patient_split_strategy_requested == "unseen_patients") min_unseen_train_units_per_condition else NA_integer_,
    train_patients_condition0 = if (patient_split_strategy_used == "unseen_patients") split_plan_res$n_train_patients[[1]] else NA_integer_,
    train_patients_condition1 = if (patient_split_strategy_used == "unseen_patients") split_plan_res$n_train_patients[[2]] else NA_integer_,
    test_patients_condition0 = if (patient_split_strategy_used == "unseen_patients") split_plan_res$n_test_patients[[1]] else NA_integer_,
    test_patients_condition1 = if (patient_split_strategy_used == "unseen_patients") split_plan_res$n_test_patients[[2]] else NA_integer_,
    min_patient_test_iterations = if (patient_split_strategy_used == "unseen_patients") {
      min(patient_coverage_summary$test_iterations[patient_coverage_summary$eligible])
    } else {
      NA_real_
    },
    median_patient_test_iterations = if (patient_split_strategy_used == "unseen_patients") {
      stats::median(patient_coverage_summary$test_iterations[patient_coverage_summary$eligible])
    } else {
      NA_real_
    },
    max_patient_test_iterations = if (patient_split_strategy_used == "unseen_patients") {
      max(patient_coverage_summary$test_iterations[patient_coverage_summary$eligible])
    } else {
      NA_real_
    },
    memory_mode = memory_mode,
    iteration_no = iteration_no,
    feature_min_pct = feature_min_pct,
    feature_min_log2FC = feature_min_log2FC,
    shap_expr_cor_cutoff = shap_expr_cor_cutoff,
    feature_bagging_fraction = feature_bagging_fraction,
    score_scaling_method = "median_anchor_unbounded_natural_bounded",
    legacy_score_scaling_method_requested = score_scaling_method,
    adaptive_weight_winsor_probs = paste(score_scaling_probs, collapse = ", "),
    shap_importance_median_anchor = shap_importance_median_anchor,
    expression_log2FC_abs_median_anchor = expression_log2FC_abs_median_anchor,
    if_adaptive_weights = if_adaptive_weights,
    shap_expr_cor_direction_filter = shap_expr_cor_direction_filter,
    visualization_inputs_saved = TRUE,
    n_rows_task_singlecell = n_rows_task_singlecell,
    n_rows_task_model_initial = n_rows_task_model_initial,
    n_rows_task_model_after_balance = nrow(dataset_task),
    n_rows_condition0_after_balance = get_condition_row_count(0),
    n_rows_condition1_after_balance = get_condition_row_count(1),
    test_coverage_percent_condition0 = get_coverage_value(0, "percent_test_units"),
    test_coverage_percent_condition1 = get_coverage_value(1, "percent_test_units"),
    test_coverage_percent_min = min(
      get_coverage_value(0, "percent_test_units"),
      get_coverage_value(1, "percent_test_units"),
      na.rm = TRUE
    ),
    min_test_appearances_condition0 = get_coverage_value(0, "min_test_appearances"),
    min_test_appearances_condition1 = get_coverage_value(1, "min_test_appearances"),
    n_candidate_features = n_candidate_features,
    n_features_after_feature_prefilter = n_features_after_prefilter,
    n_features_after_nonzero_filter = length(nonzero_feature_cols),
    candidate_genes_available_after_prefiltering = n_features_after_prefilter,
    candidate_genes_available_for_feature_bagging = feature_bagging_summary$total_candidate_features[1],
    genes_used_per_iteration = feature_bagging_summary$features_per_iteration[1],
    min_gene_inclusion_count = feature_bagging_summary$min_gene_inclusion_count[1],
    median_gene_inclusion_count = feature_bagging_summary$median_gene_inclusion_count[1],
    max_gene_inclusion_count = feature_bagging_summary$max_gene_inclusion_count[1],
    unique_features_included_by_feature_bagging = feature_bagging_summary$unique_features_included[1],
    percent_features_included_by_feature_bagging = feature_bagging_summary$percent_features_included[1],
    min_patient_condition_coverage_per_gene = if (is.null(patient_feature_coverage_summary)) {
      NA_real_
    } else {
      min(patient_feature_coverage_summary$patient_condition_units_with_shap)
    },
    median_patient_condition_coverage_per_gene = if (is.null(patient_feature_coverage_summary)) {
      NA_real_
    } else {
      stats::median(patient_feature_coverage_summary$patient_condition_units_with_shap)
    },
    max_patient_condition_coverage_per_gene = if (is.null(patient_feature_coverage_summary)) {
      NA_real_
    } else {
      max(patient_feature_coverage_summary$patient_condition_units_with_shap)
    },
    n_features_before_shap_expr_cor_direction_filter = n_features_before_direction_filter,
    n_features_after_shap_expr_cor_direction_filter = n_features_after_direction_filter,
    n_features_before_model_expression_pattern_filter = unique(combined_results$n_features_before_model_expression_pattern_filter)[1],
    n_features_after_model_expression_pattern_filter = unique(combined_results$n_features_after_model_expression_pattern_filter)[1],
    n_final_ranked_features = nrow(combined_results),
    n_top_features_requested = top_feature_show,
    n_top_features_plotted = top_feature_show_actual,
    n_patients_condition0 = patient_count_condition0,
    n_patients_condition1 = patient_count_condition1,
    min_patients = min_patients,
    cross_patient_eligible_patients_condition0 = cross_patient_result$n_patients_condition0,
    cross_patient_eligible_patients_condition1 = cross_patient_result$n_patients_condition1,
    cross_patient_eligible_patients_total = cross_patient_result$n_patients_total,
    cross_patient_support_available = cross_patient_support_available,
    cross_patient_support_unavailable_reason = cross_patient_result$reason,
    cross_patient_support_method = "patient_block_bootstrap",
    cross_patient_support_bootstraps_requested = cross_patient_result$n_bootstraps_requested,
    cross_patient_support_bootstraps_valid = cross_patient_result$n_bootstraps_valid_task,
    cross_patient_bootstrap_backend = cross_patient_result$runtime$backend,
    cross_patient_bootstrap_workers_requested = cross_patient_result$runtime$requested_workers,
    cross_patient_bootstrap_workers_selected = cross_patient_result$runtime$selected_workers,
    cross_patient_bootstrap_memory_safe_workers = cross_patient_result$runtime$memory_safe_workers,
    cross_patient_bootstrap_fallback_used = cross_patient_result$runtime$fallback_used,
    cross_patient_bootstrap_fallback_reason = cross_patient_result$runtime$fallback_reason,
    cross_patient_bootstrap_total_seconds = cross_patient_result$runtime$total_seconds,
    pooled_auc = as.numeric(roc_plot_all_iteration$auc),
    mean_iteration_auc = mean_auc,
    sd_iteration_auc = sd_auc,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  write.csv(
    run_summary,
    "Integrated_adaptive_scoring_results/run_summary.csv",
    row.names = FALSE
  )
  
  if (memory_mode == "full") {
    shp_S <- as.matrix(shp_S)
    shp_X <- as.data.frame(shp_X, check.names = FALSE)
    shp_testdata <- as.data.frame(shp_testdata, check.names = FALSE)
    shp_S[!is.finite(shp_S)] <- 0
    if (is.null(colnames(shp_S)) || any(colnames(shp_S) == "")) {
      colnames(shp_S) <- colnames(shp_X)
    }
    shp_X <- shp_X[, colnames(shp_S), drop = FALSE]
    retained_feature_names <- intersect(retained_feature_names, colnames(shp_S))
    if (length(retained_feature_names) == 0) {
      stop("No retained ranked features are present in the SHAP matrices for ", temp_split, ".")
    }
    shp_S <- shp_S[, retained_feature_names, drop = FALSE]
    shp_X <- shp_X[, retained_feature_names, drop = FALSE]
    beeswarm_source_S <- shp_S
    beeswarm_source_X <- shp_X
    beeswarm_baseline <- mean(shp_baseline)
    condition_keys <- sort(unique(as.character(shp_testdata$condition)))
    condition_aggregate_data <- setNames(vector("list", length(condition_keys)), condition_keys)
    for (condition_key in condition_keys) {
      condition_idx <- as.character(shp_testdata$condition) == condition_key
      condition_aggregate_data[[condition_key]] <- list(
        S = colMeans(beeswarm_source_S[condition_idx, , drop = FALSE], na.rm = TRUE),
        X = colMeans(beeswarm_source_X[condition_idx, , drop = FALSE], na.rm = TRUE)
      )
    }
    attr(condition_aggregate_data, "baseline") <- beeswarm_baseline
  } else {
    condition_aggregate_data <- build_condition_shap_aggregates_from_chunks(
      chunk_files = chunk_files,
      feature_names = retained_feature_names,
      condition_col = "condition"
    )
    beeswarm_plot_data <- collect_shap_plot_data_from_chunks(
      chunk_files = chunk_files,
      feature_names = unique(as.character(head(combined_results$feature, top_feature_show_actual)))
    )
    beeswarm_source_S <- as.matrix(beeswarm_plot_data$S)
    beeswarm_source_X <- as.data.frame(beeswarm_plot_data$X, check.names = FALSE)
    beeswarm_baseline <- beeswarm_plot_data$baseline
    rm(beeswarm_plot_data)
    gc(verbose = FALSE)
  }
  
  for (what_condition in sort(as.numeric(names(condition_aggregate_data)))) {
    condition_label <- ifelse(what_condition == 0, task_what_0, task_what_1)
    condition_summary <- condition_aggregate_data[[as.character(what_condition)]]
    condition_X_values <- condition_summary$X
    condition_X_values[!is.finite(condition_X_values)] <- 0
    write.csv(
      data.frame(average = round(condition_X_values, 2)),
      paste0("docs/average.exp_aggregated_condition.", condition_label, ".csv")
    )
  }
  
  saveRDS(list(fit_list, pred_list, confusionMatrix_all, ROC_all, model_performance_all, importance_data_all), file = paste0("docs/", sub_folder, "_all_fits.Rds"))
  
  density_features_to_save <- intersect(
    as.character(head(combined_results$feature, top_feature_show_actual)),
    colnames(dataset_task)
  )
  density_source <- dataset_task[, c("condition", density_features_to_save), drop = FALSE]
  
  write.csv(combined_results, "docs/combined_results.csv", row.names = FALSE)
  visualization_context <- list(
    temp_split = temp_split,
    comparison_mode = comparison_mode,
    split_by_col = ifelse(is.null(effective_split_by_col), NA_character_, effective_split_by_col),
    task_what_0 = task_what_0,
    task_what_1 = task_what_1,
    top_feature_label_actual = top_feature_label_actual,
    top_feature_show_actual = top_feature_show_actual,
    min_patients = min_patients_value,
    feature_min_log2FC = feature_min_log2FC,
    shap_expr_cor_cutoff = shap_expr_cor_cutoff,
    cross_patient_support_available = cross_patient_support_available,
    cross_patient_support_unavailable_reason = cross_patient_result$reason,
    cross_patient_support_note = cross_patient_support_note
  )
  save_shapture_visualization_inputs(
    output_dir = getwd(),
    combined_results = combined_results,
    ROC_all = ROC_all,
    AUC_all = AUC_all,
    roc_summary = list(
      pooled_auc = as.numeric(roc_plot_all_iteration$auc),
      mean_iteration_auc = mean_auc,
      sd_iteration_auc = sd_auc,
      ci_lower_iteration_auc = auc_ci[1],
      ci_upper_iteration_auc = auc_ci[2]
    ),
    model_performance_all = model_performance_all,
    importance_data_all = importance_data_all,
    condition_aggregate_data = condition_aggregate_data,
    beeswarm_source_S = beeswarm_source_S,
    beeswarm_source_X = beeswarm_source_X,
    beeswarm_baseline = beeswarm_baseline,
    density_source = density_source,
    context = visualization_context,
    group_colors = group_colors,
    continuous_colors = continuous_colors
  )
  
  rm(
    iter_results, shp_S, shp_X, shp_testdata, beeswarm_source_S,
    beeswarm_source_X, condition_aggregate_data, density_source,
    unit_shap_matrix, unit_feature_matrix, unit_shap_metadata,
    cross_patient_result
  )
  gc(verbose = FALSE)
  
  list(
    temp_split = temp_split,
    analyzed = TRUE,
    each_condition_sample_no = each.condition_total_sample_No,
    min_n_cells = min_n_cells,
    output_dir = file.path(main_dir, sub_folder)
  )
}


#' Run the SHAPTURE pipeline
#'
#' @param dataset_all Input data frame containing `condition`, `sample_id`, any
#'   split or one-vs-others columns, and numeric features.
#' @param what_genes Feature list to analyze.
#' @param main_dir Output directory.
#' @param split_by_col Optional column used for subset-specific condition
#'   comparison. Use `NULL` for no split in `condition_vs_condition` mode. In
#'   `one_vs_others` mode, this should be the column whose levels are compared
#'   one-vs-others.
#' @param subset_levels Optional `split_by_col` levels to keep in the input data
#'   before task construction.
#' @param split_levels Optional task or target levels to analyze after
#'   `subset_levels` filtering.
#' @param comparison_mode One of `"condition_vs_condition"` or
#'   `"one_vs_others"`.
#' @param one_vs_others_col Column used for one-vs-others comparisons.
#' @param one_vs_others_levels Optional target levels for one-vs-others mode.
#' @param one_vs_others_negative_label Label for the pooled negative class in
#'   one-vs-others mode.
#' @param what_0 Negative class label for condition-vs-condition mode.
#' @param what_1 Positive class label for condition-vs-condition mode.
#' @param project_name Optional output prefix.
#' @param train_ratio Training fraction.
#' @param test_ratio Test fraction.
#' @param iteration_no Number of repeated iterations.
#' @param gene_fraction_per_iteration Fraction of task features included in
#'   each coverage-aware model iteration. The default is `0.30`.
#' @param cross_patient_bootstraps Number of patient-block bootstrap replicates
#'   used for Cross-patient Support. The default is `500L`.
#' @param analysis_seed Positive integer reproducibility root. The default
#'   `123L` deterministically controls task sampling, feature selection, model
#'   fitting, metacell construction, and Cross-patient bootstrapping.
#' @param top_feature_label Number of labeled features.
#' @param top_feature_show Number of top features shown in plots and saved
#'   for post-analysis visualization.
#' @param if_adaptive_weights Whether to use adaptive weighting.
#' @param default_prior_scores_weights Named prior weights for the five final
#'   score components.
#' @param feature_min_pct Minimum fraction of cells with feature value greater
#'   than zero in at least one task group. Set to `0` to disable this filter.
#' @param feature_min_log2FC Minimum absolute expression log2FC in original
#'   single-cell task data. Set to `0` to disable this filter.
#' @param shap_expr_cor_cutoff Minimum absolute value-SHAP correlation needed
#'   to assign a model-expression pattern label.
#' @param score_scaling_method Backward-compatible configuration field. The
#'   default `"median_anchor"` uses fixed task-level median anchors for
#'   unbounded evidence and preserves naturally bounded 0-1 evidence. Legacy
#'   `"winsorized"` and `"minmax"` values remain accepted for existing scripts
#'   but do not replace the evidence-score strategy.
#' @param score_scaling_probs Two quantiles used only when winsorizing
#'   component values for optional adaptive-weight estimation.
#' @param preprocessing_method Either `"none"` or `"scale"`.
#' @param model_type Model backend.
#' @param analysis_unit_type Either `"single_cell"` or `"metacell"`.
#' @param memory_mode Either `"full"` or `"safe"`.
#' @param patient_balance_method Patient balancing method.
#' @param patient_split_strategy `"shared_patients"` keeps the existing
#'   coverage-aware cell/metacell split and permits patient overlap between
#'   training and testing. `"unseen_patients"` uses globally disjoint patient
#'   folds with a rounded 60%/40% train/test allocation and task-level fallback.
#' @param shap_expr_cor_direction_filter Direction filter for value-to-SHAP
#'   correlations.
#' @param group_colors Two unnamed colors used for left/right group displays.
#' @param continuous_colors Continuous color palette used for SHAP/expression
#'   gradients.
#' @param n_cores Requested maximum workers for model fitting and Cross-patient
#'   bootstrapping. Runtime CPU and RAM safeguards may reduce the bootstrap
#'   worker count or select serial bootstrapping.
#' @param metacell_min_cells_to_build Minimum rows before metacell building.
#' @param metacell_min_cells_per_metacell Minimum rows per metacell.
#' @param metacell_max_cells_per_metacell Maximum rows per metacell.
#' @param metacell_n_variable_features Number of variable features for metacells.
#' @param metacell_n_pcs Number of PCs used during metacell construction.
#' @param metacell_aggregation_fun Aggregation function.
#' @param parallel_metacells Whether to parallelize metacell building.
#' @param metacell_n_cores Number of metacell worker cores.
#' @param metacell_parallel_backend Metacell parallel backend.
#' @param ... Reserved for compatibility with older scripts.
#'
#' @return Invisibly returns a list with task summaries and output folders.
#' @export
run_shapture_pipeline <- function(
    dataset_all,
    what_genes,
    main_dir,
    split_by_col = NULL,
    subset_levels = NULL,
    split_levels = NULL,
    comparison_mode = c("condition_vs_condition", "one_vs_others"),
    one_vs_others_col = NULL,
    one_vs_others_levels = NULL,
    one_vs_others_negative_label = "others",
    what_0 = "low",
    what_1 = "high",
    project_name = NULL,
    train_ratio = 0.7,
    test_ratio = 0.3,
    iteration_no = 50,
    gene_fraction_per_iteration = 0.30,
    cross_patient_bootstraps = 500L,
    analysis_seed = 123L,
    top_feature_label = 30,
    top_feature_show = length(what_genes),
    if_adaptive_weights = FALSE,
    feature_min_pct = 0.25,
    feature_min_log2FC = 0.25,
    shap_expr_cor_cutoff = 0.25,
    score_scaling_method = c("median_anchor", "winsorized", "minmax"),
    score_scaling_probs = c(0.05, 0.95),
    default_prior_scores_weights = c(
      SHAP_importance_prior_weight = 0.25,
      expression_signal_prior_weight = 0.25,
      value_SHAP_relationship_prior_weight = 0.30,
      cross_patient_support_prior_weight = 0.10,
      model_reproducibility_prior_weight = 0.10
    ),
    preprocessing_method = "scale",
    model_type = "lightgbm",
    analysis_unit_type = "metacell",
    memory_mode = "full",
    patient_balance_method = "model_weights",
    patient_split_strategy = "shared_patients",
    shap_expr_cor_direction_filter = "both",
    group_colors = c("darkblue", "darkred"),
    continuous_colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
    n_cores = 14,
    metacell_min_cells_to_build = 40,
    metacell_min_cells_per_metacell = 5,
    metacell_max_cells_per_metacell = 25,
    metacell_n_variable_features = 2000,
    metacell_n_pcs = 20,
    metacell_aggregation_fun = "mean",
    parallel_metacells = TRUE,
    metacell_n_cores = 10,
    metacell_parallel_backend = "psock",
    ...
) {
  legacy_args <- list(...)
  legacy_names <- c("one_vs_rest_col", "one_vs_rest_levels", "one_vs_rest_negative_label")
  unknown_args <- setdiff(names(legacy_args), legacy_names)
  if (length(unknown_args) > 0) {
    stop("Unused arguments: ", paste(unknown_args, collapse = ", "))
  }
  if (!is.null(legacy_args$one_vs_rest_col) && is.null(one_vs_others_col)) {
    one_vs_others_col <- legacy_args$one_vs_rest_col
  }
  if (!is.null(legacy_args$one_vs_rest_levels) && is.null(one_vs_others_levels)) {
    one_vs_others_levels <- legacy_args$one_vs_rest_levels
  }
  if (!is.null(legacy_args$one_vs_rest_negative_label) &&
      identical(one_vs_others_negative_label, "others")) {
    one_vs_others_negative_label <- legacy_args$one_vs_rest_negative_label
  }

  split_by_col <- resolve_mapped_prepared_column(
    column_name = split_by_col,
    dataset_all = dataset_all
  )
  one_vs_others_col <- resolve_mapped_prepared_column(
    column_name = one_vs_others_col,
    dataset_all = dataset_all
  )
  
  comparison_mode <- normalize_comparison_mode(comparison_mode)
  comparison_mode <- match.arg(comparison_mode, c("condition_vs_condition", "one_vs_others"))
  memory_mode <- match.arg(memory_mode, c("full", "safe"))
  patient_split_strategy <- match.arg(
    patient_split_strategy,
    c("shared_patients", "unseen_patients")
  )
  analysis_seed <- validate_analysis_seed(analysis_seed)
  score_scaling_method <- match.arg(score_scaling_method)
  if (!is.numeric(gene_fraction_per_iteration) ||
      length(gene_fraction_per_iteration) != 1L ||
      !is.finite(gene_fraction_per_iteration) ||
      gene_fraction_per_iteration <= 0 ||
      gene_fraction_per_iteration > 1) {
    stop("gene_fraction_per_iteration must be one finite value greater than 0 and at most 1.")
  }
  if (!is.numeric(cross_patient_bootstraps) ||
      length(cross_patient_bootstraps) != 1L ||
      !is.finite(cross_patient_bootstraps) ||
      cross_patient_bootstraps < 1 ||
      cross_patient_bootstraps != floor(cross_patient_bootstraps)) {
    stop("cross_patient_bootstraps must be one positive whole number.")
  }
  cross_patient_bootstraps <- as.integer(cross_patient_bootstraps)
  if (!is.numeric(n_cores) || length(n_cores) != 1L ||
      !is.finite(n_cores) || n_cores < 1L || n_cores != floor(n_cores)) {
    stop("n_cores must be one positive whole number.")
  }
  n_cores <- as.integer(n_cores)
  if (!is.numeric(feature_min_pct) ||
      length(feature_min_pct) != 1 ||
      is.na(feature_min_pct) ||
      feature_min_pct < 0 ||
      feature_min_pct >= 1) {
    stop("feature_min_pct must be one numeric value >= 0 and < 1.")
  }
  if (!is.numeric(feature_min_log2FC) ||
      length(feature_min_log2FC) != 1 ||
      is.na(feature_min_log2FC) ||
      feature_min_log2FC < 0) {
    stop("feature_min_log2FC must be one numeric value >= 0.")
  }
  if (!is.numeric(shap_expr_cor_cutoff) ||
      length(shap_expr_cor_cutoff) != 1 ||
      is.na(shap_expr_cor_cutoff) ||
      shap_expr_cor_cutoff < 0 ||
      shap_expr_cor_cutoff > 1) {
    stop("shap_expr_cor_cutoff must be one numeric value between 0 and 1.")
  }
  if (!is.numeric(score_scaling_probs) ||
      length(score_scaling_probs) != 2 ||
      any(!is.finite(score_scaling_probs)) ||
      any(score_scaling_probs < 0 | score_scaling_probs > 1) ||
      score_scaling_probs[1] >= score_scaling_probs[2]) {
    stop("score_scaling_probs must be two increasing numeric probabilities between 0 and 1.")
  }
  default_prior_scores_weights <- normalize_score_weights(default_prior_scores_weights)
  load_shapture_runtime()

  if (comparison_mode == "one_vs_others") {
    one_vs_others_col <- one_vs_others_col %||% split_by_col
    one_vs_others_levels <- one_vs_others_levels %||% split_levels
  }
  
  if (!is.null(subset_levels)) {
    if (is.null(split_by_col)) {
      stop("subset_levels requires split_by_col.")
    }
    if (!split_by_col %in% colnames(dataset_all)) {
      stop("split_by_col is not present in dataset_all: ", split_by_col)
    }
    dataset_all <- dataset_all[
      as.character(dataset_all[[split_by_col]]) %in% as.character(subset_levels),
      ,
      drop = FALSE
    ]
    if (nrow(dataset_all) == 0) {
      stop("No rows remain after applying subset_levels.")
    }
  }
  
  if (!all(c("condition", "sample_id") %in% colnames(dataset_all))) {
    stop("dataset_all must contain at least 'condition' and 'sample_id'.")
  }
  
  if (is.null(project_name)) {
    project_name <- if (comparison_mode == "condition_vs_condition") {
      paste0(what_0, "_", what_1)
    } else {
      "one_vs_others"
    }
  }

  resolved_config_for_docs <- list(
    split_by_col = split_by_col,
    subset_levels = subset_levels,
    split_levels = split_levels,
    n_rows_after_subset_levels = nrow(dataset_all),
    comparison_mode = comparison_mode,
    one_vs_others_col = one_vs_others_col,
    one_vs_others_levels = one_vs_others_levels,
    one_vs_others_negative_label = one_vs_others_negative_label,
    what_0 = what_0,
    what_1 = what_1,
    project_name = project_name,
    train_ratio = train_ratio,
    test_ratio = test_ratio,
    iteration_no = iteration_no,
    gene_fraction_per_iteration = gene_fraction_per_iteration,
    cross_patient_bootstraps = cross_patient_bootstraps,
    analysis_seed = analysis_seed,
    top_feature_label = top_feature_label,
    top_feature_show = top_feature_show,
    if_adaptive_weights = if_adaptive_weights,
    feature_min_pct = feature_min_pct,
    feature_min_log2FC = feature_min_log2FC,
    shap_expr_cor_cutoff = shap_expr_cor_cutoff,
    score_scaling_method = score_scaling_method,
    score_scaling_probs = score_scaling_probs,
    default_prior_scores_weights = default_prior_scores_weights,
    preprocessing_method = preprocessing_method,
    model_type = model_type,
    analysis_unit_type = analysis_unit_type,
    memory_mode = memory_mode,
    patient_balance_method = patient_balance_method,
    patient_split_strategy = patient_split_strategy,
    shap_expr_cor_direction_filter = shap_expr_cor_direction_filter,
    group_colors = group_colors,
    continuous_colors = continuous_colors,
    n_cores = n_cores,
    metacell_min_cells_to_build = metacell_min_cells_to_build,
    metacell_min_cells_per_metacell = metacell_min_cells_per_metacell,
    metacell_max_cells_per_metacell = metacell_max_cells_per_metacell,
    metacell_n_variable_features = metacell_n_variable_features,
    metacell_n_pcs = metacell_n_pcs,
    metacell_aggregation_fun = metacell_aggregation_fun,
    parallel_metacells = parallel_metacells,
    metacell_n_cores = metacell_n_cores,
    metacell_parallel_backend = metacell_parallel_backend
  )
  
  feature_cols <- intersect(what_genes, colnames(dataset_all))
  if (length(feature_cols) == 0) {
    stop("No requested features from what_genes are present in dataset_all.")
  }
  what_genes <- feature_cols
  
  task_context <- build_shapture_tasks(
    dataset_all = dataset_all,
    comparison_mode = comparison_mode,
    split_by_col = split_by_col,
    split_levels = split_levels,
    one_vs_others_col = one_vs_others_col,
    one_vs_others_levels = one_vs_others_levels,
    what_0 = what_0,
    what_1 = what_1,
    one_vs_others_negative_label = one_vs_others_negative_label
  )
  
  dataset_all_work <- task_context$dataset
  metacell_split_col <- task_context$task_metacell_split_col
  task_table <- task_context$task_table
  
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  
  foreach::registerDoSEQ()
  
  dataset_all_metacell <- dataset_all_work
  metacell_seed <- if (analysis_unit_type == "metacell") {
    derive_shapture_seed(
      analysis_seed,
      "metacell_construction",
      comparison_mode,
      split_by_col,
      sort(as.character(subset_levels %||% "<ALL>"))
    )
  } else {
    NA_integer_
  }
  if (analysis_unit_type == "metacell") {
    metacell_res <- build_grouped_metacells(
      dataset = dataset_all_metacell,
      feature_cols = what_genes,
      sample_col = "sample_id",
      split_by_col = metacell_split_col,
      condition_col = "condition",
      use_metacells = TRUE,
      min_cells_to_build = metacell_min_cells_to_build,
      min_cells_per_metacell = metacell_min_cells_per_metacell,
      max_cells_per_metacell = metacell_max_cells_per_metacell,
      n_variable_features = metacell_n_variable_features,
      n_pcs = metacell_n_pcs,
      aggregation_fun = metacell_aggregation_fun,
      scale_features = TRUE,
      seed = metacell_seed,
      verbose = TRUE,
      parallel = parallel_metacells,
      n_cores = metacell_n_cores,
      parallel_backend = metacell_parallel_backend
    )
    
    dataset_all_analysis <- metacell_res$data
    metacell_summary <- metacell_res$summary
  } else {
    dataset_all_analysis <- dataset_all_work
    metacell_summary <- NULL
  }
  
  setwd(main_dir)
  if (comparison_mode == "condition_vs_condition" && !is.null(split_by_col)) {
    write.csv(
      t(table(dataset_all_analysis$condition, dataset_all_analysis[[split_by_col]])),
      "cell.number.split.condition_before.loop.csv"
    )
  } else if (comparison_mode == "condition_vs_condition") {
    write.csv(
      data.frame(condition = names(table(dataset_all_analysis$condition)), n = as.integer(table(dataset_all_analysis$condition))),
      "cell.number.by.condition_before.loop.csv",
      row.names = FALSE
    )
  } else {
    write.csv(
      data.frame(split_level = names(table(dataset_all_analysis[[metacell_split_col]])), n = as.integer(table(dataset_all_analysis[[metacell_split_col]]))),
      "cell.number.by.one_vs_others_group_before.loop.csv",
      row.names = FALSE
    )
  }
  
  task_results <- vector("list", nrow(task_table))
  
  for (i in seq_len(nrow(task_table))) {
    cat(paste0(date(), "\n"))
    task_results[[i]] <- run_shapture_task(
      dataset_all = dataset_all_work,
      dataset_all_analysis = dataset_all_analysis,
      what_genes = what_genes,
      temp_split = task_table$temp_split[i],
      task_what_0 = task_table$task_what_0[i],
      task_what_1 = task_table$task_what_1[i],
      main_dir = main_dir,
      project_name = project_name,
      comparison_mode = comparison_mode,
      effective_split_by_col = task_context$effective_split_by_col,
      metacell_split_col = metacell_split_col,
      one_vs_others_col = one_vs_others_col,
      metacell_summary = metacell_summary,
      continuous_colors = continuous_colors,
      group_colors = group_colors,
      train_ratio = train_ratio,
      test_ratio = test_ratio,
      iteration_no = iteration_no,
      gene_fraction_per_iteration = gene_fraction_per_iteration,
      cross_patient_bootstraps = cross_patient_bootstraps,
      n_cores = n_cores,
      analysis_seed = analysis_seed,
      metacell_seed = metacell_seed,
      top_feature_label = top_feature_label,
      top_feature_show = top_feature_show,
      if_adaptive_weights = if_adaptive_weights,
      feature_min_pct = feature_min_pct,
      feature_min_log2FC = feature_min_log2FC,
      shap_expr_cor_cutoff = shap_expr_cor_cutoff,
      score_scaling_method = score_scaling_method,
      score_scaling_probs = score_scaling_probs,
      preprocessing_method = preprocessing_method,
      model_type = model_type,
      analysis_unit_type = analysis_unit_type,
      memory_mode = memory_mode,
      patient_balance_method = patient_balance_method,
      patient_split_strategy = patient_split_strategy,
      shap_expr_cor_direction_filter = shap_expr_cor_direction_filter,
      default_prior_scores_weights = default_prior_scores_weights,
      resolved_config_for_docs = resolved_config_for_docs
    )
  }
  
  cell_number_split <- do.call(
    rbind,
    lapply(task_results, function(x) {
      data.frame(
        split_level = x$temp_split,
        each_condition_sample_No = x$each_condition_sample_no,
        analyzed = ifelse(isTRUE(x$analyzed), "Yes", "No"),
        stringsAsFactors = FALSE
      )
    })
  )
  write.csv(cell_number_split, "cell.used.per.condition.csv", row.names = FALSE)
  
  invisible(list(
    task_results = task_results,
    output_dirs = vapply(task_results, function(x) ifelse(isTRUE(x$analyzed), x$output_dir, NA_character_), character(1))
  ))
}
