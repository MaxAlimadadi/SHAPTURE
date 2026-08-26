#' Compute the patient-balanced SHAP effect size
#'
#' @param shap_matrix Unit-averaged SHAP matrix.
#' @param metadata Metadata aligned to `shap_matrix`.
#' @param sample_col Sample identifier column.
#' @param condition_col Binary condition column.
#'
#' @return Data frame with one row per feature and a patient-balanced signed
#'   SHAP difference (`condition_1 - condition_0`).
#' @keywords internal
compute_patient_balanced_shap_effect <- function(
    shap_matrix,
    metadata,
    sample_col = "sample_id",
    condition_col = "condition"
) {
  shap_matrix <- as.matrix(shap_matrix)
  metadata <- as.data.frame(metadata, check.names = FALSE)
  if (nrow(shap_matrix) != nrow(metadata)) {
    stop("shap_matrix and metadata must have the same number of rows.")
  }
  if (!all(c(sample_col, condition_col) %in% colnames(metadata))) {
    stop("metadata must contain sample and condition columns.")
  }

  condition_values <- as.character(metadata[[condition_col]])
  if (!all(stats::na.omit(unique(condition_values)) %in% c("0", "1"))) {
    stop("Cross-patient calculations require condition values coded as 0 and 1.")
  }
  valid_rows <- !is.na(metadata[[sample_col]]) & !is.na(condition_values)
  patient_keys <- paste(
    as.character(metadata[[sample_col]][valid_rows]),
    condition_values[valid_rows],
    sep = "___"
  )
  patient_means <- rowsum(
    shap_matrix[valid_rows, , drop = FALSE],
    group = patient_keys,
    reorder = FALSE,
    na.rm = TRUE
  )
  observed_counts <- rowsum(
    is.finite(shap_matrix[valid_rows, , drop = FALSE]) * 1,
    group = patient_keys,
    reorder = FALSE
  )
  patient_means[observed_counts > 0] <-
    patient_means[observed_counts > 0] / observed_counts[observed_counts > 0]
  patient_means[observed_counts == 0] <- NA_real_
  patient_conditions <- sub("^.*___", "", rownames(patient_means))

  effect_size <- colMeans(
    patient_means[patient_conditions == "1", , drop = FALSE],
    na.rm = TRUE
  ) - colMeans(
    patient_means[patient_conditions == "0", , drop = FALSE],
    na.rm = TRUE
  )

  data.frame(
    feature = colnames(shap_matrix),
    effect_size = effect_size,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


#' Calculate conservative bootstrap retention
#'
#' @param values Numeric bootstrap values.
#' @param zero_tolerance Values with a median at or below this threshold have no
#'   measurable signal and receive support zero.
#'
#' @return Named numeric vector containing support, median, fifth percentile,
#'   and number of valid bootstrap values.
#' @keywords internal
compute_bootstrap_retention <- function(values, zero_tolerance = 1e-8) {
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(c(support = 0, median = NA_real_, q05 = NA_real_, n_valid = 0))
  }
  bootstrap_median <- stats::median(values)
  bootstrap_q05 <- as.numeric(stats::quantile(
    values,
    probs = 0.05,
    names = FALSE,
    na.rm = TRUE
  ))
  support <- if (!is.finite(bootstrap_median) || bootstrap_median <= zero_tolerance) {
    0
  } else {
    min(1, max(0, bootstrap_q05 / bootstrap_median))
  }
  c(
    support = support,
    median = bootstrap_median,
    q05 = bootstrap_q05,
    n_valid = length(values)
  )
}


#' Calculate smoother SHAP-importance bootstrap retention
#'
#' Combines the frequency of positive scaled importance with conservative
#' magnitude retention among positive bootstrap values. This avoids forcing
#' support to zero solely because the lower tail contains some zero values.
#'
#' @param values Numeric bootstrap importance values on a nonnegative scale.
#' @param lower_probability Lower quantile used among positive values.
#' @param zero_tolerance Values at or below this threshold are treated as zero.
#'
#' @return Named numeric vector containing support, nonzero frequency,
#'   positive-value retention diagnostics, and valid bootstrap counts.
#' @keywords internal
compute_smoother_importance_retention <- function(
    values,
    lower_probability = 0.05,
    zero_tolerance = 1e-8
) {
  if (!is.numeric(lower_probability) ||
      length(lower_probability) != 1L ||
      !is.finite(lower_probability) ||
      lower_probability <= 0 ||
      lower_probability >= 0.5) {
    stop("lower_probability must be one finite value greater than 0 and less than 0.5.")
  }

  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(c(
      support = 0,
      nonzero_frequency = 0,
      positive_median = NA_real_,
      positive_lower = NA_real_,
      positive_retention = 0,
      n_valid = 0,
      n_positive = 0
    ))
  }
  if (any(values < 0)) {
    stop("Smoother SHAP-importance retention requires nonnegative bootstrap values.")
  }

  positive_values <- values[values > zero_tolerance]
  n_valid <- length(values)
  n_positive <- length(positive_values)
  nonzero_frequency <- n_positive / n_valid
  if (n_positive == 0L) {
    return(c(
      support = 0,
      nonzero_frequency = nonzero_frequency,
      positive_median = NA_real_,
      positive_lower = NA_real_,
      positive_retention = 0,
      n_valid = n_valid,
      n_positive = n_positive
    ))
  }

  positive_median <- stats::median(positive_values)
  positive_lower <- as.numeric(stats::quantile(
    positive_values,
    probs = lower_probability,
    names = FALSE,
    na.rm = TRUE,
    type = 7
  ))
  positive_retention <- if (!is.finite(positive_median) ||
      positive_median <= zero_tolerance) {
    0
  } else {
    min(1, max(0, positive_lower / positive_median))
  }

  c(
    support = nonzero_frequency * positive_retention,
    nonzero_frequency = nonzero_frequency,
    positive_median = positive_median,
    positive_lower = positive_lower,
    positive_retention = positive_retention,
    n_valid = n_valid,
    n_positive = n_positive
  )
}


#' Compute Model reproducibility across repeated model fits
#'
#' Summarizes how consistently a feature retains raw mean absolute SHAP
#' importance, signed
#' condition-level SHAP separation, and its value-SHAP relationship across
#' feature-bagged model iterations.
#'
#' @param iteration_data Long data frame containing one row per included
#'   feature and model iteration.
#' @param full_stats Full-analysis feature statistics containing `feature`,
#'   `effect_size`, and `shap_expr_cor`.
#' @param shap_expr_cor_cutoff Minimum absolute per-iteration value-SHAP
#'   correlation counted as retained.
#' @param lower_probability Lower quantile used for magnitude retention.
#' @param zero_tolerance Positive values at or below this tolerance are treated
#'   as zero importance.
#'
#' @return Data frame with one bounded Model reproducibility score and its
#'   three subcomponent diagnostics per feature.
#' @keywords internal
compute_model_reproducibility <- function(
    iteration_data,
    full_stats,
    shap_expr_cor_cutoff,
    lower_probability = 0.10,
    zero_tolerance = 1e-8
) {
  required_iteration_columns <- c(
    "iteration", "feature", "shap_importance", "relative_delta_shap",
    "delta_shap", "value_shap_cor"
  )
  missing_iteration_columns <- setdiff(
    required_iteration_columns,
    colnames(iteration_data)
  )
  if (length(missing_iteration_columns) > 0L) {
    stop(
      "Missing Model reproducibility iteration columns: ",
      paste(missing_iteration_columns, collapse = ", ")
    )
  }
  required_full_columns <- c("feature", "effect_size", "shap_expr_cor")
  missing_full_columns <- setdiff(required_full_columns, colnames(full_stats))
  if (length(missing_full_columns) > 0L) {
    stop(
      "Missing full-analysis Model reproducibility columns: ",
      paste(missing_full_columns, collapse = ", ")
    )
  }
  if (!is.numeric(shap_expr_cor_cutoff) ||
      length(shap_expr_cor_cutoff) != 1L ||
      !is.finite(shap_expr_cor_cutoff) ||
      shap_expr_cor_cutoff < 0 ||
      shap_expr_cor_cutoff > 1) {
    stop("shap_expr_cor_cutoff must be one finite value between 0 and 1.")
  }
  if (!is.numeric(lower_probability) ||
      length(lower_probability) != 1L ||
      !is.finite(lower_probability) ||
      lower_probability <= 0 ||
      lower_probability >= 0.5) {
    stop("lower_probability must be one finite value greater than 0 and less than 0.5.")
  }

  clamp01 <- function(x) pmin(1, pmax(0, x))

  magnitude_retention <- function(values) {
    values <- values[is.finite(values) & values >= 0]
    if (length(values) == 0L) {
      return(c(support = 0, median = NA_real_, lower = NA_real_))
    }
    value_median <- stats::median(values)
    value_lower <- as.numeric(stats::quantile(
      values,
      probs = lower_probability,
      names = FALSE,
      na.rm = TRUE,
      type = 7
    ))
    support <- if (!is.finite(value_median) || value_median <= zero_tolerance) {
      0
    } else {
      clamp01(value_lower / value_median)
    }
    c(support = support, median = value_median, lower = value_lower)
  }

  smoother_importance_retention <- function(values) {
    values <- values[is.finite(values) & values >= 0]
    n_valid <- length(values)
    positive_values <- values[values > zero_tolerance]
    n_positive <- length(positive_values)
    nonzero_frequency <- if (n_valid > 0L) n_positive / n_valid else 0
    positive_retention <- magnitude_retention(positive_values)
    c(
      support = clamp01(nonzero_frequency * positive_retention["support"]),
      nonzero_frequency = nonzero_frequency,
      positive_median = positive_retention["median"],
      positive_lower = positive_retention["lower"],
      positive_retention = positive_retention["support"],
      n_positive = n_positive
    )
  }

  split_iterations <- split(
    iteration_data,
    as.character(iteration_data$feature),
    drop = TRUE
  )
  feature_results <- lapply(seq_len(nrow(full_stats)), function(feature_index) {
    feature_name <- as.character(full_stats$feature[feature_index])
    feature_iterations <- split_iterations[[feature_name]]
    if (is.null(feature_iterations) || nrow(feature_iterations) == 0L) {
      return(data.frame(
        feature = feature_name,
        model_reproducibility = 0,
        model_shap_importance_retention = 0,
        model_shap_importance_nonzero_frequency = 0,
        model_shap_importance_positive_median = NA_real_,
        model_shap_importance_positive_q10 = NA_real_,
        model_shap_importance_positive_retention = 0,
        positive_importance_iterations = 0L,
        delta_shap_direction_retention = 0,
        delta_shap_magnitude_retention = 0,
        delta_shap_magnitude_median = NA_real_,
        delta_shap_magnitude_q10 = NA_real_,
        delta_shap_retention = 0,
        value_shap_retention = 0,
        valid_model_iterations = 0L,
        direction_matching_iterations = 0L,
        valid_delta_shap_iterations = 0L,
        value_shap_supporting_iterations = 0L,
        valid_value_shap_iterations = 0L,
        invalid_value_shap_iterations = 0L,
        full_effect_direction = sign(full_stats$effect_size[feature_index]),
        full_value_shap_direction = sign(full_stats$shap_expr_cor[feature_index]),
        stringsAsFactors = FALSE
      ))
    }

    importance_retention <- smoother_importance_retention(
      feature_iterations$shap_importance
    )

    reference_delta_sign <- sign(full_stats$effect_size[feature_index])
    valid_delta <- is.finite(feature_iterations$delta_shap)
    delta_matches <- valid_delta &
      reference_delta_sign != 0 &
      sign(feature_iterations$delta_shap) == reference_delta_sign
    delta_direction_retention <- if (any(valid_delta)) {
      sum(delta_matches) / sum(valid_delta)
    } else {
      0
    }
    delta_magnitude_retention <- magnitude_retention(
      feature_iterations$relative_delta_shap[delta_matches]
    )
    delta_retention <- clamp01(
      delta_direction_retention * delta_magnitude_retention["support"]
    )

    reference_cor_sign <- sign(full_stats$shap_expr_cor[feature_index])
    valid_cor <- is.finite(feature_iterations$value_shap_cor)
    value_matches <- valid_cor &
      reference_cor_sign != 0 &
      sign(feature_iterations$value_shap_cor) == reference_cor_sign &
      abs(feature_iterations$value_shap_cor) >= shap_expr_cor_cutoff
    value_retention <- if (any(valid_cor)) {
      sum(value_matches) / sum(valid_cor)
    } else {
      0
    }

    model_reproducibility <- mean(c(
      importance_retention["support"],
      delta_retention,
      value_retention
    ))

    data.frame(
      feature = feature_name,
      model_reproducibility = unname(clamp01(model_reproducibility)),
      model_shap_importance_retention = unname(importance_retention["support"]),
      model_shap_importance_nonzero_frequency = unname(importance_retention["nonzero_frequency"]),
      model_shap_importance_positive_median = unname(importance_retention["positive_median"]),
      model_shap_importance_positive_q10 = unname(importance_retention["positive_lower"]),
      model_shap_importance_positive_retention = unname(importance_retention["positive_retention"]),
      positive_importance_iterations = as.integer(importance_retention["n_positive"]),
      delta_shap_direction_retention = unname(delta_direction_retention),
      delta_shap_magnitude_retention = unname(delta_magnitude_retention["support"]),
      delta_shap_magnitude_median = unname(delta_magnitude_retention["median"]),
      delta_shap_magnitude_q10 = unname(delta_magnitude_retention["lower"]),
      delta_shap_retention = unname(delta_retention),
      value_shap_retention = unname(clamp01(value_retention)),
      valid_model_iterations = nrow(feature_iterations),
      direction_matching_iterations = sum(delta_matches),
      valid_delta_shap_iterations = sum(valid_delta),
      value_shap_supporting_iterations = sum(value_matches),
      valid_value_shap_iterations = sum(valid_cor),
      invalid_value_shap_iterations = sum(!valid_cor),
      full_effect_direction = reference_delta_sign,
      full_value_shap_direction = reference_cor_sign,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  do.call(rbind, feature_results)
}


#' Compute Cross-patient support by patient-block bootstrap
#'
#' Patient blocks are resampled with replacement. SHAP importance, Expression
#' signal, and the Value-SHAP relationship are recalculated for the fixed final
#' feature set. SHAP importance combines nonzero frequency with positive-value
#' fifth-percentile-to-median retention. Expression signal uses the original
#' all-bootstrap fifth-percentile-to-median retention, and the Value-SHAP
#' relationship uses direction-and-threshold retention frequency.
#'
#' @param shap_matrix Unit-averaged held-out SHAP matrix.
#' @param shap_feature_values Unit-averaged model-input feature matrix.
#' @param shap_metadata Metadata aligned to the SHAP matrices.
#' @param model_expression_data Model-unit data used for full-data log2FC.
#' @param single_cell_expression_data Original-cell data used for detection and
#'   within-condition consistency.
#' @param full_stats Full-data feature statistics.
#' @param average_importance Full-data mean absolute SHAP importance table.
#' @param feature_names Fixed retained feature list.
#' @param comparison_mode Either `"condition_vs_condition"` or
#'   `"one_vs_others"`.
#' @param feature_min_log2FC Minimum bootstrap log2FC magnitude.
#' @param shap_expr_cor_cutoff Minimum bootstrap Value-SHAP correlation.
#' @param score_scaling_method Full-analysis component scaling method.
#' @param score_scaling_probs Full-analysis winsorization probabilities.
#' @param n_bootstraps Number of patient-block bootstrap replicates.
#' @param seed Random seed.
#'
#' @return List containing per-feature support, diagnostics, and task metadata.
#' @keywords internal
compute_cross_patient_support <- function(
    shap_matrix,
    shap_feature_values,
    shap_metadata,
    model_expression_data,
    single_cell_expression_data,
    full_stats,
    average_importance,
    feature_names,
    comparison_mode = c("condition_vs_condition", "one_vs_others"),
    feature_min_log2FC = 0.25,
    shap_expr_cor_cutoff = 0.25,
    score_scaling_method = c("winsorized", "minmax"),
    score_scaling_probs = c(0.05, 0.95),
    n_bootstraps = 1000L,
    seed = 42L
) {
  comparison_mode <- match.arg(comparison_mode)
  score_scaling_method <- match.arg(score_scaling_method)
  n_bootstraps <- as.integer(n_bootstraps)
  if (!is.finite(n_bootstraps) || n_bootstraps < 1L) {
    stop("n_bootstraps must be a positive integer.")
  }

  feature_names <- unique(as.character(feature_names))
  if (length(feature_names) == 0) {
    stop("feature_names cannot be empty for Cross-patient support.")
  }
  shap_matrix <- as.matrix(shap_matrix[, feature_names, drop = FALSE])
  shap_feature_values <- as.matrix(shap_feature_values[, feature_names, drop = FALSE])
  shap_metadata <- as.data.frame(shap_metadata, check.names = FALSE)
  model_expression_data <- as.data.frame(model_expression_data, check.names = FALSE)
  single_cell_expression_data <- as.data.frame(single_cell_expression_data, check.names = FALSE)

  if (nrow(shap_matrix) != nrow(shap_feature_values) ||
      nrow(shap_matrix) != nrow(shap_metadata)) {
    stop("SHAP values, feature values, and SHAP metadata must align by row.")
  }
  required_metadata <- c("sample_id", "condition")
  data_sources <- list(
    shap = shap_metadata,
    model = model_expression_data,
    single_cell = single_cell_expression_data
  )
  for (source_name in names(data_sources)) {
    source_data <- data_sources[[source_name]]
    required_columns <- if (identical(source_name, "shap")) {
      required_metadata
    } else {
      c(required_metadata, feature_names)
    }
    missing_columns <- setdiff(required_columns, colnames(source_data))
    if (length(missing_columns) > 0) {
      stop(
        "Missing Cross-patient support columns in ", source_name, " data: ",
        paste(head(missing_columns, 10), collapse = ", ")
      )
    }
    source_conditions <- as.character(source_data$condition)
    if (!all(stats::na.omit(unique(source_conditions)) %in% c("0", "1"))) {
      stop(source_name, " data must use condition values 0 and 1.")
    }
  }

  full_match <- match(feature_names, full_stats$feature)
  importance_match <- match(feature_names, average_importance$feature)
  if (anyNA(full_match) || anyNA(importance_match)) {
    stop("Cross-patient support inputs do not cover every retained feature.")
  }
  full_stats <- full_stats[full_match, , drop = FALSE]
  full_importance <- average_importance$`average|SHAP value|`[importance_match]

  make_block_table <- function(data) {
    unique(data.frame(
      sample_id = as.character(data$sample_id),
      condition = as.character(data$condition),
      stringsAsFactors = FALSE
    ))
  }
  source_blocks <- lapply(data_sources, make_block_table)
  common_patient_ids <- Reduce(
    intersect,
    lapply(source_blocks, function(x) unique(x$sample_id[!is.na(x$sample_id)]))
  )
  source_blocks <- lapply(source_blocks, function(x) {
    x[x$sample_id %in% common_patient_ids & x$condition %in% c("0", "1"), , drop = FALSE]
  })

  if (comparison_mode == "condition_vs_condition") {
    common_block_keys <- Reduce(
      intersect,
      lapply(source_blocks, function(x) paste(x$sample_id, x$condition, sep = "___"))
    )
    common_blocks <- source_blocks$shap[
      paste(source_blocks$shap$sample_id, source_blocks$shap$condition, sep = "___") %in% common_block_keys,
      ,
      drop = FALSE
    ]
    patient_ids_by_condition <- lapply(c("0", "1"), function(level) {
      unique(common_blocks$sample_id[common_blocks$condition == level])
    })
    patient_counts <- lengths(patient_ids_by_condition)
    support_available <- all(patient_counts >= 2L)
    unavailable_reason <- if (support_available) NA_character_ else "insufficient patients"
  } else {
    patient_ids_by_condition <- NULL
    patient_counts <- vapply(c("0", "1"), function(level) {
      length(unique(source_blocks$shap$sample_id[source_blocks$shap$condition == level]))
    }, integer(1))
    support_available <- length(common_patient_ids) >= 2L &&
      all(vapply(data_sources, function(x) all(c("0", "1") %in% as.character(x$condition)), logical(1)))
    unavailable_reason <- if (support_available) NA_character_ else "insufficient patients"
  }

  unavailable_result <- function(reason) {
    support_table <- data.frame(
      feature = feature_names,
      cross_patient_support = NA_real_,
      cross_patient_support_available = FALSE,
      cross_patient_support_unavailable_reason = reason,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    diagnostics <- data.frame(
      feature = feature_names,
      cross_patient_support = NA_real_,
      shap_importance_support = NA_real_,
      expression_signal_support = NA_real_,
      value_shap_relationship_support = NA_real_,
      n_bootstraps_requested = n_bootstraps,
      n_bootstraps_valid_task = 0L,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    list(
      support = support_table,
      diagnostics = diagnostics,
      available = FALSE,
      reason = reason,
      n_patients_condition0 = patient_counts[1],
      n_patients_condition1 = patient_counts[2],
      n_patients_total = length(common_patient_ids),
      n_bootstraps_requested = n_bootstraps,
      n_bootstraps_valid_task = 0L
    )
  }

  if (!support_available) {
    warning(
      "Cross-patient support unavailable: insufficient patients. ",
      "Final scores use four renormalized components.",
      call. = FALSE
    )
    return(unavailable_result(unavailable_reason))
  }

  scalers <- list(
    importance = fit_shapture_score_scaler(
      full_importance,
      method = score_scaling_method,
      probs = score_scaling_probs
    ),
    log2fc = fit_shapture_score_scaler(
      abs(full_stats$expression_log2FC),
      method = score_scaling_method,
      probs = score_scaling_probs
    ),
    detection = fit_shapture_score_scaler(
      abs(full_stats$det_diff),
      method = score_scaling_method,
      probs = score_scaling_probs
    ),
    consistency = fit_shapture_score_scaler(
      full_stats$within_condition_consistency_raw,
      method = score_scaling_method,
      probs = score_scaling_probs
    )
  )

  n_features <- length(feature_names)
  importance_bootstrap <- matrix(NA_real_, nrow = n_bootstraps, ncol = n_features)
  expression_bootstrap <- matrix(NA_real_, nrow = n_bootstraps, ncol = n_features)
  value_shap_bootstrap <- matrix(NA_real_, nrow = n_bootstraps, ncol = n_features)
  colnames(importance_bootstrap) <- feature_names
  colnames(expression_bootstrap) <- feature_names
  colnames(value_shap_bootstrap) <- feature_names

  build_source_indices <- function(source_data, sampled_blocks) {
    unlist(lapply(seq_len(nrow(sampled_blocks)), function(i) {
      which(
        as.character(source_data$sample_id) == sampled_blocks$sample_id[i] &
          if (comparison_mode == "condition_vs_condition") {
            as.character(source_data$condition) == sampled_blocks$condition[i]
          } else {
            TRUE
          }
      )
    }), use.names = FALSE)
  }

  draw_valid_bootstrap <- function(max_attempts = 100L) {
    for (attempt in seq_len(max_attempts)) {
      if (comparison_mode == "condition_vs_condition") {
        sampled_blocks <- do.call(rbind, lapply(seq_along(patient_ids_by_condition), function(i) {
          level <- as.character(i - 1L)
          ids <- patient_ids_by_condition[[i]]
          data.frame(
            sample_id = sample(ids, size = length(ids), replace = TRUE),
            condition = level,
            stringsAsFactors = FALSE
          )
        }))
      } else {
        sampled_blocks <- data.frame(
          sample_id = sample(
            common_patient_ids,
            size = length(common_patient_ids),
            replace = TRUE
          ),
          condition = NA_character_,
          stringsAsFactors = FALSE
        )
      }

      source_indices <- lapply(data_sources, build_source_indices, sampled_blocks = sampled_blocks)
      source_valid <- vapply(seq_along(data_sources), function(i) {
        idx <- source_indices[[i]]
        length(idx) > 0 && all(c("0", "1") %in% as.character(data_sources[[i]]$condition[idx]))
      }, logical(1))
      if (all(source_valid)) {
        return(source_indices)
      }
    }
    NULL
  }

  set.seed(seed)
  task_bootstrap_valid <- logical(n_bootstraps)
  for (bootstrap_index in seq_len(n_bootstraps)) {
    source_indices <- draw_valid_bootstrap()
    if (is.null(source_indices)) {
      next
    }
    shap_index <- source_indices$shap
    model_index <- source_indices$model
    single_cell_index <- source_indices$single_cell

    bootstrap_importance <- colMeans(abs(shap_matrix[shap_index, , drop = FALSE]), na.rm = TRUE)
    importance_bootstrap[bootstrap_index, ] <- apply_shapture_score_scaler(
      bootstrap_importance,
      scalers$importance
    )

    value_shap_bootstrap[bootstrap_index, ] <- vapply(seq_len(n_features), function(j) {
      x <- shap_feature_values[shap_index, j]
      s <- shap_matrix[shap_index, j]
      valid <- is.finite(x) & is.finite(s)
      if (sum(valid) < 3L || length(unique(x[valid])) < 2L || length(unique(s[valid])) < 2L) {
        return(NA_real_)
      }
      suppressWarnings(stats::cor(x[valid], s[valid], method = "spearman"))
    }, numeric(1))

    model_matrix <- as.matrix(model_expression_data[model_index, feature_names, drop = FALSE])
    model_condition <- as.character(model_expression_data$condition[model_index])
    mean_0 <- colMeans(model_matrix[model_condition == "0", , drop = FALSE], na.rm = TRUE)
    mean_1 <- colMeans(model_matrix[model_condition == "1", , drop = FALSE], na.rm = TRUE)
    pseudocount <- compute_dynamic_pseudocount(model_matrix, fraction = 0.01, fallback = 1e-6)
    bootstrap_log2fc <- log2((mean_1 + pseudocount) / (mean_0 + pseudocount))
    log2fc_retained <- ifelse(
      is.finite(bootstrap_log2fc) &
        abs(bootstrap_log2fc) >= feature_min_log2FC &
        sign(bootstrap_log2fc) == sign(full_stats$expression_log2FC),
      abs(bootstrap_log2fc),
      0
    )
    log2fc_score <- apply_shapture_score_scaler(log2fc_retained, scalers$log2fc)
    # A direction reversal or sub-threshold effect must remain unsupported even
    # when the full-data reference component is constant and scales to neutral.
    log2fc_score[log2fc_retained == 0] <- 0

    single_cell_matrix <- as.matrix(
      single_cell_expression_data[single_cell_index, feature_names, drop = FALSE]
    )
    single_cell_condition <- as.character(single_cell_expression_data$condition[single_cell_index])
    detection_0 <- colMeans(single_cell_matrix[single_cell_condition == "0", , drop = FALSE] > 0)
    detection_1 <- colMeans(single_cell_matrix[single_cell_condition == "1", , drop = FALSE] > 0)
    bootstrap_detection_diff <- detection_1 - detection_0
    full_detection_sign <- sign(full_stats$det_diff)
    detection_retained <- ifelse(
      full_detection_sign != 0 &
        is.finite(bootstrap_detection_diff) &
        sign(bootstrap_detection_diff) == full_detection_sign,
      abs(bootstrap_detection_diff),
      0
    )
    detection_score <- apply_shapture_score_scaler(detection_retained, scalers$detection)
    detection_score[detection_retained == 0] <- 0

    consistency_result <- tryCatch(
      compute_within_condition_expression_consistency(
        expression_data = single_cell_matrix,
        condition = single_cell_condition,
        condition_levels = c("0", "1"),
        upper_quantile = 0.99,
        min_positive_cells = 20L,
        pseudocount_fraction = 0.01,
        pseudocount_fallback = 1e-6
      ),
      error = function(e) NULL
    )
    if (is.null(consistency_result)) {
      importance_bootstrap[bootstrap_index, ] <- NA_real_
      value_shap_bootstrap[bootstrap_index, ] <- NA_real_
      next
    }
    consistency_match <- match(feature_names, consistency_result$feature)
    consistency_score <- apply_shapture_score_scaler(
      consistency_result$within_condition_consistency_raw[consistency_match],
      scalers$consistency
    )
    bootstrap_expression <- (log2fc_score + detection_score + consistency_score) / 3
    bootstrap_expression[consistency_result$both_conditions_all_zero[consistency_match]] <- 0
    expression_bootstrap[bootstrap_index, ] <- bootstrap_expression
    task_bootstrap_valid[bootstrap_index] <- TRUE
  }

  if (!any(task_bootstrap_valid)) {
    warning(
      "Cross-patient support unavailable: no valid patient-block bootstraps. ",
      "Final scores use four renormalized components.",
      call. = FALSE
    )
    return(unavailable_result("no valid patient-block bootstraps"))
  }

  importance_retention <- t(vapply(seq_len(n_features), function(j) {
    compute_smoother_importance_retention(
      importance_bootstrap[, j],
      lower_probability = 0.05
    )
  }, numeric(7)))
  importance_all_value_summary <- t(vapply(seq_len(n_features), function(j) {
    compute_bootstrap_retention(importance_bootstrap[, j])
  }, numeric(4)))
  expression_retention <- t(vapply(seq_len(n_features), function(j) {
    compute_bootstrap_retention(expression_bootstrap[, j])
  }, numeric(4)))
  colnames(importance_retention) <- c(
    "support", "nonzero_frequency", "positive_median", "positive_q05",
    "positive_retention", "n_valid", "n_positive"
  )
  colnames(importance_all_value_summary) <- c("support", "median", "q05", "n_valid")
  colnames(expression_retention) <- c("support", "median", "q05", "n_valid")

  value_shap_support <- vapply(seq_len(n_features), function(j) {
    bootstrap_correlations <- value_shap_bootstrap[, j]
    valid <- is.finite(bootstrap_correlations)
    if (!any(valid)) {
      return(0)
    }
    mean(
      sign(bootstrap_correlations[valid]) == sign(full_stats$shap_expr_cor[j]) &
        abs(bootstrap_correlations[valid]) >= shap_expr_cor_cutoff
    )
  }, numeric(1))
  value_shap_valid <- colSums(is.finite(value_shap_bootstrap))

  cross_patient_support <- rowMeans(cbind(
    importance_retention[, "support"],
    expression_retention[, "support"],
    value_shap_support
  ))
  support_table <- data.frame(
    feature = feature_names,
    cross_patient_support = cross_patient_support,
    cross_patient_support_available = TRUE,
    cross_patient_support_unavailable_reason = NA_character_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  diagnostics <- data.frame(
    feature = feature_names,
    cross_patient_support = cross_patient_support,
    shap_importance_support = importance_retention[, "support"],
    expression_signal_support = expression_retention[, "support"],
    value_shap_relationship_support = value_shap_support,
    shap_importance_bootstrap_median = importance_all_value_summary[, "median"],
    shap_importance_bootstrap_q05 = importance_all_value_summary[, "q05"],
    shap_importance_nonzero_frequency = importance_retention[, "nonzero_frequency"],
    shap_importance_positive_median = importance_retention[, "positive_median"],
    shap_importance_positive_q05 = importance_retention[, "positive_q05"],
    shap_importance_positive_retention = importance_retention[, "positive_retention"],
    expression_signal_bootstrap_median = expression_retention[, "median"],
    expression_signal_bootstrap_q05 = expression_retention[, "q05"],
    shap_importance_valid_bootstraps = as.integer(importance_retention[, "n_valid"]),
    shap_importance_positive_bootstraps = as.integer(importance_retention[, "n_positive"]),
    expression_signal_valid_bootstraps = as.integer(expression_retention[, "n_valid"]),
    value_shap_valid_bootstraps = as.integer(value_shap_valid),
    n_bootstraps_requested = n_bootstraps,
    n_bootstraps_valid_task = sum(task_bootstrap_valid),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  list(
    support = support_table,
    diagnostics = diagnostics,
    available = TRUE,
    reason = NA_character_,
    n_patients_condition0 = patient_counts[1],
    n_patients_condition1 = patient_counts[2],
    n_patients_total = length(common_patient_ids),
    n_bootstraps_requested = n_bootstraps,
    n_bootstraps_valid_task = sum(task_bootstrap_valid)
  )
}


#' Compute legacy patient-level SHAP summary statistics
#'
#' Collapses observation-level SHAP values to patient-level units and computes
#' descriptive support plus permutation/bootstrap summaries.
#'
#' @param shap_matrix Numeric SHAP matrix with rows aligned to `metadata`.
#' @param metadata Data frame aligned to the SHAP matrix rows.
#' @param sample_col Sample identifier column.
#' @param condition_col Condition column.
#' @param n_permutations Requested number of patient-label permutations.
#' @param n_bootstraps Requested number of bootstrap replicates.
#' @param seed Random seed.
#'
#' @return Data frame with one row per feature.
#' @keywords internal
#' @noRd
compute_patient_level_shap_stats <- function(
    shap_matrix,
    metadata,
    sample_col = "sample_id",
    condition_col = "condition",
    n_permutations = 1000,
    n_bootstraps = 1000,
    seed = 42
) {
  if (!is.matrix(shap_matrix)) {
    shap_matrix <- as.matrix(shap_matrix)
  }
  
  if (nrow(shap_matrix) != nrow(metadata)) {
    stop("shap_matrix and metadata must have the same number of rows.")
  }
  
  if (!sample_col %in% colnames(metadata)) {
    stop("sample_col is missing from metadata.")
  }
  
  if (!condition_col %in% colnames(metadata)) {
    stop("condition_col is missing from metadata.")
  }
  
  patient_df <- data.frame(
    sample_id = metadata[[sample_col]],
    condition = metadata[[condition_col]],
    stringsAsFactors = FALSE
  )
  
  patient_units <- unique(patient_df)
  patient_units <- patient_units[
    !is.na(patient_units$sample_id) & !is.na(patient_units$condition),
    ,
    drop = FALSE
  ]
  
  if (nrow(patient_units) < 4 || length(unique(patient_units$condition)) != 2) {
    observed <- colMeans(shap_matrix[metadata[[condition_col]] == 1, , drop = FALSE], na.rm = TRUE) -
      colMeans(shap_matrix[metadata[[condition_col]] == 0, , drop = FALSE], na.rm = TRUE)
    
    return(data.frame(
      feature = colnames(shap_matrix),
      effect_size = observed,
      pval = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      bootstrap_support = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  unit_keys <- paste(patient_df$sample_id, patient_df$condition, sep = "___")
  patient_unit_keys <- paste(patient_units$sample_id, patient_units$condition, sep = "___")
  
  patient_means <- t(vapply(patient_unit_keys, function(key) {
    idx <- which(unit_keys == key)
    colMeans(shap_matrix[idx, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(shap_matrix))))
  
  rownames(patient_means) <- patient_unit_keys
  colnames(patient_means) <- colnames(shap_matrix)
  
  unit_conditions <- patient_units$condition
  .compute_patient_level_shap_stats_from_patient_means(
    patient_means = patient_means,
    unit_conditions = unit_conditions,
    n_permutations = n_permutations,
    n_bootstraps = n_bootstraps,
    seed = seed
  )
}


#' Compute patient-level SHAP statistics from patient means
#'
#' @param patient_means Numeric matrix of patient-by-feature means.
#' @param unit_conditions Condition vector aligned to `patient_means`.
#' @param n_permutations Requested number of patient-label permutations.
#' @param n_bootstraps Requested number of bootstrap replicates.
#' @param seed Random seed.
#'
#' @return Data frame with one row per feature.
#' @keywords internal
#' @noRd
.compute_patient_level_shap_stats_from_patient_means <- function(
    patient_means,
    unit_conditions,
    n_permutations = 1000,
    n_bootstraps = 1000,
    seed = 42
) {
  if (!is.matrix(patient_means)) {
    patient_means <- as.matrix(patient_means)
  }
  
  if (nrow(patient_means) != length(unit_conditions)) {
    stop("patient_means and unit_conditions must align.")
  }
  
  cond_levels <- sort(unique(unit_conditions))
  idx_0 <- which(unit_conditions == cond_levels[1])
  idx_1 <- which(unit_conditions == cond_levels[2])
  
  observed <- colMeans(patient_means[idx_1, , drop = FALSE], na.rm = TRUE) -
    colMeans(patient_means[idx_0, , drop = FALSE], na.rm = TRUE)
  
  n_distinct_perms <- choose(nrow(patient_means), length(idx_1))
  if (!is.finite(n_distinct_perms) || n_distinct_perms < 1) {
    n_distinct_perms <- 1
  }
  
  n_permutations <- as.integer(min(n_permutations, max(n_distinct_perms * 5L, 99L)))
  n_bootstraps <- as.integer(min(n_bootstraps, 500L))
  
  set.seed(seed)
  perm_null <- matrix(NA_real_, nrow = n_permutations, ncol = ncol(patient_means))
  n_cond1 <- length(idx_1)
  all_idx <- seq_len(nrow(patient_means))
  
  for (b in seq_len(n_permutations)) {
    perm_idx_1 <- sample(all_idx, size = n_cond1, replace = FALSE)
    perm_idx_0 <- setdiff(all_idx, perm_idx_1)
    perm_null[b, ] <- colMeans(patient_means[perm_idx_1, , drop = FALSE], na.rm = TRUE) -
      colMeans(patient_means[perm_idx_0, , drop = FALSE], na.rm = TRUE)
  }
  
  pvals <- (colSums(
    abs(perm_null) >= matrix(abs(observed), nrow = n_permutations, ncol = ncol(patient_means), byrow = TRUE),
    na.rm = TRUE
  ) + 1) / (n_permutations + 1)
  
  set.seed(seed + 1L)
  boot_dist <- matrix(NA_real_, nrow = n_bootstraps, ncol = ncol(patient_means))
  
  for (b in seq_len(n_bootstraps)) {
    boot_idx_0 <- sample(idx_0, size = length(idx_0), replace = TRUE)
    boot_idx_1 <- sample(idx_1, size = length(idx_1), replace = TRUE)
    boot_dist[b, ] <- colMeans(patient_means[boot_idx_1, , drop = FALSE], na.rm = TRUE) -
      colMeans(patient_means[boot_idx_0, , drop = FALSE], na.rm = TRUE)
  }
  
  ci_lower <- apply(boot_dist, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
  ci_upper <- apply(boot_dist, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  bootstrap_support <- ifelse(
    observed >= 0,
    colMeans(boot_dist > 0, na.rm = TRUE),
    colMeans(boot_dist < 0, na.rm = TRUE)
  )
  
  data.frame(
    feature = colnames(patient_means),
    effect_size = observed,
    pval = pvals,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    bootstrap_support = bootstrap_support,
    stringsAsFactors = FALSE
  )
}


#' Compute the legacy minimum achievable permutation p-value
#'
#' @param n_condition0 Number of patients in condition 0.
#' @param n_condition1 Number of patients in condition 1.
#'
#' @return Minimum non-zero p-value under exact patient-label permutations.
#' @keywords internal
#' @noRd
compute_min_achievable_permutation_p <- function(n_condition0, n_condition1) {
  if (!is.finite(n_condition0) || !is.finite(n_condition1) ||
      n_condition0 < 1 || n_condition1 < 1) {
    return(NA_real_)
  }
  
  n_distinct <- choose(n_condition0 + n_condition1, n_condition1)
  if (!is.finite(n_distinct) || n_distinct < 1) {
    return(NA_real_)
  }
  
  1 / (n_distinct + 1)
}


#' Select the legacy patient-level support component
#'
#' Retained internally so older saved workflows can still be interpreted. The
#' current ranking uses `compute_cross_patient_support()` instead.
#'
#' @param stats_results Data frame containing support and adjusted p-values.
#' @param min_patients Minimum number of patients across the two comparison arms.
#' @param bootstrap_support_col Name of the support column.
#' @param pval_col Name of the adjusted p-value column.
#'
#' @return A list describing the score, label, and interpretation metadata.
#' @keywords internal
#' @noRd
select_patient_level_support_component <- function(
    stats_results,
    min_patients,
    bootstrap_support_col = "bootstrap_support",
    pval_col = "padj"
) {
  bootstrap_support <- stats_results[[bootstrap_support_col]]
  bootstrap_support[!is.finite(bootstrap_support)] <- 0.5
  
  if (is.na(min_patients) || min_patients < 3) {
    return(list(
      score = bootstrap_support,
      label = "bootstrap_support",
      reliable = FALSE,
      plot_label = "Patient-level support",
      note = paste0(
        "min_patients=", min_patients,
        "; descriptive-only mode; bootstrap_support used as the patient-level support component."
      )
    ))
  }
  
  pval <- stats_results[[pval_col]]
  pval[!is.finite(pval)] <- 1
  score <- -log10(pmax(pval, 1e-10))
  
  list(
    score = score,
    label = "permutation_padj",
    reliable = TRUE,
    plot_label = "Patient-level support",
    note = paste0(
      "min_patients=", min_patients,
      "; permutation-adjusted p-values used as the patient-level support component."
    )
  )
}

# Backward-compatible non-exported alias for older scripts.
select_significance_component <- function(...) {
  select_patient_level_support_component(...)
}
