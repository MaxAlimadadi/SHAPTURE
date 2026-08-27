#' Normalize SHAPTURE score prior weights
#'
#' @param prior_weights Named numeric vector of score prior weights.
#'
#' @return Named numeric vector with internal component names.
#' @export
normalize_score_weights <- function(prior_weights = NULL) {
  default_weights <- c(
    importance = 0.25,
    cross_patient_support = 0.10,
    model_reproducibility = 0.10,
    expression = 0.25,
    value_shap = 0.30
  )
  
  if (is.null(prior_weights)) {
    return(default_weights)
  }
  if (!is.numeric(prior_weights) || is.null(names(prior_weights))) {
    stop("default_prior_scores_weights must be a named numeric vector.")
  }
  
  name_map <- c(
    SHAP_importance_prior_weight = "importance",
    importance = "importance",
    cross_patient_support_prior_weight = "cross_patient_support",
    cross_patient_support = "cross_patient_support",
    patient_level_support_prior_weight = "cross_patient_support",
    patient_level_support = "cross_patient_support",
    model_reproducibility_prior_weight = "model_reproducibility",
    model_reproducibility = "model_reproducibility",
    stability = "model_reproducibility",
    expression_signal_prior_weight = "expression",
    expression = "expression",
    value_SHAP_relationship_prior_weight = "value_shap",
    value_shap = "value_shap"
  )
  
  mapped_names <- unname(name_map[names(prior_weights)])
  if (any(is.na(mapped_names))) {
    stop(
      "Unrecognized default_prior_scores_weights name(s): ",
      paste(names(prior_weights)[is.na(mapped_names)], collapse = ", ")
    )
  }
  
  out <- default_weights
  for (i in seq_along(prior_weights)) {
    out[mapped_names[i]] <- prior_weights[i]
  }
  
  if (any(!is.finite(out)) || any(out < 0) || sum(out) <= 0) {
    stop("default_prior_scores_weights must contain finite non-negative values with a positive sum.")
  }
  
  out / sum(out)
}


#' Calculate a positive median anchor
#'
#' @param x Nonnegative numeric reference values.
#' @param zero_tolerance Values at or below this threshold are treated as zero.
#' @param component_name Component name used in validation messages.
#'
#' @return The median finite value above `zero_tolerance`, or `NA` when no
#'   positive reference value is available.
#' @keywords internal
compute_positive_median_anchor <- function(
    x,
    zero_tolerance = 1e-8,
    component_name = "component"
) {
  finite_values <- as.numeric(x[is.finite(x)])
  if (any(finite_values < -zero_tolerance)) {
    stop(component_name, " anchor values must be nonnegative.")
  }
  positive_values <- finite_values[finite_values > zero_tolerance]
  if (length(positive_values) == 0L) {
    return(NA_real_)
  }
  stats::median(positive_values)
}


#' Apply a median-anchor transformation
#'
#' @param x Nonnegative numeric values.
#' @param anchor Positive median anchor calculated from the fixed full-data
#'   candidate-feature set.
#' @param zero_tolerance Values at or below this threshold are treated as zero.
#' @param component_name Component name used in validation messages.
#'
#' @return Numeric vector on a 0-1 scale. Finite positive values are transformed
#'   as `x / (x + anchor)`; finite zeros remain zero and non-finite inputs remain
#'   `NA`.
#' @keywords internal
apply_median_anchor_transform <- function(
    x,
    anchor,
    zero_tolerance = 1e-8,
    component_name = "component"
) {
  out <- rep(NA_real_, length(x))
  finite_index <- is.finite(x)
  if (!any(finite_index)) {
    return(out)
  }

  finite_values <- as.numeric(x[finite_index])
  if (any(finite_values < -zero_tolerance)) {
    stop(component_name, " values must be nonnegative before transformation.")
  }
  finite_values[abs(finite_values) <= zero_tolerance] <- 0

  if (!is.numeric(anchor) || length(anchor) != 1L || !is.finite(anchor) || anchor <= 0) {
    if (any(finite_values > zero_tolerance)) {
      stop(component_name, " requires one finite positive median anchor.")
    }
    out[finite_index] <- 0
    return(out)
  }

  transformed <- finite_values / (finite_values + anchor)
  out[finite_index] <- pmin(1, pmax(0, transformed))
  out
}


#' Compute adaptive component weights
#'
#' @param combined Combined feature table containing scaled component columns.
#' @param min_weight Minimum allowed weight for any component.
#' @param use_mad Whether to use MAD instead of SD.
#' @param prior_weights Named prior weights.
#' @param alpha Shrinkage factor toward adaptive weights.
#' @param winsorize Whether to winsorize component values before dispersion.
#' @param winsor_probs Winsorization quantiles.
#'
#' @return Named numeric vector of weights.
#' @export
compute_adaptive_weights <- function(
    combined,
    min_weight = 0.03,
    use_mad = TRUE,
    prior_weights = NULL,
    alpha = 0.5,
    winsorize = TRUE,
    winsor_probs = c(0.05, 0.95)
) {
  support_available <- "cross_patient_support" %in% colnames(combined) &&
    any(is.finite(combined$cross_patient_support))
  component_cols <- c(
    importance = "importance_scaled",
    cross_patient_support = "cross_patient_support",
    model_reproducibility = "model_reproducibility_scaled",
    expression = "expression_score",
    value_shap = "shap_expr_cor_abs_scaled"
  )
  if (!support_available) {
    component_cols <- component_cols[names(component_cols) != "cross_patient_support"]
  }
  
  missing_cols <- setdiff(component_cols, colnames(combined))
  if (length(missing_cols) > 0) {
    stop(
      paste(
        "Missing columns for adaptive weighting:",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  prior_weights <- normalize_score_weights(prior_weights)
  prior_weights <- prior_weights[names(component_cols)]
  prior_weights <- prior_weights / sum(prior_weights)
  
  winsorize_vec <- function(v, probs = c(0.05, 0.95)) {
    v_finite <- v[is.finite(v)]
    if (length(v_finite) <= 2) {
      return(v)
    }
    qs <- stats::quantile(v_finite, probs = probs, na.rm = TRUE, names = FALSE)
    pmin(pmax(v, qs[1]), qs[2])
  }
  
  raw_scores <- sapply(component_cols, function(col_name) {
    v <- combined[[col_name]]
    v <- v[is.finite(v)]
    if (length(v) <= 1) {
      return(0)
    }
    if (winsorize) {
      v <- winsorize_vec(v, probs = winsor_probs)
    }
    if (use_mad) {
      stats::mad(v, center = stats::median(v), constant = 1, na.rm = TRUE)
    } else {
      stats::sd(v, na.rm = TRUE)
    }
  })
  
  if (sum(raw_scores, na.rm = TRUE) == 0) {
    raw_scores[] <- 1
  }
  
  w_adaptive <- raw_scores / sum(raw_scores, na.rm = TRUE)
  w <- alpha * w_adaptive + (1 - alpha) * prior_weights
  w <- pmax(w, min_weight)
  w <- w / sum(w)
  
  names(w) <- paste0("w_", names(component_cols))
  
  w
}


#' Build model-expression pattern labels
#'
#' @param what_0 Label for condition 0.
#' @param what_1 Label for condition 1.
#'
#' @return Character vector containing the three possible pattern labels.
#' @keywords internal
model_expression_pattern_labels <- function(what_0, what_1) {
  c(
    condition_1 = paste0("Higher feature values push toward ", what_1),
    condition_0 = paste0("Higher feature values push toward ", what_0),
    insufficient = "Insufficient directional evidence"
  )
}


#' Build display labels for model-expression pattern legends
#'
#' @param what_0 Label for condition 0.
#' @param what_1 Label for condition 1.
#'
#' @return Named character vector mapping stored labels to display labels.
#' @keywords internal
model_expression_pattern_display_labels <- function(what_0, what_1) {
  labels <- model_expression_pattern_labels(what_0 = what_0, what_1 = what_1)
  stats::setNames(
    c(
      paste0("Higher feature values push\ntoward ", what_1),
      paste0("Higher feature values push\ntoward ", what_0),
      "Insufficient directional\nevidence"
    ),
    unname(labels)
  )
}


#' Build model-expression pattern colors
#'
#' @param what_0 Label for condition 0.
#' @param what_1 Label for condition 1.
#' @param group_colors Two colors: condition 0 then condition 1.
#'
#' @return Named character vector of colors.
#' @keywords internal
model_expression_pattern_colors <- function(what_0, what_1, group_colors) {
  labels <- model_expression_pattern_labels(what_0 = what_0, what_1 = what_1)
  stats::setNames(
    c(unname(group_colors[2]), unname(group_colors[1]), "grey55"),
    unname(labels)
  )
}


#' Ordered model-expression pattern labels used in plots
#'
#' @param what_0 Label for condition 0.
#' @param what_1 Label for condition 1.
#'
#' @return Named character vector for condition 1 first, then condition 0.
#' @keywords internal
model_expression_pattern_plot_labels <- function(what_0, what_1) {
  labels <- model_expression_pattern_labels(what_0 = what_0, what_1 = what_1)
  display_labels <- model_expression_pattern_display_labels(what_0 = what_0, what_1 = what_1)
  display_labels[unname(labels[c("condition_1", "condition_0")])]
}


#' Compute model-expression pattern labels
#'
#' @param expression_log2FC Expression log2 fold-change values.
#' @param shap_expr_cor Feature-value to SHAP correlation values.
#' @param what_0 Label for condition 0.
#' @param what_1 Label for condition 1.
#' @param expression_log2FC_cutoff Minimum absolute expression log2FC needed to
#'   assign a model-expression pattern label.
#' @param shap_expr_cor_cutoff Minimum absolute value-SHAP correlation needed
#'   to assign a model-expression pattern label.
#'
#' @return Character vector of model-expression pattern labels.
#' @export
compute_model_expression_pattern <- function(
    expression_log2FC,
    shap_expr_cor,
    what_0 = "condition_0",
    what_1 = "condition_1",
    expression_log2FC_cutoff = 0.25,
    shap_expr_cor_cutoff = 0.25
) {
  labels <- model_expression_pattern_labels(what_0 = what_0, what_1 = what_1)
  dplyr::case_when(
    !is.finite(expression_log2FC) | !is.finite(shap_expr_cor) ~ unname(labels["insufficient"]),
    abs(expression_log2FC) < expression_log2FC_cutoff | abs(shap_expr_cor) < shap_expr_cor_cutoff ~ unname(labels["insufficient"]),
    sign(expression_log2FC) != sign(shap_expr_cor) ~ unname(labels["insufficient"]),
    expression_log2FC > 0 & shap_expr_cor > 0 ~ unname(labels["condition_1"]),
    expression_log2FC < 0 & shap_expr_cor < 0 ~ unname(labels["condition_0"]),
    TRUE ~ unname(labels["insufficient"])
  )
}


#' Add model-expression pattern annotations to a feature table
#'
#' @param df Data frame containing effect size, expression log2FC, and detection
#'   summaries.
#' @param what_0 Label for condition 0.
#' @param what_1 Label for condition 1.
#' @param expression_log2FC_cutoff Minimum absolute expression log2FC needed to
#'   assign a model-expression pattern label.
#' @param shap_expr_cor_cutoff Minimum absolute value-SHAP correlation needed
#'   to assign a model-expression pattern label.
#'
#' @return Annotated data frame.
#' @export
add_direction_info <- function(
    df,
    what_0,
    what_1,
    expression_log2FC_cutoff = 0.25,
    shap_expr_cor_cutoff = 0.25
) {
  if (!"expression_log2FC" %in% colnames(df)) {
    stop("df must contain expression_log2FC.")
  }
  if (!"shap_expr_cor" %in% colnames(df)) {
    df$shap_expr_cor <- NA_real_
  }
  
  df %>%
    dplyr::mutate(
      Model_expression_pattern = compute_model_expression_pattern(
        expression_log2FC = expression_log2FC,
        shap_expr_cor = shap_expr_cor,
        what_0 = what_0,
        what_1 = what_1,
        expression_log2FC_cutoff = expression_log2FC_cutoff,
        shap_expr_cor_cutoff = shap_expr_cor_cutoff
      )
    )
}


#' Compute the integrated SHAPTURE ranking
#'
#' @param average_importance Data frame with `feature` and SHAP importance.
#' @param stats_results Data frame containing SHAP statistics and expression
#'   summaries, including `within_condition_consistency_raw`.
#' @param w_importance Weight for SHAP importance.
#' @param w_cross_patient_support Weight for Cross-patient support.
#' @param w_model_reproducibility Weight for model reproducibility.
#' @param w_expression Weight for expression signal.
#' @param w_value_shap Weight for value-to-SHAP association.
#' @param adaptive_weights Whether to compute adaptive weights.
#' @param min_adaptive_weight Minimum adaptive weight.
#' @param default_prior_scores_weights Named prior weights for the five score
#'   components.
#' @param what_0 Label for condition 0.
#' @param what_1 Label for condition 1.
#' @param expression_log2FC_cutoff Minimum absolute expression log2FC needed to
#'   keep a model-expression pattern.
#' @param shap_expr_cor_cutoff Minimum absolute value-SHAP correlation needed to
#'   keep a model-expression pattern.
#' @param shap_importance_median_anchor Optional positive median anchor fitted
#'   from the fixed pre-pattern candidate-feature SHAP importances.
#' @param expression_log2FC_abs_median_anchor Optional positive median anchor
#'   fitted from the fixed pre-pattern candidate-feature absolute log2FC values.
#' @param score_scaling_method Backward-compatible configuration field. Main
#'   evidence scores use median-anchor transformation for unbounded inputs and
#'   retain naturally bounded inputs without across-feature scaling.
#' @param score_scaling_probs Quantiles used only to winsorize component values
#'   when estimating optional adaptive-weight spread.
#'
#' @return Ranked data frame.
#' @export
compute_combined_shap_ranking <- function(
    average_importance,
    stats_results,
    w_importance = 0.25,
    w_cross_patient_support = 0.10,
    w_model_reproducibility = 0.10,
    w_expression = 0.25,
    w_value_shap = 0.30,
    adaptive_weights = FALSE,
    min_adaptive_weight = 0.03,
    default_prior_scores_weights = NULL,
    what_0 = "condition_0",
    what_1 = "condition_1",
    expression_log2FC_cutoff = 0.25,
    shap_expr_cor_cutoff = 0.25,
    shap_importance_median_anchor = NULL,
    expression_log2FC_abs_median_anchor = NULL,
    score_scaling_method = c("median_anchor", "winsorized", "minmax"),
    score_scaling_probs = c(0.05, 0.95)
) {
  score_scaling_method <- match.arg(score_scaling_method)
  if (!is.numeric(score_scaling_probs) ||
      length(score_scaling_probs) != 2 ||
      any(!is.finite(score_scaling_probs)) ||
      any(score_scaling_probs < 0 | score_scaling_probs > 1) ||
      score_scaling_probs[1] >= score_scaling_probs[2]) {
    stop("score_scaling_probs must be two increasing numeric probabilities between 0 and 1.")
  }
  
  combined <- merge(average_importance, stats_results, by = "feature")
  
  if (!"expression_log2FC" %in% colnames(combined)) {
    stop("stats_results must contain expression_log2FC.")
  }

  if (is.null(shap_importance_median_anchor)) {
    shap_importance_median_anchor <- compute_positive_median_anchor(
      combined$`average|SHAP value|`,
      component_name = "SHAP importance"
    )
  }
  if (is.null(expression_log2FC_abs_median_anchor)) {
    expression_log2FC_abs_median_anchor <- compute_positive_median_anchor(
      abs(combined$expression_log2FC),
      component_name = "absolute expression log2FC"
    )
  }
  
  combined$Model_expression_pattern <- compute_model_expression_pattern(
    expression_log2FC = combined$expression_log2FC,
    shap_expr_cor = combined$shap_expr_cor,
    what_0 = what_0,
    what_1 = what_1,
    expression_log2FC_cutoff = expression_log2FC_cutoff,
    shap_expr_cor_cutoff = shap_expr_cor_cutoff
  )
  
  insufficient_label <- unname(model_expression_pattern_labels(what_0 = what_0, what_1 = what_1)["insufficient"])
  n_features_before_pattern_filter <- nrow(combined)
  combined <- combined[combined$Model_expression_pattern != insufficient_label, , drop = FALSE]
  n_features_after_pattern_filter <- nrow(combined)
  if (n_features_after_pattern_filter == 0) {
    stop("No features remain after removing features with insufficient model-expression directional evidence.")
  }
  
  combined$importance_scaled <- apply_median_anchor_transform(
    combined$`average|SHAP value|`,
    anchor = shap_importance_median_anchor,
    component_name = "SHAP importance"
  )
  
  if (!"cross_patient_support" %in% colnames(combined)) {
    stop("stats_results must contain cross_patient_support.")
  }
  if (!"cross_patient_support_available" %in% colnames(combined)) {
    stop("stats_results must contain cross_patient_support_available.")
  }
  support_available_values <- unique(combined$cross_patient_support_available)
  if (length(support_available_values) != 1 || is.na(support_available_values)) {
    stop("cross_patient_support_available must be one constant logical value within a ranking task.")
  }
  cross_patient_support_available <- isTRUE(support_available_values[1])
  if (cross_patient_support_available) {
    if (any(!is.finite(combined$cross_patient_support)) ||
        any(combined$cross_patient_support < 0 | combined$cross_patient_support > 1)) {
      stop("Available cross_patient_support values must be finite and between 0 and 1.")
    }
  } else {
    combined$cross_patient_support[] <- NA_real_
  }
  
  model_reproducibility <- if ("model_reproducibility" %in% colnames(combined)) {
    combined$model_reproducibility
  } else {
    combined$stability
  }
  if (any(!is.finite(model_reproducibility)) ||
      any(model_reproducibility < 0 | model_reproducibility > 1)) {
    stop("Model reproducibility must contain finite values from 0 to 1.")
  }
  combined$model_reproducibility <- model_reproducibility
  combined$model_reproducibility_scaled <- model_reproducibility
  combined$effect_abs <- abs(combined$effect_size)
  combined$expression_log2FC_abs <- abs(combined$expression_log2FC)
  combined$expression_log2FC_abs_scaled <- apply_median_anchor_transform(
    combined$expression_log2FC_abs,
    anchor = expression_log2FC_abs_median_anchor,
    component_name = "absolute expression log2FC"
  )
  det_diff_abs <- abs(combined$det_diff)
  if (any(!is.finite(det_diff_abs)) || any(det_diff_abs < 0 | det_diff_abs > 1)) {
    stop("Absolute detection-rate differences must contain finite values from 0 to 1.")
  }
  combined$det_diff_scaled <- det_diff_abs
  if (!"within_condition_consistency_raw" %in% colnames(combined)) {
    stop("stats_results must contain within_condition_consistency_raw.")
  }
  if (any(!is.finite(combined$within_condition_consistency_raw)) ||
      any(combined$within_condition_consistency_raw < 0 | combined$within_condition_consistency_raw > 1)) {
    stop("within_condition_consistency_raw must contain finite values from 0 to 1.")
  }
  combined$within_condition_consistency_scaled <- combined$within_condition_consistency_raw
  combined$expression_score <- (
    combined$expression_log2FC_abs_scaled +
      combined$det_diff_scaled +
      combined$within_condition_consistency_scaled
  ) / 3
  
  shap_expr_cor_abs <- abs(combined$shap_expr_cor)
  if (any(!is.finite(shap_expr_cor_abs)) || any(shap_expr_cor_abs < 0 | shap_expr_cor_abs > 1)) {
    stop("Absolute Value-SHAP correlations must contain finite values from 0 to 1.")
  }
  combined$shap_expr_cor_abs_scaled <- shap_expr_cor_abs
  combined$shap_importance_median_anchor <- shap_importance_median_anchor
  combined$expression_log2FC_abs_median_anchor <- expression_log2FC_abs_median_anchor
  combined$score_scaling_strategy <- "median_anchor_unbounded_natural_bounded"
  combined$n_features_before_model_expression_pattern_filter <- n_features_before_pattern_filter
  combined$n_features_after_model_expression_pattern_filter <- n_features_after_pattern_filter
  
  prior_weights <- normalize_score_weights(
    if (is.null(default_prior_scores_weights)) {
      c(
        importance = w_importance,
        cross_patient_support = w_cross_patient_support,
        model_reproducibility = w_model_reproducibility,
        expression = w_expression,
        value_shap = w_value_shap
      )
    } else {
      default_prior_scores_weights
    }
  )
  
  if (!cross_patient_support_available) {
    prior_weights["cross_patient_support"] <- 0
    prior_weights <- prior_weights / sum(prior_weights)
  }

  if (adaptive_weights) {
    adaptive_weight_values <- compute_adaptive_weights(
      combined = combined,
      min_weight = min_adaptive_weight,
      use_mad = TRUE,
      prior_weights = prior_weights,
      alpha = 0.5,
      winsorize = TRUE,
      winsor_probs = score_scaling_probs
    )
    w_vec <- c(
      w_importance = 0,
      w_cross_patient_support = 0,
      w_model_reproducibility = 0,
      w_expression = 0,
      w_value_shap = 0
    )
    w_vec[names(adaptive_weight_values)] <- adaptive_weight_values
  } else {
    w_vec <- prior_weights
    names(w_vec) <- c(
      "w_importance",
      "w_cross_patient_support",
      "w_model_reproducibility",
      "w_expression",
      "w_value_shap"
    )
  }
  
  combined$weight_importance <- w_vec["w_importance"]
  combined$weight_cross_patient_support <- w_vec["w_cross_patient_support"]
  combined$weight_model_reproducibility <- w_vec["w_model_reproducibility"]
  combined$weight_expression <- w_vec["w_expression"]
  combined$weight_value_shap <- w_vec["w_value_shap"]
  
  cross_patient_contribution <- rep(0, nrow(combined))
  if (cross_patient_support_available) {
    cross_patient_contribution <-
      w_vec["w_cross_patient_support"] * combined$cross_patient_support
  }

  combined$final_score <-
    w_vec["w_importance"] * combined$importance_scaled +
    cross_patient_contribution +
    w_vec["w_model_reproducibility"] * combined$model_reproducibility_scaled +
    w_vec["w_expression"] * combined$expression_score +
    w_vec["w_value_shap"] * combined$shap_expr_cor_abs_scaled
  
  combined$final_score[!is.finite(combined$final_score)] <- 0
  combined <- combined[order(combined$final_score, decreasing = TRUE), ]
  combined$rank <- seq_len(nrow(combined))
  
  combined
}
