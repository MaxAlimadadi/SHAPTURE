#' Normalize SHAPTURE score prior weights
#'
#' @param prior_weights Named numeric vector of score prior weights.
#'
#' @return Named numeric vector with internal component names.
#' @export
normalize_score_weights <- function(prior_weights = NULL) {
  default_weights <- c(
    importance = 0.25,
    patient_level_support = 0.10,
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
    patient_level_support_prior_weight = "patient_level_support",
    patient_level_support = "patient_level_support",
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
  component_cols <- c(
    importance = "importance_scaled",
    patient_level_support = "patient_level_support_scaled",
    model_reproducibility = "model_reproducibility_scaled",
    expression = "expression_score",
    value_shap = "shap_expr_cor_abs_scaled"
  )
  
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
  
  names(w) <- c(
    "w_importance",
    "w_patient_level_support",
    "w_model_reproducibility",
    "w_expression",
    "w_value_shap"
  )
  
  w
}


#' Add a descriptive-only caption to a plot when needed
#'
#' @param p A ggplot or plot-grid object.
#' @param inference_reliable Logical flag indicating whether inference is
#'   considered reliable.
#' @param min_patients Minimum patients per comparison arm.
#'
#' @return Plot object with an optional caption.
#' @export
add_inference_warning <- function(p, inference_reliable, min_patients) {
  if (!isTRUE(inference_reliable)) {
    p <- p +
      ggplot2::labs(
        caption = paste0("DESCRIPTIVE ONLY - min patients/group = ", min_patients)
      ) +
      ggplot2::theme(
        plot.caption = ggplot2::element_text(
          colour = "firebrick",
          size = 7,
          face = "italic"
        )
      )
  }
  p
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
#' @param w_patient_level_support Weight for patient-level support.
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
#' @param score_scaling_method Component scaling method. `"winsorized"` caps
#'   finite values at `score_scaling_probs` before min-max scaling; `"minmax"`
#'   uses ordinary min-max scaling.
#' @param score_scaling_probs Two quantiles used when
#'   `score_scaling_method = "winsorized"`.
#'
#' @return Ranked data frame.
#' @export
compute_combined_shap_ranking <- function(
    average_importance,
    stats_results,
    w_importance = 0.25,
    w_patient_level_support = 0.10,
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
    score_scaling_method = c("winsorized", "minmax"),
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
  
  safe_rescale <- function(x, neutral = 0.5) {
    finite_idx <- is.finite(x)
    out <- rep(neutral, length(x))
    
    if (!any(finite_idx)) {
      return(out)
    }
    
    finite_vals <- x[finite_idx]
    if (length(unique(finite_vals)) <= 1) {
      out[finite_idx] <- neutral
      return(out)
    }
    
    if (identical(score_scaling_method, "winsorized") && length(finite_vals) > 2) {
      qs <- stats::quantile(finite_vals, probs = score_scaling_probs, na.rm = TRUE, names = FALSE)
      finite_vals <- pmin(pmax(finite_vals, qs[1]), qs[2])
      if (length(unique(finite_vals)) <= 1) {
        out[finite_idx] <- neutral
        return(out)
      }
    }
    
    out[finite_idx] <- scales::rescale(finite_vals)
    out
  }
  
  combined <- merge(average_importance, stats_results, by = "feature")
  
  if (!"expression_log2FC" %in% colnames(combined)) {
    stop("stats_results must contain expression_log2FC.")
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
  
  combined$importance_scaled <- safe_rescale(combined$`average|SHAP value|`)
  
  min_patients_values <- unique(combined$min_patients)
  if (length(min_patients_values) != 1) {
    stop("stats_results$min_patients must be constant within one ranking task.")
  }
  
  support_component <- select_patient_level_support_component(
    stats_results = combined,
    min_patients = min_patients_values[1],
    bootstrap_support_col = "bootstrap_support",
    pval_col = "padj"
  )
  
  combined$bootstrap_support[!is.finite(combined$bootstrap_support)] <- 0.5
  combined$bootstrap_support_scaled <- safe_rescale(combined$bootstrap_support)
  combined$patient_level_support_score <- support_component$score
  combined$patient_level_support_scaled <- safe_rescale(support_component$score)
  combined$patient_level_support_method <- support_component$label
  combined$patient_level_support_plot_label <- support_component$plot_label
  combined$patient_level_support_note <- support_component$note
  combined$inference_reliable <- support_component$reliable
  
  combined$model_reproducibility_scaled <- safe_rescale(combined$stability)
  combined$effect_abs <- abs(combined$effect_size)
  combined$expression_log2FC_abs <- abs(combined$expression_log2FC)
  combined$expression_log2FC_abs_scaled <- safe_rescale(combined$expression_log2FC_abs)
  combined$det_diff_scaled <- safe_rescale(abs(combined$det_diff))
  if (!"within_condition_consistency_raw" %in% colnames(combined)) {
    stop("stats_results must contain within_condition_consistency_raw.")
  }
  if (any(!is.finite(combined$within_condition_consistency_raw)) ||
      any(combined$within_condition_consistency_raw < 0 | combined$within_condition_consistency_raw > 1)) {
    stop("within_condition_consistency_raw must contain finite values from 0 to 1.")
  }
  combined$within_condition_consistency_scaled <- safe_rescale(
    combined$within_condition_consistency_raw
  )
  combined$expression_score <- (
    combined$expression_log2FC_abs_scaled +
      combined$det_diff_scaled +
      combined$within_condition_consistency_scaled
  ) / 3
  
  combined$shap_expr_cor_abs_scaled <- safe_rescale(abs(combined$shap_expr_cor))
  combined$n_features_before_model_expression_pattern_filter <- n_features_before_pattern_filter
  combined$n_features_after_model_expression_pattern_filter <- n_features_after_pattern_filter
  
  prior_weights <- normalize_score_weights(
    if (is.null(default_prior_scores_weights)) {
      c(
        importance = w_importance,
        patient_level_support = w_patient_level_support,
        model_reproducibility = w_model_reproducibility,
        expression = w_expression,
        value_shap = w_value_shap
      )
    } else {
      default_prior_scores_weights
    }
  )
  
  if (adaptive_weights) {
    w_vec <- compute_adaptive_weights(
      combined = combined,
      min_weight = min_adaptive_weight,
      use_mad = TRUE,
      prior_weights = prior_weights,
      alpha = 0.5,
      winsorize = TRUE,
      winsor_probs = score_scaling_probs
    )
  } else {
    w_vec <- prior_weights
    names(w_vec) <- c(
      "w_importance",
      "w_patient_level_support",
      "w_model_reproducibility",
      "w_expression",
      "w_value_shap"
    )
  }
  
  combined$weight_importance <- w_vec["w_importance"]
  combined$weight_patient_level_support <- w_vec["w_patient_level_support"]
  combined$weight_model_reproducibility <- w_vec["w_model_reproducibility"]
  combined$weight_expression <- w_vec["w_expression"]
  combined$weight_value_shap <- w_vec["w_value_shap"]
  
  combined$final_score <-
    w_vec["w_importance"] * combined$importance_scaled +
    w_vec["w_patient_level_support"] * combined$patient_level_support_scaled +
    w_vec["w_model_reproducibility"] * combined$model_reproducibility_scaled +
    w_vec["w_expression"] * combined$expression_score +
    w_vec["w_value_shap"] * combined$shap_expr_cor_abs_scaled
  
  combined$final_score[!is.finite(combined$final_score)] <- 0
  combined <- combined[order(combined$final_score, decreasing = TRUE), ]
  combined$rank <- seq_len(nrow(combined))
  
  combined
}
