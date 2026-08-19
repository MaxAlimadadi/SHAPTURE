#' Plot score decomposition
#'
#' @param df Ranked feature table.
#' @param top_n_to_show Number of features to show.
#'
#' @return A ggplot object.
#' @export
plot_score_decomposition <- function(df, top_n_to_show) {
  required_cols <- c(
    "feature", "final_score", "importance_scaled",
    "patient_level_support_scaled", "model_reproducibility_scaled",
    "expression_score", "shap_expr_cor_abs_scaled", "weight_importance",
    "weight_patient_level_support", "weight_model_reproducibility",
    "weight_expression", "weight_value_shap"
  )
  
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }
  
  df_plot <- df %>%
    dplyr::arrange(dplyr::desc(final_score)) %>%
    head(top_n_to_show) %>%
    dplyr::mutate(
      contribution_importance = weight_importance * importance_scaled,
      contribution_patient_level_support = weight_patient_level_support * patient_level_support_scaled,
      contribution_model_reproducibility = weight_model_reproducibility * model_reproducibility_scaled,
      contribution_expression = weight_expression * expression_score,
      contribution_value_shap = weight_value_shap * shap_expr_cor_abs_scaled
    ) %>%
    dplyr::select(
      feature, final_score, contribution_importance,
      contribution_patient_level_support, contribution_model_reproducibility,
      contribution_expression, contribution_value_shap
    )
  
  df_long <- df_plot %>%
    tidyr::pivot_longer(
      cols = tidyselect::starts_with("contribution_"),
      names_to = "component",
      values_to = "value"
    )
  
  df_long$component <- dplyr::recode(
    df_long$component,
    contribution_importance = "SHAP importance",
    contribution_patient_level_support = "Patient-level support",
    contribution_model_reproducibility = "Model reproducibility",
    contribution_expression = "Expression signal",
    contribution_value_shap = "Value-SHAP relationship"
  )
  component_order <- c(
    "SHAP importance",
    "Expression signal",
    "Value-SHAP relationship",
    "Patient-level support",
    "Model reproducibility"
  )
  df_long$component <- factor(df_long$component, levels = component_order)
  
  feature_order <- df_plot$feature
  df_long$feature <- factor(df_long$feature, levels = rev(feature_order))
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = feature, y = value, fill = component)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "",
      y = "Contribution to final score",
      fill = "Component"
    ) +
    ggplot2::scale_fill_manual(values = c(
      "SHAP importance" = "#FEAF16",
      "Expression signal" = "#004CFF",
      "Value-SHAP relationship" = "#1b9e77",
      "Patient-level support" = "#d95f02",
      "Model reproducibility" = "#7570b3"
    ), breaks = component_order)
}


#' Plot the five final-score components as a heatmap
#'
#' @param df Ranked feature table.
#' @param top_n_to_show Number of features to show. Use `nrow(df)` for all
#'   retained features.
#' @param continuous_colors Continuous color palette.
#' @param title Plot title.
#' @param tile_color Tile border color. Use `NA` to remove borders.
#' @param tile_linewidth Tile border line width.
#'
#' @return A ggplot object.
#' @export
plot_five_component_heatmap <- function(
    df,
    top_n_to_show,
    continuous_colors,
    title = "Five-component support matrix",
    tile_color = "white",
    tile_linewidth = 0.25
) {
  if (is.null(tile_color)) {
    tile_color <- NA_character_
  }
  if (length(tile_color) != 1) {
    stop("tile_color must be one color value or NA.")
  }
  if (!is.na(tile_color)) {
    tryCatch(
      grDevices::col2rgb(tile_color),
      error = function(e) stop("tile_color must be a valid R color or NA.")
    )
  }
  if (!is.numeric(tile_linewidth) ||
      length(tile_linewidth) != 1 ||
      !is.finite(tile_linewidth) ||
      tile_linewidth < 0) {
    stop("tile_linewidth must be one finite numeric value >= 0.")
  }
  
  component_cols <- c(
    "SHAP importance" = "importance_scaled",
    "Expression signal" = "expression_score",
    "Value-SHAP relationship" = "shap_expr_cor_abs_scaled",
    "Patient-level support" = "patient_level_support_scaled",
    "Model reproducibility" = "model_reproducibility_scaled"
  )
  required_cols <- c("feature", unname(component_cols))
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }
  
  df_plot <- df %>%
    dplyr::arrange(dplyr::desc(final_score)) %>%
    head(top_n_to_show) %>%
    dplyr::select(feature, dplyr::all_of(unname(component_cols)))
  colnames(df_plot) <- c("feature", names(component_cols))
  df_plot$feature <- factor(df_plot$feature, levels = rev(df_plot$feature))
  
  df_long <- reshape2::melt(
    df_plot,
    id.vars = "feature",
    variable.name = "component",
    value.name = "score"
  )
  df_long$component <- factor(df_long$component, levels = names(component_cols))
  y_text_size <- if (nrow(df_plot) > 75) 4.5 else 7
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = component, y = feature, fill = score)) +
    ggplot2::geom_tile(color = tile_color, linewidth = tile_linewidth) +
    ggplot2::scale_fill_gradientn(
      colors = continuous_colors,
      limits = c(0, 1),
      oob = scales::squish,
      name = "Score"
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 70, hjust = 1, vjust = 1, size = 9, colour = "black"),
      axis.text.y = ggplot2::element_text(size = y_text_size, colour = "black"),
      axis.title = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.5),
      plot.title = ggplot2::element_text(size = 12, face = "bold", lineheight = 1.05)
    ) +
    ggplot2::labs(title = title)
}


#' Plot the top-feature support matrix
#'
#' @param df Ranked feature table.
#' @param top_n_to_show Number of features to show.
#' @param continuous_colors Continuous color palette.
#'
#' @return A ggplot object.
#' @export
plot_evidence_heatmap <- function(df, top_n_to_show, continuous_colors) {
  required_cols <- c("feature", "final_score", "stability", "patient_level_support_scaled")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }
  
  df_plot <- df %>%
    dplyr::arrange(dplyr::desc(final_score)) %>%
    head(top_n_to_show) %>%
    dplyr::transmute(
      feature = factor(feature, levels = rev(feature)),
      `Final score` = final_score,
      `Model reproducibility` = stability,
      `Patient-level support` = patient_level_support_scaled
    )
  
  df_long <- reshape2::melt(
    df_plot,
    id.vars = "feature",
    variable.name = "metric",
    value.name = "score"
  )
  df_long$metric <- factor(
    df_long$metric,
    levels = c("Final score", "Model reproducibility", "Patient-level support")
  )
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = metric, y = feature, fill = score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.45) +
    ggplot2::scale_fill_gradientn(
      colors = continuous_colors,
      limits = c(0, 1),
      oob = scales::squish,
      name = "Score"
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 70, hjust = 1, vjust = 1, size = 9, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 7, colour = "black"),
      axis.title = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.5),
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 7.5, lineheight = 1.05)
    ) +
    ggplot2::labs(
      title = "Support matrix",
      subtitle = "Direct comparisons of final scores,\nmodel reproducibility, and patient-level support"
    )
}


#' Plot absolute SHAP effect versus biological expression effect
#'
#' @param df Ranked feature table.
#' @param top_n_to_show Number of features to show.
#' @param top_n_to_name Number of features to label.
#' @param what_0 Condition 0 label.
#' @param what_1 Condition 1 label.
#' @param group_colors Two colors: condition 0 then condition 1.
#' @param expression_log2FC_cutoff Minimum absolute expression log2FC needed
#'   to assign a model-expression pattern label.
#' @param shap_expr_cor_cutoff Minimum absolute value-SHAP correlation needed
#'   to assign a model-expression pattern label.
#' @param dot_label_size Size of feature labels drawn next to selected dots.
#'
#' @return A ggplot object.
#' @export
plot_shap_vs_expression <- function(
    df,
    top_n_to_show,
    top_n_to_name,
    what_0,
    what_1,
    group_colors = c("darkblue", "darkred"),
    expression_log2FC_cutoff = 0.25,
    shap_expr_cor_cutoff = 0.25,
    dot_label_size = 3
) {
  df <- add_direction_info(
    df,
    what_0 = what_0,
    what_1 = what_1,
    expression_log2FC_cutoff = expression_log2FC_cutoff,
    shap_expr_cor_cutoff = shap_expr_cor_cutoff
  ) %>%
    dplyr::mutate(effect_abs = abs(effect_size)) %>%
    dplyr::arrange(dplyr::desc(final_score)) %>%
    head(top_n_to_show)
  pattern_colors <- model_expression_pattern_colors(
    what_0 = what_0,
    what_1 = what_1,
    group_colors = group_colors
  )
  pattern_labels <- model_expression_pattern_plot_labels(what_0 = what_0, what_1 = what_1)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = expression_log2FC, y = .data$effect_abs)) +
    ggplot2::geom_point(
      ggplot2::aes(
        color = Model_expression_pattern,
        size = `average|SHAP value|`
      ),
      alpha = 0.8
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Expression log2FC",
      y = "Absolute SHAP effect size",
      color = "Model\nexpression\npattern",
      size = "SHAP importance"
    ) +
    ggplot2::scale_color_manual(
      values = pattern_colors,
      breaks = names(pattern_labels),
      labels = unname(pattern_labels),
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.32, 0.08)))
  
  top_df <- df %>% dplyr::arrange(dplyr::desc(final_score)) %>% head(top_n_to_name)
  p + ggrepel::geom_text_repel(
    data = top_df,
    ggplot2::aes(label = feature),
    size = dot_label_size,
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.35,
    force = 1.5
  )
}


#' Plot per-feature value-to-SHAP associations
#'
#' @param df Ranked feature table.
#' @param top_n_to_show Number of features to show.
#' @param top_n_to_name Number of features to label.
#' @param what_0 Condition 0 label.
#' @param what_1 Condition 1 label.
#' @param continuous_colors Continuous color palette.
#' @param group_colors Two colors: condition 0 then condition 1.
#' @param expression_log2FC_cutoff Minimum absolute expression log2FC needed
#'   to assign a model-expression pattern label.
#' @param shap_expr_cor_cutoff Minimum absolute value-SHAP correlation needed
#'   to assign a model-expression pattern label.
#' @param dot_label_size Size of feature labels drawn next to selected dots.
#'
#' @return A ggplot object.
#' @export
plot_valueSHAP_by_feature <- function(
    df,
    top_n_to_show,
    top_n_to_name,
    what_0,
    what_1,
    continuous_colors,
    group_colors = c("darkblue", "darkred"),
    expression_log2FC_cutoff = 0.25,
    shap_expr_cor_cutoff = 0.25,
    dot_label_size = 3
) {
  df <- df %>%
    add_direction_info(
      what_0 = what_0,
      what_1 = what_1,
      expression_log2FC_cutoff = expression_log2FC_cutoff,
      shap_expr_cor_cutoff = shap_expr_cor_cutoff
    ) %>%
    dplyr::arrange(dplyr::desc(final_score)) %>%
    head(top_n_to_show) %>%
    dplyr::mutate(feature = factor(feature, levels = rev(feature)))
  pattern_colors <- model_expression_pattern_colors(
    what_0 = what_0,
    what_1 = what_1,
    group_colors = group_colors
  )
  pattern_labels <- model_expression_pattern_plot_labels(what_0 = what_0, what_1 = what_1)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = shap_expr_cor, y = feature)) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = `average|SHAP value|`,
        fill = patient_level_support_scaled,
        shape = Model_expression_pattern
      ),
      color = "black",
      alpha = 0.9,
      stroke = 0.8
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::scale_fill_gradientn(
      colors = continuous_colors,
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      oob = scales::squish
    ) +
    ggplot2::scale_shape_manual(
      values = stats::setNames(c(21, 24), names(pattern_labels)),
      breaks = names(pattern_labels),
      labels = unname(pattern_labels),
      drop = FALSE
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = paste0("Value -> SHAP correlation\n(+ = high expr drives ", what_1, ", - = high expr drives ", what_0, ")"),
      y = "",
      fill = "Patient-level\nsupport",
      shape = "Model\nexpression\npattern",
      size = "SHAP importance"
    )
  
  top_df <- df %>% dplyr::arrange(dplyr::desc(final_score)) %>% head(top_n_to_name)
  p + ggrepel::geom_text_repel(
    data = top_df,
    ggplot2::aes(label = feature),
    size = dot_label_size,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3
  )
}
