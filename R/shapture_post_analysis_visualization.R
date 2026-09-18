#' Save reusable inputs for post-analysis SHAPTURE visualization
#'
#' @param output_dir Completed SHAPTURE task folder.
#' @param combined_results Ranked feature table.
#' @param ROC_all Held-out prediction table pooled across iterations.
#' @param AUC_all Iteration-level AUC table.
#' @param roc_summary Named list with pooled and iteration AUC summaries.
#' @param model_performance_all Iteration-level model performance table.
#' @param importance_data_all Iteration-level feature-importance table.
#' @param condition_aggregate_data Condition-level mean SHAP and feature values.
#' @param beeswarm_source_S Unit-averaged SHAP matrix for saved beeswarm features.
#' @param beeswarm_source_X Unit-averaged model-input feature-value matrix/data frame.
#' @param beeswarm_baseline Mean model baseline used for SHAP-only summaries.
#' @param density_source Data frame containing condition and saved density features.
#' @param context Named list of task and plotting context.
#' @param group_colors Two group colors: condition 0 then condition 1.
#' @param continuous_colors Continuous color palette.
#'
#' @return Invisibly returns the written RDS path.
#' @keywords internal
save_shapture_visualization_inputs <- function(
    output_dir,
    combined_results,
    ROC_all,
    AUC_all,
    roc_summary,
    model_performance_all,
    importance_data_all,
    condition_aggregate_data,
    beeswarm_source_S,
    beeswarm_source_X,
    beeswarm_baseline,
    density_source,
    context,
    group_colors,
    continuous_colors
) {
  docs_dir <- file.path(output_dir, "docs")
  dir.create(docs_dir, recursive = TRUE, showWarnings = FALSE)
  
  inputs <- list(
    version = "1.1.0",
    generated_at = as.character(Sys.time()),
    context = context,
    group_colors = unname(group_colors),
    continuous_colors = unname(continuous_colors),
    roc = list(
      ROC_all = ROC_all,
      AUC_all = AUC_all,
      summary = roc_summary
    ),
    model_performance_all = model_performance_all,
    importance_data_all = importance_data_all,
    shap = list(
      condition_aggregate_data = condition_aggregate_data,
      beeswarm_source_S = beeswarm_source_S,
      beeswarm_source_X = as.data.frame(beeswarm_source_X, check.names = FALSE),
      beeswarm_baseline = beeswarm_baseline
    ),
    density = list(
      source = density_source
    ),
    combined_results_snapshot = combined_results
  )
  
  output_path <- file.path(docs_dir, "visualization_inputs.rds")
  saveRDS(inputs, output_path)
  invisible(output_path)
}


#' Find completed SHAPTURE result folders
#'
#' Searches the immediate subfolders of an analysis directory and retains only
#' folders containing the saved table and visualization-input bundle required
#' by the post-analysis visualization pipeline. Folder names and comparison
#' modes do not affect discovery.
#'
#' @param analysis_root Parent directory containing completed SHAPTURE result
#'   folders.
#'
#' @return A sorted character vector of valid result-folder paths.
#' @export
find_shapture_result_folders <- function(analysis_root) {
  if (!is.character(analysis_root) || length(analysis_root) != 1L ||
      is.na(analysis_root) || !nzchar(analysis_root) || !dir.exists(analysis_root)) {
    stop("analysis_root must be an existing directory.")
  }

  candidate_folders <- list.dirs(
    path = analysis_root,
    full.names = TRUE,
    recursive = FALSE
  )
  required_files <- c(
    file.path("docs", "combined_results.csv"),
    file.path("docs", "visualization_inputs.rds")
  )
  valid_folders <- candidate_folders[vapply(
    candidate_folders,
    function(folder) all(file.exists(file.path(folder, required_files))),
    logical(1)
  )]

  if (length(valid_folders) == 0L) {
    stop(
      "No completed SHAPTURE result folders were found directly inside analysis_root."
    )
  }

  sort(valid_folders)
}


#' Calculate an automatic SHAPTURE figure size
#'
#' Uses the actual number of displayed features to increase figure dimensions
#' within plot-specific minimum and maximum bounds. `width_per_feature` and
#' `height_per_feature` are expressed in inches per feature. The returned values
#' are converted to centimeters for the beeswarm PNG, whose production script
#' uses `units = "cm"`.
#'
#' @param figure_type Figure layout type. See Details for accepted values.
#' @param n_features Number of features displayed in the figure.
#' @param width_per_feature Global horizontal growth in inches per feature.
#' @param height_per_feature Global vertical growth in inches per feature.
#'
#' @return A list containing `width`, `height`, and `units` for the graphics
#'   device.
#'
#' @details
#' Accepted figure types are `"roc"`, `"feature_importance"`, `"waterfall"`,
#' `"force"`, `"beeswarm"`, `"density"`, `"score_decomposition"`,
#' `"component_heatmap"`, `"evidence_heatmap"`, `"shap_vs_expression"`,
#' `"value_shap"`, and `"dependency"`. ROC, dependency, and scatter layouts
#' remain fixed. Density plots use a
#' near-square facet grid with approximately two inches allocated to each
#' panel; their width includes additional room for the legend. Force-plot width
#' uses a larger internal growth factor suited to its layout.
#'
#' @export
calculate_shapture_figure_size <- function(
    figure_type,
    n_features = 0L,
    width_per_feature = 0.01,
    height_per_feature = 0.13
) {
  valid_types <- c(
    "roc", "feature_importance", "waterfall", "force", "beeswarm",
    "density", "score_decomposition", "component_heatmap",
    "evidence_heatmap", "shap_vs_expression", "value_shap", "dependency"
  )
  if (!is.character(figure_type) || length(figure_type) != 1L ||
      is.na(figure_type) || !figure_type %in% valid_types) {
    stop("figure_type must be one of: ", paste(valid_types, collapse = ", "), ".")
  }
  if (!is.numeric(n_features) || length(n_features) != 1L ||
      !is.finite(n_features) || n_features < 0 || n_features != floor(n_features)) {
    stop("n_features must be a single non-negative integer.")
  }
  for (value_name in c("width_per_feature", "height_per_feature")) {
    value <- get(value_name)
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value < 0) {
      stop(value_name, " must be a single non-negative number.")
    }
  }

  if (identical(figure_type, "density")) {
    density_feature_count <- max(1L, as.integer(n_features))
    density_ncol <- max(1L, ceiling(sqrt(density_feature_count)))
    density_nrow <- max(1L, ceiling(density_feature_count / density_ncol))

    return(list(
      width = round(1.5 + 2 * density_ncol, 3),
      height = round(1 + 2 * density_nrow, 3),
      units = "in"
    ))
  }

  spec <- switch(
    figure_type,
    roc = list(7.5, 4.0, 7.5, 7.5, 4.0, 4.0, 0, 0, "in"),
    feature_importance = list(3.0, 3.0, 3.1, 5.0, 4.0, 30.0, 1, 0.14 / 0.13, "in"),
    waterfall = list(7.0, 4.0, 7.0, 9.0, 5.0, 30.0, 0.8, 0.14 / 0.13, "in"),
    force = list(8.0, 4.0, 8.0, 15.0, 4.0, 4.0, 3.5, 0, "in"),
    beeswarm = list(15.0, 7.0, 15.0, 20.0, 9.0, 60.0, 0.02 / (0.01 * 2.54), 0.36 / (0.13 * 2.54), "cm"),
    score_decomposition = list(5.5, 3.5, 5.5, 8.0, 4.5, 30.0, 1, 1, "in"),
    component_heatmap = list(5.2, 3.5, 5.2, 8.0, 4.5, 30.0, 1, 1, "in"),
    evidence_heatmap = list(4.0, 3.5, 4.0, 6.0, 4.5, 30.0, 0.8, 1, "in"),
    shap_vs_expression = list(8.0, 7.0, 8.0, 8.0, 7.0, 7.0, 0, 0, "in"),
    value_shap = list(8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 0, 0, "in"),
    dependency = list(8.6, 6.2, 8.6, 8.6, 6.2, 6.2, 0, 0, "in")
  )
  names(spec) <- c(
    "base_width", "base_height", "min_width", "max_width",
    "min_height", "max_height", "width_multiplier",
    "height_multiplier", "units"
  )
  unit_conversion <- if (identical(spec$units, "cm")) 2.54 else 1

  width <- spec$base_width +
    width_per_feature * n_features * spec$width_multiplier * unit_conversion
  height <- spec$base_height +
    height_per_feature * n_features * spec$height_multiplier * unit_conversion
  width <- max(spec$min_width, min(spec$max_width, width))
  height <- max(spec$min_height, min(spec$max_height, height))

  list(
    width = round(width, 3),
    height = round(height, 3),
    units = spec$units
  )
}


#' Build the pooled ROC plus iteration-AUC panel
#'
#' @param ROC_all Held-out prediction table pooled across iterations.
#' @param AUC_all Iteration-level AUC table.
#' @param mean_auc Mean iteration AUC.
#'
#' @return A ggplot/cowplot object.
#' @export
build_shapture_roc_plot <- function(ROC_all, AUC_all, mean_auc) {
  roc_plot_all_iteration <- pROC::roc(
    response = ROC_all$condition,
    predictor = ROC_all$pred,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  roc_plot_all_iteration_ggroc <- pROC::ggroc(
    roc_plot_all_iteration,
    legacy.axes = TRUE,
    colour = "steelblue",
    linewidth = 1
  ) +
    ggplot2::theme_classic() +
    ggplot2::annotate("text", x = 0.62, y = 0.16, label = paste0("Mean iteration AUC = ", round(mean_auc, 2)), size = 4) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "darkgray") +
    ggplot2::theme(panel.border = ggplot2::element_rect(linetype = "solid", fill = NA, linewidth = 1)) +
    ggplot2::xlab("False Positive Rate") +
    ggplot2::ylab("True Positive Rate") +
    ggplot2::ggtitle("Pooled ROC from all held-out predictions")
  
  set.seed(42)
  auc_boxplot <- ggplot2::ggplot(AUC_all, ggplot2::aes(x = "", y = auc)) +
    ggplot2::geom_boxplot(width = 0.22, fill = "grey85", color = "black", outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.08, size = 1.8, alpha = 0.75, color = "steelblue") +
    ggplot2::geom_hline(yintercept = mean_auc, linetype = "dashed", color = "firebrick") +
    ggplot2::annotate("text", x = 1, y = min(0.20, mean_auc), label = paste0("Mean = ", round(mean_auc, 2)), size = 3.8, color = "firebrick") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(linetype = "solid", fill = NA, linewidth = 1),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::ylab("AUC") +
    ggplot2::ggtitle("AUC across iterations")
  
  cowplot::plot_grid(
    roc_plot_all_iteration_ggroc,
    auc_boxplot,
    ncol = 2,
    rel_widths = c(1.6, 0.9),
    align = "h"
  )
}


#' Build the top-feature importance plot
#'
#' Feature names use the same model-expression pattern colors as the beeswarm
#' when both condition labels are supplied. Scores, feature order, bars, and
#' points are unchanged. Calls without condition labels retain plain axis text.
#'
#' @param combined_results Ranked feature table.
#' @param top_feature_show Number of top features to show.
#' @param task_what_0 Condition 0 label. Supply both condition labels to color
#'   feature names using the saved `Model_expression_pattern` column.
#' @param task_what_1 Condition 1 label.
#' @param group_colors Two group colors: condition 0 then condition 1.
#'
#' @return A ggplot object.
#' @export
build_shapture_feature_importance_plot <- function(
    combined_results,
    top_feature_show,
    task_what_0 = NULL,
    task_what_1 = NULL,
    group_colors = c("darkblue", "darkred")
) {
  top_feature_show <- min(top_feature_show, nrow(combined_results))
  if (!is.finite(top_feature_show) || top_feature_show < 1) {
    stop("top_feature_show must be at least 1.")
  }
  combined_results_top <- combined_results[seq_len(top_feature_show), , drop = FALSE]
  combined_results_top$feature <- factor(combined_results_top$feature, levels = rev(combined_results_top$feature))
  
  plot <- ggplot2::ggplot(combined_results_top, ggplot2::aes(x = final_score, y = feature)) +
    ggplot2::geom_col(width = 0.15) +
    ggplot2::geom_point(size = 1, color = "blue") +
    ggplot2::xlab("integrated importance score") +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(paste0("Top ", top_feature_show, " features")) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.text.x = ggplot2::element_text(colour = "black", size = 10),
      axis.text.y = ggplot2::element_text(colour = "black", size = 8),
      plot.title = ggplot2::element_text(size = 6)
    )

  if (is.null(task_what_0) && is.null(task_what_1)) {
    return(plot)
  }
  valid_condition_label <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }
  if (!valid_condition_label(task_what_0) ||
      !valid_condition_label(task_what_1) || identical(task_what_0, task_what_1)) {
    stop("Supply two distinct condition labels in task_what_0 and task_what_1.")
  }
  if (!"Model_expression_pattern" %in% colnames(combined_results_top)) {
    stop("combined_results must contain Model_expression_pattern to color feature names.")
  }
  group_colors <- validate_color_vector(group_colors, "group_colors", exact_length = 2)
  pattern_colors <- model_expression_pattern_colors(task_what_0, task_what_1, group_colors)
  pattern_labels <- model_expression_pattern_plot_labels(task_what_0, task_what_1)
  feature_names <- as.character(combined_results_top$feature)
  feature_colors <- unname(pattern_colors[as.character(combined_results_top$Model_expression_pattern)])
  feature_colors[is.na(feature_colors)] <- "grey55"
  feature_labels <- stats::setNames(
    paste0("<span style='color:", feature_colors, "'>", escape_shapture_html(feature_names), "</span>"),
    feature_names
  )
  legend_data <- data.frame(
    feature = factor(rep(feature_names[1], length(pattern_labels)), levels = levels(combined_results_top$feature)),
    final_score = combined_results_top$final_score[1],
    Model_expression_pattern = names(pattern_labels),
    stringsAsFactors = FALSE
  )

  plot +
    ggplot2::scale_y_discrete(labels = feature_labels) +
    ggplot2::geom_point(
      data = legend_data,
      ggplot2::aes(x = final_score, y = feature, fill = Model_expression_pattern),
      inherit.aes = FALSE,
      shape = 22,
      alpha = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = pattern_colors,
      breaks = names(pattern_labels),
      labels = unname(pattern_labels),
      drop = FALSE
    ) +
    ggplot2::labs(fill = "Model\nexpression\npattern") +
    ggplot2::guides(fill = ggplot2::guide_legend(
      ncol = 1,
      override.aes = list(alpha = 1, size = 3)
    )) +
    ggplot2::theme(
      axis.text.y = ggtext::element_markdown(size = 8),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 7),
      legend.key.size = grid::unit(0.3, "cm")
    )
}


#' Escape feature labels for optional markdown-colored axis text
#'
#' @param x Character vector.
#'
#' @return Escaped character vector.
#' @keywords internal
escape_shapture_html <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}


#' Build the SHAPTURE beeswarm plot from saved unit-level inputs
#'
#' @param combined_results Ranked feature table.
#' @param beeswarm_source_S Unit-averaged SHAP matrix.
#' @param beeswarm_source_X Unit-averaged feature-value table.
#' @param top_feature_show Number of top features.
#' @param task_what_0 Condition 0 label.
#' @param task_what_1 Condition 1 label.
#' @param group_colors Two group colors.
#' @param continuous_colors Continuous color palette.
#'
#' @return A ggplot object.
#' @export
build_shapture_beeswarm_plot <- function(
    combined_results,
    beeswarm_source_S,
    beeswarm_source_X,
    top_feature_show,
    task_what_0,
    task_what_1,
    group_colors,
    continuous_colors
) {
  top_feature_show <- min(top_feature_show, nrow(combined_results))
  if (!is.finite(top_feature_show) || top_feature_show < 1) {
    stop("top_feature_show must be at least 1.")
  }
  top_features_final <- head(combined_results$feature, top_feature_show)
  top_features_chr <- as.character(top_features_final)
  top_features_chr <- top_features_chr[
    top_features_chr %in% colnames(beeswarm_source_S) &
      top_features_chr %in% colnames(beeswarm_source_X)
  ]
  if (length(top_features_chr) == 0) {
    stop("No saved beeswarm features are available for the requested top features.")
  }
  if (length(top_features_chr) < length(top_features_final)) {
    warning(
      "Only ", length(top_features_chr), " of ", length(top_features_final),
      " requested top beeswarm features were saved in visualization_inputs.rds."
    )
  }
  
  beeswarm_pattern_colors <- model_expression_pattern_colors(
    what_0 = task_what_0,
    what_1 = task_what_1,
    group_colors = group_colors
  )
  beeswarm_pattern_labels <- model_expression_pattern_plot_labels(
    what_0 = task_what_0,
    what_1 = task_what_1
  )
  beeswarm_feature_patterns <- combined_results$Model_expression_pattern[
    match(top_features_chr, combined_results$feature)
  ]
  beeswarm_feature_colors <- unname(beeswarm_pattern_colors[beeswarm_feature_patterns])
  beeswarm_feature_colors[is.na(beeswarm_feature_colors)] <- "grey55"
  beeswarm_feature_labels <- stats::setNames(
    paste0(
      "<span style='color:",
      beeswarm_feature_colors,
      "'>",
      escape_shapture_html(top_features_chr),
      "</span>"
    ),
    top_features_chr
  )
  beeswarm_axis_text_y <- if (requireNamespace("ggtext", quietly = TRUE)) {
    ggtext::element_markdown(size = 10)
  } else {
    ggplot2::element_text(colour = "black", size = 10)
  }
  
  shap_sub <- beeswarm_source_S[, top_features_chr, drop = FALSE]
  shap_long <- reshape2::melt(shap_sub)
  colnames(shap_long) <- c("cell", "feature", "shap_value")
  shap_long$feature <- factor(as.character(shap_long$feature), levels = rev(top_features_chr))
  feature_values <- reshape2::melt(as.matrix(beeswarm_source_X[, top_features_chr, drop = FALSE]))
  shap_long$feature_value <- feature_values$value
  
  beeswarm_legend_df <- data.frame(
    feature = factor(rep(top_features_chr[1], length(beeswarm_pattern_labels)), levels = rev(top_features_chr)),
    shap_value = 0,
    Model_expression_pattern = names(beeswarm_pattern_labels),
    stringsAsFactors = FALSE
  )
  
  ggplot2::ggplot(shap_long, ggplot2::aes(x = shap_value, y = feature)) +
    ggbeeswarm::geom_quasirandom(
      ggplot2::aes(color = feature_value),
      method = "tukey",
      width = 0.3,
      size = 0.4,
      alpha = 0.7,
      orientation = "y",
      na.rm = TRUE
    ) +
    ggplot2::scale_color_gradientn(colors = continuous_colors) +
    ggplot2::geom_point(
      data = beeswarm_legend_df,
      ggplot2::aes(x = shap_value, y = feature, fill = Model_expression_pattern),
      inherit.aes = FALSE,
      shape = 22,
      size = 3,
      alpha = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = beeswarm_pattern_colors,
      breaks = names(beeswarm_pattern_labels),
      labels = unname(beeswarm_pattern_labels),
      drop = FALSE
    ) +
    ggplot2::scale_y_discrete(labels = beeswarm_feature_labels) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y = beeswarm_axis_text_y,
      axis.text.x = ggplot2::element_text(colour = "black"),
      plot.title = ggplot2::element_text(size = 10)
    ) +
    ggplot2::xlab("SHAP value") +
    ggplot2::labs(
      color = paste0("Scaled", "\n", "feature value"),
      fill = "Model\nexpression\npattern"
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 4))) +
    ggplot2::ylab("Features ranked by integrated scores") +
    ggplot2::ggtitle(paste0("Left = ", task_what_0, " | Right = ", task_what_1)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black")
}


#' Build shared equal-count groups for binned SHAP visualizations
#'
#' Creates deterministic, approximately equal-count groups jointly from SHAP
#' and model-seen feature values. Both dimensions are mapped to 0-1 only while
#' the groups are formed; the returned means remain on their original plotted
#' scales. The same returned groups can be passed to the binned beeswarm and
#' dependency builders. This function changes visualization only.
#'
#' @param beeswarm_source_S Unit-averaged SHAP matrix.
#' @param beeswarm_source_X Unit-averaged model-input feature-value matrix or
#'   data frame.
#' @param features Ranked feature names to group, in display order.
#' @param target_groups Desired maximum number of groups per feature.
#' @param min_units_per_group Minimum unit support for grouped displays.
#' @param raw_unit_threshold Below this number of valid units, retain one group
#'   per original unit instead of averaging.
#'
#' @return A data frame with one row per feature group and its unit count, mean
#'   SHAP value, and mean model-seen feature value.
#' @export
build_shapture_equal_count_groups <- function(
    beeswarm_source_S,
    beeswarm_source_X,
    features,
    target_groups = 500L,
    min_units_per_group = 10L,
    raw_unit_threshold = 500L
) {
  validate_positive_integer <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value)) {
      stop(name, " must be a single positive integer.")
    }
    as.integer(value)
  }
  target_groups <- validate_positive_integer(target_groups, "target_groups")
  min_units_per_group <- validate_positive_integer(
    min_units_per_group,
    "min_units_per_group"
  )
  raw_unit_threshold <- validate_positive_integer(
    raw_unit_threshold,
    "raw_unit_threshold"
  )

  beeswarm_source_S <- as.matrix(beeswarm_source_S)
  beeswarm_source_X <- as.data.frame(beeswarm_source_X, check.names = FALSE)
  if (nrow(beeswarm_source_S) != nrow(beeswarm_source_X)) {
    stop("beeswarm_source_S and beeswarm_source_X must have the same number of rows.")
  }
  if (is.null(colnames(beeswarm_source_S)) || is.null(colnames(beeswarm_source_X))) {
    stop("Both beeswarm input objects must have feature column names.")
  }
  if (!is.character(features) || length(features) < 1L || anyNA(features)) {
    stop("features must be a non-empty character vector without missing values.")
  }
  features <- unique(features[nzchar(features)])
  available_features <- features[
    features %in% colnames(beeswarm_source_S) &
      features %in% colnames(beeswarm_source_X)
  ]
  unavailable_features <- setdiff(features, available_features)
  if (length(unavailable_features) > 0L) {
    warning(
      "Omitting ", length(unavailable_features),
      " requested feature(s) missing from one or both saved beeswarm inputs: ",
      paste(unavailable_features, collapse = ", "), "."
    )
  }
  if (length(available_features) < 1L) {
    stop("None of the requested features are available in both saved beeswarm inputs.")
  }

  scale_for_grouping <- function(values) {
    value_range <- range(values, finite = TRUE)
    if (!all(is.finite(value_range)) || diff(value_range) == 0) {
      return(rep(0, length(values)))
    }
    (values - value_range[1]) / diff(value_range)
  }
  balanced_group_sizes <- function(n_units, n_groups) {
    sizes <- rep.int(n_units %/% n_groups, n_groups)
    remainder <- n_units %% n_groups
    if (remainder > 0L) {
      positions <- floor(
        ((seq_len(remainder) - 0.5) * n_groups) / remainder
      ) + 1L
      sizes[positions] <- sizes[positions] + 1L
    }
    if (sum(sizes) != n_units || diff(range(sizes)) > 1L) {
      stop("Internal equal-count group allocation failed.")
    }
    sizes
  }
  make_groups <- function(feature_values, shap_values, unit_ids, n_groups) {
    capacities <- balanced_group_sizes(length(feature_values), n_groups)
    feature_scaled <- scale_for_grouping(feature_values)
    shap_scaled <- scale_for_grouping(shap_values)
    assignment <- integer(length(feature_values))

    partition_group <- function(indices, group_ids) {
      if (length(group_ids) == 1L) {
        assignment[indices] <<- group_ids
        return(invisible(NULL))
      }

      n_left_groups <- floor(length(group_ids) / 2)
      left_group_ids <- group_ids[seq_len(n_left_groups)]
      right_group_ids <- group_ids[-seq_len(n_left_groups)]
      n_left_units <- sum(capacities[left_group_ids])
      feature_spread <- diff(stats::quantile(
        feature_scaled[indices],
        probs = c(0.05, 0.95),
        names = FALSE,
        na.rm = TRUE
      ))
      shap_spread <- diff(stats::quantile(
        shap_scaled[indices],
        probs = c(0.05, 0.95),
        names = FALSE,
        na.rm = TRUE
      ))
      if (!is.finite(feature_spread)) feature_spread <- 0
      if (!is.finite(shap_spread)) shap_spread <- 0

      if (feature_spread >= shap_spread) {
        ordered <- indices[order(
          feature_scaled[indices],
          shap_scaled[indices],
          unit_ids[indices],
          method = "radix"
        )]
      } else {
        ordered <- indices[order(
          shap_scaled[indices],
          feature_scaled[indices],
          unit_ids[indices],
          method = "radix"
        )]
      }

      partition_group(ordered[seq_len(n_left_units)], left_group_ids)
      partition_group(ordered[-seq_len(n_left_units)], right_group_ids)
      invisible(NULL)
    }

    partition_group(seq_along(feature_values), seq_len(n_groups))
    assignment
  }

  grouped_features <- lapply(available_features, function(feature) {
    feature_values <- suppressWarnings(as.numeric(beeswarm_source_X[[feature]]))
    shap_values <- suppressWarnings(as.numeric(beeswarm_source_S[, feature]))
    valid <- is.finite(feature_values) & is.finite(shap_values)
    n_valid <- sum(valid)
    if (n_valid < 1L) {
      warning("Feature ", feature, " has no paired finite SHAP and feature values; omitting it.")
      return(NULL)
    }

    feature_values <- feature_values[valid]
    shap_values <- shap_values[valid]
    unit_ids <- which(valid)
    n_groups <- if (n_valid < raw_unit_threshold) {
      n_valid
    } else {
      max(1L, min(target_groups, n_valid %/% min_units_per_group))
    }
    assignment <- make_groups(
      feature_values = feature_values,
      shap_values = shap_values,
      unit_ids = unit_ids,
      n_groups = n_groups
    )
    group_sizes <- tabulate(assignment, nbins = n_groups)
    mean_shap <- as.numeric(rowsum(shap_values, assignment, reorder = TRUE)) / group_sizes
    mean_feature_value <- as.numeric(
      rowsum(feature_values, assignment, reorder = TRUE)
    ) / group_sizes

    data.frame(
      feature = feature,
      group_id = seq_len(n_groups),
      n_units = group_sizes,
      mean_shap = mean_shap,
      mean_feature_value = mean_feature_value,
      n_valid_units = n_valid,
      n_groups = n_groups,
      grouping_mode = if (n_valid < raw_unit_threshold) {
        "original_units"
      } else {
        "joint_SHAP_feature_equal_count"
      },
      stringsAsFactors = FALSE
    )
  })
  grouped_features <- Filter(Negate(is.null), grouped_features)
  if (length(grouped_features) < 1L) {
    stop("No valid feature groups could be constructed.")
  }

  result <- do.call(rbind, grouped_features)
  rownames(result) <- NULL
  attr(result, "displayed_features") <- unique(result$feature)
  attr(result, "target_groups") <- target_groups
  attr(result, "min_units_per_group") <- min_units_per_group
  attr(result, "raw_unit_threshold") <- raw_unit_threshold
  result
}


# Draw filled hexagons as grid polygons rather than font glyphs. This keeps the
# summarized-unit symbols identical across PNG/PDF devices and operating systems.
GeomShaptureHexPoint <- ggplot2::ggproto(
  "GeomShaptureHexPoint",
  ggplot2::Geom,
  required_aes = c("x", "y"),
  non_missing_aes = c("size", "colour"),
  default_aes = ggplot2::aes(
    colour = "black",
    size = 1.5,
    alpha = NA
  ),
  draw_key = function(data, params, size) {
    angles <- seq(0, 2 * pi, length.out = 7L)[-7L] + pi / 6
    hex_color <- scales::alpha(data$colour %||% "black", data$alpha %||% NA_real_)
    grid::polygonGrob(
      x = grid::unit(0.5 + 0.38 * cos(angles), "npc"),
      y = grid::unit(0.5 + 0.38 * sin(angles), "npc"),
      gp = grid::gpar(fill = hex_color, col = hex_color, lwd = 0.35)
    )
  },
  draw_panel = function(data, panel_params, coord, na.rm = FALSE) {
    coords <- coord$transform(data, panel_params)
    if (nrow(coords) < 1L) {
      return(grid::nullGrob())
    }

    angles <- seq(0, 2 * pi, length.out = 7L)[-7L] + pi / 6
    n_points <- nrow(coords)
    radius_mm <- rep(coords$size / 2, each = 6L)
    x_offsets <- radius_mm * rep(cos(angles), times = n_points)
    y_offsets <- radius_mm * rep(sin(angles), times = n_points)
    polygon_ids <- rep(seq_len(n_points), each = 6L)
    hex_colors <- scales::alpha(coords$colour, coords$alpha)

    grid::polygonGrob(
      x = grid::unit(rep(coords$x, each = 6L), "native") +
        grid::unit(x_offsets, "mm"),
      y = grid::unit(rep(coords$y, each = 6L), "native") +
        grid::unit(y_offsets, "mm"),
      id = polygon_ids,
      gp = grid::gpar(
        fill = hex_colors,
        col = hex_colors,
        lwd = 0.35
      )
    )
  }
)


geom_shapture_hex_point <- function(
    mapping = NULL,
    data = NULL,
    stat = "identity",
    position = "identity",
    ...,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE
) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomShaptureHexPoint,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


#' Build an equal-count binned SHAP beeswarm
#'
#' @param combined_results Ranked feature table.
#' @param grouped_data Output from [build_shapture_equal_count_groups()].
#' @param task_what_0 Condition 0 label.
#' @param task_what_1 Condition 1 label.
#' @param group_colors Two group colors: condition 0 then condition 1.
#' @param continuous_colors Continuous feature-value palette.
#' @param hexagon_size Positive numeric size for the filled beeswarm hexagons.
#'
#' @return A ggplot object.
#' @export
build_shapture_binned_beeswarm_plot <- function(
    combined_results,
    grouped_data,
    task_what_0,
    task_what_1,
    group_colors,
    continuous_colors,
    hexagon_size = 0.98
) {
  if (!is.numeric(hexagon_size) || length(hexagon_size) != 1L ||
      !is.finite(hexagon_size) || hexagon_size <= 0) {
    stop("hexagon_size must be one positive finite number.")
  }
  required_columns <- c("feature", "mean_shap", "mean_feature_value", "n_units")
  missing_columns <- setdiff(required_columns, colnames(grouped_data))
  if (length(missing_columns) > 0L) {
    stop("grouped_data is missing required columns: ", paste(missing_columns, collapse = ", "), ".")
  }
  features <- unique(as.character(grouped_data$feature))
  features <- features[features %in% as.character(combined_results$feature)]
  if (length(features) < 1L) {
    stop("grouped_data does not contain ranked features from combined_results.")
  }
  grouped_data <- grouped_data[
    as.character(grouped_data$feature) %in% features,
    ,
    drop = FALSE
  ]
  grouped_data$feature <- factor(as.character(grouped_data$feature), levels = rev(features))

  pattern_colors <- model_expression_pattern_colors(
    what_0 = task_what_0,
    what_1 = task_what_1,
    group_colors = group_colors
  )
  pattern_labels <- model_expression_pattern_plot_labels(
    what_0 = task_what_0,
    what_1 = task_what_1
  )
  feature_patterns <- combined_results$Model_expression_pattern[
    match(features, combined_results$feature)
  ]
  feature_colors <- unname(pattern_colors[feature_patterns])
  feature_colors[is.na(feature_colors)] <- "grey55"
  feature_labels <- stats::setNames(
    paste0(
      "<span style='color:", feature_colors, "'>",
      escape_shapture_html(features),
      "</span>"
    ),
    features
  )
  axis_text_y <- if (requireNamespace("ggtext", quietly = TRUE)) {
    ggtext::element_markdown(size = 10)
  } else {
    ggplot2::element_text(colour = "black", size = 10)
  }
  legend_data <- data.frame(
    feature = factor(rep(features[1], length(pattern_labels)), levels = rev(features)),
    mean_shap = 0,
    Model_expression_pattern = names(pattern_labels),
    stringsAsFactors = FALSE
  )

  plot <- ggplot2::ggplot(
    grouped_data,
    ggplot2::aes(x = mean_shap, y = feature)
  ) +
    geom_shapture_hex_point(
      ggplot2::aes(color = mean_feature_value),
      position = ggbeeswarm::position_quasirandom(
        method = "tukey",
        width = 0.3,
        orientation = "y",
        na.rm = TRUE
      ),
      size = hexagon_size,
      alpha = 0.98,
      na.rm = TRUE
    ) +
    ggplot2::scale_color_gradientn(colors = continuous_colors) +
    ggplot2::geom_point(
      data = legend_data,
      ggplot2::aes(x = mean_shap, y = feature, fill = Model_expression_pattern),
      inherit.aes = FALSE,
      shape = 22,
      size = 3,
      alpha = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = pattern_colors,
      breaks = names(pattern_labels),
      labels = unname(pattern_labels),
      drop = FALSE
    ) +
    ggplot2::scale_y_discrete(labels = feature_labels) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::labs(
      title = paste0("Left = ", task_what_0, " | Right = ", task_what_1),
      subtitle = paste0(
        "Each hexagon summarizes an equal-count group of cells formed jointly ",
        "from SHAP and feature values."
      ),
      x = "SHAP value",
      y = "Features ranked by integrated scores",
      color = "Mean feature value",
      fill = "Model\nexpression\npattern"
    ) +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(order = 1),
      fill = ggplot2::guide_legend(
        order = 2,
        override.aes = list(alpha = 1, size = 4)
      )
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y = axis_text_y,
      axis.text.x = ggplot2::element_text(colour = "black"),
      plot.title = ggplot2::element_text(size = 10),
      plot.subtitle = ggplot2::element_text(size = 7.5, color = "#4A5658")
    )
  attr(plot, "displayed_feature_count") <- length(features)
  plot
}


#' Build a SHAP-value dependency plot for one feature
#'
#' Hexagons use shared equal-count groups, while the smooth is fit to the
#' original unbinned unit-level values.
#'
#' @param grouped_data Output from [build_shapture_equal_count_groups()].
#' @param feature Feature to display.
#' @param beeswarm_source_S Unit-averaged SHAP matrix.
#' @param beeswarm_source_X Unit-averaged model-input feature-value matrix or
#'   data frame.
#' @param task_what_0 Condition 0 label.
#' @param task_what_1 Condition 1 label.
#' @param continuous_colors Continuous feature-value palette.
#'
#' @return A ggplot object.
#' @importFrom mgcv s
#' @export
build_shapture_dependency_plot <- function(
    grouped_data,
    feature,
    beeswarm_source_S,
    beeswarm_source_X,
    task_what_0,
    task_what_1,
    continuous_colors
) {
  if (!is.character(feature) || length(feature) != 1L || is.na(feature) || !nzchar(feature)) {
    stop("feature must be one non-empty feature name.")
  }
  required_columns <- c("feature", "mean_shap", "mean_feature_value")
  missing_columns <- setdiff(required_columns, colnames(grouped_data))
  if (length(missing_columns) > 0L) {
    stop("grouped_data is missing required columns: ", paste(missing_columns, collapse = ", "), ".")
  }
  feature_groups <- grouped_data[
    as.character(grouped_data$feature) == feature,
    ,
    drop = FALSE
  ]
  if (nrow(feature_groups) < 1L) {
    stop("grouped_data does not contain feature: ", feature, ".")
  }

  beeswarm_source_S <- as.matrix(beeswarm_source_S)
  beeswarm_source_X <- as.data.frame(beeswarm_source_X, check.names = FALSE)
  if (!feature %in% colnames(beeswarm_source_S) ||
      !feature %in% colnames(beeswarm_source_X)) {
    stop("The requested feature is missing from one or both saved beeswarm inputs.")
  }
  if (nrow(beeswarm_source_S) != nrow(beeswarm_source_X)) {
    stop("beeswarm_source_S and beeswarm_source_X must have the same number of rows.")
  }
  raw_data <- data.frame(
    scaled_feature_value = suppressWarnings(as.numeric(beeswarm_source_X[[feature]])),
    shap_value = suppressWarnings(as.numeric(beeswarm_source_S[, feature]))
  )
  raw_data <- raw_data[
    is.finite(raw_data$scaled_feature_value) & is.finite(raw_data$shap_value),
    ,
    drop = FALSE
  ]
  if (nrow(raw_data) < 1L) {
    stop("No paired finite SHAP and feature values are available for ", feature, ".")
  }

  padded_dependency_range <- function(values, padding_fraction = 0.03) {
    value_range <- range(values, finite = TRUE)
    value_span <- diff(value_range)
    if (!is.finite(value_span) || value_span <= 0) {
      value_span <- max(abs(value_range), 1)
    }
    padding <- value_span * padding_fraction
    c(value_range[1] - padding, value_range[2] + padding)
  }
  dependency_xlim <- padded_dependency_range(feature_groups$mean_shap)
  dependency_ylim <- padded_dependency_range(feature_groups$mean_feature_value)

  prediction_grid <- NULL
  n_unique_values <- length(unique(raw_data$scaled_feature_value))
  if (nrow(raw_data) >= 10L && n_unique_values >= 4L) {
    gam_k <- max(3L, min(12L, n_unique_values - 1L))
    gam_fit <- tryCatch(
      mgcv::gam(
        shap_value ~ s(scaled_feature_value, k = gam_k, bs = "cs"),
        data = raw_data,
        method = "REML"
      ),
      error = function(e) NULL
    )
    if (!is.null(gam_fit)) {
      prediction_grid <- data.frame(
        scaled_feature_value = seq(
          min(raw_data$scaled_feature_value),
          max(raw_data$scaled_feature_value),
          length.out = 401L
        )
      )
      prediction_grid$fit <- as.numeric(stats::predict(
        gam_fit,
        newdata = prediction_grid
      ))
    }
  }

  plot <- ggplot2::ggplot(
    feature_groups,
    ggplot2::aes(x = mean_shap, y = mean_feature_value)
  ) +
    geom_shapture_hex_point(
      ggplot2::aes(color = mean_feature_value),
      size = 3,
      alpha = 0.97
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey25",
      linewidth = 0.45
    ) +
    ggplot2::scale_color_gradientn(
      colors = continuous_colors
    ) +
    ggplot2::coord_cartesian(
      xlim = dependency_xlim,
      ylim = dependency_ylim,
      expand = FALSE
    ) +
    ggplot2::labs(
      title = paste0("Left = ", task_what_0, " | Right = ", task_what_1),
      subtitle = paste0(
        "Each hexagon summarizes an equal-count group of cells formed jointly ",
        "from SHAP and feature values."
      ),
      x = "SHAP value",
      y = paste0("Scaled ", feature, " value"),
      color = "Mean feature value"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      aspect.ratio = 1 / 1.3,
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(
        size = 10,
        margin = ggplot2::margin(b = 2)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 7.5,
        color = "#4A5658",
        margin = ggplot2::margin(b = 12)
      ),
      axis.text.x = ggplot2::element_text(colour = "black"),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8.5)
    )
  if (!is.null(prediction_grid)) {
    plot <- plot + ggplot2::geom_path(
      data = prediction_grid,
      ggplot2::aes(x = fit, y = scaled_feature_value),
      inherit.aes = FALSE,
      color = "#111111",
      linewidth = 1.05
    )
  }
  attr(plot, "feature") <- feature
  attr(plot, "n_groups") <- nrow(feature_groups)
  plot
}


#' Build the top-feature expression density plot
#'
#' @param density_source Saved condition plus expression data.
#' @param combined_results Ranked feature table.
#' @param top_feature_show Number of top features.
#' @param task_what_0 Condition 0 label.
#' @param task_what_1 Condition 1 label.
#' @param group_color_map Named group colors.
#'
#' @return A ggplot object.
#'
#' @details
#' The displayed features are arranged in a near-square grid using
#' `ceiling(sqrt(n_features))` columns. Facet axes remain feature-specific.
#' @export
build_shapture_density_plot <- function(
    density_source,
    combined_results,
    top_feature_show,
    task_what_0,
    task_what_1,
    group_color_map
) {
  data_density <- density_source
  data_density$condition <- ifelse(data_density$condition == 0, task_what_0, task_what_1)
  top_feature_show <- min(top_feature_show, nrow(combined_results))
  if (!is.finite(top_feature_show) || top_feature_show < 1) {
    stop("top_feature_show must be at least 1.")
  }
  top_features <- intersect(as.character(head(combined_results$feature, top_feature_show)), colnames(data_density))
  if (length(top_features) == 0) {
    stop("No saved density features are available for the requested top features.")
  }
  density_ncol <- max(1L, ceiling(sqrt(length(top_features))))
  
  data_density <- reshape2::melt(
    data = data_density,
    id.vars = c("condition"),
    measure.vars = top_features,
    variable.name = "Feat",
    value.name = "Expr"
  )
  data_density$condition <- factor(data_density$condition, levels = c(task_what_0, task_what_1))
  
  ggplot2::ggplot(data_density, ggplot2::aes(x = Expr, fill = condition)) +
    ggplot2::geom_density(ggplot2::aes(colour = condition), alpha = 0.05, linewidth = 1) +
    ggplot2::scale_fill_manual(values = group_color_map, breaks = c(task_what_1, task_what_0)) +
    ggplot2::scale_color_manual(values = group_color_map, breaks = c(task_what_1, task_what_0)) +
    ggplot2::theme_classic() +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 1))) +
    ggplot2::xlab(NULL) +
    ggplot2::facet_wrap(~Feat, ncol = density_ncol, scales = "free") +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = 10, colour = "black", angle = 0, face = "bold"),
      strip.background.x = ggplot2::element_blank()
    )
}


#' Load saved SHAPTURE inputs for editable post-analysis plotting
#'
#' @param result_dir Completed SHAPTURE task folder containing `docs`.
#' @param output_dir Optional folder where figures should be written. Defaults
#'   to `result_dir`.
#' @param group_colors Optional two-color vector overriding the saved group
#'   colors. First color is condition 0/left/others; second is condition
#'   1/right/target.
#' @param continuous_colors Optional continuous color palette overriding the
#'   saved palette.
#'
#' @return A list of tables, saved plotting inputs, context, colors, and output
#'   folders.
#' @export
load_shapture_visualization_inputs <- function(
    result_dir,
    output_dir = NULL,
    group_colors = NULL,
    continuous_colors = NULL
) {
  load_shapture_runtime()
  
  result_dir <- normalizePath(result_dir, winslash = "/", mustWork = TRUE)
  output_root <- if (is.null(output_dir)) {
    result_dir
  } else {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  }
  
  docs_dir <- file.path(result_dir, "docs")
  output_docs_dir <- file.path(output_root, "docs")
  integrated_dir <- file.path(output_root, "Integrated_adaptive_scoring_results")
  shap_only_dir <- file.path(output_root, "SHAP_only")
  dir.create(output_docs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(integrated_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(shap_only_dir, recursive = TRUE, showWarnings = FALSE)
  
  inputs_path <- file.path(docs_dir, "visualization_inputs.rds")
  if (!file.exists(inputs_path)) {
    stop(
      "Missing docs/visualization_inputs.rds in result_dir. ",
      "Rerun the analysis with the updated SHAPTURE pipeline before using post-analysis visualization."
    )
  }
  
  combined_path <- file.path(docs_dir, "combined_results.csv")
  if (!file.exists(combined_path)) {
    stop("Missing docs/combined_results.csv in result_dir.")
  }
  
  inputs <- readRDS(inputs_path)
  combined_results <- utils::read.csv(combined_path, stringsAsFactors = FALSE, check.names = FALSE)
  context <- inputs$context
  config_path <- file.path(docs_dir, "shapture_config_used.csv")
  config_used <- if (file.exists(config_path)) {
    utils::read.csv(config_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    NULL
  }
  get_numeric_context_value <- function(name, default) {
    value <- context[[name]]
    if (length(value) == 1 && !is.null(value) && !is.na(value)) {
      numeric_value <- suppressWarnings(as.numeric(value))
      if (is.finite(numeric_value)) {
        return(numeric_value)
      }
    }
    
    if (!is.null(config_used) && all(c("parameter_name", "parameter_used") %in% colnames(config_used))) {
      idx <- which(config_used$parameter_name == name)
      if (length(idx) > 0) {
        numeric_value <- suppressWarnings(as.numeric(config_used$parameter_used[idx[1]]))
        if (is.finite(numeric_value)) {
          return(numeric_value)
        }
      }
    }
    
    default
  }
  feature_min_log2FC <- get_numeric_context_value("feature_min_log2FC", 0.25)
  shap_expr_cor_cutoff <- get_numeric_context_value("shap_expr_cor_cutoff", 0.25)
  
  group_colors <- validate_color_vector(
    group_colors %||% inputs$group_colors,
    arg_name = "group_colors",
    exact_length = 2
  )
  continuous_colors <- validate_color_vector(
    continuous_colors %||% inputs$continuous_colors,
    arg_name = "continuous_colors",
    min_length = 2
  )
  
  task_what_0 <- context$task_what_0
  task_what_1 <- context$task_what_1
  group_color_map <- build_group_color_map(task_what_0, task_what_1, group_colors)
  cross_patient_support_available <- cross_patient_support_is_available(combined_results)
  if (!is.null(context$cross_patient_support_available) &&
      !identical(isTRUE(context$cross_patient_support_available), cross_patient_support_available)) {
    warning(
      "Saved Cross-patient support context disagrees with combined_results; ",
      "the visualization uses combined_results as the source of truth.",
      call. = FALSE
    )
  }
  cross_patient_support_unavailable_reason <- if (cross_patient_support_available) {
    NA_character_
  } else {
    context$cross_patient_support_unavailable_reason %||% "insufficient patients"
  }
  cross_patient_support_note <- if (cross_patient_support_available) {
    context$cross_patient_support_note %||% "Cross-patient support is available."
  } else {
    "Cross-patient support unavailable: insufficient patients. Final scores use four renormalized components."
  }
  
  list(
    result_dir = result_dir,
    output_root = output_root,
    docs_dir = docs_dir,
    output_docs_dir = output_docs_dir,
    integrated_dir = integrated_dir,
    shap_only_dir = shap_only_dir,
    inputs = inputs,
    combined_results = combined_results,
    context = context,
    task_what_0 = task_what_0,
    task_what_1 = task_what_1,
    group_colors = group_colors,
    continuous_colors = continuous_colors,
    group_color_map = group_color_map,
    cross_patient_support_available = cross_patient_support_available,
    cross_patient_support_unavailable_reason = cross_patient_support_unavailable_reason,
    cross_patient_support_note = cross_patient_support_note,
    min_patients = context$min_patients,
    feature_min_log2FC = feature_min_log2FC,
    shap_expr_cor_cutoff = shap_expr_cor_cutoff
  )
}


#' Build a one-row shapviz object for a saved condition-level SHAP summary
#'
#' @param condition_summary One element of `condition_aggregate_data`.
#' @param baseline SHAP baseline value.
#' @param condition_label Label for the condition.
#'
#' @return A `shapviz` object.
#' @export
build_shapture_condition_shapviz <- function(condition_summary, baseline, condition_label) {
  condition_S_values <- condition_summary$S
  condition_X_values <- condition_summary$X
  aggregate_feature_names <- names(condition_S_values)
  
  condition_S_values[!is.finite(condition_S_values)] <- 0
  condition_X_values[!is.finite(condition_X_values)] <- 0
  
  structure(
    list(
      S = matrix(condition_S_values, nrow = 1, dimnames = list(paste0("condition_", condition_label), aggregate_feature_names)),
      X = as.data.frame(
        matrix(condition_X_values, nrow = 1, dimnames = list(paste0("condition_", condition_label), aggregate_feature_names)),
        check.names = FALSE
      ),
      baseline = baseline,
      S_inter = NULL
    ),
    class = "shapviz"
  )
}


#' Extract feature names from shapviz waterfall/force plot data
#'
#' @param plot A shapviz ggplot object.
#' @param feature_names Valid feature names from the source SHAP object.
#'
#' @return Character vector of feature names present in `feature_names`.
#' @noRd
extract_shapture_shapviz_plot_features <- function(plot, feature_names, keep_alignment = FALSE) {
  row_features <- row.names(plot$data)
  row_features <- sub("^[^.]+\\.", "", row_features)
  row_features[!row_features %in% feature_names] <- NA_character_
  
  label_features <- sub("\\s*=.*$", "", as.character(plot$data$label))
  label_features[!label_features %in% feature_names] <- NA_character_
  resolved_features <- ifelse(!is.na(row_features), row_features, label_features)
  
  if (isTRUE(keep_alignment)) {
    return(resolved_features)
  }
  resolved_features[!is.na(resolved_features)]
}


#' Add average feature values to selected shapviz waterfall/force labels
#'
#' @param plot A shapviz ggplot object.
#' @param condition_shapviz One-row shapviz object.
#' @param top_feature_label Number of highest absolute-SHAP feature labels to keep.
#'
#' @return The modified ggplot object.
#' @noRd
add_shapture_average_labels <- function(plot, condition_shapviz, top_feature_label = Inf) {
  row_features <- extract_shapture_shapviz_plot_features(
    plot,
    colnames(condition_shapviz$X),
    keep_alignment = TRUE
  )
  feature_row_idx <- which(!is.na(row_features))
  if (length(feature_row_idx) == 0) {
    return(plot)
  }
  
  label_count <- if (is.null(top_feature_label) || !is.finite(top_feature_label)) {
    length(feature_row_idx)
  } else {
    min(max(0L, as.integer(floor(top_feature_label))), length(feature_row_idx))
  }
  if ("S" %in% colnames(plot$data)) {
    feature_rank <- order(abs(plot$data$S[feature_row_idx]), decreasing = TRUE, na.last = TRUE)
  } else {
    feature_rank <- rev(seq_along(feature_row_idx))
  }
  label_idx <- feature_row_idx[head(feature_rank, label_count)]
  blank_idx <- setdiff(feature_row_idx, label_idx)
  if (length(blank_idx) > 0) {
    plot$data$label[blank_idx] <- ""
  }
  if (length(label_idx) > 0) {
    features_to_label <- row_features[label_idx]
    average_values <- round(as.numeric(condition_shapviz$X[1, features_to_label, drop = TRUE]), 2)
    plot$data$label[label_idx] <- paste0(features_to_label, " (", average_values, ")")
  }
  
  plot
}


#' Build a condition-level SHAP waterfall plot
#'
#' @param condition_shapviz One-row shapviz object.
#' @param condition_label Label for the condition.
#' @param max_display Number of features to display.
#' @param top_feature_label Number of highest absolute-SHAP feature labels to keep.
#'
#' @return A ggplot object.
#' @export
build_shapture_condition_waterfall_plot <- function(condition_shapviz, condition_label, max_display, top_feature_label = max_display) {
  waterfall_plot <- shapviz::sv_waterfall(condition_shapviz, row_id = 1, max_display = max_display) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 11)) +
    ggplot2::ggtitle(paste0(
      "Condition-level mean SHAP waterfall summary", "\n",
      "Average feature profile and average SHAP contribution for condition ", condition_label
    )) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11))
  
  waterfall_plot <- add_shapture_average_labels(waterfall_plot, condition_shapviz, top_feature_label)
  
  sanitize_shapviz_plot(waterfall_plot)
}


#' Build a condition-level SHAP force plot
#'
#' @param condition_shapviz One-row shapviz object.
#' @param condition_label Label for the condition.
#' @param max_display Number of features to display.
#' @param top_feature_label Number of highest absolute-SHAP feature labels to keep.
#'
#' @return A ggplot object.
#' @export
build_shapture_condition_force_plot <- function(condition_shapviz, condition_label, max_display, top_feature_label = max_display) {
  force_plot <- shapviz::sv_force(condition_shapviz, row_id = 1, max_display = max_display) +
    ggplot2::ggtitle(paste0(
      "Condition-level mean SHAP force summary", "\n",
      "Average feature profile and average SHAP contribution for condition ", condition_label
    )) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 8))
  
  force_plot <- add_shapture_average_labels(force_plot, condition_shapviz, top_feature_label)
  
  sanitize_shapviz_plot(force_plot)
}


#' Write average expression values used in SHAP-only condition plots
#'
#' @param condition_shapviz One-row shapviz object.
#' @param plot Plot returned by `build_shapture_condition_waterfall_plot()`.
#' @param output_path CSV output path.
#'
#' @return Invisibly returns `output_path`.
#' @export
write_shapture_condition_average_expression <- function(condition_shapviz, plot, output_path) {
  features <- extract_shapture_shapviz_plot_features(plot, colnames(condition_shapviz$X))
  average_values <- round(as.numeric(condition_shapviz$X[1, features, drop = TRUE]), 2)
  utils::write.csv(
    data.frame(average = rev(average_values)),
    output_path
  )
  invisible(output_path)
}
