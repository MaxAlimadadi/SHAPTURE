#####
# Editable SHAPTURE multiple-folder visualization pipeline
#
# This script does not rerun the model. It finds every immediate subfolder with
# completed SHAPTURE visualization inputs, then recreates the figures in each
# folder. It works with condition-vs-condition and one-vs-others results.
#####

library(SHAPTURE)
library(RColorBrewer)

options(ggrepel.max.overlaps = Inf)

#####
# Select the parent directory containing completed SHAPTURE result folders.
#####
analysis_root <- "path/to/folder_containing_shapture_result_folders"

# First color = condition_0/left/others; second color = condition_1/right/target.
group_colors <- c("darkblue", "darkred")

# Any valid R color vector with at least 2 colors.
continuous_colors <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))

# Save generated ggplot objects to docs/plots.rds.
save_plot_objects <- FALSE

#####
# Feature display settings
#####

# Number of top-ranked features used in the main figures.
top_feature_show <- 50

# Number of displayed features labeled in plots that support selective labels.
top_feature_label <- 15

# Number of features used in the larger heatmap and decomposition figures.
top_feature_show_extended <- 100

#####
# Automatic figure-size growth
#####

# Global increase in width and height for each displayed feature, in inches.
# Plot-specific base sizes and safety limits are handled by the package helper.
# Density plots instead use an automatic near-square grid with fixed-size panels.
width_per_feature <- 0.01
height_per_feature <- 0.13

# To use a fixed size for one figure, replace calculated_size$width or
# calculated_size$height with a number directly in that figure's pdf()/png().

if (!is.numeric(top_feature_show) || length(top_feature_show) != 1L ||
    !is.finite(top_feature_show) || top_feature_show < 1L ||
    top_feature_show != floor(top_feature_show)) {
  stop("top_feature_show must be a positive integer.")
}
if (!is.numeric(top_feature_label) || length(top_feature_label) != 1L ||
    !is.finite(top_feature_label) || top_feature_label < 0L ||
    top_feature_label != floor(top_feature_label)) {
  stop("top_feature_label must be a non-negative integer.")
}
if (!is.numeric(top_feature_show_extended) || length(top_feature_show_extended) != 1L ||
    !is.finite(top_feature_show_extended) || top_feature_show_extended < 1L ||
    top_feature_show_extended != floor(top_feature_show_extended)) {
  stop("top_feature_show_extended must be a positive integer.")
}

result_folders <- find_shapture_result_folders(analysis_root)

for (result_dir in result_folders) {
  message("Creating figures for ", basename(result_dir))

  #####
  # Load saved tables and plot inputs.
  #####
  viz <- load_shapture_visualization_inputs(
    result_dir = result_dir,
    group_colors = group_colors,
    continuous_colors = continuous_colors
  )

  combined_results <- viz$combined_results
  inputs <- viz$inputs

  if (nrow(combined_results) < 1L) {
    stop("combined_results does not contain any retained features in ", result_dir, ".")
  }
  top_feature_show_actual <- min(top_feature_show, nrow(combined_results))
  top_feature_label_actual <- min(top_feature_label, top_feature_show_actual)
  top_feature_show_extended_actual <- min(
    top_feature_show_extended,
    nrow(combined_results)
  )

  plot_objects <- list()

  #####
  # ROC plot
  #####
  ROC_plot_size <- calculate_shapture_figure_size(
    figure_type = "roc",
    width_per_feature = width_per_feature,
    height_per_feature = height_per_feature
  )

  ROC_plot <- build_shapture_roc_plot(
    ROC_all = inputs$roc$ROC_all,
    AUC_all = inputs$roc$AUC_all,
    mean_auc = inputs$roc$summary$mean_iteration_auc
  )
  plot_objects$ROC_plot <- ROC_plot

  pdf(
    file.path(viz$integrated_dir, "ROC_plot.pdf"),
    width = ROC_plot_size$width,
    height = ROC_plot_size$height
  )
  print(ROC_plot)
  dev.off()

  #####
  # Feature importance plot
  #####
  feature_importance_size <- calculate_shapture_figure_size(
    figure_type = "feature_importance",
    n_features = top_feature_show_actual,
    width_per_feature = width_per_feature,
    height_per_feature = height_per_feature
  )

  feature_importance_plot <- build_shapture_feature_importance_plot(
    combined_results = combined_results,
    top_feature_show = top_feature_show_actual,
    task_what_0 = viz$task_what_0,
    task_what_1 = viz$task_what_1,
    group_colors = viz$group_colors
  )
  plot_objects$feature_importance <- feature_importance_plot

  pdf(
    file.path(viz$integrated_dir, paste0("feature_importance_top.", top_feature_show_actual, ".pdf")),
    width = feature_importance_size$width,
    height = feature_importance_size$height
  )
  print(feature_importance_plot)
  dev.off()

  #####
  # SHAP-only condition waterfall and force plots
  #####
  waterfall_plot_size <- calculate_shapture_figure_size(
    figure_type = "waterfall",
    n_features = top_feature_show_actual,
    width_per_feature = width_per_feature,
    height_per_feature = height_per_feature
  )
  force_plot_size <- calculate_shapture_figure_size(
    figure_type = "force",
    n_features = top_feature_show_actual,
    width_per_feature = width_per_feature,
    height_per_feature = height_per_feature
  )

  condition_aggregate_data <- inputs$shap$condition_aggregate_data
  aggregate_baseline <- attr(condition_aggregate_data, "baseline")
  if (is.null(aggregate_baseline) || !is.finite(aggregate_baseline)) {
    aggregate_baseline <- inputs$shap$beeswarm_baseline
  }

  for (what_condition in sort(as.numeric(names(condition_aggregate_data)))) {
    condition_key <- as.character(what_condition)
    condition_label <- ifelse(what_condition == 0, viz$task_what_0, viz$task_what_1)

    condition_shapviz <- build_shapture_condition_shapviz(
      condition_summary = condition_aggregate_data[[condition_key]],
      baseline = aggregate_baseline,
      condition_label = condition_label
    )

    waterfall_plot <- build_shapture_condition_waterfall_plot(
      condition_shapviz = condition_shapviz,
      condition_label = condition_label,
      max_display = top_feature_show_actual,
      top_feature_label = top_feature_label_actual
    )
    plot_objects[[paste0("waterfall_", condition_label)]] <- waterfall_plot

    pdf(
      file.path(viz$shap_only_dir, paste0("waterfall.plot_aggregated_condition.", condition_label, ".pdf")),
      width = waterfall_plot_size$width,
      height = waterfall_plot_size$height
    )
    print(waterfall_plot)
    dev.off()

    force_plot <- build_shapture_condition_force_plot(
      condition_shapviz = condition_shapviz,
      condition_label = condition_label,
      max_display = top_feature_show_actual,
      top_feature_label = top_feature_label_actual
    )
    plot_objects[[paste0("force_", condition_label)]] <- force_plot

    pdf(
      file.path(viz$shap_only_dir, paste0("force_plot_aggregated_condition.", condition_label, ".pdf")),
      width = force_plot_size$width,
      height = force_plot_size$height
    )
    print(force_plot)
    dev.off()

    write_shapture_condition_average_expression(
      condition_shapviz = condition_shapviz,
      plot = waterfall_plot,
      output_path = file.path(viz$output_docs_dir, paste0("average.exp_aggregated_condition.", condition_label, ".csv"))
    )
  }

  #####
  # Beeswarm plot
  #####
  beeswarm_source_S <- as.matrix(inputs$shap$beeswarm_source_S)
  beeswarm_source_X <- as.data.frame(
    inputs$shap$beeswarm_source_X,
    check.names = FALSE
  )
  requested_beeswarm_features <- as.character(utils::head(
    combined_results$feature,
    top_feature_show_actual
  ))
  beeswarm_features_actual <- requested_beeswarm_features[
    requested_beeswarm_features %in% colnames(beeswarm_source_S) &
      requested_beeswarm_features %in% colnames(beeswarm_source_X)
  ]
  if (length(beeswarm_features_actual) < 1L) {
    stop(
      "No requested top features are available in both saved beeswarm inputs for ",
      result_dir,
      "."
    )
  }
  beeswarm_feature_count <- length(beeswarm_features_actual)

  beeswarm_plot_size <- calculate_shapture_figure_size(
    figure_type = "beeswarm",
    n_features = beeswarm_feature_count,
    width_per_feature = width_per_feature,
    height_per_feature = height_per_feature
  )

  beeswarm_plot <- build_shapture_beeswarm_plot(
    combined_results = combined_results,
    beeswarm_source_S = beeswarm_source_S,
    beeswarm_source_X = beeswarm_source_X,
    top_feature_show = top_feature_show_actual,
    task_what_0 = viz$task_what_0,
    task_what_1 = viz$task_what_1,
    group_colors = viz$group_colors,
    continuous_colors = viz$continuous_colors
  )
  plot_objects$beeswarm <- beeswarm_plot

  png(
    filename = file.path(viz$integrated_dir, paste0("beeswarm_plot_top.", beeswarm_feature_count, ".features.png")),
    width = beeswarm_plot_size$width,
    height = beeswarm_plot_size$height,
    units = beeswarm_plot_size$units,
    res = 300
  )
  print(beeswarm_plot)
  dev.off()

  #####
  # Equal-count binned beeswarm and SHAP-value dependency plots
  #####
  equal_count_groups <- build_shapture_equal_count_groups(
    beeswarm_source_S = beeswarm_source_S,
    beeswarm_source_X = beeswarm_source_X,
    features = beeswarm_features_actual,
    target_groups = 500L,
    min_units_per_group = 10L,
    raw_unit_threshold = 500L
  )
  equal_count_features <- attr(equal_count_groups, "displayed_features")
  equal_count_feature_count <- length(equal_count_features)

  write.csv(
    equal_count_groups,
    file.path(viz$output_docs_dir, "equal_count_shap_group_summary.csv"),
    row.names = FALSE
  )

  binned_beeswarm_plot <- build_shapture_binned_beeswarm_plot(
    combined_results = combined_results,
    grouped_data = equal_count_groups,
    task_what_0 = viz$task_what_0,
    task_what_1 = viz$task_what_1,
    group_colors = viz$group_colors,
    continuous_colors = viz$continuous_colors
  )
  plot_objects$binned_beeswarm <- binned_beeswarm_plot

  binned_beeswarm_plot_size <- calculate_shapture_figure_size(
    figure_type = "beeswarm",
    n_features = equal_count_feature_count,
    width_per_feature = width_per_feature,
    height_per_feature = height_per_feature
  )
  png(
    filename = file.path(
      viz$integrated_dir,
      paste0("binned_beeswarm_plot_top.", equal_count_feature_count, ".features.png")
    ),
    width = binned_beeswarm_plot_size$width,
    height = binned_beeswarm_plot_size$height,
    units = binned_beeswarm_plot_size$units,
    res = 300
  )
  print(binned_beeswarm_plot)
  dev.off()

  dependency_dir <- file.path(
    viz$integrated_dir,
    "shap_value_dependency_plots"
  )
  dir.create(dependency_dir, recursive = TRUE, showWarnings = FALSE)
  old_dependency_files <- list.files(
    dependency_dir,
    pattern = "^shap_value_dependency_.*[.]pdf$",
    full.names = TRUE
  )
  if (length(old_dependency_files) > 0L) {
    unlink(old_dependency_files)
  }
  dependency_plot_size <- calculate_shapture_figure_size(
    figure_type = "dependency",
    n_features = 1L,
    width_per_feature = width_per_feature,
    height_per_feature = height_per_feature
  )
  safe_dependency_names <- gsub("[^A-Za-z0-9._-]+", "_", equal_count_features)
  empty_dependency_names <- !nzchar(safe_dependency_names)
  safe_dependency_names[empty_dependency_names] <- paste0(
    "feature_",
    which(empty_dependency_names)
  )
  safe_dependency_names <- make.unique(safe_dependency_names, sep = "_")
  plot_objects$shap_value_dependency <- list()
  for (feature_index in seq_along(equal_count_features)) {
    dependency_feature <- equal_count_features[feature_index]
    dependency_plot <- build_shapture_dependency_plot(
      grouped_data = equal_count_groups,
      feature = dependency_feature,
      beeswarm_source_S = beeswarm_source_S,
      beeswarm_source_X = beeswarm_source_X,
      task_what_0 = viz$task_what_0,
      task_what_1 = viz$task_what_1,
      continuous_colors = viz$continuous_colors
    )
    plot_objects$shap_value_dependency[[dependency_feature]] <- dependency_plot
    dependency_path <- file.path(
      dependency_dir,
      paste0("shap_value_dependency_", safe_dependency_names[feature_index], ".pdf")
    )
    if (isTRUE(capabilities("cairo"))) {
      grDevices::cairo_pdf(
        filename = dependency_path,
        width = dependency_plot_size$width,
        height = dependency_plot_size$height
      )
    } else {
      grDevices::pdf(
        file = dependency_path,
        width = dependency_plot_size$width,
        height = dependency_plot_size$height
      )
    }
    print(dependency_plot)
    dev.off()
  }

  #####
  # Density plot
  #####
  density_feature_count <- length(intersect(
    as.character(utils::head(combined_results$feature, top_feature_show_actual)),
    colnames(inputs$density$source)
  ))
  density_plot_size <- calculate_shapture_figure_size(
    figure_type = "density",
    n_features = density_feature_count
  )

  density_plot <- build_shapture_density_plot(
    density_source = inputs$density$source,
    combined_results = combined_results,
    top_feature_show = top_feature_show_actual,
    task_what_0 = viz$task_what_0,
    task_what_1 = viz$task_what_1,
    group_color_map = viz$group_color_map
  )
  plot_objects$density <- density_plot

  pdf(
    file.path(viz$integrated_dir, "density_plot.pdf"),
    width = density_plot_size$width,
    height = density_plot_size$height
  )
  print(density_plot)
  dev.off()

  #####
  # Score decomposition, top features
  #####
  score_decomposition_size <- calculate_shapture_figure_size(
    figure_type = "score_decomposition",
    n_features = top_feature_show_actual,
    width_per_feature = width_per_feature,
    height_per_feature = height_per_feature
  )

  score_decomposition_plot <- plot_score_decomposition(
    df = combined_results,
    top_n_to_show = top_feature_show_actual
  )
  plot_objects$score_decomposition <- score_decomposition_plot

  pdf(
    file.path(viz$integrated_dir, paste0("score_decomposition_top.", top_feature_show_actual, ".pdf")),
    width = score_decomposition_size$width,
    height = score_decomposition_size$height
  )
  print(score_decomposition_plot)
  dev.off()

  #####
  # Score decomposition, top selected features
  #####
  if (top_feature_show_extended_actual > top_feature_show_actual) {
    score_decomposition_extended_size <- calculate_shapture_figure_size(
      figure_type = "score_decomposition",
      n_features = top_feature_show_extended_actual,
      width_per_feature = width_per_feature,
      height_per_feature = height_per_feature
    )

    score_decomposition_all_plot <- plot_score_decomposition(
      df = combined_results,
      top_n_to_show = top_feature_show_extended_actual
    ) +
      ggplot2::ggtitle(paste0("Score decomposition\nTop ", top_feature_show_extended_actual, " features")) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5))
    plot_objects$score_decomposition_all <- score_decomposition_all_plot

    pdf(
      file.path(viz$integrated_dir, paste0("score_decomposition_top.", top_feature_show_extended_actual, ".pdf")),
      width = score_decomposition_extended_size$width,
      height = score_decomposition_extended_size$height
    )
    print(score_decomposition_all_plot)
    dev.off()
  }

  #####
  # Component heatmap, top features
  #####
  component_heatmap_size <- calculate_shapture_figure_size(
    figure_type = "component_heatmap",
    n_features = top_feature_show_actual,
    width_per_feature = width_per_feature,
    height_per_feature = height_per_feature
  )
  five_component_heatmap_top_tile_color <- "white" # use NA to remove tile borders
  five_component_heatmap_top_tile_linewidth <- 0.25
  component_count_label <- if (viz$cross_patient_support_available) "Five" else "Four"
  component_heatmap_file_prefix <- if (viz$cross_patient_support_available) {
    "five_component_heatmap_top."
  } else {
    "four_component_heatmap_top."
  }

  five_component_heatmap_top_plot <- plot_five_component_heatmap(
    df = combined_results,
    top_n_to_show = top_feature_show_actual,
    continuous_colors = viz$continuous_colors,
    title = paste0(component_count_label, "-component support matrix\nTop ", top_feature_show_actual, " retained features"),
    tile_color = five_component_heatmap_top_tile_color,
    tile_linewidth = five_component_heatmap_top_tile_linewidth
  )
  plot_objects$five_component_heatmap_top <- five_component_heatmap_top_plot

  pdf(
    file.path(viz$integrated_dir, paste0(component_heatmap_file_prefix, top_feature_show_actual, ".pdf")),
    width = component_heatmap_size$width,
    height = component_heatmap_size$height
  )
  print(five_component_heatmap_top_plot)
  dev.off()

  #####
  # Component heatmap, top selected features
  #####
  if (top_feature_show_extended_actual > top_feature_show_actual) {
    component_heatmap_extended_size <- calculate_shapture_figure_size(
      figure_type = "component_heatmap",
      n_features = top_feature_show_extended_actual,
      width_per_feature = width_per_feature,
      height_per_feature = height_per_feature
    )
    five_component_heatmap_all_tile_color <- NA # use "white" for visible tile borders
    five_component_heatmap_all_tile_linewidth <- 0.05

    five_component_heatmap_all_plot <- plot_five_component_heatmap(
      df = combined_results,
      top_n_to_show = top_feature_show_extended_actual,
      continuous_colors = viz$continuous_colors,
      title = paste0(component_count_label, "-component support matrix\nTop ", top_feature_show_extended_actual, " features"),
      tile_color = five_component_heatmap_all_tile_color,
      tile_linewidth = five_component_heatmap_all_tile_linewidth
    )
    plot_objects$five_component_heatmap_all <- five_component_heatmap_all_plot

    pdf(
      file.path(viz$integrated_dir, paste0(component_heatmap_file_prefix, top_feature_show_extended_actual, ".pdf")),
      width = component_heatmap_extended_size$width,
      height = component_heatmap_extended_size$height
    )
    print(five_component_heatmap_all_plot)
    dev.off()
  }

  #####
  # Evidence heatmap
  #####
  evidence_heatmap_size <- calculate_shapture_figure_size(
    figure_type = "evidence_heatmap",
    n_features = top_feature_show_actual,
    width_per_feature = width_per_feature,
    height_per_feature = height_per_feature
  )

  evidence_heatmap_plot <- plot_evidence_heatmap(
    df = combined_results,
    top_n_to_show = top_feature_show_actual,
    continuous_colors = viz$continuous_colors
  )
  plot_objects$evidence_heatmap <- evidence_heatmap_plot

  pdf(
    file.path(viz$integrated_dir, paste0("evidence_heatmap_top.", top_feature_show_actual, ".pdf")),
    width = evidence_heatmap_size$width,
    height = evidence_heatmap_size$height
  )
  print(evidence_heatmap_plot)
  dev.off()

  #####
  # SHAP effect versus expression log2FC
  #####
  shap_vs_expression_size <- calculate_shapture_figure_size(
    figure_type = "shap_vs_expression",
    n_features = top_feature_show_actual,
    width_per_feature = width_per_feature,
    height_per_feature = height_per_feature
  )
  shap_vs_expression_dot_label_size <- 3

  shap_vs_expression_plot <- plot_shap_vs_expression(
    df = combined_results,
    top_n_to_show = top_feature_show_actual,
    top_n_to_name = top_feature_label_actual,
    what_0 = viz$task_what_0,
    what_1 = viz$task_what_1,
    group_colors = viz$group_colors,
    expression_log2FC_cutoff = viz$feature_min_log2FC,
    shap_expr_cor_cutoff = viz$shap_expr_cor_cutoff,
    dot_label_size = shap_vs_expression_dot_label_size
  )
  plot_objects$shap_vs_expression <- shap_vs_expression_plot

  pdf(
    file.path(viz$integrated_dir, "shap_vs_expression.pdf"),
    width = shap_vs_expression_size$width,
    height = shap_vs_expression_size$height
  )
  print(shap_vs_expression_plot)
  dev.off()

  #####
  # Value-SHAP relationship by feature
  #####
  valueSHAP_by_feature_size <- calculate_shapture_figure_size(
    figure_type = "value_shap",
    n_features = top_feature_show_actual,
    width_per_feature = width_per_feature,
    height_per_feature = height_per_feature
  )
  valueSHAP_by_feature_dot_label_size <- 3

  valueSHAP_by_feature_plot <- plot_valueSHAP_by_feature(
    df = combined_results,
    top_n_to_show = top_feature_show_actual,
    top_n_to_name = top_feature_label_actual,
    what_0 = viz$task_what_0,
    what_1 = viz$task_what_1,
    continuous_colors = viz$continuous_colors,
    group_colors = viz$group_colors,
    expression_log2FC_cutoff = viz$feature_min_log2FC,
    shap_expr_cor_cutoff = viz$shap_expr_cor_cutoff,
    dot_label_size = valueSHAP_by_feature_dot_label_size
  )
  plot_objects$valueSHAP_by_feature <- valueSHAP_by_feature_plot

  pdf(
    file.path(viz$integrated_dir, "valueSHAP_by_feature.pdf"),
    width = valueSHAP_by_feature_size$width,
    height = valueSHAP_by_feature_size$height
  )
  print(valueSHAP_by_feature_plot)
  dev.off()

  if (isTRUE(save_plot_objects)) {
    saveRDS(plot_objects, file.path(viz$output_docs_dir, "plots.rds"))
  }


}
