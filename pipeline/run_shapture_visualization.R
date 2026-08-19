#####
# Editable SHAPTURE visualization pipeline
#
# This script does not rerun the model. It reads a completed analysis folder,
# calls the individual plot functions, and saves each figure with visible
# width/height settings that you can edit directly here.
#####

library(SHAPTURE)
library(RColorBrewer)

options(ggrepel.max.overlaps = Inf)

#####
# Select the completed result folder.
# This folder must contain docs/combined_results.csv and docs/visualization_inputs.rds.
#####
result_dir <- "path/to/completed_shapture_result_folder"

# First color = condition_0/left/others; second color = condition_1/right/target.
group_colors <- c("darkblue", "darkred")

# Any valid R color vector with at least 2 colors.
continuous_colors <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))

# Save generated ggplot objects to docs/plots.rds.
save_plot_objects <- FALSE

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

plot_objects <- list()

#####
# ROC plot
#####
ROC_plot_width <- 7.5
ROC_plot_height <- 4

ROC_plot <- build_shapture_roc_plot(
  ROC_all = inputs$roc$ROC_all,
  AUC_all = inputs$roc$AUC_all,
  mean_auc = inputs$roc$summary$mean_iteration_auc
)
ROC_plot <- add_inference_warning(ROC_plot, viz$inference_reliable, viz$min_patients)
plot_objects$ROC_plot <- ROC_plot

pdf(file.path(viz$integrated_dir, "ROC_plot.pdf"), width = ROC_plot_width, height = ROC_plot_height)
print(ROC_plot)
dev.off()

#####
# Feature importance plot
#####
top_feature_show <- 50
feature_importance_width <- 3
feature_importance_height <- 6

feature_importance_plot <- build_shapture_feature_importance_plot(
  combined_results = combined_results,
  top_feature_show = top_feature_show
)
feature_importance_plot <- add_inference_warning(feature_importance_plot, viz$inference_reliable, viz$min_patients)
plot_objects$feature_importance <- feature_importance_plot

pdf(
  file.path(viz$integrated_dir, paste0("feature_importance_top.", top_feature_show, ".pdf")),
  width = feature_importance_width,
  height = feature_importance_height
)
print(feature_importance_plot)
dev.off()

#####
# SHAP-only condition waterfall and force plots
#####
top_feature_show <- 50
top_feature_label <- 15
waterfall_plot_width <- 7
waterfall_plot_height <- 10
force_plot_width <- 10
force_plot_height <- 4

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
    max_display = top_feature_show,
    top_feature_label = top_feature_show
  )
  plot_objects[[paste0("waterfall_", condition_label)]] <- waterfall_plot
  
  pdf(
    file.path(viz$shap_only_dir, paste0("waterfall.plot_aggregated_condition.", condition_label, ".pdf")),
    width = waterfall_plot_width,
    height = waterfall_plot_height
  )
  print(waterfall_plot)
  dev.off()
  
  force_plot <- build_shapture_condition_force_plot(
    condition_shapviz = condition_shapviz,
    condition_label = condition_label,
    max_display = top_feature_show,
    top_feature_label = top_feature_label
  )
  plot_objects[[paste0("force_", condition_label)]] <- force_plot
  
  pdf(
    file.path(viz$shap_only_dir, paste0("force_plot_aggregated_condition.", condition_label, ".pdf")),
    width = force_plot_width,
    height = force_plot_height
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
top_feature_show <- 50
beeswarm_plot_width <- 17
beeswarm_plot_height <- 25

beeswarm_plot <- build_shapture_beeswarm_plot(
  combined_results = combined_results,
  beeswarm_source_S = as.matrix(inputs$shap$beeswarm_source_S),
  beeswarm_source_X = as.data.frame(inputs$shap$beeswarm_source_X, check.names = FALSE),
  top_feature_show = top_feature_show,
  task_what_0 = viz$task_what_0,
  task_what_1 = viz$task_what_1,
  group_colors = viz$group_colors,
  continuous_colors = viz$continuous_colors
)
beeswarm_plot <- add_inference_warning(beeswarm_plot, viz$inference_reliable, viz$min_patients)
plot_objects$beeswarm <- beeswarm_plot

png(
  filename = file.path(viz$integrated_dir, paste0("beeswarm_plot_top.", top_feature_show, ".features.png")),
  width = beeswarm_plot_width,
  height = beeswarm_plot_height,
  units = "cm",
  res = 300
)
print(beeswarm_plot)
dev.off()

#####
# Density plot
#####
top_feature_show <- 50
density_plot_width <- 15
density_plot_height <- 18

density_plot <- build_shapture_density_plot(
  density_source = inputs$density$source,
  combined_results = combined_results,
  top_feature_show = top_feature_show,
  task_what_0 = viz$task_what_0,
  task_what_1 = viz$task_what_1,
  group_color_map = viz$group_color_map
)
density_plot <- add_inference_warning(density_plot, viz$inference_reliable, viz$min_patients)
plot_objects$density <- density_plot

pdf(file.path(viz$integrated_dir, "density_plot.pdf"), width = density_plot_width, height = density_plot_height)
print(density_plot)
dev.off()

#####
# Score decomposition, top features
#####
top_feature_show <- 50
score_decomposition_width <- 6
score_decomposition_height <- 7

score_decomposition_plot <- plot_score_decomposition(
  df = combined_results,
  top_n_to_show = top_feature_show
)
score_decomposition_plot <- add_inference_warning(score_decomposition_plot, viz$inference_reliable, viz$min_patients)
plot_objects$score_decomposition <- score_decomposition_plot

pdf(
  file.path(viz$integrated_dir, paste0("score_decomposition_top.", top_feature_show, ".pdf")),
  width = score_decomposition_width,
  height = score_decomposition_height
)
print(score_decomposition_plot)
dev.off()

#####
# Score decomposition, top selected features
#####
top_feature_show <- 100 # to select all nrow(combined_results)
score_decomposition_all_width <- 6
score_decomposition_all_height <- 10

score_decomposition_all_plot <- plot_score_decomposition(
  df = combined_results,
  top_n_to_show = top_feature_show 
) +
  ggplot2::ggtitle(paste0("Score decomposition\nTop ", top_feature_show, " features")) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5))
score_decomposition_all_plot <- add_inference_warning(score_decomposition_all_plot, viz$inference_reliable, viz$min_patients)
plot_objects$score_decomposition_all <- score_decomposition_all_plot

pdf(
  file.path(viz$integrated_dir, paste0("score_decomposition_top.", top_feature_show, ".pdf")),
  width = score_decomposition_all_width,
  height = score_decomposition_all_height
)
print(score_decomposition_all_plot)
dev.off()

#####
# Five-component heatmap, top features
#####
top_feature_show <- 50
five_component_heatmap_top_width <- 5.8
five_component_heatmap_top_height <- 10
five_component_heatmap_top_tile_color <- "white" # use NA to remove tile borders
five_component_heatmap_top_tile_linewidth <- 0.25

five_component_heatmap_top_plot <- plot_five_component_heatmap(
  df = combined_results,
  top_n_to_show = top_feature_show,
  continuous_colors = viz$continuous_colors,
  title = paste0("Five-component support matrix\nTop ", top_feature_show, " retained features"),
  tile_color = five_component_heatmap_top_tile_color,
  tile_linewidth = five_component_heatmap_top_tile_linewidth
)
five_component_heatmap_top_plot <- add_inference_warning(five_component_heatmap_top_plot, viz$inference_reliable, viz$min_patients)
plot_objects$five_component_heatmap_top <- five_component_heatmap_top_plot

pdf(
  file.path(viz$integrated_dir, paste0("five_component_heatmap_top.", top_feature_show, ".pdf")),
  width = five_component_heatmap_top_width,
  height = five_component_heatmap_top_height
)
print(five_component_heatmap_top_plot)
dev.off()

#####
# Five-component heatmap, top selected features
#####
top_feature_show <- 100 # to select all nrow(combined_results)
five_component_heatmap_all_width <- 6
five_component_heatmap_all_height <- 10
five_component_heatmap_all_tile_color <- NA # use "white" for visible tile borders
five_component_heatmap_all_tile_linewidth <- 0.05

five_component_heatmap_all_plot <- plot_five_component_heatmap(
  df = combined_results,
  top_n_to_show = top_feature_show,
  continuous_colors = viz$continuous_colors,
  title = paste0("Five-component support matrix\nTop ", top_feature_show, " features"),
  tile_color = five_component_heatmap_all_tile_color,
  tile_linewidth = five_component_heatmap_all_tile_linewidth
)
five_component_heatmap_all_plot <- add_inference_warning(five_component_heatmap_all_plot, viz$inference_reliable, viz$min_patients)
plot_objects$five_component_heatmap_all <- five_component_heatmap_all_plot

pdf(
  file.path(viz$integrated_dir, paste0("five_component_heatmap_top.", top_feature_show, ".pdf")),
  width = five_component_heatmap_all_width,
  height = five_component_heatmap_all_height
)
print(five_component_heatmap_all_plot)
dev.off()

#####
# Evidence heatmap
#####
top_feature_show <- 50
evidence_heatmap_width <- 4.05
evidence_heatmap_height <- 7

evidence_heatmap_plot <- plot_evidence_heatmap(
  df = combined_results,
  top_n_to_show = top_feature_show,
  continuous_colors = viz$continuous_colors
)
evidence_heatmap_plot <- add_inference_warning(evidence_heatmap_plot, viz$inference_reliable, viz$min_patients)
plot_objects$evidence_heatmap <- evidence_heatmap_plot

pdf(
  file.path(viz$integrated_dir, paste0("evidence_heatmap_top.", top_feature_show, ".pdf")),
  width = evidence_heatmap_width,
  height = evidence_heatmap_height
)
print(evidence_heatmap_plot)
dev.off()

#####
# SHAP effect versus expression log2FC
#####
top_feature_show <- 50
top_feature_label <- 50
shap_vs_expression_width <- 8
shap_vs_expression_height <- 7
shap_vs_expression_dot_label_size <- 3

shap_vs_expression_plot <- plot_shap_vs_expression(
  df = combined_results,
  top_n_to_show = top_feature_show,
  top_n_to_name = top_feature_label,
  what_0 = viz$task_what_0,
  what_1 = viz$task_what_1,
  group_colors = viz$group_colors,
  expression_log2FC_cutoff = viz$feature_min_log2FC,
  shap_expr_cor_cutoff = viz$shap_expr_cor_cutoff,
  dot_label_size = shap_vs_expression_dot_label_size
)
shap_vs_expression_plot <- add_inference_warning(shap_vs_expression_plot, viz$inference_reliable, viz$min_patients)
plot_objects$shap_vs_expression <- shap_vs_expression_plot

pdf(
  file.path(viz$integrated_dir, "shap_vs_expression.pdf"),
  width = shap_vs_expression_width,
  height = shap_vs_expression_height
)
print(shap_vs_expression_plot)
dev.off()

#####
# Value-SHAP relationship by feature
#####
top_feature_show <- 50
top_feature_label <- 50
valueSHAP_by_feature_width <- 8
valueSHAP_by_feature_height <- 8
valueSHAP_by_feature_dot_label_size <- 3

if (isTRUE(viz$inference_reliable)) {
  valueSHAP_by_feature_plot <- plot_valueSHAP_by_feature(
    df = combined_results,
    top_n_to_show = top_feature_show,
    top_n_to_name = top_feature_label,
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
    width = valueSHAP_by_feature_width,
    height = valueSHAP_by_feature_height
  )
  print(valueSHAP_by_feature_plot)
  dev.off()
}

if (isTRUE(save_plot_objects)) {
  saveRDS(plot_objects, file.path(viz$output_docs_dir, "plots.rds"))
}
