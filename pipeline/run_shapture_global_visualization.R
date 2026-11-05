#####
# Editable SHAPTURE global visualization pipeline
#
# This script does not rerun models. It reads multiple completed SHAPTURE
# result folders and connects their ranked feature tables in two global views.
#####

library(SHAPTURE)
library(RColorBrewer)

#####
# Select the folder containing multiple completed SHAPTURE result folders.
# Each completed folder must contain docs/combined_results.csv and either
# docs/visualization_inputs.rds or Integrated_adaptive_scoring_results/run_summary.csv.
#####
results_root <- "path/to/folder_containing_multiple_shapture_result_folders"

# Optional explicit completed-result folders. NULL discovers them recursively
# under results_root.
result_dirs <- NULL

# Select this many top-ranked features independently from every result folder.
# Their union is displayed in the landscape and compared in the similarity map.
top_feature_show <- 50

# Any valid R color vector with at least 2 colors.
continuous_colors <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))

# Both figures use average-linkage hierarchical clustering.
clustering_method <- "average"

# Figures are saved here. Change this if you want a separate output folder.
output_dir <- results_root
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#####
# Load and combine the completed result folders.
#####
global_results <- load_shapture_global_results(
  results_root = results_root,
  top_feature_show = top_feature_show,
  result_dirs = result_dirs
)

message(
  "Loaded ", length(global_results$split_levels_to_analyze), " split levels and ",
  nrow(global_results$score_matrix), " unique selected features."
)

#####
# Cross-level feature landscape
#####
# NULL sizes the figure automatically from the actual feature and split-level
# counts. Replace either NULL with an exact size in inches to override it.
cross_level_landscape_width <- NULL
cross_level_landscape_height <- NULL

# Set FALSE when hundreds or thousands of feature names would be unreadable.
cross_level_landscape_show_feature_names <- TRUE

# Set FALSE to hide the target/reference directional triangles.
cross_level_landscape_show_direction_symbols <- TRUE
cross_level_landscape_direction_symbol_size_mm <- 2.5

# Use NA to remove tile borders.
cross_level_landscape_tile_color <- "white"
cross_level_landscape_tile_linewidth <- 0.7

# NULL automatically rasterizes only a very large heatmap body while keeping
# labels and dendrograms as vector elements. Use FALSE to force all-vector PDF.
cross_level_landscape_use_raster <- NULL

cross_level_feature_landscape <- build_shapture_cross_level_landscape(
  global_results = global_results,
  continuous_colors = continuous_colors,
  clustering_method = clustering_method,
  show_row_names = cross_level_landscape_show_feature_names,
  show_direction_symbols = cross_level_landscape_show_direction_symbols,
  direction_symbol_size_mm = cross_level_landscape_direction_symbol_size_mm,
  tile_color = cross_level_landscape_tile_color,
  tile_linewidth = cross_level_landscape_tile_linewidth,
  use_raster = cross_level_landscape_use_raster
)

cross_level_landscape_size <- calculate_shapture_global_figure_size(
  figure_type = "cross_level_landscape",
  n_features = nrow(global_results$score_matrix),
  n_split_levels = ncol(global_results$score_matrix),
  width = cross_level_landscape_width,
  height = cross_level_landscape_height
)

cross_level_landscape_path <- file.path(
  output_dir,
  paste0(
    "global_cross_level_feature_landscape_top.",
    top_feature_show,
    "_clustered.pdf"
  )
)

pdf(
  cross_level_landscape_path,
  width = cross_level_landscape_size$width,
  height = cross_level_landscape_size$height
)
ComplexHeatmap::draw(
  cross_level_feature_landscape$heatmap,
  annotation_legend_list = if (
    is.null(cross_level_feature_landscape$direction_legend)
  ) {
    list()
  } else {
    list(cross_level_feature_landscape$direction_legend)
  },
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = FALSE,
  padding = grid::unit(c(10, 15, 15, 10), "mm")
)
grid::grid.text(
  paste0(
    "Gray = feature not selected among the top ",
    top_feature_show,
    " features in that split level."
  ),
  x = grid::unit(10, "mm"),
  y = grid::unit(4, "mm"),
  just = c("left", "bottom"),
  gp = grid::gpar(fontsize = 8.5, col = "#667174")
)
dev.off()

#####
# Direction-aware Jaccard similarity
#####
# NULL sizes the figure automatically from the actual split-level count.
# Replace either NULL with an exact size in inches to override it.
directional_jaccard_width <- NULL
directional_jaccard_height <- NULL

# Value labels are useful for a modest number of split levels. Set FALSE when
# many levels make the matrix cells too small.
directional_jaccard_show_values <- TRUE
directional_jaccard_value_digits <- 2

# Use NA to remove tile borders.
directional_jaccard_tile_color <- "white"
directional_jaccard_tile_linewidth <- 0.8
directional_jaccard_use_raster <- NULL

directional_jaccard_similarity <- build_shapture_directional_jaccard_heatmap(
  global_results = global_results,
  continuous_colors = continuous_colors,
  clustering_method = clustering_method,
  show_values = directional_jaccard_show_values,
  value_digits = directional_jaccard_value_digits,
  tile_color = directional_jaccard_tile_color,
  tile_linewidth = directional_jaccard_tile_linewidth,
  use_raster = directional_jaccard_use_raster
)

directional_jaccard_size <- calculate_shapture_global_figure_size(
  figure_type = "directional_jaccard",
  n_split_levels = nrow(global_results$directional_jaccard),
  width = directional_jaccard_width,
  height = directional_jaccard_height
)

directional_jaccard_path <- file.path(
  output_dir,
  paste0(
    "global_direction_aware_jaccard_top.",
    top_feature_show,
    "_clustered.pdf"
  )
)

pdf(
  directional_jaccard_path,
  width = directional_jaccard_size$width,
  height = directional_jaccard_size$height
)
ComplexHeatmap::draw(
  directional_jaccard_similarity$heatmap,
  heatmap_legend_side = "right",
  padding = grid::unit(c(10, 14, 14, 10), "mm")
)
grid::grid.text(
  paste0(
    "Direction-aware Jaccard similarity compares the retained feature sets; ",
    "1 = identical and 0 = no overlap"
  ),
  x = grid::unit(10, "mm"),
  y = grid::unit(4, "mm"),
  just = c("left", "bottom"),
  gp = grid::gpar(fontsize = 8.5, col = "#667174")
)
dev.off()

message("Saved: ", cross_level_landscape_path)
message("Saved: ", directional_jaccard_path)
message(
  "Cross-level landscape size: ", cross_level_landscape_size$width, " x ",
  cross_level_landscape_size$height, " inches."
)
message(
  "Direction-aware Jaccard size: ", directional_jaccard_size$width, " x ",
  directional_jaccard_size$height, " inches."
)
