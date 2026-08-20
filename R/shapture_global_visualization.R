#' Discover completed SHAPTURE result folders
#'
#' @param results_root Folder containing multiple completed SHAPTURE result
#'   folders. Discovery searches recursively for `docs/combined_results.csv`.
#'
#' @return Normalized paths to completed result folders.
#' @export
discover_shapture_result_dirs <- function(results_root) {
  if (length(results_root) != 1L || is.na(results_root) || !dir.exists(results_root)) {
    stop("results_root must be one existing directory.")
  }

  results_root <- normalizePath(results_root, winslash = "/", mustWork = TRUE)
  result_files <- list.files(
    results_root,
    pattern = "^combined_results[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  result_files <- result_files[basename(dirname(result_files)) == "docs"]
  result_dirs <- unique(dirname(dirname(result_files)))

  if (length(result_dirs) < 2L) {
    stop(
      "At least two completed SHAPTURE result folders are required under results_root."
    )
  }

  normalizePath(result_dirs, winslash = "/", mustWork = TRUE)
}


#' Normalize a comparison mode for global visualization
#'
#' @param comparison_mode Comparison-mode label saved by SHAPTURE.
#'
#' @return `condition_vs_condition` or `one_vs_others`.
#' @keywords internal
normalize_global_comparison_mode <- function(comparison_mode) {
  comparison_mode <- tolower(trimws(as.character(comparison_mode)[1]))
  if (comparison_mode %in% c("condition_vs_condition", "condition-vs-condition")) {
    return("condition_vs_condition")
  }
  if (comparison_mode %in% c(
    "one_vs_others", "one-vs-others", "one_vs_rest", "one-vs-rest"
  )) {
    return("one_vs_others")
  }
  stop("Unsupported comparison_mode in completed results: ", comparison_mode)
}


#' Read task context needed for global visualization
#'
#' @param result_dir Completed SHAPTURE result folder.
#'
#' @return A task-context list.
#' @keywords internal
read_shapture_global_context <- function(result_dir) {
  visualization_path <- file.path(result_dir, "docs", "visualization_inputs.rds")
  summary_path <- file.path(
    result_dir,
    "Integrated_adaptive_scoring_results",
    "run_summary.csv"
  )

  if (file.exists(visualization_path)) {
    inputs <- readRDS(visualization_path)
    context <- inputs$context
    required <- c(
      "temp_split", "comparison_mode", "task_what_0", "task_what_1"
    )
    if (is.list(context) && all(required %in% names(context))) {
      return(list(
        split_level = as.character(context$temp_split)[1],
        comparison_mode = normalize_global_comparison_mode(context$comparison_mode),
        condition_0 = as.character(context$task_what_0)[1],
        condition_1 = as.character(context$task_what_1)[1]
      ))
    }
  }

  if (file.exists(summary_path)) {
    summary <- utils::read.csv(
      summary_path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    required <- c(
      "split_level", "comparison_mode", "condition_0", "condition_1"
    )
    if (nrow(summary) >= 1L && all(required %in% names(summary))) {
      return(list(
        split_level = as.character(summary$split_level)[1],
        comparison_mode = normalize_global_comparison_mode(summary$comparison_mode),
        condition_0 = as.character(summary$condition_0)[1],
        condition_1 = as.character(summary$condition_1)[1]
      ))
    }
  }

  stop(
    "Could not recover split level, comparison mode, and condition labels from ",
    result_dir,
    ". Keep docs/visualization_inputs.rds or the task run_summary.csv."
  )
}


#' Normalize a saved model-expression pattern
#'
#' @param pattern Saved `Model_expression_pattern` value.
#' @param condition_0 Negative/reference label for this task.
#' @param condition_1 Positive/target label for this task.
#'
#' @return `condition_0`, `condition_1`, or `NA_character_`.
#' @keywords internal
normalize_global_direction <- function(pattern, condition_0, condition_1) {
  pattern <- trimws(gsub("[\r\n]+", " ", as.character(pattern)))
  label_0 <- paste0("Higher feature values push toward ", condition_0)
  label_1 <- paste0("Higher feature values push toward ", condition_1)

  direction <- rep(NA_character_, length(pattern))
  direction[pattern == label_0] <- "condition_0"
  direction[pattern == label_1] <- "condition_1"
  direction
}


#' Format a saved label for figure display
#'
#' @param x Character labels.
#'
#' @return Labels with underscores replaced by spaces.
#' @keywords internal
format_shapture_global_label <- function(x) {
  gsub("_", " ", as.character(x), fixed = TRUE)
}


#' Validate a continuous color palette
#'
#' @param continuous_colors R color vector.
#'
#' @return Unnamed validated colors.
#' @keywords internal
validate_global_palette <- function(continuous_colors) {
  continuous_colors <- unname(as.character(continuous_colors))
  if (length(continuous_colors) < 2L || anyNA(continuous_colors)) {
    stop("continuous_colors must contain at least two valid R colors.")
  }
  tryCatch(
    grDevices::col2rgb(continuous_colors),
    error = function(e) stop("continuous_colors contains an invalid R color.")
  )
  continuous_colors
}


#' Load multiple completed SHAPTURE comparisons for global visualization
#'
#' Each result folder contributes its top `top_feature_show` ranked features.
#' The union of those per-folder lists defines the global feature landscape and
#' direction-aware Jaccard comparison.
#'
#' @param results_root Folder containing multiple completed result folders.
#' @param top_feature_show Number of top-ranked features selected independently
#'   from each result folder.
#' @param result_dirs Optional explicit result-folder paths. When `NULL`, they
#'   are discovered under `results_root`.
#'
#' @return A `shapture_global_results` object containing selected rows, score
#'   and direction matrices, and a direction-aware Jaccard matrix.
#' @export
load_shapture_global_results <- function(
    results_root,
    top_feature_show = 50L,
    result_dirs = NULL
) {
  if (
    length(top_feature_show) != 1L ||
      is.na(top_feature_show) ||
      !is.finite(top_feature_show) ||
      top_feature_show < 1
  ) {
    stop("top_feature_show must be one positive finite number.")
  }
  top_feature_show <- as.integer(round(top_feature_show))

  if (is.null(result_dirs)) {
    result_dirs <- discover_shapture_result_dirs(results_root)
  } else {
    missing_dirs <- result_dirs[!dir.exists(result_dirs)]
    if (length(missing_dirs) > 0L) {
      stop("These result_dirs do not exist: ", paste(missing_dirs, collapse = ", "))
    }
    result_dirs <- unique(normalizePath(
      result_dirs,
      winslash = "/",
      mustWork = TRUE
    ))
    if (length(result_dirs) < 2L) {
      stop("At least two result_dirs are required.")
    }
  }

  tasks <- lapply(result_dirs, function(result_dir) {
    result_path <- file.path(result_dir, "docs", "combined_results.csv")
    if (!file.exists(result_path)) {
      stop("Missing docs/combined_results.csv in ", result_dir)
    }

    context <- read_shapture_global_context(result_dir)
    result <- utils::read.csv(
      result_path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    required <- c("feature", "final_score", "Model_expression_pattern")
    missing_cols <- setdiff(required, names(result))
    if (length(missing_cols) > 0L) {
      stop(
        "Missing required columns in ", result_path, ": ",
        paste(missing_cols, collapse = ", ")
      )
    }
    if (nrow(result) < 1L) {
      stop("No retained features were found in ", result_path)
    }

    result$feature <- as.character(result$feature)
    result$final_score <- suppressWarnings(as.numeric(result$final_score))
    if (any(!is.finite(result$final_score))) {
      stop("final_score must be finite in ", result_path)
    }

    if ("rank" %in% names(result)) {
      result$rank <- suppressWarnings(as.numeric(result$rank))
      rank_order <- order(
        ifelse(is.finite(result$rank), result$rank, Inf),
        -result$final_score,
        result$feature
      )
    } else {
      rank_order <- order(-result$final_score, result$feature)
    }
    result <- result[rank_order, , drop = FALSE]
    result <- result[!duplicated(result$feature), , drop = FALSE]
    result <- utils::head(result, min(top_feature_show, nrow(result)))

    result$direction_code <- normalize_global_direction(
      pattern = result$Model_expression_pattern,
      condition_0 = context$condition_0,
      condition_1 = context$condition_1
    )
    result$split_level <- context$split_level
    result$result_dir <- result_dir

    list(context = context, selected = result)
  })

  split_levels <- vapply(tasks, function(task) task$context$split_level, character(1))
  if (anyNA(split_levels) || any(!nzchar(split_levels))) {
    stop("Every result folder must have a non-empty split-level label.")
  }
  if (anyDuplicated(split_levels)) {
    duplicated_levels <- unique(split_levels[duplicated(split_levels)])
    stop(
      "Duplicate split-level labels were found: ",
      paste(duplicated_levels, collapse = ", "),
      ". Use a narrower results_root or explicit result_dirs."
    )
  }

  comparison_modes <- unique(vapply(
    tasks,
    function(task) task$context$comparison_mode,
    character(1)
  ))
  if (length(comparison_modes) != 1L) {
    stop("Do not mix condition-vs-condition and one-vs-others folders.")
  }
  comparison_mode <- comparison_modes[[1]]

  if (comparison_mode == "condition_vs_condition") {
    condition_0 <- unique(vapply(
      tasks,
      function(task) task$context$condition_0,
      character(1)
    ))
    condition_1 <- unique(vapply(
      tasks,
      function(task) task$context$condition_1,
      character(1)
    ))
    if (length(condition_0) != 1L || length(condition_1) != 1L) {
      stop(
        "Condition-vs-condition folders must use the same condition labels ",
        "for direction-aware global visualization."
      )
    }
    direction_labels <- c(
      condition_1 = format_shapture_global_label(condition_1[[1]]),
      condition_0 = format_shapture_global_label(condition_0[[1]])
    )
    comparison_label <- paste(
      format_shapture_global_label(condition_0[[1]]),
      "versus",
      format_shapture_global_label(condition_1[[1]])
    )
  } else {
    direction_labels <- c(condition_1 = "Target", condition_0 = "Others")
    comparison_label <- "Target versus others"
  }

  task_order <- order(split_levels)
  tasks <- tasks[task_order]
  split_levels <- split_levels[task_order]
  selected_rows <- do.call(rbind, lapply(tasks, `[[`, "selected"))
  rownames(selected_rows) <- NULL

  unknown_direction <- sum(is.na(selected_rows$direction_code))
  if (unknown_direction > 0L) {
    warning(
      unknown_direction,
      " selected feature records had an unrecognized Model_expression_pattern ",
      "and cannot contribute a directional match."
    )
  }

  feature_summary <- stats::aggregate(
    final_score ~ feature,
    data = selected_rows,
    FUN = mean
  )
  feature_breadth <- table(selected_rows$feature)
  feature_summary$n_split_levels <- as.integer(feature_breadth[feature_summary$feature])
  feature_summary <- feature_summary[
    order(-feature_summary$n_split_levels, -feature_summary$final_score, feature_summary$feature),
    ,
    drop = FALSE
  ]
  features <- feature_summary$feature

  score_matrix <- matrix(
    NA_real_,
    nrow = length(features),
    ncol = length(split_levels),
    dimnames = list(features, split_levels)
  )
  direction_matrix <- matrix(
    NA_character_,
    nrow = length(features),
    ncol = length(split_levels),
    dimnames = list(features, split_levels)
  )

  for (i in seq_len(nrow(selected_rows))) {
    feature <- selected_rows$feature[[i]]
    split_level <- selected_rows$split_level[[i]]
    score_matrix[feature, split_level] <- selected_rows$final_score[[i]]
    direction_matrix[feature, split_level] <- selected_rows$direction_code[[i]]
  }

  directional_jaccard <- matrix(
    0,
    nrow = length(split_levels),
    ncol = length(split_levels),
    dimnames = list(split_levels, split_levels)
  )

  selected_by_level <- setNames(lapply(split_levels, function(split_level) {
    selected_rows[selected_rows$split_level == split_level, , drop = FALSE]
  }), split_levels)

  for (i in seq_along(split_levels)) {
    for (j in seq_along(split_levels)) {
      if (i == j) {
        directional_jaccard[i, j] <- 1
        next
      }

      task_i <- selected_by_level[[i]]
      task_j <- selected_by_level[[j]]
      union_features <- union(task_i$feature, task_j$feature)
      shared_features <- intersect(task_i$feature, task_j$feature)

      direction_i <- stats::setNames(task_i$direction_code, task_i$feature)
      direction_j <- stats::setNames(task_j$direction_code, task_j$feature)
      same_direction <- !is.na(direction_i[shared_features]) &
        !is.na(direction_j[shared_features]) &
        direction_i[shared_features] == direction_j[shared_features]

      directional_jaccard[i, j] <- if (length(union_features) == 0L) {
        NA_real_
      } else {
        sum(same_direction) / length(union_features)
      }
    }
  }

  result <- list(
    results_root = normalizePath(results_root, winslash = "/", mustWork = TRUE),
    result_dirs = result_dirs[task_order],
    comparison_mode = comparison_mode,
    comparison_label = comparison_label,
    direction_labels = direction_labels,
    top_feature_show = top_feature_show,
    split_levels = split_levels,
    selected_rows = selected_rows,
    feature_summary = feature_summary,
    score_matrix = score_matrix,
    direction_matrix = direction_matrix,
    directional_jaccard = directional_jaccard,
    n_selected_by_level = vapply(selected_by_level, nrow, integer(1))
  )
  class(result) <- c("shapture_global_results", "list")
  result
}


#' Build a clustered cross-level SHAPTURE feature landscape
#'
#' @param global_results Output from [load_shapture_global_results()].
#' @param continuous_colors Continuous color palette for final scores.
#' @param clustering_method Linkage method passed to [stats::hclust()].
#' @param show_row_names Whether feature names are displayed.
#' @param show_direction_symbols Whether directional triangles are displayed.
#' @param direction_symbol_size_mm Direction-symbol size in millimeters.
#' @param tile_color Tile-border color; use `NA` for no border.
#' @param tile_linewidth Tile-border width.
#' @param use_raster Whether to rasterize the heatmap body. `NULL` enables it
#'   automatically for more than 20,000 cells while preserving vector labels.
#' @param title Optional figure title.
#'
#' @return A list containing a ComplexHeatmap object, direction legend, and the
#'   plotted matrices.
#' @export
build_shapture_cross_level_landscape <- function(
    global_results,
    continuous_colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
    clustering_method = "average",
    show_row_names = TRUE,
    show_direction_symbols = TRUE,
    direction_symbol_size_mm = 2.5,
    tile_color = "white",
    tile_linewidth = 0.7,
    use_raster = NULL,
    title = NULL
) {
  if (!inherits(global_results, "shapture_global_results")) {
    stop("global_results must come from load_shapture_global_results().")
  }
  continuous_colors <- validate_global_palette(continuous_colors)
  score_matrix <- global_results$score_matrix
  direction_matrix <- global_results$direction_matrix
  clustering_matrix <- score_matrix
  clustering_matrix[is.na(clustering_matrix)] <- 0

  if (is.null(use_raster)) {
    use_raster <- length(score_matrix) > 20000L
  }
  if (length(use_raster) != 1L || is.na(use_raster)) {
    stop("use_raster must be TRUE, FALSE, or NULL.")
  }

  row_clustering <- if (nrow(clustering_matrix) > 1L) {
    stats::hclust(
      stats::dist(clustering_matrix),
      method = clustering_method
    )
  } else {
    FALSE
  }
  column_clustering <- if (ncol(clustering_matrix) > 1L) {
    stats::hclust(
      stats::dist(t(clustering_matrix)),
      method = clustering_method
    )
  } else {
    FALSE
  }

  score_colors <- circlize::colorRamp2(
    seq(0, 1, length.out = length(continuous_colors)),
    continuous_colors
  )
  if (is.null(title)) {
    title <- paste0(
      "Cross-level feature landscape\n",
      global_results$comparison_label
    )
  }

  layer_fun <- if (isTRUE(show_direction_symbols)) {
    function(j, i, x, y, width, height, fill) {
      matrix_index <- cbind(i, j)
      direction <- direction_matrix[matrix_index]
      score <- score_matrix[matrix_index]
      valid <- !is.na(direction) & is.finite(score)
      if (!any(valid)) {
        return(invisible(NULL))
      }
      symbol_color <- ifelse(
        score[valid] < 0.23 | score[valid] > 0.74,
        "white",
        "#263437"
      )
      grid::grid.points(
        x = x[valid],
        y = y[valid],
        pch = ifelse(direction[valid] == "condition_1", 24, 25),
        size = grid::unit(direction_symbol_size_mm, "mm"),
        gp = grid::gpar(
          col = symbol_color,
          fill = symbol_color,
          lwd = 0.5
        )
      )
      invisible(NULL)
    }
  } else {
    NULL
  }

  heatmap <- ComplexHeatmap::Heatmap(
    score_matrix,
    name = "Final score",
    col = score_colors,
    na_col = "#E3E4DF",
    cluster_rows = row_clustering,
    cluster_columns = column_clustering,
    show_row_dend = !identical(row_clustering, FALSE),
    show_column_dend = !identical(column_clustering, FALSE),
    row_dend_width = grid::unit(21, "mm"),
    column_dend_height = grid::unit(18, "mm"),
    show_row_names = show_row_names,
    row_names_gp = grid::gpar(
      fontsize = if (nrow(score_matrix) > 100L) 6 else 10,
      fontface = "bold",
      col = "#263437"
    ),
    column_names_gp = grid::gpar(fontsize = 9, col = "#425154"),
    column_names_rot = 58,
    rect_gp = grid::gpar(col = tile_color, lwd = tile_linewidth),
    column_title = title,
    column_title_gp = grid::gpar(
      fontsize = 17,
      fontface = "bold",
      col = "#172326"
    ),
    heatmap_legend_param = list(
      at = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.50", "0.75", "1"),
      title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 9),
      legend_height = grid::unit(43, "mm")
    ),
    use_raster = isTRUE(use_raster),
    raster_by_magick = FALSE,
    raster_quality = 2,
    layer_fun = layer_fun
  )

  direction_legend <- if (isTRUE(show_direction_symbols)) {
    ComplexHeatmap::Legend(
      title = "Higher feature\nvalues push toward",
      labels = unname(global_results$direction_labels),
      type = "points",
      pch = c(24, 25),
      legend_gp = grid::gpar(col = "#263437", fill = "#263437"),
      size = grid::unit(3, "mm"),
      title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 9)
    )
  } else {
    NULL
  }

  structure(
    list(
      heatmap = heatmap,
      direction_legend = direction_legend,
      score_matrix = score_matrix,
      direction_matrix = direction_matrix,
      clustering_method = clustering_method,
      use_raster = isTRUE(use_raster)
    ),
    class = c("shapture_global_landscape", "list")
  )
}


#' Build a direction-aware Jaccard similarity heatmap
#'
#' @param global_results Output from [load_shapture_global_results()].
#' @param continuous_colors Continuous color palette.
#' @param clustering_method Linkage method passed to [stats::hclust()].
#' @param show_values Whether similarity values are printed in heatmap cells.
#' @param value_digits Number of digits printed for similarity values.
#' @param tile_color Tile-border color; use `NA` for no border.
#' @param tile_linewidth Tile-border width.
#' @param use_raster Whether to rasterize the heatmap body. `NULL` enables it
#'   automatically for more than 20,000 cells while preserving vector labels.
#' @param title Optional figure title.
#'
#' @return A list containing a ComplexHeatmap object and similarity matrix.
#' @export
build_shapture_directional_jaccard_heatmap <- function(
    global_results,
    continuous_colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
    clustering_method = "average",
    show_values = TRUE,
    value_digits = 2L,
    tile_color = "white",
    tile_linewidth = 0.8,
    use_raster = NULL,
    title = "Split-level direction-aware similarity"
) {
  if (!inherits(global_results, "shapture_global_results")) {
    stop("global_results must come from load_shapture_global_results().")
  }
  continuous_colors <- validate_global_palette(continuous_colors)
  similarity_matrix <- global_results$directional_jaccard
  if (any(!is.finite(similarity_matrix))) {
    stop("The direction-aware Jaccard matrix contains non-finite values.")
  }

  if (is.null(use_raster)) {
    use_raster <- length(similarity_matrix) > 20000L
  }
  if (length(use_raster) != 1L || is.na(use_raster)) {
    stop("use_raster must be TRUE, FALSE, or NULL.")
  }
  if (
    length(value_digits) != 1L ||
      is.na(value_digits) ||
      value_digits < 0 ||
      value_digits > 6
  ) {
    stop("value_digits must be one integer between 0 and 6.")
  }
  value_digits <- as.integer(value_digits)

  similarity_clustering <- if (nrow(similarity_matrix) > 1L) {
    stats::hclust(
      stats::as.dist(1 - similarity_matrix),
      method = clustering_method
    )
  } else {
    FALSE
  }
  score_colors <- circlize::colorRamp2(
    seq(0, 1, length.out = length(continuous_colors)),
    continuous_colors
  )

  layer_fun <- if (isTRUE(show_values)) {
    function(j, i, x, y, width, height, fill) {
      values <- similarity_matrix[cbind(i, j)]
      label_color <- ifelse(values < 0.16 | values > 0.74, "white", "#263437")
      grid::grid.text(
        sprintf(paste0("%.", value_digits, "f"), values),
        x = x,
        y = y,
        gp = grid::gpar(
          fontsize = if (nrow(similarity_matrix) > 20L) 6 else 9,
          fontface = "bold",
          col = label_color
        )
      )
      invisible(NULL)
    }
  } else {
    NULL
  }

  heatmap <- ComplexHeatmap::Heatmap(
    similarity_matrix,
    name = "Direction-aware\nJaccard similarity",
    col = score_colors,
    cluster_rows = similarity_clustering,
    cluster_columns = similarity_clustering,
    show_row_dend = !identical(similarity_clustering, FALSE),
    show_column_dend = !identical(similarity_clustering, FALSE),
    row_dend_width = grid::unit(21, "mm"),
    column_dend_height = grid::unit(18, "mm"),
    row_names_gp = grid::gpar(fontsize = 9.5, col = "#344346"),
    column_names_gp = grid::gpar(fontsize = 9, col = "#425154"),
    column_names_rot = 58,
    rect_gp = grid::gpar(col = tile_color, lwd = tile_linewidth),
    column_title = title,
    column_title_gp = grid::gpar(
      fontsize = 17,
      fontface = "bold",
      col = "#172326"
    ),
    heatmap_legend_param = list(
      at = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.50", "0.75", "1"),
      title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 9),
      legend_height = grid::unit(43, "mm")
    ),
    use_raster = isTRUE(use_raster),
    raster_by_magick = FALSE,
    raster_quality = 2,
    layer_fun = layer_fun
  )

  structure(
    list(
      heatmap = heatmap,
      similarity_matrix = similarity_matrix,
      clustering_method = clustering_method,
      use_raster = isTRUE(use_raster)
    ),
    class = c("shapture_directional_jaccard_heatmap", "list")
  )
}
