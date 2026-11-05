#' Choose an adaptive metacell size
#'
#' @param n_cells Number of rows in the current sample-by-split-by-condition
#'   group.
#' @param min_cells_to_build Minimum rows required before metacells are built.
#' @param min_cells_per_metacell Minimum target rows per metacell.
#' @param max_cells_per_metacell Maximum target rows per metacell.
#'
#' @return Integer target metacell size.
#' @export
choose_adaptive_metacell_size <- function(
    n_cells,
    min_cells_to_build = 40,
    min_cells_per_metacell = 5,
    max_cells_per_metacell = 25
) {
  if (!is.numeric(n_cells) || length(n_cells) != 1 || is.na(n_cells) || n_cells < 1) {
    stop("n_cells must be a single positive number.")
  }
  
  n_cells <- as.integer(n_cells)
  
  if (n_cells < min_cells_to_build) {
    return(1L)
  }
  
  k <- if (n_cells < 100) {
    5L
  } else if (n_cells < 500) {
    10L
  } else if (n_cells < 2000) {
    15L
  } else if (n_cells < 10000) {
    20L
  } else {
    25L
  }
  
  k <- max(min_cells_per_metacell, min(max_cells_per_metacell, k))
  
  if (floor(n_cells / k) < 2) {
    return(max(1L, floor(n_cells / 2)))
  }
  
  as.integer(k)
}


#' Compute KNN indices for metacell construction
#'
#' @param embedding Numeric embedding matrix.
#' @param k Number of nearest neighbors.
#'
#' @return List with neighbor indices and distances.
#' @export
compute_knn_indices <- function(embedding, k) {
  n <- nrow(embedding)
  
  if (n < 2) {
    return(list(
      nn.idx = matrix(integer(0), nrow = n, ncol = 0),
      nn.dists = matrix(numeric(0), nrow = n, ncol = 0)
    ))
  }
  
  k <- min(as.integer(k), n - 1L)
  
  if (k < 1) {
    return(list(
      nn.idx = matrix(integer(0), nrow = n, ncol = 0),
      nn.dists = matrix(numeric(0), nrow = n, ncol = 0)
    ))
  }
  
  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(data = embedding, query = embedding, k = k + 1L)
    return(list(
      nn.idx = nn$nn.idx[, -1, drop = FALSE],
      nn.dists = nn$nn.dists[, -1, drop = FALSE]
    ))
  }
  
  if (n > 4000) {
    stop("Please install the 'RANN' package for efficient metacell construction on large groups.")
  }
  
  dist_mat <- as.matrix(stats::dist(embedding))
  ord <- t(apply(dist_mat, 1, order))
  
  list(
    nn.idx = ord[, 2:(k + 1), drop = FALSE],
    nn.dists = t(apply(dist_mat, 1, sort))[, 2:(k + 1), drop = FALSE]
  )
}


#' Build metacells for one sample-by-split-by-condition group
#'
#' @param idx Row indices for the current group.
#' @param dataset Full analysis data set.
#' @param feature_cols Numeric feature columns.
#' @param sample_col Sample identifier column.
#' @param split_level_column Column defining the analysis subset used for within-group
#'   metacell construction. Metacells are also built within condition.
#' @param condition_col Condition column.
#' @param min_cells_to_build Minimum rows required before metacells are built.
#' @param min_cells_per_metacell Minimum target rows per metacell.
#' @param max_cells_per_metacell Maximum target rows per metacell.
#' @param n_variable_features Maximum variable features used for KNN.
#' @param n_pcs Number of PCs used for embedding.
#' @param aggregation_fun Aggregation function: `"mean"` or `"sum"`.
#' @param scale_features Whether to scale features before PCA.
#' @param verbose Whether to print progress messages.
#'
#' @return List containing metacell data and a summary row.
#' @keywords internal
.build_one_split_metacell_group <- function(
    idx,
    dataset,
    feature_cols,
    sample_col,
    split_level_column,
    condition_col,
    min_cells_to_build,
    min_cells_per_metacell,
    max_cells_per_metacell,
    n_variable_features,
    n_pcs,
    aggregation_fun,
    scale_features,
    verbose
) {
  sub <- dataset[idx, c(condition_col, sample_col, split_level_column, feature_cols), drop = FALSE]
  n_group <- nrow(sub)
  
  this_sample <- as.character(sub[[sample_col]][1])
  this_split <- as.character(sub[[split_level_column]][1])
  this_condition <- unique(as.character(sub[[condition_col]]))
  
  if (length(this_condition) != 1) {
    stop(
      paste0(
        "Group ", this_sample, " / ", this_split,
        " contains multiple conditions. Build metacells within condition as well, or fix metadata."
      )
    )
  }
  
  k_target <- choose_adaptive_metacell_size(
    n_cells = n_group,
    min_cells_to_build = min_cells_to_build,
    min_cells_per_metacell = min_cells_per_metacell,
    max_cells_per_metacell = max_cells_per_metacell
  )
  
  if (verbose) {
    message(
      "Metacell grouping: sample=", this_sample,
      " split=", this_split,
      " condition=", this_condition,
      " n=", n_group,
      " target_k=", k_target
    )
  }
  
  if (k_target <= 1L || n_group < min_cells_to_build) {
    out_sub <- sub
    out_sub$metacell_id <- paste0(this_sample, "__", this_split, "__", this_condition, "__cell_", seq_len(n_group))
    out_sub$metacell_size <- 1L
    out_sub$analysis_unit <- "single_cell_fallback"
    
    return(list(
      data = out_sub[, c(condition_col, sample_col, split_level_column, "metacell_id", "metacell_size", "analysis_unit", feature_cols), drop = FALSE],
      summary = data.frame(
        sample_id = this_sample,
        split_level = this_split,
        condition = this_condition,
        input_cells = n_group,
        chosen_k = 1L,
        output_metacells = n_group,
        stringsAsFactors = FALSE
      )
    ))
  }
  
  expr_mat <- as.matrix(sub[, feature_cols, drop = FALSE])
  storage.mode(expr_mat) <- "double"
  
  feature_var <- apply(expr_mat, 2, stats::var, na.rm = TRUE)
  feature_var[is.na(feature_var)] <- 0
  keep_var <- feature_var > 0
  expr_use <- expr_mat[, keep_var, drop = FALSE]
  
  if (ncol(expr_use) == 0) {
    out_sub <- sub
    out_sub$metacell_id <- paste0(this_sample, "__", this_split, "__", this_condition, "__cell_", seq_len(n_group))
    out_sub$metacell_size <- 1L
    out_sub$analysis_unit <- "single_cell_zero_var_fallback"
    
    return(list(
      data = out_sub[, c(condition_col, sample_col, split_level_column, "metacell_id", "metacell_size", "analysis_unit", feature_cols), drop = FALSE],
      summary = data.frame(
        sample_id = this_sample,
        split_level = this_split,
        condition = this_condition,
        input_cells = n_group,
        chosen_k = 1L,
        output_metacells = n_group,
        stringsAsFactors = FALSE
      )
    ))
  }
  
  if (ncol(expr_use) > n_variable_features) {
    top_features <- names(sort(feature_var[keep_var], decreasing = TRUE))[seq_len(n_variable_features)]
    expr_use <- expr_use[, top_features, drop = FALSE]
  }
  
  if (scale_features) {
    expr_use <- scale(expr_use)
    expr_use[is.na(expr_use)] <- 0
  }
  
  n_pcs_use <- min(n_pcs, nrow(expr_use) - 1L, ncol(expr_use))
  
  if (n_pcs_use >= 2) {
    pca_res <- stats::prcomp(expr_use, center = FALSE, scale. = FALSE, rank. = n_pcs_use)
    embedding <- pca_res$x[, seq_len(n_pcs_use), drop = FALSE]
  } else {
    embedding <- matrix(expr_use[, 1], ncol = 1)
    row.names(embedding) <- row.names(sub)
  }
  
  knn_k <- min(k_target - 1L, nrow(embedding) - 1L)
  knn <- compute_knn_indices(embedding = embedding, k = knn_k)
  
  assignments <- rep(NA_integer_, n_group)
  unassigned <- rep(TRUE, n_group)
  metacell_counter <- 0L
  
  if (ncol(knn$nn.dists) > 0) {
    seed_order <- order(rowMeans(knn$nn.dists, na.rm = TRUE))
  } else {
    seed_order <- seq_len(n_group)
  }
  
  for (seed in seed_order) {
    if (!unassigned[seed]) {
      next
    }
    
    neighbor_pool <- if (ncol(knn$nn.idx) > 0) knn$nn.idx[seed, ] else integer(0)
    neighbor_pool <- neighbor_pool[neighbor_pool > 0]
    neighbor_pool <- neighbor_pool[unassigned[neighbor_pool]]
    member_idx <- unique(c(seed, neighbor_pool))
    
    if (length(member_idx) < min_cells_per_metacell) {
      next
    }
    
    if (length(member_idx) > k_target) {
      member_idx <- member_idx[seq_len(k_target)]
    }
    
    metacell_counter <- metacell_counter + 1L
    assignments[member_idx] <- metacell_counter
    unassigned[member_idx] <- FALSE
  }
  
  leftover_idx <- which(is.na(assignments))
  
  if (metacell_counter == 0L) {
    n_centers <- max(2L, floor(n_group / k_target))
    km_res <- stats::kmeans(embedding, centers = n_centers, iter.max = 50)
    assignments <- km_res$cluster
  } else if (length(leftover_idx) > 0) {
    centroids <- do.call(
      rbind,
      lapply(seq_len(metacell_counter), function(mc_id) {
        members <- which(assignments == mc_id)
        if (length(members) == 0) {
          return(rep(NA_real_, ncol(embedding)))
        }
        colMeans(embedding[members, , drop = FALSE], na.rm = TRUE)
      })
    )
    
    centroids <- as.matrix(centroids)
    if (nrow(centroids) == 0 || ncol(centroids) == 0) {
      stop("Failed to compute metacell centroids for leftover assignment.")
    }
    
    for (j in leftover_idx) {
      cell_vec <- as.numeric(embedding[j, , drop = TRUE])
      d <- apply(centroids, 1, function(center_vec) {
        center_vec <- as.numeric(center_vec)
        if (length(center_vec) != length(cell_vec)) {
          return(Inf)
        }
        if (any(!is.finite(center_vec)) || any(!is.finite(cell_vec))) {
          return(Inf)
        }
        sum((center_vec - cell_vec)^2)
      })
      
      assignments[j] <- if (all(!is.finite(d))) 1L else which.min(d)
    }
  }
  
  metacell_sizes <- table(assignments)
  agg_mat <- rowsum(expr_mat, group = assignments, reorder = TRUE)
  
  if (aggregation_fun == "mean") {
    agg_mat <- agg_mat / as.numeric(metacell_sizes[rownames(agg_mat)])
  }
  
  agg_df <- as.data.frame(agg_mat, check.names = FALSE)
  agg_df$metacell_id <- paste0(
    this_sample, "__", this_split, "__", this_condition, "__MC_",
    sprintf("%03d", seq_len(nrow(agg_df)))
  )
  agg_df$metacell_size <- as.integer(metacell_sizes[rownames(agg_df)])
  agg_df$analysis_unit <- "metacell"
  agg_df[[condition_col]] <- this_condition
  agg_df[[sample_col]] <- this_sample
  agg_df[[split_level_column]] <- this_split
  
  row.names(agg_df) <- agg_df$metacell_id
  
  list(
    data = agg_df[, c(condition_col, sample_col, split_level_column, "metacell_id", "metacell_size", "analysis_unit", feature_cols), drop = FALSE],
    summary = data.frame(
      sample_id = this_sample,
      split_level = this_split,
      condition = this_condition,
      input_cells = n_group,
      chosen_k = k_target,
      output_metacells = nrow(agg_df),
      stringsAsFactors = FALSE
    )
  )
}


#' Build grouped metacells without hdWGCNA
#'
#' @param dataset Input data frame containing condition, sample, split, and
#'   feature columns.
#' @param feature_cols Feature columns to aggregate.
#' @param sample_col Sample identifier column.
#' @param split_level_column Column used to keep metacells within an analysis subset.
#'   Metacells are built within sample, split, and condition to avoid mixing
#'   biological groups before modeling.
#' @param condition_col Condition column.
#' @param use_metacells Whether to build metacells or keep single rows.
#' @param min_cells_to_build Minimum rows required before metacell building.
#' @param min_cells_per_metacell Minimum target rows per metacell.
#' @param max_cells_per_metacell Maximum target rows per metacell.
#' @param n_variable_features Number of variable features used for embedding.
#' @param n_pcs Number of principal components for embedding.
#' @param aggregation_fun Aggregation function: `"mean"` or `"sum"`.
#' @param scale_features Whether to scale features before PCA.
#' @param seed Random seed.
#' @param verbose Whether to print progress.
#' @param parallel Whether to build groups in parallel.
#' @param n_cores Number of worker cores.
#' @param parallel_backend One of `"serial"`, `"mclapply"`, or `"psock"`.
#'
#' @return List with aggregated data and a summary table.
#' @export
build_grouped_metacells <- function(
    dataset,
    feature_cols,
    sample_col = "sample_id",
    split_level_column,
    condition_col = "condition",
    use_metacells = TRUE,
    min_cells_to_build = 40,
    min_cells_per_metacell = 5,
    max_cells_per_metacell = 25,
    n_variable_features = 2000,
    n_pcs = 20,
    aggregation_fun = c("mean", "sum"),
    scale_features = TRUE,
    seed = 42,
    verbose = TRUE,
    parallel = FALSE,
    n_cores = 2,
    parallel_backend = c("serial", "mclapply", "psock")
) {
  aggregation_fun <- match.arg(aggregation_fun)
  parallel_backend <- match.arg(parallel_backend)
  
  required_cols <- c(condition_col, sample_col, split_level_column, feature_cols)
  missing_cols <- setdiff(required_cols, colnames(dataset))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  if (!use_metacells) {
    out <- dataset[, required_cols, drop = FALSE]
    out$metacell_id <- row.names(out)
    out$metacell_size <- 1L
    out$analysis_unit <- "single_cell"
    
    summary_table <- data.frame(
      sample_id = out[[sample_col]],
      split_level = out[[split_level_column]],
      condition = out[[condition_col]],
      input_cells = 1L,
      chosen_k = 1L,
      output_metacells = 1L,
      stringsAsFactors = FALSE
    )
    
    return(list(
      data = out[, c(condition_col, sample_col, split_level_column, "metacell_id", "metacell_size", "analysis_unit", feature_cols), drop = FALSE],
      summary = summary_table
    ))
  }
  
  group_key <- interaction(
    dataset[[sample_col]],
    dataset[[split_level_column]],
    dataset[[condition_col]],
    drop = TRUE,
    lex.order = TRUE
  )
  
  split_idx <- split(seq_len(nrow(dataset)), group_key)
  group_seeds <- vapply(
    names(split_idx),
    function(group_name) derive_shapture_seed(seed, "metacell_group", group_name),
    integer(1)
  )
  work_items <- seq_along(split_idx)
  
  worker_fun <- function(group_position) {
    set.seed(group_seeds[[group_position]])
    idx <- split_idx[[group_position]]
    .build_one_split_metacell_group(
      idx = idx,
      dataset = dataset,
      feature_cols = feature_cols,
      sample_col = sample_col,
      split_level_column = split_level_column,
      condition_col = condition_col,
      min_cells_to_build = min_cells_to_build,
      min_cells_per_metacell = min_cells_per_metacell,
      max_cells_per_metacell = max_cells_per_metacell,
      n_variable_features = n_variable_features,
      n_pcs = n_pcs,
      aggregation_fun = aggregation_fun,
      scale_features = scale_features,
      verbose = verbose
    )
  }
  
  if (!parallel || length(split_idx) <= 1 || n_cores <= 1) {
    res_list <- lapply(work_items, worker_fun)
  } else if (parallel_backend == "mclapply" && .Platform$OS.type != "windows") {
    res_list <- parallel::mclapply(work_items, worker_fun, mc.cores = n_cores)
  } else {
    cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterExport(
      cl,
      varlist = c(
        ".build_one_split_metacell_group",
        "choose_adaptive_metacell_size",
        "compute_knn_indices",
        "dataset",
        "feature_cols",
        "sample_col",
        "split_level_column",
        "condition_col",
        "min_cells_to_build",
        "min_cells_per_metacell",
        "max_cells_per_metacell",
        "n_variable_features",
        "n_pcs",
        "aggregation_fun",
        "scale_features",
        "verbose",
        "split_idx",
        "group_seeds"
      ),
      envir = environment()
    )
    
    res_list <- parallel::parLapply(cl, work_items, worker_fun)
  }
  
  list(
    data = do.call(rbind, lapply(res_list, `[[`, "data")),
    summary = do.call(rbind, lapply(res_list, `[[`, "summary"))
  )
}
