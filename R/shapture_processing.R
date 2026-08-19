#' Prepare balanced train/test splits for SHAPTURE
#'
#' Creates balanced binary train/test splits while preserving optional metadata
#' columns and ensuring that any scaling is learned from the training data only.
#'
#' @param dataset A data frame containing the condition column, optional metadata
#'   columns, and numeric feature columns.
#' @param temp_random_seed Integer seed used for reproducible sampling.
#' @param train_ratio Fraction of sampled rows per condition assigned to the
#'   training split.
#' @param test_ratio Fraction of sampled rows per condition assigned to the test
#'   split.
#' @param each_condition_total_sample_no Target number of rows per condition to
#'   sample before train/test partitioning.
#' @param condition_col Name of the binary condition column.
#' @param preprocess_method Either `"none"` or `"scale"`.
#' @param metadata_cols Optional metadata columns to carry through the split but
#'   exclude from the feature matrix.
#' @param sampling_prob_col Optional column of per-row sampling probabilities.
#' @param precomputed_split Optional list with `unit_col`, `train_unit_ids`,
#'   and `test_unit_ids`. When supplied, these rows are used directly instead
#'   of drawing a new random split.
#'
#' @return A list containing processed and raw train/test data plus feature
#'   metadata.
#' @export
prepare_train_test_data <- function(
    dataset,
    temp_random_seed,
    train_ratio,
    test_ratio,
    each_condition_total_sample_no,
    condition_col = "condition",
    preprocess_method = c("none", "scale"),
    metadata_cols = NULL,
    sampling_prob_col = NULL,
    precomputed_split = NULL
) {
  preprocess_method <- match.arg(preprocess_method)
  metadata_cols <- unique(metadata_cols)
  metadata_cols <- metadata_cols[!is.na(metadata_cols) & nzchar(metadata_cols)]
  
  if (length(metadata_cols) > 0 && !all(metadata_cols %in% colnames(dataset))) {
    stop("Some metadata_cols are missing from dataset.")
  }
  
  if (!is.null(sampling_prob_col) && !sampling_prob_col %in% colnames(dataset)) {
    stop("sampling_prob_col is missing from dataset.")
  }
  
  if (is.null(precomputed_split)) {
    set.seed(temp_random_seed)
    dataset <- dataset[sample(nrow(dataset)), , drop = FALSE]
  }
  
  excluded_cols <- c(condition_col, metadata_cols, sampling_prob_col)
  feature <- names(dataset)[!names(dataset) %in% excluded_cols]
  
  if (length(feature) == 0) {
    stop("No feature columns found after excluding condition column.")
  }
  
  condition_levels <- levels(factor(dataset[[condition_col]]))
  if (length(condition_levels) != 2) {
    stop("This function currently expects exactly 2 condition levels.")
  }
  
  n_available_per_group <- min(
    sum(dataset[[condition_col]] == condition_levels[1]),
    sum(dataset[[condition_col]] == condition_levels[2])
  )
  
  requested_per_group <- min(
    as.integer(each_condition_total_sample_no),
    as.integer(n_available_per_group)
  )
  
  n_train_per_group <- floor(train_ratio * requested_per_group)
  n_test_per_group <- floor(test_ratio * requested_per_group)
  
  if (n_train_per_group < 1) {
    stop("Requested training set size per condition is < 1.")
  }
  
  if (n_test_per_group < 1) {
    n_test_per_group <- 1L
  }
  
  if (!is.null(precomputed_split)) {
    required_split_names <- c("unit_col", "train_unit_ids", "test_unit_ids")
    missing_split_names <- setdiff(required_split_names, names(precomputed_split))
    if (length(missing_split_names) > 0) {
      stop("precomputed_split is missing: ", paste(missing_split_names, collapse = ", "))
    }
    unit_col <- precomputed_split$unit_col
    if (!unit_col %in% colnames(dataset)) {
      stop("precomputed_split$unit_col is missing from dataset: ", unit_col)
    }
    
    train_indices <- match(precomputed_split$train_unit_ids, as.character(dataset[[unit_col]]))
    test_indices <- match(precomputed_split$test_unit_ids, as.character(dataset[[unit_col]]))
    if (anyNA(train_indices) || anyNA(test_indices)) {
      stop("Some precomputed train/test unit IDs were not found in dataset.")
    }
    if (length(intersect(train_indices, test_indices)) > 0) {
      stop("Precomputed train/test split has overlapping rows.")
    }
    
    train_raw <- dataset[train_indices, , drop = FALSE]
    test_raw <- dataset[test_indices, , drop = FALSE]
  } else {
    idx_cond1 <- which(dataset[[condition_col]] == condition_levels[1])
    idx_cond2 <- which(dataset[[condition_col]] == condition_levels[2])
    
    if (length(idx_cond1) < n_train_per_group || length(idx_cond2) < n_train_per_group) {
      stop("Not enough rows to sample the requested training set size.")
    }
    
    train_prob_condition1 <- if (is.null(sampling_prob_col)) {
      NULL
    } else {
      dataset[[sampling_prob_col]][idx_cond1]
    }
    
    train_prob_condition2 <- if (is.null(sampling_prob_col)) {
      NULL
    } else {
      dataset[[sampling_prob_col]][idx_cond2]
    }
    
    set.seed(temp_random_seed + 1L)
    train_indices_condition1 <- sample(
      idx_cond1,
      size = n_train_per_group,
      prob = train_prob_condition1
    )
    
    train_indices_condition2 <- sample(
      idx_cond2,
      size = n_train_per_group,
      prob = train_prob_condition2
    )
    
    train_indices <- c(train_indices_condition1, train_indices_condition2)
    train_raw <- dataset[train_indices, , drop = FALSE]
    
    remainder <- dataset[-train_indices, , drop = FALSE]
    rem_idx_cond1 <- which(remainder[[condition_col]] == condition_levels[1])
    rem_idx_cond2 <- which(remainder[[condition_col]] == condition_levels[2])
    
    n_test_per_group <- min(
      n_test_per_group,
      length(rem_idx_cond1),
      length(rem_idx_cond2)
    )
    
    if (n_test_per_group < 1) {
      stop("Not enough rows to sample the requested test set size from the remainder.")
    }
    
    test_prob_condition1 <- if (is.null(sampling_prob_col)) {
      NULL
    } else {
      remainder[[sampling_prob_col]][rem_idx_cond1]
    }
    
    test_prob_condition2 <- if (is.null(sampling_prob_col)) {
      NULL
    } else {
      remainder[[sampling_prob_col]][rem_idx_cond2]
    }
    
    set.seed(temp_random_seed + 3L)
    test_indices_condition1 <- sample(
      rem_idx_cond1,
      size = n_test_per_group,
      prob = test_prob_condition1
    )
    
    test_indices_condition2 <- sample(
      rem_idx_cond2,
      size = n_test_per_group,
      prob = test_prob_condition2
    )
    
    test_indices <- c(test_indices_condition1, test_indices_condition2)
    test_raw <- remainder[test_indices, , drop = FALSE]
  }
  
  set.seed(temp_random_seed + 2L)
  train_raw <- train_raw[sample(nrow(train_raw)), , drop = FALSE]
  
  set.seed(temp_random_seed + 4L)
  test_raw <- test_raw[sample(nrow(test_raw)), , drop = FALSE]
  
  keep_cols <- c(condition_col, metadata_cols, feature)
  train_raw <- train_raw[, keep_cols, drop = FALSE]
  test_raw <- test_raw[, keep_cols, drop = FALSE]
  
  if (preprocess_method == "none") {
    train_data <- train_raw
    test_data <- test_raw
    preprocessing_info <- list(method = "none")
  }
  
  if (preprocess_method == "scale") {
    x_train <- as.matrix(train_raw[, feature, drop = FALSE])
    x_test <- as.matrix(test_raw[, feature, drop = FALSE])
    
    train_means <- colMeans(x_train, na.rm = TRUE)
    train_sds <- apply(x_train, 2, stats::sd, na.rm = TRUE)
    train_sds[is.na(train_sds) | train_sds == 0] <- 1
    
    x_train_z <- sweep(x_train, 2, train_means, "-")
    x_train_z <- sweep(x_train_z, 2, train_sds, "/")
    x_test_z <- sweep(x_test, 2, train_means, "-")
    x_test_z <- sweep(x_test_z, 2, train_sds, "/")
    
    x_train_z[is.na(x_train_z)] <- 0
    x_test_z[is.na(x_test_z)] <- 0
    
    train_min <- apply(x_train_z, 2, min, na.rm = TRUE)
    train_max <- apply(x_train_z, 2, max, na.rm = TRUE)
    train_range <- train_max - train_min
    train_range[is.na(train_range) | train_range == 0] <- 1
    
    x_train_scaled <- sweep(x_train_z, 2, train_min, "-")
    x_train_scaled <- sweep(x_train_scaled, 2, train_range, "/")
    x_test_scaled <- sweep(x_test_z, 2, train_min, "-")
    x_test_scaled <- sweep(x_test_scaled, 2, train_range, "/")
    
    x_test_scaled[x_test_scaled < 0] <- 0
    x_test_scaled[x_test_scaled > 1] <- 1
    
    x_train_scaled <- as.data.frame(x_train_scaled, check.names = FALSE)
    x_test_scaled <- as.data.frame(x_test_scaled, check.names = FALSE)
    
    row.names(x_train_scaled) <- row.names(train_raw)
    row.names(x_test_scaled) <- row.names(test_raw)
    
    train_data <- data.frame(
      condition = train_raw[[condition_col]],
      train_raw[, metadata_cols, drop = FALSE],
      x_train_scaled,
      check.names = FALSE
    )
    
    test_data <- data.frame(
      condition = test_raw[[condition_col]],
      test_raw[, metadata_cols, drop = FALSE],
      x_test_scaled,
      check.names = FALSE
    )
    
    preprocessing_info <- list(
      method = "scale",
      train_means = train_means,
      train_sds = train_sds,
      train_min = train_min,
      train_max = train_max
    )
  }
  
  if (!all(feature %in% colnames(train_data))) {
    stop("Some features are missing from train_data after preprocessing.")
  }
  
  if (!all(feature %in% colnames(test_data))) {
    stop("Some features are missing from test_data after preprocessing.")
  }
  
  if (anyNA(train_data[, feature, drop = FALSE])) {
    stop("NA values found in train_data after preprocessing.")
  }
  
  if (anyNA(test_data[, feature, drop = FALSE])) {
    stop("NA values found in test_data after preprocessing.")
  }
  
  list(
    train_data = train_data,
    test_data = test_data,
    train_raw = train_raw,
    test_raw = test_raw,
    feature = feature,
    condition_levels = condition_levels,
    preprocessing_info = preprocessing_info
  )
}


#' Build coverage-aware train/test split plans
#'
#' Plans repeated binary train/test splits so cells/metacells with the fewest
#' previous held-out appearances are prioritized for the next test set.
#'
#' @param dataset Data frame containing condition, unit IDs, optional sampling
#'   probabilities, and feature columns.
#' @param random_seeds Integer vector of per-iteration seeds.
#' @param train_ratio Fraction of sampled rows per condition assigned to train.
#' @param test_ratio Fraction of sampled rows per condition assigned to test.
#' @param each_condition_total_sample_no Target sampled rows per condition.
#' @param condition_col Binary condition column.
#' @param unit_col Stable cell/metacell identifier column.
#' @param sampling_prob_col Optional sampling probability column used for train
#'   sampling and for weighted tie-breaking during test selection.
#'
#' @return List containing split plans and a coverage summary.
#' @export
build_coverage_aware_split_plan <- function(
    dataset,
    random_seeds,
    train_ratio,
    test_ratio,
    each_condition_total_sample_no,
    condition_col = "condition",
    unit_col = ".shapture_unit_id",
    sampling_prob_col = NULL
) {
  if (!condition_col %in% colnames(dataset)) {
    stop("condition_col is missing from dataset.")
  }
  if (!unit_col %in% colnames(dataset)) {
    stop("unit_col is missing from dataset.")
  }
  if (!is.null(sampling_prob_col) && !sampling_prob_col %in% colnames(dataset)) {
    stop("sampling_prob_col is missing from dataset.")
  }
  
  unit_ids <- as.character(dataset[[unit_col]])
  if (anyNA(unit_ids) || any(!nzchar(unit_ids))) {
    stop("unit_col must contain non-empty stable identifiers for all rows.")
  }
  if (anyDuplicated(unit_ids)) {
    stop("unit_col must uniquely identify rows for coverage-aware sampling.")
  }
  
  condition_levels <- levels(factor(dataset[[condition_col]]))
  if (length(condition_levels) != 2) {
    stop("Coverage-aware split planning currently expects exactly 2 condition levels.")
  }
  
  n_available_per_group <- min(
    sum(dataset[[condition_col]] == condition_levels[1]),
    sum(dataset[[condition_col]] == condition_levels[2])
  )
  requested_per_group <- min(
    as.integer(each_condition_total_sample_no),
    as.integer(n_available_per_group)
  )
  n_train_per_group <- floor(train_ratio * requested_per_group)
  n_test_per_group <- floor(test_ratio * requested_per_group)
  if (n_train_per_group < 1) {
    stop("Requested training set size per condition is < 1.")
  }
  if (n_test_per_group < 1) {
    n_test_per_group <- 1L
  }
  
  select_lowest_count <- function(idx_pool, n, counts, probs = NULL, seed = 1L) {
    if (n > length(idx_pool)) {
      stop("Requested more test rows than available in coverage-aware split planning.")
    }
    set.seed(seed)
    selected <- integer(0)
    remaining <- idx_pool
    
    while (length(selected) < n) {
      remaining_counts <- counts[unit_ids[remaining]]
      min_count <- min(remaining_counts, na.rm = TRUE)
      candidates <- remaining[remaining_counts == min_count]
      need <- n - length(selected)
      take_n <- min(need, length(candidates))
      candidate_prob <- NULL
      if (!is.null(probs)) {
        candidate_prob <- probs[match(candidates, idx_pool)]
        if (all(!is.finite(candidate_prob)) || sum(candidate_prob, na.rm = TRUE) <= 0) {
          candidate_prob <- NULL
        }
      }
      chosen <- sample(candidates, size = take_n, replace = FALSE, prob = candidate_prob)
      selected <- c(selected, chosen)
      remaining <- setdiff(remaining, chosen)
    }
    
    selected
  }
  
  test_counts <- stats::setNames(rep(0L, nrow(dataset)), unit_ids)
  train_counts <- stats::setNames(rep(0L, nrow(dataset)), unit_ids)
  plans <- vector("list", length(random_seeds))
  
  for (iter_idx in seq_along(random_seeds)) {
    seed <- as.integer(random_seeds[iter_idx])
    train_indices <- integer(0)
    test_indices <- integer(0)
    
    for (level_idx in seq_along(condition_levels)) {
      condition_level <- condition_levels[level_idx]
      idx_cond <- which(dataset[[condition_col]] == condition_level)
      if (length(idx_cond) < n_train_per_group + n_test_per_group) {
        stop("Not enough rows in one condition to build non-overlapping train/test sets.")
      }
      
      test_probs <- if (is.null(sampling_prob_col)) {
        NULL
      } else {
        dataset[[sampling_prob_col]][idx_cond]
      }
      test_idx <- select_lowest_count(
        idx_pool = idx_cond,
        n = n_test_per_group,
        counts = test_counts,
        probs = test_probs,
        seed = seed + 1000L * level_idx
      )
      
      remaining_idx <- setdiff(idx_cond, test_idx)
      train_probs <- if (is.null(sampling_prob_col)) {
        NULL
      } else {
        dataset[[sampling_prob_col]][remaining_idx]
      }
      if (!is.null(train_probs) &&
          (all(!is.finite(train_probs)) || sum(train_probs, na.rm = TRUE) <= 0)) {
        train_probs <- NULL
      }
      set.seed(seed + 2000L * level_idx)
      train_idx <- sample(
        remaining_idx,
        size = n_train_per_group,
        replace = FALSE,
        prob = train_probs
      )
      
      test_counts[unit_ids[test_idx]] <- test_counts[unit_ids[test_idx]] + 1L
      train_counts[unit_ids[train_idx]] <- train_counts[unit_ids[train_idx]] + 1L
      test_indices <- c(test_indices, test_idx)
      train_indices <- c(train_indices, train_idx)
    }
    
    plans[[iter_idx]] <- list(
      unit_col = unit_col,
      train_unit_ids = unit_ids[train_indices],
      test_unit_ids = unit_ids[test_indices]
    )
  }
  
  coverage_summary <- do.call(rbind, lapply(condition_levels, function(condition_level) {
    idx <- which(dataset[[condition_col]] == condition_level)
    ids <- unit_ids[idx]
    data.frame(
      condition = condition_level,
      total_units = length(ids),
      unique_test_units = sum(test_counts[ids] > 0),
      percent_test_units = round(100 * mean(test_counts[ids] > 0), 2),
      min_test_appearances = min(test_counts[ids]),
      median_test_appearances = stats::median(test_counts[ids]),
      max_test_appearances = max(test_counts[ids]),
      unique_train_units = sum(train_counts[ids] > 0),
      percent_train_units = round(100 * mean(train_counts[ids] > 0), 2),
      min_train_appearances = min(train_counts[ids]),
      median_train_appearances = stats::median(train_counts[ids]),
      max_train_appearances = max(train_counts[ids]),
      train_rows_across_iterations = sum(train_counts[ids]),
      test_rows_across_iterations = sum(test_counts[ids]),
      stringsAsFactors = FALSE
    )
  }))
  
  list(
    plans = plans,
    coverage_summary = coverage_summary,
    test_counts = test_counts,
    train_counts = train_counts
  )
}


#' Build coverage-aware feature subset plans
#'
#' Plans repeated feature subsets so features with the fewest previous model
#' inclusions are prioritized for the next iteration.
#'
#' @param feature_names Character vector of candidate features.
#' @param random_seeds Integer vector of per-iteration seeds.
#' @param fraction Fraction of candidate features to include in each iteration.
#'
#' @return List containing per-iteration feature plans and an aggregate summary.
#' @export
build_coverage_aware_feature_plan <- function(
    feature_names,
    random_seeds,
    fraction = 0.20
) {
  feature_names <- unique(as.character(feature_names))
  feature_names <- feature_names[!is.na(feature_names) & nzchar(feature_names)]
  if (length(feature_names) == 0) {
    stop("feature_names must contain at least one feature.")
  }
  if (!is.numeric(fraction) ||
      length(fraction) != 1 ||
      !is.finite(fraction) ||
      fraction <= 0 ||
      fraction > 1) {
    stop("fraction must be one finite numeric value in (0, 1].")
  }
  
  n_per_iteration <- max(1L, ceiling(length(feature_names) * fraction))
  n_per_iteration <- min(n_per_iteration, length(feature_names))
  inclusion_counts <- stats::setNames(rep(0L, length(feature_names)), feature_names)
  plans <- vector("list", length(random_seeds))
  
  for (iter_idx in seq_along(random_seeds)) {
    seed <- as.integer(random_seeds[iter_idx])
    set.seed(seed + 910000L)
    selected <- character(0)
    remaining <- feature_names
    
    while (length(selected) < n_per_iteration) {
      remaining_counts <- inclusion_counts[remaining]
      min_count <- min(remaining_counts, na.rm = TRUE)
      candidates <- remaining[remaining_counts == min_count]
      need <- n_per_iteration - length(selected)
      take_n <- min(need, length(candidates))
      chosen <- sample(candidates, size = take_n, replace = FALSE)
      selected <- c(selected, chosen)
      remaining <- setdiff(remaining, chosen)
    }
    
    inclusion_counts[selected] <- inclusion_counts[selected] + 1L
    plans[[iter_idx]] <- selected
  }
  
  summary <- data.frame(
    total_candidate_features = length(feature_names),
    feature_fraction_per_iteration = fraction,
    features_per_iteration = n_per_iteration,
    min_gene_inclusion_count = min(inclusion_counts),
    median_gene_inclusion_count = stats::median(inclusion_counts),
    max_gene_inclusion_count = max(inclusion_counts),
    unique_features_included = sum(inclusion_counts > 0),
    percent_features_included = round(100 * mean(inclusion_counts > 0), 2),
    stringsAsFactors = FALSE
  )
  
  list(
    plans = plans,
    summary = summary,
    inclusion_counts = inclusion_counts
  )
}


#' Build balance groups for SHAPTURE sampling and weighting
#'
#' @param dataset Data frame to augment.
#' @param comparison_mode One of `"condition_vs_condition"` or
#'   `"one_vs_others"`.
#' @param sample_col Sample identifier column.
#' @param condition_col Condition column.
#' @param split_by_col Split-defining column used in one-vs-others mode.
#' @param positive_label Positive class label.
#' @param balance_col Name of the output balance-group column.
#'
#' @return The input data frame with an added balance-group column.
#' @export
build_balance_groups <- function(
    dataset,
    comparison_mode = c("condition_vs_condition", "one_vs_others"),
    sample_col = "sample_id",
    condition_col = "condition",
    split_by_col = NULL,
    positive_label = NULL,
    balance_col = "balance_group"
) {
  if (identical(comparison_mode, "one_vs_rest")) {
    comparison_mode <- "one_vs_others"
  }
  comparison_mode <- match.arg(comparison_mode)

  if (!sample_col %in% colnames(dataset)) {
    stop("sample_col is missing from dataset.")
  }

  if (!condition_col %in% colnames(dataset)) {
    stop("condition_col is missing from dataset.")
  }

  if (comparison_mode == "condition_vs_condition") {
    dataset[[balance_col]] <- as.character(dataset[[sample_col]])
    return(dataset)
  }

  if (is.null(split_by_col) || !split_by_col %in% colnames(dataset)) {
    stop("split_by_col must be present in dataset for one_vs_others balancing.")
  }

  if (is.null(positive_label)) {
    stop("positive_label must be provided for one_vs_others balancing.")
  }

  sample_ids <- as.character(dataset[[sample_col]])
  split_levels <- as.character(dataset[[split_by_col]])
  class_labels <- as.character(dataset[[condition_col]])

  is_positive <- class_labels == as.character(positive_label)
  dataset[[balance_col]] <- ifelse(
    is_positive,
    sample_ids,
    paste(sample_ids, split_levels, sep = "__")
  )

  dataset
}


#' Downsample rows within each condition using balance groups
#'
#' @param dataset Data frame to downsample.
#' @param condition_col Condition column.
#' @param balance_col Balance-group column.
#' @param seed Integer seed for reproducibility.
#'
#' @return Downsampled data frame.
#' @export
downsample_by_condition_and_group <- function(
    dataset,
    condition_col = "condition",
    balance_col = "balance_group",
    seed = 42
) {
  if (!condition_col %in% colnames(dataset)) {
    stop("condition_col is missing from dataset.")
  }

  if (!balance_col %in% colnames(dataset)) {
    stop("balance_col is missing from dataset.")
  }

  set.seed(seed)

  condition_levels <- unique(dataset[[condition_col]])
  downsampled <- lapply(condition_levels, function(this_condition) {
    sub <- dataset[dataset[[condition_col]] == this_condition, , drop = FALSE]
    group_counts <- table(sub[[balance_col]])
    n_target <- as.integer(stats::median(group_counts))

    if (is.na(n_target) || n_target < 1) {
      return(sub)
    }

    group_indices <- split(seq_len(nrow(sub)), sub[[balance_col]])
    keep_idx <- unlist(
      lapply(group_indices, function(idx) {
        if (length(idx) <= n_target) {
          idx
        } else {
          sample(idx, size = n_target)
        }
      }),
      use.names = FALSE
    )
    sub[keep_idx, , drop = FALSE]
  })

  out <- dplyr::bind_rows(downsampled)
  row.names(out) <- NULL
  out
}


#' Add inverse-group sampling probabilities
#'
#' @param dataset Data frame to augment.
#' @param group_col Balance-group column.
#' @param condition_col Condition column used for per-condition balancing.
#' @param prob_col Name of the output probability column.
#'
#' @return The input data frame with an added probability column.
#' @export
add_inverse_group_probabilities <- function(
    dataset,
    group_col = "balance_group",
    condition_col = "condition",
    prob_col = "sampling_prob"
) {
  if (!group_col %in% colnames(dataset)) {
    stop("group_col is missing from dataset.")
  }
  
  if (!condition_col %in% colnames(dataset)) {
    stop("condition_col is missing from dataset.")
  }
  
  dataset[[prob_col]] <- NA_real_
  condition_levels <- unique(dataset[[condition_col]])
  
  for (this_condition in condition_levels) {
    idx <- which(dataset[[condition_col]] == this_condition)
    group_counts <- table(dataset[[group_col]][idx])
    inv_counts <- 1 / as.numeric(group_counts[dataset[[group_col]][idx]])
    dataset[[prob_col]][idx] <- inv_counts / mean(inv_counts, na.rm = TRUE)
  }
  
  dataset
}


#' Add inverse-sample sampling probabilities
#'
#' @param dataset Data frame to augment.
#' @param sample_col Sample identifier column.
#' @param condition_col Condition column used for per-condition balancing.
#' @param prob_col Name of the output probability column.
#'
#' @return The input data frame with an added probability column.
#' @export
add_inverse_sample_probabilities <- function(
    dataset,
    sample_col = "sample_id",
    condition_col = "condition",
    prob_col = "sampling_prob"
) {
  add_inverse_group_probabilities(
    dataset = dataset,
    group_col = sample_col,
    condition_col = condition_col,
    prob_col = prob_col
  )
}


#' Compute inverse-group model weights
#'
#' @param group_ids Vector of balance-group identifiers.
#'
#' @return Numeric weights normalized to mean 1.
#' @export
compute_model_weights_by_group <- function(group_ids) {
  group_counts <- table(group_ids)
  weights <- 1 / as.numeric(group_counts[group_ids])
  weights / mean(weights, na.rm = TRUE)
}


#' Compute inverse-sample model weights
#'
#' @param sample_ids Vector of sample identifiers.
#'
#' @return Numeric weights normalized to mean 1.
#' @export
compute_model_weights_by_sample <- function(sample_ids) {
  compute_model_weights_by_group(sample_ids)
}


#' Compute a dynamic pseudocount from expression scale
#'
#' @param x Numeric vector or matrix of expression values.
#' @param fraction Fraction of the median positive value to use.
#' @param fallback Minimum fallback value.
#'
#' @return A positive scalar pseudocount.
#' @export
compute_dynamic_pseudocount <- function(
    x,
    fraction = 0.01,
    fallback = 1e-6
) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x > 0]
  
  if (length(x) == 0) {
    return(fallback)
  }
  
  max(fraction * stats::median(x, na.rm = TRUE), fallback)
}
