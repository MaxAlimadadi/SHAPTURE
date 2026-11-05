#' Validate the SHAPTURE analysis seed
#'
#' @param analysis_seed Positive integer used as the reproducibility root.
#'
#' @return The validated integer seed.
#' @keywords internal
#' @noRd
validate_analysis_seed <- function(analysis_seed) {
  if (!is.numeric(analysis_seed) || length(analysis_seed) != 1L ||
      !is.finite(analysis_seed) || analysis_seed < 1 ||
      analysis_seed != floor(analysis_seed) || analysis_seed > 999999999) {
    stop("analysis_seed must be one whole number between 1 and 999999999.")
  }

  as.integer(analysis_seed)
}


#' Derive a stable child seed from an analysis seed
#'
#' Uses a small platform-independent string hash so seed assignment depends on
#' stable task and stage identifiers, not task order or parallel scheduling.
#'
#' @param analysis_seed Positive integer reproducibility root.
#' @param ... Values that uniquely identify the task or stochastic stage.
#'
#' @return A positive integer suitable for `set.seed()` and model backends.
#' @keywords internal
#' @noRd
derive_shapture_seed <- function(analysis_seed, ...) {
  analysis_seed <- validate_analysis_seed(analysis_seed)
  seed_parts <- list(...)
  serialize_part <- function(x) {
    if (is.null(x)) {
      return("<NULL>")
    }
    x <- as.character(x)
    x[is.na(x)] <- "<NA>"
    paste(x, collapse = "\u001e")
  }
  key <- paste(
    c(as.character(analysis_seed), vapply(seed_parts, serialize_part, character(1))),
    collapse = "\u001f"
  )

  # Keep child seeds below one billion so existing stage-specific offsets
  # cannot overflow R's integer seed range.
  modulus <- 999999937
  hash <- as.double(analysis_seed %% modulus)
  for (code in utf8ToInt(enc2utf8(key))) {
    hash <- (hash * 131 + code) %% modulus
  }
  for (mix_round in seq_len(3L)) {
    hash <- (hash * 48271 + 12820163 + mix_round) %% modulus
  }

  as.integer(hash + 1)
}


#' Derive a stable non-repeating seed sequence
#'
#' @param root_seed Positive integer seed for the task or analysis.
#' @param stage Stable stochastic-stage label.
#' @param n Number of seeds to generate.
#'
#' @return Integer vector of length `n`.
#' @keywords internal
#' @noRd
derive_shapture_seed_sequence <- function(root_seed, stage, n) {
  root_seed <- validate_analysis_seed(root_seed)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 1 || n != floor(n) || n > 999999937) {
    stop("n must be one whole number between 1 and 999999937.")
  }
  n <- as.integer(n)
  modulus <- 999999937
  current <- as.double(derive_shapture_seed(root_seed, stage, "start") - 1L)
  stride <- as.double(
    (derive_shapture_seed(root_seed, stage, "stride") - 1L) %%
      (modulus - 1) + 1L
  )
  out <- integer(n)
  for (i in seq_len(n)) {
    out[[i]] <- as.integer(current + 1)
    current <- (current + stride) %% modulus
  }
  out
}


#' Build the deterministic seed plan for one SHAPTURE task
#'
#' @param analysis_seed Positive integer reproducibility root.
#' @param comparison_mode Comparison mode.
#' @param split_level_column Split-defining column or `NULL`.
#' @param temp_split Current task level.
#' @param task_what_0 Negative/reference label.
#' @param task_what_1 Positive/target label.
#' @param iteration_no Number of model iterations.
#'
#' @return A list containing task-stage seeds and an iteration seed table.
#' @keywords internal
#' @noRd
build_shapture_task_seed_plan <- function(
    analysis_seed,
    comparison_mode,
    split_level_column,
    temp_split,
    task_what_0,
    task_what_1,
    iteration_no
) {
  analysis_seed <- validate_analysis_seed(analysis_seed)
  if (!is.numeric(iteration_no) || length(iteration_no) != 1L ||
      !is.finite(iteration_no) || iteration_no < 1 ||
      iteration_no != floor(iteration_no)) {
    stop("iteration_no must be one positive whole number.")
  }
  iteration_no <- as.integer(iteration_no)

  task_seed <- derive_shapture_seed(
    analysis_seed,
    "task",
    comparison_mode,
    split_level_column,
    temp_split,
    task_what_0,
    task_what_1
  )
  iteration <- seq_len(iteration_no)
  iteration_seeds <- data.frame(
    iteration = iteration,
    cell_sampling_seed = derive_shapture_seed_sequence(
      task_seed, "cell_sampling", iteration_no
    ),
    feature_selection_seed = derive_shapture_seed_sequence(
      task_seed, "feature_selection", iteration_no
    ),
    model_seed = derive_shapture_seed_sequence(
      task_seed, "model", iteration_no
    ),
    stringsAsFactors = FALSE
  )

  list(
    analysis_seed = analysis_seed,
    task_seed = task_seed,
    balance_seed = derive_shapture_seed(task_seed, "patient_balance"),
    cross_patient_bootstrap_seed = derive_shapture_seed(
      task_seed,
      "cross_patient_bootstrap"
    ),
    iteration_seeds = iteration_seeds
  )
}


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
        candidate_prob[!is.finite(candidate_prob)] <- 0
        if (sum(candidate_prob) <= 0) {
          candidate_prob <- NULL
        }
      }
      chosen <- candidates[sample.int(
        length(candidates),
        size = take_n,
        replace = FALSE,
        prob = candidate_prob
      )]
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
          any(!is.finite(train_probs))) {
        train_probs[!is.finite(train_probs)] <- 0
      }
      if (!is.null(train_probs) && sum(train_probs) <= 0) {
        train_probs <- NULL
      }
      set.seed(seed + 2000L * level_idx)
      train_idx <- remaining_idx[sample.int(
        length(remaining_idx),
        size = n_train_per_group,
        replace = FALSE,
        prob = train_probs
      )]
      
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


#' Select patients with balanced held-out coverage
#'
#' @param pool Candidate patient IDs.
#' @param n Number to select.
#' @param test_counts Named patient-level test counts.
#' @param already_selected Patient IDs already selected in the current fold.
#' @param pair_counts Symmetric matrix of prior test-pair counts.
#' @param seed Random seed for deterministic tie-breaking.
#'
#' @return Character vector of selected patient IDs.
#' @noRd
.select_coverage_aware_patients <- function(
    pool,
    n,
    test_counts,
    already_selected = character(0),
    pair_counts,
    seed
) {
  pool <- unique(as.character(pool))
  pool <- pool[!is.na(pool) & nzchar(pool)]
  if (n == 0) {
    return(character(0))
  }
  if (length(pool) < n) {
    stop("Not enough patients in a required eligibility pool.")
  }

  set.seed(seed)
  selected <- character(0)
  for (slot in seq_len(n)) {
    candidates <- setdiff(pool, selected)
    candidate_counts <- test_counts[candidates]
    candidates <- candidates[candidate_counts == min(candidate_counts)]

    current_selected <- unique(c(already_selected, selected))
    if (length(current_selected) > 0 && length(candidates) > 1) {
      pair_score <- rowSums(
        pair_counts[candidates, current_selected, drop = FALSE]
      )
      candidates <- candidates[pair_score == min(pair_score)]
    }

    selected <- c(selected, sample(candidates, size = 1L))
  }
  selected
}


#' Allocate a balanced row quota across selected patients
#'
#' @param available Named available-row counts by patient.
#' @param total Total rows to allocate.
#' @param minimum_each Minimum rows required from every selected patient.
#' @param seed Random seed for deterministic tie-breaking.
#'
#' @return Named integer quota vector.
#' @noRd
.allocate_patient_row_quotas <- function(available, total, minimum_each, seed) {
  available_names <- names(available)
  available <- as.integer(available)
  names(available) <- available_names
  n_patients <- length(available)
  if (n_patients == 0) {
    stop("No patients are available for row allocation.")
  }
  if (any(available < minimum_each)) {
    stop("A selected patient has fewer rows than the required minimum.")
  }
  if (total < n_patients * minimum_each || total > sum(available)) {
    stop("Requested row total is incompatible with patient-level row limits.")
  }

  quotas <- stats::setNames(rep(as.integer(minimum_each), n_patients), names(available))
  remaining <- as.integer(total - sum(quotas))
  set.seed(seed)
  while (remaining > 0L) {
    candidates <- names(quotas)[quotas < available]
    if (length(candidates) == 0) {
      stop("Unable to complete patient-balanced row allocation.")
    }
    min_quota <- min(quotas[candidates])
    candidates <- candidates[quotas[candidates] == min_quota]
    chosen <- sample(candidates, size = 1L)
    quotas[chosen] <- quotas[chosen] + 1L
    remaining <- remaining - 1L
  }
  quotas
}


#' Select rows within patients using unit-level coverage
#'
#' @param dataset Task data.
#' @param condition_level Current condition level.
#' @param patient_ids Selected patient IDs.
#' @param total_rows Total rows to select.
#' @param minimum_each Minimum rows required from every selected patient.
#' @param balance_across_patients Whether to allocate nearly equal row quotas
#'   across patients. This is used for held-out test data; training keeps the
#'   existing sampling and model-balancing behavior after representing each
#'   selected patient.
#' @param unit_counts Named prior unit-level appearance counts.
#' @param condition_col Condition column.
#' @param sample_col Patient column.
#' @param unit_col Stable unit-ID column.
#' @param sampling_prob_col Optional row sampling-probability column.
#' @param seed Random seed.
#'
#' @return Integer row indices.
#' @noRd
.select_patient_balanced_rows <- function(
    dataset,
    condition_level,
    patient_ids,
    total_rows,
    minimum_each,
    balance_across_patients = TRUE,
    unit_counts,
    condition_col,
    sample_col,
    unit_col,
    sampling_prob_col,
    seed
) {
  patient_ids <- unique(as.character(patient_ids))
  condition_values <- as.character(dataset[[condition_col]])
  sample_values <- as.character(dataset[[sample_col]])
  unit_ids <- as.character(dataset[[unit_col]])

  pools <- lapply(patient_ids, function(patient_id) {
    which(
      condition_values == as.character(condition_level) &
        sample_values == patient_id
    )
  })
  names(pools) <- patient_ids
  available <- vapply(pools, length, integer(1))
  if (!balance_across_patients) {
    if (any(available < minimum_each) ||
        total_rows < length(patient_ids) * minimum_each ||
        total_rows > sum(available)) {
      stop("Requested training rows are incompatible with the selected patients.")
    }

    set.seed(seed)
    selected <- integer(0)
    for (patient_idx in seq_along(patient_ids)) {
      pool <- pools[[patient_ids[patient_idx]]]
      patient_prob <- if (is.null(sampling_prob_col)) {
        NULL
      } else {
        dataset[[sampling_prob_col]][pool]
      }
      if (!is.null(patient_prob)) {
        patient_prob[!is.finite(patient_prob)] <- 0
        if (sum(patient_prob) <= 0) {
          patient_prob <- NULL
        }
      }
      chosen <- pool[sample.int(
        length(pool),
        size = minimum_each,
        replace = FALSE,
        prob = patient_prob
      )]
      selected <- c(selected, chosen)
    }

    remaining_needed <- total_rows - length(selected)
    if (remaining_needed > 0L) {
      remaining <- setdiff(unlist(pools, use.names = FALSE), selected)
      remaining_prob <- if (is.null(sampling_prob_col)) {
        NULL
      } else {
        dataset[[sampling_prob_col]][remaining]
      }
      if (!is.null(remaining_prob)) {
        remaining_prob[!is.finite(remaining_prob)] <- 0
        if (sum(remaining_prob) <= 0) {
          remaining_prob <- NULL
        }
      }
      selected <- c(
        selected,
        remaining[sample.int(
          length(remaining),
          size = remaining_needed,
          replace = FALSE,
          prob = remaining_prob
        )]
      )
    }
    return(selected)
  }

  quotas <- .allocate_patient_row_quotas(
    available = available,
    total = total_rows,
    minimum_each = minimum_each,
    seed = seed
  )

  selected <- integer(0)
  for (patient_idx in seq_along(patient_ids)) {
    patient_id <- patient_ids[patient_idx]
    pool <- pools[[patient_id]]
    quota <- quotas[[patient_id]]
    set.seed(seed + patient_idx * 1009L)
    patient_selected <- integer(0)
    remaining <- pool

    while (length(patient_selected) < quota) {
      remaining_counts <- unit_counts[unit_ids[remaining]]
      candidates <- remaining[remaining_counts == min(remaining_counts)]
      need <- quota - length(patient_selected)
      take_n <- min(need, length(candidates))
      candidate_prob <- NULL
      if (!is.null(sampling_prob_col)) {
        candidate_prob <- dataset[[sampling_prob_col]][candidates]
        candidate_prob[!is.finite(candidate_prob)] <- 0
        if (sum(candidate_prob) <= 0) {
          candidate_prob <- NULL
        }
      }
      chosen <- candidates[sample.int(
        length(candidates),
        size = take_n,
        replace = FALSE,
        prob = candidate_prob
      )]
      patient_selected <- c(patient_selected, chosen)
      remaining <- setdiff(remaining, chosen)
    }
    selected <- c(selected, patient_selected)
  }
  selected
}


#' Build coverage-aware unseen-patient split plans
#'
#' Constructs repeated train/test plans with global patient separation. Patients
#' are eligible only when they meet the task-specific minimum row count in a
#' condition. Sixty percent of eligible patients are assigned to training and
#' the remaining forty percent to testing, with rounded counts and a minimum of
#' three training and two testing patients per condition.
#'
#' @param dataset Task data containing condition, patient, and stable unit IDs.
#' @param random_seeds Integer vector of per-iteration seeds.
#' @param train_ratio Existing row-level training fraction used to retain the
#'   current SHAPTURE row budget.
#' @param test_ratio Existing row-level testing fraction used to retain the
#'   current SHAPTURE row budget.
#' @param each_condition_total_sample_no Current target row budget per condition.
#' @param condition_col Binary condition column.
#' @param sample_col Patient/sample identifier column.
#' @param unit_col Stable cell/metacell identifier column.
#' @param sampling_prob_col Optional row sampling-probability column.
#' @param min_patients_per_condition Minimum eligible patients per condition.
#' @param patient_train_fraction Fraction of eligible patients assigned to train.
#' @param min_units_per_patient Minimum available and selected test units per
#'   patient and condition.
#' @param min_train_units_per_condition Minimum training rows per condition.
#'
#' @return A list with `success`, split plans, patient and unit coverage, and
#'   audit metadata. Eligibility failures return `success = FALSE` and a reason.
#' @export
build_unseen_patient_split_plan <- function(
    dataset,
    random_seeds,
    train_ratio,
    test_ratio,
    each_condition_total_sample_no,
    condition_col = "condition",
    sample_col = "sample_id",
    unit_col = ".shapture_unit_id",
    sampling_prob_col = NULL,
    min_patients_per_condition = 5L,
    patient_train_fraction = 0.60,
    min_units_per_patient = 20L,
    min_train_units_per_condition = 50L
) {
  if (!is.data.frame(dataset)) {
    stop("dataset must be a data frame.")
  }
  if (!is.numeric(random_seeds) || length(random_seeds) < 1L ||
      any(!is.finite(random_seeds))) {
    stop("random_seeds must contain at least one finite numeric seed.")
  }
  validate_positive_integer <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < 1 || x != as.integer(x)) {
      stop(name, " must be one positive integer.")
    }
  }
  validate_fraction <- function(x, name, allow_one = TRUE) {
    upper_ok <- if (allow_one) x <= 1 else x < 1
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x <= 0 || !upper_ok) {
      stop(name, " must be one value in (0, ", if (allow_one) "1]." else "1).")
    }
  }
  validate_fraction(train_ratio, "train_ratio")
  validate_fraction(test_ratio, "test_ratio")
  validate_fraction(patient_train_fraction, "patient_train_fraction", allow_one = FALSE)
  validate_positive_integer(each_condition_total_sample_no, "each_condition_total_sample_no")
  validate_positive_integer(min_patients_per_condition, "min_patients_per_condition")
  validate_positive_integer(min_units_per_patient, "min_units_per_patient")
  validate_positive_integer(min_train_units_per_condition, "min_train_units_per_condition")

  required_cols <- c(condition_col, sample_col, unit_col)
  if (!all(required_cols %in% colnames(dataset))) {
    stop("Unseen-patient split planning is missing required columns: ",
         paste(setdiff(required_cols, colnames(dataset)), collapse = ", "))
  }
  if (!is.null(sampling_prob_col) && !sampling_prob_col %in% colnames(dataset)) {
    stop("sampling_prob_col is missing from dataset.")
  }
  unit_ids <- as.character(dataset[[unit_col]])
  if (anyNA(unit_ids) || any(!nzchar(unit_ids)) || anyDuplicated(unit_ids)) {
    stop("unit_col must uniquely identify every row for unseen-patient planning.")
  }
  sample_values <- as.character(dataset[[sample_col]])
  condition_values <- as.character(dataset[[condition_col]])
  if (anyNA(sample_values) || any(!nzchar(sample_values))) {
    stop("sample_col must contain non-empty patient IDs for unseen-patient planning.")
  }
  if (anyNA(condition_values) || any(!nzchar(condition_values))) {
    stop("condition_col must contain non-empty values for unseen-patient planning.")
  }

  condition_levels <- levels(factor(dataset[[condition_col]]))
  if (length(condition_levels) != 2L) {
    stop("Unseen-patient split planning currently expects exactly 2 condition levels.")
  }

  count_table <- as.data.frame(
    table(
      patient_id = sample_values,
      condition = condition_values
    ),
    stringsAsFactors = FALSE
  )
  names(count_table)[3] <- "available_units"
  count_table <- count_table[count_table$available_units > 0, , drop = FALSE]
  count_table$eligible <- count_table$available_units >= min_units_per_patient

  eligible_by_condition <- stats::setNames(
    lapply(condition_levels, function(level) {
      as.character(count_table$patient_id[
        count_table$condition == level & count_table$eligible
      ])
    }),
    condition_levels
  )
  eligible_counts <- vapply(eligible_by_condition, length, integer(1))
  if (any(eligible_counts < min_patients_per_condition)) {
    reason <- paste0(
      "unseen_patients requires at least ", min_patients_per_condition,
      " eligible patients per condition after requiring at least ",
      min_units_per_patient, " units per patient; observed ",
      paste(paste0(names(eligible_counts), "=", eligible_counts), collapse = ", "),
      "."
    )
    return(list(
      success = FALSE,
      reason = reason,
      patient_coverage_summary = count_table
    ))
  }

  n_train_patients <- as.integer(round(patient_train_fraction * eligible_counts))
  n_test_patients <- eligible_counts - n_train_patients
  if (any(n_train_patients < 3L) || any(n_test_patients < 2L)) {
    reason <- paste0(
      "The rounded 60/40 patient split does not leave at least three training ",
      "and two testing patients per condition."
    )
    return(list(
      success = FALSE,
      reason = reason,
      patient_coverage_summary = count_table
    ))
  }

  eligible_0 <- eligible_by_condition[[1]]
  eligible_1 <- eligible_by_condition[[2]]
  shared <- intersect(eligible_0, eligible_1)
  only_0 <- setdiff(eligible_0, shared)
  only_1 <- setdiff(eligible_1, shared)
  n_test_0 <- n_test_patients[[1]]
  n_test_1 <- n_test_patients[[2]]
  min_shared_test <- max(0L, n_test_0 - length(only_0), n_test_1 - length(only_1))
  max_shared_test <- min(n_test_0, n_test_1, length(shared))
  if (min_shared_test > max_shared_test) {
    return(list(
      success = FALSE,
      reason = "No globally separated patient partition can satisfy the rounded 60/40 split in both conditions.",
      patient_coverage_summary = count_table
    ))
  }
  n_shared_test <- max_shared_test
  n_only_0_test <- n_test_0 - n_shared_test
  n_only_1_test <- n_test_1 - n_shared_test

  all_eligible_patients <- sort(unique(c(eligible_0, eligible_1)))
  patient_test_counts <- stats::setNames(
    rep(0L, length(all_eligible_patients)),
    all_eligible_patients
  )
  pair_counts <- matrix(
    0L,
    nrow = length(all_eligible_patients),
    ncol = length(all_eligible_patients),
    dimnames = list(all_eligible_patients, all_eligible_patients)
  )
  patient_plans <- vector("list", length(random_seeds))

  for (iter_idx in seq_along(random_seeds)) {
    seed <- as.integer(random_seeds[iter_idx])
    selected_shared <- .select_coverage_aware_patients(
      pool = shared,
      n = n_shared_test,
      test_counts = patient_test_counts,
      pair_counts = pair_counts,
      seed = seed + 300001L
    )
    selected_only_0 <- .select_coverage_aware_patients(
      pool = only_0,
      n = n_only_0_test,
      test_counts = patient_test_counts,
      already_selected = selected_shared,
      pair_counts = pair_counts,
      seed = seed + 300101L
    )
    selected_only_1 <- .select_coverage_aware_patients(
      pool = only_1,
      n = n_only_1_test,
      test_counts = patient_test_counts,
      already_selected = c(selected_shared, selected_only_0),
      pair_counts = pair_counts,
      seed = seed + 300201L
    )
    test_patients_global <- unique(c(
      selected_shared,
      selected_only_0,
      selected_only_1
    ))
    test_by_condition <- list(
      intersect(test_patients_global, eligible_0),
      intersect(test_patients_global, eligible_1)
    )
    names(test_by_condition) <- condition_levels
    train_by_condition <- list(
      setdiff(eligible_0, test_patients_global),
      setdiff(eligible_1, test_patients_global)
    )
    names(train_by_condition) <- condition_levels

    if (any(vapply(test_by_condition, length, integer(1)) != n_test_patients) ||
        any(vapply(train_by_condition, length, integer(1)) != n_train_patients)) {
      stop("Internal patient-partition counts do not match the requested 60/40 split.")
    }

    patient_plans[[iter_idx]] <- list(
      test_patients_global = test_patients_global,
      test_by_condition = test_by_condition,
      train_by_condition = train_by_condition
    )
    patient_test_counts[test_patients_global] <-
      patient_test_counts[test_patients_global] + 1L
    if (length(test_patients_global) > 1L) {
      pair_grid <- utils::combn(test_patients_global, 2L)
      for (pair_idx in seq_len(ncol(pair_grid))) {
        p1 <- pair_grid[1, pair_idx]
        p2 <- pair_grid[2, pair_idx]
        pair_counts[p1, p2] <- pair_counts[p1, p2] + 1L
        pair_counts[p2, p1] <- pair_counts[p2, p1] + 1L
      }
    }
  }

  available_train <- matrix(0L, nrow = length(patient_plans), ncol = 2L)
  available_test <- matrix(0L, nrow = length(patient_plans), ncol = 2L)
  for (iter_idx in seq_along(patient_plans)) {
    for (condition_idx in seq_along(condition_levels)) {
      level <- condition_levels[condition_idx]
      available_train[iter_idx, condition_idx] <- sum(
        condition_values == level &
          sample_values %in% patient_plans[[iter_idx]]$train_by_condition[[level]]
      )
      available_test[iter_idx, condition_idx] <- sum(
        condition_values == level &
          sample_values %in% patient_plans[[iter_idx]]$test_by_condition[[level]]
      )
    }
  }

  n_available_per_group <- min(table(dataset[[condition_col]]))
  requested_per_group <- min(
    as.integer(each_condition_total_sample_no),
    as.integer(n_available_per_group)
  )
  base_train_rows <- max(1L, floor(train_ratio * requested_per_group))
  base_test_rows <- max(1L, floor(test_ratio * requested_per_group))
  minimum_test_rows <- max(n_test_patients * min_units_per_patient)
  train_rows_per_condition <- min(
    max(base_train_rows, min_train_units_per_condition),
    min(available_train)
  )
  test_rows_per_condition <- min(
    max(base_test_rows, minimum_test_rows),
    min(available_test)
  )
  if (train_rows_per_condition < min_train_units_per_condition ||
      test_rows_per_condition < minimum_test_rows) {
    reason <- paste0(
      "No fixed patient-held-out row budget can satisfy at least ",
      min_train_units_per_condition, " training units per condition and ",
      min_units_per_patient, " testing units per selected patient."
    )
    return(list(
      success = FALSE,
      reason = reason,
      patient_coverage_summary = count_table
    ))
  }

  test_unit_counts <- stats::setNames(rep(0L, nrow(dataset)), unit_ids)
  train_unit_counts <- stats::setNames(rep(0L, nrow(dataset)), unit_ids)
  plans <- vector("list", length(random_seeds))
  test_patient_condition_keys <- vector("list", length(random_seeds))

  for (iter_idx in seq_along(random_seeds)) {
    seed <- as.integer(random_seeds[iter_idx])
    train_indices <- integer(0)
    test_indices <- integer(0)
    iteration_keys <- character(0)
    for (condition_idx in seq_along(condition_levels)) {
      level <- condition_levels[condition_idx]
      test_patients <- patient_plans[[iter_idx]]$test_by_condition[[level]]
      train_patients <- patient_plans[[iter_idx]]$train_by_condition[[level]]
      test_idx <- .select_patient_balanced_rows(
        dataset = dataset,
        condition_level = level,
        patient_ids = test_patients,
        total_rows = test_rows_per_condition,
        minimum_each = min_units_per_patient,
        balance_across_patients = TRUE,
        unit_counts = test_unit_counts,
        condition_col = condition_col,
        sample_col = sample_col,
        unit_col = unit_col,
        sampling_prob_col = sampling_prob_col,
        seed = seed + condition_idx * 400001L
      )
      train_idx <- .select_patient_balanced_rows(
        dataset = dataset,
        condition_level = level,
        patient_ids = train_patients,
        total_rows = train_rows_per_condition,
        minimum_each = 1L,
        balance_across_patients = FALSE,
        unit_counts = train_unit_counts,
        condition_col = condition_col,
        sample_col = sample_col,
        unit_col = unit_col,
        sampling_prob_col = sampling_prob_col,
        seed = seed + condition_idx * 500009L
      )
      test_indices <- c(test_indices, test_idx)
      train_indices <- c(train_indices, train_idx)
      iteration_keys <- c(
        iteration_keys,
        paste(level, test_patients, sep = "::")
      )
    }

    train_patients_used <- unique(sample_values[train_indices])
    test_patients_used <- unique(sample_values[test_indices])
    overlapping_patients <- intersect(train_patients_used, test_patients_used)
    if (length(overlapping_patients) > 0L) {
      stop(
        "Patient leakage detected in unseen-patient split iteration ",
        iter_idx, ": ", paste(overlapping_patients, collapse = ", "),
        ". Planned test patients: ",
        paste(patient_plans[[iter_idx]]$test_patients_global, collapse = ", "),
        ". Planned condition-0 train patients: ",
        paste(patient_plans[[iter_idx]]$train_by_condition[[condition_levels[1]]], collapse = ", "),
        ". Actual train patients: ", paste(train_patients_used, collapse = ", "),
        ". Actual test patients: ", paste(test_patients_used, collapse = ", "),
        "."
      )
    }
    test_unit_counts[unit_ids[test_indices]] <-
      test_unit_counts[unit_ids[test_indices]] + 1L
    train_unit_counts[unit_ids[train_indices]] <-
      train_unit_counts[unit_ids[train_indices]] + 1L
    plans[[iter_idx]] <- list(
      unit_col = unit_col,
      train_unit_ids = unit_ids[train_indices],
      test_unit_ids = unit_ids[test_indices],
      train_patient_ids = train_patients_used,
      test_patient_ids = test_patients_used
    )
    test_patient_condition_keys[[iter_idx]] <- unique(iteration_keys)
  }

  unit_coverage_summary <- do.call(rbind, lapply(condition_levels, function(level) {
    idx <- which(condition_values == level)
    ids <- unit_ids[idx]
    data.frame(
      condition = level,
      total_units = length(ids),
      unique_test_units = sum(test_unit_counts[ids] > 0),
      percent_test_units = round(100 * mean(test_unit_counts[ids] > 0), 2),
      min_test_appearances = min(test_unit_counts[ids]),
      median_test_appearances = stats::median(test_unit_counts[ids]),
      max_test_appearances = max(test_unit_counts[ids]),
      unique_train_units = sum(train_unit_counts[ids] > 0),
      percent_train_units = round(100 * mean(train_unit_counts[ids] > 0), 2),
      min_train_appearances = min(train_unit_counts[ids]),
      median_train_appearances = stats::median(train_unit_counts[ids]),
      max_train_appearances = max(train_unit_counts[ids]),
      train_rows_across_iterations = sum(train_unit_counts[ids]),
      test_rows_across_iterations = sum(test_unit_counts[ids]),
      stringsAsFactors = FALSE
    )
  }))

  count_table$test_iterations <- vapply(seq_len(nrow(count_table)), function(i) {
    if (!count_table$eligible[i]) {
      return(0L)
    }
    patient_id <- as.character(count_table$patient_id[i])
    condition <- as.character(count_table$condition[i])
    sum(vapply(patient_plans, function(plan) {
      patient_id %in% plan$test_by_condition[[condition]]
    }, logical(1)))
  }, integer(1))
  count_table$train_iterations <- ifelse(
    count_table$eligible,
    length(random_seeds) - count_table$test_iterations,
    0L
  )

  list(
    success = TRUE,
    reason = NA_character_,
    plans = plans,
    coverage_summary = unit_coverage_summary,
    patient_coverage_summary = count_table,
    test_counts = test_unit_counts,
    train_counts = train_unit_counts,
    test_patient_condition_keys = test_patient_condition_keys,
    patient_train_fraction = patient_train_fraction,
    patient_test_fraction = 1 - patient_train_fraction,
    n_train_patients = n_train_patients,
    n_test_patients = n_test_patients,
    train_rows_per_condition = as.integer(train_rows_per_condition),
    test_rows_per_condition = as.integer(test_rows_per_condition),
    min_units_per_patient = as.integer(min_units_per_patient),
    min_train_units_per_condition = as.integer(min_train_units_per_condition)
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


#' Build feature plans coordinated with held-out patient coverage
#'
#' Prioritizes features with the fewest previous observations among the patient-
#' condition units held out in the current iteration, then uses overall feature
#' inclusion counts and seeded random values to break ties.
#'
#' @param feature_names Character vector of candidate features.
#' @param random_seeds Integer vector of per-iteration seeds.
#' @param test_patient_condition_keys List of held-out patient-condition keys by
#'   iteration.
#' @param eligible_patient_condition_keys Optional complete vector of eligible
#'   patient-condition keys. Supplying this keeps never-held-out patients in the
#'   coverage denominator with zero observations.
#' @param fraction Fraction of candidate features included per iteration.
#'
#' @return Feature plans, standard feature-bagging summary, and patient-feature
#'   coverage table.
#' @export
build_patient_coverage_aware_feature_plan <- function(
    feature_names,
    random_seeds,
    test_patient_condition_keys,
    eligible_patient_condition_keys = NULL,
    fraction = 0.20
) {
  feature_names <- unique(as.character(feature_names))
  feature_names <- feature_names[!is.na(feature_names) & nzchar(feature_names)]
  if (length(feature_names) == 0L) {
    stop("feature_names must contain at least one feature.")
  }
  if (length(test_patient_condition_keys) != length(random_seeds)) {
    stop("test_patient_condition_keys must align with random_seeds.")
  }
  planned_patient_keys <- sort(unique(unlist(
    test_patient_condition_keys,
    use.names = FALSE
  )))
  all_patient_keys <- if (is.null(eligible_patient_condition_keys)) {
    planned_patient_keys
  } else {
    sort(unique(as.character(eligible_patient_condition_keys)))
  }
  all_patient_keys <- all_patient_keys[!is.na(all_patient_keys) & nzchar(all_patient_keys)]
  if (length(all_patient_keys) == 0L) {
    stop("No eligible patient-condition keys were supplied for feature planning.")
  }
  if (!all(planned_patient_keys %in% all_patient_keys)) {
    stop("Every held-out patient-condition key must be included in the eligible key set.")
  }
  if (!is.numeric(fraction) || length(fraction) != 1L ||
      !is.finite(fraction) || fraction <= 0 || fraction > 1) {
    stop("fraction must be one finite numeric value in (0, 1].")
  }

  n_per_iteration <- min(
    length(feature_names),
    max(1L, ceiling(length(feature_names) * fraction))
  )
  inclusion_counts <- stats::setNames(rep(0L, length(feature_names)), feature_names)
  patient_feature_counts <- matrix(
    0L,
    nrow = length(feature_names),
    ncol = length(all_patient_keys),
    dimnames = list(feature_names, all_patient_keys)
  )
  plans <- vector("list", length(random_seeds))

  for (iter_idx in seq_along(random_seeds)) {
    current_keys <- unique(test_patient_condition_keys[[iter_idx]])
    current_keys <- intersect(current_keys, all_patient_keys)
    current_patient_units_covered <- rowSums(
      patient_feature_counts[, current_keys, drop = FALSE] > 0L
    )
    patient_coverage_score <- rowSums(
      patient_feature_counts[, current_keys, drop = FALSE]
    )
    set.seed(as.integer(random_seeds[iter_idx]) + 910000L)
    tie_break <- stats::runif(length(feature_names))
    order_idx <- order(
      current_patient_units_covered,
      patient_coverage_score,
      inclusion_counts[feature_names],
      tie_break
    )
    selected <- feature_names[order_idx[seq_len(n_per_iteration)]]
    inclusion_counts[selected] <- inclusion_counts[selected] + 1L
    patient_feature_counts[selected, current_keys] <-
      patient_feature_counts[selected, current_keys, drop = FALSE] + 1L
    plans[[iter_idx]] <- selected
  }

  patients_covered <- rowSums(patient_feature_counts > 0L)
  patient_feature_coverage <- data.frame(
    feature = feature_names,
    gene_inclusion_count = inclusion_counts[feature_names],
    eligible_patient_condition_units = length(all_patient_keys),
    patient_condition_units_with_shap = patients_covered[feature_names],
    percent_patient_condition_coverage = round(
      100 * patients_covered[feature_names] / length(all_patient_keys),
      2
    ),
    min_patient_feature_test_iterations = apply(patient_feature_counts, 1, min),
    median_patient_feature_test_iterations = apply(patient_feature_counts, 1, stats::median),
    max_patient_feature_test_iterations = apply(patient_feature_counts, 1, max),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  summary <- data.frame(
    total_candidate_features = length(feature_names),
    feature_fraction_per_iteration = fraction,
    features_per_iteration = n_per_iteration,
    min_gene_inclusion_count = min(inclusion_counts),
    median_gene_inclusion_count = stats::median(inclusion_counts),
    max_gene_inclusion_count = max(inclusion_counts),
    unique_features_included = sum(inclusion_counts > 0),
    percent_features_included = round(100 * mean(inclusion_counts > 0), 2),
    min_patient_condition_coverage_per_gene = min(patients_covered),
    median_patient_condition_coverage_per_gene = stats::median(patients_covered),
    max_patient_condition_coverage_per_gene = max(patients_covered),
    stringsAsFactors = FALSE
  )

  list(
    plans = plans,
    summary = summary,
    inclusion_counts = inclusion_counts,
    patient_feature_coverage = patient_feature_coverage,
    patient_feature_counts = patient_feature_counts
  )
}


#' Build balance groups for SHAPTURE sampling and weighting
#'
#' @param dataset Data frame to augment.
#' @param comparison_mode One of `"condition_vs_condition"` or
#'   `"one_vs_others"`.
#' @param sample_col Sample identifier column.
#' @param condition_col Condition column.
#' @param split_level_column Split-defining column used in one-vs-others mode.
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
    split_level_column = NULL,
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

  if (is.null(split_level_column) || !split_level_column %in% colnames(dataset)) {
    stop("split_level_column must be present in dataset for one_vs_others balancing.")
  }

  if (is.null(positive_label)) {
    stop("positive_label must be provided for one_vs_others balancing.")
  }

  sample_ids <- as.character(dataset[[sample_col]])
  split_levels_to_analyze <- as.character(dataset[[split_level_column]])
  class_labels <- as.character(dataset[[condition_col]])

  is_positive <- class_labels == as.character(positive_label)
  dataset[[balance_col]] <- ifelse(
    is_positive,
    sample_ids,
    paste(sample_ids, split_levels_to_analyze, sep = "__")
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


#' Compute within-condition expression consistency
#'
#' Calculates a bounded, condition-balanced expression-consistency score for
#' each feature. Zeros are preserved, upper-tail winsorization is estimated
#' from positive values only, and each condition receives equal weight.
#'
#' @param expression_data Numeric matrix or data frame with rows as cells and
#'   columns as features.
#' @param condition Binary condition vector with one value per row.
#' @param condition_levels Optional length-two vector defining condition 0 and
#'   condition 1. If `NULL`, the two observed values are used in order of first
#'   appearance.
#' @param upper_quantile Positive-value quantile used for upper-tail
#'   winsorization.
#' @param min_positive_cells Minimum number of positive cells required before
#'   winsorization is applied within a condition.
#' @param pseudocount_fraction Fraction of the condition-specific median
#'   positive value used as the pseudocount.
#' @param pseudocount_fallback Minimum pseudocount.
#'
#' @return A data frame with condition-specific diagnostics and the equally
#'   weighted `within_condition_consistency_raw` score for every feature.
#' @keywords internal
compute_within_condition_expression_consistency <- function(
    expression_data,
    condition,
    condition_levels = NULL,
    upper_quantile = 0.99,
    min_positive_cells = 20L,
    pseudocount_fraction = 0.01,
    pseudocount_fallback = 1e-6
) {
  if (!is.data.frame(expression_data) && !is.matrix(expression_data)) {
    stop("expression_data must be a numeric data frame or matrix.")
  }
  if (nrow(expression_data) != length(condition)) {
    stop("condition must contain one value for each row of expression_data.")
  }
  if (nrow(expression_data) < 2 || ncol(expression_data) < 1) {
    stop("expression_data must contain at least two cells and one feature.")
  }
  if (is.null(colnames(expression_data)) || any(!nzchar(colnames(expression_data)))) {
    stop("expression_data must have non-empty feature names.")
  }
  if (anyDuplicated(colnames(expression_data))) {
    stop("expression_data feature names must be unique.")
  }
  numeric_features <- if (is.matrix(expression_data)) {
    rep(is.numeric(expression_data), ncol(expression_data))
  } else {
    vapply(expression_data, is.numeric, logical(1))
  }
  if (any(!numeric_features)) {
    stop(
      "expression_data contains non-numeric feature(s): ",
      paste(head(colnames(expression_data)[!numeric_features], 10), collapse = ", ")
    )
  }
  if (anyNA(condition)) {
    stop("condition cannot contain missing values.")
  }
  if (!is.numeric(upper_quantile) || length(upper_quantile) != 1 ||
      !is.finite(upper_quantile) || upper_quantile <= 0 || upper_quantile >= 1) {
    stop("upper_quantile must be one number strictly between 0 and 1.")
  }
  if (!is.numeric(min_positive_cells) || length(min_positive_cells) != 1 ||
      !is.finite(min_positive_cells) || min_positive_cells < 1 ||
      min_positive_cells != as.integer(min_positive_cells)) {
    stop("min_positive_cells must be a positive integer.")
  }
  if (!is.numeric(pseudocount_fraction) || length(pseudocount_fraction) != 1 ||
      !is.finite(pseudocount_fraction) || pseudocount_fraction <= 0) {
    stop("pseudocount_fraction must be a positive finite number.")
  }
  if (!is.numeric(pseudocount_fallback) || length(pseudocount_fallback) != 1 ||
      !is.finite(pseudocount_fallback) || pseudocount_fallback <= 0) {
    stop("pseudocount_fallback must be a positive finite number.")
  }

  observed_levels <- unique(condition)
  if (is.null(condition_levels)) {
    condition_levels <- observed_levels
  }
  if (length(condition_levels) != 2 || anyNA(condition_levels)) {
    stop("condition_levels must contain exactly two non-missing values.")
  }
  if (anyDuplicated(condition_levels)) {
    stop("condition_levels must contain two distinct values.")
  }
  if (!all(observed_levels %in% condition_levels) ||
      !all(condition_levels %in% observed_levels)) {
    stop("condition and condition_levels must describe exactly the same two conditions.")
  }

  condition_indices <- lapply(condition_levels, function(level) which(condition == level))
  if (any(lengths(condition_indices) < 2)) {
    stop("Each condition must contain at least two cells for consistency calculation.")
  }

  summarize_condition <- function(x, row_index) {
    values <- as.numeric(x[row_index])
    if (any(!is.finite(values))) {
      stop("Expression values must be finite.")
    }
    if (any(values < 0)) {
      stop(
        "Within-condition expression consistency requires non-negative values; ",
        "use unscaled normalized expression rather than centered/scaled values."
      )
    }

    positive_values <- values[values > 0]
    n_positive <- length(positive_values)
    if (n_positive == 0) {
      return(c(
        n_positive = 0,
        all_zero = 1,
        winsorization_used = 0,
        upper_cap = NA_real_,
        pseudocount = NA_real_,
        cv = 0,
        consistency = 1
      ))
    }

    winsorization_used <- n_positive >= min_positive_cells
    upper_cap <- NA_real_
    winsorized_values <- values
    if (winsorization_used) {
      upper_cap <- as.numeric(stats::quantile(
        positive_values,
        probs = upper_quantile,
        names = FALSE,
        na.rm = FALSE,
        type = 7
      ))
      positive_index <- winsorized_values > 0
      winsorized_values[positive_index] <- pmin(
        winsorized_values[positive_index],
        upper_cap
      )
    }

    pseudocount <- max(
      pseudocount_fraction * stats::median(positive_values),
      pseudocount_fallback
    )
    condition_mean <- mean(winsorized_values)
    condition_sd <- stats::sd(winsorized_values)
    condition_cv <- condition_sd / (abs(condition_mean) + pseudocount)
    condition_consistency <- 1 / (1 + condition_cv)

    c(
      n_positive = n_positive,
      all_zero = 0,
      winsorization_used = as.numeric(winsorization_used),
      upper_cap = upper_cap,
      pseudocount = pseudocount,
      cv = condition_cv,
      consistency = condition_consistency
    )
  }

  feature_summaries <- lapply(seq_len(ncol(expression_data)), function(feature_index) {
    feature_values <- if (is.matrix(expression_data)) {
      expression_data[, feature_index]
    } else {
      expression_data[[feature_index]]
    }
    condition0 <- summarize_condition(feature_values, condition_indices[[1]])
    condition1 <- summarize_condition(feature_values, condition_indices[[2]])

    data.frame(
      feature = colnames(expression_data)[feature_index],
      within_condition_positive_n_condition0 = as.integer(condition0["n_positive"]),
      within_condition_positive_n_condition1 = as.integer(condition1["n_positive"]),
      within_condition_all_zero_condition0 = as.logical(condition0["all_zero"]),
      within_condition_all_zero_condition1 = as.logical(condition1["all_zero"]),
      within_condition_winsorization_used_condition0 = as.logical(condition0["winsorization_used"]),
      within_condition_winsorization_used_condition1 = as.logical(condition1["winsorization_used"]),
      within_condition_upper_cap_condition0 = unname(condition0["upper_cap"]),
      within_condition_upper_cap_condition1 = unname(condition1["upper_cap"]),
      within_condition_pseudocount_condition0 = unname(condition0["pseudocount"]),
      within_condition_pseudocount_condition1 = unname(condition1["pseudocount"]),
      within_condition_cv_condition0 = unname(condition0["cv"]),
      within_condition_cv_condition1 = unname(condition1["cv"]),
      within_condition_consistency_condition0 = unname(condition0["consistency"]),
      within_condition_consistency_condition1 = unname(condition1["consistency"]),
      within_condition_consistency_raw = mean(c(
        condition0["consistency"],
        condition1["consistency"]
      )),
      both_conditions_all_zero = as.logical(
        condition0["all_zero"] == 1 && condition1["all_zero"] == 1
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  do.call(rbind, feature_summaries)
}
