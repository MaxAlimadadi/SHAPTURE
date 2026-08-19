#' Compute patient-level SHAP summary statistics
#'
#' Collapses observation-level SHAP values to patient-level units and computes
#' descriptive support plus permutation/bootstrap summaries.
#'
#' @param shap_matrix Numeric SHAP matrix with rows aligned to `metadata`.
#' @param metadata Data frame aligned to the SHAP matrix rows.
#' @param sample_col Sample identifier column.
#' @param condition_col Condition column.
#' @param n_permutations Requested number of patient-label permutations.
#' @param n_bootstraps Requested number of bootstrap replicates.
#' @param seed Random seed.
#'
#' @return Data frame with one row per feature.
#' @export
compute_patient_level_shap_stats <- function(
    shap_matrix,
    metadata,
    sample_col = "sample_id",
    condition_col = "condition",
    n_permutations = 1000,
    n_bootstraps = 1000,
    seed = 42
) {
  if (!is.matrix(shap_matrix)) {
    shap_matrix <- as.matrix(shap_matrix)
  }
  
  if (nrow(shap_matrix) != nrow(metadata)) {
    stop("shap_matrix and metadata must have the same number of rows.")
  }
  
  if (!sample_col %in% colnames(metadata)) {
    stop("sample_col is missing from metadata.")
  }
  
  if (!condition_col %in% colnames(metadata)) {
    stop("condition_col is missing from metadata.")
  }
  
  patient_df <- data.frame(
    sample_id = metadata[[sample_col]],
    condition = metadata[[condition_col]],
    stringsAsFactors = FALSE
  )
  
  patient_units <- unique(patient_df)
  patient_units <- patient_units[
    !is.na(patient_units$sample_id) & !is.na(patient_units$condition),
    ,
    drop = FALSE
  ]
  
  if (nrow(patient_units) < 4 || length(unique(patient_units$condition)) != 2) {
    observed <- colMeans(shap_matrix[metadata[[condition_col]] == 1, , drop = FALSE], na.rm = TRUE) -
      colMeans(shap_matrix[metadata[[condition_col]] == 0, , drop = FALSE], na.rm = TRUE)
    
    return(data.frame(
      feature = colnames(shap_matrix),
      effect_size = observed,
      pval = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      bootstrap_support = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  unit_keys <- paste(patient_df$sample_id, patient_df$condition, sep = "___")
  patient_unit_keys <- paste(patient_units$sample_id, patient_units$condition, sep = "___")
  
  patient_means <- t(vapply(patient_unit_keys, function(key) {
    idx <- which(unit_keys == key)
    colMeans(shap_matrix[idx, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(shap_matrix))))
  
  rownames(patient_means) <- patient_unit_keys
  colnames(patient_means) <- colnames(shap_matrix)
  
  unit_conditions <- patient_units$condition
  .compute_patient_level_shap_stats_from_patient_means(
    patient_means = patient_means,
    unit_conditions = unit_conditions,
    n_permutations = n_permutations,
    n_bootstraps = n_bootstraps,
    seed = seed
  )
}


#' Compute patient-level SHAP statistics from patient means
#'
#' @param patient_means Numeric matrix of patient-by-feature means.
#' @param unit_conditions Condition vector aligned to `patient_means`.
#' @param n_permutations Requested number of patient-label permutations.
#' @param n_bootstraps Requested number of bootstrap replicates.
#' @param seed Random seed.
#'
#' @return Data frame with one row per feature.
#' @keywords internal
.compute_patient_level_shap_stats_from_patient_means <- function(
    patient_means,
    unit_conditions,
    n_permutations = 1000,
    n_bootstraps = 1000,
    seed = 42
) {
  if (!is.matrix(patient_means)) {
    patient_means <- as.matrix(patient_means)
  }
  
  if (nrow(patient_means) != length(unit_conditions)) {
    stop("patient_means and unit_conditions must align.")
  }
  
  cond_levels <- sort(unique(unit_conditions))
  idx_0 <- which(unit_conditions == cond_levels[1])
  idx_1 <- which(unit_conditions == cond_levels[2])
  
  observed <- colMeans(patient_means[idx_1, , drop = FALSE], na.rm = TRUE) -
    colMeans(patient_means[idx_0, , drop = FALSE], na.rm = TRUE)
  
  n_distinct_perms <- choose(nrow(patient_means), length(idx_1))
  if (!is.finite(n_distinct_perms) || n_distinct_perms < 1) {
    n_distinct_perms <- 1
  }
  
  n_permutations <- as.integer(min(n_permutations, max(n_distinct_perms * 5L, 99L)))
  n_bootstraps <- as.integer(min(n_bootstraps, 500L))
  
  set.seed(seed)
  perm_null <- matrix(NA_real_, nrow = n_permutations, ncol = ncol(patient_means))
  n_cond1 <- length(idx_1)
  all_idx <- seq_len(nrow(patient_means))
  
  for (b in seq_len(n_permutations)) {
    perm_idx_1 <- sample(all_idx, size = n_cond1, replace = FALSE)
    perm_idx_0 <- setdiff(all_idx, perm_idx_1)
    perm_null[b, ] <- colMeans(patient_means[perm_idx_1, , drop = FALSE], na.rm = TRUE) -
      colMeans(patient_means[perm_idx_0, , drop = FALSE], na.rm = TRUE)
  }
  
  pvals <- (colSums(
    abs(perm_null) >= matrix(abs(observed), nrow = n_permutations, ncol = ncol(patient_means), byrow = TRUE),
    na.rm = TRUE
  ) + 1) / (n_permutations + 1)
  
  set.seed(seed + 1L)
  boot_dist <- matrix(NA_real_, nrow = n_bootstraps, ncol = ncol(patient_means))
  
  for (b in seq_len(n_bootstraps)) {
    boot_idx_0 <- sample(idx_0, size = length(idx_0), replace = TRUE)
    boot_idx_1 <- sample(idx_1, size = length(idx_1), replace = TRUE)
    boot_dist[b, ] <- colMeans(patient_means[boot_idx_1, , drop = FALSE], na.rm = TRUE) -
      colMeans(patient_means[boot_idx_0, , drop = FALSE], na.rm = TRUE)
  }
  
  ci_lower <- apply(boot_dist, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
  ci_upper <- apply(boot_dist, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  bootstrap_support <- ifelse(
    observed >= 0,
    colMeans(boot_dist > 0, na.rm = TRUE),
    colMeans(boot_dist < 0, na.rm = TRUE)
  )
  
  data.frame(
    feature = colnames(patient_means),
    effect_size = observed,
    pval = pvals,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    bootstrap_support = bootstrap_support,
    stringsAsFactors = FALSE
  )
}


#' Compute the minimum achievable permutation p-value
#'
#' @param n_condition0 Number of patients in condition 0.
#' @param n_condition1 Number of patients in condition 1.
#'
#' @return Minimum non-zero p-value under exact patient-label permutations.
#' @export
compute_min_achievable_permutation_p <- function(n_condition0, n_condition1) {
  if (!is.finite(n_condition0) || !is.finite(n_condition1) ||
      n_condition0 < 1 || n_condition1 < 1) {
    return(NA_real_)
  }
  
  n_distinct <- choose(n_condition0 + n_condition1, n_condition1)
  if (!is.finite(n_distinct) || n_distinct < 1) {
    return(NA_real_)
  }
  
  1 / (n_distinct + 1)
}


#' Select the patient-level support component
#'
#' Uses descriptive patient support when patient counts are too small for
#' reliable permutation-based inference, and adjusted permutation p-values once
#' enough patients are available.
#'
#' @param stats_results Data frame containing support and adjusted p-values.
#' @param min_patients Minimum number of patients across the two comparison arms.
#' @param bootstrap_support_col Name of the support column.
#' @param pval_col Name of the adjusted p-value column.
#'
#' @return A list describing the score, label, and interpretation metadata.
#' @export
select_patient_level_support_component <- function(
    stats_results,
    min_patients,
    bootstrap_support_col = "bootstrap_support",
    pval_col = "padj"
) {
  bootstrap_support <- stats_results[[bootstrap_support_col]]
  bootstrap_support[!is.finite(bootstrap_support)] <- 0.5
  
  if (is.na(min_patients) || min_patients < 3) {
    return(list(
      score = bootstrap_support,
      label = "bootstrap_support",
      reliable = FALSE,
      plot_label = "Patient-level support",
      note = paste0(
        "min_patients=", min_patients,
        "; descriptive-only mode; bootstrap_support used as the patient-level support component."
      )
    ))
  }
  
  pval <- stats_results[[pval_col]]
  pval[!is.finite(pval)] <- 1
  score <- -log10(pmax(pval, 1e-10))
  
  list(
    score = score,
    label = "permutation_padj",
    reliable = TRUE,
    plot_label = "Patient-level support",
    note = paste0(
      "min_patients=", min_patients,
      "; permutation-adjusted p-values used as the patient-level support component."
    )
  )
}

# Backward-compatible non-exported alias for older scripts.
select_significance_component <- function(...) {
  select_patient_level_support_component(...)
}
