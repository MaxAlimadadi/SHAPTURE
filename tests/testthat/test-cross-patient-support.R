make_cross_patient_test_data <- function(patients_per_condition = 3L) {
  patient_ids <- c(
    paste0("c0_p", seq_len(patients_per_condition)),
    paste0("c1_p", seq_len(patients_per_condition))
  )
  conditions <- rep(c(0, 1), each = patients_per_condition)
  metadata <- do.call(rbind, lapply(seq_along(patient_ids), function(i) {
    data.frame(
      sample_id = patient_ids[i],
      condition = conditions[i],
      cell_index = seq_len(6),
      stringsAsFactors = FALSE
    )
  }))
  row.names(metadata) <- NULL

  set.seed(101)
  f1 <- pmax(0, metadata$condition * 2 + stats::rnorm(nrow(metadata), sd = 0.2))
  f2 <- pmax(0, (1 - metadata$condition) * 1.5 + stats::rnorm(nrow(metadata), sd = 0.25))
  expression <- data.frame(f1 = f1, f2 = f2)
  shap <- cbind(
    f1 = scale(f1)[, 1] + stats::rnorm(length(f1), sd = 0.05),
    f2 = scale(f2)[, 1] + stats::rnorm(length(f2), sd = 0.05)
  )
  model_data <- cbind(metadata[, c("sample_id", "condition")], expression)
  consistency <- compute_within_condition_expression_consistency(
    expression_data = expression,
    condition = metadata$condition,
    condition_levels = c(0, 1)
  )
  means0 <- colMeans(expression[metadata$condition == 0, , drop = FALSE])
  means1 <- colMeans(expression[metadata$condition == 1, , drop = FALSE])
  pseudocount <- compute_dynamic_pseudocount(as.matrix(expression))
  log2fc <- log2((means1 + pseudocount) / (means0 + pseudocount))
  det0 <- colMeans(expression[metadata$condition == 0, , drop = FALSE] > 0)
  det1 <- colMeans(expression[metadata$condition == 1, , drop = FALSE] > 0)
  shap_cor <- vapply(colnames(expression), function(feature) {
    stats::cor(expression[[feature]], shap[, feature], method = "spearman")
  }, numeric(1))

  full_stats <- data.frame(
    feature = colnames(expression),
    expression_log2FC = log2fc,
    det_diff = det1 - det0,
    shap_expr_cor = shap_cor,
    within_condition_consistency_raw = consistency$within_condition_consistency_raw,
    check.names = FALSE
  )
  average_importance <- data.frame(
    feature = colnames(shap),
    `average|SHAP value|` = colMeans(abs(shap)),
    check.names = FALSE
  )

  list(
    shap = shap,
    expression = expression,
    metadata = metadata,
    model_data = model_data,
    single_cell_data = model_data,
    full_stats = full_stats,
    average_importance = average_importance
  )
}


test_that("Cross-patient support is bounded and reproducible", {
  test_data <- make_cross_patient_test_data(3L)
  result1 <- compute_cross_patient_support(
    shap_matrix = test_data$shap,
    shap_feature_values = test_data$expression,
    shap_metadata = test_data$metadata,
    model_expression_data = test_data$model_data,
    single_cell_expression_data = test_data$single_cell_data,
    full_stats = test_data$full_stats,
    average_importance = test_data$average_importance,
    feature_names = colnames(test_data$shap),
    comparison_mode = "condition_vs_condition",
    n_bootstraps = 100L,
    seed = 44L
  )
  result2 <- compute_cross_patient_support(
    shap_matrix = test_data$shap,
    shap_feature_values = test_data$expression,
    shap_metadata = test_data$metadata,
    model_expression_data = test_data$model_data,
    single_cell_expression_data = test_data$single_cell_data,
    full_stats = test_data$full_stats,
    average_importance = test_data$average_importance,
    feature_names = colnames(test_data$shap),
    comparison_mode = "condition_vs_condition",
    n_bootstraps = 100L,
    seed = 44L
  )

  expect_true(result1$available)
  expect_equal(result1$n_bootstraps_valid_task, 100L)
  expect_true(all(result1$support$cross_patient_support >= 0))
  expect_true(all(result1$support$cross_patient_support <= 1))
  expect_equal(result1$support, result2$support)
  expect_equal(result1$diagnostics, result2$diagnostics)
})


test_that("Two patients per condition are supported without a low-patient warning", {
  test_data <- make_cross_patient_test_data(2L)

  expect_no_warning(
    result <- compute_cross_patient_support(
      shap_matrix = test_data$shap,
      shap_feature_values = test_data$expression,
      shap_metadata = test_data$metadata,
      model_expression_data = test_data$model_data,
      single_cell_expression_data = test_data$single_cell_data,
      full_stats = test_data$full_stats,
      average_importance = test_data$average_importance,
      feature_names = colnames(test_data$shap),
      comparison_mode = "condition_vs_condition",
      n_bootstraps = 25L,
      seed = 45L
    )
  )

  expect_true(result$available)
  expect_equal(result$n_patients_condition0, 2L)
  expect_equal(result$n_patients_condition1, 2L)
})


test_that("one-vs-others resamples patients globally", {
  test_data <- make_cross_patient_test_data(3L)
  shared_patient_ids <- rep(paste0("p", seq_len(3)), each = 6)
  test_data$metadata$sample_id <- rep(shared_patient_ids, 2)
  test_data$model_data$sample_id <- test_data$metadata$sample_id
  test_data$single_cell_data$sample_id <- test_data$metadata$sample_id

  result <- compute_cross_patient_support(
    shap_matrix = test_data$shap,
    shap_feature_values = test_data$expression,
    shap_metadata = test_data$metadata,
    model_expression_data = test_data$model_data,
    single_cell_expression_data = test_data$single_cell_data,
    full_stats = test_data$full_stats,
    average_importance = test_data$average_importance,
    feature_names = colnames(test_data$shap),
    comparison_mode = "one_vs_others",
    n_bootstraps = 50L,
    seed = 55L
  )

  expect_true(result$available)
  expect_equal(result$n_patients_total, 3L)
  expect_equal(result$n_bootstraps_valid_task, 50L)
  expect_true(all(result$support$cross_patient_support >= 0))
  expect_true(all(result$support$cross_patient_support <= 1))
})


test_that("Unavailable support is missing and remaining weights are renormalized", {
  stats_results <- data.frame(
    feature = c("f1", "f2"),
    effect_size = c(0.4, 0.3),
    stability = c(0.8, 0.9),
    min_patients = c(1, 1),
    expression_log2FC = c(1, -1),
    det_diff = c(0.4, -0.4),
    shap_expr_cor = c(0.8, -0.8),
    within_condition_consistency_raw = c(0.7, 0.8),
    cross_patient_support = NA_real_,
    cross_patient_support_available = FALSE,
    cross_patient_support_unavailable_reason = "insufficient patients",
    check.names = FALSE
  )
  average_importance <- data.frame(
    feature = c("f1", "f2"),
    `average|SHAP value|` = c(0.5, 0.4),
    check.names = FALSE
  )
  ranked <- compute_combined_shap_ranking(
    average_importance = average_importance,
    stats_results = stats_results,
    what_0 = "left",
    what_1 = "right",
    score_scaling_method = "minmax"
  )

  expect_true(all(is.na(ranked$cross_patient_support)))
  expect_equal(unique(ranked$weight_cross_patient_support), 0)
  weight_sum <- with(
    ranked[1, ],
    weight_importance + weight_cross_patient_support +
      weight_model_reproducibility + weight_expression + weight_value_shap
  )
  expect_equal(weight_sum, 1)
})
