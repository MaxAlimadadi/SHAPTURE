test_that("Smoother importance retention combines recurrence and positive magnitude", {
  importance_values <- c(rep(0, 75), rep(0.45, 23), rep(0.60, 402))

  result <- compute_smoother_importance_retention(importance_values)

  expect_equal(unname(result["nonzero_frequency"]), 425 / 500)
  expect_equal(unname(result["positive_median"]), 0.60)
  expect_equal(unname(result["positive_lower"]), 0.45)
  expect_equal(unname(result["positive_retention"]), 0.45 / 0.60)
  expect_equal(unname(result["support"]), 0.85 * 0.75)
  expect_equal(unname(result["n_valid"]), 500)
  expect_equal(unname(result["n_positive"]), 425)
})


test_that("Smoother importance retention handles absent or invalid importance", {
  all_zero <- compute_smoother_importance_retention(rep(0, 50))
  no_valid <- compute_smoother_importance_retention(c(NA_real_, Inf, -Inf))

  expect_equal(unname(all_zero["support"]), 0)
  expect_equal(unname(all_zero["nonzero_frequency"]), 0)
  expect_equal(unname(all_zero["n_valid"]), 50)
  expect_equal(unname(no_valid["support"]), 0)
  expect_equal(unname(no_valid["n_valid"]), 0)
})


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
  expect_equal(
    result1$diagnostics$shap_importance_support,
    result1$diagnostics$shap_importance_nonzero_frequency *
      result1$diagnostics$shap_importance_positive_retention
  )
  expect_true(all(
    result1$diagnostics$shap_importance_positive_bootstraps <=
      result1$diagnostics$shap_importance_valid_bootstraps
  ))
})


test_that("Cross-patient SHAP-importance support is invariant to SHAP units", {
  test_data <- make_cross_patient_test_data(3L)
  scaled_test_data <- test_data
  scaled_test_data$shap <- test_data$shap * 10
  scaled_test_data$average_importance[["average|SHAP value|"]] <-
    test_data$average_importance[["average|SHAP value|"]] * 10

  base_result <- compute_cross_patient_support(
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
    seed = 64L
  )
  scaled_result <- compute_cross_patient_support(
    shap_matrix = scaled_test_data$shap,
    shap_feature_values = scaled_test_data$expression,
    shap_metadata = scaled_test_data$metadata,
    model_expression_data = scaled_test_data$model_data,
    single_cell_expression_data = scaled_test_data$single_cell_data,
    full_stats = scaled_test_data$full_stats,
    average_importance = scaled_test_data$average_importance,
    feature_names = colnames(scaled_test_data$shap),
    comparison_mode = "condition_vs_condition",
    n_bootstraps = 100L,
    seed = 64L
  )

  expect_equal(
    base_result$diagnostics$shap_importance_support,
    scaled_result$diagnostics$shap_importance_support
  )
  expect_equal(
    scaled_result$diagnostics$shap_importance_positive_median,
    base_result$diagnostics$shap_importance_positive_median * 10
  )
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


test_that("Feature-specific Cross-patient Support contributes without recycling", {
  stats_results <- data.frame(
    feature = c("A1BG", "feature_b", "feature_c"),
    effect_size = c(0.4, 0.4, 0.4),
    stability = c(0.5, 0.5, 0.5),
    model_reproducibility = c(0.5, 0.5, 0.5),
    min_patients = c(3, 3, 3),
    expression_log2FC = c(1, 1, 1),
    det_diff = c(0.4, 0.4, 0.4),
    shap_expr_cor = c(0.8, 0.8, 0.8),
    within_condition_consistency_raw = c(0.7, 0.7, 0.7),
    cross_patient_support = c(0.2, 0.8, 0.5),
    cross_patient_support_available = TRUE,
    cross_patient_support_unavailable_reason = NA_character_,
    check.names = FALSE
  )
  average_importance <- data.frame(
    feature = stats_results$feature,
    `average|SHAP value|` = c(0.5, 0.5, 0.5),
    check.names = FALSE
  )

  ranked <- compute_combined_shap_ranking(
    average_importance = average_importance,
    stats_results = stats_results,
    default_prior_scores_weights = c(
      importance = 0,
      cross_patient_support = 1,
      model_reproducibility = 0,
      expression = 0,
      value_shap = 0
    ),
    what_0 = "left",
    what_1 = "right",
    score_scaling_method = "minmax"
  )
  ranked <- ranked[match(stats_results$feature, ranked$feature), ]

  expect_equal(ranked$final_score, stats_results$cross_patient_support)
  expect_equal(ranked$rank, c(3L, 1L, 2L))
})
