test_that("homogeneous absence is treated as perfect within-condition consistency", {
  expression_data <- data.frame(
    absent_in_condition0 = c(rep(0, 25), rep(2, 25)),
    absent_in_both = rep(0, 50),
    check.names = FALSE
  )
  condition <- rep(c("condition0", "condition1"), each = 25)

  result <- compute_within_condition_expression_consistency(
    expression_data = expression_data,
    condition = condition,
    condition_levels = c("condition0", "condition1")
  )

  absent0 <- result[result$feature == "absent_in_condition0", ]
  expect_true(absent0$within_condition_all_zero_condition0)
  expect_equal(absent0$within_condition_consistency_condition0, 1)
  expect_equal(absent0$within_condition_consistency_condition1, 1)
  expect_equal(absent0$within_condition_consistency_raw, 1)
  expect_false(absent0$both_conditions_all_zero)

  absent_both <- result[result$feature == "absent_in_both", ]
  expect_true(absent_both$both_conditions_all_zero)
})


test_that("upper-tail winsorization preserves zeros and requires 20 positive cells", {
  expression_data <- data.frame(
    eligible = c(rep(0, 4), rep(1, 20), 100, rep(1, 25)),
    sparse = c(rep(0, 21), 1, 2, 3, 4, rep(1, 25)),
    check.names = FALSE
  )
  condition <- rep(c(0, 1), each = 25)

  result <- compute_within_condition_expression_consistency(
    expression_data = expression_data,
    condition = condition,
    condition_levels = c(0, 1)
  )

  eligible <- result[result$feature == "eligible", ]
  expect_true(eligible$within_condition_winsorization_used_condition0)
  expect_equal(eligible$within_condition_positive_n_condition0, 21)
  expect_equal(eligible$within_condition_upper_cap_condition0, 80.2, tolerance = 1e-12)

  sparse <- result[result$feature == "sparse", ]
  expect_false(sparse$within_condition_winsorization_used_condition0)
  expect_equal(sparse$within_condition_positive_n_condition0, 4)
  expect_true(is.na(sparse$within_condition_upper_cap_condition0))
  expect_equal(sparse$within_condition_pseudocount_condition0, 0.025)
})


test_that("condition consistency is equally weighted despite unequal cell counts", {
  expression_data <- data.frame(
    feature_a = c(1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1),
    check.names = FALSE
  )
  condition <- c(rep("condition0", 3), rep("condition1", 9))

  result <- compute_within_condition_expression_consistency(
    expression_data = expression_data,
    condition = condition,
    condition_levels = c("condition0", "condition1")
  )

  expected <- mean(c(
    result$within_condition_consistency_condition0,
    result$within_condition_consistency_condition1
  ))
  expect_equal(result$within_condition_consistency_raw, expected)
})


test_that("within-condition consistency rejects values outside its expression domain", {
  expect_error(
    compute_within_condition_expression_consistency(
      expression_data = data.frame(feature_a = c(0, 1, -1, 2)),
      condition = c(0, 0, 1, 1),
      condition_levels = c(0, 1)
    ),
    "non-negative values"
  )

  expect_error(
    compute_within_condition_expression_consistency(
      expression_data = data.frame(feature_a = c(0, 1, NA, 2)),
      condition = c(0, 0, 1, 1),
      condition_levels = c(0, 1)
    ),
    "finite"
  )
})


test_that("Expression signal uses the new consistency score", {
  average_importance <- data.frame(
    feature = c("feature_a", "feature_b", "feature_c"),
    `average|SHAP value|` = c(1, 2, 3),
    check.names = FALSE
  )
  stats_results <- data.frame(
    feature = c("feature_a", "feature_b", "feature_c"),
    effect_size = c(0.2, 0.3, 0.4),
    cross_patient_support = c(0.7, 0.8, 0.9),
    cross_patient_support_available = TRUE,
    stability = c(0.7, 0.8, 0.9),
    min_patients = rep(3, 3),
    expression_log2FC = c(1, 2, 3),
    det_diff = c(0.1, 0.2, 0.3),
    shap_expr_cor = c(0.4, 0.6, 0.8),
    within_condition_consistency_raw = c(0.2, 0.5, 0.8),
    check.names = FALSE
  )

  result <- compute_combined_shap_ranking(
    average_importance = average_importance,
    stats_results = stats_results,
    default_prior_scores_weights = c(
      importance = 0,
      cross_patient_support = 0,
      model_reproducibility = 0,
      expression = 1,
      value_shap = 0
    ),
    score_scaling_method = "minmax"
  )
  result <- result[match(stats_results$feature, result$feature), ]

  expect_equal(result$within_condition_consistency_scaled, c(0, 0.5, 1))
  expect_equal(result$expression_score, c(0, 0.5, 1))
  expect_false(any(c(
    "expr_variance",
    "within_condition_variance_inv",
    "within_condition_variance_scaled",
    "var_inv",
    "var_scaled"
  ) %in% colnames(result)))
})
