test_that("Smoother importance retention combines recurrence and positive magnitude", {
  iteration_data <- data.frame(
    iteration = rep(seq_len(5), 2),
    feature = rep(c("feature_a", "feature_b"), each = 5),
    relative_importance = c(0, 0.1, 0.2, 0.3, 0.4, rep(0, 5)),
    delta_shap = c(rep(0.2, 5), rep(-0.2, 5)),
    relative_delta_shap = rep(0.2, 10),
    value_shap_cor = c(rep(0.8, 5), rep(-0.8, 5)),
    stringsAsFactors = FALSE
  )
  full_stats <- data.frame(
    feature = c("feature_a", "feature_b"),
    effect_size = c(0.3, -0.3),
    shap_expr_cor = c(0.7, -0.7),
    stringsAsFactors = FALSE
  )

  result <- SHAPTURE:::compute_model_reproducibility(
    iteration_data = iteration_data,
    full_stats = full_stats,
    shap_expr_cor_cutoff = 0.25,
    lower_probability = 0.10
  )

  feature_a <- result[result$feature == "feature_a", ]
  positive_values <- c(0.1, 0.2, 0.3, 0.4)
  expected_positive_retention <- as.numeric(
    quantile(positive_values, 0.10, names = FALSE, type = 7)
  ) / median(positive_values)
  expect_equal(feature_a$model_shap_importance_nonzero_frequency, 4 / 5)
  expect_equal(
    feature_a$model_shap_importance_retention,
    (4 / 5) * expected_positive_retention
  )
  expect_equal(feature_a$delta_shap_retention, 1)
  expect_equal(feature_a$value_shap_retention, 1)

  feature_b <- result[result$feature == "feature_b", ]
  expect_equal(feature_b$model_shap_importance_retention, 0)
  expect_equal(feature_b$model_reproducibility, 2 / 3)
  expect_true(all(result$model_reproducibility >= 0 & result$model_reproducibility <= 1))
})


test_that("Value-SHAP retention uses valid iteration correlations as denominator", {
  iteration_data <- data.frame(
    iteration = seq_len(4),
    feature = "feature_a",
    relative_importance = rep(0.2, 4),
    delta_shap = rep(0.3, 4),
    relative_delta_shap = rep(0.2, 4),
    value_shap_cor = c(0.8, NA, 0.1, 0.7),
    stringsAsFactors = FALSE
  )
  full_stats <- data.frame(
    feature = "feature_a",
    effect_size = 0.3,
    shap_expr_cor = 0.7,
    stringsAsFactors = FALSE
  )

  result <- SHAPTURE:::compute_model_reproducibility(
    iteration_data = iteration_data,
    full_stats = full_stats,
    shap_expr_cor_cutoff = 0.25
  )

  expect_equal(result$valid_value_shap_iterations, 3L)
  expect_equal(result$value_shap_supporting_iterations, 2L)
  expect_equal(result$value_shap_retention, 2 / 3)
})


test_that("Runtime controls retain backward-compatible defaults", {
  config <- SHAPTURE:::default_shapture_script_config(c("feature_a", "feature_b"))

  expect_identical(config$gene_fraction_per_iteration, 0.20)
  expect_identical(config$cross_patient_bootstraps, 1000L)
})
