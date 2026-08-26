test_that("Smoother importance retention combines recurrence and positive magnitude", {
  iteration_data <- data.frame(
    iteration = rep(seq_len(5), 2),
    feature = rep(c("feature_a", "feature_b"), each = 5),
    shap_importance = c(0, 0.1, 0.2, 0.3, 0.4, rep(0, 5)),
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
  expect_equal(feature_a$valid_model_iterations, 5L)
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
    shap_importance = rep(0.2, 4),
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


test_that("Model reproducibility follows the approved three-part calculation", {
  iteration_data <- data.frame(
    iteration = c(seq_len(5), seq_len(5)),
    feature = rep(c("feature_a", "feature_b"), each = 5),
    shap_importance = c(rep(0.2, 5), 1, 2, 4, 8, 16),
    delta_shap = c(0.2, 0.3, -0.4, 0.5, 0.6, rep(0.2, 5)),
    relative_delta_shap = c(0.2, 0.3, 0.4, 0.5, 0.6, rep(0.2, 5)),
    value_shap_cor = c(0.8, 0.1, -0.8, NA, 0.5, rep(0.8, 5)),
    stringsAsFactors = FALSE
  )
  full_stats <- data.frame(
    feature = c("feature_a", "feature_b"),
    effect_size = c(0.3, 0.3),
    shap_expr_cor = c(0.7, 0.7),
    stringsAsFactors = FALSE
  )

  result <- SHAPTURE:::compute_model_reproducibility(
    iteration_data = iteration_data,
    full_stats = full_stats,
    shap_expr_cor_cutoff = 0.25,
    lower_probability = 0.10
  )
  feature_a <- result[result$feature == "feature_a", ]

  matching_delta_magnitudes <- c(0.2, 0.3, 0.5, 0.6)
  expected_delta_magnitude <- as.numeric(
    quantile(matching_delta_magnitudes, 0.10, names = FALSE, type = 7)
  ) / median(matching_delta_magnitudes)
  expected_delta <- (4 / 5) * expected_delta_magnitude
  expected_value_shap <- 2 / 4

  expect_equal(feature_a$model_shap_importance_retention, 1)
  expect_equal(feature_a$delta_shap_direction_retention, 4 / 5)
  expect_equal(feature_a$delta_shap_magnitude_retention, expected_delta_magnitude)
  expect_equal(feature_a$delta_shap_retention, expected_delta)
  expect_equal(feature_a$value_shap_retention, expected_value_shap)
  expect_equal(
    feature_a$model_reproducibility,
    mean(c(1, expected_delta, expected_value_shap))
  )
})


test_that("Patient-balanced SHAP effects give patients equal condition weight", {
  shap_matrix <- matrix(
    c(rep(0, 100), 10, 2, 2),
    ncol = 1,
    dimnames = list(NULL, "feature_a")
  )
  metadata <- data.frame(
    sample_id = c(rep("patient_0a", 100), "patient_0b", "patient_1a", "patient_1a"),
    condition = c(rep(0, 101), 1, 1)
  )

  result <- SHAPTURE:::compute_patient_balanced_shap_effect(
    shap_matrix = shap_matrix,
    metadata = metadata
  )

  expect_equal(result$effect_size, 2 - mean(c(0, 10)))
})


test_that("Runtime controls use the current package defaults", {
  config <- SHAPTURE:::default_shapture_script_config(c("feature_a", "feature_b"))

  expect_identical(config$iteration_no, 50)
  expect_identical(config$gene_fraction_per_iteration, 0.30)
  expect_identical(config$cross_patient_bootstraps, 500L)
})
