test_that("Positive median anchors exclude zero and non-finite references", {
  expect_equal(compute_positive_median_anchor(c(0, 0.2, 0.6, NA, Inf)), 0.4)
  expect_true(is.na(compute_positive_median_anchor(c(0, NA, Inf))))
  expect_equal(compute_positive_median_anchor(c(1e-10, 0.2, 0.6)), 0.4)
  expect_error(compute_positive_median_anchor(c(-0.1, 0.5)), "nonnegative")
})

test_that("Median-anchor transforms preserve zero and do not stretch finite maxima", {
  expect_equal(
    apply_median_anchor_transform(c(0, 0.2, 0.4, 1.2, NA, Inf), anchor = 0.4),
    c(0, 1 / 3, 0.5, 0.75, NA_real_, NA_real_)
  )
  expect_equal(apply_median_anchor_transform(rep(0.8, 3), 0.8), rep(0.5, 3))
  expect_equal(apply_median_anchor_transform(c(0, NA), NA_real_), c(0, NA_real_))
  expect_error(apply_median_anchor_transform(0.1, NA_real_), "positive median anchor")
  expect_error(apply_median_anchor_transform(-0.1, 0.4), "nonnegative")
})

test_that("Naturally bounded evidence retains absolute support levels", {
  importance <- data.frame(
    feature = c("a", "b", "c"),
    `average|SHAP value|` = c(0.2, 0.4, 0.8),
    check.names = FALSE
  )
  stats_results <- data.frame(
    feature = importance$feature,
    effect_size = c(0.2, 0.3, 0.4),
    expression_log2FC = c(0.5, 1, 2),
    det_diff = c(0.7, 0.8, 0.9),
    within_condition_consistency_raw = c(0.7, 0.8, 0.9),
    shap_expr_cor = c(0.7, 0.8, 0.9),
    cross_patient_support = c(0.7, 0.8, 0.9),
    cross_patient_support_available = TRUE,
    model_reproducibility = rep(1, 3)
  )
  result <- compute_combined_shap_ranking(importance, stats_results)
  result <- result[match(importance$feature, result$feature), ]

  expect_equal(result$importance_scaled, c(1 / 3, 0.5, 2 / 3))
  expect_equal(result$expression_log2FC_abs_scaled, c(1 / 3, 0.5, 2 / 3))
  expect_equal(result$det_diff_scaled, c(0.7, 0.8, 0.9))
  expect_equal(result$within_condition_consistency_scaled, c(0.7, 0.8, 0.9))
  expect_equal(result$shap_expr_cor_abs_scaled, c(0.7, 0.8, 0.9))
  expect_equal(result$cross_patient_support, c(0.7, 0.8, 0.9))
  expect_equal(result$model_reproducibility_scaled, rep(1, 3))

  legacy <- compute_combined_shap_ranking(
    importance, stats_results, score_scaling_method = "winsorized"
  )
  expect_equal(legacy[order(legacy$feature), ], result[order(result$feature), ])
})
