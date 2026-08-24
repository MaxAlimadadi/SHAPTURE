make_cross_patient_visualization_data <- function(available = TRUE) {
  data.frame(
    feature = c("f1", "f2", "f3"),
    final_score = c(0.85, 0.70, 0.55),
    `average|SHAP value|` = c(1.2, 0.9, 0.6),
    importance_scaled = c(1, 0.6, 0.2),
    expression_score = c(0.8, 0.5, 0.3),
    shap_expr_cor_abs_scaled = c(0.9, 0.7, 0.4),
    model_reproducibility_scaled = c(0.8, 0.7, 0.6),
    stability = c(0.9, 0.8, 0.7),
    expression_log2FC = c(1.0, -0.8, 0.7),
    shap_expr_cor = c(0.8, -0.7, 0.6),
    effect_size = c(0.5, -0.4, 0.3),
    weight_importance = if (available) 0.25 else 0.2777778,
    weight_expression = if (available) 0.25 else 0.2777778,
    weight_value_shap = if (available) 0.30 else 0.3333333,
    weight_model_reproducibility = if (available) 0.10 else 0.1111111,
    weight_cross_patient_support = if (available) 0.10 else 0,
    cross_patient_support = if (available) c(0.9, 0.7, 0.5) else NA_real_,
    cross_patient_support_available = available,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


test_that("Cross-patient visualizations use one consistent 0-1 support scale", {
  df <- make_cross_patient_visualization_data(TRUE)
  palette <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))

  expect_true(cross_patient_support_is_available(df))

  decomposition <- plot_score_decomposition(df, top_n_to_show = 3)
  expect_true("Cross-patient support" %in% as.character(decomposition$data$component))

  component_heatmap <- plot_five_component_heatmap(
    df,
    top_n_to_show = 3,
    continuous_colors = palette
  )
  expect_true("Cross-patient support" %in% as.character(component_heatmap$data$component))
  expect_equal(component_heatmap$scales$get_scales("fill")$limits, c(0, 1))

  evidence_heatmap <- plot_evidence_heatmap(
    df,
    top_n_to_show = 3,
    continuous_colors = palette
  )
  expect_true("Cross-patient support" %in% as.character(evidence_heatmap$data$metric))
  expect_equal(evidence_heatmap$scales$get_scales("fill")$limits, c(0, 1))

  value_shap <- plot_valueSHAP_by_feature(
    df,
    top_n_to_show = 3,
    top_n_to_name = 2,
    what_0 = "control",
    what_1 = "case",
    continuous_colors = palette
  )
  expect_equal(value_shap$labels$fill, "Cross-patient\nsupport")
  expect_equal(value_shap$scales$get_scales("fill")$limits, c(0, 1))
})


test_that("Unavailable Cross-patient support is omitted rather than shown as zero", {
  df <- make_cross_patient_visualization_data(FALSE)
  palette <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  expected_note <- paste0(
    "Cross-patient support unavailable: insufficient patients.\n",
    "Final scores use four renormalized components."
  )

  expect_false(cross_patient_support_is_available(df))

  decomposition <- plot_score_decomposition(df, top_n_to_show = 3)
  expect_false("Cross-patient support" %in% as.character(decomposition$data$component))
  expect_equal(decomposition$labels$caption, expected_note)

  component_heatmap <- plot_five_component_heatmap(
    df,
    top_n_to_show = 3,
    continuous_colors = palette
  )
  expect_false("Cross-patient support" %in% as.character(component_heatmap$data$component))
  expect_match(component_heatmap$labels$title, "Four-component")
  expect_equal(component_heatmap$labels$caption, expected_note)

  evidence_heatmap <- plot_evidence_heatmap(
    df,
    top_n_to_show = 3,
    continuous_colors = palette
  )
  expect_false("Cross-patient support" %in% as.character(evidence_heatmap$data$metric))
  expect_equal(evidence_heatmap$labels$caption, expected_note)

  value_shap <- plot_valueSHAP_by_feature(
    df,
    top_n_to_show = 3,
    top_n_to_name = 2,
    what_0 = "control",
    what_1 = "case",
    continuous_colors = palette
  )
  expect_null(value_shap$labels$fill)
  expect_equal(value_shap$labels$caption, expected_note)
})


test_that("Visualization loading uses combined_results as the support source of truth", {
  result_dir <- file.path(tempdir(), paste0("shapture_viz_", sample.int(1e7, 1)))
  docs_dir <- file.path(result_dir, "docs")
  dir.create(docs_dir, recursive = TRUE)

  df <- make_cross_patient_visualization_data(TRUE)
  utils::write.csv(df, file.path(docs_dir, "combined_results.csv"), row.names = FALSE)
  saveRDS(
    list(
      context = list(
        task_what_0 = "control",
        task_what_1 = "case",
        cross_patient_support_available = FALSE,
        cross_patient_support_unavailable_reason = "insufficient patients"
      ),
      group_colors = c("darkblue", "darkred"),
      continuous_colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
    ),
    file.path(docs_dir, "visualization_inputs.rds")
  )

  expect_warning(
    viz <- load_shapture_visualization_inputs(result_dir),
    "uses combined_results as the source of truth"
  )
  expect_true(viz$cross_patient_support_available)
  expect_true(cross_patient_support_is_available(viz$combined_results))
})
