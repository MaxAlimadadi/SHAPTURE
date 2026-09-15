test_that("automatic visualization sizes follow the documented defaults", {
  expect_equal(
    calculate_shapture_figure_size("feature_importance", 50),
    list(width = 3.5, height = 10, units = "in")
  )
  expect_equal(
    calculate_shapture_figure_size("beeswarm", 50),
    list(width = 16, height = 25, units = "cm")
  )
  expect_equal(
    calculate_shapture_figure_size("density", 50),
    list(width = 17.5, height = 15, units = "in")
  )
  expect_equal(
    calculate_shapture_figure_size("score_decomposition", 50),
    list(width = 6, height = 10, units = "in")
  )
  expect_equal(
    calculate_shapture_figure_size("component_heatmap", 50),
    list(width = 5.7, height = 10, units = "in")
  )
  expect_equal(
    calculate_shapture_figure_size("evidence_heatmap", 50),
    list(width = 4.4, height = 10, units = "in")
  )
  expect_equal(
    calculate_shapture_figure_size("waterfall", 50),
    list(width = 7.4, height = 11, units = "in")
  )
  expect_equal(
    calculate_shapture_figure_size("force", 50),
    list(width = 9.75, height = 4, units = "in")
  )
})

test_that("density layouts stay balanced across feature counts", {
  expected_sizes <- list(
    `10` = list(width = 9.5, height = 7, units = "in"),
    `20` = list(width = 11.5, height = 9, units = "in"),
    `30` = list(width = 13.5, height = 11, units = "in"),
    `50` = list(width = 17.5, height = 15, units = "in")
  )

  for (feature_count in names(expected_sizes)) {
    expect_equal(
      calculate_shapture_figure_size("density", as.integer(feature_count)),
      expected_sizes[[feature_count]]
    )
  }
})

test_that("density facets use the actual available feature count", {
  set.seed(123)
  feature_names <- paste0("feature_", seq_len(10))
  density_source <- data.frame(
    condition = rep(c(0, 1), each = 20),
    matrix(
      stats::rnorm(40 * length(feature_names)),
      nrow = 40,
      dimnames = list(NULL, feature_names)
    ),
    check.names = FALSE
  )
  combined_results <- data.frame(feature = feature_names)

  density_plot <- build_shapture_density_plot(
    density_source = density_source,
    combined_results = combined_results,
    top_feature_show = 50,
    task_what_0 = "condition_0",
    task_what_1 = "condition_1",
    group_color_map = c(condition_0 = "darkblue", condition_1 = "darkred")
  )

  expect_equal(density_plot$facet$params$ncol, 4)
  expect_true(density_plot$facet$params$free$x)
  expect_true(density_plot$facet$params$free$y)
})

test_that("fixed layouts do not grow with the feature count", {
  expect_equal(
    calculate_shapture_figure_size("roc", 500),
    list(width = 7.5, height = 4, units = "in")
  )
  expect_equal(
    calculate_shapture_figure_size("shap_vs_expression", 500),
    list(width = 8, height = 7, units = "in")
  )
  expect_equal(
    calculate_shapture_figure_size("value_shap", 500),
    list(width = 8, height = 8, units = "in")
  )
})

test_that("automatic sizes respect minimum and maximum safeguards", {
  expect_equal(
    calculate_shapture_figure_size("component_heatmap", 5),
    list(width = 5.25, height = 4.5, units = "in")
  )
  expect_equal(
    calculate_shapture_figure_size("component_heatmap", 1000),
    list(width = 8, height = 30, units = "in")
  )
  expect_equal(
    calculate_shapture_figure_size("beeswarm", 1000),
    list(width = 20, height = 60, units = "cm")
  )
})

test_that("users can adjust both global growth rates", {
  expect_equal(
    calculate_shapture_figure_size(
      "component_heatmap",
      50,
      width_per_feature = 0.02,
      height_per_feature = 0.20
    ),
    list(width = 6.2, height = 13.5, units = "in")
  )
})

test_that("automatic size inputs are validated", {
  expect_error(
    calculate_shapture_figure_size("unknown", 10),
    "figure_type"
  )
  expect_error(
    calculate_shapture_figure_size("component_heatmap", -1),
    "n_features"
  )
  expect_error(
    calculate_shapture_figure_size("component_heatmap", 1.5),
    "n_features"
  )
  expect_error(
    calculate_shapture_figure_size("component_heatmap", 10, width_per_feature = NA_real_),
    "width_per_feature"
  )
  expect_error(
    calculate_shapture_figure_size("component_heatmap", 10, height_per_feature = -0.1),
    "height_per_feature"
  )
})

test_that("result-folder discovery is independent of folder naming", {
  analysis_root <- tempfile("shapture-results-")
  dir.create(analysis_root)
  valid_names <- c(
    "condition_A_condition_B_subset_1",
    "one_vs_others_target_2"
  )

  for (folder_name in valid_names) {
    docs_dir <- file.path(analysis_root, folder_name, "docs")
    dir.create(docs_dir, recursive = TRUE)
    file.create(file.path(docs_dir, "combined_results.csv"))
    file.create(file.path(docs_dir, "visualization_inputs.rds"))
  }
  incomplete_dir <- file.path(analysis_root, "incomplete", "docs")
  dir.create(incomplete_dir, recursive = TRUE)
  file.create(file.path(incomplete_dir, "combined_results.csv"))

  discovered <- find_shapture_result_folders(analysis_root)

  expect_equal(basename(discovered), sort(valid_names))
  expect_false(any(grepl("incomplete", discovered, fixed = TRUE)))
})

test_that("result-folder discovery reports invalid roots and empty roots", {
  expect_error(
    find_shapture_result_folders(tempfile("missing-root-")),
    "existing directory"
  )

  empty_root <- tempfile("empty-root-")
  dir.create(empty_root)
  expect_error(
    find_shapture_result_folders(empty_root),
    "No completed SHAPTURE result folders"
  )
})
