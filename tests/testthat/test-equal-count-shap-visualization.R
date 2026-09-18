make_equal_count_visualization_data <- function(n = 1003L) {
  unit_index <- seq_len(n)
  feature_values <- data.frame(
    gene_A = (unit_index %% 101) / 100,
    gene_B = ((unit_index * 7) %% 89) / 88,
    check.names = FALSE
  )
  shap_values <- cbind(
    gene_A = sin(feature_values$gene_A * pi * 2) + unit_index / (n * 20),
    gene_B = feature_values$gene_B - 0.5
  )
  combined_results <- data.frame(
    feature = c("gene_A", "gene_B"),
    final_score = c(0.9, 0.8),
    Model_expression_pattern = c(
      "Higher feature values push toward case",
      "Higher feature values push toward control"
    ),
    stringsAsFactors = FALSE
  )
  list(
    S = shap_values,
    X = feature_values,
    combined_results = combined_results
  )
}


test_that("Joint SHAP-feature grouping is deterministic and balanced", {
  inputs <- make_equal_count_visualization_data()
  groups <- build_shapture_equal_count_groups(
    inputs$S,
    inputs$X,
    features = c("gene_A", "gene_B"),
    target_groups = 50,
    min_units_per_group = 10,
    raw_unit_threshold = 500
  )
  repeated <- build_shapture_equal_count_groups(
    inputs$S,
    inputs$X,
    features = c("gene_A", "gene_B"),
    target_groups = 50,
    min_units_per_group = 10,
    raw_unit_threshold = 500
  )

  expect_identical(groups, repeated)
  expect_identical(attr(groups, "displayed_features"), c("gene_A", "gene_B"))
  expect_equal(as.integer(table(groups$feature)), c(50, 50))
  expect_equal(
    as.numeric(tapply(groups$n_units, groups$feature, sum)),
    c(1003, 1003)
  )
  expect_lte(max(groups$n_units) - min(groups$n_units), 1)
  expect_true(all(groups$grouping_mode == "joint_SHAP_feature_equal_count"))
})


test_that("Small datasets retain their original units", {
  inputs <- make_equal_count_visualization_data(40)
  groups <- build_shapture_equal_count_groups(
    inputs$S,
    inputs$X,
    features = "gene_A",
    target_groups = 500,
    min_units_per_group = 10,
    raw_unit_threshold = 500
  )

  expect_equal(nrow(groups), 40)
  expect_true(all(groups$n_units == 1L))
  expect_true(all(groups$grouping_mode == "original_units"))
})


test_that("Binned beeswarm and dependency plots share clear labels", {
  inputs <- make_equal_count_visualization_data()
  groups <- build_shapture_equal_count_groups(
    inputs$S,
    inputs$X,
    features = c("gene_A", "gene_B"),
    target_groups = 50,
    min_units_per_group = 10,
    raw_unit_threshold = 500
  )
  binned_beeswarm <- build_shapture_binned_beeswarm_plot(
    combined_results = inputs$combined_results,
    grouped_data = groups,
    task_what_0 = "control",
    task_what_1 = "case",
    group_colors = c("darkblue", "darkred"),
    continuous_colors = c("blue", "white", "red")
  )
  dependency <- build_shapture_dependency_plot(
    grouped_data = groups,
    feature = "gene_A",
    beeswarm_source_S = inputs$S,
    beeswarm_source_X = inputs$X,
    task_what_0 = "control",
    task_what_1 = "case",
    continuous_colors = c("blue", "white", "red")
  )
  expected_subtitle <- paste0(
    "Each hexagon summarizes an equal-count group of cells formed jointly ",
    "from SHAP and feature values."
  )

  expect_s3_class(binned_beeswarm, "ggplot")
  expect_s3_class(dependency, "ggplot")
  expect_equal(attr(binned_beeswarm, "displayed_feature_count"), 2)
  expect_equal(attr(dependency, "feature"), "gene_A")
  expect_equal(attr(dependency, "n_groups"), 50)
  expect_equal(binned_beeswarm$labels$title, "Left = control | Right = case")
  expect_equal(dependency$labels$title, "Left = control | Right = case")
  expect_equal(binned_beeswarm$labels$subtitle, expected_subtitle)
  expect_equal(dependency$labels$subtitle, expected_subtitle)
  expect_equal(binned_beeswarm$labels$colour, "Mean feature value")
  expect_equal(dependency$labels$colour, "Mean feature value")
  expect_equal(dependency$labels$y, "Scaled gene_A value")

  expect_s3_class(binned_beeswarm$layers[[1]]$geom, "GeomShaptureHexPoint")
  expect_s3_class(dependency$layers[[1]]$geom, "GeomShaptureHexPoint")
  beeswarm_hexagons <- ggplot2::layer_grob(binned_beeswarm, 1L)[[1]]
  dependency_hexagons <- ggplot2::layer_grob(dependency, 1L)[[1]]
  expect_s3_class(beeswarm_hexagons, "polygon")
  expect_s3_class(dependency_hexagons, "polygon")
  expect_true(all(!is.na(beeswarm_hexagons$gp$fill)))
  expect_true(all(!is.na(dependency_hexagons$gp$fill)))
  expect_equal(length(beeswarm_hexagons$id), nrow(groups) * 6L)
  expect_equal(length(dependency_hexagons$id), sum(groups$feature == "gene_A") * 6L)
})


test_that("Equal-count grouping validates inputs and reports missing features", {
  inputs <- make_equal_count_visualization_data(100)
  expect_error(
    build_shapture_equal_count_groups(inputs$S[-1, ], inputs$X, "gene_A"),
    "same number of rows"
  )
  expect_warning(
    groups <- build_shapture_equal_count_groups(
      inputs$S,
      inputs$X,
      c("gene_A", "missing_gene")
    ),
    "missing_gene"
  )
  expect_identical(attr(groups, "displayed_features"), "gene_A")
  expect_error(
    build_shapture_dependency_plot(
      groups,
      "missing_gene",
      inputs$S,
      inputs$X,
      "control",
      "case",
      c("blue", "red")
    ),
    "does not contain feature"
  )
})
