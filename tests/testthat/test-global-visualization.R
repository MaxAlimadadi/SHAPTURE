make_global_result_folder <- function(
    root,
    folder_name,
    split_level,
    comparison_mode,
    condition_0,
    condition_1,
    features,
    scores,
    directions
) {
  result_dir <- file.path(root, folder_name)
  docs_dir <- file.path(result_dir, "docs")
  integrated_dir <- file.path(result_dir, "Integrated_adaptive_scoring_results")
  dir.create(docs_dir, recursive = TRUE)
  dir.create(integrated_dir, recursive = TRUE)

  direction_labels <- ifelse(
    directions == "condition_1",
    condition_1,
    condition_0
  )
  combined <- data.frame(
    feature = features,
    final_score = scores,
    rank = seq_along(features),
    Model_expression_pattern = paste(
      "Higher feature values push toward",
      direction_labels
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    combined,
    file.path(docs_dir, "combined_results.csv"),
    row.names = FALSE
  )

  run_summary <- data.frame(
    split_level = split_level,
    comparison_mode = comparison_mode,
    condition_0 = condition_0,
    condition_1 = condition_1,
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    run_summary,
    file.path(integrated_dir, "run_summary.csv"),
    row.names = FALSE
  )

  result_dir
}

test_that("global loading selects the requested top features per folder", {
  root <- tempfile("shapture-global-")
  dir.create(root)

  make_global_result_folder(
    root, "comparison_A", "A", "condition_vs_condition", "Control", "Case",
    c("GeneA", "GeneB", "GeneC"), c(0.9, 0.8, 0.1),
    c("condition_1", "condition_0", "condition_1")
  )
  make_global_result_folder(
    root, "comparison_B", "B", "condition_vs_condition", "Control", "Case",
    c("GeneA", "GeneB", "GeneD"), c(0.85, 0.7, 0.2),
    c("condition_1", "condition_1", "condition_0")
  )

  result <- load_shapture_global_results(root, top_feature_show = 2)

  expect_setequal(rownames(result$score_matrix), c("GeneA", "GeneB"))
  expect_equal(unname(result$n_selected_by_level), c(2L, 2L))
  expect_equal(result$directional_jaccard["A", "B"], 0.5)
  expect_equal(unname(diag(result$directional_jaccard)), c(1, 1))
})

test_that("direction-aware Jaccard distinguishes opposite feature patterns", {
  root <- tempfile("shapture-global-")
  dir.create(root)

  make_global_result_folder(
    root, "comparison_A", "A", "condition_vs_condition", "Control", "Case",
    c("GeneA", "GeneB"), c(0.9, 0.8),
    c("condition_1", "condition_0")
  )
  make_global_result_folder(
    root, "comparison_B", "B", "condition_vs_condition", "Control", "Case",
    c("GeneA", "GeneB"), c(0.9, 0.8),
    c("condition_0", "condition_0")
  )

  result <- load_shapture_global_results(root, top_feature_show = 2)

  expect_equal(result$directional_jaccard["A", "B"], 0.5)
})

test_that("one-vs-others target labels normalize across folders", {
  root <- tempfile("shapture-global-")
  dir.create(root)

  make_global_result_folder(
    root, "one_vs_others_A", "A", "one_vs_others", "others", "A",
    c("GeneA", "GeneB"), c(0.9, 0.8),
    c("condition_1", "condition_0")
  )
  make_global_result_folder(
    root, "one_vs_others_B", "B", "one_vs_others", "others", "B",
    c("GeneA", "GeneB"), c(0.9, 0.8),
    c("condition_1", "condition_0")
  )

  result <- load_shapture_global_results(root, top_feature_show = 2)

  expect_equal(result$comparison_mode, "one_vs_others")
  expect_equal(result$directional_jaccard["A", "B"], 1)
  expect_equal(unname(result$direction_labels), c("Target", "Others"))
})

test_that("global heatmaps use average linkage and preserve matrices", {
  root <- tempfile("shapture-global-")
  dir.create(root)

  make_global_result_folder(
    root, "comparison_A", "A", "condition_vs_condition", "Control", "Case",
    c("GeneA", "GeneB", "GeneC"), c(0.9, 0.8, 0.2),
    c("condition_1", "condition_0", "condition_1")
  )
  make_global_result_folder(
    root, "comparison_B", "B", "condition_vs_condition", "Control", "Case",
    c("GeneA", "GeneB", "GeneD"), c(0.85, 0.7, 0.3),
    c("condition_1", "condition_1", "condition_0")
  )

  result <- load_shapture_global_results(root, top_feature_show = 3)
  landscape <- build_shapture_cross_level_landscape(result)
  similarity <- build_shapture_directional_jaccard_heatmap(result)

  expect_equal(landscape$clustering_method, "average")
  expect_equal(similarity$clustering_method, "average")
  expect_equal(landscape$score_matrix, result$score_matrix)
  expect_equal(similarity$similarity_matrix, result$directional_jaccard)
})
