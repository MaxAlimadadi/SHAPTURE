test_that("safe and full modes preserve the same observed and missing SHAP plot values", {
  skip_on_cran()
  root <- tempfile("memory-mode-regression-")
  dir.create(root)
  original_directory <- getwd()
  on.exit({
    setwd(original_directory)
    unlink(root, recursive = TRUE)
  }, add = TRUE)

  n <- 480L
  condition <- rep(c("left", "right"), each = n / 2L)
  expression <- outer(seq_len(n), seq_len(20L), function(i, j) {
    0.7 + (sin(i * (j + 1)) + 1) / 4
  })
  expression[condition == "right", 1:6] <- expression[condition == "right", 1:6] * 4
  expression[condition == "left", 7:12] <- expression[condition == "left", 7:12] * 4
  colnames(expression) <- paste0("Gene", seq_len(ncol(expression)))
  rownames(expression) <- paste0("Cell", seq_len(n))
  dataset <- data.frame(
    condition = condition, sample_id = rep(paste0("Patient", seq_len(8L)), each = 60L),
    expression, check.names = FALSE
  )
  arguments <- list(
    dataset_all = dataset, what_genes = colnames(expression),
    comparison_mode = "condition_vs_condition", what_0 = "left", what_1 = "right",
    project_name = "memory_regression", iteration_no = 6L,
    gene_fraction_per_iteration = 0.5, cross_patient_bootstraps = 20L,
    analysis_unit_type = "single_cell", model_type = "xgboost", n_cores = 1L,
    analysis_seed = 123L, top_feature_show = 5L, top_feature_label = 5L
  )
  safe_root <- file.path(root, "safe")
  full_root <- file.path(root, "full")
  dir.create(safe_root)
  dir.create(full_root)
  run_comparison <- function(settings) {
    withCallingHandlers(do.call(run_shapture_pipeline, settings), warning = function(warning) {
      # Ignore dependency build-version notices, but expose analytical warnings.
      if (grepl("^package .* was built under R version", conditionMessage(warning))) {
        invokeRestart("muffleWarning")
      }
    })
  }
  safe <- run_comparison(c(arguments, list(main_dir = safe_root)))
  full <- run_comparison(c(arguments, list(main_dir = full_root, memory_mode = "full")))
  safe_dir <- safe$output_dirs[[1]]
  full_dir <- full$output_dirs[[1]]
  safe_results <- utils::read.csv(file.path(safe_dir, "docs", "combined_results.csv"), check.names = FALSE)
  full_results <- utils::read.csv(file.path(full_dir, "docs", "combined_results.csv"), check.names = FALSE)
  expect_gt(nrow(safe_results), 0L)
  expect_equal(safe_results, full_results, tolerance = 1e-12)

  safe_inputs <- readRDS(file.path(safe_dir, "docs", "visualization_inputs.rds"))
  full_inputs <- readRDS(file.path(full_dir, "docs", "visualization_inputs.rds"))
  features <- colnames(safe_inputs$shap$beeswarm_source_S)
  expect_gt(sum(is.na(safe_inputs$shap$beeswarm_source_S)), 0L)
  expect_equal(
    safe_inputs$shap$beeswarm_source_S,
    full_inputs$shap$beeswarm_source_S[, features, drop = FALSE], tolerance = 1e-12
  )
  expect_equal(
    safe_inputs$shap$beeswarm_source_X,
    full_inputs$shap$beeswarm_source_X[, features, drop = FALSE], tolerance = 1e-12
  )
  expect_identical(safe_inputs$roc, full_inputs$roc)
  for (key in names(safe_inputs$shap$condition_aggregate_data)) {
    expect_equal(
      safe_inputs$shap$condition_aggregate_data[[key]],
      full_inputs$shap$condition_aggregate_data[[key]], tolerance = 1e-12
    )
  }
  expect_false(dir.exists(file.path(safe_dir, "docs", ".tmp_iteration_chunks")))
})
