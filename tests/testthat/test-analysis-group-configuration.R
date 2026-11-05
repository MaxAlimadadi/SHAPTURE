make_group_config_data <- function() {
  data.frame(
    condition = rep(c("control", "case"), 3),
    sample_id = paste0("patient_", seq_len(6)),
    grouping = rep(c("A", "B", "C"), each = 2),
    feature = seq_len(6),
    stringsAsFactors = FALSE
  )
}

test_that("analysis-group configuration uses only the new names", {
  config <- default_shapture_script_config("feature")
  expect_true(all(c("split_level_column", "split_levels_to_analyze") %in% names(config)))
  expect_null(config$split_level_column)
  expect_null(config$split_levels_to_analyze)
  removed_names <- c("subset_levels", "split_by_col", "split_levels")
  expect_false(any(removed_names %in% names(config)))
  expect_false(any(removed_names %in% names(formals(run_shapture_pipeline))))
  expect_true("split_level_column" %in% names(formals(build_dataset_from_seurat)))
  expect_true("split_level_column" %in% names(formals(build_dataset_from_sce)))
})

test_that("removed config names fail clearly, including NULL values", {
  for (name in c("subset_levels", "split_by_col", "split_levels")) {
    env <- list2env(list(shapture_config = setNames(list(NULL), name)))
    expect_error(
      resolve_shapture_script_config(env, make_group_config_data(), "feature", tempdir()),
      "Unknown shapture_config parameters"
    )
    expect_error(
      do.call(run_shapture_pipeline, setNames(list(NULL), name)),
      "Unused arguments"
    )
  }
})

test_that("script resolution selects tasks without filtering input rows", {
  input <- make_group_config_data()
  for (mode in c("condition_vs_condition", "one_vs_others")) {
    env <- list2env(list(shapture_config = list(
      split_level_column = "grouping",
      split_levels_to_analyze = "A",
      comparison_mode = mode,
      what_0 = "control", what_1 = "case"
    )))
    config <- resolve_shapture_script_config(env, input, "feature", tempdir())
    expect_identical(config$split_level_column, "grouping")
    expect_identical(config$split_levels_to_analyze, "A")
    expect_identical(row.names(config$dataset_all), row.names(input))
    expect_identical(config$dataset_all$grouping, input$grouping)
    if (mode == "one_vs_others") {
      expect_identical(config$one_vs_others_col, "grouping")
      expect_identical(config$one_vs_others_levels, "A")
    }
  }
})

test_that("task selection and manual input filtering have distinct roles", {
  input <- make_group_config_data()
  make_tasks <- function(data, mode, levels = NULL, column = "grouping") {
    build_shapture_tasks(
      dataset_all = data, comparison_mode = mode,
      split_level_column = column, split_levels_to_analyze = levels,
      one_vs_others_col = NULL, one_vs_others_levels = NULL,
      what_0 = "control", what_1 = "case", one_vs_others_negative_label = "others"
    )
  }
  conditions <- make_tasks(input, "condition_vs_condition", "A")
  expect_identical(conditions$task_table$temp_split, "A")
  expect_identical(conditions$dataset, input)
  expect_identical(make_tasks(input, "condition_vs_condition")$task_table$temp_split, c("A", "B", "C"))
  expect_identical(make_tasks(input, "condition_vs_condition", column = NULL)$task_table$temp_split, "all_data")

  others <- make_tasks(input, "one_vs_others", "A")
  expect_identical(others$task_table$temp_split, "A")
  task <- materialize_task_dataset(others$dataset, "one_vs_others", "A", NULL, "grouping", "others", "A")
  expect_identical(task$grouping[task$condition == "others"], c("B", "B", "C", "C"))
  filtered <- input[input$grouping %in% c("A", "B"), , drop = FALSE]
  restricted <- make_tasks(filtered, "one_vs_others", "A")
  restricted_task <- materialize_task_dataset(restricted$dataset, "one_vs_others", "A", NULL, "grouping", "others", "A")
  expect_identical(restricted_task$grouping[restricted_task$condition == "others"], c("B", "B"))
  expect_identical(make_tasks(filtered[filtered$grouping == "A", ], "condition_vs_condition")$task_table$temp_split, "A")
})
