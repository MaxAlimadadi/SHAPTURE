test_that("analysis_seed defaults to a fixed reproducible value", {
  config <- SHAPTURE:::default_shapture_script_config(c("gene_a", "gene_b"))

  expect_identical(config$analysis_seed, 123L)
  expect_identical(SHAPTURE:::validate_analysis_seed(123), 123L)
  expect_error(
    SHAPTURE:::validate_analysis_seed(0),
    "between 1 and 999999999"
  )
  expect_error(
    SHAPTURE:::validate_analysis_seed(1.5),
    "whole number"
  )
})


test_that("task seed plans are stable and independent of iteration count", {
  make_plan <- function(iteration_no) {
    SHAPTURE:::build_shapture_task_seed_plan(
      analysis_seed = 123L,
      comparison_mode = "condition_vs_condition",
      split_by_col = "cell_type",
      temp_split = "subset_A",
      task_what_0 = "control",
      task_what_1 = "case",
      iteration_no = iteration_no
    )
  }

  plan_30_a <- make_plan(30L)
  plan_30_b <- make_plan(30L)
  plan_60 <- make_plan(60L)

  expect_identical(plan_30_a, plan_30_b)
  expect_identical(
    plan_30_a$iteration_seeds,
    plan_60$iteration_seeds[seq_len(30), , drop = FALSE]
  )
  expect_equal(nrow(plan_60$iteration_seeds), 60L)
  for (seed_col in c(
    "cell_sampling_seed",
    "feature_selection_seed",
    "model_seed"
  )) {
    expect_false(anyDuplicated(plan_30_a$iteration_seeds[[seed_col]]) > 0L)
  }
})


test_that("task seeds depend on stable task identity rather than task order", {
  plan_a_first <- SHAPTURE:::build_shapture_task_seed_plan(
    123L, "condition_vs_condition", "cell_type", "subset_A",
    "control", "case", 5L
  )
  invisible(SHAPTURE:::build_shapture_task_seed_plan(
    123L, "condition_vs_condition", "cell_type", "subset_B",
    "control", "case", 5L
  ))
  plan_a_after <- SHAPTURE:::build_shapture_task_seed_plan(
    123L, "condition_vs_condition", "cell_type", "subset_A",
    "control", "case", 5L
  )
  plan_b <- SHAPTURE:::build_shapture_task_seed_plan(
    123L, "condition_vs_condition", "cell_type", "subset_B",
    "control", "case", 5L
  )

  expect_identical(plan_a_first, plan_a_after)
  expect_false(identical(plan_a_first$task_seed, plan_b$task_seed))
  expect_false(identical(
    plan_a_first$iteration_seeds$model_seed,
    plan_b$iteration_seeds$model_seed
  ))
})


test_that("changing analysis_seed creates a distinct reproducible plan", {
  plan_123 <- SHAPTURE:::build_shapture_task_seed_plan(
    123L, "one_vs_others", "cluster", "cluster_1",
    "others", "cluster_1", 10L
  )
  plan_124 <- SHAPTURE:::build_shapture_task_seed_plan(
    124L, "one_vs_others", "cluster", "cluster_1",
    "others", "cluster_1", 10L
  )

  expect_false(identical(plan_123$task_seed, plan_124$task_seed))
  expect_false(identical(plan_123$iteration_seeds, plan_124$iteration_seeds))
  expect_false(identical(
    plan_123$cross_patient_bootstrap_seed,
    plan_124$cross_patient_bootstrap_seed
  ))
})
