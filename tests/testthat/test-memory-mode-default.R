make_memory_default_data <- function() {
  data.frame(
    condition = rep(c("control", "case"), each = 6L),
    sample_id = rep(paste0("Sample", seq_len(4L)), each = 3L),
    grouping = rep(c("A", "B"), 6L),
    GeneA = seq_len(12L)
  )
}

test_that("the package and script configuration default to safe memory mode", {
  expect_identical(formals(run_shapture_pipeline)$memory_mode, "safe")
  expect_identical(default_shapture_script_config("GeneA")$memory_mode, "safe")
})

test_that("both comparison modes resolve to safe unless full is requested", {
  data <- make_memory_default_data()
  for (comparison in c("condition_vs_condition", "one_vs_others")) {
    config <- list(
      comparison_mode = comparison, split_level_column = "grouping",
      what_0 = "control", what_1 = "case"
    )
    scope <- list2env(list(shapture_config = config))
    resolved <- resolve_shapture_script_config(scope, data, "GeneA", tempdir())
    expect_identical(resolved$memory_mode, "safe")

    scope$shapture_config$memory_mode <- "full"
    resolved <- resolve_shapture_script_config(scope, data, "GeneA", tempdir())
    expect_identical(resolved$memory_mode, "full")

    scope$memory_mode <- "safe"
    resolved <- resolve_shapture_script_config(scope, data, "GeneA", tempdir())
    expect_identical(resolved$memory_mode, "safe")
  }
})

test_that("invalid memory modes still fail clearly", {
  scope <- list2env(list(shapture_config = list(
    what_0 = "control", what_1 = "case", memory_mode = "invalid"
  )))
  expect_error(
    resolve_shapture_script_config(scope, make_memory_default_data(), "GeneA", tempdir()),
    "memory_mode must be one of"
  )
})

test_that("saved default configuration records safe memory mode", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  write_shapture_config_csv(default_shapture_script_config("GeneA"), path)
  saved <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_identical(saved$parameter_used[saved$parameter_name == "memory_mode"], "safe")
})
