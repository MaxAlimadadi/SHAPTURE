make_project_layout_data <- function() {
  data.frame(
    condition = rep(c("control", "case"), 6L),
    sample_id = rep(paste0("Patient", seq_len(6L)), each = 2L),
    grouping = rep(c("A", "B"), each = 6L), GeneA = seq_len(12L)
  )
}

test_that("project names are single folder names or NULL", {
  expect_null(validate_shapture_project_name(NULL))
  for (name in c("analysis_01", "run with spaces", "2026-11-06")) {
    expect_identical(validate_shapture_project_name(name), name)
  }
  for (name in list("", " ", ".", "..", "a/b", "a\\b", "/absolute", "a\nb",
                    NA_character_, character(), c("A", "B"), 1L, list("A"))) {
    expect_error(validate_shapture_project_name(name), "project_name must be NULL or one")
  }
})

test_that("a project creates one parent while NULL uses the base directly", {
  root <- tempfile("project-layout-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  layout <- prepare_shapture_output_layout(root, "analysis_01", "condition_vs_condition", "control", "case")
  expect_identical(layout$main_dir, normalizePath(root, winslash = "/"))
  expect_identical(layout$project_dir, normalizePath(file.path(root, "analysis_01"), winslash = "/"))
  expect_identical(layout$comparison_prefix, "control_case")
  direct <- prepare_shapture_output_layout(root, NULL, "one_vs_others", "control", "case")
  expect_identical(direct$project_dir, direct$main_dir)
  expect_identical(direct$comparison_prefix, "one_vs_others")
  expect_identical(prepare_shapture_output_layout(root, "analysis_01", "condition_vs_condition", "control", "case"), layout)
})

test_that("output layout handles relative directories without changing working directory", {
  root <- tempfile("project-relative-")
  dir.create(root)
  previous <- getwd()
  on.exit({setwd(previous); unlink(root, recursive = TRUE)}, add = TRUE)
  setwd(root)
  layout <- prepare_shapture_output_layout("relative/output", "run 01", "condition_vs_condition", "control", "case")
  expect_identical(getwd(), normalizePath(root, winslash = "/"))
  expect_identical(layout$project_dir, normalizePath(file.path(root, "relative/output/run 01"), winslash = "/"))
})

test_that("output layout reports invalid paths and file collisions", {
  root <- tempfile("project-files-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  for (value in list(NULL, "", NA_character_, c("A", "B"), 1L)) {
    expect_error(prepare_shapture_output_layout(value, NULL, "condition_vs_condition", "control", "case"), "main_dir must")
  }
  existing_file <- file.path(root, "existing_file")
  writeLines("keep", existing_file)
  expect_error(prepare_shapture_output_layout(existing_file, NULL, "condition_vs_condition", "control", "case"), "Cannot create the base")
  expect_error(prepare_shapture_output_layout(root, "existing_file", "condition_vs_condition", "control", "case"), "Cannot create the project")
  expect_identical(readLines(existing_file), "keep")
})

test_that("configuration resolution preserves NULL or the chosen project name", {
  data <- make_project_layout_data()
  for (mode in c("condition_vs_condition", "one_vs_others")) {
    for (project in list(NULL, "analysis_01")) {
      scope <- list2env(list(shapture_config = list(
        comparison_mode = mode, split_level_column = "grouping",
        what_0 = "control", what_1 = "case", project_name = project
      )))
      resolved <- resolve_shapture_script_config(scope, data, "GeneA", tempdir())
      expect_true("project_name" %in% names(resolved))
      expect_identical(resolved$project_name, project)
    }
  }
})

test_that("pipeline routes comparisons, audits and saved configurations into the project", {
  root <- tempfile("project-routing-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  original_directory <- getwd()
  sentinel <- file.path(root, "keep.txt")
  writeLines("unrelated file", sentinel)
  calls <- list()
  testthat::local_mocked_bindings(
    load_shapture_runtime = function() invisible(NULL),
    run_shapture_task = function(...) {
      arguments <- list(...)
      calls[[length(calls) + 1L]] <<- arguments
      output <- file.path(arguments$main_dir, paste0(arguments$comparison_prefix, "_", arguments$temp_split))
      dir.create(file.path(output, "docs"), recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(data.frame(feature = "GeneA", final_score = 0.5), file.path(output, "docs", "combined_results.csv"), row.names = FALSE)
      saveRDS(list(), file.path(output, "docs", "visualization_inputs.rds"))
      write_shapture_config_csv(arguments$resolved_config_for_docs, file.path(output, "docs", "shapture_config_used.csv"))
      list(temp_split = arguments$temp_split, analyzed = TRUE,
        each_condition_sample_no = 2L, output_dir = output)
    },
    .package = "SHAPTURE"
  )
  for (mode in c("condition_vs_condition", "one_vs_others")) {
    for (project in list("analysis_01", NULL)) {
      base <- file.path(root, mode, if (is.null(project)) "direct" else "named")
      result <- run_shapture_pipeline(
        make_project_layout_data(), "GeneA", base,
        split_level_column = "grouping", comparison_mode = mode,
        what_0 = "control", what_1 = "case", project_name = project,
        analysis_unit_type = "single_cell", n_cores = 1L
      )
      prefix <- if (mode == "condition_vs_condition") "control_case" else "one_vs_others"
      expected_parent <- normalizePath(if (is.null(project)) base else file.path(base, project), winslash = "/")
      expect_identical(result$project_dir, expected_parent)
      expect_identical(unname(result$output_dirs), file.path(expected_parent, paste0(prefix, "_", c("A", "B"))))
      expect_identical(find_shapture_result_folders(result$project_dir), sort(unname(result$output_dirs)))
      expect_true(file.exists(file.path(expected_parent, "cell.used.per.condition.csv")))
      expect_length(list.files(expected_parent, "before[.]loop[.]csv$"), 1L)
      if (!is.null(project)) {
        expect_false(file.exists(file.path(base, "cell.used.per.condition.csv")))
      }
      saved <- utils::read.csv(file.path(result$output_dirs[[1]], "docs", "shapture_config_used.csv"), stringsAsFactors = FALSE)
      read_value <- function(name) saved$parameter_used[saved$parameter_name == name]
      expect_identical(read_value("project_name"), if (is.null(project)) "NULL" else project)
      expect_identical(read_value("project_dir"), expected_parent)
      expect_identical(read_value("main_dir"), normalizePath(base, winslash = "/"))
      expect_identical(read_value("comparison_prefix"), prefix)
      expect_identical(getwd(), original_directory)
    }
  }
  expect_identical(readLines(sentinel), "unrelated file")
  expect_length(calls, 8L)
})

test_that("metacell construction receives identical scientific settings across project names", {
  root <- tempfile("project-metacells-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  metacell_calls <- list()
  testthat::local_mocked_bindings(
    load_shapture_runtime = function() invisible(NULL),
    build_grouped_metacells = function(...) {
      settings <- list(...)
      metacell_calls[[length(metacell_calls) + 1L]] <<- settings
      list(data = settings$dataset, summary = NULL)
    },
    run_shapture_task = function(...) {
      settings <- list(...)
      list(temp_split = settings$temp_split, analyzed = FALSE, each_condition_sample_no = 0L)
    },
    .package = "SHAPTURE"
  )
  for (project in list(NULL, "analysis_01")) {
    run_shapture_pipeline(make_project_layout_data(), "GeneA", root,
      project_name = project, analysis_unit_type = "metacell", what_0 = "control", what_1 = "case")
  }
  expect_length(metacell_calls, 2L)
  expect_identical(metacell_calls[[1]], metacell_calls[[2]])
})
