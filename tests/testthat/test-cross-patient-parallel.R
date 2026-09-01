make_parallel_controller_fixture <- function(n_plans = 12L) {
  worker <- new.env(parent = baseenv())
  worker$n_features <- 2L
  kernel <- function(plan) {
    plan <- as.numeric(plan)
    list(
      importance = c(plan, plan + 1),
      expression = c(0.25, 0.75),
      value_shap = c(-0.5, 0.5),
      valid = TRUE
    )
  }
  environment(kernel) <- worker

  list(
    plans = as.list(seq_len(n_plans)),
    kernel = kernel,
    memory_estimate = data.frame(
      features = 2L,
      max_shap_rows = 10L,
      max_model_rows = 10L,
      max_single_cell_rows = 10L,
      exported_input_gib = 0.001,
      working_set_gib = 0.001,
      static_worker_estimate_gib = 0.05,
      stringsAsFactors = FALSE
    ),
    resources = list(
      cores_override = 4L,
      total_ram_gib_override = 128,
      available_ram_gib_override = 120,
      platform_override = "test"
    )
  )
}


test_that("Cross-patient workers respect CPU, RAM, and user limits", {
  abundant <- data.frame(
    total_ram_gib = 64,
    available_ram_gib = 56,
    available_cores = 8L
  )
  constrained <- data.frame(
    total_ram_gib = 16,
    available_ram_gib = 4,
    available_cores = 8L
  )

  cpu_limited <- select_cross_patient_workers(
    requested_workers = 12L,
    resources = abundant,
    master_rss_gib = 2,
    worker_estimate_gib = 1
  )
  memory_limited <- select_cross_patient_workers(
    requested_workers = 8L,
    resources = constrained,
    master_rss_gib = 2,
    worker_estimate_gib = 1.5
  )
  unavailable <- select_cross_patient_workers(
    requested_workers = 8L,
    resources = transform(abundant, total_ram_gib = NA_real_),
    master_rss_gib = 2,
    worker_estimate_gib = 1
  )

  expect_equal(cpu_limited$selected_workers, 8L)
  expect_equal(cpu_limited$backend, "psock_recycled")
  expect_equal(memory_limited$selected_workers, 1L)
  expect_equal(memory_limited$backend, "serial")
  expect_equal(unavailable$selected_workers, 1L)
  expect_equal(unavailable$selection_reason, "resource information unavailable")
})


test_that("Recycled PSOCK execution is exactly equal to serial execution", {
  skip_on_cran()
  fixture <- make_parallel_controller_fixture()
  serial_values <- lapply(fixture$plans, fixture$kernel)

  result <- run_cross_patient_bootstrap_controller(
    plans = fixture$plans,
    kernel = fixture$kernel,
    requested_workers = 2L,
    memory_estimate = fixture$memory_estimate,
    worker_log_file = tempfile("shapture-worker-", fileext = ".log"),
    recycle_every = 3L,
    resource_overrides = fixture$resources
  )

  expect_identical(result$values, serial_values)
  expect_equal(result$runtime$backend, "psock_recycled")
  expect_equal(result$runtime$selected_workers, 2L)
  expect_true(result$runtime$health_check_passed)
  expect_true(result$runtime$health_result_retained)
  expect_false(result$runtime$fallback_used)
  expect_equal(result$runtime$worker_recycling_interval, 3L)
})


test_that("Full Cross-patient statistics are identical in serial and parallel", {
  skip_on_cran()
  test_data <- make_cross_patient_test_data(3L)
  arguments <- list(
    shap_matrix = test_data$shap,
    shap_feature_values = test_data$expression,
    shap_metadata = test_data$metadata,
    model_expression_data = test_data$model_data,
    single_cell_expression_data = test_data$single_cell_data,
    full_stats = test_data$full_stats,
    average_importance = test_data$average_importance,
    feature_names = colnames(test_data$shap),
    comparison_mode = "condition_vs_condition",
    n_bootstraps = 50L,
    seed = 887L
  )
  serial_result <- do.call(
    compute_cross_patient_support,
    c(arguments, list(n_cores = 1L))
  )

  testthat::local_mocked_bindings(
    detect_shapture_bootstrap_resources = function(...) {
      data.frame(
        platform = "test",
        available_cores = 4L,
        total_ram_gib = 128,
        available_ram_gib = 120,
        cpu_source = "test override",
        ram_source = "test override",
        stringsAsFactors = FALSE
      )
    },
    .package = "SHAPTURE"
  )
  parallel_result <- do.call(
    compute_cross_patient_support,
    c(arguments, list(n_cores = 2L))
  )

  serial_runtime <- serial_result$runtime
  parallel_runtime <- parallel_result$runtime
  serial_result$runtime <- NULL
  parallel_result$runtime <- NULL
  expect_identical(parallel_result, serial_result)
  expect_equal(serial_runtime$backend, "serial")
  expect_equal(parallel_runtime$backend, "psock_recycled")
  expect_false(parallel_runtime$fallback_used)
})


test_that("Cluster setup failure replays every bootstrap serially", {
  fixture <- make_parallel_controller_fixture(8L)
  serial_values <- lapply(fixture$plans, fixture$kernel)
  failing_cluster_factory <- function(...) {
    stop("test cluster setup failure")
  }

  result <- run_cross_patient_bootstrap_controller(
    plans = fixture$plans,
    kernel = fixture$kernel,
    requested_workers = 2L,
    memory_estimate = fixture$memory_estimate,
    resource_overrides = fixture$resources,
    cluster_factory = failing_cluster_factory
  )

  expect_identical(result$values, serial_values)
  expect_equal(result$runtime$backend, "serial")
  expect_true(result$runtime$fallback_used)
  expect_match(result$runtime$fallback_reason, "test cluster setup failure")
})


test_that("Mid-run worker failure replays every bootstrap serially", {
  skip_on_cran()
  fixture <- make_parallel_controller_fixture(8L)
  serial_values <- lapply(fixture$plans, fixture$kernel)
  failing_dispatch <- function(item, gc_interval = 1L) {
    stop("test mid-run worker failure")
  }

  result <- run_cross_patient_bootstrap_controller(
    plans = fixture$plans,
    kernel = fixture$kernel,
    requested_workers = 2L,
    memory_estimate = fixture$memory_estimate,
    worker_log_file = tempfile("shapture-worker-", fileext = ".log"),
    recycle_every = 2L,
    resource_overrides = fixture$resources,
    batch_dispatch = failing_dispatch
  )

  expect_identical(result$values, serial_values)
  expect_true(result$runtime$fallback_used)
  expect_equal(result$runtime$backend, "serial")
  expect_equal(result$runtime$selected_workers, 1L)
  expect_match(result$runtime$fallback_reason, "test mid-run worker failure")
  expect_false(result$runtime$health_result_retained)
})


test_that("Invalid parallel output triggers a complete serial replay", {
  skip_on_cran()
  fixture <- make_parallel_controller_fixture(8L)
  serial_values <- lapply(fixture$plans, fixture$kernel)
  invalid_dispatch <- function(item, gc_interval = 1L) {
    lapply(item$plans, function(plan) {
      list(
        importance = 1,
        expression = 2,
        value_shap = 3,
        valid = TRUE
      )
    })
  }

  result <- run_cross_patient_bootstrap_controller(
    plans = fixture$plans,
    kernel = fixture$kernel,
    requested_workers = 2L,
    memory_estimate = fixture$memory_estimate,
    worker_log_file = tempfile("shapture-worker-", fileext = ".log"),
    recycle_every = 2L,
    resource_overrides = fixture$resources,
    batch_dispatch = invalid_dispatch
  )

  expect_identical(result$values, serial_values)
  expect_true(result$runtime$fallback_used)
  expect_equal(result$runtime$backend, "serial")
  expect_equal(result$runtime$selected_workers, 1L)
  expect_match(result$runtime$fallback_reason, "invalid|outside", ignore.case = TRUE)
})


test_that("Unavailable Cross-patient work has an explicit runtime state", {
  runtime <- unavailable_cross_patient_runtime(6L, "insufficient patients")

  expect_equal(runtime$requested_workers, 6L)
  expect_equal(runtime$selected_workers, 0L)
  expect_equal(runtime$backend, "not_applicable")
  expect_equal(runtime$selection_reason, "insufficient patients")
  expect_false(runtime$fallback_used)
})
