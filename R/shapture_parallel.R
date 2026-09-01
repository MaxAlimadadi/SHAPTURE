utils::globalVariables(c(
  "feature_chunks",
  "full_detection_sign",
  "full_log2fc",
  "model_condition_all",
  "model_expression_matrix",
  "n_features",
  "shared_expression_data",
  "single_cell_condition_all",
  "single_cell_expression_matrix"
))


#' Detect CPU and memory resources for Cross-patient bootstrapping
#'
#' @param cores_override Optional test-only CPU override.
#' @param total_ram_gib_override Optional test-only installed-RAM override.
#' @param available_ram_gib_override Optional test-only available-RAM override.
#' @param platform_override Optional test-only platform override.
#'
#' @return One-row resource summary.
#' @noRd
detect_shapture_bootstrap_resources <- function(
    cores_override = NULL,
    total_ram_gib_override = NULL,
    available_ram_gib_override = NULL,
    platform_override = NULL
) {
  gib <- 1024^3
  cpu_source <- "parallelly::availableCores"
  available_cores <- cores_override
  if (is.null(available_cores)) {
    available_cores <- tryCatch(
      suppressWarnings(as.integer(parallelly::availableCores(omit = 0L)[1L])),
      error = function(e) NA_integer_
    )
  } else {
    cpu_source <- "test override"
  }
  if (!is.finite(available_cores) || available_cores < 1L) {
    available_cores <- 1L
    cpu_source <- paste0(cpu_source, "; conservative fallback")
  }

  total_ram_gib <- total_ram_gib_override
  available_ram_gib <- available_ram_gib_override
  ram_source <- "ps::ps_system_memory"
  if (is.null(total_ram_gib) || is.null(available_ram_gib)) {
    memory <- tryCatch(ps::ps_system_memory(), error = function(e) NULL)
    if (!is.null(memory)) {
      if (is.null(total_ram_gib)) {
        total_ram_gib <- as.numeric(memory$total) / gib
      }
      if (is.null(available_ram_gib)) {
        available_ram_gib <- as.numeric(memory$avail) / gib
      }
    }
  } else {
    ram_source <- "test override"
  }
  if (!is.finite(total_ram_gib) || total_ram_gib <= 0 ||
      !is.finite(available_ram_gib) || available_ram_gib < 0) {
    total_ram_gib <- NA_real_
    available_ram_gib <- NA_real_
    ram_source <- paste0(ram_source, "; unavailable")
  }

  platform <- platform_override
  if (is.null(platform)) {
    sysname <- tolower(as.character(Sys.info()[["sysname"]]))
    platform <- switch(
      sysname,
      darwin = "macos",
      windows = "windows",
      linux = "linux",
      "unknown"
    )
  }

  data.frame(
    platform = platform,
    available_cores = as.integer(available_cores),
    total_ram_gib = total_ram_gib,
    available_ram_gib = available_ram_gib,
    cpu_source = cpu_source,
    ram_source = ram_source,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


#' Read current R-process memory use
#'
#' @return Resident memory in GiB, or `NA_real_` when unavailable.
#' @noRd
current_shapture_process_rss_gib <- function() {
  value <- tryCatch(
    as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]]) / 1024^3,
    error = function(e) NA_real_
  )
  if (is.finite(value) && value >= 0) value else NA_real_
}


#' Estimate memory required by one Cross-patient bootstrap worker
#'
#' @param kernel Bootstrap kernel with its compact worker environment.
#' @param plans Precomputed deterministic patient-block plans.
#' @param chunk_size Feature chunk size used by the kernel.
#'
#' @return One-row memory summary.
#' @noRd
estimate_cross_patient_worker_memory <- function(kernel, plans, chunk_size = 64L) {
  gib <- 1024^3
  worker <- environment(kernel)
  n_features <- get("n_features", envir = worker, inherits = FALSE)
  input_names <- c(
    "shap_matrix", "shap_feature_values", "model_expression_matrix",
    "model_condition_all", "single_cell_expression_matrix",
    "single_cell_condition_all", "full_log2fc", "full_detection_sign"
  )
  input_names <- input_names[vapply(
    input_names,
    exists,
    logical(1),
    envir = worker,
    inherits = FALSE
  )]
  exported_input_bytes <- sum(vapply(input_names, function(name) {
    as.numeric(utils::object.size(get(name, envir = worker, inherits = FALSE)))
  }, numeric(1)))

  valid_plan_index <- which(!vapply(plans, is.null, logical(1)))[1L]
  source_names <- if (is.na(valid_plan_index)) character(0) else names(plans[[valid_plan_index]])
  max_rows <- stats::setNames(numeric(length(source_names)), source_names)
  for (source_name in source_names) {
    max_rows[source_name] <- max(vapply(plans, function(plan) {
      if (is.null(plan)) 0L else length(plan[[source_name]])
    }, integer(1)))
  }
  get_max_rows <- function(name) {
    if (name %in% names(max_rows)) unname(max_rows[[name]]) else 0
  }

  model_bytes <- get_max_rows("model") * n_features * 8
  single_cell_bytes <- get_max_rows("single_cell") * n_features * 8
  shap_chunk_bytes <- get_max_rows("shap") * min(n_features, chunk_size) * 8
  working_bytes <- model_bytes + single_cell_bytes + 4 * shap_chunk_bytes

  # This multiplier covers ordinary R copies and temporary allocations. A real
  # health-worker RSS measurement replaces it whenever parallel work is viable.
  static_worker_bytes <- exported_input_bytes + 3 * working_bytes + 256 * 1024^2
  data.frame(
    features = n_features,
    max_shap_rows = get_max_rows("shap"),
    max_model_rows = get_max_rows("model"),
    max_single_cell_rows = get_max_rows("single_cell"),
    exported_input_gib = exported_input_bytes / gib,
    working_set_gib = working_bytes / gib,
    static_worker_estimate_gib = static_worker_bytes / gib,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


#' Select a memory-safe Cross-patient bootstrap worker count
#'
#' @param requested_workers User-requested maximum workers.
#' @param resources Resource summary.
#' @param master_rss_gib Current main-process memory use.
#' @param worker_estimate_gib Estimated memory per worker.
#' @param reserve_fraction Fraction of installed RAM kept in reserve.
#' @param minimum_reserve_gib Minimum RAM reserve.
#'
#' @return One-row worker-selection summary.
#' @noRd
select_cross_patient_workers <- function(
    requested_workers,
    resources,
    master_rss_gib,
    worker_estimate_gib,
    reserve_fraction = 0.25,
    minimum_reserve_gib = 4
) {
  if (any(!is.finite(c(
    resources$total_ram_gib,
    resources$available_ram_gib,
    master_rss_gib,
    worker_estimate_gib
  ))) || worker_estimate_gib <= 0) {
    return(data.frame(
      requested_workers = requested_workers,
      reserve_gib = NA_real_,
      memory_safe_workers = 1L,
      selected_workers = 1L,
      backend = "serial",
      selection_reason = "resource information unavailable",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  reserve_gib <- max(minimum_reserve_gib, reserve_fraction * resources$total_ram_gib)
  usable_gib <- max(
    0,
    min(resources$available_ram_gib, resources$total_ram_gib - reserve_gib) - master_rss_gib
  )
  memory_safe_workers <- floor(usable_gib / worker_estimate_gib)
  selected_workers <- min(
    as.integer(requested_workers),
    as.integer(resources$available_cores),
    as.integer(memory_safe_workers)
  )
  if (selected_workers < 2L) {
    selected_workers <- 1L
    backend <- "serial"
    selection_reason <- "fewer than two workers are safe"
  } else {
    backend <- "psock_recycled"
    selection_reason <- "parallel resources available"
  }

  data.frame(
    requested_workers = requested_workers,
    reserve_gib = reserve_gib,
    memory_safe_workers = as.integer(memory_safe_workers),
    selected_workers = as.integer(selected_workers),
    backend = backend,
    selection_reason = selection_reason,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


#' Validate one Cross-patient bootstrap result
#'
#' @param value One bootstrap result.
#' @param n_features Expected feature count.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
validate_cross_patient_bootstrap_value <- function(value, n_features) {
  expected_names <- c("importance", "expression", "value_shap", "valid")
  if (!is.list(value) || !identical(sort(names(value)), sort(expected_names))) {
    stop("Cross-patient bootstrap output has an invalid structure.")
  }
  if (length(value$importance) != n_features ||
      length(value$expression) != n_features ||
      length(value$value_shap) != n_features ||
      length(value$valid) != 1L ||
      !is.logical(value$valid)) {
    stop("Cross-patient bootstrap output dimensions are invalid.")
  }
  if (any(is.finite(value$importance) & value$importance < 0)) {
    stop("Cross-patient SHAP-importance output contains negative values.")
  }
  if (any(is.finite(value$expression) & (value$expression < 0 | value$expression > 1))) {
    stop("Cross-patient Expression-signal output is outside 0-1.")
  }
  if (any(is.finite(value$value_shap) & abs(value$value_shap) > 1 + 1e-12)) {
    stop("Cross-patient Value-SHAP output is outside -1 to 1.")
  }
  invisible(TRUE)
}


#' Validate all Cross-patient bootstrap results
#'
#' @param values Ordered bootstrap results.
#' @param n_plans Expected plan count.
#' @param n_features Expected feature count.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
validate_cross_patient_bootstrap_values <- function(values, n_plans, n_features) {
  if (length(values) != n_plans) {
    stop("Cross-patient bootstrap result count does not match the plan count.")
  }
  invisible(lapply(values, validate_cross_patient_bootstrap_value, n_features = n_features))
  invisible(TRUE)
}


cross_patient_health_dispatch <- function(plan) {
  value <- get("cross_patient_bootstrap_kernel", envir = .GlobalEnv)(plan)
  rss_gib <- tryCatch(
    as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]]) / 1024^3,
    error = function(e) NA_real_
  )
  list(value = value, worker_rss_gib = rss_gib)
}
environment(cross_patient_health_dispatch) <- baseenv()


cross_patient_batch_dispatch <- function(item, gc_interval = 1L) {
  kernel <- get("cross_patient_bootstrap_kernel", envir = .GlobalEnv)
  values <- vector("list", length(item$plans))
  for (i in seq_along(item$plans)) {
    values[[i]] <- kernel(item$plans[[i]])
    if (gc_interval > 0L && i %% gc_interval == 0L) {
      gc(verbose = FALSE)
    }
  }
  values
}
environment(cross_patient_batch_dispatch) <- baseenv()


#' Stop a PSOCK cluster without masking an earlier error
#'
#' @param cluster Cluster object or `NULL`.
#'
#' @return Invisibly returns `TRUE` when no cluster is supplied or the cluster
#'   stops cleanly, and `FALSE` otherwise.
#' @noRd
safe_stop_shapture_cluster <- function(cluster) {
  if (is.null(cluster)) {
    return(invisible(TRUE))
  }
  stopped <- !inherits(
    try(parallel::stopCluster(cluster), silent = TRUE),
    "try-error"
  )
  invisible(stopped)
}


#' Run one PSOCK worker health check
#'
#' @param plan First deterministic bootstrap plan.
#' @param kernel Bootstrap kernel.
#' @param worker_log_file Worker output file.
#' @param cluster_factory Injectable cluster constructor for tests.
#' @param health_dispatch Injectable health dispatcher for tests.
#'
#' @return Health-check result and measured worker RSS.
#' @noRd
run_cross_patient_health_check <- function(
    plan,
    kernel,
    worker_log_file,
    cluster_factory = parallel::makePSOCKcluster,
    health_dispatch = cross_patient_health_dispatch
) {
  cluster <- cluster_factory(1L, outfile = worker_log_file)
  on.exit(safe_stop_shapture_cluster(cluster), add = TRUE)
  export_env <- new.env(parent = emptyenv())
  export_env$cross_patient_bootstrap_kernel <- kernel
  parallel::clusterExport(
    cluster,
    "cross_patient_bootstrap_kernel",
    envir = export_env
  )
  parallel::clusterCall(cluster, health_dispatch, plan = plan)[[1L]]
}


#' Run one recycled PSOCK wave
#'
#' @param items Bootstrap batches assigned to the current workers.
#' @param kernel Bootstrap kernel.
#' @param worker_log_file Worker output file.
#' @param cluster_factory Injectable cluster constructor for tests.
#' @param batch_dispatch Injectable batch dispatcher for tests.
#'
#' @return List of batch results.
#' @noRd
run_cross_patient_psock_wave <- function(
    items,
    kernel,
    worker_log_file,
    cluster_factory = parallel::makePSOCKcluster,
    batch_dispatch = cross_patient_batch_dispatch
) {
  cluster <- cluster_factory(length(items), outfile = worker_log_file)
  on.exit(safe_stop_shapture_cluster(cluster), add = TRUE)
  export_env <- new.env(parent = emptyenv())
  export_env$cross_patient_bootstrap_kernel <- kernel
  parallel::clusterExport(
    cluster,
    "cross_patient_bootstrap_kernel",
    envir = export_env
  )
  parallel::parLapply(
    cluster,
    items,
    batch_dispatch,
    gc_interval = 1L
  )
}


#' Run deterministic bootstrap plans through recycled PSOCK workers
#'
#' @param plans Ordered bootstrap plans.
#' @param kernel Bootstrap kernel.
#' @param workers Number of workers.
#' @param retained_health_value Result from the first health-check plan.
#' @param worker_log_file Worker output file.
#' @param recycle_every Maximum plans handled by one worker process.
#' @param cluster_factory Injectable cluster constructor for tests.
#' @param batch_dispatch Injectable batch dispatcher for tests.
#'
#' @return Ordered list of bootstrap results.
#' @noRd
run_recycled_cross_patient_psock <- function(
    plans,
    kernel,
    workers,
    retained_health_value,
    worker_log_file,
    recycle_every = 5L,
    cluster_factory = parallel::makePSOCKcluster,
    batch_dispatch = cross_patient_batch_dispatch
) {
  remaining_indices <- seq_along(plans)[-1L]
  if (!length(remaining_indices)) {
    return(list(retained_health_value))
  }

  batch_id <- ceiling(seq_along(remaining_indices) / recycle_every)
  index_batches <- split(remaining_indices, batch_id)
  wave_id <- ceiling(seq_along(index_batches) / workers)
  waves <- split(seq_along(index_batches), wave_id)
  values <- vector("list", length(plans))
  values[[1L]] <- retained_health_value

  for (wave in waves) {
    items <- lapply(wave, function(batch_index) {
      list(id = batch_index, plans = plans[index_batches[[batch_index]]])
    })
    wave_values <- run_cross_patient_psock_wave(
      items = items,
      kernel = kernel,
      worker_log_file = worker_log_file,
      cluster_factory = cluster_factory,
      batch_dispatch = batch_dispatch
    )
    for (j in seq_along(wave)) {
      values[index_batches[[wave[[j]]]]] <- wave_values[[j]]
    }
  }
  values
}


#' Build an unavailable Cross-patient runtime summary
#'
#' @param requested_workers User-requested worker maximum.
#' @param reason Why execution was not applicable.
#'
#' @return One-row runtime summary.
#' @noRd
unavailable_cross_patient_runtime <- function(requested_workers, reason) {
  data.frame(
    platform = NA_character_,
    requested_workers = as.integer(requested_workers),
    available_cores = NA_integer_,
    total_ram_gib = NA_real_,
    available_ram_gib = NA_real_,
    master_rss_gib = NA_real_,
    features = NA_integer_,
    max_shap_rows = NA_integer_,
    max_model_rows = NA_integer_,
    max_single_cell_rows = NA_integer_,
    exported_input_gib = NA_real_,
    static_worker_estimate_gib = NA_real_,
    health_worker_rss_gib = NA_real_,
    calibrated_worker_estimate_gib = NA_real_,
    reserved_ram_gib = NA_real_,
    memory_safe_workers = NA_integer_,
    selected_workers = 0L,
    backend = "not_applicable",
    selection_reason = reason,
    health_check_passed = NA,
    health_result_retained = FALSE,
    worker_recycling_interval = NA_integer_,
    fallback_used = FALSE,
    fallback_reason = NA_character_,
    patient_plan_seconds = 0,
    compute_seconds = 0,
    total_seconds = 0,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


#' Execute Cross-patient bootstrap plans with safe parallel fallback
#'
#' @param plans Ordered deterministic patient-block plans.
#' @param kernel Scientific bootstrap kernel.
#' @param requested_workers User-requested maximum workers.
#' @param memory_estimate Static worker-memory estimate.
#' @param patient_plan_seconds Time used to create deterministic plans.
#' @param worker_log_file Worker output file, or `NULL`.
#' @param recycle_every Maximum bootstraps per worker process.
#' @param resource_overrides Optional test-only resource overrides.
#' @param cluster_factory Injectable cluster constructor for tests.
#' @param health_dispatch Injectable health dispatcher for tests.
#' @param batch_dispatch Injectable batch dispatcher for tests.
#' @param output_validator Injectable output validator for tests.
#'
#' @return List containing ordered values and a runtime summary.
#' @noRd
run_cross_patient_bootstrap_controller <- function(
    plans,
    kernel,
    requested_workers,
    memory_estimate,
    patient_plan_seconds = 0,
    worker_log_file = NULL,
    recycle_every = 5L,
    resource_overrides = list(),
    cluster_factory = parallel::makePSOCKcluster,
    health_dispatch = cross_patient_health_dispatch,
    batch_dispatch = cross_patient_batch_dispatch,
    output_validator = validate_cross_patient_bootstrap_values
) {
  requested_workers <- as.integer(requested_workers)
  if (!is.finite(requested_workers) || requested_workers < 1L) {
    stop("requested_workers must be a positive integer.")
  }
  if (!is.numeric(recycle_every) || length(recycle_every) != 1L ||
      !is.finite(recycle_every) || recycle_every < 1L ||
      recycle_every != floor(recycle_every)) {
    stop("recycle_every must be a positive integer.")
  }
  recycle_every <- as.integer(recycle_every)
  worker_log_file <- worker_log_file %||% ""
  n_features <- get("n_features", envir = environment(kernel), inherits = FALSE)
  started <- unname(proc.time()[["elapsed"]])
  resource_args <- c(
    list(),
    resource_overrides[c(
      "cores_override", "total_ram_gib_override",
      "available_ram_gib_override", "platform_override"
    )]
  )
  resources <- do.call(detect_shapture_bootstrap_resources, resource_args)
  master_rss_gib <- current_shapture_process_rss_gib()
  if (!is.finite(master_rss_gib)) {
    master_rss_gib <- as.numeric(utils::object.size(list(plans = plans, kernel = kernel))) / 1024^3
  }
  selection <- select_cross_patient_workers(
    requested_workers = requested_workers,
    resources = resources,
    master_rss_gib = master_rss_gib,
    worker_estimate_gib = memory_estimate$static_worker_estimate_gib
  )

  health_passed <- FALSE
  health_result <- NULL
  health_worker_rss_gib <- NA_real_
  calibrated_worker_estimate_gib <- memory_estimate$static_worker_estimate_gib
  fallback_used <- FALSE
  fallback_reason <- NA_character_

  if (selection$backend == "psock_recycled") {
    health <- tryCatch({
      answer <- run_cross_patient_health_check(
        plan = plans[[1L]],
        kernel = kernel,
        worker_log_file = worker_log_file,
        cluster_factory = cluster_factory,
        health_dispatch = health_dispatch
      )
      validate_cross_patient_bootstrap_value(answer$value, n_features)
      answer
    }, error = function(e) {
      fallback_reason <<- conditionMessage(e)
      NULL
    })

    if (is.null(health)) {
      selection$backend <- "serial"
      selection$selected_workers <- 1L
      selection$selection_reason <- "health check failed"
      fallback_used <- TRUE
    } else {
      health_passed <- TRUE
      health_result <- health$value
      health_worker_rss_gib <- health$worker_rss_gib
      measured_worker_gib <- if (is.finite(health_worker_rss_gib)) {
        1.25 * health_worker_rss_gib
      } else {
        NA_real_
      }
      calibrated_worker_estimate_gib <- max(
        memory_estimate$static_worker_estimate_gib,
        measured_worker_gib,
        na.rm = TRUE
      )
      resources <- do.call(detect_shapture_bootstrap_resources, resource_args)
      master_after_health_gib <- current_shapture_process_rss_gib()
      if (!is.finite(master_after_health_gib)) {
        master_after_health_gib <- master_rss_gib
      }
      selection <- select_cross_patient_workers(
        requested_workers = requested_workers,
        resources = resources,
        master_rss_gib = master_after_health_gib,
        worker_estimate_gib = calibrated_worker_estimate_gib
      )
    }
  }

  compute_started <- unname(proc.time()[["elapsed"]])
  if (selection$backend == "serial") {
    values <- lapply(plans, kernel)
  } else {
    values <- tryCatch({
      answer <- run_recycled_cross_patient_psock(
        plans = plans,
        kernel = kernel,
        workers = selection$selected_workers,
        retained_health_value = health_result,
        worker_log_file = worker_log_file,
        recycle_every = recycle_every,
        cluster_factory = cluster_factory,
        batch_dispatch = batch_dispatch
      )
      output_validator(answer, length(plans), n_features)
      answer
    }, error = function(e) {
      fallback_used <<- TRUE
      fallback_reason <<- conditionMessage(e)
      gc(verbose = FALSE)
      lapply(plans, kernel)
    })
  }
  output_validator(values, length(plans), n_features)
  if (fallback_used && selection$backend == "psock_recycled") {
    selection$backend <- "serial"
    selection$selected_workers <- 1L
    selection$selection_reason <- "parallel execution failed; complete serial replay used"
  }
  compute_seconds <- unname(proc.time()[["elapsed"]]) - compute_started

  runtime <- data.frame(
    platform = resources$platform,
    requested_workers = requested_workers,
    available_cores = resources$available_cores,
    total_ram_gib = resources$total_ram_gib,
    available_ram_gib = resources$available_ram_gib,
    master_rss_gib = master_rss_gib,
    features = memory_estimate$features,
    max_shap_rows = memory_estimate$max_shap_rows,
    max_model_rows = memory_estimate$max_model_rows,
    max_single_cell_rows = memory_estimate$max_single_cell_rows,
    exported_input_gib = memory_estimate$exported_input_gib,
    static_worker_estimate_gib = memory_estimate$static_worker_estimate_gib,
    health_worker_rss_gib = health_worker_rss_gib,
    calibrated_worker_estimate_gib = calibrated_worker_estimate_gib,
    reserved_ram_gib = selection$reserve_gib,
    memory_safe_workers = selection$memory_safe_workers,
    selected_workers = selection$selected_workers,
    backend = selection$backend,
    selection_reason = selection$selection_reason,
    health_check_passed = health_passed,
    health_result_retained = health_passed &&
      selection$backend == "psock_recycled" && !fallback_used,
    worker_recycling_interval = if (selection$backend == "psock_recycled") {
      recycle_every
    } else {
      NA_integer_
    },
    fallback_used = fallback_used,
    fallback_reason = fallback_reason,
    patient_plan_seconds = patient_plan_seconds,
    compute_seconds = compute_seconds,
    total_seconds = unname(proc.time()[["elapsed"]]) - started + patient_plan_seconds,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  list(values = values, runtime = runtime)
}
