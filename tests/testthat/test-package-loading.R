test_that("package workflows expose their functions without a separate loader", {
  workflow_functions <- c(
    "build_dataset_from_seurat", "build_dataset_from_sce",
    "resolve_shapture_script_config", "run_shapture_pipeline",
    "load_shapture_visualization_inputs", "load_shapture_global_results"
  )
  for (name in workflow_functions) {
    expect_true(is.function(getExportedValue("SHAPTURE", name)), info = name)
  }
})

test_that("parallel support is available inside the package namespace", {
  namespace <- asNamespace("SHAPTURE")
  helpers <- c(
    "safe_stop_shapture_cluster", "select_cross_patient_workers",
    "run_cross_patient_bootstrap_controller"
  )
  for (name in helpers) {
    expect_true(exists(name, envir = namespace, inherits = FALSE), info = name)
    expect_true(is.function(get(name, envir = namespace)), info = name)
  }
})
