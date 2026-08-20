source_shapture_modules <- function() {
  this_file <- tryCatch(
    normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
    error = function(e) normalizePath("pipeline/helper_function.R", winslash = "/", mustWork = TRUE)
  )
  
  repo_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
  module_files <- c(
    "R/shapture_processing.R",
    "R/shapture_inference.R",
    "R/shapture_metacells.R",
    "R/shapture_data_prep.R",
    "R/shapture_scoring.R",
    "R/shapture_visualization.R",
    "R/shapture_post_analysis_visualization.R",
    "R/shapture_global_visualization.R",
    "R/shapture_pipeline.R"
  )
  
  for (module_file in module_files) {
    source(file.path(repo_root, module_file), local = globalenv())
  }
  
  required_operator_packages <- c("magrittr", "foreach")
  missing_operator_packages <- required_operator_packages[
    !vapply(required_operator_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_operator_packages) > 0) {
    stop(
      "Missing packages required for sourced SHAPTURE workflows: ",
      paste(missing_operator_packages, collapse = ", ")
    )
  }
  
  assign("%>%", get("%>%", envir = asNamespace("magrittr")), envir = globalenv())
  assign("%dopar%", get("%dopar%", envir = asNamespace("foreach")), envir = globalenv())
  
  invisible(TRUE)
}


source_shapture_modules()
