#####
# Run SHAPTURE analysis
#
# Before running this file, create these objects in your R session:
# - dataset_all
# - what_genes
# - shapture_config
#####

library(SHAPTURE)

# Base output directory. A non-NULL project_name creates a parent folder here.
analysis_dir <- getwd()
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("dataset_all")) {
  stop("dataset_all must exist before running SHAPTURE.")
}

if (!exists("what_genes")) {
  stop("what_genes must exist before running SHAPTURE.")
}

if (!exists("shapture_config")) {
  stop("shapture_config must exist before running SHAPTURE.")
}

pipeline_config <- resolve_shapture_script_config(
  envir = globalenv(),
  dataset_all = dataset_all,
  what_genes = what_genes,
  main_dir = analysis_dir
)

analysis_result <- do.call(run_shapture_pipeline, pipeline_config)
message("Results saved in ", analysis_result$project_dir)
