#' Validate a project folder name
#'
#' @param project_name A single folder name, or `NULL` for no project folder.
#' @return The unchanged valid project name.
#' @keywords internal
validate_shapture_project_name <- function(project_name) {
  if (is.null(project_name)) {
    return(NULL)
  }
  if (!is.character(project_name) || length(project_name) != 1L ||
      is.na(project_name) || !nzchar(trimws(project_name)) ||
      project_name %in% c(".", "..") ||
      grepl("[/\\\\]", project_name) || grepl("[[:cntrl:]]", project_name)) {
    stop("project_name must be NULL or one non-empty folder name, not a path.")
  }
  project_name
}


#' Prepare the project output directory
#'
#' @param main_dir Base output directory.
#' @param project_name Optional single project-folder name.
#' @param comparison_mode Comparison mode.
#' @param what_0 Negative condition label.
#' @param what_1 Positive condition label.
#' @return A list with the base directory, project directory, and standard
#'   comparison-folder prefix. Returned paths are absolute.
#' @keywords internal
prepare_shapture_output_layout <- function(
    main_dir, project_name, comparison_mode, what_0, what_1
) {
  project_name <- validate_shapture_project_name(project_name)
  if (!is.character(main_dir) || length(main_dir) != 1L || is.na(main_dir) ||
      !nzchar(trimws(main_dir))) {
    stop("main_dir must be one non-empty directory path.")
  }
  main_dir <- path.expand(main_dir)
  if (!dir.exists(main_dir) && !dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create the base output directory: ", main_dir)
  }
  main_dir <- normalizePath(main_dir, winslash = "/", mustWork = TRUE)
  project_dir <- if (is.null(project_name)) main_dir else file.path(main_dir, project_name)
  if (!dir.exists(project_dir) && !dir.create(project_dir, showWarnings = FALSE)) {
    stop("Cannot create the project output directory: ", project_dir)
  }
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  comparison_prefix <- if (comparison_mode == "condition_vs_condition") {
    paste0(what_0, "_", what_1)
  } else {
    "one_vs_others"
  }
  list(main_dir = main_dir, project_dir = project_dir, comparison_prefix = comparison_prefix)
}
