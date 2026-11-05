#' Build SHAPTURE input from a Seurat object
#'
#' @param seurat_obj A Seurat object.
#' @param condition_col Metadata column containing the binary condition labels.
#' @param sample_col Metadata column containing sample identifiers.
#' @param split_level_column Optional metadata column to carry for subset-specific or
#'   one-vs-others analyses.
#' @param assay Assay to read.
#' @param layer Assay layer to read. Defaults to `"data"`.
#' @param genes Optional feature list. If `NULL`, all features in the requested
#'   assay/layer are used.
#' @param condition_levels Optional levels to retain from `condition_col`.
#'
#' @return A list with `dataset_all`, `what_genes`, the prepared-data
#'   `split_level_column` name, and a `column_mapping` table. Source metadata names
#'   may be arbitrary. If the requested split name conflicts with SHAPTURE's
#'   prepared-data columns (`condition` or `sample_id`) or a feature name, the
#'   split is stored under a deterministic internal alias returned in
#'   `split_level_column`.
#' @export
build_dataset_from_seurat <- function(
    seurat_obj,
    condition_col,
    sample_col,
    split_level_column = NULL,
    assay = "RNA",
    layer = "data",
    genes = NULL,
    condition_levels = NULL
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  
  meta_cols <- c(condition_col, sample_col, split_level_column)
  meta_cols <- unique(meta_cols[!is.na(meta_cols) & nzchar(meta_cols)])
  missing_cols <- setdiff(meta_cols, colnames(seurat_obj@meta.data))
  if (length(missing_cols) > 0) {
    stop("Missing metadata columns in Seurat object: ", paste(missing_cols, collapse = ", "))
  }
  
  Seurat::DefaultAssay(seurat_obj) <- assay
  if (exists("JoinLayers", where = asNamespace("SeuratObject"), inherits = FALSE)) {
    seurat_obj <- SeuratObject::JoinLayers(seurat_obj)
  }
  seurat_obj <- Seurat::NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  
  if (is.null(condition_levels)) {
    condition_levels <- unique(as.character(seurat_obj@meta.data[[condition_col]]))
  }
  
  seurat_obj <- seurat_obj[, seurat_obj@meta.data[[condition_col]] %in% condition_levels]
  seurat_obj@meta.data[[condition_col]] <- factor(seurat_obj@meta.data[[condition_col]])
  
  assay_mat <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  
  if (is.null(genes)) {
    genes <- rownames(assay_mat)
  } else {
    genes <- intersect(as.character(genes), rownames(assay_mat))
  }
  
  if (length(genes) == 0) {
    stop("No requested genes were found in the selected Seurat assay/layer.")
  }

  prepared_split_level_column <- resolve_prepared_split_column(
    split_level_column = split_level_column,
    feature_names = genes
  )
  
  meta_df <- data.frame(
    condition = seurat_obj@meta.data[[condition_col]],
    sample_id = seurat_obj@meta.data[[sample_col]],
    check.names = FALSE
  )
  
  if (!is.null(prepared_split_level_column)) {
    meta_df[[prepared_split_level_column]] <- seurat_obj@meta.data[[split_level_column]]
  }
  
  dataset_all <- data.frame(
    meta_df,
    Matrix::t(assay_mat[genes, , drop = FALSE]),
    check.names = FALSE
  )

  column_mapping <- build_prepared_column_mapping(
    condition_col = condition_col,
    sample_col = sample_col,
    split_level_column = split_level_column,
    prepared_split_level_column = prepared_split_level_column
  )
  if (!identical(split_level_column, prepared_split_level_column)) {
    attr(dataset_all, "shapture_column_mapping") <- column_mapping
  }
  
  list(
    dataset_all = dataset_all,
    what_genes = genes,
    split_level_column = prepared_split_level_column,
    column_mapping = column_mapping
  )
}


#' Build SHAPTURE input from a SingleCellExperiment / CATALYST object
#'
#' @param sce A SingleCellExperiment object.
#' @param condition_col Column in `colData(sce)` containing condition labels.
#' @param sample_col Column in `colData(sce)` containing sample identifiers.
#' @param split_level_column Optional split column. If `NULL`, no subset column is
#'   added. If `cluster_resolution` is supplied, the resulting column will use
#'   the name stored in `split_level_column`.
#' @param cluster_resolution Optional CATALYST clustering resolution key used
#'   with `cluster_ids()`.
#' @param feature_cols Optional feature list. If `NULL`, CATALYST type/state
#'   markers are used.
#'
#' @return A list with `dataset_all`, `what_genes`, the prepared-data
#'   `split_level_column` name, and a `column_mapping` table. Source metadata names
#'   may be arbitrary. If the requested split name conflicts with SHAPTURE's
#'   prepared-data columns (`condition` or `sample_id`) or a feature name, the
#'   split is stored under a deterministic internal alias returned in
#'   `split_level_column`.
#' @export
build_dataset_from_sce <- function(
    sce,
    condition_col,
    sample_col,
    split_level_column = NULL,
    cluster_resolution = NULL,
    feature_cols = NULL
) {
  if (!methods::is(sce, "SingleCellExperiment")) {
    stop("sce must be a SingleCellExperiment.")
  }
  
  meta_cols <- c(condition_col, sample_col)
  missing_cols <- setdiff(meta_cols, colnames(SummarizedExperiment::colData(sce)))
  if (length(missing_cols) > 0) {
    stop("Missing metadata columns in sce: ", paste(missing_cols, collapse = ", "))
  }
  
  if (is.null(feature_cols)) {
    feature_cols <- c(CATALYST::type_markers(sce), CATALYST::state_markers(sce))
  }

  prepared_split_level_column <- resolve_prepared_split_column(
    split_level_column = split_level_column,
    feature_names = feature_cols
  )
  
  meta_df <- data.frame(
    condition = SummarizedExperiment::colData(sce)[[condition_col]],
    sample_id = SummarizedExperiment::colData(sce)[[sample_col]],
    check.names = FALSE
  )
  
  if (!is.null(prepared_split_level_column)) {
    if (!is.null(cluster_resolution)) {
      meta_df[[prepared_split_level_column]] <- CATALYST::cluster_ids(sce, k = cluster_resolution)
    } else {
      if (!split_level_column %in% colnames(SummarizedExperiment::colData(sce))) {
        stop("split_level_column is not present in colData(sce) and no cluster_resolution was supplied.")
      }
      meta_df[[prepared_split_level_column]] <- SummarizedExperiment::colData(sce)[[split_level_column]]
    }
  }
  
  dataset_all <- data.frame(
    meta_df,
    t(sce@assays@data$exprs[feature_cols, , drop = FALSE]),
    check.names = FALSE
  )

  column_mapping <- build_prepared_column_mapping(
    condition_col = condition_col,
    sample_col = sample_col,
    split_level_column = split_level_column,
    prepared_split_level_column = prepared_split_level_column
  )
  if (!identical(split_level_column, prepared_split_level_column)) {
    attr(dataset_all, "shapture_column_mapping") <- column_mapping
  }
  
  list(
    dataset_all = dataset_all,
    what_genes = feature_cols,
    split_level_column = prepared_split_level_column,
    column_mapping = column_mapping
  )
}


#' Resolve a collision-safe split column for prepared data
#'
#' @param split_level_column Source split column name or `NULL`.
#' @param feature_names Feature names that will be added to the prepared table.
#'
#' @return The prepared-data split column name or `NULL`.
#' @keywords internal
resolve_prepared_split_column <- function(split_level_column, feature_names = character()) {
  if (is.null(split_level_column)) {
    return(NULL)
  }

  if (length(split_level_column) != 1L || is.na(split_level_column) || !nzchar(split_level_column)) {
    stop("split_level_column must be NULL or one non-empty column name.")
  }

  occupied_names <- unique(c("condition", "sample_id", as.character(feature_names)))
  if (!split_level_column %in% occupied_names) {
    return(split_level_column)
  }

  alias_base <- paste0(".shapture_split_", make.names(split_level_column))
  alias <- alias_base
  suffix <- 2L
  while (alias %in% occupied_names) {
    alias <- paste0(alias_base, "_", suffix)
    suffix <- suffix + 1L
  }

  alias
}


#' Record source-to-prepared metadata column names
#'
#' @param condition_col Source condition column.
#' @param sample_col Source sample column.
#' @param split_level_column Optional source split column.
#' @param prepared_split_level_column Optional prepared-data split column.
#'
#' @return A data frame describing the adapter's metadata-column mapping.
#' @keywords internal
build_prepared_column_mapping <- function(
    condition_col,
    sample_col,
    split_level_column = NULL,
    prepared_split_level_column = NULL
) {
  roles <- c("condition", "sample_id")
  source_columns <- c(condition_col, sample_col)
  prepared_columns <- c("condition", "sample_id")

  if (!is.null(split_level_column)) {
    roles <- c(roles, "split_level_column")
    source_columns <- c(source_columns, split_level_column)
    prepared_columns <- c(prepared_columns, prepared_split_level_column)
  }

  data.frame(
    role = roles,
    source_column = source_columns,
    prepared_column = prepared_columns,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


#' Resolve a source split name recorded by an object adapter
#'
#' @param column_name Configured split column name or `NULL`.
#' @param dataset_all Prepared SHAPTURE data.
#'
#' @return The corresponding prepared-data column name, when a collision map
#'   is available, or `column_name` unchanged.
#' @keywords internal
resolve_mapped_prepared_column <- function(column_name, dataset_all) {
  if (is.null(column_name)) {
    return(NULL)
  }

  column_mapping <- attr(dataset_all, "shapture_column_mapping", exact = TRUE)
  required_mapping_cols <- c("role", "source_column", "prepared_column")
  if (!is.data.frame(column_mapping) ||
      !all(required_mapping_cols %in% colnames(column_mapping))) {
    return(column_name)
  }

  mapped_rows <- column_mapping$role == "split_level_column" &
    column_mapping$source_column == column_name
  mapped_names <- unique(column_mapping$prepared_column[mapped_rows])
  mapped_names <- mapped_names[!is.na(mapped_names) & nzchar(mapped_names)]

  if (length(mapped_names) == 1L && mapped_names %in% colnames(dataset_all)) {
    return(mapped_names)
  }

  column_name
}
