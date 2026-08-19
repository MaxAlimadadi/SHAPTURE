#' Build SHAPTURE input from a Seurat object
#'
#' @param seurat_obj A Seurat object.
#' @param condition_col Metadata column containing the binary condition labels.
#' @param sample_col Metadata column containing sample identifiers.
#' @param split_by_col Optional metadata column to carry for subset-specific or
#'   one-vs-others analyses.
#' @param assay Assay to read.
#' @param layer Assay layer to read. Defaults to `"data"`.
#' @param genes Optional feature list. If `NULL`, all features in the requested
#'   assay/layer are used.
#' @param condition_levels Optional levels to retain from `condition_col`.
#'
#' @return A list with `dataset_all` and `what_genes`.
#' @export
build_dataset_from_seurat <- function(
    seurat_obj,
    condition_col,
    sample_col,
    split_by_col = NULL,
    assay = "RNA",
    layer = "data",
    genes = NULL,
    condition_levels = NULL
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  
  meta_cols <- c(condition_col, sample_col, split_by_col)
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
  
  meta_df <- data.frame(
    condition = seurat_obj@meta.data[[condition_col]],
    sample_id = seurat_obj@meta.data[[sample_col]],
    check.names = FALSE
  )
  
  if (!is.null(split_by_col)) {
    meta_df[[split_by_col]] <- seurat_obj@meta.data[[split_by_col]]
  }
  
  dataset_all <- data.frame(
    meta_df,
    t(assay_mat[genes, , drop = FALSE]),
    check.names = FALSE
  )
  
  list(
    dataset_all = dataset_all,
    what_genes = genes
  )
}


#' Build SHAPTURE input from a SingleCellExperiment / CATALYST object
#'
#' @param sce A SingleCellExperiment object.
#' @param condition_col Column in `colData(sce)` containing condition labels.
#' @param sample_col Column in `colData(sce)` containing sample identifiers.
#' @param split_by_col Optional split column. If `NULL`, no subset column is
#'   added. If `cluster_resolution` is supplied, the resulting column will use
#'   the name stored in `split_by_col`.
#' @param cluster_resolution Optional CATALYST clustering resolution key used
#'   with `cluster_ids()`.
#' @param feature_cols Optional feature list. If `NULL`, CATALYST type/state
#'   markers are used.
#'
#' @return A list with `dataset_all` and `what_genes`.
#' @export
build_dataset_from_sce <- function(
    sce,
    condition_col,
    sample_col,
    split_by_col = NULL,
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
  
  meta_df <- data.frame(
    condition = SummarizedExperiment::colData(sce)[[condition_col]],
    sample_id = SummarizedExperiment::colData(sce)[[sample_col]],
    check.names = FALSE
  )
  
  if (!is.null(split_by_col)) {
    if (!is.null(cluster_resolution)) {
      meta_df[[split_by_col]] <- CATALYST::cluster_ids(sce, k = cluster_resolution)
    } else {
      if (!split_by_col %in% colnames(SummarizedExperiment::colData(sce))) {
        stop("split_by_col is not present in colData(sce) and no cluster_resolution was supplied.")
      }
      meta_df[[split_by_col]] <- SummarizedExperiment::colData(sce)[[split_by_col]]
    }
  }
  
  dataset_all <- data.frame(
    meta_df,
    t(sce@assays@data$exprs[feature_cols, , drop = FALSE]),
    check.names = FALSE
  )
  
  list(
    dataset_all = dataset_all,
    what_genes = feature_cols
  )
}
