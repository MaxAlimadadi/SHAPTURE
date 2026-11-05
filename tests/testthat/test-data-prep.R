test_that("Seurat data preparation transposes sparse assay layers", {
  counts <- Matrix::Matrix(
    matrix(
      c(1, 0, 2, 3, 0, 4),
      nrow = 2,
      dimnames = list(c("geneA", "geneB"), paste0("cell", 1:3))
    ),
    sparse = TRUE
  )
  object <- Seurat::CreateSeuratObject(counts = counts)
  object$condition <- c("left", "left", "right")
  object$sample_id <- c("patient_1", "patient_1", "patient_2")

  result <- build_dataset_from_seurat(
    seurat_obj = object,
    condition_col = "condition",
    sample_col = "sample_id",
    genes = c("geneA", "geneB"),
    condition_levels = c("left", "right")
  )

  expect_equal(nrow(result$dataset_all), 3L)
  expect_equal(result$what_genes, c("geneA", "geneB"))
  expect_true(all(c("condition", "sample_id", "geneA", "geneB") %in%
    colnames(result$dataset_all)))
  expect_true(all(vapply(result$dataset_all[c("geneA", "geneB")], is.numeric, logical(1))))
  expect_null(result$split_level_column)
  expect_equal(
    result$column_mapping$prepared_column,
    c("condition", "sample_id")
  )
})


test_that("Seurat adapter preserves nonconflicting split names", {
  counts <- Matrix::Matrix(
    matrix(
      c(1, 0, 2, 3, 0, 4),
      nrow = 2,
      dimnames = list(c("geneA", "geneB"), paste0("cell", 1:3))
    ),
    sparse = TRUE
  )
  object <- Seurat::CreateSeuratObject(counts = counts)
  object$outcome <- c("left", "left", "right")
  object$donor <- c("patient_1", "patient_1", "patient_2")
  object$subset_label <- c("A", "B", "B")

  result <- build_dataset_from_seurat(
    seurat_obj = object,
    condition_col = "outcome",
    sample_col = "donor",
    split_level_column = "subset_label",
    genes = c("geneA", "geneB"),
    condition_levels = c("left", "right")
  )

  expect_identical(result$split_level_column, "subset_label")
  expect_identical(
    colnames(result$dataset_all),
    c("condition", "sample_id", "subset_label", "geneA", "geneB")
  )
  expect_equal(as.character(result$dataset_all$condition), c("left", "left", "right"))
  expect_equal(as.character(result$dataset_all$sample_id), c("patient_1", "patient_1", "patient_2"))
  expect_equal(as.character(result$dataset_all$subset_label), c("A", "B", "B"))
})


test_that("Seurat adapter aliases split names that collide with prepared metadata", {
  counts <- Matrix::Matrix(
    matrix(
      c(1, 0, 2, 3, 0, 4),
      nrow = 2,
      dimnames = list(c("geneA", "geneB"), paste0("cell", 1:3))
    ),
    sparse = TRUE
  )
  object <- Seurat::CreateSeuratObject(counts = counts)
  object$outcome <- c("left", "left", "right")
  object$donor <- c("patient_1", "patient_1", "patient_2")
  object$condition <- c("subset_A", "subset_B", "subset_B")

  result <- build_dataset_from_seurat(
    seurat_obj = object,
    condition_col = "outcome",
    sample_col = "donor",
    split_level_column = "condition",
    genes = c("geneA", "geneB"),
    condition_levels = c("left", "right")
  )

  expect_identical(result$split_level_column, ".shapture_split_condition")
  expect_false(anyDuplicated(colnames(result$dataset_all)) > 0L)
  expect_equal(as.character(result$dataset_all$condition), c("left", "left", "right"))
  expect_equal(
    as.character(result$dataset_all[[result$split_level_column]]),
    c("subset_A", "subset_B", "subset_B")
  )
  expect_equal(
    result$column_mapping$prepared_column,
    c("condition", "sample_id", ".shapture_split_condition")
  )
  expect_identical(
    resolve_mapped_prepared_column("condition", result$dataset_all),
    ".shapture_split_condition"
  )
})


test_that("SCE adapter aliases split names that collide with prepared metadata", {
  exprs <- matrix(
    c(1, 0, 2, 3, 0, 4),
    nrow = 2,
    dimnames = list(c("markerA", "markerB"), paste0("cell", 1:3))
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(exprs = exprs),
    colData = data.frame(
      outcome = c("left", "left", "right"),
      donor = c("patient_1", "patient_1", "patient_2"),
      condition = c("subset_A", "subset_B", "subset_B")
    )
  )

  result <- build_dataset_from_sce(
    sce = sce,
    condition_col = "outcome",
    sample_col = "donor",
    split_level_column = "condition",
    cluster_resolution = NULL,
    feature_cols = c("markerA", "markerB")
  )

  expect_identical(result$split_level_column, ".shapture_split_condition")
  expect_false(anyDuplicated(colnames(result$dataset_all)) > 0L)
  expect_equal(as.character(result$dataset_all$condition), c("left", "left", "right"))
  expect_equal(as.character(result$dataset_all$sample_id), c("patient_1", "patient_1", "patient_2"))
  expect_equal(
    as.character(result$dataset_all[[result$split_level_column]]),
    c("subset_A", "subset_B", "subset_B")
  )
})


test_that("SCE adapter leaves nonconflicting prepared data unchanged", {
  exprs <- matrix(
    c(1, 0, 2, 3, 0, 4),
    nrow = 2,
    dimnames = list(c("markerA", "markerB"), paste0("cell", 1:3))
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(exprs = exprs),
    colData = data.frame(
      outcome = c("left", "left", "right"),
      donor = c("patient_1", "patient_1", "patient_2"),
      subset_label = c("A", "B", "B")
    )
  )

  result <- build_dataset_from_sce(
    sce = sce,
    condition_col = "outcome",
    sample_col = "donor",
    split_level_column = "subset_label",
    cluster_resolution = NULL,
    feature_cols = c("markerA", "markerB")
  )
  expected <- data.frame(
    condition = SummarizedExperiment::colData(sce)[["outcome"]],
    sample_id = SummarizedExperiment::colData(sce)[["donor"]],
    subset_label = SummarizedExperiment::colData(sce)[["subset_label"]],
    t(exprs[c("markerA", "markerB"), , drop = FALSE]),
    check.names = FALSE
  )

  expect_identical(result$split_level_column, "subset_label")
  expect_identical(result$dataset_all, expected)
})


test_that("script config resolves source split names from collision metadata", {
  dataset_all <- data.frame(
    condition = factor(c("left", "right", "left", "right")),
    sample_id = c("P1", "P2", "P3", "P4"),
    .shapture_split_condition = c("A", "A", "B", "B"),
    markerA = c(0, 1, 2, 3),
    check.names = FALSE
  )
  attr(dataset_all, "shapture_column_mapping") <- data.frame(
    role = c("condition", "sample_id", "split_level_column"),
    source_column = c("outcome", "donor", "condition"),
    prepared_column = c("condition", "sample_id", ".shapture_split_condition"),
    stringsAsFactors = FALSE
  )
  config_env <- list2env(list(
    condition_levels = c("left", "right"),
    shapture_config = list(
      split_level_column = "condition",
      comparison_mode = "condition_vs_condition",
      what_0 = "left",
      what_1 = "right"
    )
  ))

  resolved <- resolve_shapture_script_config(
    envir = config_env,
    dataset_all = dataset_all,
    what_genes = "markerA",
    main_dir = tempdir()
  )

  expect_identical(resolved$split_level_column, ".shapture_split_condition")
  expect_identical(
    as.character(resolved$dataset_all[[resolved$split_level_column]]),
    c("A", "A", "B", "B")
  )
})


test_that("split aliases avoid feature-name collisions deterministically", {
  expect_identical(
    resolve_prepared_split_column(
      split_level_column = "condition",
      feature_names = c("geneA", ".shapture_split_condition")
    ),
    ".shapture_split_condition_2"
  )
  expect_identical(
    resolve_prepared_split_column(
      split_level_column = "sample_id",
      feature_names = "geneA"
    ),
    ".shapture_split_sample_id"
  )
})
