make_seurat_assay_compatibility_object <- function(type = "Assay5", split = FALSE) {
  counts <- Matrix::Matrix(
    matrix(
      (seq_len(6L * 12L) * 7L) %% 19L,
      nrow = 6L,
      dimnames = list(paste0("Gene", seq_len(6L)), paste0("Cell", seq_len(12L)))
    ),
    sparse = TRUE
  )
  assay_counts <- if (split) {
    list(batch_A = counts[, 1:6], batch_B = counts[, 7:12])
  } else {
    counts
  }
  assay <- switch(type,
    Assay = SeuratObject::CreateAssayObject(counts = counts),
    Assay5 = SeuratObject::CreateAssay5Object(counts = assay_counts),
    SCTAssay = methods::new("SCTAssay", SeuratObject::CreateAssayObject(counts = counts))
  )
  object <- SeuratObject::CreateSeuratObject(counts = assay)
  object$outcome <- rep(c("left", "right", "excluded"), c(5L, 5L, 2L))
  object$donor <- rep(paste0("Patient", seq_len(4L)), each = 3L)
  object$grouping <- rep(c("A", "B"), 6L)
  SeuratObject::VariableFeatures(object) <- c("Gene2", "Gene4", "Gene6")
  object
}

test_that("legacy and already-joined assays are left unchanged", {
  for (type in c("Assay", "SCTAssay", "Assay5")) {
    object <- make_seurat_assay_compatibility_object(type)
    expect_identical(join_seurat_layers_if_supported(object), object, info = type)
    expect_identical(join_seurat_layers_if_supported(object, "RNA"), object, info = type)
    expect_identical(class(object[["RNA"]])[[1]], type)
  }
})

test_that("Seurat preparation preserves the existing normalized values for each assay class", {
  for (type in c("Assay", "SCTAssay", "Assay5")) {
    object <- make_seurat_assay_compatibility_object(type)
    snapshot <- serialize(object, NULL)
    genes <- c("Gene6", "Gene2", "Gene4")
    expected_object <- Seurat::NormalizeData(
      object, assay = "RNA", normalization.method = "LogNormalize",
      scale.factor = 10000, verbose = FALSE
    )
    expected_object <- expected_object[, expected_object$outcome %in% c("left", "right")]
    expected_meta <- data.frame(
      condition = factor(expected_object$outcome),
      sample_id = expected_object$donor,
      grouping = expected_object$grouping
    )
    expected_matrix <- SeuratObject::GetAssayData(expected_object, assay = "RNA", layer = "data")
    expected <- data.frame(
      expected_meta, Matrix::t(expected_matrix[genes, , drop = FALSE]), check.names = FALSE
    )

    result <- build_dataset_from_seurat(
      object, condition_col = "outcome", sample_col = "donor",
      split_level_column = "grouping", assay = "RNA", genes = genes,
      condition_levels = c("left", "right")
    )

    expect_identical(result$dataset_all, expected, info = type)
    expect_identical(result$what_genes, genes)
    expect_identical(result$split_level_column, "grouping")
    expect_identical(serialize(object, NULL), snapshot)
  }
})

test_that("split v5 layers join automatically and match the previous working path", {
  object <- make_seurat_assay_compatibility_object(split = TRUE)
  snapshot <- serialize(object, NULL)
  expected <- SeuratObject::JoinLayers(object, assay = "RNA")
  joined <- join_seurat_layers_if_supported(object, assay = "RNA")

  expect_identical(joined, expected)
  expect_identical(SeuratObject::Layers(joined[["RNA"]]), "counts")
  expect_identical(SeuratObject::VariableFeatures(joined), SeuratObject::VariableFeatures(object))
  expect_identical(serialize(object, NULL), snapshot)

  arguments <- list(
    condition_col = "outcome", sample_col = "donor", split_level_column = "grouping",
    genes = c("Gene6", "Gene2", "Gene4"), condition_levels = c("left", "right")
  )
  actual <- do.call(build_dataset_from_seurat, c(list(seurat_obj = object), arguments))
  baseline <- do.call(build_dataset_from_seurat, c(list(seurat_obj = expected), arguments))
  expect_identical(actual, baseline)
})

test_that("normalized split layers are joined without changing their values", {
  object <- make_seurat_assay_compatibility_object(split = TRUE)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  expected <- SeuratObject::JoinLayers(object, assay = "RNA")
  joined <- join_seurat_layers_if_supported(object)

  expect_identical(joined, expected)
  expect_setequal(SeuratObject::Layers(joined[["RNA"]]), c("counts", "data"))
  expect_identical(
    SeuratObject::GetAssayData(joined, layer = "data"),
    SeuratObject::GetAssayData(expected, layer = "data")
  )
})

test_that("split data and scale.data layers are detected independently", {
  for (family in c("data", "scale.data")) {
    object <- make_seurat_assay_compatibility_object()
    assay <- object[["RNA"]]
    counts <- SeuratObject::GetAssayData(object, layer = "counts")
    SeuratObject::LayerData(assay, layer = paste0(family, ".batch_A")) <-
      as.matrix(counts[, 1:6])
    SeuratObject::LayerData(assay, layer = paste0(family, ".batch_B")) <-
      as.matrix(counts[, 7:12])
    object[["RNA"]] <- assay
    joined <- join_seurat_layers_if_supported(object)

    expect_identical(SeuratObject::GetAssayData(joined, layer = "counts"), counts)
    expect_setequal(SeuratObject::Layers(joined[["RNA"]]), c("counts", family))
    expected_values <- as.matrix(counts)
    expect_equal(
      as.matrix(SeuratObject::GetAssayData(joined, layer = family)), expected_values
    )
  }
})

test_that("inherited joining methods are recognized without changing assay class", {
  check_subclass <- function() {
    class_name <- "SHAPTURECompatibilityTestAssay5"
    methods::setClass(class_name, contains = "Assay5", where = globalenv())
    on.exit(methods::removeClass(class_name, where = globalenv()), add = TRUE)
    object <- make_seurat_assay_compatibility_object(split = TRUE)
    custom_assay <- methods::new(class_name, object[["RNA"]])
    object <- SeuratObject::CreateSeuratObject(counts = custom_assay)
    joined <- join_seurat_layers_if_supported(object)

    expect_identical(class(joined[["RNA"]])[[1]], class_name)
    expect_identical(SeuratObject::Layers(joined[["RNA"]]), "counts")
    expect_identical(joined, SeuratObject::JoinLayers(object, assay = "RNA"))
  }
  check_subclass()
})

test_that("joining targets the requested assay rather than the default assay", {
  object <- make_seurat_assay_compatibility_object("Assay")
  split_assay <- make_seurat_assay_compatibility_object(split = TRUE)[["RNA"]]
  SeuratObject::Key(split_assay) <- "adt_"
  object[["ADT"]] <- split_assay
  original_rna <- object[["RNA"]]
  original_metadata <- object[[]]
  joined <- join_seurat_layers_if_supported(object, assay = "ADT")

  expect_identical(SeuratObject::DefaultAssay(joined), "RNA")
  expect_identical(joined[["RNA"]], original_rna)
  expect_identical(joined[[]], original_metadata)
  expect_identical(SeuratObject::Layers(joined[["ADT"]]), "counts")
})

test_that("RNA SCT ADT and custom assay names use the same preparation path", {
  original <- make_seurat_assay_compatibility_object(split = TRUE)
  arguments <- list(condition_col = "outcome", sample_col = "donor", genes = c("Gene2", "Gene4"))
  baseline <- do.call(build_dataset_from_seurat, c(list(seurat_obj = original), arguments))

  for (assay_name in c("RNA", "SCT", "ADT", "CustomProtein")) {
    assay <- original[["RNA"]]
    SeuratObject::Key(assay) <- paste0(tolower(assay_name), "_")
    object <- SeuratObject::CreateSeuratObject(
      counts = assay, assay = assay_name, meta.data = original[[]]
    )
    joined <- join_seurat_layers_if_supported(object, assay = assay_name)
    result <- do.call(
      build_dataset_from_seurat,
      c(list(seurat_obj = object, assay = assay_name), arguments)
    )

    expect_identical(SeuratObject::Layers(joined[[assay_name]]), "counts")
    expect_identical(SeuratObject::DefaultAssay(joined), assay_name)
    expect_identical(result, baseline)
  }
})

test_that("a single remaining split layer is renamed by native joining", {
  object <- make_seurat_assay_compatibility_object()
  counts <- SeuratObject::GetAssayData(object, layer = "counts")
  split_object <- SeuratObject::CreateSeuratObject(
    counts = SeuratObject::CreateAssay5Object(counts = list(batch_A = counts))
  )
  joined <- join_seurat_layers_if_supported(split_object)

  expect_identical(SeuratObject::Layers(joined[["RNA"]]), "counts")
  expect_identical(SeuratObject::GetAssayData(joined, layer = "counts"), counts)
})

test_that("joining falls back unchanged when the generic is unavailable", {
  object <- make_seurat_assay_compatibility_object(split = TRUE)
  result <- testthat::with_mocked_bindings(
    join_seurat_layers_if_supported(object),
    JoinLayers = NULL, .package = "SeuratObject"
  )
  expect_identical(result, object)
})

test_that("actual layer-joining errors are not swallowed", {
  object <- make_seurat_assay_compatibility_object(split = TRUE)
  # An invalid layer must retain Seurat's native failure rather than be skipped.
  object@assays$RNA@layers[["counts.batch_B"]] <- "invalid matrix"
  native_error <- tryCatch(
    SeuratObject::JoinLayers(object, assay = "RNA"), error = identity
  )
  expect_s3_class(native_error, "error")
  expect_error(
    join_seurat_layers_if_supported(object), conditionMessage(native_error), fixed = TRUE
  )
})

test_that("assay selection errors are reported clearly", {
  object <- make_seurat_assay_compatibility_object()
  expect_error(join_seurat_layers_if_supported(list()), "must be a Seurat object")
  for (assay in list("missing", "", NA_character_, character(), c("RNA", "ADT"), 1L)) {
    expect_error(join_seurat_layers_if_supported(object, assay), "assay must name one assay")
  }
})
