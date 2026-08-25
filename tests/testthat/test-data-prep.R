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
})
