make_cross_patient_test_data <- function(patients_per_condition = 3L) {
  patient_ids <- c(
    paste0("c0_p", seq_len(patients_per_condition)),
    paste0("c1_p", seq_len(patients_per_condition))
  )
  conditions <- rep(c(0, 1), each = patients_per_condition)
  metadata <- do.call(rbind, lapply(seq_along(patient_ids), function(i) {
    data.frame(
      sample_id = patient_ids[i],
      condition = conditions[i],
      cell_index = seq_len(6),
      stringsAsFactors = FALSE
    )
  }))
  row.names(metadata) <- NULL

  set.seed(101)
  f1 <- pmax(0, metadata$condition * 2 + stats::rnorm(nrow(metadata), sd = 0.2))
  f2 <- pmax(0, (1 - metadata$condition) * 1.5 + stats::rnorm(nrow(metadata), sd = 0.25))
  expression <- data.frame(f1 = f1, f2 = f2)
  shap <- cbind(
    f1 = scale(f1)[, 1] + stats::rnorm(length(f1), sd = 0.05),
    f2 = scale(f2)[, 1] + stats::rnorm(length(f2), sd = 0.05)
  )
  model_data <- cbind(metadata[, c("sample_id", "condition")], expression)
  consistency <- compute_within_condition_expression_consistency(
    expression_data = expression,
    condition = metadata$condition,
    condition_levels = c(0, 1)
  )
  means0 <- colMeans(expression[metadata$condition == 0, , drop = FALSE])
  means1 <- colMeans(expression[metadata$condition == 1, , drop = FALSE])
  pseudocount <- compute_dynamic_pseudocount(as.matrix(expression))
  log2fc <- log2((means1 + pseudocount) / (means0 + pseudocount))
  det0 <- colMeans(expression[metadata$condition == 0, , drop = FALSE] > 0)
  det1 <- colMeans(expression[metadata$condition == 1, , drop = FALSE] > 0)
  shap_cor <- vapply(colnames(expression), function(feature) {
    stats::cor(expression[[feature]], shap[, feature], method = "spearman")
  }, numeric(1))

  full_stats <- data.frame(
    feature = colnames(expression),
    expression_log2FC = log2fc,
    det_diff = det1 - det0,
    shap_expr_cor = shap_cor,
    within_condition_consistency_raw = consistency$within_condition_consistency_raw,
    check.names = FALSE
  )
  average_importance <- data.frame(
    feature = colnames(shap),
    `average|SHAP value|` = colMeans(abs(shap)),
    check.names = FALSE
  )

  list(
    shap = shap,
    expression = expression,
    metadata = metadata,
    model_data = model_data,
    single_cell_data = model_data,
    full_stats = full_stats,
    average_importance = average_importance
  )
}
