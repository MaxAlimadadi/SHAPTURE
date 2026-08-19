#' SHAPTURE: Sample-Aware SHAP-Based Feature Discovery
#'
#' SHAPTURE provides a sample-aware workflow for subset-specific and
#' one-vs-others feature discovery from single-cell or metacell data using
#' repeated machine-learning fits, SHAP-based interpretation, patient-aware
#' balancing, and integrated marker ranking.
#'
#' @docType package
#' @name SHAPTURE
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @importFrom grDevices dev.off pdf png
#' @importFrom rlang .data
#' @importFrom stats aggregate predict setNames
#' @importFrom utils capture.output head str write.csv
"_PACKAGE"
