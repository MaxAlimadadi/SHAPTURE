make_condition_label_shapviz <- function(n_features) {
  features <- paste0("feature_", seq_len(n_features))
  shap <- matrix(
    seq_len(n_features) * rep(c(-1, 1), length.out = n_features) / 100,
    nrow = 1,
    dimnames = list(NULL, features)
  )
  values <- as.data.frame(
    matrix(seq_len(n_features) / n_features, nrow = 1,
           dimnames = list(NULL, features))
  )
  shapviz::shapviz(shap, X = values, baseline = 0.25)
}

test_that("waterfalls keep all displayed rows independently of selective labels", {
  for (n_features in c(10L, 50L, 60L)) {
    object <- make_condition_label_shapviz(n_features)
    original <- object
    displayed_count <- min(50L, n_features)

    for (label_count in c(0L, 15L, 50L)) {
      plot <- build_shapture_condition_waterfall_plot(
        object, "condition_A", max_display = 50L,
        top_feature_label = label_count
      )
      built <- ggplot2::ggplot_build(plot)

      expect_equal(nrow(plot$data), displayed_count)
      expect_true(all(nzchar(plot$data$label)))
      expect_equal(length(unique(built$data[[1]]$y)), displayed_count)
      expect_equal(sum(plot$data$S), sum(object$S))
      expect_equal(tail(plot$data$to, 1), object$baseline + sum(object$S))
      expect_identical(object, original)

      features <- extract_shapture_shapviz_plot_features(
        plot, colnames(object$X), keep_alignment = TRUE
      )
      feature_rows <- which(!is.na(features))
      expected_labels <- paste0(
        features[feature_rows], " (",
        round(as.numeric(object$X[1, features[feature_rows], drop = TRUE]), 2),
        ")"
      )
      expect_equal(plot$data$label[feature_rows], expected_labels)
      expect_equal(sum(grepl("other features", plot$data$label)),
                   as.integer(n_features > 50L))
    }
  }
})

test_that("waterfall defaults label every selected feature and render without warnings", {
  object <- make_condition_label_shapviz(60L)
  for (display_count in c(1L, 5L, 50L, 100L)) {
    plot <- build_shapture_condition_waterfall_plot(object, "condition_A", display_count)
    expect_equal(nrow(plot$data), min(display_count, ncol(object$S)))
    expect_true(all(nzchar(plot$data$label)))
    expect_warning(ggplot2::ggplotGrob(plot), NA)
  }
})

test_that("force plots retain selective top-feature labels and other-feature summaries", {
  for (n_features in c(10L, 50L, 60L)) {
    object <- make_condition_label_shapviz(n_features)
    for (label_count in c(0L, 15L, 50L)) {
      plot <- build_shapture_condition_force_plot(
        object, "condition_A", max_display = 50L,
        top_feature_label = label_count
      )
      features <- extract_shapture_shapviz_plot_features(
        plot, colnames(object$X), keep_alignment = TRUE
      )
      feature_rows <- which(!is.na(features))
      labeled_rows <- feature_rows[nzchar(plot$data$label[feature_rows])]
      expected_rows <- head(
        feature_rows[order(abs(plot$data$S[feature_rows]), decreasing = TRUE)],
        label_count
      )

      expect_equal(nrow(plot$data), min(50L, n_features))
      expect_setequal(labeled_rows, expected_rows)
      expect_length(labeled_rows, min(label_count, length(feature_rows)))
      expect_true(all(nzchar(plot$data$label[is.na(features)])))
      expect_equal(sum(plot$data$S), sum(object$S))
    }
    expect_warning(ggplot2::ggplotGrob(plot), NA)
  }
})

test_that("waterfall average-expression exports include every displayed feature", {
  for (n_features in c(10L, 60L)) {
    object <- make_condition_label_shapviz(n_features)
    plot <- build_shapture_condition_waterfall_plot(object, "condition_A", 50L, 15L)
    features <- extract_shapture_shapviz_plot_features(plot, colnames(object$X))
    output <- tempfile(fileext = ".csv")
    write_shapture_condition_average_expression(object, plot, output)
    values <- utils::read.csv(output, row.names = 1)$average

    expect_length(values, min(n_features, 49L))
    expect_equal(values, rev(round(as.numeric(object$X[1, features, drop = TRUE]), 2)))
    unlink(output)
  }
})
