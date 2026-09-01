make_feature_label_test_data <- function(what_0 = "control", what_1 = "case") {
  data.frame(
    feature = c("gene<&>", "gene_B", "gene_C"),
    final_score = c(0.85, 0.7, 0.4),
    Model_expression_pattern = c(
      paste("Higher feature values push toward", what_1),
      paste("Higher feature values push toward", what_0),
      "Insufficient directional evidence"
    ),
    stringsAsFactors = FALSE
  )
}


test_that("Feature-importance labels match beeswarm colors in both comparison modes", {
  for (condition_labels in list(c("control", "case"), c("others", "target"))) {
    df <- make_feature_label_test_data(condition_labels[1], condition_labels[2])
    original <- df
    colors <- c("#124E78", "#A23B34")
    plot <- build_shapture_feature_importance_plot(
      df, 3,
      task_what_0 = condition_labels[1],
      task_what_1 = condition_labels[2],
      group_colors = colors
    )
    values <- matrix(seq_len(24), nrow = 8, dimnames = list(NULL, df$feature))
    beeswarm <- build_shapture_beeswarm_plot(
      df,
      beeswarm_source_S = values / 10 - 1,
      beeswarm_source_X = as.data.frame(values, check.names = FALSE),
      top_feature_show = 3,
      task_what_0 = condition_labels[1],
      task_what_1 = condition_labels[2],
      group_colors = colors,
      continuous_colors = c("blue", "red")
    )

    expect_identical(df, original)
    expect_identical(as.character(plot$data$feature), df$feature)
    expect_identical(plot$data$final_score, df$final_score)
    expect_identical(levels(plot$data$feature), rev(df$feature))
    expect_identical(
      plot$scales$get_scales("y")$labels,
      beeswarm$scales$get_scales("y")$labels
    )
    expect_identical(
      unname(plot$scales$get_scales("y")$labels),
      c(
        "<span style='color:#A23B34'>gene&lt;&amp;&gt;</span>",
        "<span style='color:#124E78'>gene_B</span>",
        "<span style='color:grey55'>gene_C</span>"
      )
    )
    expect_s3_class(plot$theme$axis.text.y, "element_markdown")
    expect_equal(
      plot$scales$get_scales("fill")$breaks,
      beeswarm$scales$get_scales("fill")$breaks
    )
    expect_equal(
      plot$scales$get_scales("fill")$labels,
      beeswarm$scales$get_scales("fill")$labels
    )
    expect_equal(plot$labels$fill, "Model\nexpression\npattern")
    expect_false(any(grepl("Insufficient", plot$scales$get_scales("fill")$labels)))

    baseline <- build_shapture_feature_importance_plot(df, 3)
    baseline_built <- ggplot2::ggplot_build(baseline)
    colored_built <- ggplot2::ggplot_build(plot)
    expect_equal(colored_built$data[1:2], baseline_built$data)
    expect_equal(
      colored_built$layout$panel_params[[1]]$x.range,
      baseline_built$layout$panel_params[[1]]$x.range
    )
    expect_equal(
      colored_built$layout$panel_params[[1]]$y.range,
      baseline_built$layout$panel_params[[1]]$y.range
    )
    expect_warning(ggplot2::ggplotGrob(plot), NA)
  }
})


test_that("Feature-importance top-N selection and old two-argument calls still work", {
  df <- make_feature_label_test_data()
  plot <- build_shapture_feature_importance_plot(df, 1, "control", "case")
  expect_equal(as.character(plot$data$feature), df$feature[1])
  expect_equal(plot$data$final_score, df$final_score[1])
  expect_match(unname(plot$scales$get_scales("y")$labels), "darkred", fixed = TRUE)
  expect_length(plot$scales$get_scales("fill")$breaks, 2)
  expect_warning(ggplot2::ggplotGrob(plot), NA)

  all_features <- build_shapture_feature_importance_plot(df, 100, "control", "case")
  expect_equal(nrow(all_features$data), nrow(df))
  expect_warning(ggplot2::ggplotGrob(all_features), NA)

  legacy <- build_shapture_feature_importance_plot(df[, c("feature", "final_score")], 3)
  expect_null(legacy$scales$get_scales("fill"))
  expect_equal(legacy$theme$axis.text.y$colour, "black")
  expect_warning(ggplot2::ggplotGrob(legacy), NA)
})


test_that("Feature-importance color inputs are validated", {
  df <- make_feature_label_test_data()
  expect_error(
    build_shapture_feature_importance_plot(df, 3, task_what_0 = "control"),
    "two distinct condition labels"
  )
  expect_error(
    build_shapture_feature_importance_plot(df, 3, "case", "case"),
    "two distinct condition labels"
  )
  expect_error(
    build_shapture_feature_importance_plot(df[, c("feature", "final_score")], 3, "control", "case"),
    "must contain Model_expression_pattern"
  )
  expect_error(
    build_shapture_feature_importance_plot(df, 3, "control", "case", "red"),
    "exactly 2 colors"
  )
  expect_error(
    build_shapture_feature_importance_plot(df, 3, "control", "case", c("bad_color", "red")),
    "invalid R color"
  )
})
