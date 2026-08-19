make_patient_dataset <- function(
    patients_condition0,
    patients_condition1,
    units_per_patient = 25L
) {
  make_condition <- function(patient_ids, condition) {
    do.call(rbind, lapply(patient_ids, function(patient_id) {
      data.frame(
        condition = condition,
        sample_id = patient_id,
        .shapture_unit_id = paste(condition, patient_id, seq_len(units_per_patient), sep = "_"),
        feature_1 = seq_len(units_per_patient),
        stringsAsFactors = FALSE
      )
    }))
  }

  rbind(
    make_condition(patients_condition0, "0"),
    make_condition(patients_condition1, "1")
  )
}

expect_globally_disjoint_plans <- function(result) {
  for (plan in result$plans) {
    expect_length(intersect(plan$train_patient_ids, plan$test_patient_ids), 0L)
  }
}

test_that("shared-patient planning handles singleton numeric candidate pools", {
  dataset <- data.frame(
    condition = c("0", "0", "1", "1"),
    .shapture_unit_id = paste0("unit_", seq_len(4)),
    stringsAsFactors = FALSE
  )

  result <- build_coverage_aware_split_plan(
    dataset = dataset,
    random_seeds = 17,
    train_ratio = 0.5,
    test_ratio = 0.5,
    each_condition_total_sample_no = 2
  )

  plan <- result$plans[[1]]
  expect_setequal(
    c(plan$train_unit_ids, plan$test_unit_ids),
    dataset$.shapture_unit_id
  )
  expect_length(intersect(plan$train_unit_ids, plan$test_unit_ids), 0L)
})

test_that("five patients per condition produce exact 3/2 patient folds", {
  dataset <- make_patient_dataset(
    patients_condition0 = paste0("control_", seq_len(5)),
    patients_condition1 = paste0("case_", seq_len(5))
  )

  result <- build_unseen_patient_split_plan(
    dataset = dataset,
    random_seeds = seq_len(10),
    train_ratio = 0.7,
    test_ratio = 0.3,
    each_condition_total_sample_no = 62,
    min_units_per_patient = 20,
    min_train_units_per_condition = 50
  )

  expect_true(result$success)
  expect_equal(unname(result$n_train_patients), c(3L, 3L))
  expect_equal(unname(result$n_test_patients), c(2L, 2L))
  expect_equal(result$train_rows_per_condition, 50L)
  expect_equal(result$test_rows_per_condition, 40L)
  expect_globally_disjoint_plans(result)

  eligible <- result$patient_coverage_summary$eligible
  expect_equal(result$patient_coverage_summary$test_iterations[eligible], rep(4L, 10))
  expect_equal(result$patient_coverage_summary$train_iterations[eligible], rep(6L, 10))
})

test_that("eight patients use rounded 60/40 patient allocation", {
  dataset <- make_patient_dataset(
    patients_condition0 = paste0("control_", seq_len(8)),
    patients_condition1 = paste0("case_", seq_len(8))
  )

  result <- build_unseen_patient_split_plan(
    dataset = dataset,
    random_seeds = seq_len(8),
    train_ratio = 0.7,
    test_ratio = 0.3,
    each_condition_total_sample_no = 100,
    min_units_per_patient = 20,
    min_train_units_per_condition = 50
  )

  expect_true(result$success)
  expect_equal(unname(result$n_train_patients), c(5L, 5L))
  expect_equal(unname(result$n_test_patients), c(3L, 3L))
  expect_globally_disjoint_plans(result)
})

test_that("paired one-vs-others donors remain globally separated", {
  patients <- paste0("donor_", seq_len(5))
  dataset <- make_patient_dataset(
    patients_condition0 = patients,
    patients_condition1 = patients
  )

  result <- build_unseen_patient_split_plan(
    dataset = dataset,
    random_seeds = seq_len(10),
    train_ratio = 0.7,
    test_ratio = 0.3,
    each_condition_total_sample_no = 62,
    min_units_per_patient = 20,
    min_train_units_per_condition = 50
  )

  expect_true(result$success)
  expect_globally_disjoint_plans(result)
  for (plan in result$plans) {
    expect_equal(length(plan$test_patient_ids), 2L)
    expect_equal(length(plan$train_patient_ids), 3L)
  }
})

test_that("insufficient eligible patients return a fallback reason", {
  dataset <- make_patient_dataset(
    patients_condition0 = paste0("control_", seq_len(5)),
    patients_condition1 = paste0("case_", seq_len(5))
  )
  dataset <- dataset[!(dataset$sample_id == "case_5" &
                         ave(seq_len(nrow(dataset)), dataset$sample_id, FUN = seq_along) > 10), ]

  result <- build_unseen_patient_split_plan(
    dataset = dataset,
    random_seeds = seq_len(3),
    train_ratio = 0.7,
    test_ratio = 0.3,
    each_condition_total_sample_no = 50,
    min_units_per_patient = 20,
    min_train_units_per_condition = 50
  )

  expect_false(result$success)
  expect_match(result$reason, "at least 5 eligible patients per condition")
  expect_match(result$reason, "1=4", fixed = TRUE)
})

test_that("metacell thresholds support eligible patient-held-out folds", {
  dataset <- make_patient_dataset(
    patients_condition0 = paste0("control_", seq_len(5)),
    patients_condition1 = paste0("case_", seq_len(5)),
    units_per_patient = 7
  )

  result <- build_unseen_patient_split_plan(
    dataset = dataset,
    random_seeds = seq_len(5),
    train_ratio = 0.7,
    test_ratio = 0.3,
    each_condition_total_sample_no = 12,
    min_units_per_patient = 5,
    min_train_units_per_condition = 20
  )

  expect_true(result$success)
  expect_equal(result$train_rows_per_condition, 20L)
  expect_equal(result$test_rows_per_condition, 10L)
  expect_globally_disjoint_plans(result)
})

test_that("rows are balanced across held-out patients", {
  make_unequal_condition <- function(condition, prefix) {
    do.call(rbind, lapply(seq_len(5), function(i) {
      n <- if (i == 1L) 100L else 20L
      data.frame(
        condition = condition,
        sample_id = paste0(prefix, i),
        .shapture_unit_id = paste(condition, prefix, i, seq_len(n), sep = "_"),
        stringsAsFactors = FALSE
      )
    }))
  }
  dataset <- rbind(
    make_unequal_condition("0", "control_"),
    make_unequal_condition("1", "case_")
  )

  result <- build_unseen_patient_split_plan(
    dataset = dataset,
    random_seeds = 101,
    train_ratio = 0.7,
    test_ratio = 0.3,
    each_condition_total_sample_no = 70,
    min_units_per_patient = 20,
    min_train_units_per_condition = 50
  )

  expect_true(result$success)
  plan <- result$plans[[1]]
  selected <- dataset[dataset$.shapture_unit_id %in% plan$test_unit_ids, ]
  counts <- table(selected$condition, selected$sample_id)
  nonzero_counts <- counts[counts > 0]
  expect_equal(length(unique(as.integer(nonzero_counts))), 1L)

  training <- dataset[dataset$.shapture_unit_id %in% plan$train_unit_ids, ]
  represented_training_patients <- split(
    training$sample_id,
    training$condition
  )
  expect_equal(
    vapply(represented_training_patients, function(x) length(unique(x)), integer(1)),
    c(`0` = 3L, `1` = 3L)
  )
})

test_that("patient-aware feature plans cover genes and held-out patients", {
  patient_keys <- list(
    c("0::p1", "1::p6"),
    c("0::p2", "1::p7"),
    c("0::p3", "1::p8"),
    c("0::p4", "1::p9"),
    c("0::p5", "1::p10")
  )

  result <- build_patient_coverage_aware_feature_plan(
    feature_names = paste0("gene_", seq_len(10)),
    random_seeds = seq_len(5),
    test_patient_condition_keys = patient_keys,
    eligible_patient_condition_keys = c(
      unlist(patient_keys, use.names = FALSE),
      "0::never_held_out"
    ),
    fraction = 0.20
  )

  expect_equal(vapply(result$plans, length, integer(1)), rep(2L, 5))
  expect_equal(result$summary$unique_features_included, 10L)
  expect_equal(result$summary$min_gene_inclusion_count, 1L)
  expect_equal(result$patient_feature_coverage$eligible_patient_condition_units, rep(11L, 10))
  expect_equal(result$patient_feature_coverage$gene_inclusion_count, rep(1L, 10))
  expect_true(all(result$patient_feature_coverage$patient_condition_units_with_shap >= 2L))
  expect_true(all(result$patient_feature_counts[, "0::never_held_out"] == 0L))
})
