test_that("center_feature_batch_medians", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    rows <- which(example_proteome$peptide_group_label == "10062_NVGVSFYADKPEVTQEQK_3")

    proteome <- example_proteome[rows, ]
    median_proteome <- center_feature_batch(
        proteome, example_sample_annotation,
        no_fit_imputed = FALSE,
        stat = "median", format = "long"
    )

    n_batch <- length(unique(median_proteome$MS_batch))
    expect_equal(length(unique(median_proteome$diff_medians)), n_batch)
    expect_equal(length(unique(median_proteome$median_batch)), n_batch)
})


test_that("adjust_batch_trend", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    short_df <- example_proteome[example_proteome[["peptide_group_label"]]
    %in% c("10062_NVGVSFYADKPEVTQEQK_3", "101233_QGFNVVVESGAGEASK_2"), ]

    adjusted <- adjust_batch_trend_df(short_df, example_sample_annotation,
        span = 0.7,
        abs_threshold = 5, pct_threshold = 0.20, keep_all = "all",
        no_fit_imputed = FALSE
    )

    n_batch <- length(unique(example_sample_annotation$MS_batch))

    expect_equal(adjusted[["peptide_group_label"]][1], "10062_NVGVSFYADKPEVTQEQK_3")
    expect_equal(length(unique(adjusted$MS_batch)), n_batch)
    expect_equal(adjusted$fit[1], 1830358, tolerance = 1, ignore_attr = TRUE)
})


test_that("correct_with_ComBat_df", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    short_df <- example_proteome[example_proteome[["peptide_group_label"]]
    %in% c("10062_NVGVSFYADKPEVTQEQK_3", "10063_NVGVSFYADKPEVTQEQKK_3"), ]
    combat_df <- correct_with_ComBat(short_df, example_sample_annotation, format = "long")

    expect_equal(combat_df[["peptide_group_label"]][1], "10062_NVGVSFYADKPEVTQEQK_3")
    expect_equal(combat_df[["Intensity"]][1], 768661.4, tolerance = 1)

    batch_1 <- example_sample_annotation$FullRunName[example_sample_annotation$MS_batch == "Batch_1"]
    batch_2 <- example_sample_annotation$FullRunName[example_sample_annotation$MS_batch == "Batch_2"]

    matrix_batch_1 <- short_df[short_df$FullRunName %in% batch_1, ]
    matrix_batch_2 <- short_df[short_df$FullRunName %in% batch_2, ]

    combat_batch_1 <- combat_df[combat_df$FullRunName %in% batch_1, ]
    combat_batch_2 <- combat_df[combat_df$FullRunName %in% batch_2, ]

    t_test_matrix <- t.test(matrix_batch_1$Intensity, matrix_batch_2$Intensity)
    t_test_combat <- t.test(combat_batch_1$Intensity, combat_batch_2$Intensity)

    expect_lt(t_test_matrix$p.value, 0.05)
    expect_gt(t_test_combat$p.value, 0.05)
})

# test_that("center_feature_batch_means_df", {
#     data(example_proteome, package = "proBatch")
#     data(example_sample_annotation, package = "proBatch")

#     rows <- which(example_proteome$peptide_group_label == "10062_NVGVSFYADKPEVTQEQK_3")
#     proteome <- example_proteome[rows, ]
#     means_df <- center_feature_batch_means_df(proteome, example_sample_annotation, no_fit_imputed = FALSE)

#     n_batch <- length(unique(means_df$MS_batch))
#     expect_equal(length(unique(means_df$diff)), n_batch)
#     expect_equal(length(unique(means_df$mean_batch)), n_batch)
# })

test_that("correct_batch_effects_df wrapper", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    short_df <- example_proteome[example_proteome[["peptide_group_label"]] %in%
        c("10062_NVGVSFYADKPEVTQEQK_3", "101233_QGFNVVVESGAGEASK_2"), ]

    corrected <- correct_batch_effects(short_df, example_sample_annotation,
        continuous_func = "loess_regression",
        discrete_func = "MedianCentering",
        span = 0.7,
        min_measurements = 8,
        no_fit_imputed = FALSE,
        format = "long"
    )

    expect_true("fit" %in% names(corrected))
    expect_equal(nrow(corrected), nrow(short_df))
})

test_that("correct_batch_effects_dm returns matrix", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    corrected <- correct_batch_effects(example_proteome_matrix, example_sample_annotation,
        discrete_func = "MedianCentering", no_fit_imputed = FALSE,
        format = "wide"
    )

    expect_true(is.matrix(corrected))
    expect_equal(dim(corrected), dim(example_proteome_matrix))
})


test_that("adjust_batch_trend_df keeps order column", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    short_df <- example_proteome[example_proteome[["peptide_group_label"]] %in%
        c("10062_NVGVSFYADKPEVTQEQK_3", "101233_QGFNVVVESGAGEASK_2"), ]

    adjusted <- adjust_batch_trend_df(short_df, example_sample_annotation,
        order_col = "order", keep_all = "all", fit_func = "loess_regression",
        min_measurements = 8, no_fit_imputed = FALSE
    )

    expect_true("order" %in% names(adjusted))
    expect_true("fit" %in% names(adjusted))
    expect_equal(nrow(adjusted), nrow(short_df))
})

test_that("adjust_batch_trend_dm forwards arguments", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    feature_subset <- rownames(example_proteome_matrix) %in%
        c("10062_NVGVSFYADKPEVTQEQK_3", "101233_QGFNVVVESGAGEASK_2")
    sub_matrix <- example_proteome_matrix[feature_subset, , drop = FALSE]

    res <- adjust_batch_trend_dm(sub_matrix, example_sample_annotation,
        order_col = "order", fit_func = "loess_regression",
        min_measurements = 8, no_fit_imputed = FALSE
    )

    expect_true(is.list(res))
    expect_true(is.matrix(res$corrected_dm))
    expect_equal(nrow(res$corrected_dm), sum(feature_subset))
})
