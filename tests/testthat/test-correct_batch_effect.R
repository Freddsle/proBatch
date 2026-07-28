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

test_that("center_feature_batch validates wide inputs by sample identity", {
    data_matrix <- matrix(
        c(1, 3, 2, 4),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    annotation <- data.frame(
        FullRunName = c("s1", "s2", "unused"),
        MS_batch = c("b1", "b2", "b3"),
        stringsAsFactors = FALSE
    )

    expect_warning(
        centered <- center_feature_batch(
            data_matrix,
            annotation,
            format = "wide",
            stat = "means",
            no_fit_imputed = FALSE
        ),
        "will merge on intersecting IDs only",
        fixed = TRUE
    )

    expect_true(is.matrix(centered))
    expect_identical(dim(centered), dim(data_matrix))
    expect_setequal(colnames(centered), colnames(data_matrix))

    expect_error(
        center_feature_batch(
            as.data.frame(data_matrix),
            annotation,
            format = "wide",
            no_fit_imputed = FALSE
        ),
        "format='wide' requires a numeric matrix",
        fixed = TRUE
    )
    expect_error(
        center_feature_batch(
            matrix(
                as.character(data_matrix),
                nrow = nrow(data_matrix),
                dimnames = dimnames(data_matrix)
            ),
            annotation,
            format = "wide",
            no_fit_imputed = FALSE
        ),
        "format='wide' requires a numeric matrix",
        fixed = TRUE
    )
    expect_error(
        center_feature_batch(
            data_matrix,
            sample_annotation = NULL,
            format = "wide",
            no_fit_imputed = FALSE
        ),
        "requires `sample_annotation` and column names",
        fixed = TRUE
    )
    expect_error(
        center_feature_batch(
            unname(data_matrix),
            annotation,
            format = "wide",
            no_fit_imputed = FALSE
        ),
        "requires `sample_annotation` and column names",
        fixed = TRUE
    )
    expect_error(
        center_feature_batch(
            data_matrix,
            annotation[annotation$FullRunName != "s2", , drop = FALSE],
            format = "wide",
            no_fit_imputed = FALSE
        ),
        "Not all matrix column names found in sample_annotation[[sample_id_col]].",
        fixed = TRUE
    )
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
    expect_equal(adjusted$fit[[1]], 1830358, tolerance = 1, ignore_attr = TRUE)
})


test_that("correct_with_ComBat_df", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    short_df <- example_proteome[example_proteome[["peptide_group_label"]] %in% c("10062_NVGVSFYADKPEVTQEQK_3", "10063_NVGVSFYADKPEVTQEQKK_3"), ]
    combat_df <- correct_with_ComBat(short_df, example_sample_annotation, format = "long")

    # Example: using ComBat directly
    example_matrix <- reshape2::dcast(short_df, peptide_group_label ~ FullRunName, value.var = "Intensity")
    rownames(example_matrix) <- example_matrix$peptide_group_label
    example_matrix$peptide_group_label <- NULL
    combat_example <- ComBat(dat = as.matrix(example_matrix), batch = as.factor(example_sample_annotation$MS_batch))

    expect_equal(combat_df[["peptide_group_label"]][1], "10062_NVGVSFYADKPEVTQEQK_3")
    expect_equal(combat_df[1, ][["Intensity"]], combat_example[combat_df[1, ][["peptide_group_label"]], combat_df[1, ][["FullRunName"]]], tolerance = 1)

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

test_that("ComBat entry points forward supported engine arguments", {
    calls <- list()
    testthat::local_mocked_bindings(
        run_ComBat_core = function(...) {
            arguments <- list(...)
            calls[[length(calls) + 1L]] <<- arguments
            arguments$data_matrix
        },
        .package = "proBatch"
    )

    data_matrix <- matrix(
        seq_len(8),
        nrow = 2,
        dimnames = list(
            c("f1", "f2"),
            c("s1", "s2", "s3", "s4")
        )
    )
    annotation <- data.frame(
        FullRunName = colnames(data_matrix),
        MS_batch = c("b1", "b1", "b2", "b2"),
        row.names = colnames(data_matrix),
        stringsAsFactors = FALSE
    )
    long <- data.frame(
        Feature = rep(rownames(data_matrix), times = ncol(data_matrix)),
        FullRunName = rep(colnames(data_matrix), each = nrow(data_matrix)),
        Intensity = as.vector(data_matrix),
        stringsAsFactors = FALSE
    )

    direct <- correct_with_ComBat(
        data_matrix,
        annotation,
        format = "wide",
        mean.only = TRUE,
        ref.batch = "b1"
    )
    direct_long <- correct_with_ComBat(
        long,
        annotation,
        feature_id_col = "Feature",
        format = "long",
        mean.only = TRUE,
        ref.batch = "b1"
    )
    compatibility <- suppressWarnings(correct_with_ComBat_dm(
        data_matrix,
        annotation,
        mean.only = TRUE,
        ref.batch = "b1"
    ))

    expect_equal(direct, data_matrix)
    expect_s3_class(direct_long, "data.frame")
    expect_equal(compatibility, data_matrix)
    expect_length(calls, 3L)
    for (arguments in calls) {
        expect_true(arguments$mean.only)
        expect_identical(arguments$ref.batch, "b1")
    }
})

test_that("batch annotation alignment rejects duplicate identifiers", {
    duplicate <- data.frame(
        FullRunName = c("s1", "s1"),
        MS_batch = c("b1", "b2"),
        stringsAsFactors = FALSE
    )

    expect_error(
        proBatch:::.align_sample_annotation(
            duplicate,
            sample_ids = "s1",
            sample_id_col = "FullRunName"
        ),
        "duplicate identifiers"
    )
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
    pb_test_load_example_data()

    corrected <- correct_batch_effects(example_proteome_matrix, example_sample_annotation,
        discrete_func = "MedianCentering", no_fit_imputed = FALSE,
        format = "wide"
    )

    expect_true(is.matrix(corrected))
    expect_equal(dim(corrected), dim(example_proteome_matrix))
})

test_that("deprecated batch-correction wrappers select their default method", {
    forwarded <- list()
    testthat::local_mocked_bindings(
        correct_batch_effects = function(...) {
            arguments <- list(...)
            forwarded[[length(forwarded) + 1L]] <<- arguments
            arguments$x
        },
        .package = "proBatch"
    )

    data_matrix <- matrix(
        1:4,
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    annotation <- data.frame(
        FullRunName = colnames(data_matrix),
        MS_batch = c("b1", "b2"),
        stringsAsFactors = FALSE
    )
    df_long <- data.frame(
        peptide_group_label = rep(rownames(data_matrix), times = ncol(data_matrix)),
        FullRunName = rep(colnames(data_matrix), each = nrow(data_matrix)),
        Intensity = as.vector(data_matrix),
        stringsAsFactors = FALSE
    )

    expect_warning(
        df_result <- correct_batch_effects_df(
            df_long,
            annotation,
            fill_the_missing = FALSE
        ),
        "deprecated"
    )
    expect_warning(
        dm_result <- correct_batch_effects_dm(
            data_matrix,
            annotation,
            fill_the_missing = FALSE
        ),
        "deprecated"
    )
    expect_warning(
        dm_remove_result <- correct_batch_effects_dm(
            data_matrix,
            annotation,
            discrete_func = "removeBatchEffect",
            fill_the_missing = FALSE
        ),
        "deprecated"
    )

    expect_identical(df_result, df_long)
    expect_identical(dm_result, data_matrix)
    expect_identical(dm_remove_result, data_matrix)
    expect_length(forwarded, 3L)
    expect_identical(forwarded[[1L]]$format, "long")
    expect_identical(forwarded[[2L]]$format, "wide")
    expect_identical(forwarded[[3L]]$format, "wide")
    expect_identical(forwarded[[1L]]$discrete_func, "MedianCentering")
    expect_identical(forwarded[[2L]]$discrete_func, "MedianCentering")
    expect_identical(forwarded[[3L]]$discrete_func, "removeBatchEffect")
    expect_false(forwarded[[1L]]$fill_the_missing)
    expect_false(forwarded[[2L]]$fill_the_missing)
    expect_false(forwarded[[3L]]$fill_the_missing)
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
    pb_test_load_example_data()

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
