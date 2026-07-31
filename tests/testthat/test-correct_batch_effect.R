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
    expect_false("diff.na" %in% names(adjusted))
})


test_that("correct_with_ComBat_df", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    short_df <- example_proteome[example_proteome[["peptide_group_label"]] %in% c("10062_NVGVSFYADKPEVTQEQK_3", "10063_NVGVSFYADKPEVTQEQKK_3"), ]
    combat_df <- correct_with_ComBat(
        short_df,
        example_sample_annotation,
        format = "long",
        fill_the_missing = "error"
    )

    # Example: using ComBat directly
    feature_ids <- sort(unique(short_df$peptide_group_label))
    sample_ids <- sort(unique(short_df$FullRunName))
    example_matrix <- matrix(
        NA_real_,
        nrow = length(feature_ids),
        ncol = length(sample_ids),
        dimnames = list(feature_ids, sample_ids)
    )
    matrix_indices <- cbind(
        match(short_df$peptide_group_label, feature_ids),
        match(short_df$FullRunName, sample_ids)
    )
    expect_false(anyDuplicated(data.frame(matrix_indices)) > 0L)
    example_matrix[matrix_indices] <- short_df$Intensity
    annotation_indices <- match(
        colnames(example_matrix),
        example_sample_annotation$FullRunName
    )
    expect_false(anyNA(annotation_indices))
    combat_example <- ComBat(
        dat = example_matrix,
        batch = as.factor(
            example_sample_annotation$MS_batch[annotation_indices]
        )
    )

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

test_that("ComBat / removeBatchEffect tolerate unused batch factor levels", {
    engine_calls <- list()
    testthat::local_mocked_bindings(
        ComBat = function(dat, batch, mod, ...) {
            engine_calls$combat <<- list(batch = batch, design = mod)
            dat
        },
        removeBatchEffect = function(x, batch, design, ...) {
            engine_calls$remove_batch_effect <<- list(
                batch = batch,
                design = design
            )
            x
        },
        .package = "proBatch"
    )

    data_matrix <- matrix(
        seq_len(32),
        nrow = 4,
        dimnames = list(
            paste0("f", seq_len(4)),
            paste0("s", seq_len(8))
        )
    )
    annotation <- data.frame(
        FullRunName = colnames(data_matrix),
        MS_batch = factor(
            rep(c("b1", "b2"), each = 4),
            levels = c("ghost_batch", "b1", "b2")
        ),
        condition = factor(
            rep(c("A", "B"), times = 4),
            levels = c("ghost_condition", "A", "B")
        ),
        row.names = colnames(data_matrix)
    )

    combat_out <- correct_with_ComBat(
        data_matrix,
        annotation,
        format = "wide",
        covariates_cols = "condition",
        fill_the_missing = "error"
    )
    remove_batch_effect_out <- correct_with_removeBatchEffect(
        data_matrix,
        annotation,
        format = "wide",
        batch_col = "MS_batch",
        covariates_cols = "condition",
        fill_the_missing = "error"
    )

    expect_identical(dim(combat_out), dim(data_matrix))
    expect_identical(dim(remove_batch_effect_out), dim(data_matrix))
    for (engine_call in engine_calls) {
        expect_identical(levels(engine_call$batch), c("b1", "b2"))
        expect_false(any(grepl(
            "ghost",
            colnames(engine_call$design),
            fixed = TRUE
        )))
        expect_identical(qr(engine_call$design)$rank, ncol(engine_call$design))
    }
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
        fill_the_missing = "error",
        mean.only = TRUE,
        ref.batch = "b1"
    )
    direct_long <- correct_with_ComBat(
        long,
        annotation,
        feature_id_col = "Feature",
        format = "long",
        fill_the_missing = "error",
        mean.only = TRUE,
        ref.batch = "b1"
    )
    compatibility <- suppressWarnings(correct_with_ComBat_dm(
        data_matrix,
        annotation,
        fill_the_missing = "error",
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

test_that("canonical keep bypasses missing-value preprocessing", {
    data_matrix <- matrix(
        c(1, NA_real_, 3, 4),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    annotation <- data.frame(
        FullRunName = colnames(data_matrix),
        MS_batch = c("b1", "b2"),
        row.names = colnames(data_matrix),
        stringsAsFactors = FALSE
    )

    matrix_result <- expect_silent(proBatch:::.run_matrix_method(
        data_matrix,
        annotation,
        sample_id_col = "FullRunName",
        fill_the_missing = "keep",
        missing_warning = "missing values would be handled",
        method_fun = function(data_matrix, sample_annotation) data_matrix
    ))

    expect_identical(matrix_result, data_matrix)

    df_long <- data.frame(
        Feature = rep(rownames(data_matrix), times = ncol(data_matrix)),
        FullRunName = rep(colnames(data_matrix), each = nrow(data_matrix)),
        Intensity = as.vector(data_matrix),
        stringsAsFactors = FALSE
    )
    handled <- expect_silent(proBatch:::.handle_missing_for_batch_df(
        df_long = df_long,
        sample_annotation = annotation,
        feature_id_col = "Feature",
        sample_id_col = "FullRunName",
        measure_col = "Intensity",
        fill_the_missing = "keep",
        warning_message = "missing values would be handled"
    ))

    expect_identical(handled$df_long, df_long)
    expect_identical(handled$sample_annotation, annotation)
})

test_that("removeBatchEffect maps keep, drop, and fill missing outcomes", {
    engine_calls <- list()
    testthat::local_mocked_bindings(
        removeBatchEffect = function(x, batch, design, ...) {
            engine_calls[[length(engine_calls) + 1L]] <<- list(
                matrix = x,
                dots = list(...)
            )
            x
        },
        .package = "proBatch"
    )

    data_matrix <- matrix(
        c(1, 2, NA_real_, 4, 5, 6, 7, 8),
        nrow = 2,
        dimnames = list(
            c("f1", "f2"),
            c("s1", "s2", "s3", "s4")
        )
    )
    annotation <- data.frame(
        FullRunName = colnames(data_matrix),
        MS_batch = rep(c("b1", "b2"), each = 2),
        row.names = colnames(data_matrix),
        stringsAsFactors = FALSE
    )
    df_long <- data.frame(
        Feature = rep(rownames(data_matrix), times = ncol(data_matrix)),
        FullRunName = rep(colnames(data_matrix), each = nrow(data_matrix)),
        Intensity = as.vector(data_matrix),
        stringsAsFactors = FALSE
    )

    expect_error(
        correct_with_removeBatchEffect(
            data_matrix,
            annotation,
            format = "wide"
        ),
        "cannot operate with missing values"
    )
    expect_error(
        correct_with_removeBatchEffect(
            df_long,
            annotation,
            feature_id_col = "Feature",
            format = "long"
        ),
        "cannot operate with missing values"
    )

    complete_matrix <- data_matrix
    complete_matrix[is.na(complete_matrix)] <- 0
    expect_error(
        correct_with_removeBatchEffect(
            complete_matrix,
            annotation,
            format = "wide",
            fill_the_missing = NULL
        ),
        "NULL.*ambiguous"
    )

    kept_wide <- expect_silent(correct_with_removeBatchEffect(
        data_matrix,
        annotation,
        format = "wide",
        fill_the_missing = "keep"
    ))
    kept_long <- expect_silent(correct_with_removeBatchEffect(
        df_long,
        annotation,
        feature_id_col = "Feature",
        format = "long",
        fill_the_missing = "keep"
    ))

    expect_identical(kept_wide, data_matrix)
    expect_equal(sum(is.na(kept_long$Intensity)), 1L)

    dropped <- pb_test_expect_warnings(
        correct_with_removeBatchEffect(
            df_long,
            annotation,
            feature_id_col = "Feature",
            format = "long",
            fill_the_missing = "drop_features"
        ),
        "removed 1 rows and 0 columns",
        fixed = TRUE
    )

    expect_identical(unique(dropped$Feature), "f2")

    filled <- expect_silent(
        correct_with_removeBatchEffect(
            df_long,
            annotation,
            feature_id_col = "Feature",
            format = "long",
            fill_the_missing = "fill",
            fill_value = 0
        )
    )

    filled_value <- filled$Intensity[
        filled$Feature == "f1" & filled$FullRunName == "s2"
    ]
    expect_identical(filled_value, 0)
    expect_false(anyNA(filled$Intensity))
    expect_true(is.na(filled$preBatchCorr_Intensity[
        filled$Feature == "f1" & filled$FullRunName == "s2"
    ]))

    legacy_keep <- pb_test_expect_warnings(
        correct_with_removeBatchEffect(
            data_matrix,
            annotation,
            format = "wide",
            fill_the_missing = FALSE
        ),
        "deprecated.*keep"
    )
    expect_identical(legacy_keep, data_matrix)

    expect_length(engine_calls, 5L)
    for (call in engine_calls) {
        expect_false("fill_value" %in% names(call$dots))
    }
})

test_that("unified centering fills working values without rewriting provenance", {
    df_long <- data.frame(
        Feature = rep(c("f1", "f2"), times = 2),
        FullRunName = rep(c("s1", "s2"), each = 2),
        Intensity = c(1, 2, NA_real_, 4),
        stringsAsFactors = FALSE
    )
    annotation <- data.frame(
        FullRunName = c("s1", "s2"),
        MS_batch = c("b1", "b2"),
        stringsAsFactors = FALSE
    )

    for (method in c("MedianCentering", "MeanCentering")) {
        corrected <- expect_silent(correct_batch_effects(
            df_long,
            annotation,
            format = "long",
            discrete_func = method,
            feature_id_col = "Feature",
            fill_the_missing = "fill",
            fill_value = 0,
            no_fit_imputed = FALSE
        ))

        missing_row <- corrected$Feature == "f1" &
            corrected$FullRunName == "s2"
        expect_false(is.na(corrected$Intensity[missing_row]))
        expect_true(is.na(corrected$preBatchCorr_Intensity[missing_row]))

        corrected_keys <- paste(corrected$Feature, corrected$FullRunName)
        input_keys <- paste(df_long$Feature, df_long$FullRunName)
        expected_original <- df_long$Intensity[
            match(corrected_keys, input_keys)
        ]
        expect_equal(
            corrected$preBatchCorr_Intensity,
            expected_original
        )
    }
})

test_that("continuous correction provenance is policy independent", {
    testthat::local_mocked_bindings(
        adjust_batch_trend_df = function(df_long, measure_col, ...) {
            df_long[[paste0("preTrendFit_", measure_col)]] <-
                df_long[[measure_col]]
            df_long$fit <- 10
            df_long[[measure_col]] <- df_long[[measure_col]] + df_long$fit
            df_long
        },
        .package = "proBatch"
    )

    df_long <- data.frame(
        Feature = rep(c("f1", "f2"), times = 4),
        FullRunName = rep(paste0("s", 1:4), each = 2),
        Intensity = seq_len(8),
        stringsAsFactors = FALSE
    )
    annotation <- data.frame(
        FullRunName = paste0("s", 1:4),
        MS_batch = rep(c("b1", "b2"), each = 2),
        stringsAsFactors = FALSE
    )

    corrected <- lapply(c("error", "keep"), function(policy) {
        expect_silent(correct_batch_effects(
            df_long,
            annotation,
            format = "long",
            continuous_func = "mock_trend",
            discrete_func = "MedianCentering",
            feature_id_col = "Feature",
            fill_the_missing = policy,
            no_fit_imputed = FALSE
        ))
    })

    input_keys <- paste(df_long$Feature, df_long$FullRunName)
    for (result in corrected) {
        result_keys <- paste(result$Feature, result$FullRunName)
        expect_equal(
            result$preBatchCorr_Intensity,
            df_long$Intensity[match(result_keys, input_keys)]
        )
    }
    expect_equal(
        corrected[[1L]]$preBatchCorr_Intensity,
        corrected[[2L]]$preBatchCorr_Intensity
    )
})

test_that("ComBat applies canonical fill policy before engine dispatch", {
    engine_calls <- list()
    testthat::local_mocked_bindings(
        run_ComBat_core = function(...) {
            arguments <- list(...)
            engine_calls[[length(engine_calls) + 1L]] <<- arguments
            arguments$data_matrix
        },
        .package = "proBatch"
    )

    data_matrix <- matrix(
        c(1, 2, NA_real_, 4, 5, 6, 7, 8),
        nrow = 2,
        dimnames = list(
            c("f1", "f2"),
            c("s1", "s2", "s3", "s4")
        )
    )
    annotation <- data.frame(
        FullRunName = colnames(data_matrix),
        MS_batch = rep(c("b1", "b2"), each = 2),
        row.names = colnames(data_matrix),
        stringsAsFactors = FALSE
    )

    expect_error(
        correct_with_ComBat(
            data_matrix,
            annotation,
            format = "wide"
        ),
        "ComBat cannot operate with missing values"
    )

    filled <- expect_silent(correct_with_ComBat(
        data_matrix,
        annotation,
        format = "wide",
        fill_the_missing = "fill",
        fill_value = 0
    ))

    expect_false(anyNA(filled))
    expect_identical(unname(filled[1, 2]), 0)
    expect_length(engine_calls, 1L)
    expect_false("fill_value" %in% names(engine_calls[[1L]]))
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
        fill_the_missing = "keep",
        format = "long"
    )

    expect_true("fit" %in% names(corrected))
    expect_equal(nrow(corrected), nrow(short_df))
})

test_that("correct_batch_effects_dm returns matrix", {
    pb_test_load_example_data()

    corrected <- correct_batch_effects(example_proteome_matrix, example_sample_annotation,
        discrete_func = "MedianCentering", no_fit_imputed = FALSE,
        fill_the_missing = "keep",
        format = "wide"
    )

    expect_true(is.matrix(corrected))
    expect_equal(dim(corrected), dim(example_proteome_matrix))
})

test_that("deprecated batch-correction wrappers forward compatibility calls", {
    forwarded <- list()
    remove_forwarded <- list()
    testthat::local_mocked_bindings(
        correct_batch_effects = function(...) {
            arguments <- list(...)
            forwarded[[length(forwarded) + 1L]] <<- arguments
            arguments$x
        },
        correct_with_removeBatchEffect = function(...) {
            arguments <- list(...)
            remove_forwarded[[length(remove_forwarded) + 1L]] <<- arguments
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
        Condition = c("A", "B"),
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
            fill_the_missing = "keep"
        ),
        "deprecated"
    )
    expect_warning(
        dm_result <- correct_batch_effects_dm(
            data_matrix,
            annotation,
            fill_the_missing = "keep"
        ),
        "deprecated"
    )
    expect_warning(
        dm_remove_result <- correct_batch_effects_dm(
            data_matrix,
            annotation,
            discrete_func = "removeBatchEffect",
            fill_the_missing = "keep"
        ),
        "deprecated"
    )
    expect_warning(
        remove_result <- correct_with_removeBatchEffect_dm(
            data_matrix,
            annotation,
            covariates_cols = "Condition",
            fill_the_missing = "keep",
            robust = TRUE
        ),
        "deprecated"
    )

    expect_identical(df_result, df_long)
    expect_identical(dm_result, data_matrix)
    expect_identical(dm_remove_result, data_matrix)
    expect_identical(remove_result, data_matrix)
    expect_length(forwarded, 3L)
    expect_length(remove_forwarded, 1L)
    expect_identical(forwarded[[1L]]$format, "long")
    expect_identical(forwarded[[2L]]$format, "wide")
    expect_identical(forwarded[[3L]]$format, "wide")
    expect_identical(forwarded[[1L]]$discrete_func, "MedianCentering")
    expect_identical(forwarded[[2L]]$discrete_func, "MedianCentering")
    expect_identical(forwarded[[3L]]$discrete_func, "removeBatchEffect")
    expect_identical(forwarded[[1L]]$fill_the_missing, "keep")
    expect_identical(forwarded[[2L]]$fill_the_missing, "keep")
    expect_identical(forwarded[[3L]]$fill_the_missing, "keep")
    expect_identical(remove_forwarded[[1L]]$format, "wide")
    expect_identical(
        remove_forwarded[[1L]]$covariates_cols,
        "Condition"
    )
    expect_identical(remove_forwarded[[1L]]$fill_the_missing, "keep")
    expect_true(remove_forwarded[[1L]]$robust)
})

test_that("correction wrappers have one top-level definition each", {
    check_root <- normalizePath(
        file.path(testthat::test_path(), "..", ".."),
        mustWork = TRUE
    )
    candidates <- c(
        check_root,
        file.path(check_root, "00_pkg_src", "proBatch")
    )
    valid <- vapply(
        candidates,
        function(candidate) {
            file.exists(file.path(candidate, "DESCRIPTION")) &&
                file.exists(
                    file.path(candidate, "R", "correct_batch_effects.R")
                )
        },
        logical(1L)
    )
    if (!any(valid)) {
        skip("Definition checks require the proBatch source package")
    }

    source_root <- normalizePath(
        candidates[[which(valid)[[1L]]]],
        mustWork = TRUE
    )
    r_files <- sort(list.files(
        file.path(source_root, "R"),
        pattern = "[.]R$",
        full.names = TRUE
    ))
    definitions <- unlist(lapply(r_files, function(file) {
        expressions <- parse(file = file)
        vapply(expressions, function(expression) {
            is_assignment <- is.call(expression) &&
                identical(expression[[1L]], as.name("<-"))
            is_function <- is_assignment &&
                is.call(expression[[3L]]) &&
                identical(expression[[3L]][[1L]], as.name("function"))
            if (!is_function || !is.symbol(expression[[2L]])) {
                return(NA_character_)
            }
            as.character(expression[[2L]])
        }, character(1L))
    }), use.names = FALSE)
    definitions <- definitions[!is.na(definitions)]

    wrapper_symbols <- c(
        "correct_batch_effects_df",
        "correct_batch_effects_dm",
        "correct_with_removeBatchEffect_dm"
    )
    definition_counts <- vapply(
        wrapper_symbols,
        function(symbol) sum(definitions == symbol),
        integer(1L)
    )

    expect_identical(
        definition_counts,
        stats::setNames(rep(1L, length(wrapper_symbols)), wrapper_symbols)
    )
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
    observed <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        adjust_batch_trend_df = function(df_long, sample_annotation, span, ...) {
            observed$span <- span
            df_long
        },
        .package = "proBatch"
    )
    input <- matrix(
        1:4,
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    annotation <- data.frame(
        FullRunName = c("s1", "s2"),
        stringsAsFactors = FALSE
    )

    result <- adjust_batch_trend_dm(
        input,
        annotation,
        return_fit_df = FALSE,
        span = 0.37
    )

    expect_identical(observed$span, 0.37)
    expect_equal(result, input)
})

test_that("correct_batch_effects reports sample mismatches and orders survivors", {
    data_matrix <- matrix(
        seq_len(6),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s3", "s1", "s2"))
    )
    annotation <- data.frame(
        FullRunName = c("s2", "s3", "s4"),
        MS_batch = c("b1", "b2", "b3"),
        stringsAsFactors = FALSE
    )

    corrected <- pb_test_expect_warnings(
        correct_batch_effects(
            data_matrix,
            annotation,
            format = "wide",
            discrete_func = "MedianCentering",
            fill_the_missing = "keep",
            no_fit_imputed = FALSE
        ),
        paste0(
            "Mismatch between sample_annotation and df_long samples; ",
            "will merge on intersecting IDs only. ",
            "1 sample(s) only in sample_annotation: s4; ",
            "1 sample(s) only in df_long: s1."
        ),
        fixed = TRUE
    )

    expect_identical(colnames(corrected), c("s3", "s2"))
})
