test_that("fit works", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    test_annotation <- example_sample_annotation[example_sample_annotation$MS_batch == "Batch_1", ]
    selected_files <- test_annotation$FullRunName

    df_selected <- example_proteome[example_proteome$peptide_group_label == example_proteome$peptide_group_label[1], ]
    df_selected <- df_selected[df_selected$FullRunName %in% selected_files, ]
    df_selected <- merge(df_selected, test_annotation, by = "FullRunName")

    expect_warning(
        fit_values <- fit_nonlinear(df_selected),
        "Imputed-value column present; fitting only to measured (non-imputed) values.",
        fixed = TRUE
    )
    expect_length(fit_values, nrow(df_selected))
})

test_that("rle_func finds longest measured run", {
    df <- data.frame(order = 1:6, Intensity = c(1, 2, 3, NA, 4, 5))
    expect_equal(rle_func(df), 3)
})

test_that("rle_func handles custom columns and all NA", {
    df <- data.frame(idx = 6:1, value = c(NA, NA, 3, 4, 5, NA))
    expect_equal(rle_func(df, measure_col = "value", order_col = "idx"), 3)
    df_na <- data.frame(idx = 1:3, value = c(NA, NA, NA))
    expect_equal(rle_func(df_na, measure_col = "value", order_col = "idx"), 0)
})

test_that("fit_nonlinear returns fitted curve with sufficient data", {
    set.seed(1)
    df <- data.frame(
        order = 1:10,
        Intensity = sin(1:10),
        m_score = rep(0, 10)
    )
    vals <- fit_nonlinear(df, min_measurements = 5, qual_col = NULL, no_fit_imputed = FALSE)
    expect_length(vals, nrow(df))
    expect_true(any(!is.na(vals)))
})

test_that("fit_nonlinear returns NA vector when insufficient data", {
    df <- data.frame(
        order = 1:4,
        Intensity = c(1, NA, NA, 2),
        m_score = c(0, 2, 2, 0)
    )
    vals <- suppressWarnings(
        fit_nonlinear(df, min_measurements = 5, qual_col = NULL, no_fit_imputed = FALSE)
    )
    expect_true(all(is.na(vals)))
})

test_that("fit_nonlinear excludes imputed values when requested", {
    df <- data.frame(
        order = 1:10,
        Intensity = c(1, 1, 1, 100, 1, 1, 1, 1, 1, 1),
        m_score = c(0, 0, 0, 2, 0, 0, 0, 0, 0, 0)
    )

    vals_without_imputed <- suppressWarnings(
        fit_nonlinear(
            df_feature_batch = df,
            min_measurements = 3,
            no_fit_imputed = TRUE,
            qual_col = "m_score",
            qual_value = 2
        )
    )

    df_masked <- df
    df_masked$Intensity[df_masked$m_score == 2] <- NA_real_
    vals_expected <- suppressWarnings(
        fit_nonlinear(
            df_feature_batch = df_masked,
            min_measurements = 3,
            no_fit_imputed = FALSE,
            qual_col = NULL
        )
    )

    expect_equal(vals_without_imputed, vals_expected)
})

test_that("LOESS helpers suppress predictions outside fitted support", {
    x_to_fit <- 1:12
    y_to_fit <- sin(x_to_fit / 3)
    x_all <- c(0, x_to_fit, 13)
    y_all <- rep(NA_real_, length(x_all))

    direct <- loess_regression(
        x_to_fit,
        y_to_fit,
        x_all,
        y_all,
        span = 1
    )
    optimized <- loess_regression_opt(
        x_to_fit,
        y_to_fit,
        x_all,
        y_all,
        bws = c(0.5, 1, 1.5)
    )

    for (fit in list(direct, optimized)) {
        expect_true(all(is.na(fit[c(1, length(fit))])))
        expect_true(all(is.finite(fit[-c(1, length(fit))])))
    }
})

test_that("LOESS warning and error paths return NA fallbacks", {
    warning_fit <- NULL
    warning_messages <- capture.output(
        expect_warning(
            warning_fit <- loess_regression(
                rep(1, 6),
                1:6,
                1:4,
                rep(NA_real_, 4),
                feature_id = "feature",
                batch_id = "batch"
            ),
            NA
        ),
        type = "message"
    )
    expect_true(all(is.na(warning_fit)))

    direct_error <- NULL
    direct_messages <- capture.output(
        direct_error <- loess_regression(
            1:3,
            1:2,
            1:4,
            rep(NA_real_, 4),
            feature_id = "feature",
            batch_id = "batch"
        ),
        type = "message"
    )
    optimized_error <- NULL
    optimized_messages <- capture.output(
        optimized_error <- loess_regression_opt(
            1:3,
            1:2,
            1:4,
            rep(NA_real_, 4),
            feature_id = "feature",
            batch_id = "batch",
            bws = 1
        ),
        type = "message"
    )

    expect_match(
        paste(warning_messages, collapse = "\n"),
        "could not be fit with LOESS"
    )
    expect_match(
        paste(direct_messages, collapse = "\n"),
        "could not be fit with LOESS"
    )
    expect_match(
        paste(optimized_messages, collapse = "\n"),
        "could not be fit with optimised LOESS"
    )
    expect_true(all(is.na(direct_error)))
    expect_true(all(is.na(optimized_error)))
})

test_that("fit_nonlinear does not mutate caller-owned tabular inputs", {
    inputs <- list(
        data.frame(
            order = seq_len(4),
            Intensity = c(1, 100, 1, 1),
            m_score = c(0, 2, 0, 0)
        ),
        data.table::data.table(
            order = seq_len(4),
            Intensity = c(1, 100, 1, 1),
            m_score = c(0, 2, 0, 0)
        )
    )

    for (input in inputs) {
        before <- data.table::copy(input)
        suppressWarnings(fit_nonlinear(input, min_measurements = 5))
        expect_identical(input, before)
    }
})
