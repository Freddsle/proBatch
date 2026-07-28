test_that("check_sample_consistency", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")
    # TODO: check if the warnings are raised expect_warning

    df_test <- check_sample_consistency(
        sample_annotation = example_sample_annotation,
        df_long = example_proteome,
        sample_id_col = "FullRunName",
        batch_col = NULL, order_col = NULL, facet_col = NULL
    )

    expect_equal(ncol(df_test), ncol(example_sample_annotation) + ncol(example_proteome) - 1)
    expect_equal(nrow(df_test), nrow(example_proteome))

    options(warn = 0)
    expect_warning(check_sample_consistency(
        sample_annotation = NULL, df_long = example_proteome,
        sample_id_col = "FullRunName",
        batch_col = NULL, order_col = NULL, facet_col = NULL
    ))

    falsified_annotation <- example_sample_annotation
    colnames(falsified_annotation)[colnames(falsified_annotation) == "FullRunName"] <- "FalseRunName"
    expect_error(check_sample_consistency(
        sample_annotation = falsified_annotation,
        df_long = example_proteome,
        sample_id_col = "FullRunName",
        batch_col = NULL, order_col = NULL, facet_col = NULL
    ))

    duplicate_annotation <- example_sample_annotation
    duplicate_annotation$FullRunName[2] <-
        duplicate_annotation$FullRunName[1]
    expect_error(
        check_sample_consistency(
            sample_annotation = duplicate_annotation,
            df_long = example_proteome,
            sample_id_col = "FullRunName",
            batch_col = NULL,
            order_col = NULL,
            facet_col = NULL
        ),
        "duplicate identifiers"
    )
})

test_that("define_sample_order", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")
    # TODO: check if the returned class is integer expected_type
    order_col <- "order"
    sample_order <- define_sample_order(
        order_col = order_col,
        sample_annotation = example_sample_annotation,
        facet_col = NULL, batch_col = "MS_batch",
        df_long = example_proteome,
        sample_id_col = "FullRunName",
        color_by_batch = TRUE
    )
    new_order_col <- sample_order$order_col
    df_long <- sample_order$df_long

    expect_type(sample_order, "list")
    expect_equal(df_long$order_col, example_sample_annotation$order_col)
    expect_type(df_long[[new_order_col]], "double")

    order_col <- NULL
    pb_test_expect_warnings(
        sample_order <- define_sample_order(
            order_col = order_col,
            sample_annotation = example_sample_annotation,
            facet_col = NULL, batch_col = "MS_batch",
            df_long = example_proteome,
            sample_id_col = "FullRunName",
            color_by_batch = TRUE
        ),
        c("Order column is NULL", "ordering the samples by batch"),
        fixed = TRUE
    )
    new_order_col <- sample_order$order_col
    df_long <- sample_order$df_long

    expect_type(sample_order, "list")
    expect_equal(df_long$order_col, example_sample_annotation$order_col)
    expect_s3_class(df_long[[new_order_col]], "factor")
})

test_that("adjust_units converts mm to inches", {
    res <- adjust_units("mm", width = 25.4, height = 50.8)
    expect_equal(res$unit, "in")
    expect_equal(res$width, 1)
    expect_equal(res$height, 2)
})

test_that("adjust_units converts cm to inches", {
    res <- adjust_units("cm", width = 2.54, height = 5.08)
    expect_equal(res$unit, "in")
    expect_equal(res$width, 1)
    expect_equal(res$height, 2)
})

test_that("adjust_units leaves inch units unchanged", {
    res <- adjust_units("in", width = 1, height = 2)
    expect_equal(res$unit, "in")
    expect_equal(res$width, 1)
    expect_equal(res$height, 2)
})

test_that("long-matrix preparation forwards canonical missing policy", {
    long <- data.frame(
        Feature = rep(c("f1", "f2"), times = 2),
        FullRunName = rep(c("s1", "s2"), each = 2),
        Intensity = c(1, 2, NA_real_, 4),
        stringsAsFactors = FALSE
    )
    annotation <- data.frame(
        FullRunName = c("s1", "s2"),
        row.names = c("s1", "s2"),
        stringsAsFactors = FALSE
    )

    filled <- expect_silent(proBatch:::.pb_prepare_long_matrix(
        df_long = long,
        sample_annotation = annotation,
        sample_id_col = "FullRunName",
        feature_id_col = "Feature",
        measure_col = "Intensity",
        fill_the_missing = "fill",
        fill_value = 0,
        warning_message = "missing values"
    ))

    expect_false(anyNA(filled$data_matrix))
    expect_true(is.na(
        filled$df_long$Intensity[
            filled$df_long$Feature == "f1" &
                filled$df_long$FullRunName == "s2"
        ]
    ))
    expect_identical(unname(filled$data_matrix["f1", "s2"]), 0)

    qualified <- long
    qualified$Intensity[3] <- 3
    qualified$quality <- c("ok", "ok", "masked", "ok")
    filled_mask <- suppressMessages(proBatch:::.pb_prepare_long_matrix(
        df_long = qualified,
        sample_annotation = annotation,
        sample_id_col = "FullRunName",
        feature_id_col = "Feature",
        measure_col = "Intensity",
        fill_the_missing = "fill",
        fill_value = 0,
        warning_message = "missing values",
        qual_col = "quality",
        qual_value = "masked"
    ))
    expect_false(anyNA(filled_mask$data_matrix))
    expect_identical(unname(filled_mask$data_matrix["f1", "s2"]), 0)
    expect_identical(
        filled_mask$df_long$Intensity[
            filled_mask$df_long$Feature == "f1" &
                filled_mask$df_long$FullRunName == "s2"
        ],
        3
    )

    kept <- suppressMessages(proBatch:::.pb_prepare_long_matrix(
        df_long = qualified,
        sample_annotation = annotation,
        sample_id_col = "FullRunName",
        feature_id_col = "Feature",
        measure_col = "Intensity",
        fill_the_missing = "keep",
        warning_message = "missing values",
        qual_col = "quality",
        qual_value = "masked"
    ))

    expect_identical(kept$df_long$Intensity, qualified$Intensity)
    expect_true(anyNA(kept$data_matrix))
})
