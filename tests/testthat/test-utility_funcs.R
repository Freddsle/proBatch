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
    expect_warning(
        expect_warning(
            sample_order <- define_sample_order(
                order_col = order_col,
                sample_annotation = example_sample_annotation,
                facet_col = NULL, batch_col = "MS_batch",
                df_long = example_proteome,
                sample_id_col = "FullRunName",
                color_by_batch = TRUE
            ),
            "Order column is NULL"
        ),
        "ordering the samples by batch"
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
