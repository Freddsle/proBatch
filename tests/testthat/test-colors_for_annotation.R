test_that("map_factors_to_colors", {
    data(example_sample_annotation, package = "proBatch")

    factor_columns <- c("MS_batch", "Strain")
    sample_annt <- example_sample_annotation[, factor_columns]
    sample_annt[factor_columns] <- lapply(sample_annt[factor_columns], factor)

    expect_warning(
        sample_color_annt <- map_factors_to_colors(sample_annt),
        "Too many colors, consider merging some covariate values for better visualisation"
    )

    expect_equal(names(sample_color_annt$MS_batch[1]), "Batch_1")
    expect_equal(sample_color_annt$MS_batch[[1]], "#7FC97F")
    expect_equal(sample_color_annt$MS_batch[[2]], "#BEAED4")
})


test_that("map_numbers_to_colors", {
    data(example_sample_annotation, package = "proBatch")

    numeric_columns <- c("DateTime", "order")
    sample_annt <- example_sample_annotation[, numeric_columns]
    sample_color_annt <- map_numbers_to_colors(sample_annt)

    expect_equal(names(sample_color_annt), c("DateTime", "order"))
    expect_equal(sample_color_annt$DateTime[[1]], "#A6611A")
    expect_equal(sample_color_annt$order[[1]], "#D01C8B")
})


test_that("sample_annotation_to_colors", {
    data(example_sample_annotation, package = "proBatch")

    factor_columns <- c(
        "MS_batch", "EarTag", "Strain",
        "Diet", "digestion_batch", "Sex"
    )
    numeric_columns <- c("DateTime", "order")

    sample_annotation <- example_sample_annotation[, c(factor_columns, numeric_columns)]
    sample_annotation[factor_columns] <- lapply(sample_annotation[factor_columns], factor)

    expect_warning(
        color_scheme <- sample_annotation_to_colors(
            sample_annotation,
            factor_columns = factor_columns,
            numeric_columns = numeric_columns
        ),
        "Too many colors, consider merging some covariate values for better visualisation"
    )

    expect_equal(names(color_scheme), c(factor_columns, numeric_columns))

    expect_equal(color_scheme$MS_batch[[1]], "turquoise")
    expect_equal(color_scheme$MS_batch[[2]], "blue")

    expect_equal(color_scheme$DateTime[[1]], "#A6611A")
    expect_equal(color_scheme$order[[1]], "#D01C8B")
})

test_that("color_list_to_df handles empty intersection", {
    data(example_sample_annotation, package = "proBatch")

    color_list <- list(AnotherFactor = c("red", "blue"))
    expect_warning(
        color_df <- proBatch:::color_list_to_df(color_list, example_sample_annotation),
        "color list and sample annotation have different factors,
            using only intersection in color scheme!"
    )

    expect_equal(nrow(color_df), nrow(example_sample_annotation))
    expect_equal(ncol(color_df), 0)
    expect_equal(rownames(color_df), example_sample_annotation$FullRunName)
})
