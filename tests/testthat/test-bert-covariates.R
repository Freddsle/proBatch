test_that(".sanitize_covariates_for_BERT converts factors to numeric codes", {
    sa <- data.frame(
        sample = paste0("S", 1:4),
        condition = factor(c("A", "B", "A", "C")),
        stringsAsFactors = FALSE
    )

    res <- proBatch:::.sanitize_covariates_for_BERT(sa, "condition")

    expect_equal(res$covariates_cols, "condition")
    expect_true(is.numeric(res$sample_annotation$condition))
    expect_equal(
        res$sample_annotation$condition,
        as.numeric(sa$condition)
    )
})

test_that(".sanitize_covariates_for_BERT handles character covariates", {
    sa <- data.frame(
        sample = paste0("S", 1:3),
        condition = c("ctrl", "case", "ctrl"),
        stringsAsFactors = FALSE
    )

    res <- proBatch:::.sanitize_covariates_for_BERT(sa, "condition")

    expected <- as.numeric(factor(sa$condition, exclude = NULL))
    expect_true(is.numeric(res$sample_annotation$condition))
    expect_equal(res$sample_annotation$condition, expected)
})

test_that(".sanitize_covariates_for_BERT keeps numeric columns numeric", {
    sa <- data.frame(
        sample = paste0("S", 1:3),
        age = c(32.5, 41.1, 41.1)
    )

    res <- proBatch:::.sanitize_covariates_for_BERT(sa, "age")

    expect_true(is.numeric(res$sample_annotation$age))
    expect_equal(res$sample_annotation$age, sa$age)
})

test_that(".sanitize_covariates_for_BERT errors for missing columns", {
    sa <- data.frame(sample = paste0("S", 1:2))

    expect_error(
        proBatch:::.sanitize_covariates_for_BERT(sa, "age"),
        "Covariates missing in sample_annotation"
    )
})

test_that(".sanitize_covariates_for_BERT returns early for NULL covariates", {
    sa <- data.frame(sample = paste0("S", 1:2))

    res <- proBatch:::.sanitize_covariates_for_BERT(sa, NULL)

    expect_null(res$covariates_cols)
    expect_equal(res$sample_annotation, sa)
})
