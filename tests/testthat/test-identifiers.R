test_that("shared identifier validation rejects missing and duplicate IDs", {
    expect_error(
        proBatch:::.pb_validate_identifiers(NULL, "Sample axis"),
        "must be named"
    )
    expect_error(
        proBatch:::.pb_validate_identifiers(c("s1", ""), "Sample axis"),
        "NA or empty"
    )
    expect_error(
        proBatch:::.pb_validate_identifiers(c("s1", "s1"), "Sample axis"),
        "duplicate identifiers.*s1"
    )
    expect_identical(
        proBatch:::.pb_validate_identifiers(
            c("s1", "s1"),
            "Long samples",
            require_unique = FALSE
        ),
        c("s1", "s1")
    )
})

test_that("long-key validation rejects duplicate feature/sample pairs", {
    duplicate <- data.frame(
        feature = c("f1", "f1"),
        sample = c("s1", "s1"),
        value = c(1, 2)
    )

    expect_error(
        proBatch:::.pb_validate_long_keys(
            duplicate,
            feature_id_col = "feature",
            sample_id_col = "sample"
        ),
        "duplicate feature/sample keys.*f1/s1"
    )
})
