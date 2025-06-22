test_that("numeric fill replaces missing values", {
    data(example_proteome_matrix, package = "proBatch")
    mat <- example_proteome_matrix[8:10, 1:3]
    index_missing <- which(is.na(mat), arr.ind = TRUE)

    expect_warning(
        expect_warning(
            res <- handle_missing_values(mat, "warn", fill_the_missing = -1),
            "warn"
        ),
        "filling missing values with"
    )

    expect_true(!any(is.na(res)))
    expect_equal(res[index_missing], rep(-1, sum(is.na(mat))))
})


test_that("rows with NAs removed for rectangular matrix", {
    mat <- matrix(c(1, 2, NA, 3, 4, 5), nrow = 3, byrow = TRUE)
    expect_warning(
        expect_warning(
            res <- handle_missing_values(mat, "warn"),
            "warn",
            fixed = TRUE
        ),
        "removed 1 rows",
        fixed = TRUE
    )
    expect_equal(nrow(res), 2)
    expect_equal(res[, 1], c(1, 4))
})

test_that("symmetric square matrix removes rows and columns with NAs", {
    mat <- matrix(
        c(
            1, 2, 3,
            2, 1, NA,
            3, NA, 1
        ),
        nrow = 3, byrow = TRUE,
        dimnames = list(paste0("r", 1:3), paste0("c", 1:3))
    )
    expect_warning(
        expect_warning(
            expect_warning(
                res <- handle_missing_values(mat, "warn"),
                "removed 2 rows",
                fixed = TRUE
            ),
            "matrix is square, but not symmetric",
            fixed = TRUE
        ),
        "warn",
        fixed = TRUE
    )

    expect_equal(dim(res), c(1, 3))
    expect_equal(res, matrix(c(1, 2, 3), nrow = 1, dimnames = list("r1", paste0("c", 1:3))))
})


test_that("square but non-symmetric matrix removes rows only", {
    mat <- matrix(c(1, 2, 3, NA), nrow = 2, byrow = TRUE)
    expect_warning(
        expect_warning(
            expect_warning(
                res <- handle_missing_values(mat, "warn"),
                "matrix is square, but not symmetric",
                fixed = TRUE
            ),
            "removed 1 rows",
            fixed = TRUE
        ),
        "warn",
        fixed = TRUE
    )
    expect_equal(dim(res), c(1, 2))
    expect_true(all(!is.na(res)))
})


test_that("all rows incomplete leads to empty matrix", {
    mat <- matrix(
        c(
            NA, 2, 3,
            2, NA, 4,
            3, 4, NA
        ),
        nrow = 3, byrow = TRUE,
        dimnames = list(paste0("r", 1:3), paste0("c", 1:3))
    )

    expect_warning(
        expect_warning(
            expect_warning(
                res <- handle_missing_values(mat, "warn"),
                "matrix is square, but not symmetric",
                fixed = TRUE
            ),
            "warn",
            fixed = TRUE
        ),
        "removed 3 rows",
        fixed = TRUE
    )
    expect_equal(dim(res), c(0, 3))
})


test_that("non-numeric fill replaces missing values with 0", {
    data(example_proteome_matrix, package = "proBatch")
    mat <- example_proteome_matrix[8:10, 1:3]

    expect_true(any(is.na(mat)))

    expect_warning(
        expect_warning(
            expect_warning(
                res <- handle_missing_values(mat, "warn", fill_the_missing = "a"),
                "filling value is not numeric",
                fixed = TRUE
            ),
            "warn",
            fixed = TRUE
        ),
        "filling missing values with 0",
    )

    expect_true(!any(is.na(res)))
    expect_equal(res[is.na(mat)], rep(0, sum(is.na(mat))))
})
