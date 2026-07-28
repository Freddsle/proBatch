test_that("numeric fill replaces missing values", {
    data(example_proteome_matrix, package = "proBatch")
    mat <- example_proteome_matrix[8:10, 1:3]
    index_missing <- which(is.na(mat), arr.ind = TRUE)

    pb_test_expect_warnings(
        res <- handle_missing_values(mat, "warn", fill_the_missing = -1),
        c("warn", "filling missing values with"),
        fixed = TRUE
    )

    expect_true(!any(is.na(res)))
    expect_equal(res[index_missing], rep(-1, sum(is.na(mat))))
})


test_that("rows with NAs removed for rectangular matrix", {
    mat <- matrix(c(1, 2, NA, 3, 4, 5), nrow = 3, byrow = TRUE)
    pb_test_expect_warnings(
        res <- handle_missing_values(mat, "warn"),
        c("warn", "removed 1 rows and 0 columns"),
        fixed = TRUE
    )
    expect_equal(nrow(res), 2)
    expect_equal(res[, 1], c(1, 4))
})

test_that("symmetric matrix removes matching rows and columns", {
    mat <- matrix(
        c(
            1, 2, NA,
            2, 3, 4,
            NA, 4, 5
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("r", 1:3), paste0("r", 1:3))
    )

    pb_test_expect_warnings(
        res <- handle_missing_values(mat, "warn"),
        c("warn", "removed 2 rows and 2 columns"),
        fixed = TRUE
    )

    expect_identical(
        res,
        matrix(3, dimnames = list("r2", "r2"))
    )
})

test_that("square matrix with unmatched dimensions removes rows only", {
    mat <- matrix(
        c(
            1, 2, 3,
            2, 1, NA,
            3, NA, 1
        ),
        nrow = 3, byrow = TRUE,
        dimnames = list(paste0("r", 1:3), paste0("c", 1:3))
    )
    pb_test_expect_warnings(
        res <- handle_missing_values(mat, "warn"),
        c(
            "removed 2 rows and 0 columns",
            "Matrix is square but not symmetric",
            "warn"
        ),
        fixed = TRUE
    )

    expect_equal(dim(res), c(1, 3))
    expect_equal(res, matrix(c(1, 2, 3), nrow = 1, dimnames = list("r1", paste0("c", 1:3))))
})


test_that("square but non-symmetric matrix removes rows only", {
    mat <- matrix(c(1, 2, 3, NA), nrow = 2, byrow = TRUE)
    pb_test_expect_warnings(
        res <- handle_missing_values(mat, "warn"),
        c(
            "Matrix is square but not symmetric",
            "removed 1 rows and 0 columns",
            "warn"
        ),
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

    pb_test_expect_warnings(
        res <- handle_missing_values(mat, "warn"),
        c(
            "Matrix is square but not symmetric",
            "warn",
            "removed 3 rows and 0 columns"
        ),
        fixed = TRUE
    )
    expect_equal(dim(res), c(0, 3))
})


test_that("non-numeric fill replaces missing values with 0", {
    data(example_proteome_matrix, package = "proBatch")
    mat <- example_proteome_matrix[8:10, 1:3]

    expect_true(any(is.na(mat)))

    pb_test_expect_warnings(
        res <- handle_missing_values(mat, "warn", fill_the_missing = "a"),
        c(
            "filling value is not a finite numeric scalar",
            "warn",
            "filling missing values with 0"
        ),
        fixed = TRUE
    )

    expect_true(!any(is.na(res)))
    expect_equal(res[is.na(mat)], rep(0, sum(is.na(mat))))
})

test_that("false keeps missing values unchanged", {
    mat <- matrix(
        c(1, NA, 3, 4),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )

    pb_test_expect_warnings(
        res <- handle_missing_values(
            mat,
            "warn",
            fill_the_missing = FALSE
        ),
        c(
            "warn",
            "`fill_the_missing` is FALSE: keeping missing values unchanged"
        ),
        fixed = TRUE
    )

    expect_identical(res, mat)
})
