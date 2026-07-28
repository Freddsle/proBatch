test_that("missing-value policy defaults to an explicit error", {
    matrix <- matrix(
        c(1, NA_real_, 3, 4),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )

    expect_error(
        handle_missing_values(
            matrix,
            "method cannot operate with missing values"
        ),
        "method cannot operate with missing values"
    )
    expect_error(
        handle_missing_values(
            matrix,
            "unused",
            fill_the_missing = NULL
        ),
        "NULL.*ambiguous"
    )

    complete <- matrix
    complete[is.na(complete)] <- 0
    expect_error(
        handle_missing_values(
            complete,
            "unused",
            fill_the_missing = NULL
        ),
        "NULL.*ambiguous"
    )
})

test_that("canonical keep, drop_features, and fill policies are deterministic", {
    matrix <- matrix(
        c(1, NA_real_, 3, 4, 5, 6),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(c("f1", "f2", "f3"), c("s1", "s2"))
    )

    kept <- expect_silent(handle_missing_values(
        matrix,
        "unused",
        fill_the_missing = "keep"
    ))
    expect_identical(kept, matrix)

    dropped <- pb_test_expect_warnings(
        handle_missing_values(
            matrix,
            "unused",
            fill_the_missing = "drop_features"
        ),
        "removed 1 rows and 0 columns",
        fixed = TRUE
    )
    expect_identical(rownames(dropped), c("f2", "f3"))
    expect_identical(colnames(dropped), colnames(matrix))

    filled <- expect_silent(handle_missing_values(
        matrix,
        "unused",
        fill_the_missing = "fill",
        fill_value = -1
    ))
    expect_false(anyNA(filled))
    expect_identical(unname(filled[1, 2]), -1)
})

test_that("drop_features never removes matching sample columns", {
    matrix <- matrix(
        c(
            1, 2, 3,
            2, 1, NA,
            3, NA, 1
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:3))
    )

    result <- pb_test_expect_warnings(
        handle_missing_values(
            matrix,
            "unused",
            fill_the_missing = "drop_features"
        ),
        "removed 2 rows and 0 columns",
        fixed = TRUE
    )

    expect_identical(dim(result), c(1L, 3L))
    expect_identical(rownames(result), "f1")
    expect_identical(colnames(result), colnames(matrix))
})

test_that("drop_features reports only removals that occur", {
    complete <- matrix(
        1:6,
        nrow = 3,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:2))
    )
    expect_identical(
        expect_silent(handle_missing_values(
            complete,
            "unused",
            fill_the_missing = "drop_features"
        )),
        complete
    )

    incomplete <- complete
    incomplete[, 1] <- NA_real_
    result <- pb_test_expect_warnings(
        handle_missing_values(
            incomplete,
            "unused",
            fill_the_missing = "drop_features"
        ),
        "removed 3 rows and 0 columns",
        fixed = TRUE
    )
    expect_identical(dim(result), c(0L, 2L))
})

test_that("legacy missing values translate with one deprecation warning", {
    matrix <- matrix(
        c(1, NA_real_, 3, 4),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )

    kept <- pb_test_expect_warnings(
        handle_missing_values(
            matrix,
            "unused",
            fill_the_missing = FALSE
        ),
        "deprecated.*keep"
    )
    expect_true(anyNA(kept))

    filled <- pb_test_expect_warnings(
        handle_missing_values(
            matrix,
            "unused",
            fill_the_missing = 0
        ),
        "deprecated.*fill"
    )
    expect_identical(unname(filled[1, 2]), 0)

    for (legacy in c("remove", "rm", "REMOVE")) {
        dropped <- pb_test_expect_warnings(
            handle_missing_values(
                matrix,
                "unused",
                fill_the_missing = legacy
            ),
            c("deprecated.*drop_features", "removed 1 rows and 0 columns")
        )
        expect_identical(rownames(dropped), "f2")
    }
})

test_that("invalid policies and fill values are rejected", {
    matrix <- matrix(
        c(1, NA_real_, 3, 4),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )

    for (invalid in list(TRUE, "anything", "d", c("keep", "fill"))) {
        expect_error(
            handle_missing_values(
                matrix,
                "unused",
                fill_the_missing = invalid
            ),
            "must be one of"
        )
    }
    expect_error(
        handle_missing_values(
            matrix,
            "unused",
            fill_the_missing = "fill"
        ),
        "fill_value"
    )
    expect_error(
        handle_missing_values(
            matrix,
            "unused",
            fill_the_missing = "keep",
            fill_value = 0
        ),
        "only valid"
    )
    expect_error(
        handle_missing_values(
            matrix,
            "unused",
            fill_the_missing = 0,
            fill_value = 1
        ),
        "legacy numeric"
    )
    expect_error(
        handle_missing_values(
            matrix,
            "unused",
            fill_the_missing = Inf
        ),
        "finite"
    )
    expect_error(
        handle_missing_values(
            as.data.frame(matrix),
            "unused",
            fill_the_missing = "keep"
        ),
        "numeric matrix"
    )

    character_matrix <- matrix(c("1", NA_character_), nrow = 1)
    expect_error(
        handle_missing_values(
            character_matrix,
            "unused",
            fill_the_missing = "keep"
        ),
        "numeric matrix"
    )
})
