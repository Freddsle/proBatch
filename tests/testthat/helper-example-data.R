pb_test_load_example_data <- function() {
    data(
        list = c("example_proteome_matrix", "example_sample_annotation"),
        package = "proBatch",
        envir = parent.frame()
    )
    invisible(TRUE)
}

pb_test_expect_warnings <- function(object, patterns, fixed = FALSE) {
    object_expr <- substitute(object)
    caller <- parent.frame()
    warning_messages <- character()

    value <- withCallingHandlers(
        eval(object_expr, envir = caller),
        warning = function(condition) {
            warning_messages <<- c(
                warning_messages,
                conditionMessage(condition)
            )
            invokeRestart("muffleWarning")
        }
    )

    testthat::expect_length(warning_messages, length(patterns))
    for (pattern in patterns) {
        testthat::expect_true(
            any(grepl(pattern, warning_messages, fixed = fixed)),
            info = paste0(
                "Expected warning matching `",
                pattern,
                "`, got: ",
                paste(warning_messages, collapse = " | ")
            )
        )
    }

    invisible(value)
}

pb_test_make_pbf <- function(
    n_rows = 30,
    n_cols = 6,
    add_log2 = FALSE,
    complete = FALSE
) {
    pb_test_load_example_data()

    matrix_small <- example_proteome_matrix[1:n_rows, 1:n_cols]
    if (isTRUE(complete)) {
        matrix_small[is.na(matrix_small)] <- 0
    }
    sample_ids <- colnames(matrix_small)
    sample_ann <- example_sample_annotation[
        match(sample_ids, example_sample_annotation$FullRunName),
    ]

    pbf <- suppressMessages(ProBatchFeatures(
        data_matrix = matrix_small,
        sample_annotation = sample_ann,
        sample_id_col = "FullRunName",
        name = "feature::raw"
    ))

    if (isTRUE(add_log2)) {
        pbf <- suppressMessages(pb_transform(
            pbf,
            from = "feature::raw",
            steps = "log2",
            store_fast_steps = TRUE
        ))
    }

    pbf
}
