test_that("PVCA stacked-summary symbols are declared globally", {
    check_root <- normalizePath(
        file.path(testthat::test_path(), "..", ".."),
        mustWork = TRUE
    )
    candidates <- c(
        check_root,
        file.path(check_root, "00_pkg_src", "proBatch")
    )
    valid <- vapply(
        candidates,
        function(candidate) {
            file.exists(file.path(candidate, "DESCRIPTION")) &&
                file.exists(file.path(candidate, "R", "proBatch.R"))
        },
        logical(1L)
    )
    if (!any(valid)) {
        skip("Global-variable checks require the proBatch source package")
    }

    source_root <- normalizePath(
        candidates[[which(valid)[[1L]]]],
        mustWork = TRUE
    )
    expressions <- parse(file.path(source_root, "R", "proBatch.R"))
    global_calls <- list()
    collect_global_calls <- function(expression) {
        if (
            is.call(expression) &&
                identical(
                    paste(deparse(expression[[1L]]), collapse = ""),
                    "utils::globalVariables"
                )
        ) {
            global_calls[[length(global_calls) + 1L]] <<- expression
        }
        if (is.recursive(expression)) {
            lapply(as.list(expression), collect_global_calls)
        }
        invisible(NULL)
    }
    lapply(as.list(expressions), collect_global_calls)

    expect_length(global_calls, 1L)
    declared <- eval(global_calls[[1L]][[2L]], envir = baseenv())
    expect_contains(declared, "weight")
})
