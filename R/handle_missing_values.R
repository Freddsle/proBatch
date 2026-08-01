.pb_normalize_missing_policy <- function(
    missing,
    fill_value = NULL,
    argument = "missing"
) {
    choices <- c("error", "keep", "drop_features", "fill")
    legacy <- NULL

    if (is.null(missing)) {
        stop(
            "`",
            argument,
            " = NULL` is ambiguous. Use one of ",
            paste(sprintf('"%s"', choices), collapse = ", "),
            " and supply `fill_value` only with \"fill\".",
            call. = FALSE
        )
    }

    if (isFALSE(missing)) {
        legacy <- "FALSE"
        missing <- "keep"
    } else if (is.numeric(missing)) {
        if (length(missing) != 1L || is.na(missing) || !is.finite(missing)) {
            stop(
                "A legacy numeric `",
                argument,
                "` must be one finite, non-missing value.",
                call. = FALSE
            )
        }
        if (!is.null(fill_value)) {
            stop(
                "Do not supply `fill_value` together with a legacy numeric `",
                argument,
                "`; use `",
                argument,
                " = \"fill\"` instead.",
                call. = FALSE
            )
        }
        legacy <- "numeric"
        fill_value <- missing
        missing <- "fill"
    } else if (
        is.character(missing) &&
            length(missing) == 1L &&
            !is.na(missing) &&
            missing %in% c("remove", "rm", "REMOVE")
    ) {
        legacy <- missing
        missing <- "drop_features"
    } else if (
        !is.character(missing) ||
            length(missing) != 1L ||
            is.na(missing) ||
            !(missing %in% choices)
    ) {
        stop(
            "`",
            argument,
            "` must be one of ",
            paste(sprintf('"%s"', choices), collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    if (!is.null(legacy)) {
        replacement <- if (identical(legacy, "numeric")) {
            paste0(
                "`",
                argument,
                " = \"fill\"` with `fill_value = ",
                format(fill_value),
                "`"
            )
        } else {
            paste0("`", argument, " = \"", missing, "\"`")
        }
        warning(
            "Legacy `",
            argument,
            " = ",
            legacy,
            "` is deprecated; use ",
            replacement,
            ".",
            call. = FALSE
        )
    }

    if (identical(missing, "fill")) {
        if (
            !is.numeric(fill_value) ||
                length(fill_value) != 1L ||
                is.na(fill_value) ||
                !is.finite(fill_value)
        ) {
            stop(
                "`fill_value` must be one finite, non-missing numeric value ",
                "when `",
                argument,
                " = \"fill\"`.",
                call. = FALSE
            )
        }
    } else if (!is.null(fill_value)) {
        stop(
            "`fill_value` is only valid when `",
            argument,
            " = \"fill\"`.",
            call. = FALSE
        )
    }

    list(policy = missing, fill_value = fill_value)
}

#' Handle missing values in a numeric matrix
#'
#' Applies the package-wide missing-value policy used before matrix-oriented
#' batch correction. Feature dropping never removes sample columns, including
#' when the input happens to be square and symmetric.
#'
#' @param data_matrix Numeric feature-by-sample matrix.
#' @param warning_message Error message used when missing values are present
#'   under the default \code{"error"} policy.
#' @param fill_the_missing One of \code{"error"}, \code{"keep"},
#'   \code{"drop_features"}, or \code{"fill"}. For one release, \code{FALSE}, a
#'   numeric scalar, and \code{"remove"}, \code{"rm"}, or \code{"REMOVE"} are
#'   translated with a deprecation warning. Explicit \code{NULL} is an error.
#' @param fill_value Finite numeric scalar used only with
#'   \code{fill_the_missing = "fill"}.
#'
#' @return A numeric matrix with the requested missing-value policy applied.
#' @export
#' @examples
#' mat <- matrix(c(1, NA, 3, 4), nrow = 2)
#' handle_missing_values(
#'     mat,
#'     warning_message = "demo",
#'     fill_the_missing = "fill",
#'     fill_value = 0
#' )
handle_missing_values <- function(
    data_matrix,
    warning_message,
    fill_the_missing = "error",
    fill_value = NULL
) {
    if (!is.matrix(data_matrix) || !is.numeric(data_matrix)) {
        stop("`data_matrix` must be a numeric matrix.", call. = FALSE)
    }

    policy <- .pb_normalize_missing_policy(
        fill_the_missing,
        fill_value = fill_value,
        argument = "fill_the_missing"
    )
    if (!anyNA(data_matrix)) {
        return(data_matrix)
    }

    if (identical(policy$policy, "error")) {
        stop(warning_message, call. = FALSE)
    }
    if (identical(policy$policy, "keep")) {
        return(data_matrix)
    }
    if (identical(policy$policy, "drop_features")) {
        kept <- data_matrix[complete.cases(data_matrix), , drop = FALSE]
        removed_rows <- nrow(data_matrix) - nrow(kept)
        if (removed_rows > 0L) {
            warning(
                sprintf("removed %d rows and 0 columns", removed_rows),
                call. = FALSE
            )
        }
        return(kept)
    }

    data_matrix[is.na(data_matrix)] <- policy$fill_value
    data_matrix
}
