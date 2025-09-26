#' Handle missing values in a data matrix
#'
#' This function can either fill missing values with a specified value
#' or remove rows (and columns, if applicable) with missing values.
#' It is primarily intended for use prior to batch correction methods
#' that cannot handle missing values, such as ComBat or limma's
#' removeBatchEffect, or plotting functions that require complete data.
#'
#' Semantics:
#' - If there are no NAs: return input unchanged.
#' - If fill_the_missing is explicitly FALSE: do nothing (keep NAs).
#' - If fill_the_missing is missing (argument not supplied) or one of
#'   "remove","rm","REMOVE": remove rows with any NA; if the matrix is
#'   square and symmetric (na.rm=TRUE), remove matching rows AND columns
#'   using the same row keep-mask.
#' - Otherwise: if non-numeric or NA, coerce to 0 with a warning; then fill NAs.
#' @param data_matrix A numeric matrix with features in rows and samples in columns.
#' @param warning_message A character string with a warning shown if missing values are found.
#' @param fill_the_missing A control value:
#'   - FALSE: do nothing (keep NAs).
#'   - Missing (arg not supplied) or "remove"/"rm"/"REMOVE": remove rows with any NA
#'     (and matching columns if square & symmetric).
#'   - Numeric scalar: fill NAs with this value.
#'   - Non-numeric: coerced to 0 with a warning and used to fill NAs.
#' @return A matrix with missing values handled as specified.
#' @keywords internal
handle_missing_values <- function(data_matrix, warning_message, fill_the_missing = NULL) {
    # 1) Validate input and coerce to matrix
    if (!is.matrix(data_matrix)) {
        warning("Coercing input to matrix")
        data_matrix <- as.matrix(data_matrix)
    }
    orig <- data_matrix

    # 2) If no NAs, return early
    if (!anyNA(orig)) {
        return(data_matrix)
    }

    # There ARE NAs
    warning(warning_message)

    # 3) EXPLICIT no-op branch: FALSE => keep NAs, do nothing
    if (isFALSE(fill_the_missing)) {
        warning("`fill_the_missing` is FALSE: keeping missing values unchanged")
        return(data_matrix)
    }

    # 4) REMOVAL branch: old default (missing arg) OR explicit "remove" token
    removal_tokens <- c("remove", "rm", "REMOVE")
    removal_requested <- missing(fill_the_missing) || is.null(fill_the_missing) ||
        (is.character(fill_the_missing) && length(fill_the_missing) == 1L &&
            fill_the_missing %in% removal_tokens)

    if (removal_requested) {
        nr <- nrow(data_matrix)
        nc <- ncol(data_matrix)

        # Define row-wise completeness
        keep_rows <- complete.cases(data_matrix)

        if (nr == nc && isSymmetric(data_matrix, na.rm = TRUE)) {
            # Square & (NA-tolerant) symmetric: remove rows AND the corresponding columns
            if (!any(keep_rows)) {
                warning("All rows contain missing values; returning 0x0 matrix")
                data_matrix <- data_matrix[FALSE, FALSE, drop = FALSE]
            } else {
                message("Removing rows and corresponding columns with missing values (square & symmetric)")
                data_matrix <- data_matrix[keep_rows, keep_rows, drop = FALSE]
            }
        } else {
            # Non-square OR square but not symmetric: remove only rows with any NA
            if (nr == nc && !isSymmetric(data_matrix, na.rm = TRUE)) {
                warning("Matrix is square but not symmetric; removing rows with missing values")
            } else {
                message("Removing rows with missing values")
            }
            data_matrix <- data_matrix[keep_rows, , drop = FALSE]
        }

        # 5) Report removals
        post <- data_matrix
        removed_rows <- nrow(orig) - nrow(post)
        if (removed_rows > 0) warning(sprintf("removed %d rows", removed_rows))
        removed_cols <- ncol(orig) - ncol(post)
        if (removed_cols > 0) warning(sprintf("removed %d columns", removed_cols))

        # Filling doesn't remove rows/cols; still report zeros to mirror previous behavior
        message(
            paste(
                sprintf("removed %d rows", removed_rows)
            ),
            sprintf("and %d columns", removed_cols)
        )
        return(data_matrix)
    }

    # 6) any other specified value
    fill_val <- fill_the_missing
    if (!is.numeric(fill_val) || length(fill_val) != 1L || is.na(fill_val)) {
        warning("filling value is not a finite numeric scalar; coercing to 0")
        fill_val <- 0
    }
    warning(sprintf("filling missing values with %s", as.character(fill_val)))
    nas <- is.na(data_matrix)
    if (any(nas)) {
        data_matrix[nas] <- fill_val
    }

    message(paste0(
        sprintf("replaced values in %d rows", sum(rowSums(nas) > 0)),
        sprintf("and %d columns", sum(colSums(nas) > 0))
    ))
    return(data_matrix)
}
