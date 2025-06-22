handle_missing_values <- function(data_matrix, warning_message, fill_the_missing = NULL) {
    # 1. Validate input and coerce to matrix
    if (!is.matrix(data_matrix)) {
        warning("Coercing input to matrix")
        data_matrix <- as.matrix(data_matrix)
    }
    orig <- data_matrix

    # 2. Only proceed if there are any NAs
    if (any(is.na(orig))) {
        warning(warning_message)

        # 2a. Fill-missing path
        if (!is.null(fill_the_missing)) {
            if (!is.numeric(fill_the_missing)) {
                warning("filling value is not numeric, coercing to 0")
                fill_the_missing <- 0
            }
            warning(sprintf("filling missing values with %s", fill_the_missing))
            data_matrix[is.na(data_matrix)] <- fill_the_missing

            # 2b. Removal path
        } else {
            nr <- nrow(data_matrix)
            nc <- ncol(data_matrix)

            # 2b.i Square & symmetric (ignoring NAs)
            if (nr == nc && isSymmetric(data_matrix, na.rm = TRUE)) {
                message("removing rows and columns with missing values, as matrix is square")
                keep <- complete.cases(data_matrix)
                if (!any(keep)) {
                    warning("removing rows with all values missing")
                    data_matrix <- data_matrix[FALSE, FALSE]
                } else {
                    data_matrix <- data_matrix[keep, keep, drop = FALSE]
                }

                # 2b.ii Square but not symmetric
            } else if (nr == nc) {
                warning("matrix is square, but not symmetric, coincidence or error?")
                data_matrix <- data_matrix[complete.cases(data_matrix), , drop = FALSE]

                # 2b.iii Non-square: remove rows with any NA
            } else {
                data_matrix <- data_matrix[complete.cases(data_matrix), , drop = FALSE]
            }
        }

        # 3. Report removals
        post <- data_matrix
        removed_rows <- nrow(orig) - nrow(post)
        if (removed_rows > 0) {
            warning(sprintf("removed %d rows", removed_rows))
        }
        removed_cols <- ncol(orig) - ncol(post)
        if (removed_cols > 0) {
            warning(sprintf("removed %d columns", removed_cols))
        }
    }

    return(data_matrix)
}
