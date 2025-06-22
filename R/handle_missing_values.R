handle_missing_values <- function(data_matrix, warning_message, fill_the_missing = NULL) {
    # Validate input
    if (!is.matrix(data_matrix)) {
        warning("Coercing input to matrix")
        data_matrix <- as.matrix(data_matrix)
    }
    orig <- data_matrix
    # Only proceed if any NA
    if (any(is.na(orig))) {
        warning(warning_message)

        # 1. Fill path
        if (!is.null(fill_the_missing)) {
            if (!is.numeric(fill_the_missing)) {
                warning("filling value is not numeric, coercing to 0")
                fill_the_missing <- 0
            }
            warning(sprintf("filling missing values with %s", fill_the_missing))
            data_matrix[is.na(data_matrix)] <- fill_the_missing
        } else {
            # 2. Removal path
            nr <- nrow(data_matrix)
            nc <- ncol(data_matrix)
            if (nr == nc) {
                # Square: check symmetry safely
                if (all(!is.na(data_matrix)) && isSymmetric(data_matrix)) {
                    message("removing rows and columns with missing values, as matrix is square")
                    # detect rows entirely NA
                    all_na_row <- apply(data_matrix, 1, function(r) all(is.na(r)))
                    if (any(all_na_row)) {
                        warning("removing rows (and cols) with all values missing")
                        good <- which(!all_na_row)
                        data_matrix <- data_matrix[good, good]
                    }
                } else if (isSymmetric(data_matrix, na.rm = TRUE)) {
                    warning("matrix symmetric once NAs are ignored—no removal performed")
                } else {
                    warning("matrix is square but not symmetric; removing rows with any NAs")
                    data_matrix <- data_matrix[complete.cases(data_matrix), ]
                }
            } else {
                warning("non-square matrix; removing rows with any NAs")
                data_matrix <- data_matrix[complete.cases(data_matrix), ]
            }
        }

        # 3. Report removals
        post <- data_matrix
        # rows removed
        removed_rows <- nrow(orig) - nrow(post)
        if (removed_rows > 0) {
            warning(sprintf("removed %d rows", removed_rows))
        }
        # cols removed
        removed_cols <- ncol(orig) - ncol(post)
        if (removed_cols > 0) {
            warning(sprintf("removed %d columns", removed_cols))
        }
    }
    return(data_matrix)
}
