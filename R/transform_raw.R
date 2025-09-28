#' Functions to log transform raw data before normalization and batch correction
#'
#' @inheritParams proBatch
#' @param log_base base of the logarithm for transformation
#' @param offset small positive number to prevent 0 conversion to \code{-Inf}
#'
#' @return `log_transform_df()` returns \code{df_long}-size data frame, with
#' \code{measure_col} log transformed; with old value in another column
#' called "beforeLog_intensity" if "intensity" was the value of
#' \code{measure_col};
#' `log_transform_dm()` returns \code{data_matrix} format matrix
#'
#' @examples
#' data(list = c("example_proteome", "example_proteome_matrix"), package = "proBatch")
#'
#' log_transformed_df <- log_transform_df(example_proteome)
#'
#' log_transformed_matrix <- log_transform_dm(example_proteome_matrix,
#'     log_base = 10, offset = 1
#' )
#' @name transform_raw_data
NULL

#' Log transformation of the data long format.
#' @rdname transform_raw_data
#' @export
log_transform_df <- function(df_long, log_base = 2, offset = 1,
                             measure_col = "Intensity") {
    if (!is.null(log_base)) {
        df_long <- df_long %>%
            mutate(!!(paste("beforeLog", measure_col, sep = "_")) :=
                !!(sym(measure_col))) %>%
            mutate(!!(sym(measure_col)) :=
                log(!!(sym(measure_col)) + offset, base = log_base))
    } else {
        warning("Log base is NULL, returning the original data frame")
    }
    return(df_long)
}


#' "Unlog" transformation of the data to pre-log form (for quantification, forcing log-transform)
#'
#' @export
#' @rdname transform_raw_data
#'
unlog_df <- function(df_long, log_base = 2, offset = 1,
                     measure_col = "Intensity") {
    if (!is.null(log_base)) {
        df_long <- df_long %>%
            mutate(!!(paste("beforeUnLog", measure_col, sep = "_")) :=
                !!(sym(measure_col))) %>%
            mutate(!!(sym(measure_col)) :=
                log_base^(!!sym(measure_col)) - offset)
    } else {
        warning("Log base is NULL, returning the original data frame")
    }
    return(df_long)
}

#' @export
#' @rdname transform_raw_data
#'
log_transform_dm <- function(x, ...) UseMethod("log_transform_dm")

#' Log transformation of the data matrix format.
#' @method log_transform_dm default
#' @export
#' @rdname transform_raw_data
log_transform_dm.default <- function(data_matrix, log_base = 2, offset = 1) {
    if (!is.null(log_base)) {
        # Validate numeric data_matrix
        if (!is.numeric(data_matrix)) {
            stop("data_matrix must be numeric")
        }
        data_matrix_log <- log(data_matrix + offset, base = log_base)
    } else {
        warning("Log base is NULL, returning the original data matrix")
        data_matrix_log <- data_matrix
    }
    return(data_matrix_log)
}

#' @rdname transform_raw_data
#' @method log_transform_dm ProBatchFeatures
#' @export
log_transform_dm.ProBatchFeatures <- function(x, log_base = 2, offset = 1,
                                              pbf_name = NULL, final_name = NULL) {
    object <- x
    if (is.null(pbf_name)) {
        pbf_name <- pb_current_assay(object)
        message("`pbf_name` not provided, using the most recent assay: ", pbf_name)
    }
    step <- if (!is.null(log_base) && log_base == 2 && offset == 1) "log2" else "log"
    object <- pb_transform(
        object,
        from = pbf_name,
        steps = step,
        funs = list(log_transform_dm.default),
        params_list = list(list(log_base = log_base, offset = offset)),
        final_name = final_name
    )
    object
}

#' @export
#' @rdname transform_raw_data
#'
unlog_dm <- function(x, ...) UseMethod("unlog_dm")

#' @rdname transform_raw_data
#' @method unlog_dm default
#' @export
unlog_dm.default <- function(data_matrix, log_base = 2, offset = 1) {
    if (!is.null(log_base)) {
        data_matrix_unlog <- log_base^(data_matrix) - offset
    } else {
        warning("Log base is NULL, returning the original data matrix")
        data_matrix_unlog <- data_matrix
    }
    return(data_matrix_unlog)
}

#' @rdname transform_raw_data
#' @method unlog_dm ProBatchFeatures
#' @export
unlog_dm.ProBatchFeatures <- function(x, log_base = 2, offset = 1,
                                      pbf_name = NULL, final_name = NULL) {
    object <- x
    if (is.null(pbf_name)) {
        pbf_name <- pb_current_assay(object)
        message("`pbf_name` not provided, using the most recent assay: ", pbf_name)
    }
    object <- pb_transform(
        object,
        from = pbf_name,
        steps = "unlog",
        funs = list(unlog_dm.default),
        params_list = list(list(log_base = log_base, offset = offset)),
        final_name = final_name
    )
    object
}
