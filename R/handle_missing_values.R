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


####################################################
# filters for missing values
####################################################


#' Apply `QFeatures` missing-data helpers to stored assays
#'
#' These wrappers delegate to the corresponding `QFeatures` generics while
#' ensuring that the requested assays remain part of the `ProBatchFeatures`
#' object. Only assays that are already materialised can be modified.
#' If a transformation step was applied as a "fast" step (log, log2, etc.),
#' consider re-running it with `store_fast_steps = TRUE`.
#'
#' @param object A `ProBatchFeatures` object.
#' @param pbf_name Character vector of assay names. Defaults to
#'   `pb_current_assay(object)`.
#' @param inplace Logical (used by `pb_filterNA()`, `pb_infIsNA()` and `pb_zeroIsNA()`), whether to modify
#'   the object in place. Default: `FALSE`.
#'   If `FALSE`, the modified assay(s) will be added to the object with
#'   `final_name` (if provided) or the original name(s) with suffix `_filteredNA`.
#' @param final_name Character (used by  by `pb_filterNA()`, `pb_infIsNA()` and `pb_zeroIsNA()`), name for the
#'   modified assay(s) if
#'   `inplace` is `FALSE`. If `NULL` (default), the original name(s) with
#'   suffix `_filteredNA` will be used.
#' @param ... Additional parameters forwarded to the underlying
#'   `QFeatures` method where applicable.
#' @return `pb_zeroIsNA()`, `pb_infIsNA()` and `pb_filterNA()` return the
#'   updated `ProBatchFeatures` object. `pb_nNA()` returns the output of the
#'   corresponding `QFeatures::nNA()` call (a `list` of `DataFrame`s).
#' @name pb_missing_helpers
NULL

#' @rdname pb_missing_helpers
#' @export
pb_zeroIsNA <- function(
    object, 
    pbf_name = pb_current_assay(object),
    inplace = FALSE,
    final_name = NULL,
    ...) {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    stopifnot(is.logical(inplace), length(inplace) == 1L)

    assays <- .pb_require_materialised_assays(object, pbf_name)
    params <- list(...)
    if (!inplace) {
        final_name <- .get_final_name(assays, final_name, suffix = "_zeroIsNA")
    } else if (!is.null(final_name)) {
        warning("`final_name` is ignored when `inplace = TRUE`.")
    }

    for (idx in seq_along(assays)) {
        nm <- assays[[idx]]
        filtered <- do.call(QFeatures::zeroIsNA, c(list(object[[nm]]), params))
        if (inplace) {
            object[[nm]] <- filtered
            to_nm <- nm
        } else {
            new_nm <- final_name[[idx]]
            existing <- names(object)
            if (new_nm %in% existing) {
                new_nm <- make.unique(c(existing, new_nm))[length(existing) + 1L]
            }
            object <- QFeatures::addAssay(object, filtered, name = new_nm)
            to_nm <- new_nm
        }
    }
    object
}

#' @rdname pb_missing_helpers
#' @export
pb_infIsNA <- function(
    object, 
    pbf_name = pb_current_assay(object),
    inplace = FALSE,
    final_name = NULL,
    ...) {    
    stopifnot(methods::is(object, "ProBatchFeatures"))
    stopifnot(is.logical(inplace), length(inplace) == 1L)

    assays <- .pb_require_materialised_assays(object, pbf_name)
    params <- list(...)
    if (!inplace) {
        final_name <- .get_final_name(assays, final_name, suffix = "_infIsNA")
    } else if (!is.null(final_name)) {
        warning("`final_name` is ignored when `inplace = TRUE`.")
    }

    for (idx in seq_along(assays)) {
        nm <- assays[[idx]]
        filtered <- do.call(QFeatures::infIsNA, c(list(object[[nm]]), params))
        if (inplace) {
            object[[nm]] <- filtered
            to_nm <- nm
        } else {
            new_nm <- final_name[[idx]]
            existing <- names(object)
            if (new_nm %in% existing) {
                new_nm <- make.unique(c(existing, new_nm))[length(existing) + 1L]
            }
            object <- QFeatures::addAssay(object, filtered, name = new_nm)
            to_nm <- new_nm
        }
    }
    object
}

#' @rdname pb_missing_helpers
#' @export
pb_nNA <- function(object, pbf_name = pb_current_assay(object), ...) {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    assays <- .pb_require_materialised_assays(object, pbf_name)
    params <- list(...)
    res <- lapply(assays, function(nm) do.call(QFeatures::nNA, c(list(object[[nm]]), params)))
    if (length(res) == 1L) {
        return(res[[1L]])
    }
    stats::setNames(res, assays)
}

#' @rdname pb_missing_helpers
#' @export
pb_filterNA <- function(
    object,
    pbf_name = pb_current_assay(object),
    inplace = FALSE,
    final_name = NULL,
    ...) {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    stopifnot(is.logical(inplace), length(inplace) == 1L)
    assays <- .pb_require_materialised_assays(object, pbf_name)
    params <- list(...)

    if (!inplace) {
        final_name <- .get_final_name(assays, final_name, suffix = "_filteredNA")
    } else if (!is.null(final_name)) {
        warning("`final_name` is ignored when `inplace = TRUE`.")
    }

    for (idx in seq_along(assays)) {
        nm <- assays[[idx]]
        filtered <- do.call(QFeatures::filterNA, c(list(object[[nm]]), params))
        if (inplace) {
            object[[nm]] <- filtered
            to_nm <- nm
        } else {
            new_nm <- final_name[[idx]]
            existing <- names(object)
            if (new_nm %in% existing) {
                new_nm <- make.unique(c(existing, new_nm))[length(existing) + 1L]
            }
            object <- QFeatures::addAssay(object, filtered, name = new_nm)
            to_nm <- new_nm
        }
        object <- .pb_add_log_entry(
            object,
            step = "filterNA",
            fun = "filterNA",
            from = nm,
            to = to_nm,
            params = c(params, list(inplace = inplace))
        )
    }
    object
}

# Internal helper to validate that assays are materialised in the object
.pb_require_materialised_assays <- function(object, assays) {
    if (missing(assays) || is.null(assays)) {
        assays <- pb_current_assay(object)
    }
    assays <- as.character(assays)
    assays <- unique(assays[!is.na(assays) & nzchar(assays)])
    if (!length(assays)) {
        stop(
            "No assay names available. Provide `pbf_name` or ensure the object stores assays.",
            call. = FALSE
        )
    }

    missing <- setdiff(assays, names(object))
    if (length(missing)) {
        extra_msg <- ""
        log_df <- tryCatch(get_operation_log(object), error = function(e) NULL)
        if (!is.null(log_df) && nrow(log_df)) {
            recorded <- unique(as.character(log_df$to))
            logged_only <- intersect(missing, recorded)
            if (length(logged_only)) {
                extra_msg <- paste0(
                    " These assays only exist as logged fast steps. ",
                    "Re-run the originating transformation with `store_fast_steps = TRUE` ",
                    "or `store_intermediate = TRUE` to materialise them."
                )
            }
        }
        stop(
            "Assay(s) ", paste(missing, collapse = ", "),
            " are not stored in the object.", extra_msg,
            call. = FALSE
        )
    }
    assays
}


.get_final_name <- function(original, final_name, suffix) {
    if (is.null(final_name)) {
        final_name <- paste0(original, suffix)
    }
    if (length(final_name) == 1L && length(assays) > 1L) {
            final_name <- rep(final_name, length(assays))
        }
    if (length(final_name) != length(assays)) {
        stop(
            "`final_name` must be length 1 or match `pbf_name` when `inplace = FALSE`.",
            call. = FALSE
        )
    }
    final_name <- as.character(final_name)
    final_name
}