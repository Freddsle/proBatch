<<<<<<< HEAD
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
#' @examples
#' mat <- matrix(c(1, NA, 3, 4), nrow = 2)
#' suppressWarnings(handle_missing_values(mat, warning_message = "demo", fill_the_missing = 0))
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
        message(sprintf(
            "removed %d rows and %d columns",
            removed_rows,
            removed_cols
        ))
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

    message(sprintf(
        "replaced values in %d rows and %d columns",
        sum(rowSums(nas) > 0),
        sum(colSums(nas) > 0)
    ))
=======
handle_missing_values <- function(data_matrix, warning_message,
                                  fill_the_missing = NULL) {
    if (any(is.na(as.vector(data_matrix)))) {
        warning(warning_message)
        pre_corr_dim <- nrow(data_matrix)
        if (!is.null(fill_the_missing)) {
            if (!is.numeric(fill_the_missing)) {
                warning("filling value is not numeric, setting to 0")
                fill_the_missing <- 0
            } else {
                warning(sprintf("filling missing value with %s", fill_the_missing))
                data_matrix[is.na(data_matrix)] <- fill_the_missing
            }
        } else {
            complete_cases <- complete.cases(data_matrix)
            if (ncol(data_matrix) == nrow(data_matrix)) {
                pre_corr_dim <- nrow(data_matrix)
                if (isSymmetric(data_matrix)) {
                    message("removing rows and columns with missing values, as matrix is square")
                    if (all(!complete_cases)) {
                        warning("removing rows with all values missing")
                        all_missing_rows <- apply(data_matrix, 2, function(x) all(is.na(x)))
                        bad_rows <- names(all_missing_rows[all_missing_rows])
                        good_rows <- setdiff(colnames(data_matrix), bad_rows)
                        data_matrix <- (data_matrix[good_rows, good_rows])
                    }
                    data_matrix <- data_matrix[complete.cases(data_matrix), complete.cases(data_matrix)]
                } else {
                    warning("matrix is square, but not symmetric, coincidence or error?")
                    data_matrix <- data_matrix[complete.cases(data_matrix), ]
                }
            } else {
                warning("filling value is NULL, removing rows with missing values from data matrix")
                data_matrix <- data_matrix[complete.cases(data_matrix), ]
            }
        }
        after_corr_dim <- nrow(data_matrix)
        rem_rows <- pre_corr_dim - after_corr_dim
        if (rem_rows > 0) {
            warning(sprintf("removed %s rows of the matrix with missing values", rem_rows))
        }
    }
>>>>>>> 7f231190 (4 spaces (BioCheck), added test for  transform log funcs, fixed seed in colors to hex sorting)
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
#'   `names(object)` - all assays.
#' @param inplace Logical (used by `pb_filterNA()` only), whether to modify the
#'   object in place. Default: `FALSE`.
#'   If `FALSE`, the modified assay(s) will be added to the object with
#'   `final_name` (if provided) or the original name(s) with suffix `_filteredNA`.
#' @param final_name Character (used by `pb_filterNA()` only), name for the
#'   modified assay(s) if `inplace` is `FALSE`. If `NULL` (default), the
#'   original name(s) with suffix `_filteredNA` will be used.
#' @param ... Additional parameters forwarded to the underlying
#'   `QFeatures` method where applicable.
#' @return `pb_zeroIsNA()`, `pb_infIsNA()` and `pb_filterNA()` return the
#'   updated `ProBatchFeatures` object. `pb_nNA()` returns the output of the
#'   corresponding `QFeatures::nNA()` call (a `list` of `DataFrame`s).
#' @name pb_missing_helpers
NULL

#' @rdname pb_missing_helpers
#' @export
pb_zeroIsNA <- function(object, pbf_name = names(object), ...) {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    assays <- .pb_require_materialised_assays(object, pbf_name)
    params <- .pb_collect_missing_params(list(...), forbidden = "i")
    for (nm in assays) {
        prior <- object
        object <- do.call(QFeatures::zeroIsNA, c(list(object, i = nm), params))
        object <- .as_ProBatchFeatures(object, from = prior)
        object <- .pb_add_log_entry(
            object,
            step = "zeroIsNA",
            fun = "zeroIsNA",
            from = nm,
            to = nm,
            params = params
        )
    }
    object
}

#' @rdname pb_missing_helpers
#' @export
pb_infIsNA <- function(object, pbf_name = names(object), ...) {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    assays <- .pb_require_materialised_assays(object, pbf_name)
    params <- .pb_collect_missing_params(list(...), forbidden = "i")
    for (nm in assays) {
        prior <- object
        object <- do.call(QFeatures::infIsNA, c(list(object, i = nm), params))
        object <- .as_ProBatchFeatures(object, from = prior)
        object <- .pb_add_log_entry(
            object,
            step = "infIsNA",
            fun = "infIsNA",
            from = nm,
            to = nm,
            params = params
        )
    }
    object
}

#' @rdname pb_missing_helpers
#' @export
pb_nNA <- function(object, pbf_name = names(object), ...) {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    assays <- .pb_require_materialised_assays(object, pbf_name)
    params <- .pb_collect_missing_params(list(...), forbidden = "i")
    res <- lapply(
        assays,
        function(nm) do.call(QFeatures::nNA, c(list(object, i = nm), params))
    )
    if (length(res) == 1L) {
        return(res[[1L]])
    }
    res <- stats::setNames(res, assays)
    # join "nNA" from each assay into a single DataFrame and add it to the result
    # only for "nNA" and add it as a last result entry
    res$nNA <- do.call(rbind, lapply(res, `[[`, "nNA"))
    res
}

#' @rdname pb_missing_helpers
#' @export
pb_filterNA <- function(
    object,
    pbf_name = NULL,
    inplace = FALSE,
    final_name = NULL,
    ...) {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    stopifnot(is.logical(inplace), length(inplace) == 1L)

    if (is.null(pbf_name)) {
        pbf_name <- names(object)
        message("`pbf_name` not provided, using all assays: ", paste(pbf_name, collapse = ", "))
    }

    assays <- .pb_require_materialised_assays(object, pbf_name)
    params <- .pb_collect_missing_params(list(...), forbidden = c("i", "name"))

    if (!inplace) {
        final_name <- .pb_prepare_final_names(assays, final_name, suffix = "_filteredNA")
    } else if (!is.null(final_name)) {
        warning("`final_name` is ignored when `inplace = TRUE`.")
    }

    for (idx in seq_along(assays)) {
        nm <- assays[[idx]]
        message("Processing assay:", nm)
        message("  Features before filtering:\t", length(object[[nm]]))
        if (inplace) {
            prior <- object
            object <- do.call(QFeatures::filterNA, c(list(object, i = nm), params))
            object <- .as_ProBatchFeatures(object, from = prior)
            to_nm <- nm
            message("  Features after filtering:\t", length(object[[nm]]))
        } else {
            filtered_obj <- do.call(QFeatures::filterNA, c(list(object, i = nm), params))
            filtered_obj <- .as_ProBatchFeatures(filtered_obj, from = object)
            filtered <- filtered_obj[[nm]]
            new_nm <- .pb_unique_assay_name(object, final_name[[idx]])
            prior <- object
            object <- QFeatures::addAssay(object, filtered, name = new_nm)
            object <- .as_ProBatchFeatures(object, from = prior)
            to_nm <- new_nm
            message("  Features after filtering:\t", length(object[[new_nm]]))
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


.pb_collect_missing_params <- function(params, forbidden = character()) {
    if (!length(params)) {
        return(params)
    }
    nm <- names(params)
    if (!is.null(nm)) {
        bad <- intersect(nm[nzchar(nm)], forbidden)
        if (length(bad)) {
            stop(
                "Argument(s) ", paste(bad, collapse = ", "),
                " must not be supplied via `...`.",
                call. = FALSE
            )
        }
    }
    params
}

.pb_prepare_final_names <- function(assays, final_name, suffix) {
    if (is.null(final_name)) {
        final_name <- paste0(assays, suffix)
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
    as.character(final_name)
}

.pb_unique_assay_name <- function(object, proposed) {
    existing <- names(object)
    if (!length(existing) || !(proposed %in% existing)) {
        return(proposed)
    }
    make.unique(c(existing, proposed))[length(existing) + 1L]
}
