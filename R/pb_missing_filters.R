####################################################
# filters for missing values in ProBatchFeatures
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
