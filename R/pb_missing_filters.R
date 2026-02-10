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
#' @param group_cols Character vector (used by `pb_groupfilterNA()` only)
#'   specifying sample-annotation column(s) that define the groups within which
#'   missingness must be evaluated.
#' @param min_valid Integer scalar (used by `pb_groupfilterNA()` only),
#'   minimum number of non-missing values required within each group to retain a
#'   feature. Default: `2L`.
#' @param pNA Numeric scalar (used by `pb_groupfilterNA()` only), maximum
#'   proportion of missing values allowed within each group to retain a feature.
#'   Must be in the range `[0, 1]`. If both `min_valid` and `pNA` are provided,
#'   the stricter requirement is applied per group by enforcing the larger
#'   minimum number of observed values implied by either threshold.
#' @param ... Additional parameters forwarded to the underlying
#'   `QFeatures` method where applicable.
#' @return `pb_zeroIsNA()`, `pb_infIsNA()`, `pb_filterNA()` and
#'   `pb_groupfilterNA()` return the updated `ProBatchFeatures` object.
#'   `pb_nNA()` returns the output of the corresponding `QFeatures::nNA()` call
#'   (a `list` of `DataFrame`s).
#' @details For grouped filtering, features are retained if they meet the
#'   missingness criteria in at least one group defined by `group_cols`.
#' @name pb_missing_helpers
NULL

#' @rdname pb_missing_helpers
#' @export
pb_zeroIsNA <- function(object, pbf_name = names(object), ...) {
    stopifnot(is(object, "ProBatchFeatures"))
    assays <- .pb_require_materialised_assays(object, pbf_name)
    params <- .pb_collect_missing_params(list(...), forbidden = "i")
    for (nm in assays) {
        prior <- object
        object <- do.call(zeroIsNA, c(list(object, i = nm), params))
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
    stopifnot(is(object, "ProBatchFeatures"))
    assays <- .pb_require_materialised_assays(object, pbf_name)
    params <- .pb_collect_missing_params(list(...), forbidden = "i")
    for (nm in assays) {
        prior <- object
        object <- do.call(infIsNA, c(list(object, i = nm), params))
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
    stopifnot(is(object, "ProBatchFeatures"))
    assays <- .pb_require_materialised_assays(object, pbf_name)
    params <- .pb_collect_missing_params(list(...), forbidden = "i")
    res <- lapply(
        assays,
        function(nm) do.call(nNA, c(list(object, i = nm), params))
    )
    if (length(res) == 1L) {
        return(res[[1L]])
    }
    res <- setNames(res, assays)
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
  ...
) {
    stopifnot(is(object, "ProBatchFeatures"))
    stopifnot(is.logical(inplace), length(inplace) == 1L)

    assays <- .pb_missing_assays_for_input(
        object = object,
        pbf_name = pbf_name,
        default = "all",
        inform_if_default = TRUE
    )
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
            object <- do.call(filterNA, c(list(object, i = nm), params))
            object <- .as_ProBatchFeatures(object, from = prior)
            to_nm <- nm
            message("  Features after filtering:\t", length(object[[nm]]))
        } else {
            filtered_obj <- do.call(filterNA, c(list(object, i = nm), params))
            filtered_obj <- .as_ProBatchFeatures(filtered_obj, from = object)
            filtered <- filtered_obj[[nm]]
            new_nm <- .pb_unique_assay_name(object, final_name[[idx]])
            prior <- object
            object <- addAssay(object, filtered, name = new_nm)
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

#' @rdname pb_missing_helpers
#' @export
pb_groupfilterNA <- function(
  object,
  pbf_name = NULL,
  group_cols,
  min_valid = 2L,
  pNA = NULL,
  inplace = FALSE,
  final_name = NULL,
  ...
) {
    stopifnot(is(object, "ProBatchFeatures"))
    stopifnot(is.logical(inplace), length(inplace) == 1L)

    if (missing(group_cols) || is.null(group_cols) || !length(group_cols)) {
        stop("`group_cols` must be provided and non-empty.", call. = FALSE)
    }
    group_cols <- as.character(group_cols)
    if (anyNA(group_cols) || !all(nzchar(group_cols))) {
        stop("`group_cols` must contain non-missing, non-empty column names.", call. = FALSE)
    }

    if (!is.null(min_valid)) {
        min_valid <- as.integer(min_valid)
        if (length(min_valid) != 1L || is.na(min_valid) || min_valid < 0L) {
            stop("`min_valid` must be a single non-negative integer.", call. = FALSE)
        }
    }

    if (!is.null(pNA)) {
        if (!is.numeric(pNA) || length(pNA) != 1L || is.na(pNA) || pNA < 0 || pNA > 1) {
            stop("`pNA` must be a single numeric value between 0 and 1.", call. = FALSE)
        }
        pNA <- as.numeric(pNA)
    }

    if (is.null(min_valid) && is.null(pNA)) {
        stop("Provide at least one of `min_valid` or `pNA` to perform filtering.", call. = FALSE)
    }

    assays <- .pb_missing_assays_for_input(
        object = object,
        pbf_name = pbf_name,
        default = "all",
        inform_if_default = TRUE
    )
    params <- .pb_collect_missing_params(list(...), forbidden = c("i", "name", "min", "pNA"))

    if (!inplace) {
        final_name <- .pb_prepare_final_names(assays, final_name, suffix = "_groupfilteredNA")
    } else if (!is.null(final_name)) {
        warning("`final_name` is ignored when `inplace = TRUE`.")
    }

    for (idx in seq_along(assays)) {
        nm <- assays[[idx]]
        message("Processing assay:", nm)
        current <- object[[nm]]
        features_before <- nrow(current)
        message("  Features before filtering:\t", features_before)

        cd <- SummarizedExperiment::colData(current)
        cd_df <- as.data.frame(cd)

        missing_cols <- setdiff(group_cols, names(cd_df))
        if (length(missing_cols)) {
            stop(
                "Assay '", nm, "' is missing group column(s): ",
                paste(missing_cols, collapse = ", "),
                call. = FALSE
            )
        }

        group_df <- cd_df[, group_cols, drop = FALSE]
        has_na_group <- vapply(group_df, function(col) any(is.na(col)), logical(1L))
        if (any(has_na_group)) {
            bad_cols <- paste(group_cols[has_na_group], collapse = ", ")
            stop(
                "Group column(s) ", bad_cols,
                " contain NA values in assay '", nm, "'.",
                call. = FALSE
            )
        }

        group_factor <- interaction(group_df, drop = TRUE, lex.order = TRUE)
        split_indices <- split(seq_along(group_factor), group_factor, drop = TRUE)

        feature_ids <- rownames(current)
        if (is.null(feature_ids)) {
            stop(
                "Assay '", nm, "' has no rownames; cannot perform grouped filtering.",
                call. = FALSE
            )
        }
        keep_logical <- setNames(rep(FALSE, length(feature_ids)), feature_ids)

        if (length(feature_ids)) {
            for (grp_name in names(split_indices)) {
                idx_cols <- split_indices[[grp_name]]
                group_size <- length(idx_cols)
                if (!is.null(min_valid) && group_size < min_valid) {
                    stop(
                        "Assay '", nm, "' has group '", grp_name,
                        "' with ", group_size,
                        " sample(s); requires at least ", min_valid, ".",
                        call. = FALSE
                    )
                }
                sub_se <- current[, idx_cols, drop = FALSE]
                tmp_name <- "tmp_group"
                tmp_obj <- QFeatures(setNames(list(sub_se), tmp_name))
                # Derive the per-group minimum number of observed values implied by
                # `min_valid` and `pNA`, then convert to an allowed missingness proportion.
                inferred_min_valid <- 0L
                if (!is.null(min_valid)) {
                    inferred_min_valid <- max(inferred_min_valid, min_valid)
                }
                if (!is.null(pNA)) {
                    required_from_pna <- as.integer(ceiling((1 - pNA) * group_size))
                    inferred_min_valid <- max(inferred_min_valid, required_from_pna)
                }
                inferred_min_valid <- as.integer(inferred_min_valid)
                if (inferred_min_valid == 0L) {
                    p_na <- 1
                } else {
                    p_na <- 1 - (inferred_min_valid / group_size)
                }
                if (!is.finite(p_na) || p_na < 0) {
                    p_na <- 0
                } else if (p_na > 1) {
                    p_na <- 1
                }
                filtered_tmp <- do.call(
                    filterNA,
                    c(list(tmp_obj, i = tmp_name, pNA = p_na), params)
                )
                keep_group <- rownames(filtered_tmp[[tmp_name]])
                if (length(keep_group)) {
                    common <- intersect(keep_group, names(keep_logical))
                    keep_logical[common] <- TRUE
                }
            }
        }

        keep_features <- names(keep_logical)[keep_logical]
        filtered_se <- current[keep_features, , drop = FALSE]
        features_after <- nrow(filtered_se)

        if (inplace) {
            prior <- object
            object[[nm]] <- filtered_se
            object <- .as_ProBatchFeatures(object, from = prior)
            to_nm <- nm
            message("  Features after filtering:\t", features_after)
        } else {
            new_nm <- .pb_unique_assay_name(object, final_name[[idx]])
            prior <- object
            object <- addAssay(object, filtered_se, name = new_nm)
            object <- .as_ProBatchFeatures(object, from = prior)
            to_nm <- new_nm
            message("  Features after filtering:\t", features_after)
        }

        log_params <- c(
            list(
                group_cols = group_cols,
                min_valid = min_valid,
                pNA = pNA,
                inplace = inplace
            ),
            params
        )

        object <- .pb_add_log_entry(
            object,
            step = "groupfilterNA",
            fun = "pb_groupfilterNA",
            from = nm,
            to = to_nm,
            params = log_params
        )
    }

    object
}

.pb_missing_assays_for_input <- function(object,
                                         pbf_name,
                                         default = c("current", "all"),
                                         inform_if_default = FALSE) {
    default <- match.arg(default)
    using_default <- is.null(pbf_name) || !length(pbf_name)

    assays <- .pb_resolve_assays_for_input(
        object = object,
        pbf_name = pbf_name,
        default = default,
        deduplicate = TRUE,
        empty_message = "No assay names available. Provide `pbf_name` or ensure the object stores assays."
    )

    if (using_default && isTRUE(inform_if_default)) {
        if (identical(default, "all")) {
            message("`pbf_name` not provided, using all assays: ", paste(assays, collapse = ", "))
        } else {
            message("`pbf_name` not provided, using the most recent assay: ", assays[[1]])
        }
    }

    .pb_require_materialised_assays(object, assays)
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
