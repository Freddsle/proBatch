####################################################
# filters for missing values in ProBatchFeatures
####################################################
#' Apply `QFeatures` missing-data helpers to stored assays
#'
#' These wrappers delegate to the corresponding `QFeatures` generics while
#' ensuring that the requested assays remain part of the `ProBatchFeatures`
#' object. Only assays that are already materialised can be modified. If a
#' transformation step was applied as a "fast" step (log, log2, etc.), consider
#' re-running it with `store_fast_steps = TRUE`.
#'
#' @param object A `ProBatchFeatures` object.
#' @param pbf_name Character vector of assay names. Defaults to
#'   `names(object)` - all assays.
#' @param inplace Logical (used by `pb_filterNA()` and
#'   `pb_groupfilterNA()`), whether to modify the object in place. Default:
#'   `FALSE`. If `FALSE`, the modified assay(s) are added to the object with
#'   `final_name` (if provided) or a function-specific suffix.
#' @param final_name Character (used by `pb_filterNA()` and
#'   `pb_groupfilterNA()`), name for the modified assay(s) if `inplace` is
#'   `FALSE`. An explicit value must provide one unique, non-missing, non-empty
#'   name per selected assay. A collision with a stored assay or virtual
#'   operation-log target is accepted only for an exact idempotent retry with
#'   identical data, parent assay, and stable lineage origin. Explicit names
#'   are never recycled or disambiguated. If `NULL` (default), the original
#'   name(s) receive suffix `_filteredNA` or `_groupfilteredNA`, respectively;
#'   generated collisions across the stored and virtual result namespace are
#'   disambiguated with numeric suffixes.
#' @param group_cols Character vector (used by `pb_groupfilterNA()` only)
#'   specifying sample-annotation column(s) that define the groups within which
#'   missingness must be evaluated.
#' @param min_valid Integer scalar (used by `pb_groupfilterNA()` only),
#'   minimum number of non-missing values required within each group to retain
#'   a feature. Default: `2L`.
#' @param pNA Numeric scalar (used by `pb_groupfilterNA()` only), maximum
#'   proportion of missing values allowed within each group to retain a
#'   feature. Must be in the range `[0, 1]`. If both `min_valid` and `pNA` are
#'   provided, the stricter requirement is applied per group by enforcing the
#'   larger minimum number of observed values implied by either threshold.
#' @param mask_failing Logical scalar (used by `pb_groupfilterNA()` only).
#'   When `TRUE` (default), values in groups that did not pass the missingness
#'   threshold are set to `NA`. When `FALSE`, the original values in failing
#'   groups are kept as-is (the feature is still retained if it
#'   passes in at least one group).
#' @param ... Additional parameters forwarded to the underlying
#'   `QFeatures` method where applicable.
#' @return `pb_zeroIsNA()`, `pb_infIsNA()`, `pb_filterNA()` and
#'   `pb_groupfilterNA()` return the updated `ProBatchFeatures` object. For one
#'   assay, `pb_nNA()` returns the corresponding `QFeatures::nNA()` result (a
#'   list of `DataFrame`s). For multiple assays it returns a named outer list
#'   with one such result per assay and a combined `nNA` `DataFrame`.
#' @details For grouped filtering, features are retained if they meet the
#'   missingness criteria in at least one group defined by `group_cols`. When
#'   `mask_failing = TRUE` (the default), the values in groups where the
#'   feature did not pass the threshold are replaced with `NA`.
#'
#'   With in-place `pb_filterNA()`, QFeatures preserves valid assay links and
#'   prunes their hits to the retained features. In-place `pb_groupfilterNA()`
#'   replaces the assay directly and may remove affected links when feature
#'   rows change; QFeatures warns when this occurs. Use `inplace = FALSE` to
#'   keep existing links unchanged.
#' @examples
#' values <- matrix(
#'     c(0, 1, 2, Inf, 4, 5, 6, 7, 8),
#'     nrow = 3,
#'     byrow = TRUE,
#'     dimnames = list(
#'         paste0("feature", 1:3),
#'         paste0("sample", 1:3)
#'     )
#' )
#' annotation <- data.frame(sample = colnames(values))
#' pbf <- ProBatchFeatures(
#'     values,
#'     sample_annotation = annotation,
#'     sample_id_col = "sample"
#' )
#' pbf <- pb_zeroIsNA(pbf, "feature::raw")
#' pbf <- pb_infIsNA(pbf, "feature::raw")
#' pb_nNA(pbf, "feature::raw")$nNA
#'
#' filtered <- pb_filterNA(
#'     pbf,
#'     pbf_name = "feature::raw",
#'     pNA = 0
#' )
#' rownames(filtered[["feature::raw_filteredNA"]])
#' @name pb_missing_helpers
#' @md
NULL

.pb_apply_missing_qf_step <- function(
    object,
    assays,
    fun,
    step,
    params,
    fun_name = step
) {
    for (nm in assays) {
        prior <- object
        object <- do.call(fun, c(list(object, i = nm), params))
        object <- .as_ProBatchFeatures(object, from = prior)
        object <- .pb_add_log_entry(
            object,
            step = step,
            fun = fun_name,
            from = nm,
            to = nm,
            params = params
        )
    }
    object
}

#' @rdname pb_missing_helpers
#' @export
pb_zeroIsNA <- function(object, pbf_name = names(object), ...) {
    stopifnot(is(object, "ProBatchFeatures"))
    assays <- .pb_require_materialised_assays(object, pbf_name)
    params <- .pb_collect_missing_params(list(...), forbidden = "i")
    ##### changed: deduplicated with .pb_apply_missing_qf_step()
    .pb_apply_missing_qf_step(
        object = object,
        assays = assays,
        fun = zeroIsNA,
        step = "zeroIsNA",
        params = params
    )
}

#' @rdname pb_missing_helpers
#' @export
pb_infIsNA <- function(object, pbf_name = names(object), ...) {
    stopifnot(is(object, "ProBatchFeatures"))
    assays <- .pb_require_materialised_assays(object, pbf_name)
    params <- .pb_collect_missing_params(list(...), forbidden = "i")
    .pb_apply_missing_qf_step(
        object = object,
        assays = assays,
        fun = infIsNA,
        step = "infIsNA",
        params = params
    )
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
    # join "nNA" from each assay into a single DataFrame
    # and add it to the result
    # only for "nNA" and add it as a last result entry
    res$nNA <- do.call(rbind, lapply(res, `[[`, "nNA"))
    res
}

.pb_filterNA_matrix <- function(data_matrix, ...) {
    input <- SummarizedExperiment::SummarizedExperiment(
        assays = list(intensity = data_matrix)
    )
    filtered <- do.call(
        filterNA,
        c(list(input), list(...))
    )
    SummarizedExperiment::assay(filtered, "intensity")
}

.pb_groupfilterNA_matrix <- function(
    data_matrix,
    sample_annotation,
    group_cols,
    min_valid = 2L,
    pNA = NULL,
    mask_failing = TRUE,
    ...
) {
    sample_annotation <- as.data.frame(sample_annotation)
    temporary <- suppressMessages(ProBatchFeatures(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = NULL,
        name = "feature::.pb_groupfilterNA_input"
    ))
    source_name <- pb_current_assay(temporary)
    result_name <- "feature::.pb_groupfilterNA_result"
    filtered <- suppressMessages(pb_groupfilterNA(
        temporary,
        pbf_name = source_name,
        group_cols = group_cols,
        min_valid = min_valid,
        pNA = pNA,
        inplace = FALSE,
        final_name = result_name,
        mask_failing = mask_failing,
        ...
    ))
    SummarizedExperiment::assay(filtered[[result_name]], "intensity")
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

    assays <- .pb_require_materialised_assays(
        object = object,
        assays = .pb_resolve_assays_for_input(
            object = object,
            pbf_name = pbf_name,
            default = "all",
            deduplicate = TRUE,
            inform_if_default = TRUE,
            empty_message = "No assay names available. Provide `pbf_name` or ensure the object stores assays."
        )
    )
    params <- .pb_collect_missing_params(list(...), forbidden = c("i", "name"))

    if (!inplace) {
        final_name <- .pb_prepare_final_names(
            object,
            assays,
            final_name,
            suffix = "_filteredNA"
        )
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
            new_nm <- final_name[[idx]]
            log_params <- c(params, list(inplace = inplace))
            retry_status <- .pb_target_retry_status(
                object = object,
                to = new_nm,
                from = nm,
                data = SummarizedExperiment::assay(filtered, "intensity"),
                step = "filterNA",
                fun = "filterNA",
                params = log_params
            )
            if (!identical(retry_status, "stored_idempotent")) {
                prior <- object
                object <- addAssay(object, filtered, name = new_nm)
                object <- .as_ProBatchFeatures(object, from = prior)
            }
            to_nm <- new_nm
            message("  Features after filtering:\t", length(object[[new_nm]]))
        }
        log_params <- if (inplace) {
            c(params, list(inplace = inplace))
        } else {
            log_params
        }
        object <- .pb_add_log_entry(
            object,
            step = "filterNA",
            fun = "filterNA",
            from = nm,
            to = to_nm,
            params = log_params
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
    mask_failing = TRUE,
    ...
) {
    stopifnot(is(object, "ProBatchFeatures"))
    stopifnot(
        is.logical(inplace),
        length(inplace) == 1L,
        is.logical(mask_failing),
        length(mask_failing) == 1L
    )

    if (missing(group_cols) || is.null(group_cols) || !length(group_cols)) {
        stop("`group_cols` must be provided and non-empty.", call. = FALSE)
    }
    group_cols <- as.character(group_cols)
    if (anyNA(group_cols) || !all(nzchar(group_cols))) {
        stop(
            "`group_cols` must contain non-missing, non-empty column names.",
            call. = FALSE
        )
    }

    if (!is.null(min_valid)) {
        min_valid <- as.integer(min_valid)
        if (length(min_valid) != 1L || is.na(min_valid) || min_valid < 0L) {
            stop(
                "`min_valid` must be a single non-negative integer.",
                call. = FALSE
            )
        }
    }

    if (!is.null(pNA)) {
        if (
            !is.numeric(pNA) ||
                length(pNA) != 1L ||
                is.na(pNA) ||
                pNA < 0 ||
                pNA > 1
        ) {
            stop(
                "`pNA` must be a single numeric value between 0 and 1.",
                call. = FALSE
            )
        }
        pNA <- as.numeric(pNA)
    }

    if (is.null(min_valid) && is.null(pNA)) {
        stop(
            "Provide at least one of `min_valid` or `pNA` to perform filtering.",
            call. = FALSE
        )
    }

    assays <- .pb_require_materialised_assays(
        object = object,
        assays = .pb_resolve_assays_for_input(
            object = object,
            pbf_name = pbf_name,
            default = "all",
            deduplicate = TRUE,
            inform_if_default = TRUE,
            empty_message = "No assay names available. Provide `pbf_name` or ensure the object stores assays."
        )
    )
    params <- .pb_collect_missing_params(
        list(...),
        forbidden = c("i", "name", "min", "pNA")
    )

    if (!inplace) {
        final_name <- .pb_prepare_final_names(
            object,
            assays,
            final_name,
            suffix = "_groupfilteredNA"
        )
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
                "Assay '",
                nm,
                "' is missing group column(s): ",
                paste(missing_cols, collapse = ", "),
                call. = FALSE
            )
        }

        group_df <- cd_df[, group_cols, drop = FALSE]
        has_na_group <- vapply(
            group_df,
            function(col) any(is.na(col)),
            logical(1L)
        )
        if (any(has_na_group)) {
            bad_cols <- paste(group_cols[has_na_group], collapse = ", ")
            stop(
                "Group column(s) ",
                bad_cols,
                " contain NA values in assay '",
                nm,
                "'.",
                call. = FALSE
            )
        }

        group_factor <- interaction(group_df, drop = TRUE, lex.order = TRUE)
        split_indices <- split(
            seq_along(group_factor),
            group_factor,
            drop = TRUE
        )

        feature_ids <- rownames(current)
        if (is.null(feature_ids)) {
            stop(
                "Assay '",
                nm,
                "' has no rownames; cannot perform grouped filtering.",
                call. = FALSE
            )
        }
        keep_logical <- setNames(rep(FALSE, length(feature_ids)), feature_ids)
        group_pass <- list()

        if (length(feature_ids)) {
            for (grp_name in names(split_indices)) {
                idx_cols <- split_indices[[grp_name]]
                group_size <- length(idx_cols)
                if (!is.null(min_valid) && group_size < min_valid) {
                    stop(
                        "Assay '",
                        nm,
                        "' has group '",
                        grp_name,
                        "' with ",
                        group_size,
                        " sample(s); requires at least ",
                        min_valid,
                        ".",
                        call. = FALSE
                    )
                }
                sub_se <- current[, idx_cols, drop = FALSE]
                tmp_name <- "tmp_group"
                tmp_obj <- QFeatures(setNames(list(sub_se), tmp_name))
                # Derive the per-group minimum number of
                # observed values implied by
                # `min_valid` and `pNA`, then convert to an
                # allowed missingness proportion.
                inferred_min_valid <- 0L
                if (!is.null(min_valid)) {
                    inferred_min_valid <- max(inferred_min_valid, min_valid)
                }
                if (!is.null(pNA)) {
                    required_from_pna <-
                        as.integer(ceiling((1 - pNA) * group_size))
                    inferred_min_valid <-
                        max(inferred_min_valid, required_from_pna)
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
                group_pass[[grp_name]] <- keep_group
                if (length(keep_group)) {
                    common <- intersect(keep_group, names(keep_logical))
                    keep_logical[common] <- TRUE
                }
            }
        }

        keep_features <- names(keep_logical)[keep_logical]
        filtered_se <- current[keep_features, , drop = FALSE]
        if (mask_failing && length(keep_features)) {
            mat <- SummarizedExperiment::assay(filtered_se)
            for (grp_name in names(split_indices)) {
                idx_cols <- split_indices[[grp_name]]
                failed_features <- setdiff(
                    keep_features,
                    group_pass[[grp_name]]
                )
                if (length(failed_features)) {
                    mat[failed_features, idx_cols] <- NA
                }
            }
            SummarizedExperiment::assay(filtered_se) <- mat
        }
        features_after <- nrow(filtered_se)
        log_params <- c(
            list(
                group_cols = group_cols,
                min_valid = min_valid,
                pNA = pNA,
                mask_failing = mask_failing,
                inplace = inplace
            ),
            params
        )

        if (inplace) {
            prior <- object
            object[[nm]] <- filtered_se
            object <- .as_ProBatchFeatures(object, from = prior)
            to_nm <- nm
            message("  Features after filtering:\t", features_after)
        } else {
            new_nm <- final_name[[idx]]
            retry_status <- .pb_target_retry_status(
                object = object,
                to = new_nm,
                from = nm,
                data = SummarizedExperiment::assay(filtered_se, "intensity"),
                step = "groupfilterNA",
                fun = "pb_groupfilterNA",
                params = log_params
            )
            if (!identical(retry_status, "stored_idempotent")) {
                prior <- object
                object <- addAssay(object, filtered_se, name = new_nm)
                object <- .as_ProBatchFeatures(object, from = prior)
            }
            to_nm <- new_nm
            message("  Features after filtering:\t", features_after)
        }

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
            "Assay(s) ",
            paste(missing, collapse = ", "),
            " are not stored in the object.",
            extra_msg,
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
                "Argument(s) ",
                paste(bad, collapse = ", "),
                " must not be supplied via `...`.",
                call. = FALSE
            )
        }
    }
    params
}

.pb_result_namespace <- function(object) {
    log <- get_operation_log(object)
    virtual_targets <- if (nrow(log)) as.character(log$to) else character()
    result_names <- unique(c(names(object), virtual_targets))
    result_names[!is.na(result_names) & nzchar(result_names)]
}

.pb_disambiguate_generated_names <- function(proposed, reserved) {
    result <- character(length(proposed))
    for (idx in seq_along(proposed)) {
        candidate <- proposed[[idx]]
        suffix <- 0L
        while (candidate %in% c(reserved, result)) {
            suffix <- suffix + 1L
            candidate <- paste0(proposed[[idx]], ".", suffix)
        }
        result[[idx]] <- candidate
    }
    result
}

.pb_prepare_final_names <- function(object, assays, final_name, suffix) {
    if (is.null(final_name)) {
        return(.pb_disambiguate_generated_names(
            proposed = paste0(assays, suffix),
            reserved = .pb_result_namespace(object)
        ))
    }

    if (!is.character(final_name)) {
        stop("`final_name` must be a character vector.", call. = FALSE)
    }
    if (length(final_name) != length(assays)) {
        stop(
            "`final_name` must contain exactly one name per selected assay ",
            "(expected ",
            length(assays),
            ", got ",
            length(final_name),
            "); ",
            "scalar names are not recycled.",
            call. = FALSE
        )
    }
    if (anyNA(final_name) || any(!nzchar(final_name))) {
        stop(
            "`final_name` must contain only non-missing, non-empty names.",
            call. = FALSE
        )
    }
    if (anyDuplicated(final_name)) {
        duplicated_names <- unique(final_name[duplicated(final_name)])
        stop(
            "`final_name` entries must be unique; duplicated name(s): ",
            paste(duplicated_names, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    final_name
}
