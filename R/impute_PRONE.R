#' Impute missing intensities with PRONE (optional)
#'
#' Adapter around `PRONE::impute_se()` that lets proBatch users run PRONE's
#' mixed kNN/MNAR imputation on:
#'   - long data.frames (one feature per row and sample),
#'   - wide data matrices (features in rows, samples in columns).
#'
#' For `ProBatchFeatures` pipelines, use [pb_transform()] with the registered
#' `"PRONEImpute"` step.
#'
#' @param x Object to impute; one of:
#'   \itemize{
#'     \item long data.frame (proBatch long format, see [proBatch()]);
#'     \item numeric matrix/data.frame with features in rows and samples in columns.
#'   }
#' @param sample_annotation Optional sample annotation table; must contain
#'   `sample_id_col` if provided.
#' @param sample_id_col Column in `sample_annotation` that matches matrix column
#'   names or long-format sample IDs. Default `"FullRunName"`.
#' @param feature_id_col Column in long-format data holding peptide/protein IDs.
#'   Default `"feature_id"`.
#' @param measure_col Column with intensities in long-format data.
#' @param condition_col Optional column name in `sample_annotation` that is
#'   passed as the `condition` argument to `PRONE::impute_se()` (controls group
#'   separation before imputation). If `NULL`, PRONE imputes globally.
#'
#' @return Object of the same class as `x`, with imputed values.
#' @rdname imputePRONE
#' @export
imputePRONE <- function(x,
                        sample_annotation = NULL,
                        sample_id_col = "FullRunName",
                        feature_id_col = "feature_id",
                        measure_col = "Intensity",
                        condition_col = NULL) {
    .pb_requireNamespace("PRONE")

    if (is.matrix(x)) {
        return(imputePRONE_dm(
            x = x,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            condition_col = condition_col,
            assay_in = "tmp_raw"
        ))
    }

    if (is.data.frame(x)) {
        if (.pb_is_long_df_input(x, feature_id_col, sample_id_col, measure_col)) {
            return(imputePRONE_df(
                x = x,
                sample_annotation = sample_annotation,
                sample_id_col = sample_id_col,
                feature_id_col = feature_id_col,
                measure_col = measure_col,
                condition_col = condition_col,
                assay_in = "tmp_raw"
            ))
        }
        return(imputePRONE_dm(
            x = x,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            condition_col = condition_col,
            assay_in = "tmp_raw"
        ))
    }

    if (inherits(x, "ProBatchFeatures") || inherits(x, "QFeatures")) {
        stop(
            "imputePRONE(): ProBatchFeatures/QFeatures inputs are not supported. ",
            "Use pb_transform(..., steps = 'PRONEImpute') to apply PRONE within a pipeline."
        )
    }

    stop("imputePRONE(): unsupported object of class ", paste(class(x), collapse = "/"))
}

# ---- long-format front-end -------------------------------------------------

#' @rdname imputePRONE
#' @export
imputePRONE_df <- function(x,
                           sample_annotation = NULL,
                           sample_id_col = "FullRunName",
                           feature_id_col = "feature_id",
                           measure_col = "Intensity",
                           condition_col = NULL,
                           assay_in = "raw") {
    .pb_requireNamespace("PRONE")

    .pb_transform_long_via_matrix(
        df_long = x,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col,
        matrix_fun = function(data_matrix) {
            .prone_matrix_step(
                data_matrix = data_matrix,
                sample_annotation = sample_annotation,
                sample_id_col = sample_id_col,
                condition_col = condition_col,
                assay_in = assay_in
            )
        },
        sample_annotation = sample_annotation
    )
}

# ---- wide-format front-end -------------------------------------------------

#' @rdname imputePRONE
#' @export
imputePRONE_dm <- function(x,
                           sample_annotation = NULL,
                           sample_id_col = "FullRunName",
                           condition_col = NULL,
                           assay_in = "raw") {
    .pb_requireNamespace("PRONE")

    data_matrix <- if (is.matrix(x)) x else as.matrix(x)

    .prone_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        condition_col = condition_col,
        assay_in = assay_in
    )
}

# ---- shared helper --------------------------------------------------------

.pb_prone_prepare_sample_annotation <- function(sample_annotation,
                                                sample_ids,
                                                sample_id_col = "FullRunName",
                                                align = TRUE,
                                                placeholder_col = ".pb_prone_sample_id") {
    sample_annotation_local <- sample_annotation
    sample_id_col_local <- sample_id_col
    placeholder_col_local <- NULL

    if (is.null(sample_annotation_local)) {
        placeholder_col_local <- placeholder_col
        sample_annotation_local <- data.frame(sample_ids, stringsAsFactors = FALSE)
        colnames(sample_annotation_local) <- placeholder_col_local
        sample_id_col_local <- placeholder_col_local
    }

    if (isTRUE(align)) {
        sample_annotation_local <- .align_sample_annotation(
            sample_annotation = sample_annotation_local,
            sample_ids = sample_ids,
            sample_id_col = sample_id_col_local
        )
    }

    list(
        sample_annotation = sample_annotation_local,
        sample_id_col = sample_id_col_local,
        placeholder_col = placeholder_col_local
    )
}

.pb_prone_make_se <- function(data_matrix,
                              sample_annotation,
                              assay_in = "raw",
                              force_column = FALSE) {
    sample_df <- as.data.frame(sample_annotation, stringsAsFactors = FALSE)
    sample_ids <- colnames(data_matrix)
    feature_ids <- rownames(data_matrix)

    if (is.null(feature_ids)) {
        feature_ids <- sprintf("feature_%d", seq_len(nrow(data_matrix)))
    } else {
        feature_ids <- as.character(feature_ids)
        bad_ids <- is.na(feature_ids) | !nzchar(feature_ids)
        if (any(bad_ids)) {
            feature_ids[bad_ids] <- sprintf("feature_%d", which(bad_ids))
        }
        feature_ids <- make.unique(feature_ids, sep = "_dup")
    }

    rownames(sample_df) <- sample_ids
    col_data <- S4Vectors::DataFrame(sample_df)
    rownames(col_data) <- sample_ids
    if (isTRUE(force_column) || !("Column" %in% colnames(col_data))) {
        col_data$Column <- rownames(col_data)
    }
    row_data <- S4Vectors::DataFrame(
        "Protein.IDs" = feature_ids,
        "IDs" = feature_ids,
        check.names = FALSE
    )

    SummarizedExperiment::SummarizedExperiment(
        assays = setNames(list(data_matrix), assay_in),
        colData = col_data,
        rowData = row_data
    )
}

.pb_prone_prepare_condition <- function(sample_df, sample_ids, condition_col = NULL) {
    condition_arg <- NULL
    if (is.null(condition_col)) {
        return(list(sample_df = sample_df, condition_arg = condition_arg))
    }

    if (is.character(condition_col) && length(condition_col) == 1L) {
        if (!(condition_col %in% colnames(sample_df))) {
            stop(
                "PRONE imputation: condition column '", condition_col,
                "' not found in sample annotation.",
                call. = FALSE
            )
        }
        condition_arg <- condition_col
        return(list(sample_df = sample_df, condition_arg = condition_arg))
    }

    cond_values <- condition_col
    if (!is.null(names(cond_values))) {
        idx <- match(sample_ids, names(cond_values))
        if (anyNA(idx)) {
            stop(
                "PRONE imputation: condition vector is missing values for some samples.",
                call. = FALSE
            )
        }
        cond_values <- cond_values[idx]
    } else if (length(cond_values) != length(sample_ids)) {
        stop(
            "PRONE imputation: unnamed condition vector must match the number of samples.",
            call. = FALSE
        )
    }

    cond_values <- unname(as.vector(cond_values))

    existing_names <- colnames(sample_df)
    base_name <- ".pb_prone_condition"
    new_name <- base_name
    counter <- 1L
    while (!is.null(existing_names) && new_name %in% existing_names) {
        counter <- counter + 1L
        new_name <- paste0(base_name, "_", counter)
    }

    sample_df[[new_name]] <- cond_values
    condition_arg <- new_name

    list(sample_df = sample_df, condition_arg = condition_arg)
}

.prone_matrix_step <- function(data_matrix,
                               sample_annotation = NULL,
                               sample_id_col = "FullRunName",
                               condition_col = NULL,
                               assay_in = "raw") {
    .pb_requireNamespace("PRONE")

    if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
    }
    storage.mode(data_matrix) <- "double"

    if (is.null(colnames(data_matrix))) {
        stop("PRONE imputation requires matrix column names (sample identifiers).", call. = FALSE)
    }
    sample_ids <- colnames(data_matrix)
    annotation_prep <- .pb_prone_prepare_sample_annotation(
        sample_annotation = sample_annotation,
        sample_ids = sample_ids,
        sample_id_col = sample_id_col,
        align = FALSE
    )
    sample_annotation_local <- annotation_prep$sample_annotation
    sample_id_col_local <- annotation_prep$sample_id_col
    placeholder_col <- annotation_prep$placeholder_col

    .run_matrix_method(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation_local,
        sample_id_col = sample_id_col_local,
        fill_the_missing = NULL,
        missing_warning = "PRONE imputation cannot operate without sample identifiers.",
        method_fun = function(data_matrix, sample_annotation) {
            original_data_matrix <- data_matrix
            sample_ids_local <- colnames(data_matrix)
            sample_df <- as.data.frame(sample_annotation, stringsAsFactors = FALSE)

            if (!is.null(placeholder_col) && placeholder_col %in% names(sample_df)) {
                sample_df[[placeholder_col]] <- NULL
            }

            condition_prep <- .pb_prone_prepare_condition(
                sample_df = sample_df,
                sample_ids = sample_ids_local,
                condition_col = condition_col
            )
            sample_df <- condition_prep$sample_df
            condition_arg <- condition_prep$condition_arg

            se <- .pb_prone_make_se(
                data_matrix = data_matrix,
                sample_annotation = sample_df,
                assay_in = assay_in,
                force_column = TRUE
            )

            prone_impute <- .pb_prone_impute_fun()
            # TODO: switch back to PRONE::impute_se once the upstream bug is fixed.
            se_imp <- prone_impute(
                se,
                ain = assay_in,
                condition = condition_arg
            )

            # Prefer the override's "<assay>_imputed" output when available; fall back
            # to the original assay to stay compatible with upstream PRONE builds.
            imputed_assay <- paste0(assay_in, "_imputed")
            assay_names <- SummarizedExperiment::assayNames(se_imp)
            assay_to_use <- if (imputed_assay %in% assay_names) imputed_assay else assay_in

            imputed_mat <- SummarizedExperiment::assay(se_imp, assay_to_use)
            if (!is.matrix(imputed_mat)) {
                imputed_mat <- as.matrix(imputed_mat)
            }
            storage.mode(imputed_mat) <- "double"
            imputed_mat <- .pb_prone_restore_dimnames(
                original_matrix = original_data_matrix,
                imputed_matrix = imputed_mat
            )

            if (!is.null(condition_arg)) {
                imputed_mat <- .pb_prone_condition_cleanup(
                    original_matrix = original_data_matrix,
                    imputed_matrix = imputed_mat
                )
            }

            imputed_mat
        }
    )
}

.pb_prone_restore_dimnames <- function(original_matrix, imputed_matrix) {
    out <- imputed_matrix

    expected_cols <- colnames(original_matrix)
    if (!is.null(expected_cols)) {
        current_cols <- colnames(out)
        has_bad_cols <- is.null(current_cols) ||
            length(current_cols) != ncol(out) ||
            anyNA(current_cols) ||
            any(!nzchar(current_cols))

        if (has_bad_cols) {
            if (ncol(out) != length(expected_cols)) {
                stop(
                    "PRONE imputation returned a matrix with unexpected sample dimensions.",
                    call. = FALSE
                )
            }
            colnames(out) <- expected_cols
        } else if (length(expected_cols) == ncol(out) &&
            !anyDuplicated(expected_cols) &&
            !anyDuplicated(current_cols) &&
            setequal(expected_cols, current_cols)) {
            out <- out[, expected_cols, drop = FALSE]
        } else if (ncol(out) == length(expected_cols)) {
            colnames(out) <- expected_cols
        } else {
            stop(
                "PRONE imputation returned columns that cannot be aligned to the original samples.",
                call. = FALSE
            )
        }
    }

    expected_rows <- rownames(original_matrix)
    if (!is.null(expected_rows)) {
        current_rows <- rownames(out)
        has_bad_rows <- is.null(current_rows) ||
            length(current_rows) != nrow(out) ||
            anyNA(current_rows) ||
            any(!nzchar(current_rows))

        if (!has_bad_rows && !anyDuplicated(current_rows) && !anyDuplicated(expected_rows)) {
            if (nrow(out) != length(expected_rows) || !setequal(expected_rows, current_rows)) {
                message(
                    "PRONE imputation returned a different feature set than input; ",
                    "missing features will be padded with NA."
                )
            }
            aligned <- matrix(
                NA_real_,
                nrow = length(expected_rows),
                ncol = ncol(out),
                dimnames = list(expected_rows, colnames(out))
            )
            idx_expected <- match(current_rows, expected_rows)
            keep <- !is.na(idx_expected)
            if (any(keep)) {
                aligned[idx_expected[keep], ] <- out[keep, , drop = FALSE]
            }
            out <- aligned
        } else if (has_bad_rows && nrow(out) == length(expected_rows)) {
            rownames(out) <- expected_rows
        } else if (has_bad_rows && nrow(out) < length(expected_rows)) {
            message(
                "PRONE imputation returned fewer features without valid row names; ",
                "padding missing features with NA by position."
            )
            aligned <- matrix(
                NA_real_,
                nrow = length(expected_rows),
                ncol = ncol(out),
                dimnames = list(expected_rows, colnames(out))
            )
            if (nrow(out) > 0L) {
                aligned[seq_len(nrow(out)), ] <- out
            }
            out <- aligned
        } else {
            stop(
                "PRONE imputation returned rows that cannot be aligned to the original features.",
                call. = FALSE
            )
        }
    }

    out
}

.pb_prone_condition_cleanup <- function(original_matrix, imputed_matrix) {
    rows_with_na_before <- rowSums(is.na(original_matrix)) > 0L
    rows_with_na_after <- rowSums(is.na(imputed_matrix)) > 0L

    before_count <- sum(rows_with_na_before)
    after_count <- sum(rows_with_na_after)

    base_msg <- paste0(
        "PRONE condition imputation: ",
        before_count,
        " rows had missing values before; ",
        after_count,
        " remain after."
    )

    if (after_count > 0L) {
        message(base_msg, " Keeping all rows; unresolved values remain NA.")
        return(imputed_matrix)
    }

    message(base_msg, " No additional NA rows remained.")
    imputed_matrix
}

.pb_prone_impute_fun <- local({
    override_fun <- NULL
    function() {
        opt_fun <- getOption("proBatch.prone_impute_se", NULL)
        if (is.function(opt_fun)) {
            return(opt_fun)
        }
        if (is.null(override_fun)) {
            override_fun <<- .pb_load_prone_impute()
        }
        override_fun
    }
})

.pb_load_prone_impute <- function() {
    override_path <- .pb_prone_override_path()
    if (!is.null(override_path)) {
        override_env <- new.env(parent = environment())
        tryCatch(
            {
                sys.source(override_path, envir = override_env)
                fun <- override_env$impute_se
                if (is.function(fun)) {
                    return(fun)
                }
            },
            error = function(e) NULL
        )
    }
    PRONE::impute_se
}

.pb_prone_override_path <- function() {
    inst_path <- system.file("overrides/PRONE/Imputation.R", package = "proBatch")
    if (nzchar(inst_path) && file.exists(inst_path)) {
        return(inst_path)
    }
    pkg_root <- tryCatch(find.package("proBatch"), error = function(e) "")
    if (nzchar(pkg_root)) {
        dev_path <- file.path(pkg_root, "inst", "overrides", "PRONE", "Imputation.R")
        if (file.exists(dev_path)) {
            return(dev_path)
        }
    }
    NULL
}
