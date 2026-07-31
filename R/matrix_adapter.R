#' Apply a matrix-oriented method safely
#'
#' Converts supported inputs to a numeric feature-by-sample matrix, aligns
#' sample annotation by identifiers, invokes a function or registered step,
#' validates its result, and reconstructs long data without
#' joins or Cartesian expansion.
#'
#' @param x Numeric feature-by-sample matrix or long data frame.
#'   \code{ProBatchFeatures} objects must be transformed with
#'   \code{pb_transform()} so storage and lineage are preserved.
#' @param fun Function or registered step identifier. The method receives the
#'   numeric matrix as its first argument. Aligned \code{sample_annotation} is
#'   also supplied when the method declares that argument or accepts \code{...}.
#' @param sample_annotation Optional data frame aligned by
#'   \code{sample_id_col}, or by unique row names when that column is absent.
#' @param feature_id_col Feature identifier column for long input.
#' @param sample_id_col Sample identifier column for long input and annotation.
#' @param measure_col Numeric measurement column for long input.
#' @param missing Missing-value policy. \code{"error"} (the default) rejects
#'   missing values, \code{"keep"} passes them through, \code{"drop_features"}
#'   removes incomplete feature rows, and \code{"fill"} replaces missing
#'   values with \code{fill_value}.
#' @param fill_value Finite numeric scalar used only with
#'   \code{missing = "fill"}.
#' @param keep_all Long-output column policy. \code{"default"} preserves the
#'   input columns, \code{"all"} also appends non-conflicting annotation
#'   columns, and \code{"minimal"} retains only the feature, sample,
#'   and measurement columns.
#' @param ... Additional arguments forwarded to \code{fun}. The adapter-owned
#'   \code{data_matrix} and \code{sample_annotation} arguments must not be
#'   supplied through \code{...}.
#'
#' @return A validated numeric matrix for matrix input, or a long data frame
#'   whose retained rows remain in their original relative order. If the method
#'   returns a \code{pb_step_result}, the adapter returns a
#'   \code{pb_step_result} with the validated or restored data
#'   and unchanged artifacts.
#' @export
#' @examples
#' matrix_input <- matrix(
#'     1:6,
#'     nrow = 2,
#'     dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
#' )
#' annotation <- data.frame(
#'     sample = c("s3", "s1", "s2"),
#'     shift = c(3, 1, 2)
#' )
#' pb_apply_matrix_method(
#'     matrix_input,
#'     function(data_matrix, sample_annotation) {
#'         sweep(data_matrix, 2, sample_annotation$shift, "+")
#'     },
#'     sample_annotation = annotation,
#'     sample_id_col = "sample"
#' )
pb_apply_matrix_method <- function(
    x,
    fun,
    sample_annotation = NULL,
    feature_id_col = "peptide_group_label",
    sample_id_col = "FullRunName",
    measure_col = "Intensity",
    missing = "error",
    fill_value = NULL,
    keep_all = "default",
    ...
) {
    if (inherits(x, "ProBatchFeatures")) {
        stop(
            "`ProBatchFeatures` transformations must use `pb_transform()` ",
            "so assay storage, lineage, and operation logging are preserved.",
            call. = FALSE
        )
    }

    feature_id_col <- .pb_adapter_scalar_name(
        feature_id_col,
        "feature_id_col"
    )
    sample_id_col <- .pb_adapter_scalar_name(
        sample_id_col,
        "sample_id_col"
    )
    measure_col <- .pb_adapter_scalar_name(measure_col, "measure_col")
    if (anyDuplicated(c(feature_id_col, sample_id_col, measure_col))) {
        stop(
            "`feature_id_col`, `sample_id_col`, and `measure_col` must name ",
            "distinct columns.",
            call. = FALSE
        )
    }

    canonical_missing <- c("error", "keep", "drop_features", "fill")
    if (!is.character(missing) ||
        length(missing) != 1L ||
        is.na(missing) ||
        !missing %in% canonical_missing) {
        stop(
            "`missing` must be one of: ",
            paste(shQuote(canonical_missing), collapse = ", "),
            ".",
            call. = FALSE
        )
    }
    missing_policy <- .pb_normalize_missing_policy(
        missing,
        fill_value = fill_value,
        argument = "missing"
    )
    keep_all <- match.arg(
        keep_all,
        choices = c("default", "all", "minimal")
    )

    dots <- list(...)
    dot_names <- names(dots)
    if (!is.null(dot_names)) {
        named <- !is.na(dot_names) & nzchar(dot_names)
        reserved <- intersect(
            dot_names[named],
            c("data_matrix", "sample_annotation")
        )
        if (length(reserved)) {
            stop(
                "Do not supply adapter-owned argument(s) through `...`: ",
                paste(unique(reserved), collapse = ", "),
                ".",
                call. = FALSE
            )
        }
    }

    input_kind <- if (is.matrix(x)) {
        "matrix"
    } else if (is.data.frame(x)) {
        "long"
    } else {
        stop(
            "`x` must be a numeric matrix or a long data frame.",
            call. = FALSE
        )
    }

    prepared <- if (identical(input_kind, "matrix")) {
        list(
            matrix = .pb_adapter_validate_matrix(x, "Input matrix"),
            original = NULL
        )
    } else {
        .pb_adapter_long_to_matrix(
            x,
            feature_id_col = feature_id_col,
            sample_id_col = sample_id_col,
            measure_col = measure_col
        )
    }

    working <- .pb_adapter_apply_missing(
        prepared$matrix,
        missing = missing_policy$policy,
        fill_value = missing_policy$fill_value
    )
    annotation <- .pb_adapter_align_annotation(
        sample_annotation,
        sample_ids = colnames(working),
        sample_id_col = sample_id_col
    )

    method <- .pb_adapter_resolve_method(fun)
    method_formals <- tryCatch(
        names(formals(method)),
        error = function(...) character()
    )
    method_args <- dots
    if (any(c("sample_annotation", "...") %in% method_formals)) {
        method_args <- c(
            list(sample_annotation = annotation),
            method_args
        )
    }
    method_result <- .pb_step_result_parts(do.call(
        method,
        c(list(working), method_args)
    ))
    result <- .pb_adapter_validate_output(
        method_result$data,
        input_matrix = working,
        missing = missing_policy$policy
    )

    if (identical(input_kind, "long")) {
        result <- .pb_adapter_restore_long(
            original = prepared$original,
            result = result,
            annotation = annotation,
            feature_id_col = feature_id_col,
            sample_id_col = sample_id_col,
            measure_col = measure_col,
            keep_all = keep_all
        )
    }

    if (method_result$structured) {
        return(pb_step_result(result, method_result$artifacts))
    }
    result
}

.pb_adapter_scalar_name <- function(value, argument) {
    if (!is.character(value) || length(value) != 1L ||
        is.na(value) || !nzchar(value)) {
        stop(
            "`", argument, "` must be a non-empty character scalar.",
            call. = FALSE
        )
    }
    value
}

.pb_adapter_validate_matrix <- function(value, context) {
    if (!is.matrix(value) || !is.numeric(value)) {
        stop(context, " must be a numeric matrix.", call. = FALSE)
    }
    if (!nrow(value) || !ncol(value)) {
        stop(
            context,
            " must have at least one feature and one sample.",
            call. = FALSE
        )
    }

    rownames(value) <- .pb_validate_identifiers(
        rownames(value),
        paste0(context, " feature axis")
    )
    colnames(value) <- .pb_validate_identifiers(
        colnames(value),
        paste0(context, " sample axis")
    )
    value
}

.pb_adapter_long_to_matrix <- function(
    value,
    feature_id_col,
    sample_id_col,
    measure_col
) {
    if (anyDuplicated(names(value))) {
        stop("Long input must have unique column names.", call. = FALSE)
    }
    required <- c(feature_id_col, sample_id_col, measure_col)
    missing_columns <- setdiff(required, names(value))
    if (length(missing_columns)) {
        stop(
            "Long input is missing required column(s): ",
            paste(missing_columns, collapse = ", "),
            ".",
            call. = FALSE
        )
    }
    if (!is.numeric(value[[measure_col]])) {
        stop(
            "Long input measurement column `", measure_col,
            "` must be numeric.",
            call. = FALSE
        )
    }

    .pb_validate_long_keys(
        value,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col
    )
    feature_order <- unique(as.character(value[[feature_id_col]]))
    sample_order <- unique(as.character(value[[sample_id_col]]))
    converted <- long_to_matrix(
        value,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col
    )
    converted <- converted[
        feature_order,
        sample_order,
        drop = FALSE
    ]

    list(matrix = converted, original = value)
}

.pb_adapter_apply_missing <- function(value, missing, fill_value) {
    result <- handle_missing_values(
        value,
        warning_message = paste0(
            "Input contains missing values. Choose an explicit `missing` ",
            "policy."
        ),
        fill_the_missing = missing,
        fill_value = fill_value
    )
    if (!nrow(result)) {
        stop(
            "`missing = \"drop_features\"` removed every feature.",
            call. = FALSE
        )
    }
    result
}

.pb_adapter_align_annotation <- function(
    sample_annotation,
    sample_ids,
    sample_id_col
) {
    if (is.null(sample_annotation)) {
        annotation <- data.frame(
            sample_ids,
            stringsAsFactors = FALSE,
            row.names = sample_ids
        )
        names(annotation) <- sample_id_col
        return(annotation)
    }
    if (!is.data.frame(sample_annotation)) {
        stop(
            "`sample_annotation` must be a data frame or NULL.",
            call. = FALSE
        )
    }
    if (anyDuplicated(names(sample_annotation))) {
        stop(
            "`sample_annotation` must have unique column names.",
            call. = FALSE
        )
    }

    annotation <- as.data.frame(
        sample_annotation,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    annotation_ids <- if (sample_id_col %in% names(annotation)) {
        .pb_validate_identifiers(
            annotation[[sample_id_col]],
            "`sample_annotation` sample identifiers"
        )
    } else {
        .pb_validate_identifiers(
            rownames(annotation),
            "`sample_annotation` row names"
        )
    }
    matches <- match(sample_ids, annotation_ids)
    if (anyNA(matches)) {
        stop(
            "`sample_annotation` is missing sample identifier(s): ",
            paste(sample_ids[is.na(matches)], collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    annotation <- annotation[matches, , drop = FALSE]
    annotation <- .align_sample_annotation(
        annotation,
        sample_ids = sample_ids,
        sample_id_col = sample_id_col
    )
    rownames(annotation) <- sample_ids
    if (!(sample_id_col %in% names(annotation))) {
        annotation[[sample_id_col]] <- sample_ids
    }
    annotation
}

.pb_adapter_resolve_method <- function(fun) {
    if (is.function(fun)) {
        return(fun)
    }

    resolution <- .pb_resolve_step(fun)
    if (!is.list(resolution) || !is.function(resolution$fun)) {
        stop(
            "Internal step resolution did not return a callable function.",
            call. = FALSE
        )
    }
    resolution$fun
}

.pb_adapter_validate_output <- function(
    result,
    input_matrix,
    missing,
    allow_unnamed_features = FALSE
) {
    restore_unnamed_features <- isTRUE(allow_unnamed_features) &&
        is.null(rownames(input_matrix)) &&
        is.null(rownames(result))
    if (restore_unnamed_features &&
        is.matrix(result) &&
        is.numeric(result) &&
        nrow(result) &&
        ncol(result)) {
        if (nrow(result) != nrow(input_matrix)) {
            stop(
                "Method output cannot change the number of features when ",
                "the input feature axis is unnamed.",
                call. = FALSE
            )
        }
        temporary_feature_ids <- paste0(
            ".pb_unnamed_feature_",
            seq_len(nrow(input_matrix))
        )
        rownames(input_matrix) <- temporary_feature_ids
        rownames(result) <- temporary_feature_ids
    }

    result <- .pb_adapter_validate_matrix(result, "Method output")
    unexpected_features <- setdiff(rownames(result), rownames(input_matrix))
    unexpected_samples <- setdiff(colnames(result), colnames(input_matrix))
    if (length(unexpected_features) || length(unexpected_samples)) {
        stop(
            "Method output introduced unknown identifier(s). Features: ",
            paste(unexpected_features, collapse = ", "),
            "; samples: ",
            paste(unexpected_samples, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    expected_feature_order <- rownames(input_matrix)[
        rownames(input_matrix) %in% rownames(result)
    ]
    expected_sample_order <- colnames(input_matrix)[
        colnames(input_matrix) %in% colnames(result)
    ]
    if (!identical(rownames(result), expected_feature_order) ||
        !identical(colnames(result), expected_sample_order)) {
        stop(
            "Method output must preserve input feature and sample order; ",
            "ordered subsetting is allowed.",
            call. = FALSE
        )
    }
    if (!identical(missing, "keep") && anyNA(result)) {
        stop(
            "Method output contains missing values that are incompatible with ",
            "`missing = \"", missing, "\"`.",
            call. = FALSE
        )
    }
    if (restore_unnamed_features) {
        rownames(result) <- NULL
    }
    result
}

.pb_adapter_restore_long <- function(
    original,
    result,
    annotation,
    feature_id_col,
    sample_id_col,
    measure_col,
    keep_all
) {
    feature_index <- match(
        as.character(original[[feature_id_col]]),
        rownames(result)
    )
    sample_index <- match(
        as.character(original[[sample_id_col]]),
        colnames(result)
    )
    keep <- !is.na(feature_index) & !is.na(sample_index)
    restored <- original[keep, , drop = FALSE]
    restored[[measure_col]] <- result[cbind(
        feature_index[keep],
        sample_index[keep]
    )]

    if (identical(keep_all, "all")) {
        annotation_columns <- setdiff(
            names(annotation),
            c(sample_id_col, names(restored))
        )
        annotation_rows <- match(
            as.character(restored[[sample_id_col]]),
            rownames(annotation)
        )
        for (column in annotation_columns) {
            restored[[column]] <- annotation[[column]][annotation_rows]
        }
    } else if (identical(keep_all, "minimal")) {
        restored <- restored[, unique(c(
            feature_id_col,
            sample_id_col,
            measure_col
        )), drop = FALSE]
    }
    restored
}
