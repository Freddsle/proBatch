#' Data normalization methods
#'
#' @description Normalization of raw (usually log-transformed) data.
#' Normalization brings the samples to the same scale.
#' Currently the following normalization functions are implemented:
#' #' \enumerate{
#'   \item Quantile normalization: `quantile_normalize_dm()`.
#'   Quantile normalization of the data.
#'   \item Median normalization: `normalize_sample_medians_dm()`.
#'   Normalization by centering sample medians to global median of the data
#' }
#' Alternatively, one can call normalization function with `normalize_data_dm()`
#' wrapper.
#'
#'
#' @inheritParams proBatch
#' @param normalize_func global batch normalization method
#' (`quantile` or `MedianCentering`)
#' @param log_base whether to log transform data matrix
#' before normalization (e.g. `NULL`, `2` or `10`)
#' @param offset small positive number to prevent 0 conversion to \code{-Inf}
#' @param sample_annotation optional data frame with sample-level metadata
#'   used to support within-batch normalization for data frames or matrices.
#' @param group_col character(1); column in `sample_annotation` (or `df_long`)
#'   identifying batches/centers for within-batch normalization.
#' @param inside_batch logical; when `TRUE`, perform median centering within
#'   each level of `group_col` after excluding rows that are entirely `NA`
#'   within the corresponding subset.
#'
#' @return the data in the same format as input (\code{data_matrix} or
#' \code{df_long}).
#' For \code{df_long} the data frame stores the original values of
#' \code{measure_col}
#' in another column called "preNorm_intensity" if "intensity", and the
#' normalized values in \code{measure_col} column.
#'
#' @examples
#' data(list = c("example_proteome", "example_proteome_matrix"), package = "proBatch")
#'
#' # Quantile normalization:
#' quantile_normalized_matrix <- quantile_normalize_dm(example_proteome_matrix)
#'
#' # Median centering:
#' median_normalized_df <- normalize_sample_medians_df(example_proteome)
#'
#' # Transform the data in one go:
#' quantile_normalized_matrix <- normalize_data_dm(example_proteome_matrix,
#'     normalize_func = "quantile", log_base = 2, offset = 1
#' )
#'
#' @name normalize
NULL

#'
#' @export
#' @rdname normalize
#'
quantile_normalize_dm <- function(data_matrix) {
    q_norm_proteome <- normalize.quantiles(data_matrix)
    colnames(q_norm_proteome) <- colnames(data_matrix)
    rownames(q_norm_proteome) <- rownames(data_matrix)
    return(q_norm_proteome)
}

#'
#' @export
#' @rdname normalize
#'
quantile_normalize_df <- function(df_long,
                                  feature_id_col = "peptide_group_label",
                                  sample_id_col = "FullRunName",
                                  measure_col = "Intensity",
                                  no_fit_imputed = TRUE,
                                  qual_col = NULL,
                                  qual_value = 2,
                                  keep_all = "default") {
    if (is.null(qual_col) & no_fit_imputed) {
        warning("imputed value flag column is NULL, changing no_fit_imputed to FALSE")
        no_fit_imputed <- FALSE
    }

    if (no_fit_imputed) {
        if (!(qual_col %in% names(df_long))) {
            stop("imputed value flag column (qual_col) is not in the data frame!")
        }
        message("removing imputed values (requants) from the matrix")
        data_matrix <- long_to_matrix(
            df_long,
            feature_id_col = feature_id_col,
            measure_col = measure_col,
            sample_id_col = sample_id_col,
            qual_col = qual_col,
            qual_value = qual_value
        )
    } else {
        if (!is.null(qual_col) && (qual_col %in% names(df_long))) {
            warning("imputed value (requant) column is in the data, are you sure you
                want to use imputed (requant) values in quantile inference?")
        }
        data_matrix <- long_to_matrix(
            df_long,
            feature_id_col = feature_id_col,
            measure_col = measure_col,
            sample_id_col = sample_id_col,
            qual_col = NULL
        )
    }

    q_norm_proteome <- normalize.quantiles(data_matrix)
    colnames(q_norm_proteome) <- colnames(data_matrix)
    rownames(q_norm_proteome) <- rownames(data_matrix)
    normalized_df <- matrix_to_long(
        q_norm_proteome,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col
    )

    old_measure_col <- paste("preNorm", measure_col, sep = "_")

    df_long <- df_long %>%
        rename(!!(old_measure_col) := !!(sym(measure_col)))

    normalized_df <- normalized_df %>%
        merge(
            df_long %>% select(-one_of(setdiff(
                names(normalized_df),
                c(feature_id_col, sample_id_col, measure_col)
            ))),
            by = c(feature_id_col, sample_id_col)
        )

    default_cols <- names(normalized_df)
    minimal_cols <- c(sample_id_col, feature_id_col, measure_col, old_measure_col)

    if (!is.null(qual_col) && qual_col %in% names(normalized_df)) {
        minimal_cols <- c(minimal_cols, qual_col)
    }
    normalized_df <- subset_keep_cols(
        normalized_df,
        keep_all,
        default_cols = default_cols,
        minimal_cols = minimal_cols
    )

    return(normalized_df)
}

#'
#' @export
#' @rdname normalize
#'
normalize_sample_medians_dm <- function(data_matrix,
                                        sample_annotation = NULL,
                                        sample_id_col = "FullRunName",
                                        group_col = NULL,
                                        inside_batch = FALSE) {
    norm_res <- .pb_median_center_matrix(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        group_col = group_col,
        inside_batch = inside_batch
    )
    norm_res$matrix
}

.pb_median_center_matrix <- function(data_matrix,
                                     sample_annotation = NULL,
                                     sample_id_col = "FullRunName",
                                     group_col = NULL,
                                     inside_batch = FALSE) {
    if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
    }
    if (!inside_batch) {
        sample_medians <- apply(data_matrix, 2, median, na.rm = TRUE)
        global_median <- median(data_matrix, na.rm = TRUE)
        if (is.na(global_median)) {
            return(list(
                matrix = data_matrix,
                sample_medians = sample_medians,
                reference_medians = rep(NA_real_, length(sample_medians))
            ))
        }
        adjustments <- global_median - sample_medians
        adjustments[is.na(adjustments)] <- 0
        centered <- sweep(data_matrix, 2, adjustments, FUN = "+")
        reference_medians <- rep(global_median, length(sample_medians))
        return(list(
            matrix = centered,
            sample_medians = sample_medians,
            reference_medians = reference_medians
        ))
    }

    if (is.null(group_col) || !nzchar(group_col)) {
        stop("`group_col` must be provided when `inside_batch = TRUE`.")
    }
    if (is.null(sample_annotation)) {
        stop("`sample_annotation` must be provided when `inside_batch = TRUE`.")
    }
    sample_annotation <- as.data.frame(sample_annotation)
    if (!(sample_id_col %in% names(sample_annotation))) {
        stop("`sample_id_col` not found in supplied `sample_annotation`.")
    }
    if (!(group_col %in% names(sample_annotation))) {
        stop("`group_col` not found in supplied `sample_annotation`.")
    }
    sa_unique <- sample_annotation %>%
        select(dplyr::all_of(c(sample_id_col, group_col))) %>%
        distinct()
    if (anyDuplicated(sa_unique[[sample_id_col]])) {
        dup_ids <- unique(sa_unique[[sample_id_col]][duplicated(sa_unique[[sample_id_col]])])
        stop("Duplicate entries for samples in `sample_annotation`: ", paste(dup_ids, collapse = ", "))
    }
    id_match <- match(colnames(data_matrix), sa_unique[[sample_id_col]])
    if (any(is.na(id_match))) {
        missing_ids <- colnames(data_matrix)[is.na(id_match)]
        stop("Missing annotations for samples: ", paste(missing_ids, collapse = ", "))
    }
    groups <- sa_unique[[group_col]][id_match]
    if (any(is.na(groups))) {
        stop("`group_col` must not contain NA values when `inside_batch = TRUE`.")
    }

    centered <- data_matrix
    sample_medians <- rep(NA_real_, ncol(data_matrix))
    reference_medians <- rep(NA_real_, ncol(data_matrix))
    split_indices <- split(seq_along(groups), groups, drop = FALSE)

    for (idx in split_indices) {
        if (!length(idx)) {
            next
        }
        sub_matrix <- data_matrix[, idx, drop = FALSE]
        valid_rows <- rowSums(!is.na(sub_matrix)) > 0
        if (!any(valid_rows)) {
            next
        }
        medians <- apply(sub_matrix[valid_rows, , drop = FALSE], 2, median, na.rm = TRUE)
        global_median <- median(sub_matrix[valid_rows, , drop = FALSE], na.rm = TRUE)
        if (is.na(global_median)) {
            next
        }
        adjustments <- global_median - medians
        adjustments[is.na(adjustments)] <- 0
        centered[, idx] <- sweep(sub_matrix, 2, adjustments, FUN = "+")
        sample_medians[idx] <- medians
        reference_medians[idx] <- global_median
    }

    list(
        matrix = centered,
        sample_medians = sample_medians,
        reference_medians = reference_medians
    )
}

#'
#' @export
#' @rdname normalize
normalize_sample_medians_df <- function(df_long,
                                        feature_id_col = "peptide_group_label",
                                        sample_id_col = "FullRunName",
                                        measure_col = "Intensity",
                                        no_fit_imputed = FALSE,
                                        qual_col = NULL,
                                        qual_value = 2,
                                        keep_all = "default",
                                        sample_annotation = NULL,
                                        group_col = NULL,
                                        inside_batch = FALSE) {
    df_processed <- df_long
    if (no_fit_imputed) {
        if (!(qual_col %in% names(df_processed))) {
            stop("imputed value flag column (qual_col) is not in the data frame!")
        }
        message("removing imputed values (requants) from the matrix")
        df_processed <- df_processed %>%
            mutate(!!sym(measure_col) := ifelse(!!sym(qual_col) == qual_value,
                NA, !!sym(measure_col)
            ))
    } else {
        if (!is.null(qual_col) && (qual_col %in% names(df_processed))) {
            warning("imputed value (requant) column is in the data, are you sure you
              want to use imputed (requant) values in sample median inference?")
        }
    }

    old_measure_col <- paste("preNorm", measure_col, sep = "_")

    if (!inside_batch) {
        normalized_df <- df_processed %>%
            group_by_at(vars(one_of(sample_id_col))) %>%
            mutate(median_run = median(!!sym(measure_col), na.rm = TRUE)) %>%
            ungroup() %>%
            mutate(
                median_global = median(!!sym(measure_col), na.rm = TRUE),
                !!(old_measure_col) := !!(sym(measure_col))
            ) %>%
            mutate(diff_norm = median_global - median_run) %>%
            mutate(!!(sym(measure_col)) := !!(sym(measure_col)) + diff_norm)

        default_cols <- names(normalized_df)
        minimal_cols <- c(sample_id_col, feature_id_col, measure_col, old_measure_col)

        if (!is.null(qual_col) && qual_col %in% names(normalized_df)) {
            minimal_cols <- c(minimal_cols, qual_col)
        }
        normalized_df <- subset_keep_cols(
            normalized_df,
            keep_all,
            default_cols = default_cols,
            minimal_cols = minimal_cols
        )

        return(normalized_df)
    }

    if (is.null(group_col) || !nzchar(group_col)) {
        stop("`group_col` must be provided when `inside_batch = TRUE`.")
    }

    grouping_annotation <- sample_annotation
    if (!is.null(grouping_annotation)) {
        grouping_annotation <- as.data.frame(grouping_annotation)
        if (!(group_col %in% names(grouping_annotation))) {
            stop("`group_col` not found in supplied `sample_annotation`.")
        }
        grouping_annotation <- grouping_annotation %>%
            select(dplyr::all_of(c(sample_id_col, group_col))) %>%
            distinct()
    } else {
        if (!(group_col %in% names(df_processed))) {
            stop("Provide `sample_annotation` or include `group_col` in `df_long` when `inside_batch = TRUE`.")
        }
        grouping_annotation <- df_processed %>%
            select(dplyr::all_of(c(sample_id_col, group_col))) %>%
            distinct()
    }

    data_matrix <- long_to_matrix(
        df_processed,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col
    )

    norm_res <- .pb_median_center_matrix(
        data_matrix = data_matrix,
        sample_annotation = grouping_annotation,
        sample_id_col = sample_id_col,
        group_col = group_col,
        inside_batch = inside_batch
    )
    normalized_matrix <- norm_res$matrix
    sample_medians <- norm_res$sample_medians
    reference_medians <- norm_res$reference_medians

    normalized_df <- matrix_to_long(
        normalized_matrix,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col
    )

    df_processed <- df_processed %>%
        rename(!!(old_measure_col) := !!(sym(measure_col)))

    median_info <- data.frame(
        sample_id_temp = colnames(normalized_matrix),
        median_run = sample_medians,
        median_global = reference_medians,
        stringsAsFactors = FALSE
    )
    names(median_info)[1] <- sample_id_col

    normalized_df <- normalized_df %>%
        left_join(median_info, by = sample_id_col) %>%
        mutate(diff_norm = median_global - median_run)

    extra_cols <- setdiff(names(df_processed), names(normalized_df))
    if (length(extra_cols)) {
        keep_cols <- unique(c(feature_id_col, sample_id_col, extra_cols))
        normalized_df <- normalized_df %>%
            left_join(
                df_processed %>%
                    select(dplyr::all_of(keep_cols)),
                by = c(feature_id_col, sample_id_col)
            )
    }

    default_cols <- names(normalized_df)
    minimal_cols <- c(sample_id_col, feature_id_col, measure_col, old_measure_col)

    if (!is.null(qual_col) && qual_col %in% names(normalized_df)) {
        minimal_cols <- c(minimal_cols, qual_col)
    }
    normalized_df <- subset_keep_cols(
        normalized_df,
        keep_all,
        default_cols = default_cols,
        minimal_cols = minimal_cols
    )

    return(normalized_df)
}

#'
#' @export
#' @rdname normalize
normalize_data_dm <- function(data_matrix,
                              normalize_func = c("quantile", "medianCentering"),
                              log_base = NULL, offset = 1) {
    normalize_func <- match.arg(normalize_func)
    if (!is.null(log_base)) {
        data_matrix <- log_transform_dm(data_matrix, log_base = log_base, offset = offset)
    }

    if (normalize_func == "quantile") {
        normalized_matrix <- quantile_normalize_dm(data_matrix)
    } else if (normalize_func == "medianCentering") {
        normalized_matrix <- normalize_sample_medians_dm(data_matrix)
    } else {
        stop("Only quantile and median centering normalization methods implemented")
    }

    return(normalized_matrix)
}

#' @export
#' @rdname normalize
normalize_data_df <- function(df_long,
                              normalize_func = c("quantile", "medianCentering"),
                              log_base = NULL, offset = 1,
                              feature_id_col = "peptide_group_label",
                              sample_id_col = "FullRunName",
                              measure_col = "Intensity",
                              no_fit_imputed = TRUE,
                              qual_col = NULL,
                              qual_value = 2,
                              keep_all = "default") {
    normalize_func <- match.arg(normalize_func)


    if (!is.null(log_base)) {
        df_long <- log_transform_df(df_long, log_base = log_base, offset = offset)
    }

    if (normalize_func == "quantile") {
        normalized_df <- quantile_normalize_df(
            df_long,
            feature_id_col = feature_id_col,
            measure_col = measure_col,
            sample_id_col = sample_id_col,
            keep_all = keep_all,
            no_fit_imputed = no_fit_imputed,
            qual_col = qual_col,
            qual_value = qual_value
        )
    } else if (normalize_func == "medianCentering") {
        normalized_df <- normalize_sample_medians_df(
            df_long,
            sample_id_col = sample_id_col,
            measure_col = measure_col,
            keep_all = keep_all,
            no_fit_imputed = no_fit_imputed,
            qual_col = qual_col,
            qual_value = qual_value
        )
    } else {
        stop("Only quantile and median centering normalization methods implemented")
    }

    return(normalized_df)
}
