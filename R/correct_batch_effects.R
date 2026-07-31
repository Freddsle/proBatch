#' @title Batch correction methods for normalized data
#'
#' @description
#' Batch correction of normalized data to bring each feature within each batch
#' to a comparable scale. The following methods are available:
#' \enumerate{
#'   \item \strong{Per-feature centering by median/mean}:
#'     \code{\link{center_feature_batch}()} — unified API for long
#'     (\code{"long"}) and wide (\code{"wide"}) data via \code{format}, and for
#'     medians/means via \code{stat}. \emph{Deprecated wrappers}:
#'     \code{center_feature_batch_medians_df()},
#'     \code{center_feature_batch_means_df()},
#'     \code{center_feature_batch_medians_dm()},
#'     \code{center_feature_batch_means_dm()} now forward to
#'     \code{center_feature_batch()} and emit deprecation warnings.
#'
#'   \item \strong{Discrete batch correction with ComBat}:
#'     \code{\link{correct_with_ComBat}()} adjusts for discrete batch effects
#'     using parametric or non-parametric empirical Bayes (Johnson et al.,
#'     2007). \strong{Missing data:} ComBat requires an \emph{NA-free} matrix.
#'     If your data contain missing values, use
#'     \code{fill_the_missing = "drop_features"} to drop NA-containing features
#'     or \code{fill_the_missing = "fill"} with an explicit \code{fill_value}
#'     to impute before calling ComBat.
#'
#'   \item \strong{Linear batch correction with limma}:
#'     \code{\link{correct_with_removeBatchEffect}()} removes linear batch
#'     effects via \code{limma::removeBatchEffect}. \strong{Missing data:} NA
#'     values in the \emph{expression/intensity matrix} are allowed; the
#'     \emph{design matrix} (batch/covariates) must be free of missing values.
#'
#'   \item \strong{Continuous drift correction}:
#'     \code{\link{adjust_batch_trend_df}()} fits and removes within-batch
#'     trends (e.g., LOESS) and is typically followed by a discrete adjustment
#'     such as \code{center_feature_batch()}, \code{correct_with_ComBat()}, or
#'     \code{correct_with_removeBatchEffect()}.
#'
#' }
#'
#' Alternatively, use the wrapper \code{\link{correct_batch_effects}()} to
#' combine continuous and/or discrete corrections in one call.
#'
#' @inheritParams proBatch
#'
#' @param return_fit_df Logical; whether \code{\link{adjust_batch_trend_dm}()}
#'   returns the \code{fit_df} (for curve inspection)
#'   alongside the corrected matrix.
#' @param fit_func Function used for trend fitting (e.g.,
#'   \code{"loess_regression"}).
#' @param min_measurements Minimum number of samples per
#'   batch required for fitting.
#' @param par.prior Logical; use parametric prior (ComBat) or non-parametric.
#' @param continuous_func Which function to use for the
#'   continuous fit (currently only \code{"loess_regression"}); set \code{NULL}
#'   if order-associated drift correction is not required.
#' @param discrete_func Which function to use for discrete batch effects in the
#'   wrapper: one of \code{"MedianCentering"}, \code{"MeanCentering"},
#'   \code{"ComBat"}, or \code{"removeBatchEffect"}.
#' @param fill_the_missing Missing-value policy applied \emph{before} discrete
#'   correction: \code{"error"} (default), \code{"keep"},
#'   \code{"drop_features"}, or \code{"fill"}. With \code{"keep"}, missing
#'   values are passed through; ComBat still requires an NA-free matrix, while
#'   removeBatchEffect can model expression values containing NA. For one
#'   release, \code{FALSE}, a numeric scalar, and \code{"remove"}, \code{"rm"},
#'   or \code{"REMOVE"} are translated with a deprecation warning. Explicit
#'   \code{NULL} is an error.
#' @param fill_value Finite numeric scalar used only with
#'   \code{fill_the_missing = "fill"}.
#' @param ... Additional parameters. The trend-adjustment functions forward
#'   them to the chosen \code{fit_func}. \code{correct_batch_effects()}
#'   forwards them to that fitter when continuous correction is enabled and to
#'   \code{limma::removeBatchEffect()} when that discrete method is selected.
#'
#' @return
#' Returns data in the same format as input (\code{data_matrix} or
#' \code{df_long}). For long format, original values from \code{measure_col}
#' are preserved in \code{"preBatchCorr_[measure_col]"} and corrected values
#' are written to \code{measure_col}.
#'
#' The function \code{\link{adjust_batch_trend_dm}()}, if
#' \code{return_fit_df = TRUE}, returns a list with:
#' \enumerate{
#'   \item \code{corrected_dm} — corrected data matrix
#'   \item \code{fit_df} — data frame to inspect fitted curves
#' }
#'
#' @examples
#' data(
#'     list = c("example_sample_annotation", "example_proteome"),
#'     package = "proBatch"
#' )
#' test_peptides <- unique(example_proteome$peptide_group_label)[1:3]
#' test_df <- subset(example_proteome, peptide_group_label %in% test_peptides)
#'
#' # 1) Per-feature median centering
#' median_centered_df <- center_feature_batch(
#'     x = test_df,
#'     sample_annotation = example_sample_annotation,
#'     format = "long",
#'     stat = "medians",
#'     no_fit_imputed = FALSE
#' )
#'
#' # 2) ComBat (discrete) — drop NA features first if needed
#' combat_corrected_df <- correct_with_ComBat(
#'     x = test_df,
#'     sample_annotation = example_sample_annotation,
#'     format = "long",
#'     fill_the_missing = "drop_features",
#'     no_fit_imputed = FALSE
#' )
#'
#' # 3) Continuous drift correction (LOESS), then discrete centering if desired
#' adjusted_df <- adjust_batch_trend_df(
#'     df_long = test_df,
#'     sample_annotation = example_sample_annotation,
#'     span = 0.7, min_measurements = 8,
#'     no_fit_imputed = FALSE
#' )
#' plot_fit <- plot_with_fitting_curve(
#'     unique(adjusted_df$peptide_group_label),
#'     df_long = adjusted_df,
#'     measure_col = "preTrendFit_Intensity",
#'     fit_df = adjusted_df,
#'     sample_annotation = example_sample_annotation
#' )
#'
#' # 4) One-call wrapper (continuous + discrete)
#' batch_corrected_df <- correct_batch_effects(
#'     x = test_df, sample_annotation = example_sample_annotation,
#'     format = "long",
#'     continuous_func = "loess_regression",
#'     discrete_func = "MedianCentering",
#'     batch_col = "MS_batch",
#'     fill_the_missing = "keep",
#'     no_fit_imputed = FALSE,
#'     span = 0.7, min_measurements = 8
#' )
#'
#' @seealso
#' \code{\link{center_feature_batch}}, \code{\link{adjust_batch_trend_df}},
#' \code{\link{adjust_batch_trend_dm}}, \code{\link{correct_with_ComBat}},
#' \code{\link{correct_with_removeBatchEffect}},
#' \code{\link{correct_batch_effects}}
#'
#' @references
#' Johnson WE, Li C, Rabinovic A (2007). Adjusting batch effects in microarray
#' expression data using empirical Bayes methods. \emph{Biostatistics}
#' 8(1):118–127. Smyth GK (2025). \emph{limma User's Guide}, Bioconductor (see
#' removeBatchEffect / lmFit). Leek JT et al. (2024).
#' \emph{sva} vignette, Bioconductor.
#'
#' @name correct_batch_effects
NULL

#' @title Center features per-batch by median/mean (unified)
#' @md
#' @description Centers each feature *within each batch* to the global
#' location (median/mean) of that feature. Works for long-data
#' frames and wide matrices.
#'
#' @inheritParams correct_batch_effects
#' @param x Either a long data.frame or a numeric matrix (features in rows,
#'   samples in columns), depending on `format`.
#' @param format One of `"long"` or `"wide"`. `"long"` expects a data.frame
#'   with columns `sample_id_col`, `feature_id_col`, `measure_col`. `"wide"`
#'   expects a numeric matrix and uses `matrix_to_long()` / `long_to_matrix()`.
#' @param stat One of `"medians"` or `"means"`. Aliases `"median"`/`"mean"`
#'   are accepted.
#' @param keep_all Passed to `subset_keep_cols()` for long
#'   format, ignored for wide.
#' @param no_fit_imputed If `TRUE` and `qual_col` is provided, imputed values
#'   are masked during location estimation (original values remain unchanged).
#' @param qual_col,qual_value Column and value that flag imputed entries.
#'
#' @return If `format = "long"`, a data.frame; if `format = "wide"`, a matrix.
#'
#' @details
#' For `"wide"`, conversion uses `matrix_to_long()` / `long_to_matrix()`. For
#' `"long"`, the function adds `"preBatchCorr_[measure_col]"`
#' and the diagnostic columns:
#' - for medians: `median_batch`, `median_global`, `diff_medians`
#' - for means:   `mean_batch`,   `mean_global`,   `diff_means`
#'
#' @examples
#' # LONG
#' # Load necessary datasets
#' data(
#'     list = c("example_sample_annotation", "example_proteome"),
#'     package = "proBatch"
#' )
#' out_long <- center_feature_batch(
#'     x = example_proteome, sample_annotation = example_sample_annotation,
#'     format = "long", stat = "medians",
#'     sample_id_col = "FullRunName", batch_col = "MS_batch",
#'     feature_id_col = "peptide_group_label", measure_col = "Intensity",
#'     no_fit_imputed = FALSE
#' )
#'
#' # WIDE
#' data(example_proteome_matrix, package = "proBatch")
#' out_wide <- center_feature_batch(
#'     x = example_proteome_matrix,
#'     sample_annotation = example_sample_annotation,
#'     format = "wide", stat = "means",
#'     sample_id_col = "FullRunName", batch_col = "MS_batch",
#'     feature_id_col = "peptide_group_label", measure_col = "Intensity",
#'     no_fit_imputed = FALSE
#' )
#'
#' @export
center_feature_batch <- function(
    x,
    sample_annotation = NULL,
    format = c("long", "wide"),
    stat = c("medians", "means"),
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    keep_all = "default",
    no_fit_imputed = TRUE,
    qual_col = NULL,
    qual_value = NULL
) {
    format <- match.arg(format)
    stat_in <- tolower(stat[1])
    stat_in <- switch(
        stat_in,
        medians = "median",
        median = "median",
        means = "mean",
        mean = "mean",
        stop("`stat` must be one of 'medians' or 'means'.")
    )
    stat_names <- if (identical(stat_in, "median")) {
        c(
            batch = "median_batch",
            global = "median_global",
            diff = "diff_medians"
        )
    } else {
        c(batch = "mean_batch", global = "mean_global", diff = "diff_means")
    }

    if (identical(format, "wide")) {
        if (!is.matrix(x) || !is.numeric(x)) {
            stop(
                "format='wide' requires a numeric matrix with features in rows and samples in columns."
            )
        }
        if (is.null(sample_annotation) || is.null(colnames(x))) {
            stop(
                "format='wide' requires `sample_annotation` and column names on the matrix."
            )
        }
        if (!all(colnames(x) %in% sample_annotation[[sample_id_col]])) {
            stop(
                "Not all matrix column names found in sample_annotation[[sample_id_col]]."
            )
        }
        # Convert wide -> long
        df_long <- matrix_to_long(
            data_matrix = x,
            feature_id_col = feature_id_col,
            measure_col = measure_col,
            sample_id_col = sample_id_col
        )
    } else {
        df_long <- x
    }

    # format == "long"
    if (!is.data.frame(df_long) && identical(format, "long")) {
        stop("format='long' requires a long data format")
    }
    corrected_long <- .center_feature_batch_df_core(
        df_long = df_long,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        keep_all = keep_all,
        no_fit_imputed = no_fit_imputed,
        qual_col = qual_col,
        qual_value = qual_value,
        stat = stat_in,
        stat_names = stat_names
    )

    if (identical(format, "wide")) {
        # Convert long -> wide
        out_wide <- long_to_matrix(
            corrected_long,
            feature_id_col = feature_id_col,
            measure_col = measure_col,
            sample_id_col = sample_id_col
        )
        return(out_wide)
    } # else  format == "long"
    corrected_long
}

.center_feature_batch_df_core <- function(
    df_long,
    sample_annotation = NULL,
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    keep_all = "default",
    stat = c("median", "mean"),
    no_fit_imputed = TRUE,
    qual_col = NULL,
    qual_value = NULL,
    stat_names = c(
        batch = "median_batch",
        global = "median_global",
        diff = "diff_medians"
    )
) {
    stat <- match.arg(stat)
    original_cols <- names(df_long)

    # Merge/check annotations; ensure batch present post-merge
    df_long <- check_sample_consistency(
        sample_annotation,
        sample_id_col,
        df_long,
        batch_col,
        order_col = NULL,
        facet_col = NULL,
        merge = TRUE
    )

    if (no_fit_imputed && is.null(qual_col)) {
        warning(
            "`qual_col` is NULL, setting `no_fit_imputed = FALSE`; imputed flags will be ignored."
        )
        no_fit_imputed <- FALSE
    }

    # Choose location function
    summariser <- switch(stat, median = median, mean = mean)

    # Optionally mask imputed values during inference
    tmp_col <- NULL
    if (isTRUE(no_fit_imputed)) {
        df_long <- .mask_imputed_measure(
            df_long = df_long,
            measure_col = measure_col,
            qual_col = qual_col,
            qual_value = qual_value
        )
        tmp_col <- attr(df_long, "temp_measure_col")
    }
    measure_for_inference <- if (!is.null(tmp_col)) tmp_col else measure_col
    old_measure_col <- .make_pre_col("preBatchCorr", measure_col)

    # Compute per-(batch,feature) and per-feature locations, then apply shift
    corrected_df <- df_long %>%
        group_by(across(any_of(c(batch_col, feature_id_col)))) %>%
        mutate(
            !!stat_names["batch"] := summariser(
                .data[[measure_for_inference]],
                na.rm = TRUE
            )
        ) %>%
        ungroup() %>%
        group_by(across(any_of(feature_id_col))) %>%
        mutate(
            !!stat_names["global"] := summariser(
                .data[[measure_for_inference]],
                na.rm = TRUE
            )
        ) %>%
        ungroup() %>%
        mutate(
            !!stat_names["diff"] := .data[[stat_names["global"]]] -
                .data[[stat_names["batch"]]]
        ) %>%
        rename(!!old_measure_col := !!measure_col) %>%
        mutate(
            !!measure_col := .data[[old_measure_col]] +
                .data[[stat_names["diff"]]]
        )

    # Drop temporary masked column if present
    if (!is.null(tmp_col) && tmp_col %in% names(corrected_df)) {
        corrected_df <- select(corrected_df, -all_of(tmp_col))
    }

    # Column retention for long format
    default_cols <- unique(c(
        original_cols,
        batch_col,
        old_measure_col,
        unname(stat_names["batch"]),
        unname(stat_names["global"]),
        unname(stat_names["diff"])
    ))
    minimal_cols <- unique(c(
        sample_id_col,
        feature_id_col,
        measure_col,
        old_measure_col,
        batch_col,
        unname(stat_names["batch"]),
        unname(stat_names["diff"])
    ))
    if (!is.null(qual_col) && qual_col %in% names(corrected_df)) {
        default_cols <- c(default_cols, qual_col)
        minimal_cols <- c(minimal_cols, qual_col)
    }

    subset_keep_cols(
        df = corrected_df,
        keep_all = keep_all,
        default_cols = default_cols,
        minimal_cols = minimal_cols
    )
}

#' @title Adjust batch trend with custom (continuous) fit
#' @description Adjust batch signal trend with the custom (continuous) fit.
#' Should be followed by discrete corrections, e.g.
#' \code{center_feature_batch()} or \code{correct_with_ComBat()}. Available for
#' both long format data frame (\code{adjust_batch_trend_df()}) and data matrix
#' (\code{adjust_batch_trend_dm()}).
#' @export
#' @rdname correct_batch_effects
#'
#' @seealso \code{\link{fit_nonlinear}}, \code{\link{plot_with_fitting_curve}}
adjust_batch_trend_dm <- function(
    data_matrix,
    sample_annotation,
    batch_col = "MS_batch",
    feature_id_col = "peptide_group_label",
    sample_id_col = "FullRunName",
    measure_col = "Intensity",
    order_col = "order",
    fit_func = "loess_regression",
    return_fit_df = TRUE,
    no_fit_imputed = TRUE,
    qual_col = NULL,
    qual_value = NULL,
    min_measurements = 8,
    ...
) {
    df_long <- matrix_to_long(
        data_matrix,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col
    )

    corrected_data <- adjust_batch_trend_df(
        df_long,
        sample_annotation,
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        order_col = order_col,
        fit_func = fit_func,
        no_fit_imputed = no_fit_imputed,
        qual_col = qual_col,
        qual_value = qual_value,
        min_measurements = min_measurements,
        ...
    )

    corrected_df <- corrected_data
    corrected_dm <- long_to_matrix(
        corrected_df,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col
    )

    if (return_fit_df) {
        # extract only the columns relevant for inspecting the fit
        # (only non-empty columns)
        fit_columns <- c(
            sample_id_col,
            feature_id_col,
            batch_col,
            order_col,
            "fit"
        )
        # if any of the fit columns are not present in the
        # corrected_df, remove them and warn
        if (any(!fit_columns %in% names(corrected_df))) {
            missing_cols <- fit_columns[!fit_columns %in% names(corrected_df)]
            message(
                "The following columns are not present in the corrected_df and will be removed from fit_df: ",
                toString(missing_cols)
            )
            fit_columns <- fit_columns[fit_columns %in% names(corrected_df)]
        }
        fit_df <- corrected_df[, fit_columns, drop = FALSE]
        return(list(
            corrected_dm = corrected_dm,
            fit_df = fit_df
        ))
    } else {
        return(corrected_dm)
    }
}

#'
#' @export
#' @rdname correct_batch_effects
#'
#' @seealso \code{\link{fit_nonlinear}}, \code{\link{plot_with_fitting_curve}}
adjust_batch_trend_df <- function(
    df_long,
    sample_annotation = NULL,
    batch_col = "MS_batch",
    feature_id_col = "peptide_group_label",
    sample_id_col = "FullRunName",
    measure_col = "Intensity",
    order_col = "order",
    keep_all = "default",
    fit_func = "loess_regression",
    no_fit_imputed = TRUE,
    qual_col = NULL,
    qual_value = NULL,
    min_measurements = 8,
    ...
) {
    if (is(df_long, "ProBatchFeatures")) {
        object <- df_long
        assay_name <- .pb_resolve_assay_for_input(object)
        df_long <- .pb_pbf_to_long(
            object = object,
            assay_name = assay_name,
            feature_id_col = feature_id_col,
            sample_id_col = sample_id_col,
            measure_col = measure_col
        )
        sample_annotation <- .pb_default_sample_annotation(
            object = object,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            sample_ids = unique(df_long[[sample_id_col]])
        )
    }

    original_cols <- names(df_long)

    df_long <- check_sample_consistency(
        sample_annotation,
        sample_id_col,
        df_long,
        batch_col,
        order_col = order_col,
        facet_col = NULL,
        merge = TRUE
    )

    if (no_fit_imputed) {
        if (is.null(qual_col)) {
            warning(
                "`qual_col` is NULL, setting `no_fit_imputed = FALSE` so imputed flags are ignored."
            )
            no_fit_imputed <- FALSE
        } else if (!(qual_col %in% names(df_long))) {
            stop(
                "imputed value flag column (qual_col) is not in the data frame!"
            )
        }
    } else {
        if (!is.null(qual_col)) {
            # flags provided but explicitly ignored --
            # keep behavior, inform user
            warning(
                "`qual_col` provided but `no_fit_imputed = FALSE`; imputed flags will be ignored for curve fitting."
            )
        }
    }

    # If a batch column is requested, ensure it's present
    # after consistency checks.
    # If no per-batch stratification; fit per-feature across all samples
    group_vars <- c(feature_id_col, if (!is.null(batch_col)) batch_col)
    corrected_df <- df_long %>%
        nest(data = -all_of(group_vars)) %>%
        mutate(
            fit = pmap(
                list(data = .data$data),
                function(data) {
                    fit_nonlinear(
                        df_feature_batch = data,
                        measure_col = measure_col,
                        order_col = order_col,
                        fit_func = fit_func,
                        no_fit_imputed = no_fit_imputed,
                        qual_col = qual_col,
                        qual_value = qual_value,
                        min_measurements = min_measurements,
                        ...
                    )
                }
            )
        )
    old_measure_col <- .make_pre_col("preTrendFit", measure_col)

    corrected_df <- corrected_df %>%
        # only unnest 'data' (original rows) and 'fit' (vector of fitted values)
        unnest(cols = c(data, fit)) %>%
        group_by(across(any_of(c(feature_id_col, batch_col)))) %>%
        mutate(mean_fit = mean(fit, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(
            diff_fit = mean_fit - fit,
            diff.na = ifelse(is.na(diff_fit), 0, diff_fit)
        ) %>%
        rename(!!sym(old_measure_col) := !!sym(measure_col)) %>%
        # Conditional shift: use diff.na so NA in fit
        # doesn't propagate the shift
        mutate(
            !!sym(measure_col) := .data[["diff.na"]] +
                .data[[old_measure_col]]
        ) %>%
        select(-any_of("diff.na"))

    default_cols <- unique(c(
        original_cols,
        old_measure_col,
        "fit",
        if (!is.null(batch_col)) batch_col
    ))
    minimal_cols <- unique(c(
        sample_id_col,
        feature_id_col,
        measure_col,
        old_measure_col,
        "fit",
        if (!is.null(batch_col) && batch_col %in% names(corrected_df)) batch_col
    ))

    if (!is.null(qual_col) && qual_col %in% names(corrected_df)) {
        default_cols <- c(default_cols, qual_col)
        minimal_cols <- c(minimal_cols, qual_col)
    }

    corrected_df <- subset_keep_cols(
        corrected_df,
        keep_all,
        default_cols = default_cols,
        minimal_cols = minimal_cols
    )

    return(corrected_df)
}

##############################################################################
# ComBat-based batch correction

#' @title ComBat-based batch correction (unified)
#' @description Adjusts for discrete batch effects using \code{ComBat}.
#' Works for long-format data frames and wide matrices via \code{format}.
#' Optionally accepts covariates through
#' \code{covariates_cols} (passed as \code{mod}).
#'
#' @inheritParams proBatch
#' @param x Data in long (\code{data.frame}) or wide (\code{matrix}) form,
#'   controlled by \code{format}.
#' @param format One of \code{"long"} or \code{"wide"}.
#' @param par.prior Logical; ComBat parametric prior (vs non-parametric).
#' @param covariates_cols Optional character vector of
#'   \code{sample_annotation} columns included in \code{mod} for ComBat
#'   (biological or nuisance covariates).
#' @param fill_the_missing Missing-value policy prior to ComBat:
#'   \code{"error"} (default), \code{"keep"}, \code{"drop_features"}, or
#'   \code{"fill"}. ComBat still requires the resulting matrix to be NA-free.
#' @param keep_all For long format, columns to retain (see
#'   \code{subset_keep_cols()}).
#' @param no_fit_imputed If \code{TRUE} and \code{qual_col}
#'   provided, masked values are excluded when building the matrix (original
#'   values still corrected).
#' @param fill_value Finite numeric scalar used only with
#'   \code{fill_the_missing = "fill"}.
#' @param ... Further arguments passed to \code{sva::ComBat()}.
#'
#' @return Matrix if \code{format="wide"}, data.frame if \code{format="long"}.
#' @examples
#' data(
#'     list = c("example_sample_annotation", "example_proteome"),
#'     package = "proBatch"
#' )
#' feature_ids <- c(
#'     "10062_NVGVSFYADKPEVTQEQK_3",
#'     "10063_NVGVSFYADKPEVTQEQKK_3"
#' )
#' example_subset <- example_proteome[
#'     example_proteome$peptide_group_label %in% feature_ids,
#' ]
#' combat_corrected <- correct_with_ComBat(
#'     x = example_subset,
#'     sample_annotation = example_sample_annotation,
#'     format = "long",
#'     fill_the_missing = "error"
#' )
#' @references Johnson WE et al. (2007) \emph{Biostatistics}
#'   8(1):118–127; \emph{sva} vignette.
#' @export
correct_with_ComBat <- function(
    x,
    sample_annotation = NULL,
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    format = c("long", "wide"),
    par.prior = TRUE,
    covariates_cols = NULL,
    fill_the_missing = "error",
    keep_all = "default",
    no_fit_imputed = TRUE,
    qual_col = NULL,
    qual_value = NULL,
    ...,
    fill_value = NULL
) {
    format <- match.arg(format)
    missing_policy <- .pb_normalize_missing_policy(
        fill_the_missing,
        fill_value = fill_value,
        argument = "fill_the_missing"
    )
    fill_the_missing <- missing_policy$policy
    fill_value <- missing_policy$fill_value

    if (identical(format, "wide")) {
        if (!is.matrix(x)) {
            stop("format='wide' requires a numeric matrix.")
        }
        corrected_matrix <- .combat_matrix_step(
            data_matrix = x,
            sample_annotation = sample_annotation,
            batch_col = batch_col,
            sample_id_col = sample_id_col,
            par.prior = par.prior,
            fill_the_missing = fill_the_missing,
            fill_value = fill_value,
            covariates_cols = covariates_cols,
            ...
        )
        return(corrected_matrix)
    }

    qual_for_matrix <- if (no_fit_imputed) qual_col else NULL
    qual_val_for_matrix <- if (no_fit_imputed) qual_value else NULL
    prep <- .pb_prepare_long_matrix(
        df_long = x,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        batch_col = batch_col,
        fill_the_missing = fill_the_missing,
        fill_value = fill_value,
        warning_message = "ComBat cannot operate with missing values in the matrix",
        qual_col = qual_for_matrix,
        qual_value = qual_val_for_matrix,
        error_message = "format='long' requires a data.frame."
    )
    df_long <- prep$df_long
    sample_annotation <- prep$sample_annotation
    data_matrix <- prep$data_matrix
    original_cols <- prep$original_cols

    # ComBat on matrix (method ensures numeric & SA alignment)
    corrected_matrix <- .combat_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        batch_col = batch_col,
        sample_id_col = sample_id_col,
        par.prior = par.prior,
        fill_the_missing = fill_the_missing,
        fill_value = fill_value,
        covariates_cols = covariates_cols,
        ...
    )

    .post_correction_to_long(
        corrected_matrix = corrected_matrix,
        df_long = df_long,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col,
        original_cols = original_cols,
        keep_all = keep_all
    )
}

# ---- removeBatchEffect unified ------------------------------------------------

#' @title Batch effect correction with limma::removeBatchEffect (unified)
#' @description Removes batch-associated linear effects with
#'   removeBatchEffect from limma.
#' Works for long or wide via \code{format}. Use \code{covariates_cols} to keep
#' biological effects in the design (not removed).
#' @inheritParams correct_with_ComBat
#' @inheritParams proBatch
#' @param covariates_cols Optional \code{sample_annotation} columns for the
#'   design matrix (biological or nuisance covariates).
#' @param fill_the_missing Missing-value policy applied before modeling:
#'   \code{"error"} (default), \code{"keep"}, \code{"drop_features"}, or
#'   \code{"fill"}. With \code{"keep"}, NA values are passed unchanged to
#'   \code{limma::removeBatchEffect()}. The design matrix (\code{batch_col} and
#'   \code{covariates_cols}) must be NA-free for every policy.
#' @param ... Further arguments passed to
#'   \code{limma::removeBatchEffect()}.
#'
#' @return Matrix if \code{format="wide"}, data.frame if \code{format="long"}
#'   with batch effects removed
#' @examples
#' data(
#'     list = c("example_sample_annotation", "example_proteome_matrix"),
#'     package = "proBatch"
#' )
#' batch_corrected_matrix <- correct_with_removeBatchEffect(
#'     example_proteome_matrix,
#'     example_sample_annotation,
#'     format = "wide",
#'     batch_col = "MS_batch",
#'     covariates_cols = c("Diet", "Sex"),
#'     fill_the_missing = "drop_features"
#' )
#' @seealso \code{\link[limma:removeBatchEffect]{removeBatchEffect}}
#' @export
correct_with_removeBatchEffect <- function(
    x,
    sample_annotation,
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    format = c("long", "wide"),
    covariates_cols = NULL,
    fill_the_missing = "error",
    keep_all = "default",
    ...,
    fill_value = NULL
) {
    format <- match.arg(format)
    missing_policy <- .pb_normalize_missing_policy(
        fill_the_missing,
        fill_value = fill_value,
        argument = "fill_the_missing"
    )
    fill_the_missing <- missing_policy$policy
    fill_value <- missing_policy$fill_value

    if (identical(format, "wide")) {
        if (!is.matrix(x)) {
            stop("format='wide' requires a numeric matrix.")
        }
        corrected_matrix <- .removeBatchEffect_matrix_step(
            data_matrix = x,
            sample_annotation = sample_annotation,
            batch_col = batch_col,
            sample_id_col = sample_id_col,
            covariates_cols = covariates_cols,
            fill_the_missing = fill_the_missing,
            fill_value = fill_value,
            ...
        )
        return(corrected_matrix)
    }

    # LONG
    warning_message <- paste0(
        "removeBatchEffect cannot operate with missing values under ",
        "fill_the_missing = \"error\"; use \"keep\" to pass NA to limma, ",
        "or choose \"drop_features\"/\"fill\". The design matrix ",
        "(batch/covariates) must be free of NA in every case."
    )

    prep <- .pb_prepare_long_matrix(
        df_long = x,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        batch_col = batch_col,
        fill_the_missing = fill_the_missing,
        fill_value = fill_value,
        warning_message = warning_message,
        check_samples = FALSE,
        error_message = "format='long' requires a data.frame."
    )
    df_long <- prep$df_long
    sample_annotation <- prep$sample_annotation
    data_matrix <- prep$data_matrix
    original_cols <- prep$original_cols

    corrected_matrix <- .removeBatchEffect_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        batch_col = batch_col,
        sample_id_col = sample_id_col,
        covariates_cols = covariates_cols,
        fill_the_missing = fill_the_missing,
        fill_value = fill_value,
        ...
    )

    .post_correction_to_long(
        corrected_matrix = corrected_matrix,
        df_long = df_long,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col,
        original_cols = original_cols,
        keep_all = keep_all,
        join_method = "merge"
    )
}


# ---- Orchestrator (unified) --------------------------------------------------

#' @title Combine continuous and discrete corrections (unified)
#' @description Optional continuous drift removal + discrete adjustment
#' via \code{"MedianCentering"}, \code{"MeanCentering"}, \code{"ComBat"}, or
#' \code{"removeBatchEffect"}. Works for long or wide via \code{format}.
#' @inheritParams correct_with_ComBat
#' @param continuous_func e.g. \code{"loess_regression"} or \code{NULL}.
#' @param discrete_func batch method name.
#' @param ... Further arguments passed to the continuous fitter when enabled
#'   and to \code{limma::removeBatchEffect()} when that discrete method is
#'   selected. They are not forwarded to centering or ComBat.
#' @export
correct_batch_effects <- function(
    x,
    sample_annotation,
    format = c("long", "wide"),
    continuous_func = NULL,
    discrete_func = c(
        "MedianCentering",
        "MeanCentering",
        "ComBat",
        "removeBatchEffect"
    ),
    batch_col = "MS_batch",
    feature_id_col = "peptide_group_label",
    sample_id_col = "FullRunName",
    measure_col = "Intensity",
    order_col = "order",
    keep_all = "default",
    no_fit_imputed = TRUE,
    qual_col = NULL,
    qual_value = NULL,
    fill_the_missing = "error",
    par.prior = TRUE,
    covariates_cols = NULL,
    min_measurements = 8,
    ...,
    fill_value = NULL
) {
    format <- match.arg(format)
    discrete_func <- match.arg(discrete_func)
    missing_policy <- .pb_normalize_missing_policy(
        fill_the_missing,
        fill_value = fill_value,
        argument = "fill_the_missing"
    )
    fill_the_missing <- missing_policy$policy
    fill_value <- missing_policy$fill_value
    input_feature_ids <- NULL
    input_sample_ids <- NULL

    # Standardize to LONG for the pipeline, then back if needed
    if (identical(format, "wide")) {
        stopifnot(is.matrix(x))
        input_feature_ids <- rownames(x)
        input_sample_ids <- colnames(x)
        df_long <- matrix_to_long(
            data_matrix = x,
            feature_id_col = feature_id_col,
            measure_col = measure_col,
            sample_id_col = sample_id_col
        )
        back_to_wide <- TRUE
    } else {
        stopifnot(is.data.frame(x))
        df_long <- x
        back_to_wide <- FALSE
    }

    original_cols <- names(df_long)

    # Pre-handle missingness if requested (trims df_long/SA consistently)
    handled <- .handle_missing_for_batch_df(
        df_long = df_long,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col,
        fill_the_missing = fill_the_missing,
        fill_value = fill_value,
        warning_message = "Batch correction cannot operate with missing values in the matrix",
        qual_col = if (no_fit_imputed) qual_col else NULL,
        qual_value = if (no_fit_imputed) qual_value else NULL
    )
    df_long <- handled$df_long
    sample_annotation <- handled$sample_annotation
    original_measurements <- df_long[,
        c(
            feature_id_col,
            sample_id_col,
            measure_col
        ),
        drop = FALSE
    ]
    if (!is.null(handled$data_matrix)) {
        df_long <- .pb_apply_matrix_to_long(
            df_long = df_long,
            data_matrix = handled$data_matrix,
            feature_id_col = feature_id_col,
            sample_id_col = sample_id_col,
            measure_col = measure_col
        )
    }

    # Optional continuous drift removal (always long)
    if (!is.null(continuous_func)) {
        df_long <- adjust_batch_trend_df(
            df_long = df_long,
            sample_annotation = sample_annotation,
            batch_col = batch_col,
            feature_id_col = feature_id_col,
            sample_id_col = sample_id_col,
            measure_col = measure_col,
            order_col = order_col,
            keep_all = keep_all,
            no_fit_imputed = no_fit_imputed,
            qual_col = qual_col,
            qual_value = qual_value,
            fit_func = continuous_func,
            min_measurements = min_measurements,
            ...
        )
    }

    # Discrete registry (use unified APIs)
    df_long <- switch(
        discrete_func,
        MedianCentering = center_feature_batch(
            x = df_long,
            sample_annotation = sample_annotation,
            format = "long",
            stat = "medians",
            sample_id_col = sample_id_col,
            batch_col = batch_col,
            feature_id_col = feature_id_col,
            measure_col = measure_col,
            keep_all = keep_all,
            no_fit_imputed = no_fit_imputed,
            qual_col = qual_col,
            qual_value = qual_value
        ),
        MeanCentering = center_feature_batch(
            x = df_long,
            sample_annotation = sample_annotation,
            format = "long",
            stat = "means",
            sample_id_col = sample_id_col,
            batch_col = batch_col,
            feature_id_col = feature_id_col,
            measure_col = measure_col,
            keep_all = keep_all,
            no_fit_imputed = no_fit_imputed,
            qual_col = qual_col,
            qual_value = qual_value
        ),
        ComBat = correct_with_ComBat(
            x = df_long,
            sample_annotation = sample_annotation,
            feature_id_col = feature_id_col,
            measure_col = measure_col,
            sample_id_col = sample_id_col,
            batch_col = batch_col,
            format = "long",
            par.prior = par.prior,
            covariates_cols = covariates_cols,
            fill_the_missing = fill_the_missing,
            fill_value = fill_value,
            keep_all = keep_all,
            no_fit_imputed = no_fit_imputed,
            qual_col = qual_col,
            qual_value = qual_value
        ),
        removeBatchEffect = correct_with_removeBatchEffect(
            x = df_long,
            sample_annotation = sample_annotation,
            feature_id_col = feature_id_col,
            measure_col = measure_col,
            sample_id_col = sample_id_col,
            batch_col = batch_col,
            format = "long",
            covariates_cols = covariates_cols,
            fill_the_missing = fill_the_missing,
            fill_value = fill_value,
            keep_all = keep_all,
            ...
        )
    )

    df_long <- .pb_restore_pre_batch_values(
        df_long = df_long,
        original_measurements = original_measurements,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col
    )

    if (back_to_wide) {
        corrected_matrix <- long_to_matrix(
            df_long,
            feature_id_col = feature_id_col,
            measure_col = measure_col,
            sample_id_col = sample_id_col
        )

        if (
            !is.null(input_sample_ids) &&
                !is.null(colnames(corrected_matrix)) &&
                setequal(colnames(corrected_matrix), input_sample_ids)
        ) {
            corrected_matrix <-
                corrected_matrix[, input_sample_ids, drop = FALSE]
        }

        if (
            !is.null(input_feature_ids) &&
                !is.null(rownames(corrected_matrix)) &&
                setequal(rownames(corrected_matrix), input_feature_ids)
        ) {
            corrected_matrix <-
                corrected_matrix[input_feature_ids, , drop = FALSE]
        }

        return(corrected_matrix)
    }

    # Ensure provenance columns are retained for long
    old_measure_col <- .make_pre_col("preBatchCorr", measure_col)
    default_cols <- unique(c(
        original_cols,
        batch_col,
        old_measure_col,
        if (!is.null(continuous_func)) {
            c(
                .make_pre_col("preTrendFit", measure_col),
                "fit"
            )
        }
    ))
    minimal_cols <- c(
        sample_id_col,
        feature_id_col,
        measure_col,
        old_measure_col,
        batch_col
    )

    subset_keep_cols(
        df_long,
        keep_all,
        default_cols = default_cols,
        minimal_cols = minimal_cols
    )
}

############################################################################
# Internal functions
.make_pre_col <-
    function(prefix, measure_col) paste(prefix, measure_col, sep = "_")


.align_sample_annotation <- function(
    sample_annotation,
    sample_ids,
    sample_id_col = NULL
) {
    if (is.null(sample_annotation)) {
        stop("sample_annotation must be provided for batch correction")
    }

    sample_annotation <- as.data.frame(sample_annotation)
    sample_ids <- .pb_validate_identifiers(
        sample_ids,
        "`data_matrix` sample axis"
    )

    if (!is.null(sample_id_col)) {
        if (!(sample_id_col %in% names(sample_annotation))) {
            if (!is.null(rownames(sample_annotation))) {
                annotation_ids <- .pb_validate_identifiers(
                    rownames(sample_annotation),
                    "`sample_annotation` row names"
                )
                matches <- match(sample_ids, annotation_ids)
            } else {
                stop(sprintf(
                    "Sample ID column %s is not defined in sample annotation",
                    sample_id_col
                ))
            }
        } else {
            dummy_df <- data.frame(
                temp_id = sample_ids,
                stringsAsFactors = FALSE
            )
            names(dummy_df) <- sample_id_col
            check_sample_consistency(
                sample_annotation,
                sample_id_col,
                dummy_df,
                batch_col = NULL,
                order_col = NULL,
                facet_col = NULL,
                merge = FALSE
            )
            annotation_ids <- .pb_validate_identifiers(
                sample_annotation[[sample_id_col]],
                "`sample_annotation` sample identifiers"
            )
            matches <- match(sample_ids, annotation_ids)
        }
    } else if (!is.null(rownames(sample_annotation))) {
        annotation_ids <- .pb_validate_identifiers(
            rownames(sample_annotation),
            "`sample_annotation` row names"
        )
        matches <- match(sample_ids, annotation_ids)
    } else {
        stop(
            "Either sample_id_col must be supplied or sample_annotation must have rownames"
        )
    }

    if (anyNA(matches)) {
        stop(
            "sample_annotation is missing entries for: ",
            paste(sample_ids[is.na(matches)], collapse = ", ")
        )
    }

    sample_annotation[matches, , drop = FALSE]
}

.handle_missing_for_batch_df <- function(
    df_long,
    sample_annotation,
    feature_id_col,
    sample_id_col,
    measure_col,
    fill_the_missing,
    warning_message,
    qual_col = NULL,
    qual_value = NULL,
    fill_value = NULL
) {
    missing_policy <- .pb_normalize_missing_policy(
        fill_the_missing,
        fill_value = fill_value,
        argument = "fill_the_missing"
    )
    fill_the_missing <- missing_policy$policy
    fill_value <- missing_policy$fill_value
    if (identical(fill_the_missing, "keep")) {
        return(list(
            df_long = df_long,
            sample_annotation = sample_annotation,
            data_matrix = NULL
        ))
    }

    data_matrix <- long_to_matrix(
        df_long,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col,
        qual_col = qual_col,
        qual_value = qual_value
    )

    if (!anyNA(data_matrix)) {
        return(list(
            df_long = df_long,
            sample_annotation = sample_annotation,
            data_matrix = data_matrix
        ))
    }

    data_matrix <- handle_missing_values(
        data_matrix,
        warning_message = warning_message,
        fill_the_missing = fill_the_missing,
        fill_value = fill_value
    )

    if (!nrow(data_matrix) || !ncol(data_matrix)) {
        stop(
            "No data remaining after handling missing values for batch correction"
        )
    }

    kept_features <- rownames(data_matrix)
    kept_samples <- colnames(data_matrix)

    keep_mask <- df_long[[feature_id_col]] %in%
        kept_features &
        df_long[[sample_id_col]] %in% kept_samples
    df_long <- df_long[keep_mask, , drop = FALSE]

    if (!is.null(sample_annotation)) {
        sample_annotation <- .align_sample_annotation(
            sample_annotation,
            sample_ids = kept_samples,
            sample_id_col = sample_id_col
        )
    }

    list(
        df_long = df_long,
        sample_annotation = sample_annotation,
        data_matrix = data_matrix
    )
}

.pb_apply_matrix_to_long <- function(
    df_long,
    data_matrix,
    feature_id_col,
    sample_id_col,
    measure_col
) {
    feature_idx <- match(df_long[[feature_id_col]], rownames(data_matrix))
    sample_idx <- match(df_long[[sample_id_col]], colnames(data_matrix))
    if (anyNA(feature_idx) || anyNA(sample_idx)) {
        stop(
            "Processed matrix is not aligned with long data.",
            call. = FALSE
        )
    }
    df_long[[measure_col]] <- data_matrix[cbind(feature_idx, sample_idx)]
    df_long
}

.pb_restore_pre_batch_values <- function(
    df_long,
    original_measurements,
    feature_id_col,
    sample_id_col,
    measure_col
) {
    old_measure_col <- .make_pre_col("preBatchCorr", measure_col)
    if (!(old_measure_col %in% names(df_long))) {
        stop(
            "Corrected long data has no provenance column.",
            call. = FALSE
        )
    }

    key_cols <- c(feature_id_col, sample_id_col)
    if (anyDuplicated(original_measurements[, key_cols, drop = FALSE])) {
        same_keys <- nrow(df_long) == nrow(original_measurements) &&
            all(vapply(
                key_cols,
                function(column) {
                    identical(
                        as.character(df_long[[column]]),
                        as.character(original_measurements[[column]])
                    )
                },
                logical(1)
            ))
        if (!same_keys) {
            stop(
                "Duplicate long-data keys changed order during correction.",
                call. = FALSE
            )
        }
        df_long[[old_measure_col]] <- original_measurements[[measure_col]]
        return(df_long)
    }

    original_value_col <- make.unique(c(
        names(df_long),
        ".pb_original_measure"
    ))[[ncol(df_long) + 1L]]
    present_col <- make.unique(c(
        names(df_long),
        original_value_col,
        ".pb_original_measure_present"
    ))[[ncol(df_long) + 2L]]
    lookup <- original_measurements
    names(lookup)[names(lookup) == measure_col] <- original_value_col
    lookup[[present_col]] <- TRUE

    df_long <- left_join(
        df_long,
        lookup,
        by = c(feature_id_col, sample_id_col)
    )
    if (anyNA(df_long[[present_col]])) {
        stop(
            "Corrected long data cannot be matched to its input.",
            call. = FALSE
        )
    }
    df_long[[old_measure_col]] <- df_long[[original_value_col]]
    df_long[[original_value_col]] <- NULL
    df_long[[present_col]] <- NULL
    df_long
}

# Core ComBat matrix call with optional covariates (mod)
run_ComBat_core <- function(
    sample_annotation,
    batch_col,
    data_matrix,
    par.prior,
    covariates_cols = NULL,
    ...
) {
    if (is.null(sample_annotation)) {
        stop("sample_annotation is required for ComBat correction")
    }

    sample_annotation <- as.data.frame(sample_annotation)
    if (!(batch_col %in% names(sample_annotation))) {
        stop("Batch column is not present in sample_annotation")
    }

    # ONE batch factor only (ComBat constraint). Remove empty levels before
    # constructing the engine model to avoid rank-deficient designs.
    batches <- droplevels(factor(sample_annotation[[batch_col]]))

    if (!is.null(covariates_cols) && length(covariates_cols)) {
        missing_cov <- setdiff(covariates_cols, names(sample_annotation))
        if (length(missing_cov)) {
            stop(
                "Covariates missing in sample_annotation: ",
                paste(missing_cov, collapse = ", ")
            )
        }
        covariates <- as.data.frame(
            sample_annotation[, covariates_cols, drop = FALSE]
        )
        covariates[] <- lapply(covariates, function(column_values) {
            if (is.factor(column_values)) {
                droplevels(column_values)
            } else {
                column_values
            }
        })
        mod <- model.matrix(~., data = covariates)
    } else {
        mod <- model.matrix(~1, data = sample_annotation)
    }

    ComBat(
        dat = data_matrix,
        batch = batches,
        mod = mod,
        par.prior = par.prior,
        ...
    )
}


.loess_limmaRBE_matrix_step <- function(
    data_matrix,
    sample_annotation,
    batch_col = "MS_batch",
    sample_id_col = "FullRunName",
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    order_col = "order",
    covariates_cols = NULL,
    fill_the_missing = "error",
    min_measurements = 8,
    no_fit_imputed = TRUE,
    qual_col = NULL,
    qual_value = NULL,
    ...,
    fill_value = NULL
) {
    if (
        !is.null(sample_annotation) &&
            !(sample_id_col %in% names(sample_annotation))
    ) {
        sample_ids <- rownames(sample_annotation)
        if (is.null(sample_ids)) {
            stop(
                "sample_annotation must contain column '",
                sample_id_col,
                "' or row names for loessLimmaRBE."
            )
        }
        sample_annotation[[sample_id_col]] <- sample_ids
    }

    correct_batch_effects(
        x = data_matrix,
        sample_annotation = sample_annotation,
        format = "wide",
        continuous_func = "loess_regression",
        discrete_func = "removeBatchEffect",
        batch_col = batch_col,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col,
        order_col = order_col,
        no_fit_imputed = no_fit_imputed,
        qual_col = qual_col,
        qual_value = qual_value,
        fill_the_missing = fill_the_missing,
        fill_value = fill_value,
        covariates_cols = covariates_cols,
        min_measurements = min_measurements,
        ...
    )
}

.combat_matrix_step <- function(
    data_matrix,
    sample_annotation,
    batch_col = "MS_batch",
    sample_id_col = NULL,
    par.prior = TRUE,
    fill_the_missing = "error",
    covariates_cols = NULL,
    ...,
    fill_value = NULL
) {
    combat_args <- list(...)
    .run_matrix_method(
        data_matrix,
        sample_annotation,
        sample_id_col = sample_id_col,
        fill_the_missing = fill_the_missing,
        fill_value = fill_value,
        missing_warning = "ComBat cannot operate with missing values in the matrix",
        method_fun = function(data_matrix, sample_annotation) {
            do.call(
                run_ComBat_core,
                c(
                    list(
                        sample_annotation = sample_annotation,
                        batch_col = batch_col,
                        data_matrix = data_matrix,
                        par.prior = par.prior,
                        covariates_cols = covariates_cols
                    ),
                    combat_args
                )
            )
        }
    )
}

.removeBatchEffect_matrix_step <- function(
    data_matrix,
    sample_annotation,
    batch_col = "MS_batch",
    sample_id_col = NULL,
    covariates_cols = NULL,
    fill_the_missing = "error",
    ...,
    fill_value = NULL
) {
    .run_matrix_method(
        data_matrix,
        sample_annotation,
        sample_id_col = sample_id_col,
        fill_the_missing = fill_the_missing,
        fill_value = fill_value,
        missing_warning = paste0(
            "removeBatchEffect cannot operate with missing values under ",
            "fill_the_missing = \"error\"."
        ),
        method_fun = function(data_matrix, sample_annotation) {
            if (!(batch_col %in% names(sample_annotation))) {
                stop("Batch column is not present in sample_annotation")
            }
            batches <- droplevels(as.factor(sample_annotation[[batch_col]]))

            # design matrix (covariates optional, never include batch twice)
            if (!is.null(covariates_cols)) {
                missing_cov <- setdiff(
                    covariates_cols,
                    names(sample_annotation)
                )
                if (length(missing_cov)) {
                    stop(
                        "Covariate columns missing in sample_annotation: ",
                        paste(missing_cov, collapse = ", ")
                    )
                }
                if (batch_col %in% covariates_cols) {
                    stop(
                        "`covariates_cols` must not include `batch_col` when using removeBatchEffect."
                    )
                }
                covariates <- as.data.frame(
                    sample_annotation[, covariates_cols, drop = FALSE]
                )
                covariates[] <- lapply(covariates, function(column_values) {
                    if (is.factor(column_values)) {
                        droplevels(column_values)
                    } else {
                        column_values
                    }
                })
                degenerate_cov <- names(covariates)[vapply(
                    covariates,
                    function(column_values) {
                        values <- unique(column_values[!is.na(column_values)])
                        length(values) < 2L
                    },
                    logical(1)
                )]
                if (length(degenerate_cov)) {
                    warning(
                        "Dropping covariates with <2 observed values for removeBatchEffect: ",
                        paste(degenerate_cov, collapse = ", ")
                    )
                    covariates <-
                        covariates[,
                            setdiff(names(covariates), degenerate_cov),
                            drop = FALSE
                        ]
                }
                if (ncol(covariates)) {
                    mod <- model.matrix(~., data = covariates)
                } else {
                    mod <- model.matrix(~1, data = sample_annotation)
                }
            } else {
                mod <- model.matrix(~1, data = sample_annotation)
            }

            removeBatchEffect(data_matrix, batch = batches, design = mod, ...)
        }
    )
}

.mask_imputed_measure <- function(
    df,
    measure_col,
    qual_col,
    qual_value,
    temp_suffix = "temp",
    df_long = NULL
) {
    # Backward/forward compatibility: accept either `df` or `df_long` argument
    if (is.null(df) && !is.null(df_long)) {
        df <- df_long
    } else if (!is.null(df_long) && !identical(df, df_long)) {
        stop("Provide either `df` or `df_long` (or ensure they are identical).")
    }

    if (is.null(qual_col)) {
        return(df)
    }
    if (!(qual_col %in% names(df))) {
        stop("imputed value flag column (qual_col) is not in the data frame!")
    }
    temp_measure_col <- paste0(temp_suffix, "_", measure_col)
    df[[temp_measure_col]] <- ifelse(
        df[[qual_col]] == qual_value,
        NA,
        df[[measure_col]]
    )
    attr(df, "temp_measure_col") <- temp_measure_col
    df
}

.run_matrix_method <- function(
    data_matrix,
    sample_annotation,
    sample_id_col = NULL,
    fill_the_missing = "error",
    missing_warning = "This method cannot operate with missing values in the matrix",
    method_fun,
    fill_value = NULL
) {
    # ensure numeric matrix input for downstream modeling (sva/limma)
    if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
    }
    storage.mode(data_matrix) <- "double"
    if (!is.numeric(data_matrix)) {
        stop(
            "Input must be coercible to a numeric matrix for batch correction."
        )
    }

    data_matrix <- handle_missing_values(
        data_matrix,
        warning_message = missing_warning,
        fill_the_missing = fill_the_missing,
        fill_value = fill_value
    )
    if (!nrow(data_matrix) || !ncol(data_matrix)) {
        stop(
            "No data remaining after handling missing values for batch correction"
        )
    }

    # SA alignment
    sample_annotation <- .align_sample_annotation(
        sample_annotation,
        sample_ids = colnames(data_matrix),
        sample_id_col = sample_id_col
    )

    # call method
    method_fun(data_matrix = data_matrix, sample_annotation = sample_annotation)
}

# Effective compatibility wrappers from the later C-collated source file.
# Keeping them here removes load-order ambiguity without
# retaining a duplicate file.
#' @title DEPRECATED: center_feature_batch_medians_dm
#' @md
#' @description Use [center_feature_batch()] with
#'   `format="wide", stat="medians"`.
#' @inheritParams correct_batch_effects
#' @return A numeric matrix of batch-corrected values
#'   (features \eqn{\times} samples).
#' @examples
#' # Use the unified replacement with `format = "wide"` and `stat = "medians"`.
#' args(center_feature_batch)
#' @export
center_feature_batch_medians_dm <- function(
    data_matrix,
    sample_annotation,
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity"
) {
    .Deprecated("center_feature_batch")
    center_feature_batch(
        x = data_matrix,
        sample_annotation = sample_annotation,
        format = "wide",
        stat = "medians",
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        feature_id_col = feature_id_col,
        measure_col = measure_col
    )
}

#' @title DEPRECATED: center_feature_batch_means_dm
#' @md
#' @description Use [center_feature_batch()] with `format="wide", stat="means"`.
#' @inheritParams correct_batch_effects
#' @return A numeric matrix of batch-corrected values
#'   (features \eqn{\times} samples).
#' @examples
#' # Use the unified replacement with `format = "wide"` and `stat = "means"`.
#' args(center_feature_batch)
#' @export
center_feature_batch_means_dm <- function(
    data_matrix,
    sample_annotation,
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity"
) {
    .Deprecated("center_feature_batch")
    center_feature_batch(
        x = data_matrix,
        sample_annotation = sample_annotation,
        format = "wide",
        stat = "means",
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        feature_id_col = feature_id_col,
        measure_col = measure_col
    )
}

#' @title DEPRECATED: center_feature_batch_medians_df
#' @md
#' @description Use [center_feature_batch()] with
#'   `format="long", stat="medians"`.
#' @inheritParams correct_batch_effects
#' @return A data.frame in long format with batch-corrected
#'   values in \code{measure_col} and original values preserved in
#'   \code{preBatchCorr_[measure_col]}.
#' @examples
#' # Use the unified replacement with `format = "long"` and `stat = "medians"`.
#' args(center_feature_batch)
#' @export
center_feature_batch_medians_df <- function(
    df_long,
    sample_annotation = NULL,
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    keep_all = "default",
    no_fit_imputed = TRUE,
    qual_col = NULL,
    qual_value = NULL
) {
    .Deprecated("center_feature_batch")
    center_feature_batch(
        x = df_long,
        sample_annotation = sample_annotation,
        format = "long",
        stat = "medians",
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        keep_all = keep_all,
        no_fit_imputed = no_fit_imputed,
        qual_col = qual_col,
        qual_value = qual_value
    )
}

#' @title DEPRECATED: center_feature_batch_means_df
#' @md
#' @description Use [center_feature_batch()] with `format="long", stat="means"`.
#' @inheritParams correct_batch_effects
#' @return A data.frame in long format with batch-corrected
#'   values in \code{measure_col} and original values preserved in
#'   \code{preBatchCorr_[measure_col]}.
#' @examples
#' # Use the unified replacement with `format = "long"` and `stat = "means"`.
#' args(center_feature_batch)
#' @export
center_feature_batch_means_df <- function(
    df_long,
    sample_annotation = NULL,
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    keep_all = "default",
    no_fit_imputed = TRUE,
    qual_col = NULL,
    qual_value = NULL
) {
    .Deprecated("center_feature_batch")
    center_feature_batch(
        x = df_long,
        sample_annotation = sample_annotation,
        format = "long",
        stat = "means",
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        keep_all = keep_all,
        no_fit_imputed = no_fit_imputed,
        qual_col = qual_col,
        qual_value = qual_value
    )
}


#' @title DEPRECATED: correct_with_ComBat_df
#' @md
#' @description Use [correct_with_ComBat()] with `format="long"`.
#' @inheritParams correct_batch_effects
#' @param ... Further arguments passed to \code{sva::ComBat()}.
#' @return A data.frame in long format with ComBat-corrected
#'   values in \code{measure_col} and original values preserved in
#'   \code{preBatchCorr_[measure_col]}.
#' @examples
#' # Use the unified replacement with `format = "long"`.
#' args(correct_with_ComBat)
#' @export
correct_with_ComBat_df <- function(
    df_long,
    sample_annotation = NULL,
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    par.prior = TRUE,
    fill_the_missing = "error",
    no_fit_imputed = TRUE,
    qual_col = NULL,
    qual_value = NULL,
    keep_all = "default",
    covariates_cols = NULL,
    ...,
    fill_value = NULL
) {
    .Deprecated("correct_with_ComBat")
    correct_with_ComBat(
        x = df_long,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        format = "long",
        par.prior = par.prior,
        covariates_cols = covariates_cols,
        fill_the_missing = fill_the_missing,
        fill_value = fill_value,
        keep_all = keep_all,
        no_fit_imputed = no_fit_imputed,
        qual_col = qual_col,
        qual_value = qual_value,
        ...
    )
}

#' @title DEPRECATED: correct_with_ComBat_dm
#' @md
#' @description Use [correct_with_ComBat()] with `format="wide"`.
#' @inheritParams correct_batch_effects
#' @param ... Further arguments passed to \code{sva::ComBat()}.
#' @return A numeric matrix of ComBat-corrected values
#'   (features \eqn{\times} samples).
#' @examples
#' # Use the unified replacement with `format = "wide"`.
#' args(correct_with_ComBat)
#' @export
correct_with_ComBat_dm <- function(
    data_matrix,
    sample_annotation = NULL,
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    par.prior = TRUE,
    fill_the_missing = "error",
    covariates_cols = NULL,
    ...,
    fill_value = NULL
) {
    .Deprecated("correct_with_ComBat")
    correct_with_ComBat(
        x = data_matrix,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        format = "wide",
        par.prior = par.prior,
        covariates_cols = covariates_cols,
        fill_the_missing = fill_the_missing,
        fill_value = fill_value,
        ...
    )
}

#' @title DEPRECATED: correct_batch_effects_df
#' @md
#' @description Use [correct_batch_effects()] with `format="long"`.
#' @inheritParams correct_batch_effects
#' @param ... Further arguments forwarded to [correct_batch_effects()].
#' @return A data.frame in long format with batch-corrected
#'   values in \code{measure_col}.
#' @examples
#' # Use the unified replacement with `format = "long"`.
#' args(correct_batch_effects)
#' @export
correct_batch_effects_df <- function(
    df_long,
    sample_annotation,
    continuous_func = NULL,
    discrete_func = c(
        "MedianCentering",
        "MeanCentering",
        "ComBat",
        "removeBatchEffect"
    ),
    batch_col = "MS_batch",
    feature_id_col = "peptide_group_label",
    sample_id_col = "FullRunName",
    measure_col = "Intensity",
    order_col = "order",
    keep_all = "default",
    no_fit_imputed = TRUE,
    qual_col = NULL,
    qual_value = NULL,
    fill_the_missing = "error",
    par.prior = TRUE,
    covariates_cols = NULL,
    min_measurements = 8,
    ...,
    fill_value = NULL
) {
    .Deprecated("correct_batch_effects")
    correct_batch_effects(
        x = df_long,
        sample_annotation = sample_annotation,
        format = "long",
        continuous_func = continuous_func,
        discrete_func = match.arg(discrete_func),
        batch_col = batch_col,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col,
        order_col = order_col,
        keep_all = keep_all,
        no_fit_imputed = no_fit_imputed,
        qual_col = qual_col,
        qual_value = qual_value,
        fill_the_missing = fill_the_missing,
        fill_value = fill_value,
        par.prior = par.prior,
        covariates_cols = covariates_cols,
        min_measurements = min_measurements,
        ...
    )
}

#' @title DEPRECATED: correct_batch_effects_dm
#' @md
#' @description Use [correct_batch_effects()] with `format="wide"`.
#' @inheritParams correct_batch_effects
#' @param ... Further arguments forwarded to [correct_batch_effects()].
#' @return A numeric matrix of batch-corrected values
#'   (features \eqn{\times} samples).
#' @examples
#' # Use the unified replacement with `format = "wide"`.
#' args(correct_batch_effects)
#' @export
correct_batch_effects_dm <- function(
    data_matrix,
    sample_annotation,
    continuous_func = NULL,
    discrete_func = c(
        "MedianCentering",
        "MeanCentering",
        "ComBat",
        "removeBatchEffect"
    ),
    batch_col = "MS_batch",
    feature_id_col = "peptide_group_label",
    sample_id_col = "FullRunName",
    measure_col = "Intensity",
    order_col = "order",
    min_measurements = 8,
    no_fit_imputed = TRUE,
    fill_the_missing = "error",
    par.prior = TRUE,
    covariates_cols = NULL,
    ...,
    fill_value = NULL
) {
    .Deprecated("correct_batch_effects")
    correct_batch_effects(
        x = data_matrix,
        sample_annotation = sample_annotation,
        format = "wide",
        continuous_func = continuous_func,
        discrete_func = match.arg(discrete_func),
        batch_col = batch_col,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col,
        order_col = order_col,
        min_measurements = min_measurements,
        no_fit_imputed = no_fit_imputed,
        fill_the_missing = fill_the_missing,
        fill_value = fill_value,
        par.prior = par.prior,
        covariates_cols = covariates_cols,
        ...
    )
}

#' @title DEPRECATED: correct_with_removeBatchEffect_df
#' @md
#' @description Use [correct_with_removeBatchEffect()] with `format="long"`.
#' @inheritParams correct_batch_effects
#' @param ... Further arguments passed to
#'   \code{limma::removeBatchEffect()}.
#' @return A data.frame in long format with batch effects
#'   removed in \code{measure_col}.
#' @examples
#' # Use the unified replacement with `format = "long"`.
#' args(correct_with_removeBatchEffect)
#' @export
correct_with_removeBatchEffect_df <- function(
    df_long,
    sample_annotation = NULL,
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    covariates_cols = NULL,
    fill_the_missing = "error",
    keep_all = "default",
    ...,
    fill_value = NULL
) {
    .Deprecated("correct_with_removeBatchEffect")
    correct_with_removeBatchEffect(
        x = df_long,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        format = "long",
        covariates_cols = covariates_cols,
        fill_the_missing = fill_the_missing,
        fill_value = fill_value,
        keep_all = keep_all,
        ...
    )
}

#' @title DEPRECATED: correct_with_removeBatchEffect_dm
#' @md
#' @description Use [correct_with_removeBatchEffect()] with `format="wide"`.
#' @inheritParams correct_batch_effects
#' @param ... Further arguments passed to \code{limma::removeBatchEffect()}.
#' @return A numeric matrix with batch effects removed
#'   (features \eqn{\times} samples).
#' @examples
#' data(
#'     list = c("example_sample_annotation", "example_proteome_matrix"),
#'     package = "proBatch"
#' )
#' example_proteome_small <- example_proteome_matrix[1:100, ]
#' batch_corrected_matrix <- correct_with_removeBatchEffect(
#'     example_proteome_small,
#'     example_sample_annotation,
#'     format = "wide",
#'     batch_col = "MS_batch",
#'     covariates_cols = c("Diet", "Sex"),
#'     fill_the_missing = "drop_features"
#' )
#' @seealso \code{\link[limma:removeBatchEffect]{removeBatchEffect}}
#' @export
correct_with_removeBatchEffect_dm <- function(
    data_matrix,
    sample_annotation,
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    covariates_cols = NULL,
    fill_the_missing = "error",
    ...,
    fill_value = NULL
) {
    .Deprecated("correct_with_removeBatchEffect")
    correct_with_removeBatchEffect(
        x = data_matrix,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        format = "wide",
        covariates_cols = covariates_cols,
        fill_the_missing = fill_the_missing,
        fill_value = fill_value,
        ...
    )
}
