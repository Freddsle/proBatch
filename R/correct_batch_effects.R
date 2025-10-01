#' @title Batch correction of normalized data
#'
#' @description Batch correction of normalized data. Batch correction
#' brings each feature in each batch to the comparable shape.
#' Currently the following batch correction functions are implemented:
#' \enumerate{
#'   \item Per-feature median centering:
#'   \code{center_feature_batch_medians_df()}.
#'   Median centering of the features (per batch median).
#'   \item correction with ComBat:  \code{correct_with_ComBat_df()}.
#' Adjusts for discrete batch effects using ComBat. ComBat, described in
#' Johnson et al. 2007. It uses either parametric or
#' non-parametric empirical Bayes frameworks for adjusting data for batch
#' effects. Users are returned an expression matrix that has been corrected for
#' batch effects. The input data are assumed to be free of missing values
#' and normalized before batch effect removal. Please note that missing values
#' are common in proteomics, which is why in some cases corrections like
#' \code{center_peptide_batch_medians_df} are more appropriate.
#'   \item Continuous drift correction:  \code{adjust_batch_trend_df()}.
#' Adjust batch signal trend with the custom (continuous) fit.
#' Should be followed by discrete corrections,
#' e.g. \code{center_feature_batch_medians_df()} or
#' \code{correct_with_ComBat_df()}.
#' }
#' Alternatively, one can call the correction function with
#' \code{correct_batch_effects_df()} wrapper.
#' Batch correction method allows correction of
#' continuous signal drift within batch (if required) and adjustment for
#' discrete difference across batches.
#'
#'
#' @inheritParams proBatch
#' @param return_fit_df (logical) whether to return the \code{fit_df} from
#' \code{adjust_batch_trend_dm} or only the data matrix
#' @param fit_func function to fit the (non)-linear trend
#' @param min_measurements the number of samples in a batch required for curve fitting.
#' @param par.prior use parametrical or non-parametrical prior
#' @param continuous_func function to use for the fit (currently
#' only \code{loess_regression} available); if order-associated fix is not
#' required, should be \code{NULL}.
#' @param discrete_func function to use for adjustment of discrete batch effects
#' (\code{MedianCentering} or \code{ComBat}).
#' @param fill_the_missing numeric value used to impute missing measurements before
#' correction. If \code{NULL} (default) missing values are left untouched, if
#' \code{FALSE} rows with missing values are removed prior to correction.
#' @param ... other parameters, usually of \code{adjust_batch_trend},
#' and \code{fit_func}.

#'
#' @return the data in the same format as input (\code{data_matrix} or
#' \code{df_long}).
#' For \code{df_long} the data frame stores the original values of
#' \code{measure_col}
#' in another column called "preBatchCorr_[measure_col]", and the normalized
#' values in \code{measure_col} column.
#'
#' The function \code{adjust_batch_trend_dm()}, if \code{return_fit_df} is
#' \code{TRUE} returns list of two items:
#' \enumerate{
#'   \item \code{data_matrix}
#'   \item \code{fit_df}, used to examine the fitting curves
#' }

#'
#' @examples
#' # Load necessary datasets
#' data(
#'     list = c("example_sample_annotation", "example_proteome"),
#'     package = "proBatch"
#' )
#'
#' # Median centering per feature per batch:
#' median_centered_df <- center_feature_batch_medians_df(
#'     example_proteome, example_sample_annotation
#' )
#'
#' # Correct with ComBat:
#' combat_corrected_df <- correct_with_ComBat_df(
#'     example_proteome,
#'     example_sample_annotation
#' )
#'
#' # Adjust the MS signal drift:
#' test_peptides <- unique(example_proteome$peptide_group_label)[1:3]
#' test_peptide_filter <- example_proteome$peptide_group_label %in% test_peptides
#' test_proteome <- example_proteome[test_peptide_filter, ]
#' adjusted_df <- adjust_batch_trend_df(test_proteome,
#'     example_sample_annotation,
#'     span = 0.7,
#'     min_measurements = 8
#' )
#' plot_fit <- plot_with_fitting_curve(unique(adjusted_df$peptide_group_label),
#'     df_long = adjusted_df, measure_col = "preTrendFit_Intensity",
#'     fit_df = adjusted_df, sample_annotation = example_sample_annotation
#' )
#'
#' # Correct the data in one go:
#' batch_corrected_matrix <- correct_batch_effects_df(example_proteome,
#'     example_sample_annotation,
#'     continuous_func = "loess_regression",
#'     discrete_func = "MedianCentering",
#'     batch_col = "MS_batch",
#'     span = 0.7, min_measurements = 8
#' )
#'
#' @seealso \code{\link{fit_nonlinear}}
#' @name correct_batch_effects
NULL

#'
#' @export
#' @rdname correct_batch_effects
#'
center_feature_batch_medians_dm <- function(data_matrix, sample_annotation,
                                            sample_id_col = "FullRunName",
                                            batch_col = "MS_batch",
                                            feature_id_col = "peptide_group_label",
                                            measure_col = "Intensity") {
    .center_feature_batch_stat_dm(
        data_matrix, sample_annotation,
        sample_id_col, batch_col, feature_id_col, measure_col,
        df_stat_fun = center_feature_batch_medians_df
    )
}

#' @export
#' @rdname correct_batch_effects
#'
center_feature_batch_means_dm <- function(data_matrix, sample_annotation,
                                          sample_id_col = "FullRunName",
                                          batch_col = "MS_batch",
                                          feature_id_col = "peptide_group_label",
                                          measure_col = "Intensity") {
    .center_feature_batch_stat_dm(
        data_matrix, sample_annotation,
        sample_id_col, batch_col, feature_id_col, measure_col,
        df_stat_fun = center_feature_batch_means_df
    )
}

.center_feature_batch_stat_dm <- function(data_matrix, sample_annotation,
                                          sample_id_col, batch_col,
                                          feature_id_col, measure_col,
                                          df_stat_fun) {
    df_long <- matrix_to_long(
        data_matrix,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col
    )
    corrected_df <- df_stat_fun(
        df_long, sample_annotation,
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        feature_id_col = feature_id_col,
        measure_col = measure_col
    )
    long_to_matrix(
        corrected_df,
        sample_id_col = sample_id_col,
        measure_col = measure_col,
        feature_id_col = feature_id_col
    )
}

#'
#' @export
#' @rdname correct_batch_effects
#'
center_feature_batch_means_df <- function(df_long, sample_annotation = NULL,
                                          sample_id_col = "FullRunName",
                                          batch_col = "MS_batch",
                                          feature_id_col = "peptide_group_label",
                                          measure_col = "Intensity",
                                          keep_all = "default",
                                          no_fit_imputed = TRUE,
                                          qual_col = NULL,
                                          qual_value = NULL) {
    .center_feature_batch_stat_df(
        df_long, sample_annotation,
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        keep_all = keep_all,
        stat = "mean",
        no_fit_imputed = no_fit_imputed,
        qual_col = qual_col,
        qual_value = qual_value,
        stat_names = c(batch = "mean_batch", global = "mean_global", diff = "diff_means")
    )
}

#'
#' @export
#' @rdname correct_batch_effects
#'
center_feature_batch_medians_df <- function(df_long, sample_annotation = NULL,
                                            sample_id_col = "FullRunName",
                                            batch_col = "MS_batch",
                                            feature_id_col = "peptide_group_label",
                                            measure_col = "Intensity",
                                            keep_all = "default",
                                            no_fit_imputed = TRUE,
                                            qual_col = NULL,
                                            qual_value = NULL) {
    .center_feature_batch_stat_df(
        df_long, sample_annotation,
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        keep_all = keep_all,
        stat = "median",
        no_fit_imputed = no_fit_imputed,
        qual_col = qual_col,
        qual_value = qual_value,
        stat_names = c(batch = "median_batch", global = "median_global", diff = "diff_medians")
    )
}

#' @export
#' @rdname correct_batch_effects
#'
center_feature_batch_means_dm <- function(data_matrix, sample_annotation,
                                          sample_id_col = "FullRunName",
                                          batch_col = "MS_batch",
                                          feature_id_col = "peptide_group_label",
                                          measure_col = "Intensity") {
    df_long <- matrix_to_long(
        data_matrix,
        feature_id_col = feature_id_col,
        measure_col = measure_col, sample_id_col = sample_id_col
    )
    corrected_df <- center_feature_batch_means_df(
        df_long, sample_annotation,
        sample_id_col = sample_id_col, batch_col = batch_col,
        feature_id_col = feature_id_col, measure_col = measure_col
    )
    long_to_matrix(corrected_df,
        sample_id_col = sample_id_col,
        measure_col = measure_col, feature_id_col = feature_id_col
    )
}

#'
#' @export
#' @rdname correct_batch_effects
#'
#' @seealso \code{\link{fit_nonlinear}}, \code{\link{plot_with_fitting_curve}}
adjust_batch_trend_df <- function(df_long, sample_annotation = NULL,
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
                                  min_measurements = 8, ...) {
    original_cols <- names(df_long)

    df_long <- check_sample_consistency(
        sample_annotation,
        sample_id_col, df_long,
        batch_col,
        order_col = order_col,
        facet_col = NULL,
        merge = TRUE
    )

    if (no_fit_imputed) {
        if (is.null(qual_col)) {
            warning("`qual_col` is NULL, setting `no_fit_imputed = FALSE` so imputed flags are ignored.")
            no_fit_imputed <- FALSE
        } else if (!(qual_col %in% names(df_long))) {
            stop("imputed value flag column (qual_col) is not in the data frame!")
        }
    } else {
        if (!is.null(qual_col)) {
            # flags provided but explicitly ignored -- keep behavior, inform user
            warning("`qual_col` provided but `no_fit_imputed = FALSE`; imputed flags will be ignored for curve fitting.")
        }
    }

    # If a batch column is requested, ensure it's present after consistency checks.
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
                        min_measurements = min_measurements, ...
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
            diff.na  = ifelse(is.na(diff_fit), 0, diff_fit)
        ) %>%
        rename(!!sym(old_measure_col) := !!sym(measure_col)) %>%
        # Conditional shift: use diff.na so NA in fit doesn't propagate the shift
        mutate(!!sym(measure_col) := diff.na + .data[[old_measure_col]]) %>%
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

#'
#' @export
#' @rdname correct_batch_effects
#'
#' @seealso \code{\link{fit_nonlinear}}, \code{\link{plot_with_fitting_curve}}
adjust_batch_trend_dm <- function(data_matrix, sample_annotation,
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
                                  min_measurements = 8, ...) {
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
        min_measurements = min_measurements
    )

    corrected_df <- corrected_data
    corrected_dm <- long_to_matrix(
        corrected_df,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col
    )
    if (return_fit_df) {
        # extract only the columns relevant for inspecting the fit (only non-empty columns)
        fit_columns <- c(sample_id_col, feature_id_col, batch_col, order_col, "fit")
        # if any of the fit columns are not present in the corrected_df, remove them and warn
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
correct_with_ComBat_df <- function(df_long, sample_annotation = NULL,
                                   feature_id_col = "peptide_group_label",
                                   measure_col = "Intensity",
                                   sample_id_col = "FullRunName",
                                   batch_col = "MS_batch",
                                   par.prior = TRUE,
                                   fill_the_missing = NULL,
                                   no_fit_imputed = TRUE,
                                   qual_col = NULL,
                                   qual_value = NULL,
                                   keep_all = "default") {
    original_cols <- names(df_long)
    df_long <- check_sample_consistency(
        sample_annotation,
        sample_id_col,
        df_long,
        batch_col,
        order_col = NULL,
        facet_col = NULL,
        merge = FALSE
    )

    if (!is.null(qual_col) && no_fit_imputed) {
        if (!(qual_col %in% names(df_long))) {
            stop("imputed value flag column is not in the data frame!")
        }
        data_matrix <- long_to_matrix(
            df_long,
            feature_id_col = feature_id_col,
            measure_col = measure_col,
            sample_id_col = sample_id_col,
            qual_col = qual_col,
            qual_value = qual_value
        )
    } else {
        data_matrix <- long_to_matrix(
            df_long,
            feature_id_col = feature_id_col,
            measure_col = measure_col,
            sample_id_col = sample_id_col,
            qual_col = NULL
        )
    }

    handled <- .handle_missing_for_batch_df(
        df_long = df_long,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col,
        fill_the_missing = fill_the_missing,
        warning_message = "ComBat cannot operate with missing values in the matrix",
        qual_col = if (no_fit_imputed) qual_col else NULL,
        qual_value = if (no_fit_imputed) qual_value else NULL
    )
    df_long <- handled$df_long
    sample_annotation <- handled$sample_annotation

    corrected_matrix <- .combat_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        batch_col = batch_col,
        sample_id_col = sample_id_col,
        par.prior = par.prior,
        fill_the_missing = fill_the_missing
    )

    corrected_df <- matrix_to_long(
        corrected_matrix,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col
    )

    old_measure_col <- paste("preBatchCorr", measure_col, sep = "_")

    df_long <- df_long %>%
        rename(!!sym(old_measure_col) := !!sym(measure_col))

    corrected_df <- corrected_df %>%
        merge(df_long, by = c(feature_id_col, sample_id_col))

    default_cols <- c(original_cols, old_measure_col)
    minimal_cols <- c(sample_id_col, feature_id_col, measure_col, old_measure_col)

    corrected_df <- subset_keep_cols(
        corrected_df,
        keep_all,
        default_cols = default_cols,
        minimal_cols = minimal_cols
    )

    return(corrected_df)
}

#'
#' @export
#' @rdname correct_batch_effects
#'
correct_with_ComBat_dm <- function(data_matrix, sample_annotation = NULL,
                                   feature_id_col = "peptide_group_label",
                                   measure_col = "Intensity",
                                   sample_id_col = "FullRunName",
                                   batch_col = "MS_batch",
                                   par.prior = TRUE,
                                   fill_the_missing = NULL) {
    corrected_matrix <- .combat_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        batch_col = batch_col,
        sample_id_col = sample_id_col,
        par.prior = par.prior,
        fill_the_missing = fill_the_missing
    )

    return(corrected_matrix)
}

run_ComBat_core <- function(sample_annotation, batch_col, data_matrix,
                            par.prior, ...) {
    # TODO: handle the case of multiple batch factors
    if (is.null(sample_annotation)) {
        stop("sample_annotation is required for ComBat correction")
    }
    sample_annotation <- as.data.frame(sample_annotation)
    if (!(batch_col %in% names(sample_annotation))) {
        stop("Batch column is not present in sample_annotation")
    }
    # Coerce to factor
    batches <- as.factor(sample_annotation[[batch_col]])
    modCombat <- model.matrix(~1, data = sample_annotation)
    corrected_matrix <- ComBat(
        dat = data_matrix, batch = batches,
        mod = modCombat, par.prior = par.prior, ...
    )
    return(corrected_matrix)
}


#'
#' @export
#' @rdname correct_batch_effects
#'
correct_batch_effects_df <- function(df_long, sample_annotation,
                                     continuous_func = NULL,
                                     discrete_func = c("MedianCentering", "MeanCentering", "ComBat", "removeBatchEffect"),
                                     batch_col = "MS_batch",
                                     feature_id_col = "peptide_group_label",
                                     sample_id_col = "FullRunName",
                                     measure_col = "Intensity",
                                     order_col = "order",
                                     keep_all = "default",
                                     no_fit_imputed = TRUE,
                                     qual_col = NULL,
                                     qual_value = NULL,
                                     fill_the_missing = NULL,
                                     min_measurements = 8, ...) {
    discrete_func <- match.arg(discrete_func)
    original_cols <- names(df_long)

    handled <- .handle_missing_for_batch_df(
        df_long = df_long,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col,
        fill_the_missing = fill_the_missing,
        warning_message = "Batch correction cannot operate with missing values in the matrix",
        qual_col = if (no_fit_imputed) qual_col else NULL,
        qual_value = if (no_fit_imputed) qual_value else NULL
    )
    df_long <- handled$df_long
    sample_annotation <- handled$sample_annotation

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
            min_measurements = min_measurements, ...
        )
    }

    # registry
    discrete_methods <- list(
        MedianCentering = function() {
            center_feature_batch_medians_df(
                df_long, sample_annotation, sample_id_col, batch_col,
                feature_id_col, measure_col, keep_all, no_fit_imputed, qual_col, qual_value
            )
        },
        MeanCentering = function() {
            center_feature_batch_means_df(
                df_long, sample_annotation, sample_id_col, batch_col,
                feature_id_col, measure_col, keep_all, no_fit_imputed, qual_col, qual_value
            )
        },
        ComBat = function() {
            correct_with_ComBat_df(
                df_long, sample_annotation, feature_id_col, measure_col,
                sample_id_col, batch_col,
                par.prior = TRUE, fill_the_missing = fill_the_missing
            )
        },
        removeBatchEffect = function() {
            correct_with_removeBatchEffect_df(
                df_long, sample_annotation, feature_id_col, measure_col,
                sample_id_col, batch_col,
                covariates_cols = NULL,
                fill_the_missing = fill_the_missing, keep_all = keep_all, ...
            )
        }
    )

    corrected_df <- discrete_methods[[discrete_func]]()

    old_measure_col <- .make_pre_col("preBatchCorr", measure_col)
    if (!is.null(continuous_func)) {
        preFit_measure_col <- .make_pre_col("preTrendFit", measure_col)
        default_cols <- c(original_cols, batch_col, old_measure_col, preFit_measure_col, "fit")
    } else {
        default_cols <- c(original_cols, batch_col, old_measure_col)
    }
    minimal_cols <- c(sample_id_col, feature_id_col, measure_col, old_measure_col, batch_col)

    subset_keep_cols(corrected_df, keep_all,
        default_cols = default_cols,
        minimal_cols = minimal_cols
    )
}

#'
#' @export
#' @rdname correct_batch_effects
#'
correct_batch_effects_dm <- function(data_matrix, sample_annotation,
                                     continuous_func = NULL,
                                     discrete_func = c(
                                         "MedianCentering",
                                         "MeanCentering",
                                         "ComBat"
                                     ),
                                     batch_col = "MS_batch",
                                     feature_id_col = "peptide_group_label",
                                     sample_id_col = "FullRunName",
                                     measure_col = "Intensity",
                                     order_col = "order",
                                     min_measurements = 8,
                                     no_fit_imputed = TRUE,
                                     fill_the_missing = NULL,
                                     ...) {
    df_long <- matrix_to_long(
        data_matrix,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col
    )

    corrected_df <- correct_batch_effects_df(
        df_long,
        sample_annotation,
        continuous_func = continuous_func,
        discrete_func = discrete_func,
        batch_col = batch_col,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col,
        order_col = order_col,
        min_measurements = min_measurements,
        no_fit_imputed = no_fit_imputed,
        fill_the_missing = fill_the_missing,
        qual_col = NULL,
        qual_value = NULL,
        keep_all = "default", ...
    )
    # Convert the corrected data frame back to matrix format
    corrected_matrix <- long_to_matrix(
        corrected_df,
        sample_id_col = sample_id_col,
        measure_col = measure_col,
        feature_id_col = feature_id_col
    )

    return(corrected_matrix)
}

#' @export
#' @rdname correct_batch_effects
correct_with_removeBatchEffect_df <- function(df_long, sample_annotation = NULL,
                                              feature_id_col = "peptide_group_label",
                                              measure_col = "Intensity",
                                              sample_id_col = "FullRunName",
                                              batch_col = "MS_batch",
                                              covariates_cols = NULL,
                                              fill_the_missing = NULL,
                                              keep_all = "default", ...) {
    original_cols <- names(df_long)

    # long -> matrix
    data_matrix <- long_to_matrix(
        df_long,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col
    )

    # run limma path (reuses your existing matrix step)
    corrected_matrix <- .removeBatchEffect_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        batch_col = batch_col,
        sample_id_col = sample_id_col,
        covariates_cols = covariates_cols,
        fill_the_missing = fill_the_missing, ...
    )

    # matrix -> long and merge provenance
    corrected_df <- matrix_to_long(
        corrected_matrix,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col
    )
    old_measure_col <- .make_pre_col("preBatchCorr", measure_col)
    df_long <- dplyr::rename(df_long, !!old_measure_col := !!measure_col)
    corrected_df <- merge(corrected_df, df_long, by = c(feature_id_col, sample_id_col))

    default_cols <- c(original_cols, old_measure_col)
    minimal_cols <- c(sample_id_col, feature_id_col, measure_col, old_measure_col)

    subset_keep_cols(corrected_df, keep_all,
        default_cols = default_cols,
        minimal_cols = minimal_cols
    )
}


#' @title Batch effect correction with removeBatchEffect from limma
#' @description Batch effect correction with removeBatchEffect.
#' @param data_matrix data matrix with features in rows and samples in columns
#' @param sample_annotation data frame with sample annotations
#' @param feature_id_col column name in \code{data_matrix} with feature IDs
#' @param measure_col column name in \code{data_matrix} with measured values
#' @param sample_id_col column name in \code{sample_annotation} with sample IDs
#' @param batch_col column name in \code{sample_annotation} with batch IDs
#' @param covariates_cols vector of column names in \code{sample_annotation}
#' with covariates to include in the model
#' @param fill_the_missing numeric value used to impute missing measurements
#' before correction. If \code{FALSE} rows with missing values are removed.
#' @param ... other parameters to pass to \code{removeBatchEffect}
#' @return data matrix with batch effects removed
#' @examples
#' data(
#'     list = c("example_sample_annotation", "example_proteome_matrix"),
#'     package = "proBatch"
#' )
#' batch_corrected_matrix <- correct_with_removeBatchEffect_dm(
#'     example_proteome_matrix,
#'     example_sample_annotation,
#'     batch_col = "MS_batch",
#'     covariates_cols = c("Condition", "Type")
#' )
#' @seealso \code{\link{removeBatchEffect}}
#' @export
correct_with_removeBatchEffect_dm <- function(data_matrix, sample_annotation,
                                              feature_id_col = "peptide_group_label",
                                              measure_col = "Intensity",
                                              sample_id_col = "FullRunName",
                                              batch_col = "MS_batch",
                                              covariates_cols = NULL,
                                              fill_the_missing = NULL, ...) {
    corrected_matrix <- .removeBatchEffect_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        batch_col = batch_col,
        sample_id_col = sample_id_col,
        covariates_cols = covariates_cols,
        fill_the_missing = fill_the_missing,
        ...
    )
    return(corrected_matrix)
}

############################################################################
# Internal functions
.make_pre_col <- function(prefix, measure_col) paste(prefix, measure_col, sep = "_")


.align_sample_annotation <- function(sample_annotation, sample_ids,
                                     sample_id_col = NULL) {
    if (is.null(sample_annotation)) {
        stop("sample_annotation must be provided for batch correction")
    }

    sample_annotation <- as.data.frame(sample_annotation)

    if (!is.null(sample_id_col)) {
        if (!(sample_id_col %in% names(sample_annotation))) {
            if (!is.null(rownames(sample_annotation))) {
                matches <- match(sample_ids, rownames(sample_annotation))
            } else {
                stop(sprintf(
                    "Sample ID column %s is not defined in sample annotation",
                    sample_id_col
                ))
            }
        } else {
            dummy_df <- data.frame(temp_id = sample_ids, stringsAsFactors = FALSE)
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
            matches <- match(sample_ids, sample_annotation[[sample_id_col]])
        }
    } else if (!is.null(rownames(sample_annotation))) {
        matches <- match(sample_ids, rownames(sample_annotation))
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

.handle_missing_for_batch_df <- function(df_long,
                                         sample_annotation,
                                         feature_id_col,
                                         sample_id_col,
                                         measure_col,
                                         fill_the_missing,
                                         warning_message,
                                         qual_col = NULL,
                                         qual_value = NULL) {
    handle_flag <- !is.null(fill_the_missing) || identical(fill_the_missing, FALSE)
    if (!handle_flag) {
        return(list(df_long = df_long, sample_annotation = sample_annotation))
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
        return(list(df_long = df_long, sample_annotation = sample_annotation))
    }

    data_matrix <- handle_missing_values(
        data_matrix,
        warning_message = warning_message,
        fill_the_missing = fill_the_missing
    )

    if (!nrow(data_matrix) || !ncol(data_matrix)) {
        stop("No data remaining after handling missing values for batch correction")
    }

    kept_features <- rownames(data_matrix)
    kept_samples <- colnames(data_matrix)

    keep_mask <- df_long[[feature_id_col]] %in% kept_features &
        df_long[[sample_id_col]] %in% kept_samples
    df_long <- df_long[keep_mask, , drop = FALSE]

    feature_idx <- match(df_long[[feature_id_col]], kept_features)
    sample_idx <- match(df_long[[sample_id_col]], kept_samples)
    df_long[[measure_col]] <- data_matrix[cbind(feature_idx, sample_idx)]

    if (!is.null(sample_annotation)) {
        sample_annotation <- .align_sample_annotation(
            sample_annotation,
            sample_ids = kept_samples,
            sample_id_col = sample_id_col
        )
    }

    list(df_long = df_long, sample_annotation = sample_annotation)
}

.combat_matrix_step <- function(data_matrix, sample_annotation,
                                batch_col = "MS_batch",
                                sample_id_col = NULL,
                                par.prior = TRUE,
                                fill_the_missing = NULL, ...) {
    .run_matrix_method(
        data_matrix, sample_annotation,
        sample_id_col = sample_id_col,
        fill_the_missing = fill_the_missing,
        missing_warning = "ComBat cannot operate with missing values in the matrix",
        method_fun = function(data_matrix, sample_annotation) {
            if (!(batch_col %in% names(sample_annotation))) {
                stop("Batch column is not present in sample_annotation")
            }
            batches <- as.factor(sample_annotation[[batch_col]])
            modCombat <- model.matrix(~1, data = sample_annotation)
            ComBat(dat = data_matrix, batch = batches, mod = modCombat, par.prior = par.prior, ...)
        }
    )
}


.removeBatchEffect_matrix_step <- function(data_matrix, sample_annotation,
                                           batch_col = "MS_batch",
                                           sample_id_col = NULL,
                                           covariates_cols = NULL,
                                           fill_the_missing = NULL, ...) {
    .run_matrix_method(
        data_matrix, sample_annotation,
        sample_id_col = sample_id_col,
        fill_the_missing = fill_the_missing,
        missing_warning = "removeBatchEffect cannot operate with missing values in the matrix",
        method_fun = function(data_matrix, sample_annotation) {
            if (!(batch_col %in% names(sample_annotation))) {
                stop("Batch column is not present in sample_annotation")
            }
            batches <- as.factor(sample_annotation[[batch_col]])

            # design matrix (covariates optional, never include batch twice)
            if (!is.null(covariates_cols)) {
                missing_cov <- setdiff(covariates_cols, names(sample_annotation))
                if (length(missing_cov)) {
                    stop("Covariate columns missing in sample_annotation: ", paste(missing_cov, collapse = ", "))
                }
                if (batch_col %in% covariates_cols) {
                    stop("`covariates_cols` must not include `batch_col` when using removeBatchEffect.")
                }
                covariates <- as.data.frame(sample_annotation[, covariates_cols, drop = FALSE])
                mod <- model.matrix(~., data = covariates)
            } else {
                mod <- model.matrix(~1, data = sample_annotation)
            }

            removeBatchEffect(data_matrix, batch = batches, design = mod, ...)
        }
    )
}

.mask_imputed_measure <- function(df, measure_col, qual_col, qual_value, temp_suffix = "temp") {
    if (is.null(qual_col)) {
        return(df)
    }
    if (!(qual_col %in% names(df))) {
        stop("imputed value flag column (qual_col) is not in the data frame!")
    }
    temp_measure_col <- paste0(temp_suffix, "_", measure_col)
    df[[temp_measure_col]] <- ifelse(df[[qual_col]] == qual_value, NA, df[[measure_col]])
    attr(df, "temp_measure_col") <- temp_measure_col
    df
}

.center_feature_batch_stat_df <- function(df_long, sample_annotation = NULL,
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
                                              batch = "median_batch", global = "median_global",
                                              diff = "diff_medians"
                                          )) {
    stat <- match.arg(stat)
    original_cols <- names(df_long)

    df_long <- check_sample_consistency(
        sample_annotation, sample_id_col, df_long,
        batch_col,
        order_col = NULL, facet_col = NULL, merge = TRUE
    )

    if (no_fit_imputed && is.null(qual_col)) {
        warning("`qual_col` is NULL, setting `no_fit_imputed = FALSE` so imputed flags are ignored.")
        no_fit_imputed <- FALSE
    }

    # choose summariser
    summariser <- switch(stat,
        median = stats::median,
        mean = base::mean
    )

    # Optionally mask imputed values for inference
    tmp_col <- NULL
    if (no_fit_imputed) {
        df_long <- .mask_imputed_measure(df_long, measure_col, qual_col, qual_value)
        tmp_col <- attr(df_long, "temp_measure_col")
    }

    measure_for_inference <- if (!is.null(tmp_col)) tmp_col else measure_col
    old_measure_col <- .make_pre_col("preBatchCorr", measure_col)

    # infer per-batch, per-feature location; then per-feature global location
    corrected_df <- df_long %>%
        dplyr::group_by(dplyr::across(tidyselect::any_of(c(batch_col, feature_id_col)))) %>%
        dplyr::mutate(!!stat_names["batch"] := summariser(.data[[measure_for_inference]], na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(dplyr::across(tidyselect::any_of(feature_id_col))) %>%
        dplyr::mutate(!!stat_names["global"] := summariser(.data[[measure_for_inference]], na.rm = TRUE)) %>%
        dplyr::ungroup()

    # compute shift and apply
    corrected_df <- corrected_df %>%
        dplyr::mutate(!!stat_names["diff"] := .data[[stat_names["global"]]] - .data[[stat_names["batch"]]]) %>%
        dplyr::rename(!!old_measure_col := !!measure_col) %>%
        dplyr::mutate(!!measure_col := .data[[old_measure_col]] + .data[[stat_names["diff"]]])

    # Drop internal temp column if it exists
    if (!is.null(tmp_col) && tmp_col %in% names(corrected_df)) {
        corrected_df <- dplyr::select(corrected_df, -tidyselect::all_of(tmp_col))
    }

    # keep columns
    default_cols <- unique(c(
        original_cols, batch_col, old_measure_col,
        unname(stat_names["batch"]), unname(stat_names["global"]), unname(stat_names["diff"])
    ))
    minimal_cols <- unique(c(
        sample_id_col, feature_id_col, measure_col, old_measure_col,
        batch_col, unname(stat_names["batch"]), unname(stat_names["diff"])
    ))
    if (!is.null(qual_col) && qual_col %in% names(corrected_df)) {
        default_cols <- c(default_cols, qual_col)
        minimal_cols <- c(minimal_cols, qual_col)
    }

    subset_keep_cols(corrected_df, keep_all,
        default_cols = default_cols,
        minimal_cols = minimal_cols
    )
}

.run_matrix_method <- function(data_matrix, sample_annotation,
                               sample_id_col = NULL,
                               fill_the_missing = NULL,
                               missing_warning = "This method cannot operate with missing values in the matrix",
                               method_fun) {
    # optional NA handling
    handle_flag <- !is.null(fill_the_missing) || identical(fill_the_missing, FALSE)
    if (handle_flag && anyNA(data_matrix)) {
        data_matrix <- handle_missing_values(
            data_matrix,
            warning_message = missing_warning,
            fill_the_missing = fill_the_missing
        )
        if (!nrow(data_matrix) || !ncol(data_matrix)) {
            stop("No data remaining after handling missing values for batch correction")
        }
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
