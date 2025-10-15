#' @title Batch correction methods for normalized data
#'
#' @description
#' Batch correction of normalized data to bring each feature within each batch
#' to a comparable scale. The following methods are available:
#' \enumerate{
#'   \item \strong{Per-feature centering by median/mean}:
#'     \code{\link{center_feature_batch}()} — unified API for long (\code{"long"})
#'     and wide (\code{"wide"}) data via \code{format}, and for medians/means via
#'     \code{stat}. \emph{Deprecated wrappers}:
#'     \code{center_feature_batch_medians_df()}, \code{center_feature_batch_means_df()},
#'     \code{center_feature_batch_medians_dm()}, \code{center_feature_batch_means_dm()}
#'     now forward to \code{center_feature_batch()} and emit deprecation warnings.
#'
#'   \item \strong{Discrete batch correction with ComBat}:
#'     \code{\link{correct_with_ComBat}()} adjusts for discrete batch effects
#'     using parametric or non-parametric empirical Bayes (Johnson et al., 2007).
#'     \strong{Missing data:} ComBat requires an \emph{NA-free} matrix. If your data
#'     contain missing values, either set \code{fill_the_missing = FALSE} to drop
#'     NA-containing features/columns or provide a numeric value to impute before
#'     calling ComBat.
#'
#'   \item \strong{Linear batch correction with limma}:
#'     \code{\link{correct_with_removeBatchEffect}()} removes linear batch effects
#'     via \code{limma::removeBatchEffect}. \strong{Missing data:} NA values in the
#'     \emph{expression/intensity matrix} are allowed; the \emph{design matrix}
#'     (batch/covariates) must be free of missing values.
#'
#'   \item \strong{Continuous drift correction}:
#'     \code{\link{adjust_batch_trend_df}()} fits and removes within-batch trends
#'     (e.g., LOESS) and is typically followed by a discrete adjustment such as
#'     \code{center_feature_batch()}, \code{correct_with_ComBat()}, or
#'     \code{correct_with_removeBatchEffect()}.
#'
#'   \item \strong{Tree-based batch correction with BERT}:
#'     \code{\link{correct_with_BERT}()} performs batch-effect correction
#'     that \emph{tolerates missing values} (no pre-imputation
#'     required) and can parallelize sub-trees. Internally calls
#'     limmaRBE or ComBat BERT-versions.
#'     (Available only if \pkg{BERT} is installed.)
#' }
#'
#' Alternatively, use the wrapper \code{\link{correct_batch_effects}()} to
#' combine continuous and/or discrete corrections in one call.
#'
#' @inheritParams proBatch
#'
#' @param return_fit_df Logical; whether \code{\link{adjust_batch_trend_dm}()}
#'   returns the \code{fit_df} (for curve inspection) alongside the corrected matrix.
#' @param fit_func Function used for trend fitting (e.g., \code{"loess_regression"}).
#' @param min_measurements Minimum number of samples per batch required for fitting.
#' @param par.prior Logical; use parametric prior (ComBat) or non-parametric.
#' @param continuous_func Which function to use for the continuous fit (currently
#'   only \code{"loess_regression"}); set \code{NULL} if order-associated drift
#'   correction is not required.
#' @param discrete_func Which function to use for discrete batch effects in the
#'   wrapper: one of \code{"MedianCentering"}, \code{"MeanCentering"},
#'   \code{"ComBat"}, or \code{"removeBatchEffect"}.
#' @param fill_the_missing Missing-value policy applied \emph{before} discrete correction.
#'   If \code{NULL} (default), missing values are left as is (no imputation/dropping).
#'   For \code{ComBat}, the matrix must be NA-free — set \code{FALSE} to drop rows with NA
#'   or provide a numeric value to impute. For \code{removeBatchEffect}, NA values in the data
#'   matrix are permitted; however, the design (batch/covariates) must not contain NA.
#' @param ... Additional parameters passed to \code{adjust_batch_trend_df()} and
#'   the chosen \code{fit_func}.
#'
#' @return
#' Returns data in the same format as input (\code{data_matrix} or \code{df_long}).
#' For long format, original values from \code{measure_col} are preserved in
#' \code{"preBatchCorr_[measure_col]"} and corrected values are written to
#' \code{measure_col}.
#'
#' The function \code{\link{adjust_batch_trend_dm}()}, if \code{return_fit_df = TRUE},
#' returns a list with:
#' \enumerate{
#'   \item \code{corrected_dm} — corrected data matrix
#'   \item \code{fit_df} — data frame to inspect fitted curves
#' }
#'
#' @examples
#' # Load example data
#' data(
#'     list = c("example_sample_annotation", "example_proteome"),
#'     package = "proBatch"
#' )
#'
#' # 1) Per-feature centering (LONG): medians
#' median_centered_df <- center_feature_batch(
#'     x = example_proteome,
#'     sample_annotation = example_sample_annotation,
#'     format = "long", stat = "medians",
#'     sample_id_col = "FullRunName", batch_col = "MS_batch",
#'     feature_id_col = "peptide_group_label", measure_col = "Intensity"
#' )
#'
#' # 2) Per-feature centering (WIDE): means
#' data(example_proteome_matrix, package = "proBatch")
#' mean_centered_mat <- center_feature_batch(
#'     x = example_proteome_matrix,
#'     sample_annotation = example_sample_annotation,
#'     format = "wide", stat = "means",
#'     sample_id_col = "FullRunName", batch_col = "MS_batch",
#'     feature_id_col = "peptide_group_label", measure_col = "Intensity"
#' )
#'
#' # 3) ComBat (discrete) — drop NA features/samples first if needed
#' combat_corrected_df <- correct_with_ComBat(
#'     x = example_proteome,
#'     sample_annotation = example_sample_annotation,
#'     format = "long",
#'     fill_the_missing = FALSE
#' )
#'
#' # 4) Continuous drift correction (LOESS), then discrete centering if desired
#' test_peptides <- unique(example_proteome$peptide_group_label)[1:3]
#' test_df <- subset(example_proteome, peptide_group_label %in% test_peptides)
#' adjusted_df <- adjust_batch_trend_df(
#'     df_long = test_df,
#'     sample_annotation = example_sample_annotation,
#'     span = 0.7, min_measurements = 8
#' )
#' plot_fit <- plot_with_fitting_curve(
#'     unique(adjusted_df$peptide_group_label),
#'     df_long = adjusted_df, measure_col = "preTrendFit_Intensity",
#'     fit_df = adjusted_df, sample_annotation = example_sample_annotation
#' )
#'
#' # 5) One-call wrapper (continuous + discrete)
#' batch_corrected_matrix <- correct_batch_effects(
#'     x = example_proteome, sample_annotation = example_sample_annotation,
#'     format = "long",
#'     continuous_func = "loess_regression",
#'     discrete_func = "MedianCentering",
#'     batch_col = "MS_batch",
#'     span = 0.7, min_measurements = 8
#' )
#'
#' @seealso
#' \code{\link{center_feature_batch}},
#' \code{\link{adjust_batch_trend_df}},
#' \code{\link{adjust_batch_trend_dm}},
#' \code{\link{correct_with_ComBat}},
#' \code{\link{correct_with_removeBatchEffect}},
#' \code{\link{correct_batch_effects}}
#'
#' @references
#' Johnson WE, Li C, Rabinovic A (2007). Adjusting batch effects in microarray
#' expression data using empirical Bayes methods. \emph{Biostatistics} 8(1):118–127.
#' Smyth GK (2025). \emph{limma User's Guide}, Bioconductor (see removeBatchEffect / lmFit).
#' Leek JT et al. (2024). \emph{sva} vignette, Bioconductor.
#'
#' @name correct_batch_effects
NULL

#' @title Center features per-batch by median/mean (unified)
#' @description Centers each feature *within each batch* to the global
#' location (median/mean) of that feature. Works for long-data frames and
#' wide matrices.
#'
#' @inheritParams correct_batch_effects
#' @param x Either a long data.frame or a numeric matrix (features in rows,
#'   samples in columns), depending on `format`.
#' @param format One of `"long"` or `"wide"`. `"long"` expects a data.frame
#'   with columns `sample_id_col`, `feature_id_col`, `measure_col`. `"wide"`
#'   expects a numeric matrix and uses `matrix_to_long()` / `long_to_matrix()`.
#' @param stat One of `"medians"` or `"means"`. Aliases `"median"`/`"mean"`
#'   are accepted.
#' @param keep_all Passed to `subset_keep_cols()` for long format, ignored for wide.
#' @param no_fit_imputed If `TRUE` and `qual_col` is provided, imputed values
#'   are masked during location estimation (original values remain unchanged).
#' @param qual_col,qual_value Column and value that flag imputed entries.
#'
#' @return If `format = "long"`, a data.frame; if `format = "wide"`, a matrix.
#'
#' @details
#' For `"wide"`, conversion uses `matrix_to_long()` / `long_to_matrix()`.
#' For `"long"`, the function adds `"preBatchCorr_[measure_col]"` and the
#' diagnostic columns:
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
#'     feature_id_col = "peptide_group_label", measure_col = "Intensity"
#' )
#'
#' # WIDE
#' data(example_proteome_matrix, package = "proBatch")
#' out_wide <- center_feature_batch(
#'     x = example_proteome_matrix, sample_annotation = example_sample_annotation,
#'     format = "wide", stat = "means",
#'     sample_id_col = "FullRunName", batch_col = "MS_batch",
#'     feature_id_col = "peptide_group_label", measure_col = "Intensity"
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
    qual_value = NULL) {
    format <- match.arg(format)
    stat_in <- tolower(stat[1])
    stat_in <- switch(stat_in,
        medians = "median",
        median  = "median",
        means   = "mean",
        mean    = "mean",
        stop("`stat` must be one of 'medians' or 'means'.")
    )
    stat_names <- if (identical(stat_in, "median")) {
        c(batch = "median_batch", global = "median_global", diff = "diff_medians")
    } else {
        c(batch = "mean_batch", global = "mean_global", diff = "diff_means")
    }

    if (identical(format, "wide")) {
        if (length(unique(sample_annotation[[sample_id_col]])) != ncol(x)) {
            stop("format='wide' requires a numeric matrix with features in rows and samples in columns.")
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
        df_long           = df_long,
        sample_annotation = sample_annotation,
        sample_id_col     = sample_id_col,
        batch_col         = batch_col,
        feature_id_col    = feature_id_col,
        measure_col       = measure_col,
        keep_all          = keep_all,
        no_fit_imputed    = no_fit_imputed,
        qual_col          = qual_col,
        qual_value        = qual_value,
        stat              = stat_in,
        stat_names        = stat_names
    )

    if (identical(format, "wide")) {
        # Convert long -> wide
        out_wide <- long_to_matrix(
            corrected_long,
            feature_id_col = feature_id_col,
            measure_col    = measure_col,
            sample_id_col  = sample_id_col
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
    stat_names = c(batch = "median_batch", global = "median_global", diff = "diff_medians")) {
    stat <- match.arg(stat)
    original_cols <- names(df_long)

    # Merge/check annotations; ensure batch present post-merge
    df_long <- check_sample_consistency(
        sample_annotation, sample_id_col, df_long,
        batch_col,
        order_col = NULL, facet_col = NULL, merge = TRUE
    )

    if (no_fit_imputed && is.null(qual_col)) {
        warning("`qual_col` is NULL, setting `no_fit_imputed = FALSE`; imputed flags will be ignored.")
        no_fit_imputed <- FALSE
    }

    # Choose location function
    summariser <- switch(stat,
        median = median,
        mean = mean
    )

    # Optionally mask imputed values during inference
    tmp_col <- NULL
    if (isTRUE(no_fit_imputed)) {
        df_long <- .mask_imputed_measure(
            df_long     = df_long,
            measure_col = measure_col,
            qual_col    = qual_col,
            qual_value  = qual_value
        )
        tmp_col <- attr(df_long, "temp_measure_col")
    }
    measure_for_inference <- if (!is.null(tmp_col)) tmp_col else measure_col
    old_measure_col <- .make_pre_col("preBatchCorr", measure_col)

    # Compute per-(batch,feature) and per-feature locations, then apply shift
    corrected_df <- df_long %>%
        group_by(across(any_of(c(batch_col, feature_id_col)))) %>%
        mutate(!!stat_names["batch"] := summariser(.data[[measure_for_inference]], na.rm = TRUE)) %>%
        ungroup() %>%
        group_by(across(any_of(feature_id_col))) %>%
        mutate(!!stat_names["global"] := summariser(.data[[measure_for_inference]], na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(!!stat_names["diff"] := .data[[stat_names["global"]]] - .data[[stat_names["batch"]]]) %>%
        rename(!!old_measure_col := !!measure_col) %>%
        mutate(!!measure_col := .data[[old_measure_col]] + .data[[stat_names["diff"]]])

    # Drop temporary masked column if present
    if (!is.null(tmp_col) && tmp_col %in% names(corrected_df)) {
        corrected_df <- select(corrected_df, -all_of(tmp_col))
    }

    # Column retention for long format
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

    subset_keep_cols(
        df           = corrected_df,
        keep_all     = keep_all,
        default_cols = default_cols,
        minimal_cols = minimal_cols
    )
}

#' @title Adjust batch trend with custom (continuous) fit
#' @description Adjust batch signal trend with the custom (continuous) fit.
#' Should be followed by discrete corrections,
#' e.g. \code{center_feature_batch_medians_df()} or \code{correct_with_ComBat_df()}.
#' Available for both long format data frame (\code{adjust_batch_trend_df()})
#' and data matrix (\code{adjust_batch_trend_dm()}).
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

##############################################################################
# ComBat-based batch correction

#' @title ComBat-based batch correction (unified)
#' @description Adjusts for discrete batch effects using \code{ComBat}.
#' Works for long-format data frames and wide matrices via \code{format}.
#' Optionally accepts covariates through \code{covariates_cols} (passed as \code{mod}).
#'
#' @inheritParams correct_batch_effects
#' @param x Data in long (\code{data.frame}) or wide (\code{matrix}) form, controlled by \code{format}.
#' @param format One of \code{"long"} or \code{"wide"}.
#' @param par.prior Logical; ComBat parametric prior (vs non-parametric).
#' @param covariates_cols Optional character vector of \code{sample_annotation} columns
#'   included in \code{mod} for ComBat (biological or nuisance covariates).
#' @param fill_the_missing Missing-value policy prior to ComBat. If \code{NULL}, no action is taken
#'   and the call will \emph{fail if the matrix contains NA}. Set \code{FALSE} to drop rows with NA,
#'   or provide a numeric value to impute (use with caution).
#' @param keep_all For long format, columns to retain (see \code{subset_keep_cols()}).
#' @param no_fit_imputed If \code{TRUE} and \code{qual_col} provided, masked values are
#'   excluded when building the matrix (original values still corrected).
#' @param mComBat_center If using \code{use_mComBat=TRUE}, the center which should be used
#'  to center the data (https://doi.org/10.1186/s12859-015-0478-3).
#' @param use_mComBat Logical; whether to use the modified ComBat (mComBat) version.
#'
#' @return Matrix if \code{format="wide"}, data.frame if \code{format="long"}.
#' @references Johnson WE et al. (2007) \emph{Biostatistics} 8(1):118–127; \emph{sva} vignette.
#' @export
correct_with_ComBat <- function(
    x, sample_annotation = NULL,
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    format = c("long", "wide"),
    par.prior = TRUE,
    covariates_cols = NULL,
    fill_the_missing = NULL,
    keep_all = "default",
    no_fit_imputed = TRUE,
    qual_col = NULL,
    qual_value = NULL,
    mComBat_center = NULL,
    use_mComBat = FALSE) {
    format <- match.arg(format)

    if (identical(format, "wide")) {
        if (!is.matrix(x)) stop("format='wide' requires a numeric matrix.")
        corrected_matrix <- .combat_matrix_step(
            data_matrix = x,
            sample_annotation = sample_annotation,
            batch_col = batch_col,
            sample_id_col = sample_id_col,
            par.prior = par.prior,
            fill_the_missing = fill_the_missing,
            covariates_cols = covariates_cols,
            mComBat_center = mComBat_center,
            use_mComBat = use_mComBat
        )
        return(corrected_matrix)
    }

    # LONG format
    if (!is.data.frame(x)) stop("format='long' requires a data.frame.")
    df_long <- x
    original_cols <- names(df_long)

    df_long <- check_sample_consistency(
        sample_annotation, sample_id_col, df_long,
        batch_col,
        order_col = NULL, facet_col = NULL, merge = FALSE
    )

    # Handle missingness first; trim df_long / sample_annotation consistently
    qual_for_matrix <- if (no_fit_imputed) qual_col else NULL
    qual_val_for_matrix <- if (no_fit_imputed) qual_value else NULL
    handled <- .handle_missing_for_batch_df(
        df_long = df_long,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col,
        fill_the_missing = fill_the_missing,
        warning_message = "ComBat cannot operate with missing values in the matrix",
        qual_col = qual_for_matrix,
        qual_value = qual_val_for_matrix
    )
    df_long <- handled$df_long
    sample_annotation <- handled$sample_annotation

    # Build matrix AFTER trimming/masking; ComBat on synchronized matrix
    data_matrix <- long_to_matrix(
        df_long,
        feature_id_col = feature_id_col,
        measure_col    = measure_col,
        sample_id_col  = sample_id_col,
        qual_col       = qual_for_matrix,
        qual_value     = qual_val_for_matrix
    )

    # ComBat on matrix (method ensures numeric & SA alignment)
    corrected_matrix <- .combat_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        batch_col = batch_col,
        sample_id_col = sample_id_col,
        par.prior = par.prior,
        fill_the_missing = fill_the_missing,
        covariates_cols = covariates_cols,
        mComBat_center = mComBat_center,
        use_mComBat = use_mComBat
    )

    corrected_df <- matrix_to_long(
        corrected_matrix,
        feature_id_col = feature_id_col,
        measure_col    = measure_col,
        sample_id_col  = sample_id_col
    )

    old_measure_col <- .make_pre_col("preBatchCorr", measure_col)
    df_long <- rename(df_long, !!old_measure_col := !!measure_col)

    corrected_df <- left_join(
        corrected_df,
        df_long,
        by = setNames(c(feature_id_col, sample_id_col), c(feature_id_col, sample_id_col))
    )

    default_cols <- c(original_cols, old_measure_col)
    minimal_cols <- c(sample_id_col, feature_id_col, measure_col, old_measure_col)

    subset_keep_cols(
        corrected_df, keep_all,
        default_cols = default_cols,
        minimal_cols = minimal_cols
    )
}

# ---- removeBatchEffect unified ------------------------------------------------

#' @title Batch effect correction with limma::removeBatchEffect (unified)
#' @description Removes batch-associated linear effects with removeBatchEffect from limma.
#' Works for long or wide via \code{format}. Use \code{covariates_cols}
#' to keep biological effects in the design (not removed).
#' @inheritParams correct_with_ComBat
#' @param covariates_cols Optional \code{sample_annotation} columns for the design matrix (biological or nuisance covariates).
#' @param fill_the_missing Missing-value policy applied before modeling. If \code{NULL} (default),
#'   \emph{NA values in the data matrix are left as is} and handled by limma's linear modeling;
#'   the design matrix (\code{batch_col} and \code{covariates_cols}) must be NA-free.
#'   Set \code{FALSE} to drop rows with NA, or provide a numeric value to impute explicitly.
#'
#' @return Matrix if \code{format="wide"}, data.frame if \code{format="long"} with batch effects removed
#' @examples
#' data(
#'     list = c("example_sample_annotation", "example_proteome_matrix"),
#'     package = "proBatch"
#' )
#' batch_corrected_matrix <- correct_with_removeBatchEffect(
#'     example_proteome_matrix,
#'     example_sample_annotation,
#'     batch_col = "MS_batch",
#'     covariates_cols = c("Condition", "Type")
#' )
#' @seealso \code{\link{removeBatchEffect}}
#' @export
correct_with_removeBatchEffect <- function(
    x, sample_annotation,
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    format = c("long", "wide"),
    covariates_cols = NULL,
    fill_the_missing = NULL,
    keep_all = "default",
    ...) {
    format <- match.arg(format)

    if (identical(format, "wide")) {
        if (!is.matrix(x)) stop("format='wide' requires a numeric matrix.")
        corrected_matrix <- .removeBatchEffect_matrix_step(
            data_matrix       = x,
            sample_annotation = sample_annotation,
            batch_col         = batch_col,
            sample_id_col     = sample_id_col,
            covariates_cols   = covariates_cols,
            fill_the_missing  = fill_the_missing,
            ...
        )
        return(corrected_matrix)
    }

    # LONG
    if (!is.data.frame(x)) stop("format='long' requires a data.frame.")
    df_long <- x
    original_cols <- names(df_long)

    # Pre-handle missingness the same way as ComBat
    handled <- .handle_missing_for_batch_df(
        df_long = df_long,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col,
        fill_the_missing = fill_the_missing,
        warning_message = if (is.null(fill_the_missing) || !fill_the_missing) {
            "removeBatchEffect will leave NA as-is in the matrix; design matrix (batch/covariates) must be free of NA."
        } else {
            "removeBatchEffect can operate with missing values; applying requested NA handling before modeling."
        },
        qual_col = NULL,
        qual_value = NULL
    )

    df_long <- handled$df_long
    sample_annotation <- handled$sample_annotation

    data_matrix <- long_to_matrix(
        df_long,
        feature_id_col = feature_id_col,
        measure_col    = measure_col,
        sample_id_col  = sample_id_col
    )

    corrected_matrix <- .removeBatchEffect_matrix_step(
        data_matrix       = data_matrix,
        sample_annotation = sample_annotation,
        batch_col         = batch_col,
        sample_id_col     = sample_id_col,
        covariates_cols   = covariates_cols,
        fill_the_missing  = fill_the_missing,
        ...
    )

    corrected_df <- matrix_to_long(
        corrected_matrix,
        feature_id_col = feature_id_col,
        measure_col    = measure_col,
        sample_id_col  = sample_id_col
    )

    old_measure_col <- .make_pre_col("preBatchCorr", measure_col)
    df_long <- rename(df_long, !!old_measure_col := !!measure_col)
    corrected_df <- merge(corrected_df, df_long, by = c(feature_id_col, sample_id_col))

    default_cols <- c(original_cols, old_measure_col)
    minimal_cols <- c(sample_id_col, feature_id_col, measure_col, old_measure_col)

    subset_keep_cols(
        corrected_df, keep_all,
        default_cols = default_cols,
        minimal_cols = minimal_cols
    )
}


# ---- Orchestrator (unified) --------------------------------------------------

#' @title Combine continuous and discrete corrections (unified)
#' @description Optional continuous drift removal + discrete adjustment
#' via \code{"MedianCentering"}, \code{"MeanCentering"}, \code{"ComBat"},
#' or \code{"removeBatchEffect"}. Works for long or wide via \code{format}.
#' @inheritParams correct_with_ComBat
#' @param continuous_func e.g. \code{"loess_regression"} or \code{NULL}.
#' @param discrete_func batch method name.
#' @export
correct_batch_effects <- function(
    x, sample_annotation,
    format = c("long", "wide"),
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
    par.prior = TRUE,
    covariates_cols = NULL,
    min_measurements = 8,
    ...) {
    format <- match.arg(format)
    discrete_func <- match.arg(discrete_func)

    # Standardize to LONG for the pipeline, then back if needed
    if (identical(format, "wide")) {
        stopifnot(is.matrix(x))
        df_long <- matrix_to_long(
            data_matrix    = x,
            feature_id_col = feature_id_col,
            measure_col    = measure_col,
            sample_id_col  = sample_id_col
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
        warning_message = "Batch correction cannot operate with missing values in the matrix",
        qual_col = if (no_fit_imputed) qual_col else NULL,
        qual_value = if (no_fit_imputed) qual_value else NULL
    )
    df_long <- handled$df_long
    sample_annotation <- handled$sample_annotation

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
    df_long <- switch(discrete_func,
        MedianCentering = center_feature_batch(
            x = df_long, sample_annotation = sample_annotation,
            format = "long", stat = "medians",
            sample_id_col = sample_id_col, batch_col = batch_col,
            feature_id_col = feature_id_col, measure_col = measure_col,
            keep_all = keep_all, no_fit_imputed = no_fit_imputed,
            qual_col = qual_col, qual_value = qual_value
        ),
        MeanCentering = center_feature_batch(
            x = df_long, sample_annotation = sample_annotation,
            format = "long", stat = "means",
            sample_id_col = sample_id_col, batch_col = batch_col,
            feature_id_col = feature_id_col, measure_col = measure_col,
            keep_all = keep_all, no_fit_imputed = no_fit_imputed,
            qual_col = qual_col, qual_value = qual_value
        ),
        ComBat = correct_with_ComBat(
            x = df_long, sample_annotation = sample_annotation,
            feature_id_col = feature_id_col, measure_col = measure_col,
            sample_id_col = sample_id_col, batch_col = batch_col,
            format = "long", par.prior = par.prior,
            covariates_cols = covariates_cols,
            fill_the_missing = fill_the_missing,
            keep_all = keep_all, no_fit_imputed = no_fit_imputed,
            qual_col = qual_col, qual_value = qual_value
        ),
        removeBatchEffect = correct_with_removeBatchEffect(
            x = df_long, sample_annotation = sample_annotation,
            feature_id_col = feature_id_col, measure_col = measure_col,
            sample_id_col = sample_id_col, batch_col = batch_col,
            format = "long", covariates_cols = covariates_cols,
            fill_the_missing = fill_the_missing, keep_all = keep_all, ...
        )
    )

    if (back_to_wide) {
        return(long_to_matrix(
            df_long,
            feature_id_col = feature_id_col,
            measure_col    = measure_col,
            sample_id_col  = sample_id_col
        ))
    }

    # Ensure provenance columns are retained for long
    old_measure_col <- .make_pre_col("preBatchCorr", measure_col)
    default_cols <- unique(c(
        original_cols, batch_col, old_measure_col,
        if (!is.null(continuous_func)) c(.make_pre_col("preTrendFit", measure_col), "fit")
    ))
    minimal_cols <- c(sample_id_col, feature_id_col, measure_col, old_measure_col, batch_col)

    subset_keep_cols(df_long, keep_all,
        default_cols = default_cols, minimal_cols = minimal_cols
    )
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

# Core ComBat matrix call with optional covariates (mod)
run_ComBat_core <- function(sample_annotation, batch_col, data_matrix,
                            par.prior, covariates_cols = NULL,
                            mComBat_center = NULL, use_mComBat = FALSE, ...) {
    if (is.null(sample_annotation)) {
        stop("sample_annotation is required for ComBat correction")
    }

    sample_annotation <- as.data.frame(sample_annotation)
    if (!(batch_col %in% names(sample_annotation))) {
        stop("Batch column is not present in sample_annotation")
    }

    # ONE batch factor only (ComBat constraint).
    batches <- factor(sample_annotation[[batch_col]])
    if (use_mComBat) {
        # Check if mComBat_center is valid and present among batches
        if (is.null(mComBat_center) || !(mComBat_center %in% levels(batches))) {
            stop("mComBat_center must be specified and present in the batch levels when use_mComBat is TRUE.")
        }
    }

    if (!is.null(covariates_cols) && length(covariates_cols)) {
        missing_cov <- setdiff(covariates_cols, names(sample_annotation))
        if (length(missing_cov)) {
            stop("Covariates missing in sample_annotation: ", paste(missing_cov, collapse = ", "))
        }
        covariates <- as.data.frame(sample_annotation[, covariates_cols, drop = FALSE])
        mod <- model.matrix(~., data = covariates)
    } else {
        mod <- model.matrix(~1, data = sample_annotation)
    }

    if (use_mComBat) {
        message(
            "Correction using M-ComBat with center batch: ", mComBat_center,
            ". \n\t\tSee https://doi.org/10.1186/s12859-015-0478-3 for details."
        )
        .m_COMBAT(
            dat = data_matrix, batch = batches, center = mComBat_center, mod = mod
        )
    } else {
        ComBat(
            dat = data_matrix, batch = batches, mod = mod,
            par.prior = par.prior, ...
        )
    }
}

.mComBat_matrix_step <- function(data_matrix, sample_annotation,
                                 batch_col = "MS_batch",
                                 sample_id_col = NULL,
                                 fill_the_missing = NULL,
                                 covariates_cols = NULL,
                                 mComBat_center = NULL,
                                 use_mComBat = TRUE,
                                 ...) {
    .combat_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        batch_col = batch_col,
        sample_id_col = sample_id_col,
        fill_the_missing = fill_the_missing,
        covariates_cols = covariates_cols,
        mComBat_center = mComBat_center,
        use_mComBat = use_mComBat,
        ...
    )
}

.combat_matrix_step <- function(data_matrix, sample_annotation,
                                batch_col = "MS_batch",
                                sample_id_col = NULL,
                                par.prior = TRUE,
                                fill_the_missing = NULL,
                                covariates_cols = NULL,
                                mComBat_center = NULL,
                                use_mComBat = FALSE,
                                ...) {
    .run_matrix_method(
        data_matrix, sample_annotation,
        sample_id_col = sample_id_col,
        fill_the_missing = fill_the_missing,
        missing_warning = "ComBat cannot operate with missing values in the matrix",
        method_fun = function(data_matrix, sample_annotation) {
            run_ComBat_core(
                sample_annotation = sample_annotation,
                batch_col = batch_col,
                data_matrix = data_matrix,
                par.prior = par.prior,
                covariates_cols = covariates_cols,
                mComBat_center = mComBat_center,
                use_mComBat = use_mComBat
            )
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

.run_matrix_method <- function(data_matrix, sample_annotation,
                               sample_id_col = NULL,
                               fill_the_missing = NULL,
                               missing_warning = "This method cannot operate with missing values in the matrix",
                               method_fun) {
    # ensure numeric matrix input for downstream modeling (sva/limma)
    if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
    }
    storage.mode(data_matrix) <- "double"
    if (!is.numeric(data_matrix)) {
        stop("Input must be coercible to a numeric matrix for batch correction.")
    }

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
