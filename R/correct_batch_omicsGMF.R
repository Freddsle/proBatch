#' Batch correction with omicsGMF
#'
#' @description
#' Removes batch structure by fitting the \pkg{omicsGMF} matrix factorisation model
#' and reconstructing intensities from the learned latent space. Works on long
#' data.frames (`format = "long"`) and wide matrices (`format = "wide"`) by
#' converting between representations internally.
#'
#' Before applying this correction, run [estimate_omicsGMF_rank()] to determine
#' the appropriate latent dimensionality (`ncomponents`). The method relies on
#' the same omicsGMF workflow used for imputation, so optional arguments can be
#' forwarded via `gmf_args` and `impute_args` lists.
#'
#' @param x Data to correct. Provide either a long-format `data.frame` or a numeric
#'   matrix in wide format.
#' @param sample_annotation Data frame with sample-level metadata referenced in
#'   `design_formula`. Must contain `sample_id_col`.
#' @param ncomponents Positive integer latent dimensionality returned by
#'   [estimate_omicsGMF_rank()].
#' @param feature_id_col Column identifying features in long-format data. Ignored for matrices.
#' @param measure_col Column containing intensity/abundance measurements in long-format data.
#' @param sample_id_col Column identifying samples.
#' @param format One of `"long"` or `"wide"`.
#' @param design_formula Formula (or character string coercible to one) used to
#'   build the omicsGMF design matrix from `sample_annotation`.
#' @param family GLM family passed to omicsGMF (default: `gaussian()`).
#' @param keep_all Columns retained in the long-format output; passed to [subset_keep_cols()].
#' @param ... Optional named arguments. Use `gmf_args` and `impute_args` lists to
#'   forward parameters to `omicsGMF::runGMF()` and `omicsGMF::imputeGMF()`.
#'
#' @return A numeric matrix (for `format = "wide"`) or a long `data.frame` with
#'   corrected measurements and preserved pre-correction values in
#'   `preBatchCorr_[measure_col]`.
#' @seealso [impute_with_omicsGMF()], [estimate_omicsGMF_rank()]
#' @export
correct_with_omicsGMF <- function(
  x,
  sample_annotation,
  ncomponents,
  feature_id_col = "peptide_group_label",
  measure_col = "Intensity",
  sample_id_col = "FullRunName",
  format = c("long", "wide"),
  design_formula = ~1,
  family = gaussian(),
  keep_all = "default",
  ...
) {
    if (missing(ncomponents) || is.null(ncomponents)) {
        stop("`ncomponents` must be supplied. Run estimate_omicsGMF_rank() first.", call. = FALSE)
    }
    if (missing(sample_annotation) || is.null(sample_annotation)) {
        stop('argument "sample_annotation" is missing, with no default', call. = FALSE)
    }

    .pb_requireNamespace("omicsGMF")
    .pb_requireNamespace("sgdGMF")
    .pb_requireNamespace("SingleCellExperiment")

    format <- match.arg(format)
    dots <- list(...)
    if ("ntop" %in% names(dots)) {
        warning("forcing ntop = NULL for omicsGMF correction")
        dots$ntop <- NULL
    }

    if (identical(format, "wide")) {
        if (!is.matrix(x)) {
            stop("format = 'wide' requires a numeric matrix.", call. = FALSE)
        }
        matrix_args <- c(
            list(
                data_matrix = x,
                sample_annotation = sample_annotation,
                sample_id_col = sample_id_col,
                design_formula = design_formula,
                family = family,
                ncomponents = ncomponents
            ),
            dots
        )
        return(do.call(.omicsgmf_correct_matrix_step, matrix_args))
    }

    prep <- .pb_prepare_long_matrix(
        df_long = x,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        batch_col = NULL,
        error_message = "format = 'long' requires a data.frame.",
        error_call = FALSE
    )
    df_long <- prep$df_long
    sample_annotation <- prep$sample_annotation
    data_matrix <- prep$data_matrix
    original_cols <- prep$original_cols

    aligned_sa <- .align_sample_annotation(
        sample_annotation = sample_annotation,
        sample_ids = colnames(data_matrix),
        sample_id_col = sample_id_col
    )

    matrix_args <- c(
        list(
            data_matrix = data_matrix,
            sample_annotation = aligned_sa,
            sample_id_col = sample_id_col,
            design_formula = design_formula,
            family = family,
            ncomponents = ncomponents
        ),
        dots
    )
    corrected_matrix <- do.call(.omicsgmf_correct_matrix_step, matrix_args)

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

.omicsgmf_correct_matrix_step <- function(
  data_matrix,
  sample_annotation,
  sample_id_col = "FullRunName",
  design_formula = ~1,
  family = gaussian(),
  ncomponents,
  gmf_args = list(),
  impute_args = list()
) {
    .pb_requireNamespace("SingleCellExperiment")
    .pb_requireNamespace("omicsGMF")
    .pb_requireNamespace("sgdGMF")

    ncomponents <- as.integer(ncomponents)
    if (is.na(ncomponents) || ncomponents < 1L) {
        stop("`ncomponents` must be a positive integer.", call. = FALSE)
    }

    .run_matrix_method(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        fill_the_missing = NULL,
        missing_warning = "omicsGMF correction removed rows/columns while handling missing values.",
        method_fun = function(data_matrix, sample_annotation) {
            fit <- .omicsgmf_fit_and_impute(
                data_matrix = data_matrix,
                sample_annotation = sample_annotation,
                design_formula = design_formula,
                family = family,
                ncomponents = ncomponents,
                gmf_args = gmf_args,
                impute_args = impute_args
            )

            gmf_results <- SingleCellExperiment::reducedDim(fit$sce, fit$dimred_name)
            if (is.null(gmf_results)) {
                stop(
                    sprintf("Reduced dimension '%s' not found after runGMF.", fit$dimred_name),
                    call. = FALSE
                )
            }
            rotation_matrix <- attr(gmf_results, "rotation")
            if (is.null(rotation_matrix)) {
                stop("Rotation matrix attribute missing from GMF results; cannot compute corrected intensities.", call. = FALSE)
            }

            reconstructed <- gmf_results %*% t(rotation_matrix)
            reconstructed <- t(reconstructed)

            if (!is.null(rownames(data_matrix))) {
                rownames(reconstructed) <- rownames(data_matrix)
            }
            if (!is.null(colnames(data_matrix))) {
                colnames(reconstructed) <- colnames(data_matrix)
            }
            storage.mode(reconstructed) <- "double"
            reconstructed
        }
    )
}
