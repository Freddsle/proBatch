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
#' @param batch_col Optional batch variable name present in `design_formula`.
#'   If provided and omicsGMF stores design attributes (`X`, `Beta`) on the
#'   reduced dimension, reconstruction preserves non-batch fixed effects by
#'   adding `X_no_batch %*% Beta` to the latent component.
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
  batch_col = NULL,
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

    .pb_require_omicsgmf_stack()

    format <- match.arg(format)
    batch_col <- .omicsgmf_normalize_batch_col(batch_col)
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
                batch_col = batch_col,
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
            batch_col = batch_col,
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
  batch_col = NULL,
  family = gaussian(),
  ncomponents,
  gmf_args = list(),
  impute_args = list()
) {
    .pb_require_omicsgmf_stack()

    ncomponents <- .pb_positive_integer(ncomponents, "ncomponents")

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
            reconstructed <- .omicsgmf_reconstruct_corrected_matrix(
                gmf_results = gmf_results,
                data_matrix = data_matrix,
                batch_col = batch_col
            )

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

.omicsgmf_reconstruct_corrected_matrix <- function(gmf_results, data_matrix, batch_col = NULL) {
    rotation_matrix <- attr(gmf_results, "rotation")
    if (is.null(rotation_matrix)) {
        stop(
            "Rotation matrix attribute missing from GMF results; cannot compute corrected intensities.",
            call. = FALSE
        )
    }

    latent_component <- t(gmf_results %*% t(rotation_matrix))
    reconstructed <- latent_component

    # If batch_col is supplied, preserve non-batch fixed-effect terms from X %*% Beta.
    if (!is.null(batch_col)) {
        design_terms <- .omicsgmf_extract_design_terms(gmf_results)
        if (is.null(design_terms)) {
            warning(
                "Could not access omicsGMF design attributes (X/Beta); returning latent-only reconstruction."
            )
        } else {
            X <- design_terms$X
            Beta <- design_terms$Beta

            dims_ok <- nrow(X) == ncol(latent_component) &&
                ncol(Beta) == nrow(latent_component)
            if (!dims_ok) {
                warning(
                    "omicsGMF design attributes are incompatible with latent dimensions; returning latent-only reconstruction."
                )
            } else {
                batch_idx <- .omicsgmf_match_batch_columns(colnames(X), batch_col)
                if (!length(batch_idx)) {
                    warning(
                        "`batch_col` was not found among omicsGMF design columns; returning latent-only reconstruction."
                    )
                } else {
                    X_keep <- X
                    X_keep[, batch_idx] <- 0
                    reconstructed <- latent_component + t(X_keep %*% Beta)
                }
            }
        }
    }

    storage.mode(reconstructed) <- "double"
    reconstructed
}

.omicsgmf_extract_design_terms <- function(gmf_results) {
    X <- attr(gmf_results, "X")
    Beta <- attr(gmf_results, "Beta")
    if (is.null(X) || is.null(Beta)) {
        return(NULL)
    }
    X <- as.matrix(X)
    Beta <- as.matrix(Beta)
    if (!nrow(X) || !ncol(X) || !nrow(Beta) || !ncol(Beta)) {
        return(NULL)
    }
    if (ncol(X) != nrow(Beta)) {
        return(NULL)
    }
    list(X = X, Beta = Beta)
}

.omicsgmf_normalize_batch_col <- function(batch_col) {
    if (is.null(batch_col)) {
        return(NULL)
    }

    if (length(batch_col) != 1L) {
        stop("`batch_col` must be NULL or a single non-empty column name.", call. = FALSE)
    }

    batch_col <- as.character(batch_col[[1]])
    if (is.na(batch_col) || !nzchar(batch_col)) {
        stop("`batch_col` must be NULL or a single non-empty column name.", call. = FALSE)
    }

    batch_col
}

.omicsgmf_match_batch_columns <- function(design_columns, batch_col) {
    if (
        is.null(design_columns) || !length(design_columns) ||
            is.null(batch_col)
    ) {
        return(integer())
    }

    batch_col <- .omicsgmf_normalize_batch_col(batch_col)
    batch_col <- gsub("`", "", batch_col, fixed = TRUE)
    design_columns <- as.character(design_columns)
    design_columns <- gsub("`", "", design_columns, fixed = TRUE)

    # Also match interaction/function terms such as DietB:MS_batchb2 or factor(MS_batch)b2.
    split_terms <- strsplit(design_columns, ":", fixed = TRUE)
    has_batch <- vapply(
        split_terms,
        function(parts) {
            any(vapply(
                parts,
                function(part) {
                    part <- trimws(part)
                    startsWith(part, batch_col) ||
                        grepl(paste0("(", batch_col, ")"), part, fixed = TRUE)
                },
                logical(1)
            ))
        },
        logical(1)
    )

    which(has_batch)
}
