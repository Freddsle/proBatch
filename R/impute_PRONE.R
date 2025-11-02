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
        cols <- names(x)
        looks_long <- all(c(feature_id_col, sample_id_col, measure_col) %in% cols)
        if (looks_long) {
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

    data_matrix <- long_to_matrix(
        df_long = x,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col
    )

    imputed_matrix <- .prone_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        condition_col = condition_col,
        assay_in = assay_in
    )

    matrix_to_long(
        data_matrix = imputed_matrix,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col
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

    if (is.null(sample_annotation)) {
        col_data <- S4Vectors::DataFrame(row.names = sample_ids)
    } else {
        aligned_sa <- .align_sample_annotation(
            sample_annotation = sample_annotation,
            sample_ids = sample_ids,
            sample_id_col = sample_id_col
        )
        col_data <- S4Vectors::DataFrame(aligned_sa)
        rownames(col_data) <- sample_ids
    }

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = setNames(list(data_matrix), assay_in),
        colData = col_data
    )

    se_imp <- PRONE::impute_se(
        se,
        ain = assay_in,
        condition = condition_col
    )

    SummarizedExperiment::assay(se_imp, assay_in)
}
