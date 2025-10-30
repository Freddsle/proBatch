#' Impute missing values with omicsGMF
#'
#' @description
#' Function [impute_with_omicsGMF()] reconstructs missing values using the \pkg{omicsGMF} 
#' probabilistic matrix factorisation workflow. Works on long data.frames (`format = "long"`), 
#' wide matrices (`format = "wide"`), and can be applied to `ProBatchFeatures` pipelines via
#' [pb_transform()] (see Details).
#'
#' Before using, run [estimate_omicsGMF_rank()] to run [`omicsGMF::runRankGMF()`] and determine
#' the latent dimensionality prior to fitting [`omicsGMF::runGMF()`] and
#' imputing with [`omicsGMF::imputeGMF()`]. Provide the resulting `ncomponents`
#' explicitly here.
#' 
#' For more details on the omicsGMF method, see the package documentation:
#' <https://bioconductor.org/packages/omicsGMF/> and the GitHub repository:
#' <https://github.com/statOmics/omicsGMF>.
#'
#' @param x Data to impute. Supply either a long-format `data.frame`, a numeric
#'   matrix in wide format, or a `ProBatchFeatures` object.
#' @param sample_annotation Data frame describing samples. Required for matrices
#'   and long data; inferred from `x` when it is a `ProBatchFeatures`.
#' @param format One of `"long"` or `"wide"`. Ignored when `x` is a
#'   `ProBatchFeatures` object.
#' @param sample_id_col Column name identifying samples.
#' @param batch_col Optional column indicating batch membership in long-format data.
#' @param design_formula Formula (or character string coercible to one) used to
#'   build the design matrix from sample_annotation and passed as `X` to omicsGMF 
#'   - sample-level covariate matrix. Defaults to `~ 1`.
#' @param family GLM family passed to omicsGMF (default: `gaussian()`).
#' @param ncomponents Positive integer with the latent dimensionality returned by
#'   [estimate_omicsGMF_rank()].
#' @param ntop Numeric scalar specifying the number of features with the highest variances to
#'   use for dimensionality reduction. Default uses all features.
#' @param feature_id_col Column name identifying features in long-format data. Ignored for wide matrices.
#' @param measure_col Column name storing intensity/abundance measurements in long-format data. 
#'  Ignored for wide matrices.
#' @param keep_all Columns retained in the long-format output; passed to
#'   [subset_keep_cols()]. Ignored for wide matrices.
#' @param pbf_name When `x` is a `ProBatchFeatures`, source assay name. Defaults
#'   to [pb_current_assay()].
#' @param final_name Optional name for the stored assay when `x` is a
#'   `ProBatchFeatures`.
#' @param max_rank Upper bound supplied to `runRankGMF()` via `maxcomp`.
#'   Applicable to [estimate_omicsGMF_rank()] only.
#' @param ... Additional arguments passed as lists to omicsGMF functions:
#'   - `omicsGMF::runGMF()`,
#'   - `omicsGMF::imputeGMF()`, or
#'   - `omicsGMF::runRankGMF()`.
#'   For the full list of parameters, see the respective function documentation.
#' 
#' @details
#' For `ProBatchFeatures` inputs, the function delegates to [pb_transform()] with
#' the registered `"omicsGMFImpute"` step. Any unspecified `sample_annotation`
#' argument is sourced from the object's `colData`.
#'
#' The helper requires the following packages to be installed:
#' \pkg{omicsGMF}, \pkg{sgdGMF}, and \pkg{SingleCellExperiment}. They are checked
#' at run time via `.pb_requireNamespace()`.
#'
#' @return
#' `impute_with_omicsGMF()` returns a matrix (for `format = "wide"`), a long
#' data.frame (for `format = "long"`), or a `ProBatchFeatures` instance when `x`
#' is already a pipeline object. `estimate_omicsGMF_rank()` returns a list with
#' the ranked `SingleCellExperiment` (`sce`) and the inferred `ncomponents`
#' (integer or `NULL` when the rank could not be determined).
#'
#' @seealso [handle_missing_values()], [pb_transform()],
#'   `omicsGMF::runGMF`, `omicsGMF::imputeGMF`
#' @export
NULL 
 

#' @rdname impute_with_omicsGMF
#' @export
impute_with_omicsGMF <- function(
  x,
  sample_annotation,
  ncomponents,
  feature_id_col = "peptide_group_label",
  measure_col = "Intensity",
  sample_id_col = "FullRunName",
  batch_col = NULL,
  format = c("long", "wide"),
  design_formula = ~1,
  family = gaussian(),
  keep_all = "default",
  pbf_name = NULL,
  final_name = NULL,
  ...
) {
    if (missing(ncomponents) || is.null(ncomponents)) {
        stop("`ncomponents` must be supplied. Run estimate_omicsGMF_rank() first.", call. = FALSE)
    }

    if (is(x, "ProBatchFeatures")) {
        return(impute_with_omicsGMF.ProBatchFeatures(
            x = x,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            design_formula = design_formula,
            family = family,
            ncomponents = ncomponents,
            pbf_name = pbf_name,
            final_name = final_name,
            ...
        ))
    }

    .pb_requireNamespace("omicsGMF")
    .pb_requireNamespace("sgdGMF")

    format <- match.arg(format)

    if (identical(format, "wide")) {
        if (!is.matrix(x)) {
            stop("format = 'wide' requires a numeric matrix.", call. = FALSE)
        }
        return(.omicsgmf_matrix_step(
            data_matrix = x,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            design_formula = design_formula,
            family = family,
            ncomponents = ncomponents,
            ...
        ))
    }

    if (!is.data.frame(x)) {
        stop("format = 'long' requires a data.frame.", call. = FALSE)
    }
    df_long <- x
    original_cols <- names(df_long)

    df_long <- check_sample_consistency(
        sample_annotation, sample_id_col, df_long,
        batch_col,
        order_col = NULL, facet_col = NULL, merge = FALSE
    )

    data_matrix <- long_to_matrix(
        df_long,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col,
        qual_col = NULL
    )

    imputed_matrix <- .omicsgmf_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        design_formula = design_formula,
        family = family,
        ncomponents = ncomponents,
        ...
    )

    imputed_df <- matrix_to_long(
        imputed_matrix,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col
    )

    old_measure_col <- .make_pre_col("preImpute", measure_col)
    df_long <- rename(df_long, !!old_measure_col := !!sym(measure_col))

    imputed_df <- left_join(
        imputed_df,
        df_long,
        by = setNames(
            c(feature_id_col, sample_id_col),
            c(feature_id_col, sample_id_col)
        )
    )

    default_cols <- unique(c(original_cols, old_measure_col))
    minimal_cols <- c(sample_id_col, feature_id_col, measure_col, old_measure_col)

    subset_keep_cols(
        imputed_df,
        keep_all,
        default_cols = default_cols,
        minimal_cols = minimal_cols
    )
}

#' @rdname impute_with_omicsGMF
#' @export
impute_with_omicsGMF.ProBatchFeatures <- function(
  x,
  sample_annotation = NULL,
  ncomponents,
  sample_id_col = "FullRunName",
  design_formula = ~1,
  family = gaussian(),
  pbf_name = NULL,
  final_name = NULL,
  ...
) {
    object <- x
    .pb_requireNamespace("omicsGMF")
    .pb_requireNamespace("sgdGMF")

    if (missing(ncomponents) || is.null(ncomponents)) {
        stop("`ncomponents` must be supplied. Run estimate_omicsGMF_rank() first.", call. = FALSE)
    }

    if (is.null(pbf_name)) {
        pbf_name <- pb_current_assay(object)
        message("`pbf_name` not provided, using the most recent assay: ", pbf_name)
    }

    params <- list(
        sample_id_col = sample_id_col,
        design_formula = design_formula,
        family = family,
        ncomponents = ncomponents
    )
    # add ... arguments
    extra_args <- list(...)
    if (length(extra_args)) {
        params <- modifyList(params, extra_args)
    }

    if (!is.null(sample_annotation)) {
        params$sample_annotation <- sample_annotation
    }

    pb_transform(
        object = object,
        from = pbf_name,
        steps = "omicsGMFImpute",
        funs = list(.omicsgmf_matrix_step),
        params_list = list(params),
        final_name = final_name
    )
}

#' @rdname impute_with_omicsGMF
#' @export
estimate_omicsGMF_rank <- function(
  x,
  sample_annotation,
  feature_id_col = "peptide_group_label",
  measure_col = "Intensity",
  sample_id_col = "FullRunName",
  batch_col = NULL,
  format = c("long", "wide"),
  design_formula = ~1,
  family = gaussian(),
  max_rank = 10,
  ...
) {
    .pb_requireNamespace("omicsGMF")
    .pb_requireNamespace("sgdGMF")

    format <- match.arg(format)

    if (identical(format, "wide")) {
        if (!is.matrix(x)) {
            stop("format = 'wide' requires a numeric matrix.", call. = FALSE)
        }
        return(.omicsgmf_rank_matrix_step(
            data_matrix = x,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            design_formula = design_formula,
            family = family,
            max_rank = max_rank,
            ...
        ))
    }

    if (!is.data.frame(x)) {
        stop("format = 'long' requires a data.frame.", call. = FALSE)
    }
    df_long <- x

    df_long <- check_sample_consistency(
        sample_annotation, sample_id_col, df_long,
        batch_col = batch_col,
        order_col = NULL, facet_col = NULL, merge = FALSE
    )

    data_matrix <- long_to_matrix(
        df_long,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col,
        qual_col = NULL
    )

    .omicsgmf_rank_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        design_formula = design_formula,
        family = family,
        max_rank = max_rank,
        rank_args = rank_args
    )
}

#################################################################################
# Internal functions
#################################################################################

.omicsgmf_rank_matrix_step <- function(
  data_matrix,
  sample_annotation,
  sample_id_col = "FullRunName",
  design_formula = ~1,
  family = gaussian(),
  max_rank = 10,
  ...
) {
    .pb_requireNamespace("omicsGMF")
    .pb_requireNamespace("sgdGMF")

    max_rank <- as.integer(max_rank)
    if (is.na(max_rank) || max_rank < 1L) {
        stop("`max_rank` must be a positive integer.", call. = FALSE)
    }

    .run_matrix_method(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        fill_the_missing = NULL,
        missing_warning = "omicsGMF rank estimation removed rows/columns while handling missing values.",
        method_fun = function(data_matrix, sample_annotation) {
            sample_df <- as.data.frame(sample_annotation)
            exprs_name <- "omicsGMF_input"

            sce <- SingleCellExperiment::SingleCellExperiment(
                assays = setNames(list(data_matrix), exprs_name),
                colData = S4Vectors::DataFrame(sample_df)
            )

            design_formula_local <- .omicsgmf_normalize_formula(design_formula)
            X <- stats::model.matrix(
                design_formula_local,
                data = as.data.frame(SummarizedExperiment::colData(sce))
            )

            rank_args_local <- rank_args %||% list()

            .omicsgmf_estimate_ncomponents(
                sce = sce,
                X = X,
                exprs_name = exprs_name,
                family = family,
                max_rank = max_rank,
                rank_args = rank_args_local
            )
        }
    )
}

.omicsgmf_matrix_step <- function(
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
        missing_warning = "omicsGMF imputation removed rows/columns while handling missing values.",
        method_fun = function(data_matrix, sample_annotation) {
            sample_df <- as.data.frame(sample_annotation)
            exprs_name <- "omicsGMF_input"

            sce <- SingleCellExperiment::SingleCellExperiment(
                assays = setNames(list(data_matrix), exprs_name),
                colData = S4Vectors::DataFrame(sample_df)
            )

            design_formula_local <- .omicsgmf_normalize_formula(design_formula)
            X <- stats::model.matrix(design_formula_local, data = as.data.frame(SummarizedExperiment::colData(sce)))

            gmf_args <- gmf_args %||% list()
            impute_args <- impute_args %||% list()

            dimred_name <- gmf_args$name %||% "omicsGMF"
            gmf_call <- modifyList(
                list(object = sce, X = X, exprs_values = exprs_name, family = family, ncomponents = ncomponents, name = dimred_name),
                gmf_args
            )
            sce <- do.call(omicsGMF::runGMF, gmf_call)

            imputed_name <- impute_args$name %||% "omicsGMF_imputed"
            impute_call <- modifyList(
                list(
                    object = sce,
                    exprs_values = exprs_name,
                    reducedDimName = dimred_name,
                    name = imputed_name
                ),
                impute_args
            )
            sce <- do.call(omicsGMF::imputeGMF, impute_call)

            final_impute_name <- impute_call$name %||% imputed_name
            imputed <- SummarizedExperiment::assay(sce, final_impute_name)
            storage.mode(imputed) <- "double"
            imputed
        }
    )
}

.omicsgmf_normalize_formula <- function(formula_or_string) {
    if (is.null(formula_or_string)) {
        return(stats::as.formula("~ 1"))
    }
    if (is.character(formula_or_string)) {
        if (!length(formula_or_string) || !nzchar(formula_or_string[[1]])) {
            return(stats::as.formula("~ 1"))
        }
        return(stats::as.formula(formula_or_string[[1]]))
    }
    if (!inherits(formula_or_string, "formula")) {
        stop("`design_formula` must be a formula or character string.", call. = FALSE)
    }
    formula_or_string
}

.omicsgmf_estimate_ncomponents <- function(
  sce,
  X,
  exprs_name,
  family,
  max_rank,
  rank_args
) {
    rank_call <- modifyList(
        list(object = sce, X = X, exprs_values = exprs_name, family = family, maxcomp = max_rank),
        rank_args
    )
    ranked_sce <- do.call(omicsGMF::runRankGMF, rank_call)
    list(
        sce = ranked_sce,
        ncomponents = .omicsgmf_infer_rank(ranked_sce)
    )
}

.omicsgmf_infer_rank <- function(sce) {
    meta <- S4Vectors::metadata(sce)
    if (is.null(meta) || !length(meta)) {
        return(NULL)
    }
    gmf_meta <- meta$GMF %||% meta$gmf
    if (is.null(gmf_meta)) {
        return(NULL)
    }
    rank_info <- gmf_meta$rank %||% gmf_meta$Rank %||% gmf_meta$rank_results
    if (is.null(rank_info)) {
        return(NULL)
    }

    # 1) Check named scalars with "opt" in the name
    if (is.list(rank_info)) {
        for (nm in names(rank_info)) {
            if (!length(nm) || !grepl("opt", nm, ignore.case = TRUE)) {
                next
            }
            candidate <- rank_info[[nm]]
            if (is.numeric(candidate) && length(candidate) == 1L && is.finite(candidate)) {
                return(as.integer(round(candidate)))
            }
        }
    }

    # 2) Look for attributes storing an optimal rank
    attrs <- attributes(rank_info)
    if (length(attrs)) {
        for (nm in names(attrs)) {
            if (!grepl("opt", nm, ignore.case = TRUE)) {
                next
            }
            candidate <- attrs[[nm]]
            if (is.numeric(candidate) && length(candidate) == 1L && is.finite(candidate)) {
                return(as.integer(round(candidate)))
            }
        }
    }

    # 3) Data-frame like summaries
    if (is.data.frame(rank_info)) {
        col_opts <- c("opt", "optimal", "best")
        for (col in col_opts) {
            mt <- grep(col, names(rank_info), ignore.case = TRUE, value = TRUE)
            if (!length(mt)) next
            candidate <- rank_info[[mt[[1]]]]
            if (is.numeric(candidate) && length(candidate) >= 1L) {
                idx <- if ("deviance" %in% names(rank_info)) {
                    which.min(rank_info$deviance)
                } else {
                    1L
                }
                return(as.integer(round(candidate[idx])))
            }
        }
        if (all(c("components", "deviance") %in% names(rank_info))) {
            idx <- which.min(rank_info$deviance)
            candidate <- rank_info$components[idx]
            if (is.numeric(candidate) && length(candidate) == 1L) {
                return(as.integer(round(candidate)))
            }
        }
    }

    # 4) Named numeric vector (e.g., deviances keyed by component)
    if (is.numeric(rank_info) && !is.null(names(rank_info))) {
        idx <- which.min(rank_info)
        candidate <- suppressWarnings(as.integer(names(rank_info)[idx]))
        if (!is.na(candidate) && candidate > 0L) {
            return(candidate)
        }
    }

    NULL
}
