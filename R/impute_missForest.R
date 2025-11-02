## ------------------------------------------------------------------------
## impute_missForest.R
##
## Optional missForest-based imputation for proBatch / ProBatchFeatures
## ------------------------------------------------------------------------

#' Impute missing intensities with missForest (optional)
#'
#' Adapter around **`missForest::missForest()`** that lets proBatch users run
#' non-parametric random-forest imputation on:
#'
#' - long data.frames (one feature per row and sample),
#' - wide data matrices (features in rows, samples in columns),
#' - and, via [pb_transform()], `ProBatchFeatures` pipelines (step name:
#'   `"missForestImpute"`).
#'
#' The original missForest implementation expects **observations in rows and
#' variables in columns**; proteomics matrices in proBatch are **features in rows
#' and samples in columns**. This adapter therefore **transposes** before
#' calling `missForest()` and **transposes back** afterward. This is required by
#' the upstream package design and cannot be avoided. See the missForest manual
#' for details. (Stekhoven & Bühlmann 2012; CRAN missForest pdf).
#'
#' @param x Object to impute; one of:
#'   \itemize{
#'     \item long data.frame (proBatch long format, e.g. output of [proBatch()]);
#'     \item numeric matrix/data.frame with features in rows and samples in columns;
#'     \item for pipelines, call through [pb_transform()] instead.
#'   }
#' @param sample_annotation Optional sample annotation table; must contain
#'   `sample_id_col` if provided. Kept for API symmetry with other methods.
#' @param sample_id_col Column in `sample_annotation` that matches matrix column
#'   names or long-format sample IDs. Default `"FullRunName"`.
#' @param feature_id_col Column in long-format data holding peptide/protein IDs.
#'   Default `"feature_id"`.
#' @param measure_col Column with intensities in long-format data.
#'   Default `"Intensity"`.
#' @param ... Extra arguments forwarded *verbatim* to
#'   [missForest::missForest()], e.g. `maxiter`, `ntree`, `mtry`, `parallelize`,
#'   `backend`, `num.threads`, `replace`, `classwt`, etc. See upstream docs.
#'
#' @return Object of the same class as `x`, with imputed values.
#'
#' @references
#' Stekhoven DJ, Bühlmann P. MissForest—non-parametric missing value imputation for mixed-type data.
#' Bioinformatics. 2012 Jan 1;28(1):112-8.
#' DOI:  https://doi.org/10.1093/bioinformatics/btr597.
#'
#' @rdname imputeMissForest
#' @export
imputeMissForest <- function(x,
                             sample_annotation = NULL,
                             sample_id_col = "FullRunName",
                             feature_id_col = "feature_id",
                             measure_col = "Intensity",
                             ...) {
    .pb_requireNamespace("missForest")

    ## --- wide matrix / data.frame path ------------------------------------
    if (is.matrix(x)) {
        return(imputeMissForest_dm(
            x = x,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            ...
        ))
    }

    if (is.data.frame(x)) {
        cols <- names(x)
        looks_long <- all(c(feature_id_col, sample_id_col, measure_col) %in% cols)

        if (looks_long) {
            return(imputeMissForest_df(
                x = x,
                sample_annotation = sample_annotation,
                sample_id_col = sample_id_col,
                feature_id_col = feature_id_col,
                measure_col = measure_col,
                ...
            ))
        }

        ## not clearly long -> treat as wide
        return(imputeMissForest_dm(
            x = x,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            ...
        ))
    }

    ## --- pipelines should go through pb_transform() -----------------------
    if (inherits(x, "ProBatchFeatures")) {
        return(imputeMissForest.ProBatchFeatures(
            x = x,
            sample_id_col = sample_id_col,
            ...
        ))
    }

    if (inherits(x, "QFeatures")) {
        stop(
            "imputeMissForest(): QFeatures inputs are not supported directly. ",
            "Use pb_transform(..., steps = 'missForestImpute') to apply missForest within a pipeline.",
            call. = FALSE
        )
    }

    stop("imputeMissForest(): unsupported object of class ",
        paste(class(x), collapse = "/"),
        call. = FALSE
    )
}

# ---------------------------------------------------------------------
# long-format front-end
# ---------------------------------------------------------------------

#' @rdname imputeMissForest
#' @export
imputeMissForest_df <- function(x,
                                sample_annotation = NULL,
                                sample_id_col = "FullRunName",
                                feature_id_col = "feature_id",
                                measure_col = "Intensity",
                                ...) {
    .pb_requireNamespace("missForest")

    ## 1) long -> wide matrix (features x samples)
    data_matrix <- long_to_matrix(
        df_long        = x,
        feature_id_col = feature_id_col,
        sample_id_col  = sample_id_col,
        measure_col    = measure_col
    )

    ## 2) run missForest on matrix
    imputed_matrix <- .missforest_matrix_step(
        data_matrix     = data_matrix,
        missforest_args = list(...)
    )

    ## 3) wide -> long again
    matrix_to_long(
        data_matrix = imputed_matrix,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col
    )
}

# ---------------------------------------------------------------------
# wide-format front-end
# ---------------------------------------------------------------------

#' @rdname imputeMissForest
#' @export
imputeMissForest_dm <- function(x,
                                sample_annotation = NULL,
                                sample_id_col = "FullRunName",
                                ...) {
    .pb_requireNamespace("missForest")

    data_matrix <- if (is.matrix(x)) x else as.matrix(x)

    .missforest_matrix_step(
        data_matrix     = data_matrix,
        missforest_args = list(...)
    )
}

# ---------------------------------------------------------------------
# ProBatchFeatures front-end (to be called via pb_transform)
# ---------------------------------------------------------------------

#' @rdname imputeMissForest
#' @export
imputeMissForest.ProBatchFeatures <- function(
  x,
  sample_id_col = "FullRunName",
  pbf_name = NULL,
  final_name = NULL,
  ...
) {
    object <- x
    .pb_requireNamespace("missForest")

    if (is.null(pbf_name)) {
        pbf_name <- pb_current_assay(object)
        message("`pbf_name` not provided, using the most recent assay: ", pbf_name)
    }

    extra <- list(...)
    mf_args <- list()
    if (length(extra) && "missforest_args" %in% names(extra)) {
        supplied <- extra[["missforest_args"]]
        if (length(supplied)) {
            mf_args <- supplied
        }
        extra[["missforest_args"]] <- NULL
    }
    if (length(extra)) {
        mf_args <- modifyList(mf_args, extra)
    }
    params <- if (length(mf_args)) list(missforest_args = mf_args) else list()

    step_name <- "missForestImpute"
    step_label <- .pb_step_label(step_name, params)

    out <- .pb_apply_step(
        object = object,
        from   = pbf_name,
        step   = step_label,
        fun    = step_name,
        params = params
    )

    object <- out$object
    assay_new <- out$assay

    if (!is.null(final_name) && !is.null(assay_new) && assay_new %in% names(object)) {
        names(object)[match(assay_new, names(object))] <- final_name
    }

    object
}

# ---------------------------------------------------------------------
# actual step function used inside pb_transform()
# ---------------------------------------------------------------------

#' ProBatchFeatures step: missForest imputation
#'
#' This is the function registered under `"missForestImpute"`. It receives an
#' assay matrix (features in rows, samples in columns) and returns the same
#' matrix shape with missForest-imputed values.
#'
#' @param x numeric matrix, features in rows, samples in columns.
#' @param ... forwarded to [missForest::missForest()].
#'
#' @return numeric matrix with imputed values.
#'
#' @keywords internal
#' @export
missForestImpute <- function(x, ...) {
    .pb_requireNamespace("missForest")

    if (!is.matrix(x)) {
        x <- as.matrix(x)
    }

    dots <- list(...)
    if (length(dots) == 1L && !is.null(dots$missforest_args)) {
        mf_args <- dots$missforest_args
    } else {
        mf_args <- dots
    }

    .missforest_matrix_step(
        data_matrix     = x,
        missforest_args = mf_args
    )
}

# ---------------------------------------------------------------------
# shared internal helper
# ---------------------------------------------------------------------

#' @keywords internal
.missforest_matrix_step <- function(data_matrix,
                                    missforest_args = list()) {
    .pb_requireNamespace("missForest")

    # 0) accept the pipeline's "list(missforest_args = ...)" shape
    if (length(missforest_args) == 1L &&
        !is.null(missforest_args$missforest_args)) {
        missforest_args <- missforest_args$missforest_args
    }

    if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
    }
    storage.mode(data_matrix) <- "double"

    if (is.null(colnames(data_matrix))) {
        stop("missForest imputation requires matrix column names (sample identifiers).",
            call. = FALSE
        )
    }

    ## preserve names
    feat_ids <- rownames(data_matrix)
    sample_ids <- colnames(data_matrix)

    ## --------------------------------------------------------------
    ## 1) protect missForest from degenerate rows / cols
    ##    missForest (per docs) needs at least one observed value per
    ##    variable; pure-NA columns/rows must be removed first.
    ## --------------------------------------------------------------
    all_na_rows <- which(rowSums(!is.na(data_matrix)) == 0L)
    all_na_cols <- which(colSums(!is.na(data_matrix)) == 0L)

    core_mat <- data_matrix
    if (length(all_na_rows)) {
        core_mat <- core_mat[-all_na_rows, , drop = FALSE]
    }
    if (length(all_na_cols)) {
        core_mat <- core_mat[, -all_na_cols, drop = FALSE]
    }

    if (!nrow(core_mat) || !ncol(core_mat)) {
        stop("missForest(): no data left to impute after removing all-NA rows/columns.",
            call. = FALSE
        )
    }

    ## --------------------------------------------------------------
    ## 2) transpose to samples x features (missForest layout)
    ## --------------------------------------------------------------
    core_df <- as.data.frame(t(core_mat), stringsAsFactors = FALSE)

    ## --------------------------------------------------------------
    ## 3) call missForest with user-supplied args
    ##    upstream returns a list with $ximp, $OOBerror, $error, ...
    ## --------------------------------------------------------------
    mf_call <- c(list(xmis = core_df), missforest_args)
    mf_res <- do.call(missForest::missForest, mf_call)

    imputed_df <- mf_res$ximp

    ## --------------------------------------------------------------
    ## 4) transpose back to features x samples
    ## --------------------------------------------------------------
    imputed_core <- t(as.matrix(imputed_df))
    storage.mode(imputed_core) <- "double"

    rownames(imputed_core) <- rownames(core_mat)
    colnames(imputed_core) <- colnames(core_mat)

    ## --------------------------------------------------------------
    ## 5) reinsert removed all-NA rows/cols
    ##    (we keep them NA – we cannot fabricate signal)
    ## --------------------------------------------------------------
    if (length(all_na_rows) || length(all_na_cols)) {
        full_mat <- matrix(
            NA_real_,
            nrow = nrow(data_matrix),
            ncol = ncol(data_matrix),
            dimnames = list(feat_ids, sample_ids)
        )

        ## place the imputed core
        r_idx <- seq_len(nrow(data_matrix))
        c_idx <- seq_len(ncol(data_matrix))
        if (length(all_na_rows)) {
            r_idx <- r_idx[-all_na_rows]
        }
        if (length(all_na_cols)) {
            c_idx <- c_idx[-all_na_cols]
        }
        full_mat[r_idx, c_idx] <- imputed_core

        imputed <- full_mat

        warning(
            "missForest(): removed and reinserted ",
            length(all_na_rows), " all-NA rows and ",
            length(all_na_cols), " all-NA columns; these positions remain NA."
        )
    } else {
        imputed <- imputed_core
    }

    ## attach OOB diagnostics (does not break API)
    attr(imputed, "missForest_OOBerror") <- mf_res$OOBerror
    attr(imputed, "missForest_error") <- mf_res$error

    imputed
}
