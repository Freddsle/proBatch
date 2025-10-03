# Additional methods for batch effect correction
# Include: BERT


#' @title BERT-based batch correction (unified; long/matrix)
#' @description Adjusts batch effects using BERT (Batch-Effect Reduction Trees),
#' which aligns batches and tolerates NAs.
#' Works on long data.frames (\code{format="long"}) and wide matrices
#' (\code{format="wide"}) by converting to/from matrices under the hood.
#'
#' @inheritParams correct_with_ComBat
#' @param format One of \code{"long"} or \code{"wide"}.
#' @param bert_method Edge adjustment inside BERT: one of \code{"ComBat"}, \code{"limma"}, \code{"ref"}.
#'   Defaults to \code{"ComBat"}.
#' @param combatmode Integer in \code{1:4} passed to BERT (parametric/non-parametric, mean.only). Default 1.
#'   See \code{?BERT::BERT}. Ignored if \code{bert_method != "ComBat"}.
#' @param cores Integer; number of workers for BERT. If \code{NULL} (default) BERT uses \code{BiocParallel::bpparam()}.
#' @param BPPARAM Optional \code{BiocParallelParam} object for BERT.
#' @param keep_all For long output, columns retained (as in other correct_* functions).
#' @details
#' BERT itself requires a \emph{samples × features} table; proBatch's wrapper builds it from the input
#' \emph{features × samples} matrix. Missing values (NAs) are allowed.
#' @return If \code{format="wide"}, a corrected numeric matrix (features×samples).
#'   If \code{format="long"}, the original long data with corrected \code{measure_col} and
#'   a \code{preBatchCorr_[measure_col]} column preserved.
#' @references
#'   BERT Bioconductor page & manual. Schumann Y. & Schlumbohm S., 2025.
#'   \doi{10.1038/s41467-025-62237-4}. Nature Communications (2025)
#' @seealso \code{\link{correct_batch_effects}}, \code{\link{BERT::BERT}}
#' @export
correct_with_BERT <- function(
    x, sample_annotation,
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    format = c("long", "wide"),
    bert_method = c("ComBat", "limma", "ref"),
    combatmode = 1,
    covariates_cols = NULL,
    keep_all = "default",
    cores = NULL,
    BPPARAM = NULL,
    ...) {
    .pb_requireNamespace("BERT")
    format <- match.arg(format)
    bert_method <- match.arg(bert_method)

    if (identical(format, "wide")) {
        if (!is.matrix(x)) stop("format='wide' requires a numeric matrix.")
        corrected_matrix <- .bert_matrix_step(
            data_matrix       = x,
            sample_annotation = sample_annotation,
            sample_id_col     = sample_id_col,
            batch_col         = batch_col,
            covariates_cols   = covariates_cols,
            bert_method       = bert_method,
            combatmode        = combatmode,
            cores             = cores,
            BPPARAM           = BPPARAM,
            ...
        )
        return(corrected_matrix)
    }

    # LONG -> matrix
    if (!is.data.frame(x)) stop("format='long' requires a data.frame.")
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
        measure_col    = measure_col,
        sample_id_col  = sample_id_col
    )

    corrected_matrix <- .bert_matrix_step(
        data_matrix       = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col     = sample_id_col,
        batch_col         = batch_col,
        covariates_cols   = covariates_cols,
        bert_method       = bert_method,
        combatmode        = combatmode,
        cores             = cores,
        BPPARAM           = BPPARAM,
        ...
    )

    corrected_df <- matrix_to_long(
        corrected_matrix,
        feature_id_col = feature_id_col,
        measure_col    = measure_col,
        sample_id_col  = sample_id_col
    )

    # preserve original values in preBatchCorr_*
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

.run_BERT_core <- function(
    data_matrix, # features × samples (numeric)
    sample_annotation, # data.frame
    sample_id_col,
    batch_col,
    covariates_cols = NULL,
    bert_method = c("ComBat", "limma", "ref"),
    combatmode = 1,
    cores = NULL,
    BPPARAM = NULL,
    ...) {
    bert_method <- match.arg(bert_method)

    # Align SA to matrix columns (samples)
    sample_annotation <- .align_sample_annotation(
        sample_annotation,
        sample_ids    = colnames(data_matrix),
        sample_id_col = sample_id_col
    )

    # Build samples×features df for BERT, appending SA columns
    # Keep a strict list of feature columns to extract back.
    feature_cols <- rownames(data_matrix) # features (rows of input matrix)
    stopifnot(!is.null(feature_cols))
    df_t <- t(data_matrix) # samples × features
    storage.mode(df_t) <- "double"

    # Ensure SA contains required columns
    if (!(batch_col %in% names(sample_annotation))) {
        stop("Batch column is not present in sample_annotation")
    }

    # Tranform covariates_cols to design format (BERT requires numeric values in covariates)
    if (!is.null(covariates_cols) && length(covariates_cols)) {
        missing_cov <- setdiff(covariates_cols, names(sample_annotation))
        if (length(missing_cov)) {
            stop("Covariates missing in sample_annotation: ", paste(missing_cov, collapse = ", "))
        }
        covariates <- as.data.frame(sample_annotation[, covariates_cols, drop = FALSE])
        mod <- model.matrix(~., data = covariates)
        # replace original covariates with design matrix columns (excluding intercept)
        orig_covariates <- covariates_cols
        if (ncol(mod) > 1) {
            new_covariates <- colnames(mod)[-1]
            sample_annotation <- cbind(
                sample_annotation[, setdiff(names(sample_annotation), orig_covariates), drop = FALSE],
                as.data.frame(mod[, -1, drop = FALSE])
            )
            covariates_cols <- new_covariates
        }
    }

    # Minimal BERT metadata columns; rely on 'batchname', 'samplename', 'covariatename' args
    meta_cols <- c(
        sample_id_col, batch_col,
        if (!is.null(covariates_cols)) covariates_cols else character(0)
    )
    meta_df <- as.data.frame(sample_annotation[, meta_cols, drop = FALSE])

    # Preserve row order (samples) and make a clean data.frame for BERT
    stopifnot(nrow(meta_df) == nrow(df_t))
    df_bert <- cbind(meta_df, as.data.frame(df_t, check.names = TRUE))
    rownames(df_bert) <- df_bert[[sample_id_col]]
    rownames_order <- rownames(df_bert)
    df_bert[[sample_id_col]] <- NULL

    # Call BERT; allow user to pass parallel settings & extra args
    res <- BERT::BERT(
        data = df_bert,
        method = bert_method,
        combatmode = combatmode,
        batchname = batch_col,
        samplename = "Sample", # sample IDs were removed
        covariatename = covariates_cols,
        cores = cores,
        BPPARAM = BPPARAM,
        verify = TRUE,
        ...
    )
    # Return the sample_id_col if it was removed
    if (!is.null(sample_id_col) && !(sample_id_col %in% names(res))) {
        res[[sample_id_col]] <- rownames(res)
        res <- res[rownames_order, , drop = FALSE]
    }

    # Result mirrors input: samples×features (+ meta); select only the original feature columns
    res_df <- as.data.frame(res, check.names = FALSE)
    missing_back <- setdiff(feature_cols, colnames(res_df))
    if (length(missing_back)) {
        stop("BERT output is missing expected features: ", paste(missing_back, collapse = ", "))
    }

    corrected_t <- as.matrix(res_df[, feature_cols, drop = FALSE]) # samples×features
    storage.mode(corrected_t) <- "double"
    t(corrected_t) # features×samples
}

.bert_matrix_step <- function(
    data_matrix, sample_annotation,
    sample_id_col, batch_col,
    covariates_cols = NULL,
    bert_method = c("ComBat", "limma", "ref"),
    combatmode = 1,
    cores = NULL,
    BPPARAM = NULL,
    ...) {
    # Coerce to numeric matrix (do NOT impute; BERT tolerates NA)
    # But check that after NA removal it's still numeric and has >=2 features
    if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
    }
    if (!is.numeric(data_matrix)) {
        stop("Input must be coercible to a numeric matrix for BERT correction.")
    }
    storage.mode(data_matrix) <- "double"
    if (sum(rowSums(!is.na(data_matrix)) > 0L) < 2L) {
        stop("BERT requires at least two not-NA features.")
    }

    .run_BERT_core(
        data_matrix       = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col     = sample_id_col,
        batch_col         = batch_col,
        covariates_cols   = covariates_cols,
        bert_method       = match.arg(bert_method),
        combatmode        = combatmode,
        cores             = cores,
        BPPARAM           = BPPARAM,
        ...
    )
}
