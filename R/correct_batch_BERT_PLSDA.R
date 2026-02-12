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
#' @param combatmode Integer in \code{1:4} passed to BERT (parametric/non-parametric, mean.only).
#'   If \code{NULL}, the wrapper does not pass it and BERT uses its own default.
#'   Ignored if \code{bert_method != "ComBat"}.
#' @param cores Integer; number of workers for BERT. If \code{NULL} and \code{BPPARAM} is also \code{NULL},
#'   the wrapper leaves worker selection to \code{BERT::BERT()}. Provide \code{cores} (and optionally a matching
#'   \code{BPPARAM}) to control parallel execution (see \code{?BiocParallel}).
#' @param BPPARAM Optional \code{BiocParallelParam} object for BERT. When \code{NULL}, the wrapper does not
#'   set a default and \code{BERT::BERT()} uses its own default.
#' @param ... Additional arguments passed to \code{BERT::BERT()}. Note that by default the wrapper sets
#'   \code{referencename = " "} (a single space) unless you override it via \code{...}.
#' @param keep_all For long output, columns retained (as in other correct_* functions).
#'
#' @details
#' BERT itself requires a \emph{samples × features} table; proBatch's wrapper builds it from the input
#' \emph{features × samples} matrix. Missing values (NAs) are allowed.
#'
#' When invoked through \code{pb_transform()}, the resulting pipeline step is tagged with the chosen
#' \code{bert_method}: \code{BERTc} (ComBat), \code{BERTl} (limma), or \code{BERTr} (ref). This keeps
#' method-specific corrections distinct when multiple BERT variants are stored on the same object.
#'
#' @return If \code{format="wide"}, a corrected numeric matrix (features×samples).
#'   If \code{format="long"}, the original long data with corrected \code{measure_col} and
#'   a \code{preBatchCorr_[measure_col]} column preserved.
#'
#' @references
#'   BERT Bioconductor page & manual. Schumann Y. & Schlumbohm S., 2025.
#'   \doi{10.1038/s41467-025-62237-4}. Nature Communications (2025)
#'
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
  combatmode = NULL,
  covariates_cols = NULL,
  keep_all = "default",
  cores = NULL,
  BPPARAM = NULL,
  ...
) {
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
    prep <- .pb_prepare_long_matrix(
        df_long = x,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        batch_col = batch_col,
        error_message = "format='long' requires a data.frame."
    )
    df_long <- prep$df_long
    sample_annotation <- prep$sample_annotation
    data_matrix <- prep$data_matrix
    original_cols <- prep$original_cols

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

    .post_correction_to_long(
        corrected_matrix, df_long,
        feature_id_col, measure_col, sample_id_col,
        original_cols, keep_all
    )
}

.coerce_covariate_to_numeric <- function(cov_values) {
    if (is.null(cov_values)) {
        return(cov_values)
    }
    if (is.logical(cov_values)) {
        return(as.numeric(cov_values))
    }
    if (is.factor(cov_values)) {
        return(as.numeric(cov_values))
    }
    if (inherits(cov_values, c("Date", "POSIXct", "POSIXlt"))) {
        return(as.numeric(cov_values))
    }
    if (!is.numeric(cov_values)) {
        fac <- factor(cov_values, exclude = NULL)
        return(as.numeric(fac))
    }
    as.numeric(cov_values)
}

.coerce_batch_label_to_numeric <- function(values) {
    if (is.null(values)) {
        return(values)
    }
    as.numeric(as.factor(values))
}

.sanitize_covariates_for_BERT <- function(sample_annotation, covariates_cols) {
    if (is.null(covariates_cols) || !length(covariates_cols)) {
        return(list(
            sample_annotation = sample_annotation,
            covariates_cols = covariates_cols
        ))
    }

    if (!is.character(covariates_cols)) {
        stop("`covariates_cols` must be a character vector of column names.")
    }

    missing_cov <- setdiff(covariates_cols, names(sample_annotation))
    if (length(missing_cov)) {
        stop(
            "Covariates missing in sample_annotation: ",
            paste(missing_cov, collapse = ", ")
        )
    }

    for (covariate in covariates_cols) {
        cov_values <- sample_annotation[[covariate]]

        if (is.null(cov_values)) {
            stop("Covariate '", covariate, "' unexpectedly NULL.")
        }

        sample_annotation[[covariate]] <- .coerce_covariate_to_numeric(cov_values)
    }

    list(
        sample_annotation = sample_annotation,
        covariates_cols = covariates_cols
    )
}

.run_BERT_core <- function(
  data_matrix, # features × samples (numeric)
  sample_annotation, # data.frame
  sample_id_col,
  batch_col,
  covariates_cols = NULL,
  bert_method = c("ComBat", "limma", "ref"),
  combatmode = NULL,
  cores = NULL,
  BPPARAM = NULL,
  reference_name = " ",
  ...
) {
    bert_method <- match.arg(bert_method)

    if (is.null(sample_annotation)) {
        stop("sample_annotation must be provided for batch correction")
    }
    sample_annotation <- as.data.frame(sample_annotation)

    if (!(batch_col %in% names(sample_annotation))) {
        stop("Batch column is not present in sample_annotation")
    }

    # Align SA to matrix columns (samples)
    sample_annotation <- .align_sample_annotation(
        sample_annotation,
        sample_ids    = colnames(data_matrix),
        sample_id_col = sample_id_col
    )

    dots <- list(...)
    if (!is.null(dots$referencename)) {
        if (missing(reference_name) || identical(reference_name, " ")) {
            reference_name <- dots$referencename
        }
        dots$referencename <- NULL
    }
    labelname <- dots$labelname
    samplename <- dots$samplename
    if (!is.null(dots$samplename)) {
        dots$samplename <- NULL
    }
    if (is.null(samplename) || !nzchar(samplename)) {
        samplename <- "Sample"
    }
    if (identical(samplename, "Sample")) {
        sample_annotation[[samplename]] <- colnames(data_matrix)
    }
    if (!is.null(samplename) && !(samplename %in% names(sample_annotation))) {
        stop("samplename column is not present in sample_annotation: ", samplename)
    }
    if (!is.null(labelname) && !(labelname %in% names(sample_annotation))) {
        stop("labelname column is not present in sample_annotation: ", labelname)
    }
    sample_annotation[[batch_col]] <- .coerce_batch_label_to_numeric(
        sample_annotation[[batch_col]]
    )
    if (!is.null(labelname)) {
        sample_annotation[[labelname]] <- .coerce_batch_label_to_numeric(
            sample_annotation[[labelname]]
        )
    }

    # Build samples×features df for BERT, appending SA columns
    # Keep a strict list of feature columns to extract back.
    feature_cols <- rownames(data_matrix) # features (rows of input matrix)
    stopifnot(!is.null(feature_cols))
    df_t <- t(data_matrix) # samples × features
    storage.mode(df_t) <- "double"

    # Tranform covariates_cols to design format (BERT requires numeric values in covariates)
    if (!is.null(covariates_cols) && length(covariates_cols) > 0L) {
        prepared <- .sanitize_covariates_for_BERT(sample_annotation, covariates_cols)
        sample_annotation <- prepared$sample_annotation
        covariates_cols <- prepared$covariates_cols
    }

    # Early guard: BERT requires >=2 samples per batch×covariate level
    if (!is.null(covariates_cols) && length(covariates_cols) > 0L) {
        check_df <- sample_annotation[, c(batch_col, covariates_cols), drop = FALSE]
        # Treat NAs as an explicit level to avoid silent drops
        for (nm in names(check_df)) {
            if (anyNA(check_df[[nm]])) {
                check_df[[nm]] <- as.character(check_df[[nm]])
                check_df[[nm]][is.na(check_df[[nm]])] <- "<NA>"
            }
        }
        tab <- table(check_df, useNA = "ifany")
        if (any(tab < 2L)) {
            bad_n <- sum(tab < 2L)
            warning(
                "BERT pre-check: found ", bad_n,
                " batch×covariate strata with <2 samples; stopping early."
            )
            stop("Not enough samples at batch/covariate level.", call. = FALSE)
        }
    }

    # Minimal BERT metadata columns; rely on 'batchname', 'samplename', 'covariatename' args
    meta_cols <- unique(c(
        samplename,
        batch_col,
        if (!is.null(covariates_cols)) covariates_cols else character(0),
        if (!is.null(labelname)) labelname else character(0)
    ))
    meta_df <- as.data.frame(sample_annotation[, meta_cols, drop = FALSE])

    # Preserve row order (samples) and make a clean data.frame for BERT
    stopifnot(nrow(meta_df) == nrow(df_t))
    df_bert <- cbind(meta_df, as.data.frame(df_t, check.names = TRUE))
    sample_order <- df_bert[[samplename]]
    rownames(df_bert) <- seq_len(nrow(df_bert))

    if (!is.null(BPPARAM) || !is.null(cores)) {
        .pb_requireNamespace("BiocParallel")
    }
    if (!is.null(BPPARAM) && is.null(cores)) {
        cores <- BiocParallel::bpworkers(BPPARAM)
        if (is.null(cores) || is.na(cores) || cores < 1L) {
            cores <- 1L
        }
    }

    # Call BERT; allow user to pass parallel settings & extra args
    bert_args <- list(
        data = df_bert,
        method = bert_method,
        batchname = batch_col,
        samplename = samplename,
        covariatename = covariates_cols
    )
    if (!is.null(combatmode)) {
        bert_args$combatmode <- combatmode
    }
    if (!is.null(reference_name)) {
        bert_args$referencename <- reference_name
    }
    if (!is.null(cores)) {
        bert_args$cores <- cores
    }
    if (!is.null(BPPARAM)) {
        bert_args$BPPARAM <- BPPARAM
    }
    res <- do.call(BERT::BERT, c(bert_args, dots))
    if (!is.null(samplename) && (samplename %in% names(res))) {
        res <- res[match(sample_order, res[[samplename]]), , drop = FALSE]
    }
    rownames(res) <- sample_order
    # Return the sample_id_col if it was removed
    if (!is.null(sample_id_col) && !(sample_id_col %in% names(res))) {
        res[[sample_id_col]] <- sample_order
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
  combatmode = NULL,
  cores = NULL,
  BPPARAM = NULL,
  reference_name = " ",
  ...
) {
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
        reference_name    = reference_name,
        ...
    )
}

#' @title Batch correction via (s)PLSDA-batch (optional; requires PLSDAbatch)
#' @description Wrapper around \code{PLSDAbatch::PLSDA_batch()} to remove batch variation
#'   while preserving treatment/biological variation. Uses PLSDA-batch by default;
#'   switches to sPLSDA-batch when \code{keepX_trt} is supplied or when
#'   \code{run_splsda=TRUE} and \code{keepX_trt} is tuned internally.
#'   Internally works on matrices; long/PBF are converted to/from matrices.
#' @inheritParams correct_with_ComBat
#' @param effect_col Column in \code{sample_annotation} with biological group (treatment).
#'   If \code{NULL}, \code{Y.trt} is omitted (PLSDA-batch cannot preserve treatment in that case).
#' @param ncomp_trt,ncomp_bat Integers: # of treatment/batch components.
#' @param keepX_trt Optional numeric vector of length \code{ncomp_trt};
#'   triggers sPLSDA-batch (variables kept per treatment component).
#' @param keepX_bat Optional numeric vector of length \code{ncomp_bat} (usually leave default = all).
#' @param balance logical or "auto". If "auto" (default), we set
#'   \code{FALSE} when any \code{table(batch, treatment)} cell is zero (unbalanced design),
#'   else \code{TRUE} (as recommended by PLSDAbatch).
#' @param near_zero_var Passed to PLSDAbatch; keep \code{TRUE} for many zeros.
#' @param format Either \code{"wide"} (matrix, features x samples) or \code{"long"} (data.frame).
#' @param fill_the_missing Missing-value policy before invoking PLSDA-batch. If \code{NULL},
#'   missing values are left as-is (and the method will fail if NAs remain). Set \code{FALSE}
#'   to keep NA entries untouched, or supply a numeric value to impute prior to correction.
#' @param max.iter,tol Passed to PLSDAbatch; max iterations and tolerance for convergence.
#' @param run_splsda Logical; if \code{TRUE} runs sPLSDA correction.
#' @return A batch-corrected matrix (features x samples) for \code{format="wide"};
#'   for \code{format="long"} returns a long data.frame with \code{measure_col} corrected.
#' @seealso \link[PLSDAbatch]{PLSDA_batch} for details.
#' @references Yiwen Wang, Kim-Anh Lê Cao, PLSDA-batch: a multivariate framework to correct for batch
#' effects in microbiome data, Briefings in Bioinformatics, Volume 24, Issue 2, March 2023, bbac622,
#' https://doi.org/10.1093/bib/bbac622
#' @export
correct_with_PLSDA_batch <- function(x,
                                     sample_annotation = NULL,
                                     sample_id_col = "FullRunName",
                                     feature_id_col = "peptide_group_label",
                                     measure_col = "Intensity",
                                     batch_col = "MS_batch",
                                     effect_col = NULL,
                                     ncomp_trt = NULL,
                                     ncomp_bat = NULL,
                                     keepX_trt = NULL,
                                     keepX_bat = NULL,
                                     balance = "auto",
                                     near_zero_var = TRUE,
                                     max.iter = 500,
                                     tol = 1e-06,
                                     format = c("wide", "long"),
                                     fill_the_missing = NULL,
                                     keep_all = "default",
                                     run_splsda = FALSE,
                                     ...) {
    .pb_requireNamespace("PLSDAbatch")
    format <- match.arg(format)

    # Clear early error & namespace hint
    if (is.null(sample_annotation)) {
        stop(
            "sample_annotation is required for PLSDA-batch (to supply ", batch_col,
            if (!is.null(effect_col)) paste0(" and ", effect_col) else "", ")."
        )
    }

    if (is.character(balance)) {
        if (!identical(balance, "auto")) {
            stop("balance must be logical or \"auto\".")
        }
    } else if (!is.logical(balance) || length(balance) != 1L) {
        stop("balance must be a single logical or \"auto\".")
    }

    if (identical(format, "wide")) {
        if (!is.matrix(x)) stop("format='wide' requires a numeric matrix.")
        corrected_matrix <- .plsda_matrix_step(
            data_matrix       = x,
            sample_annotation = sample_annotation,
            sample_id_col     = sample_id_col,
            batch_col         = batch_col,
            effect_col        = effect_col,
            ncomp_trt         = ncomp_trt,
            ncomp_bat         = ncomp_bat,
            keepX_trt         = keepX_trt,
            keepX_bat         = keepX_bat,
            balance           = balance,
            near_zero_var     = near_zero_var,
            max.iter          = max.iter,
            tol               = tol,
            run_splsda        = run_splsda,
            fill_the_missing  = fill_the_missing,
            ...
        )
        return(corrected_matrix)
    }

    # LONG -> matrix
    prep <- .pb_prepare_long_matrix(
        df_long = x,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        batch_col = batch_col,
        fill_the_missing = fill_the_missing,
        warning_message = "PLSDA-batch cannot operate with missing values in the matrix",
        error_message = "format='long' requires a data.frame."
    )
    df_long <- prep$df_long
    sample_annotation <- prep$sample_annotation
    data_matrix <- prep$data_matrix
    original_cols <- prep$original_cols

    corrected_matrix <- .plsda_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        effect_col = effect_col,
        ncomp_trt = ncomp_trt,
        ncomp_bat = ncomp_bat,
        keepX_trt = keepX_trt,
        keepX_bat = keepX_bat,
        balance = balance,
        near_zero_var = near_zero_var,
        max.iter = max.iter,
        tol = tol,
        run_splsda = run_splsda,
        fill_the_missing = fill_the_missing,
        ...
    )

    .post_correction_to_long(
        corrected_matrix, df_long,
        feature_id_col, measure_col, sample_id_col,
        original_cols, keep_all
    )
}

.splsda_matrix_step <- function(run_splsda = TRUE, ...) {
    .plsda_matrix_step(run_splsda = run_splsda, ...)
}


.plsda_matrix_step <- function(
  data_matrix, sample_annotation = NULL,
  sample_id_col = "FullRunName",
  batch_col = "MS_batch",
  effect_col = NULL,
  ncomp_trt = NULL,
  ncomp_bat = NULL,
  keepX_trt = NULL,
  keepX_bat = NULL,
  balance = "auto",
  near_zero_var = TRUE,
  max.iter = 500,
  tol = 1e-06,
  run_splsda = FALSE,
  fill_the_missing = NULL,
  ...
) {
    .run_matrix_method(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        fill_the_missing = fill_the_missing,
        missing_warning = "PLSDA-batch cannot operate with missing values in the matrix",
        method_fun = function(data_matrix, sample_annotation) {
            if (sum(rowSums(!is.na(data_matrix)) > 0L) < 2L) {
                stop("PLSDA requires at least two not-NA features.")
            }

            .run_PLSDA_core(
                data_matrix = data_matrix,
                sample_annotation = sample_annotation,
                sample_id_col = sample_id_col,
                batch_col = batch_col,
                effect_col = effect_col,
                ncomp_trt = ncomp_trt,
                ncomp_bat = ncomp_bat,
                keepX_trt = keepX_trt,
                keepX_bat = keepX_bat,
                balance = balance,
                near_zero_var = near_zero_var,
                max.iter = max.iter,
                tol = tol,
                run_splsda = run_splsda,
                ...
            )
        }
    )
}

# ---- internal core (matrix -> matrix), features x samples in/out ----
.run_PLSDA_core <- function(
  data_matrix, # features × samples (numeric)
  sample_annotation, # data.frame
  sample_id_col,
  batch_col,
  effect_col = NULL,
  ncomp_trt = NULL,
  ncomp_bat = NULL,
  keepX_trt = NULL,
  keepX_bat = NULL,
  balance = "auto",
  near_zero_var = TRUE,
  max.iter = 500,
  tol = 1e-06,
  run_splsda = FALSE,
  ...
) {
    stopifnot(is.matrix(data_matrix))
    if (is.null(sample_annotation)) stop("sample_annotation is required.")
    if (run_splsda) {
        message("Running sPLSDA-batch ...")
    }

    # Align SA to matrix columns (samples)
    sample_annotation <- .align_sample_annotation(
        sample_annotation,
        sample_ids    = colnames(data_matrix),
        sample_id_col = sample_id_col
    )

    # Build samples×features matrix for PLSDAbatch
    X <- t(data_matrix) # samples × features
    if (anyNA(X)) stop("PLSDA-batch requires no NAs; impute/filter before calling.")
    storage.mode(X) <- "double"

    # Ensure SA contains required columns
    if (!(batch_col %in% names(sample_annotation))) {
        stop("Batch column '", batch_col, "' is not present in sample_annotation.")
    }

    Y.bat <- factor(sample_annotation[[batch_col]])
    Y.trt <- if (!is.null(effect_col)) factor(sample_annotation[[effect_col]]) else NULL

    # --- choose ncomp_trt if needed (robust to Y.trt NULL) ---
    if (is.null(ncomp_trt)) {
        if (is.null(Y.trt)) {
            message("ncomp_trt not supplied and effect_col is NULL; using default ncomp_trt=2.")
            ncomp_trt <- 2L
        } else {
            message("ncomp_trt not supplied; Running `mixOmics::plsda` to choose ...")
            ncomp_trt <- .select_ncomp_trt(X, Y.trt)
        }
    } else {
        message(
            "Using user-supplied ncomp_trt = ", ncomp_trt,
            ". If unsure, leave NULL to run pre-selection."
        )
    }
    ncomp_trt <- as.integer(ncomp_trt)
    if (ncomp_trt < 1L) stop("ncomp_trt must be a positive integer.")

    if (is.null(Y.trt)) {
        message("effect_col is NULL; running PLSDA-batch with no treatment preservation.")
    }

    # --- sPLSDA tuning of keepX_trt if requested ---
    dots <- list(...)
    if (!is.null(keepX_trt) && length(keepX_trt) != ncomp_trt) {
        stop("keepX_trt must be NULL or a numeric vector of length ncomp_trt.")
    }
    if (isTRUE(run_splsda) && is.null(keepX_trt)) {
        message("keepX_trt not supplied; estimating via mixOmics::tune.splsda ...")
        seed <- if (!is.null(dots[["tune.seed"]])) dots[["tune.seed"]] else 777L
        folds <- if (!is.null(dots[["tune.folds"]])) dots[["tune.folds"]] else 4L
        nrepeat <- if (!is.null(dots[["tune.nrepeat"]])) dots[["tune.nrepeat"]] else 50L
        grid <- dots[["test.keepX"]]
        keepX_trt <- .select_keepX_trt(
            X = X, Y.trt = Y.trt, ncomp_trt = ncomp_trt,
            folds = folds, nrepeat = nrepeat, test.keepX = grid
        )
        # Strip tuning-only args; pass the rest through
        tune_keys <- c("tune.seed", "tune.folds", "tune.nrepeat", "test.keepX")
        dots <- if (length(dots)) dots[setdiff(names(dots), tune_keys)] else list()
    }

    # --- balance logic: accept logical or "auto" ---
    if (is.character(balance) && identical(balance, "auto")) {
        if (is.null(Y.trt)) {
            bal <- TRUE
        } else {
            tab <- table(Y.bat, Y.trt)
            bal <- all(tab > 0L) # balanced (complete) design if no zero cells
        }
    } else if (is.logical(balance) && length(balance) == 1L) {
        bal <- balance
    } else {
        stop("balance must be logical or \"auto\".")
    }

    # --- choose ncomp_bat if needed ---
    if (is.null(ncomp_bat)) {
        message("ncomp_bat not supplied; Running `PLSDA_batch` to choose ...")
        ncomp_bat <- .select_ncomp_bat(X, Y.trt, Y.bat, ncomp_trt)
    } else {
        message(
            "Using user-supplied ncomp_bat = ", ncomp_bat,
            ". If unsure, leave NULL to run pre-selection."
        )
    }
    ncomp_bat <- as.integer(ncomp_bat)
    if (ncomp_bat < 1L) stop("ncomp_bat must be a positive integer.")

    # --- fit PLSDA-batch (pass tol and max.iter) ---
    fit <- do.call(
        PLSDAbatch::PLSDA_batch,
        c(list(
            X = X,
            Y.trt = Y.trt,
            Y.bat = Y.bat,
            ncomp.trt = ncomp_trt,
            ncomp.bat = ncomp_bat,
            keepX.trt = if (is.null(keepX_trt)) rep(ncol(X), ncomp_trt) else keepX_trt,
            keepX.bat = if (is.null(keepX_bat)) rep(ncol(X), ncomp_bat) else keepX_bat,
            max.iter = max.iter,
            tol = tol,
            near.zero.var = near_zero_var,
            balance = bal
        ), dots)
    )

    corrected <- fit$X.nobatch # samples x features
    out <- t(corrected) # back to features x samples
    dimnames(out) <- dimnames(data_matrix)
    storage.mode(out) <- "double"
    out
}

.select_ncomp_bat <- function(x, Y.trt, Y.bat, ncomp_trt) {
    n_samples <- nrow(x)
    n_features <- ncol(x)
    n_levels_bat <- nlevels(Y.bat)
    max_allowed <- min(
        n_features,
        max(1L, n_samples - 1L),
        max(1L, n_levels_bat - 1L)
    )
    candidate <- as.integer(max_allowed)
    attempt <- NULL
    while (candidate >= 1L && is.null(attempt)) {
        attempt <- try(
            PLSDAbatch::PLSDA_batch(
                X = x,
                Y.trt = Y.trt,
                Y.bat = Y.bat,
                ncomp.trt = ncomp_trt,
                ncomp.bat = candidate
            ),
            silent = TRUE
        )
        if (inherits(attempt, "try-error")) {
            candidate <- candidate - 1L
            attempt <- NULL
        }
    }
    if (is.null(attempt)) {
        stop("Unable to estimate ncomp_bat: PLSDAbatch::PLSDA_batch failed for all candidate component counts.")
    }
    if (candidate < max_allowed) {
        message("Reducing maximum ncomp_bat to ", candidate, " after PLSDAbatch::PLSDA_batch convergence issues.")
    }
    ad.trt.tune <- attempt
    ncomp_bat <- candidate

    pe <- ad.trt.tune$explained_variance.bat
    peY <- unname(pe$Y)
    ycum <- cumsum(peY)
    tab_pct <- round(100 * rbind(X = unname(pe$X), Y = peY, Y_cumsum = pmin(1, ycum)), 1)
    colnames(tab_pct) <- paste0("comp", seq_len(ncol(tab_pct)))

    message("It is recommended (by PLSDA-batch) to select ncomp_bat that explains close to 100% of variance in $Y.")
    message(
        "Proportion of variance explained (%) for \n",
        paste(capture.output(print(tab_pct, quote = FALSE, right = TRUE)), collapse = "\n")
    )

    idx <- which(ycum >= 0.99)[1]
    value <- if (length(idx) && !is.na(idx)) idx else ncomp_bat
    message("Suggesting ncomp_bat = ", value, " (of maximum ", ncomp_bat, "). Rerun with user-supplied ncomp_bat if needed.")
    as.integer(value)
}

.select_ncomp_trt <- function(x, Y.trt) {
    if (is.null(Y.trt)) {
        warning("Cannot select ncomp_trt when effect_col is NULL. Using default 2.")
        return(2L)
    }
    .pb_requireNamespace("mixOmics")
    Y.trt <- droplevels(as.factor(Y.trt))
    n_levels <- nlevels(Y.trt)
    if (n_levels < 2L) stop("effect_col must have at least two levels to select ncomp_trt.")
    tab_trt <- table(Y.trt)
    if (any(tab_trt < 3L)) warning("Some treatment levels have fewer than 3 samples; mixOmics::plsda() may fail or be unstable.")

    n_samples <- nrow(x)
    n_features <- ncol(x)
    max_allowed <- min(n_features, max(1L, n_samples - 1L), n_levels - 1L)
    if (max_allowed < 1L) stop("Unable to select ncomp_trt: insufficient samples or treatment levels.")

    rep_limit <- sum(tab_trt > 3L)
    max_ncomp <- if (rep_limit < 1L) max_allowed else min(rep_limit, max_allowed)

    ad.trt.tune <- mixOmics::plsda(X = x, Y = Y.trt, ncomp = max_ncomp)

    pe <- ad.trt.tune$prop_expl_var

    peY <- if (is.list(pe) && !is.null(pe$Y)) unname(pe$Y) else as.numeric(pe)
    ycum <- cumsum(pmax(0, peY))
    peX <- if (is.list(pe) && !is.null(pe$X)) unname(pe$X) else rep(NA_real_, length(peY))
    tab_pct <- round(100 * rbind(X = peX, Y = peY, Y_cumsum = pmin(1, ycum)), 1)
    colnames(tab_pct) <- paste0("comp", seq_len(ncol(tab_pct)))

    message("It is recommended (by PLSDA-batch) to select ncomp_trt that explains close to 100% of variance in Y.")
    message(
        "Proportion of variance explained (%)\n",
        paste(capture.output(print(tab_pct, quote = FALSE, right = TRUE)), collapse = "\n")
    )

    idx <- which(ycum >= 0.99)[1]
    value <- if (length(idx) && !is.na(idx)) idx else max_ncomp
    message("Suggesting ncomp_trt = ", value, " (of maximum ", max_ncomp, "). Rerun with user-supplied ncomp_trt if needed.")
    as.integer(value)
}


# ---- helpers (reused by sPLSDA tuning only) ----------------------------------

.default_keepX_grid <- function(p) {
    unique(sort(pmin(c(1:10, seq(20, 100, 10), seq(150, 500, 50), p), p)))
}

# X: samples×features; accepts optional tuning args via function params
.select_keepX_trt <- function(X, Y.trt, ncomp_trt,
                              folds = 10L, nrepeat = 1L,
                              test.keepX = NULL) {
    .pb_requireNamespace("mixOmics")
    p <- ncol(X)
    if (is.null(test.keepX)) test.keepX <- .default_keepX_grid(p)
    tun <- mixOmics::tune.splsda(
        X = X, Y = Y.trt,
        ncomp = as.integer(ncomp_trt),
        test.keepX = as.integer(test.keepX),
        validation = "Mfold",
        folds = as.integer(folds),
        nrepeat = as.integer(nrepeat),
        progressBar = TRUE
    )
    as.integer(tun$choice.keepX)
}
