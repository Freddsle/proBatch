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

    .post_correction_to_long(
        corrected_matrix, df_long,
        feature_id_col, measure_col, sample_id_col,
        original_cols, keep_all
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
        ...
    )

    .post_correction_to_long(
        corrected_matrix, df_long,
        feature_id_col, measure_col, sample_id_col,
        original_cols, keep_all
    )
}

.post_correction_to_long <- function(corrected_matrix, df_long,
                                     feature_id_col, measure_col, sample_id_col,
                                     original_cols, keep_all) {
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

.splsda_matrix_step <- function(run_splsda = TRUE, ...) {
    .plsda_matrix_step(run_splsda = run_splsda, ...)
}


.plsda_matrix_step <- function(
    data_matrix, sample_annotation = NULL,
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
    keep_all = "default",
    run_splsda = FALSE,
    ...) {
    # Coerce to numeric matrix
    if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
    }
    if (!is.numeric(data_matrix)) {
        stop("Input must be coercible to a numeric matrix for PLSDA correction.")
    }
    storage.mode(data_matrix) <- "double"
    if (sum(rowSums(!is.na(data_matrix)) > 0L) < 2L) {
        stop("PLSDA requires at least two not-NA features.")
    }

    .run_PLSDA_core(
        data_matrix, # features × samples (numeric)
        sample_annotation, # data.frame
        sample_id_col,
        batch_col,
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
    ...) {
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
    value <- if (length(idx)) idx else ncomp_bat
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
