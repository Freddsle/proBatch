#' RUV-III-C batch normalization (requires RUVIIIC)
#'
#' @description
#' Adjusts batch effects using the [RUV-III-C](https://github.com/CMRI-ProCan/RUV-III-C)
#' algorithm, which extends RUV-III to handle missing values by correcting each
#' feature independently. Works on long data.frames (`format = "long"`) and wide
#' matrices (`format = "wide"`) by converting to/from matrices.
#'
#' RUV-III-C normalizes each peptide/protein on the subset of
#' samples where that peptide and all controls are observed. Input - log2-transformed.
#'
#' @inheritParams correct_with_ComBat
#' @param format One of `"long"` or `"wide"`.
#' @param replicate_col Column in `sample_annotation` that identifies technical
#'   replicate groups (used to build the design matrix `M`). Must be present and
#'   non-missing for all samples involved in the correction.
#' @param negative_control_features Character vector naming the negative-control
#'   features (columns) used by RUV-III-C. All entries must be row names of the
#'   input data matrix (i.e. feature identifiers).
#' @param k Integer; number of unwanted factors to remove. For each peptide normalized
#'   it must satisfy rows(M) − cols(M) ≥ k on the subset used for that peptide.
#'   If k is too high, some features will be dropped by necessity.
#' @param to_correct Optional character vector restricting the set of features to
#'   correct. Defaults to all features present in the matrix.
#' @param keep_all Columns retained in the long-format output; passed to
#'   [subset_keep_cols()]. Ignored for wide matrices.
#' @param version Either `"CPP"` (default; fast C++ implementation) or `"R"`
#'   (prototype R implementation). Forwarded to [RUVIIIC::RUVIII_C()].
#' @param with_extra,with_alpha,with_w,progress Optional flags forwarded to
#'   [RUVIIIC::RUVIII_C()].
#' @param use_pseudorep Logical; if `TRUE` and **no** technical replicates are provided
#'   (all replicate groups are size 1), construct **PRPS** pseudo-replicate pairs and
#'   use them as `M`. If `FALSE` (default) and no technical replicates exist, stop with
#'   an informative error.
#' @param prps_group_cols Optional character vector of column names in `sample_annotation`
#'   defining *biological* groups within which to build pseudo-replicates (e.g. `"condition"`).
#'   If `NULL`, the function will use `"condition"` when present; otherwise it errors when
#'   PRPS is requested.
#' @param rep_force Logical; if `TRUE`, proceed with RUV-III-C correction even when no
#'   technical replicates are detected (all replicate groups are size 1) and
#'   `use_pseudorep = FALSE`. This is not recommended.
#' @param ... Additional arguments forwarded to [RUVIIIC::RUVIII_C()].
#'
#' @return For `format = "wide"`, a numeric matrix (features × samples) or, if
#'   `with_extra = TRUE`, a list matching the output of
#'   [RUVIIIC::RUVIII_C()] with `newY` transposed back to features × samples.
#'   For `format = "long"`, a long data.frame with corrected values in
#'   `measure_col` and original values preserved in `preBatchCorr_[measure_col]`.
#'   When `with_extra = TRUE`, a list is returned containing `newY` (features ×
#'   samples) plus an extra `corrected_long` element with the long-format output.
#'
#' @references Poulos, R.C., Hains, P.G., Shah, R. et al. Strategies to enable
#'   large-scale proteomics for reproducible research. Nat Commun 11, 3793 (2020).
#'   https://doi.org/10.1038/s41467-020-17641-3
#'
#'   GitHub: https://github.com/CMRI-ProCan/RUV-III-C
#' @seealso [correct_batch_effects()], [RUVIIIC::RUVIII_C()]
#' @export
correct_with_RUVIII_C <- function(
    x,
    sample_annotation,
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    sample_id_col = "FullRunName",
    replicate_col,
    negative_control_features,
    k,
    format = c("long", "wide"),
    keep_all = "default",
    version = c("CPP", "R"),
    to_correct = NULL,
    with_extra = FALSE,
    with_alpha = FALSE,
    with_w = FALSE,
    progress = TRUE,
    use_pseudorep = FALSE,
    prps_group_cols = NULL,
    rep_force = FALSE,
    ...) {

    .pb_requireNamespace("RUVIIIC")
    format <- match.arg(format)
    version <- match.arg(version)

    .check_ruviiic_inputs(
        replicate_col = replicate_col,
        negative_control_features = negative_control_features,
        k = k
    )

    negative_control_features <- unique(as.character(negative_control_features))
    if (identical(format, "wide")) {
        if (!is.matrix(x)) {
            stop("format='wide' requires a numeric matrix.")
        }
        return(.ruviiic_matrix_step(
            data_matrix = x,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            replicate_col = replicate_col,
            negative_control_features = negative_control_features,
            k = k,
            version = version,
            to_correct = to_correct,
            with_extra = with_extra,
            with_alpha = with_alpha,
            with_w = with_w,
            progress = progress,
            use_pseudorep = use_pseudorep,
            prps_group_cols = prps_group_cols,
            rep_force = TRUE,
            ...
        ))
    }

    if (!is.data.frame(x)) {
        stop("format='long' requires a data.frame.")
    }
    df_long <- x
    original_cols <- names(df_long)

    df_long <- check_sample_consistency(
        sample_annotation, sample_id_col, df_long,
        batch_col = NULL,
        order_col = NULL,
        facet_col = NULL,
        merge = FALSE
    )
    data_matrix <- long_to_matrix(
        df_long,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col
    )
    corrected <- .ruviiic_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        replicate_col = replicate_col,
        negative_control_features = negative_control_features,
        k = k,
        version = version,
        to_correct = to_correct,
        with_extra = with_extra,
        with_alpha = with_alpha,
        with_w = with_w,
        progress = progress,
        use_pseudorep = use_pseudorep,
        prps_group_cols = prps_group_cols,
        rep_force = TRUE,
        ...
    )
    if (is.list(corrected)) {
        corrected$corrected_long <- .post_correction_to_long(
            corrected_matrix = corrected$newY,
            df_long = df_long,
            feature_id_col = feature_id_col,
            measure_col = measure_col,
            sample_id_col = sample_id_col,
            original_cols = original_cols,
            keep_all = keep_all
        )
        return(corrected)
    }
    .post_correction_to_long(
        corrected_matrix = corrected,
        df_long = df_long,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col,
        original_cols = original_cols,
        keep_all = keep_all
    )
}

.check_ruviiic_inputs <- function(
    replicate_col, negative_control_features, k) {
    if (missing(replicate_col) || is.null(replicate_col) || !nzchar(replicate_col)) {
        stop("replicate_col must be provided and non-empty.")
    }
    if (missing(negative_control_features) || !length(negative_control_features)) {
        stop("negative_control_features must be a non-empty character vector.")
    }
    if (missing(k) || length(k) != 1 || is.na(k) || !is.numeric(k)) {
        stop("k must be a single positive numeric value.")
    }
    k <- as.integer(k)
    if (k < 1) {
        stop("k must be >= 1 for RUV-III-C correction.")
    }
}

.ruviiic_matrix_step <- function(
    data_matrix,
    sample_annotation,
    sample_id_col = "FullRunName",
    replicate_col,
    negative_control_features,
    k,
    version = c("CPP", "R"),
    to_correct = NULL,
    with_extra = FALSE,
    with_alpha = FALSE,
    with_w = FALSE,
    progress = TRUE,
    use_pseudorep = FALSE,
    prps_group_cols = NULL,
    rep_force = TRUE,
    ...) {

    if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
    }
    storage.mode(data_matrix) <- "double"

    feature_ids <- rownames(data_matrix)
    sample_ids <- colnames(data_matrix)

    if (is.null(feature_ids) || !length(feature_ids)) {
        stop("data_matrix must have row names (feature identifiers).")
    }
    if (is.null(sample_ids) || !length(sample_ids)) {
        stop("data_matrix must have column names (sample identifiers).")
    }

    sample_annotation <- .align_sample_annotation(
        sample_annotation,
        sample_ids = sample_ids,
        sample_id_col = sample_id_col
    )
    if (!(replicate_col %in% names(sample_annotation))) {
        stop(sprintf("replicate_col '%s' is not present in sample_annotation.", replicate_col))
    }
    replicate_ids <- sample_annotation[[replicate_col]]
    if (any(is.na(replicate_ids))) {
        stop("replicate_col contains missing values; provide replicate identifiers for all samples.")
    }
    replicate_ids <- as.character(replicate_ids)

    # ---- PRPS decision: only when ALL replicate groups are singletons ----
    tab <- table(replicate_ids)
    all_singletons <- all(tab == 1L)

    prps_used <- FALSE
    kept <- rep(TRUE, length(sample_ids))
    if (all_singletons) {
        if (!isTRUE(use_pseudorep)) {
            stop("No technical replicates detected (all replicate groups are size 1). ",
                 "Set use_pseudorep = TRUE and supply prps_group_cols (e.g., 'condition') to enable PRPS."
                 "Alternatively, force correction without replicates - using rep_force = TRUE - with caution.",
                 "See: https://www.nature.com/articles/s41587-022-01440-w")
        }
        prps_ids <- .prps_make_ids(
            data_matrix = data_matrix,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            group_cols = prps_group_cols,
            control_features = negative_control_features
        )
        if (!anyDuplicated(prps_ids[!is.na(prps_ids)])) {
            stop("PRPS could not form any replicate pairs; not enough within-group similarity or too few samples.")
        }
        kept <- !is.na(prps_ids)
        if (sum(kept) < 2L) {
            stop("PRPS formed fewer than 2 samples; cannot proceed.")
        }
        if (any(!kept)) {
            warning(sprintf(
                "PRPS: %d sample(s) could not be paired and will be passed through unchanged.",
                sum(!kept)
            ))
        }
        replicate_ids <- prps_ids
        prps_used <- TRUE
    }

    # build design matrix M (only over kept samples if PRPS used)
    replicate_levels <- unique(replicate_ids[kept])
    design_matrix <- matrix(
        0,
        nrow = sum(kept),
        ncol = length(replicate_levels),
        dimnames = list(sample_ids[kept], replicate_levels)
    )
    design_matrix[
        cbind(seq_len(sum(kept)), match(replicate_ids[kept], replicate_levels))
    ] <- 1
    storage.mode(design_matrix) <- "double"

    # to_correct + controls checks (unchanged)
    if (is.null(to_correct)) {
        to_correct <- feature_ids
    } else {
        to_correct <- unique(as.character(to_correct))
        missing_features <- setdiff(to_correct, feature_ids)
        if (length(missing_features)) {
            stop(
                "to_correct contains features absent from data_matrix: ",
                paste(missing_features, collapse = ", ")
            )
        }
    }

    missing_controls <- setdiff(negative_control_features, feature_ids)
    if (length(missing_controls)) {
        stop(
            "negative_control_features missing from data_matrix: ",
            paste(missing_controls, collapse = ", ")
        )
    }

    # core run (possibly subset on kept if PRPS used)
    result <- .run_RUVIIIC_core(
        k = k,
        Y = t(if (prps_used) data_matrix[, kept, drop = FALSE] else data_matrix),
        M = design_matrix,
        to_correct = to_correct,
        controls = negative_control_features,
        version = match.arg(version),
        with_extra = with_extra,
        with_alpha = with_alpha,
        with_w = with_w,
        progress = progress,
        ...
    )

    # transpose back and reassemble if PRPS dropped samples
    if (is.list(result) && !is.null(result$newY)) {
        corrected <- as.matrix(result$newY)  # samples x features
        storage.mode(corrected) <- "double"
        corrected <- t(corrected)            # features x samples (kept only)
        if (prps_used && any(!kept)) {
            full <- data_matrix
            full[, kept] <- corrected
            corrected <- full
        }
        result$newY <- corrected
        return(result)
    }

    corrected <- as.matrix(result)  # samples x features
    storage.mode(corrected) <- "double"
    corrected <- t(corrected)       # features x samples (kept only)
    if (prps_used && any(!kept)) {
        full <- data_matrix
        full[, kept] <- corrected
        corrected <- full
    }
    corrected
}

.run_RUVIIIC_core <- function(
    k,
    Y,
    M,
    to_correct,
    controls,
    version = c("CPP", "R"),
    with_extra = FALSE,
    with_alpha = FALSE,
    with_w = FALSE,
    progress = TRUE,
    ...) {
    version <- match.arg(version)
    RUVIIIC::RUVIII_C(
        k = k,
        Y = Y,
        M = M,
        toCorrect = to_correct,
        controls = controls,
        withExtra = with_extra,
        withAlpha = with_alpha,
        withW = with_w,
        version = version,
        progress = progress,
        ...
    )
}

# --- internal helper to build PRPS replicate ids (greedy DN pairing) ---
# TODO: NOT TESTED YET! Use with caution.
.prps_make_ids <- function(
    data_matrix,                      # features x samples (double; with NAs)
    sample_annotation,                # data.frame aligned to samples
    sample_id_col = "FullRunName",
    group_cols = NULL,                # biological grouping (e.g., "condition")
    control_features = character(0)   # features to compute similarity on; fallback -> all
) {
    message("PRPS: constructing pseudo-replicates based on within-group similarity",
            "This functionality is experimental; use with caution.")
    sample_ids <- colnames(data_matrix)
    if (is.null(group_cols) || !length(group_cols)) {
        if ("condition" %in% names(sample_annotation)) {
            group_cols <- "condition"
        } else {
            stop("PRPS requested but no prps_group_cols provided and no 'condition' column found in sample_annotation.")
        }
    }
    if (!all(group_cols %in% names(sample_annotation))) {
        stop("All prps_group_cols must be present in sample_annotation.")
    }
    # group key per sample
    key <- do.call(interaction, c(sample_annotation[group_cols], drop = TRUE))
    key <- as.character(key)
    prps_id <- rep(NA_character_, length(sample_ids))
    names(prps_id) <- sample_ids

    # restrict features for similarity
    features <- rownames(data_matrix)
    use_feats <- intersect(control_features, features)
    if (!length(use_feats)) {
        use_feats <- features
        warning("PRPS: none of the negative_control_features are present among rows; using all features to build pseudo-replicates.")
    }

    # greedy, within-group disjoint pairing by maximum correlation (pairwise complete)
    for (g in unique(key)) {
        idx <- which(key == g)
        if (length(idx) < 2) next
        X <- t(data_matrix[use_feats, idx, drop = FALSE])  # samples x features
        suppressWarnings({
            C <- stats::cor(t(X), use = "pairwise.complete.obs")  # pairwise, robust to NAs
        })
        if (anyNA(C)) C[is.na(C)] <- 0
        diag(C) <- -Inf
        used <- rep(FALSE, length(idx))
        pair_no <- 0L

        while (sum(!used) >= 2L) {
            avail <- which(!used)
            Csub <- C[avail, avail, drop = FALSE]
            if (!length(Csub)) break
            # best pair among available
            kmax <- which(Csub == max(Csub), arr.ind = TRUE)[1, , drop = TRUE]
            i <- avail[kmax[1]]; j <- avail[kmax[2]]
            if (i == j) break
            pair_no <- pair_no + 1L
            rid <- paste0("PRPS.", g, ".p", pair_no)
            prps_id[idx[c(i, j)]] <- rid
            used[c(i, j)] <- TRUE
        }
        # leftover (if any) stays NA -> will be passed through unchanged
    }
    prps_id
}
