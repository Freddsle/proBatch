#' Calculate clustering agreement metrics (ARI, MCC, silhouette, Jaccard)
#'
#' Compute clustering agreement metrics using known class labels and either
#' provided cluster labels or k-means clustering when no cluster column is
#' supplied. The function returns one row per requested `known_col` with ARI,
#' MCC, mean silhouette width, and two Jaccard summaries. Sample-level
#' assignments and per-cluster silhouette summaries are attached as attributes.
#'
#' @inheritParams proBatch
#' @param data_matrix Input object: matrix-like data or a `ProBatchFeatures` instance.
#' @param sample_annotation Data frame with sample-level metadata. Must include
#'   `known_col` (and optionally `cluster_col`).
#' @param known_col Column name(s) in `sample_annotation` with known class
#'   labels.
#' @param cluster_col Optional column in `sample_annotation` with provided
#'   cluster labels. When `NULL`, k-means is used to derive clusters.
#' @param k Integer; number of clusters for k-means. If `NULL`, the number of
#'   distinct classes in `known_col` is used.
#' @param raw_matrix Optional uncorrected reference matrix used to compute
#'   neighbor-preservation Jaccard (`jaccard_knn`). When `NULL`, `jaccard_knn`
#'   is returned as `NA`.
#' @param k_neighbors Integer; number of neighbors used by the kNN-based
#'   Jaccard metric.
#' @param fill_the_missing Missing-value policy applied before k-means and
#'   silhouette calculations (see [handle_missing_values()]).
#' @param dist_method Distance method passed to [stats::dist()] for silhouette
#'   computation.
#' @param pbf_name Assay name(s) used when `data_matrix` is a `ProBatchFeatures`.
#' @param raw_pbf_name Assay name(s) providing the uncorrected reference matrix
#'   when `data_matrix` is a `ProBatchFeatures` object. Use length 1 to apply
#'   one reference assay to all evaluated assays, or provide one name per assay.
#' @param ... Additional arguments forwarded to [stats::kmeans()] when k-means
#'   clustering is used.
#'
#' @return A data frame with one row per evaluated `known_col`, including
#'   classification summary metrics and metadata. In addition to `ARI`, `MCC`,
#'   and `silhouette`, the output includes `jaccard_knn` and `jaccard_label`.
#'   The result carries two attributes: `pb_assignments` and
#'   `pb_silhouette_clusters`. Each is a single object for one `known_col`, or a
#'   named list when multiple `known_col` values are evaluated.
#' @export
#'
#' @examples
#' data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#' metrics <- calculate_classification_metrics(
#'     example_proteome_matrix,
#'     sample_annotation = example_sample_annotation,
#'     raw_matrix = example_proteome_matrix,
#'     known_col = "MS_batch",
#'     k = 2
#' )
calculate_classification_metrics <- function(data_matrix, ...) {
    UseMethod("calculate_classification_metrics")
}

#' @rdname calculate_classification_metrics
#' @method calculate_classification_metrics default
#' @export
calculate_classification_metrics.default <- function(data_matrix,
                                                     sample_annotation,
                                                     known_col,
                                                     cluster_col = NULL,
                                                     sample_id_col = "FullRunName",
                                                     k = NULL,
                                                     fill_the_missing = -1,
                                                     dist_method = "euclidean",
                                                     raw_matrix = NULL,
                                                     k_neighbors = 15,
                                                     ...) {
    if (missing(known_col) || is.null(known_col)) {
        stop("`known_col` must be provided.")
    }
    known_cols <- .pb_normalize_known_cols(known_col)

    if (is(data_matrix, "SummarizedExperiment")) {
        data_matrix <- assay(data_matrix)
    }

    alignment <- .pb_align_matrix_and_annotation(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        check_args = list(
            batch_col = NULL,
            order_col = NULL,
            facet_col = NULL
        )
    )
    data_matrix <- alignment$data_matrix
    sample_annotation <- alignment$sample_annotation
    sample_ids <- alignment$sample_ids

    raw_matrix_aligned <- NULL
    if (!is.null(raw_matrix)) {
        if (is(raw_matrix, "SummarizedExperiment")) {
            raw_matrix <- assay(raw_matrix)
        }
        raw_alignment <- .pb_align_matrix_and_annotation(
            data_matrix = raw_matrix,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            check_args = list(
                batch_col = NULL,
                order_col = NULL,
                facet_col = NULL
            )
        )
        raw_matrix_aligned <- raw_alignment$data_matrix
    }

    if (is.null(sample_annotation)) {
        stop("`sample_annotation` must be provided.")
    }

    if (length(known_cols) > 1L) {
        missing_cols <- setdiff(known_cols, names(sample_annotation))
        if (length(missing_cols)) {
            stop(
                "Known label column(s) not found in sample_annotation: ",
                paste(missing_cols, collapse = ", "),
                "."
            )
        }

        metrics_list <- vector("list", length(known_cols))
        assignments_list <- vector("list", length(known_cols))
        silhouette_clusters_list <- vector("list", length(known_cols))
        names(metrics_list) <- known_cols
        names(assignments_list) <- known_cols
        names(silhouette_clusters_list) <- known_cols

        for (i in seq_along(known_cols)) {
            known_col_i <- known_cols[[i]]
            res <- .pb_calculate_classification_metrics_single(
                data_matrix = data_matrix,
                sample_annotation = sample_annotation,
                sample_ids = sample_ids,
                known_col = known_col_i,
                cluster_col = cluster_col,
                k = k,
                raw_matrix = raw_matrix_aligned,
                k_neighbors = k_neighbors,
                fill_the_missing = fill_the_missing,
                dist_method = dist_method,
                ...
            )
            metrics_list[[i]] <- res
            assignments_list[[i]] <- attr(res, "pb_assignments")
            silhouette_clusters_list[[i]] <- attr(res, "pb_silhouette_clusters")
        }

        metrics <- do.call(rbind, metrics_list)
        rownames(metrics) <- NULL
        attr(metrics, "pb_assignments") <- assignments_list
        attr(metrics, "pb_silhouette_clusters") <- silhouette_clusters_list
        return(metrics)
    }

    .pb_calculate_classification_metrics_single(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_ids = sample_ids,
        known_col = known_cols[[1]],
        cluster_col = cluster_col,
        k = k,
        raw_matrix = raw_matrix_aligned,
        k_neighbors = k_neighbors,
        fill_the_missing = fill_the_missing,
        dist_method = dist_method,
        ...
    )
}

#' @rdname calculate_classification_metrics
#' @method calculate_classification_metrics ProBatchFeatures
#' @export
calculate_classification_metrics.ProBatchFeatures <- function(data_matrix,
                                                              pbf_name = NULL,
                                                              sample_annotation = NULL,
                                                              known_col,
                                                              cluster_col = NULL,
                                                              sample_id_col = "FullRunName",
                                                              k = NULL,
                                                              fill_the_missing = -1,
                                                              dist_method = "euclidean",
                                                              raw_pbf_name = NULL,
                                                              k_neighbors = 15,
                                                              ...) {
    object <- data_matrix
    prep <- .pb_prepare_multi_assay(
        object = object,
        pbf_name = pbf_name,
        dots = list(...),
        plot_title = NULL,
        default_title_fun = function(x) x
    )
    assays <- prep$assays
    split_arg <- prep$split_arg

    default_sample_annotation <- .pb_default_sample_annotation(
        object = object,
        sample_id_col = sample_id_col
    )
    sample_ann_list <- split_arg(sample_annotation)

    metrics_list <- vector("list", length(assays))
    assignments_list <- vector("list", length(assays))
    silhouette_clusters_list <- vector("list", length(assays))
    names(metrics_list) <- assays
    names(assignments_list) <- assays
    names(silhouette_clusters_list) <- assays

    raw_assay_map <- rep(NA_character_, length(assays))
    names(raw_assay_map) <- assays
    if (!is.null(raw_pbf_name)) {
        raw_assays <- as.character(raw_pbf_name)
        if (length(raw_assays) == 1L) {
            raw_assay_map[] <- raw_assays
        } else if (length(raw_assays) == length(assays)) {
            raw_assay_map[] <- raw_assays
        } else {
            stop("`raw_pbf_name` must be length 1 or match the number of selected assays.")
        }
    }

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        dm <- pb_assay_matrix(object, assay_nm)
        raw_dm <- NULL
        raw_assay_nm <- raw_assay_map[[i]]
        if (!is.na(raw_assay_nm) && nzchar(raw_assay_nm)) {
            raw_dm <- pb_assay_matrix(object, raw_assay_nm)
        }
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- default_sample_annotation
        }

        res <- calculate_classification_metrics.default(
            data_matrix = dm,
            sample_annotation = sample_ann,
            known_col = known_col,
            cluster_col = cluster_col,
            sample_id_col = sample_id_col,
            k = k,
            raw_matrix = raw_dm,
            k_neighbors = k_neighbors,
            fill_the_missing = fill_the_missing,
            dist_method = dist_method,
            ...
        )
        res$assay <- assay_nm
        metrics_list[[i]] <- res
        assignments_list[[i]] <- attr(res, "pb_assignments")
        silhouette_clusters_list[[i]] <- attr(res, "pb_silhouette_clusters")
    }

    metrics <- do.call(rbind, metrics_list)
    rownames(metrics) <- NULL
    attr(metrics, "pb_assignments") <- assignments_list
    attr(metrics, "pb_silhouette_clusters") <- silhouette_clusters_list

    metrics
}

.pb_normalize_known_cols <- function(known_col) {
    vals <- unique(as.character(known_col))
    vals <- vals[!is.na(vals)]
    vals <- trimws(vals)
    vals <- vals[nzchar(vals)]
    if (!length(vals)) {
        stop("`known_col` must be provided.")
    }
    vals
}

.pb_calculate_classification_metrics_single <- function(data_matrix,
                                                        sample_annotation,
                                                        sample_ids,
                                                        known_col,
                                                        cluster_col = NULL,
                                                        k = NULL,
                                                        raw_matrix = NULL,
                                                        k_neighbors = 15,
                                                        fill_the_missing = -1,
                                                        dist_method = "euclidean",
                                                        ...) {
    if (!known_col %in% names(sample_annotation)) {
        stop("Known label column '", known_col, "' not found in sample_annotation.")
    }

    known <- sample_annotation[[known_col]]
    keep <- !is.na(known)
    if (!all(keep)) {
        warning(
            "Removing samples with missing known labels: ",
            paste(sample_ids[!keep], collapse = ", ")
        )
        data_matrix <- data_matrix[, keep, drop = FALSE]
        sample_annotation <- sample_annotation[keep, , drop = FALSE]
        sample_ids <- sample_ids[keep]
        known <- known[keep]
    }

    if (!ncol(data_matrix)) {
        stop("No samples remain after filtering missing known labels.")
    }

    known <- .pb_as_factor(known)

    if (is.null(cluster_col)) {
        if (!is.null(k) && length(k) != 1L) {
            stop("`k` must be a single integer or NULL.")
        }
        if (is.null(k)) {
            k <- nlevels(known)
        }
        k <- as.integer(k)
        if (is.na(k) || k < 2L) {
            stop("`k` must be >= 2.")
        }
        if (ncol(data_matrix) < k) {
            stop("`k` cannot be larger than the number of samples.")
        }

        data_matrix_k <- .pb_handle_missing_wrapper(
            data_matrix = data_matrix,
            warning_message = "k-means cannot operate with missing values in the matrix",
            fill_the_missing = fill_the_missing,
            drop_on_false = TRUE
        )

        if (!nrow(data_matrix_k) || !ncol(data_matrix_k)) {
            stop("No data remaining after handling missing values.")
        }
        if (anyNA(data_matrix_k)) {
            stop("Missing values remain after handling; set `fill_the_missing` to a numeric value or remove NAs.")
        }

        km_args <- list(...)
        if ("centers" %in% names(km_args)) {
            stop("Do not pass `centers` in `...`; use the `k` argument instead.")
        }

        km_call <- c(list(
            x = t(as.matrix(data_matrix_k)),
            centers = k
        ), km_args)
        km_res <- do.call(stats::kmeans, km_call)

        predicted <- factor(km_res$cluster)
        cluster_source <- "kmeans"
        k_used <- k
    } else {
        if (!cluster_col %in% names(sample_annotation)) {
            stop("Cluster label column '", cluster_col, "' not found in sample_annotation.")
        }
        if (!is.null(k)) {
            warning("`k` is ignored because `cluster_col` was supplied.")
        }
        predicted <- sample_annotation[[cluster_col]]
        keep <- !is.na(predicted)
        if (!all(keep)) {
            warning(
                "Removing samples with missing cluster labels: ",
                paste(sample_ids[!keep], collapse = ", ")
            )
            data_matrix <- data_matrix[, keep, drop = FALSE]
            sample_annotation <- sample_annotation[keep, , drop = FALSE]
            sample_ids <- sample_ids[keep]
            known <- known[keep]
            predicted <- predicted[keep]
        }

        if (!ncol(data_matrix)) {
            stop("No samples remain after filtering missing cluster labels.")
        }
        predicted <- .pb_as_factor(predicted)
        cluster_source <- "provided"
        k_used <- nlevels(predicted)
    }

    if (length(known) != length(predicted)) {
        stop("Known and cluster labels must have the same length.")
    }

    tab <- table(known, predicted)
    ari <- .pb_adjusted_rand_index(tab)
    mcc <- .pb_multiclass_mcc(tab)
    silhouette_mean <- .pb_silhouette_mean(
        data_matrix = data_matrix,
        clusters = predicted,
        dist_method = dist_method,
        fill_the_missing = fill_the_missing
    )
    silhouette_clusters <- .pb_silhouette_per_cluster(
        data_matrix = data_matrix,
        clusters = predicted,
        dist_method = dist_method,
        fill_the_missing = fill_the_missing
    )
    jaccard_knn <- .pb_knn_jaccard(
        raw_matrix = raw_matrix,
        corrected_matrix = data_matrix,
        k = k_neighbors,
        dist_method = dist_method
    )
    jaccard_label <- .pb_label_jaccard(
        known = known,
        predicted = predicted
    )

    metrics <- data.frame(
        known_col = known_col,
        ARI = ari,
        MCC = mcc,
        silhouette = silhouette_mean,
        jaccard_knn = jaccard_knn,
        jaccard_label = jaccard_label,
        n_samples = length(known),
        n_classes = nlevels(known),
        n_clusters = nlevels(predicted),
        cluster_source = cluster_source,
        k = k_used,
        stringsAsFactors = FALSE
    )

    assignments <- data.frame(
        sample_id = sample_ids,
        known = as.character(known),
        predicted = as.character(predicted),
        stringsAsFactors = FALSE
    )
    attr(metrics, "pb_assignments") <- assignments
    attr(metrics, "pb_silhouette_clusters") <- silhouette_clusters

    metrics
}

.pb_as_factor <- function(x) {
    if (is.factor(x)) {
        return(x)
    }
    if (is.character(x)) {
        return(factor(x))
    }
    factor(x, exclude = NULL)
}

.pb_adjusted_rand_index <- function(confusion) {
    if (length(confusion) == 0L) {
        return(NA_real_)
    }
    n <- sum(confusion)
    if (n < 2L) {
        return(NA_real_)
    }
    n_choose_2 <- function(x) x * (x - 1) / 2
    sum_ij <- sum(n_choose_2(confusion))
    sum_i <- sum(n_choose_2(rowSums(confusion)))
    sum_j <- sum(n_choose_2(colSums(confusion)))
    total <- n_choose_2(n)
    if (total == 0) {
        return(NA_real_)
    }
    expected <- (sum_i * sum_j) / total
    max_index <- 0.5 * (sum_i + sum_j)
    denom <- max_index - expected
    if (denom == 0) {
        return(NA_real_)
    }
    (sum_ij - expected) / denom
}

.pb_multiclass_mcc <- function(confusion) {
    if (length(confusion) == 0L) {
        return(NA_real_)
    }

    row_nms <- rownames(confusion)
    col_nms <- colnames(confusion)
    if (!is.null(row_nms) &&
        !is.null(col_nms) &&
        length(row_nms) == length(col_nms) &&
        setequal(row_nms, col_nms) &&
        !identical(row_nms, col_nms)) {
        confusion <- confusion[, row_nms, drop = FALSE]
    }

    t_k <- rowSums(confusion)
    p_k <- colSums(confusion)
    if (length(t_k) < length(p_k)) {
        t_k <- c(t_k, rep(0, length(p_k) - length(t_k)))
    } else if (length(p_k) < length(t_k)) {
        p_k <- c(p_k, rep(0, length(t_k) - length(p_k)))
    }

    c <- sum(diag(confusion))
    s <- sum(confusion)
    numerator <- (c * s) - sum(p_k * t_k)
    denom <- sqrt((s^2 - sum(p_k^2)) * (s^2 - sum(t_k^2)))
    if (denom == 0) {
        return(NA_real_)
    }
    numerator / denom
}

.pb_knn_jaccard <- function(raw_matrix,
                            corrected_matrix,
                            k = 15,
                            dist_method = "euclidean") {
    if (is.null(raw_matrix)) {
        return(NA_real_)
    }
    if (is.null(colnames(raw_matrix)) || is.null(colnames(corrected_matrix))) {
        warning("kNN Jaccard requires sample IDs in both matrices; returning NA.")
        return(NA_real_)
    }

    if (length(k) != 1L || is.na(k)) {
        stop("`k_neighbors` must be a single integer.")
    }
    k <- as.integer(k)
    if (k < 1L) {
        stop("`k_neighbors` must be >= 1.")
    }

    sample_ids <- intersect(colnames(corrected_matrix), colnames(raw_matrix))
    if (length(sample_ids) < 2L) {
        warning("kNN Jaccard requires at least two overlapping samples; returning NA.")
        return(NA_real_)
    }
    corrected_matrix <- corrected_matrix[, sample_ids, drop = FALSE]
    raw_matrix <- raw_matrix[, sample_ids, drop = FALSE]

    corrected_matrix <- .pb_handle_missing_wrapper(
        data_matrix = corrected_matrix,
        warning_message = "kNN Jaccard cannot operate with missing values in the corrected matrix",
        fill_the_missing = FALSE,
        drop_on_false = TRUE
    )
    raw_matrix <- .pb_handle_missing_wrapper(
        data_matrix = raw_matrix,
        warning_message = "kNN Jaccard cannot operate with missing values in the raw matrix",
        fill_the_missing = FALSE,
        drop_on_false = TRUE
    )

    if (!nrow(corrected_matrix) || !ncol(corrected_matrix) || !nrow(raw_matrix) || !ncol(raw_matrix)) {
        warning("No data remaining for kNN Jaccard after missing-value handling; returning NA.")
        return(NA_real_)
    }
    if (anyNA(corrected_matrix) || anyNA(raw_matrix)) {
        warning("Missing values remain after handling; kNN Jaccard not computed.")
        return(NA_real_)
    }

    n_samples <- length(sample_ids)
    k_eff <- min(k, n_samples - 1L)
    if (k_eff < 1L) {
        return(NA_real_)
    }

    corrected_dist <- stats::dist(t(as.matrix(corrected_matrix)), method = dist_method)
    raw_dist <- stats::dist(t(as.matrix(raw_matrix)), method = dist_method)
    corrected_neighbors <- .pb_knn_indices(as.matrix(corrected_dist), k = k_eff)
    raw_neighbors <- .pb_knn_indices(as.matrix(raw_dist), k = k_eff)

    scores <- vapply(seq_len(n_samples), function(i) {
        neigh_raw <- raw_neighbors[[i]]
        neigh_corr <- corrected_neighbors[[i]]
        union_size <- length(union(neigh_raw, neigh_corr))
        if (!union_size) {
            return(NA_real_)
        }
        length(intersect(neigh_raw, neigh_corr)) / union_size
    }, numeric(1))

    mean(scores, na.rm = TRUE)
}

.pb_knn_indices <- function(dist_matrix, k) {
    n <- nrow(dist_matrix)
    lapply(seq_len(n), function(i) {
        idx <- order(dist_matrix[i, ], na.last = NA)
        idx <- idx[idx != i]
        if (!length(idx)) {
            return(integer(0))
        }
        idx[seq_len(min(k, length(idx)))]
    })
}

.pb_label_jaccard <- function(known, predicted) {
    known <- .pb_as_factor(known)
    predicted <- .pb_as_factor(predicted)
    tab <- table(known, predicted)
    if (length(tab) == 0L) {
        return(NA_real_)
    }

    row_totals <- rowSums(tab)
    col_totals <- colSums(tab)
    if (!length(row_totals) || !length(col_totals)) {
        return(NA_real_)
    }

    jaccard_matrix <- matrix(0,
        nrow = nrow(tab),
        ncol = ncol(tab),
        dimnames = dimnames(tab)
    )
    for (i in seq_len(nrow(tab))) {
        for (j in seq_len(ncol(tab))) {
            intersection <- tab[i, j]
            union_size <- row_totals[[i]] + col_totals[[j]] - intersection
            jaccard_matrix[i, j] <- if (union_size > 0) intersection / union_size else 0
        }
    }

    candidates <- expand.grid(
        row = seq_len(nrow(jaccard_matrix)),
        col = seq_len(ncol(jaccard_matrix)),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    candidates$score <- jaccard_matrix[cbind(candidates$row, candidates$col)]
    candidates <- candidates[order(-candidates$score, candidates$row, candidates$col), , drop = FALSE]

    used_rows <- rep(FALSE, nrow(jaccard_matrix))
    used_cols <- rep(FALSE, ncol(jaccard_matrix))
    selected_scores <- numeric(0)

    for (idx in seq_len(nrow(candidates))) {
        row_i <- candidates$row[[idx]]
        col_i <- candidates$col[[idx]]
        if (used_rows[[row_i]] || used_cols[[col_i]]) {
            next
        }
        used_rows[[row_i]] <- TRUE
        used_cols[[col_i]] <- TRUE
        selected_scores <- c(selected_scores, candidates$score[[idx]])
    }

    if (!length(selected_scores)) {
        return(NA_real_)
    }
    mean(selected_scores)
}

.pb_silhouette_mean <- function(data_matrix,
                                clusters,
                                dist_method = "euclidean",
                                fill_the_missing = -1) {
    clusters <- .pb_as_factor(clusters)
    if (nlevels(clusters) < 2L) {
        warning("Silhouette requires at least two clusters; returning NA.")
        return(NA_real_)
    }
    if (ncol(data_matrix) < 2L) {
        warning("Silhouette requires at least two samples; returning NA.")
        return(NA_real_)
    }

    data_matrix <- .pb_handle_missing_wrapper(
        data_matrix = data_matrix,
        warning_message = "Silhouette cannot operate with missing values in the matrix",
        fill_the_missing = fill_the_missing,
        drop_on_false = TRUE
    )
    if (!nrow(data_matrix) || !ncol(data_matrix)) {
        warning("No data remaining after handling missing values; silhouette not computed.")
        return(NA_real_)
    }
    if (anyNA(data_matrix)) {
        warning("Missing values remain after handling; silhouette not computed.")
        return(NA_real_)
    }

    dist_matrix <- stats::dist(t(as.matrix(data_matrix)), method = dist_method)
    sil <- silhouette(as.integer(clusters), dist_matrix)
    if (!nrow(sil)) {
        return(NA_real_)
    }
    mean(sil[, "sil_width"], na.rm = TRUE)
}

.pb_silhouette_per_cluster <- function(data_matrix,
                                       clusters,
                                       dist_method = "euclidean",
                                       fill_the_missing = -1) {
    clusters <- .pb_as_factor(clusters)
    cluster_levels <- levels(clusters)
    empty <- setNames(rep(NA_real_, length(cluster_levels)), cluster_levels)

    if (nlevels(clusters) < 2L) {
        warning("Silhouette requires at least two clusters; returning NA per cluster.")
        return(empty)
    }
    if (ncol(data_matrix) < 2L) {
        warning("Silhouette requires at least two samples; returning NA per cluster.")
        return(empty)
    }

    data_matrix <- .pb_handle_missing_wrapper(
        data_matrix = data_matrix,
        warning_message = "Silhouette cannot operate with missing values in the matrix",
        fill_the_missing = fill_the_missing,
        drop_on_false = TRUE
    )
    if (!nrow(data_matrix) || !ncol(data_matrix)) {
        warning("No data remaining after handling missing values; per-cluster silhouette not computed.")
        return(empty)
    }
    if (anyNA(data_matrix)) {
        warning("Missing values remain after handling; per-cluster silhouette not computed.")
        return(empty)
    }

    dist_matrix <- stats::dist(t(as.matrix(data_matrix)), method = dist_method)
    sil <- silhouette(as.integer(clusters), dist_matrix)
    if (!nrow(sil)) {
        return(empty)
    }

    sil_clusters <- as.integer(sil[, "cluster"])
    sil_width <- sil[, "sil_width"]
    mean_by_cluster <- tapply(sil_width, sil_clusters, mean, na.rm = TRUE)

    out <- empty
    cluster_idx <- suppressWarnings(as.integer(names(mean_by_cluster)))
    valid <- !is.na(cluster_idx) &
        cluster_idx >= 1L &
        cluster_idx <= length(cluster_levels)
    if (any(valid)) {
        out[cluster_levels[cluster_idx[valid]]] <- as.numeric(mean_by_cluster[valid])
    }

    out
}
