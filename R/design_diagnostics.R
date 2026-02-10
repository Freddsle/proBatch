#' Validate batch/design annotations
#'
#' @inheritParams proBatch
#' @param condition_col optional column in \code{sample_annotation} indicating
#'   the primary biological condition.
#' @param covariates optional character vector of additional covariate columns.
#' @param strict logical; when \code{TRUE}, stop on validation errors.
#'
#' @return A list with \code{errors}, \code{warnings}, and \code{summary}. The
#'   object has class \code{"pb_design_check"}.
#' @export
#'
#' @examples
#' data("example_ecoli_data", package = "proBatch")
#' meta <- example_ecoli_data$all_metadata
#' validate_batch_design(
#'     meta,
#'     batch_col = "Lab",
#'     condition_col = "Condition",
#'     sample_id_col = "Run",
#'     strict = FALSE
#' )
validate_batch_design <- function(sample_annotation,
                                  batch_col,
                                  condition_col = NULL,
                                  covariates = NULL,
                                  sample_id_col = NULL,
                                  strict = TRUE) {
    if (is.null(sample_annotation) || !is.data.frame(sample_annotation)) {
        stop("sample_annotation must be a data.frame.")
    }
    if (is.null(batch_col) || !nzchar(batch_col)) {
        stop("`batch_col` must be provided.")
    }

    sample_annotation <- as.data.frame(sample_annotation)

    errors <- character()
    warnings <- character()

    design_cols <- .pb_design_collect_cols(
        batch_col = batch_col,
        condition_col = condition_col,
        covariates = covariates
    )

    missing_cols <- setdiff(design_cols, names(sample_annotation))
    if (length(missing_cols)) {
        errors <- c(errors, sprintf(
            "Missing columns in sample_annotation: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    id_source <- NULL
    if (!is.null(sample_id_col)) {
        if (!(sample_id_col %in% names(sample_annotation))) {
            errors <- c(errors, sprintf(
                "Sample ID column '%s' not found in sample_annotation.",
                sample_id_col
            ))
        } else {
            id_source <- sample_id_col
            ids <- sample_annotation[[sample_id_col]]
            if (anyNA(ids) || any(!nzchar(as.character(ids)))) {
                errors <- c(errors, "Sample ID column contains NA/empty values.")
            }
            if (any(duplicated(ids))) {
                errors <- c(errors, "Duplicated sample IDs detected in sample_id_col.")
            }
        }
    } else if (!is.null(rownames(sample_annotation))) {
        if (.pb_has_default_rownames(sample_annotation)) {
            warnings <- c(
                warnings,
                "Sample annotation rownames look like defaults (1..n); provide sample_id_col for ID checks."
            )
        } else {
            ids <- rownames(sample_annotation)
            if (anyNA(ids) || any(!nzchar(as.character(ids)))) {
                errors <- c(errors, "Sample annotation rownames contain NA/empty values.")
            }
            if (any(duplicated(ids))) {
                errors <- c(errors, "Duplicated sample IDs detected in sample_annotation rownames.")
            }
            id_source <- "rownames"
        }
    } else {
        warnings <- c(warnings, "No sample_id_col or rownames provided; skipping sample ID duplication checks.")
    }

    present_cols <- intersect(design_cols, names(sample_annotation))
    for (col in present_cols) {
        values <- sample_annotation[[col]]
        if (anyNA(values)) {
            errors <- c(errors, sprintf("Column '%s' contains missing values.", col))
        }

        if (is.character(values)) {
            warnings <- c(warnings, sprintf(
                "Column '%s' is character; consider converting to factor for consistent coding.",
                col
            ))
        }

        if (is.numeric(values)) {
            n_unique <- length(unique(values[!is.na(values)]))
            if (n_unique <= 10 || n_unique < 0.1 * nrow(sample_annotation)) {
                warnings <- c(warnings, sprintf(
                    "%s column has very few values, but is numeric-like; should it be treated as factor?",
                    col
                ))
            }
        }

        if (is.factor(values)) {
            level_counts <- table(values, useNA = "no")
            empty_levels <- names(level_counts)[level_counts == 0]
            if (length(empty_levels)) {
                warnings <- c(warnings, sprintf(
                    "Column '%s' has unused/empty factor levels: %s",
                    col,
                    paste(empty_levels, collapse = ", ")
                ))
            }
            trimmed <- trimws(levels(values))
            if (any(trimmed != levels(values))) {
                warnings <- c(warnings, sprintf(
                    "Column '%s' has factor levels with leading/trailing spaces.",
                    col
                ))
            }
        } else if (is.character(values)) {
            trimmed_vals <- trimws(values)
            if (any(trimmed_vals != values, na.rm = TRUE)) {
                warnings <- c(warnings, sprintf(
                    "Column '%s' has values with leading/trailing spaces.",
                    col
                ))
            }
        }
    }

    if (batch_col %in% names(sample_annotation)) {
        batch_values <- sample_annotation[[batch_col]]
        batch_counts <- table(batch_values, useNA = "no")
        singletons <- names(batch_counts)[batch_counts == 1]
        if (length(singletons)) {
            warnings <- c(warnings, sprintf(
                "Single-sample batches detected in '%s': %s",
                batch_col,
                paste(singletons, collapse = ", ")
            ))
        }

        if (!is.null(condition_col) && condition_col %in% names(sample_annotation)) {
            idx <- !is.na(batch_values) & !is.na(sample_annotation[[condition_col]])
            if (any(idx)) {
                condition_values <- sample_annotation[[condition_col]][idx]
                batch_values_idx <- batch_values[idx]

                pairs <- split(condition_values, batch_values_idx)
                nested <- all(vapply(pairs, function(x) length(unique(x)) == 1, logical(1)))
                if (nested) {
                    warnings <- c(warnings, sprintf(
                        "Batches appear nested within condition '%s'.",
                        condition_col
                    ))
                }

                tab <- table(batch_values_idx, condition_values, useNA = "no")
                batch_conditions <- rowSums(tab > 0)
                condition_batches <- colSums(tab > 0)

                if (any(condition_batches <= 1)) {
                    warnings <- c(warnings, sprintf(
                        "Some condition levels occur in only one batch (possible confounding of '%s' with '%s').",
                        condition_col,
                        batch_col
                    ))
                }
                if (any(batch_conditions <= 1)) {
                    warnings <- c(warnings, sprintf(
                        "Some batches contain only one condition level (possible nesting of '%s' within '%s').",
                        batch_col,
                        condition_col
                    ))
                }
                if (all(condition_batches == 1) && all(batch_conditions == 1)) {
                    errors <- c(errors, sprintf(
                        "Condition '%s' and batch '%s' are completely confounded (one-to-one mapping).",
                        condition_col,
                        batch_col
                    ))
                }
            }
        }
    }

    summary <- list(
        n_samples = nrow(sample_annotation),
        batch_col = batch_col,
        condition_col = condition_col,
        covariates = covariates,
        n_batches = if (batch_col %in% names(sample_annotation)) {
            length(unique(sample_annotation[[batch_col]][!is.na(sample_annotation[[batch_col]])]))
        } else {
            NA_integer_
        },
        n_conditions = if (!is.null(condition_col) && condition_col %in% names(sample_annotation)) {
            length(unique(sample_annotation[[condition_col]][!is.na(sample_annotation[[condition_col]])]))
        } else {
            NA_integer_
        },
        sample_id_source = id_source,
        missing_columns = missing_cols
    )

    out <- list(
        errors = unique(errors),
        warnings = unique(warnings),
        summary = summary
    )
    class(out) <- "pb_design_check"

    if (isTRUE(strict) && length(out$errors)) {
        msg <- paste(
            "Design validation failed with the following errors:",
            paste(sprintf("- %s", out$errors), collapse = "\n"),
            sep = "\n"
        )
        stop(msg, call. = FALSE)
    }

    out
}

#' Summarise design matrix structure
#'
#' @inheritParams proBatch
#' @param condition_col optional column in \code{sample_annotation} indicating
#'   the primary biological condition.
#' @param covariates optional character vector of additional covariate columns.
#' @param add_intercept logical; include intercept in the design matrix.
#'
#' @return A list with design diagnostics: \code{design_matrix_dim},
#'   \code{rank}, \code{deficiency}, \code{kappa}, \code{aliased_terms},
#'   \code{flagged_covariates}, and \code{notes}. The object has class
#'   \code{"pb_design_summary"}.
#' @export
#'
#' @examples
#' data("example_ecoli_data", package = "proBatch")
#' meta <- example_ecoli_data$all_metadata
#' summarize_design(
#'     meta,
#'     batch_col = "Lab",
#'     condition_col = "Condition",
#'     add_intercept = TRUE
#' )
summarize_design <- function(sample_annotation,
                             batch_col,
                             condition_col = NULL,
                             covariates = NULL,
                             add_intercept = TRUE) {
    if (is.null(sample_annotation) || !is.data.frame(sample_annotation)) {
        stop("sample_annotation must be a data.frame.")
    }
    if (is.null(batch_col) || !nzchar(batch_col)) {
        stop("`batch_col` must be provided.")
    }

    sample_annotation <- as.data.frame(sample_annotation)
    design_cols <- .pb_design_collect_cols(
        batch_col = batch_col,
        condition_col = condition_col,
        covariates = covariates
    )

    missing_cols <- setdiff(design_cols, names(sample_annotation))
    if (length(missing_cols)) {
        stop(sprintf(
            "Missing columns in sample_annotation: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    design_df <- sample_annotation[, design_cols, drop = FALSE]
    complete_idx <- stats::complete.cases(design_df)
    notes <- character()
    if (any(!complete_idx)) {
        notes <- c(notes, sprintf(
            "Dropped %d samples with missing values in design variables.",
            sum(!complete_idx)
        ))
        design_df <- design_df[complete_idx, , drop = FALSE]
    }
    if (!nrow(design_df)) {
        stop("No samples remain after removing missing values from design variables.")
    }

    design_formula <- stats::reformulate(
        termlabels = design_cols,
        response = NULL,
        intercept = isTRUE(add_intercept)
    )
    design_matrix <- stats::model.matrix(design_formula, data = design_df)

    qr_obj <- base::qr(design_matrix)
    rank <- qr_obj$rank
    deficiency <- ncol(design_matrix) - rank
    aliased_terms <- character()
    if (deficiency > 0) {
        aliased_idx <- qr_obj$pivot[(rank + 1):ncol(design_matrix)]
        aliased_terms <- colnames(design_matrix)[aliased_idx]
        notes <- c(notes, "Rank deficiency detected; aliased terms reported.")
    }

    kappa_value <- tryCatch(
        base::kappa(design_matrix),
        error = function(e) NA_real_
    )
    if (!is.na(kappa_value) && is.finite(kappa_value) && kappa_value > 30) {
        notes <- c(notes, "High condition number indicates potential collinearity.")
    }

    flagged_covariates <- list()

    numeric_cols <- names(design_df)[vapply(design_df, is.numeric, logical(1))]
    high_corr_pairs <- data.frame()
    if (length(numeric_cols) >= 2) {
        cor_mat <- stats::cor(design_df[, numeric_cols, drop = FALSE],
            use = "pairwise.complete.obs"
        )
        cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA_real_
        idx <- which(abs(cor_mat) >= 0.9, arr.ind = TRUE)
        if (nrow(idx)) {
            high_corr_pairs <- data.frame(
                var1 = rownames(cor_mat)[idx[, 1]],
                var2 = colnames(cor_mat)[idx[, 2]],
                correlation = cor_mat[idx],
                stringsAsFactors = FALSE
            )
        }
    }

    factor_cols <- names(design_df)[vapply(design_df, function(x) {
        is.factor(x) || is.character(x)
    }, logical(1))]

    rare_levels <- data.frame()
    if (length(factor_cols)) {
        for (col in factor_cols) {
            f <- as.factor(design_df[[col]])
            counts <- table(f, useNA = "no")
            rare <- counts[counts <= 1]
            if (length(rare)) {
                rare_levels <- rbind(rare_levels, data.frame(
                    variable = col,
                    level = names(rare),
                    n = as.integer(rare),
                    stringsAsFactors = FALSE
                ))
            }
        }
    }

    sparse_cells <- data.frame()
    near_separation <- data.frame()
    if (length(factor_cols) >= 2) {
        pairs <- utils::combn(factor_cols, 2, simplify = FALSE)
        for (pair in pairs) {
            f1 <- as.factor(design_df[[pair[1]]])
            f2 <- as.factor(design_df[[pair[2]]])
            tab <- table(f1, f2, useNA = "no")
            zero_cells <- sum(tab == 0)
            if (zero_cells > 0) {
                sparse_cells <- rbind(sparse_cells, data.frame(
                    var1 = pair[1],
                    var2 = pair[2],
                    zero_cells = zero_cells,
                    total_cells = length(tab),
                    zero_fraction = zero_cells / length(tab),
                    stringsAsFactors = FALSE
                ))
            }
            row_nonzero <- rowSums(tab > 0)
            col_nonzero <- colSums(tab > 0)
            if (any(row_nonzero <= 1) || any(col_nonzero <= 1)) {
                near_separation <- rbind(near_separation, data.frame(
                    var1 = pair[1],
                    var2 = pair[2],
                    stringsAsFactors = FALSE
                ))
            }
        }
    }

    flagged_covariates$high_correlations <- high_corr_pairs
    flagged_covariates$rare_levels <- rare_levels
    flagged_covariates$sparse_cells <- sparse_cells
    flagged_covariates$near_separation <- near_separation

    out <- list(
        design_matrix_dim = c(n_samples = nrow(design_matrix), n_terms = ncol(design_matrix)),
        rank = rank,
        deficiency = deficiency,
        kappa = kappa_value,
        aliased_terms = aliased_terms,
        flagged_covariates = flagged_covariates,
        notes = unique(notes)
    )
    class(out) <- "pb_design_summary"
    out
}

#' Detect nested batch variables
#'
#' @param sample_annotation data frame with sample annotations.
#' @param batch_cols character vector with two or more batch variables.
#'
#' @return A list with adjacency matrix of nesting relations (parent -> child),
#'   edge list, equivalent (bijective) partitions, topological ordering (if
#'   acyclic), and suggestions. The object has class \code{"pb_nested_batches"}.
#' @export
#'
#' @examples
#' data("example_ecoli_data", package = "proBatch")
#' meta <- example_ecoli_data$all_metadata
#' detect_nested_batches(meta, batch_cols = c("Lab", "Condition"))
detect_nested_batches <- function(sample_annotation, batch_cols) {
    if (is.null(sample_annotation) || !is.data.frame(sample_annotation)) {
        stop("sample_annotation must be a data.frame.")
    }
    if (is.null(batch_cols) || length(batch_cols) < 2) {
        stop("`batch_cols` must contain at least two column names.")
    }

    sample_annotation <- as.data.frame(sample_annotation)
    missing_cols <- setdiff(batch_cols, names(sample_annotation))
    if (length(missing_cols)) {
        stop(sprintf(
            "Missing columns in sample_annotation: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }

    batch_cols <- unique(batch_cols)
    n <- length(batch_cols)
    adjacency <- matrix(FALSE,
        nrow = n, ncol = n,
        dimnames = list(batch_cols, batch_cols)
    )

    edges <- data.frame(parent = character(), child = character(), stringsAsFactors = FALSE)
    equivalent <- list()

    for (i in seq_len(n - 1)) {
        for (j in (i + 1):n) {
            col_a <- batch_cols[i]
            col_b <- batch_cols[j]
            idx <- !is.na(sample_annotation[[col_a]]) & !is.na(sample_annotation[[col_b]])
            if (!any(idx)) {
                next
            }

            a_in_b <- .pb_is_nested(
                child = sample_annotation[[col_a]][idx],
                parent = sample_annotation[[col_b]][idx]
            )
            b_in_a <- .pb_is_nested(
                child = sample_annotation[[col_b]][idx],
                parent = sample_annotation[[col_a]][idx]
            )

            if (a_in_b && b_in_a) {
                equivalent[[length(equivalent) + 1]] <- c(col_a, col_b)
                next
            }

            if (a_in_b) {
                adjacency[col_b, col_a] <- TRUE
                edges <- rbind(edges, data.frame(
                    parent = col_b,
                    child = col_a,
                    stringsAsFactors = FALSE
                ))
            }
            if (b_in_a) {
                adjacency[col_a, col_b] <- TRUE
                edges <- rbind(edges, data.frame(
                    parent = col_a,
                    child = col_b,
                    stringsAsFactors = FALSE
                ))
            }
        }
    }

    ordering <- .pb_toposort(adjacency)
    suggestion <- .pb_nested_suggestion(adjacency, equivalent = equivalent)

    out <- list(
        adjacency = adjacency,
        edges = edges,
        equivalent = equivalent,
        ordering = ordering,
        suggestion = suggestion
    )
    class(out) <- "pb_nested_batches"
    out
}

#' Detect outlier samples in PCA space
#'
#' @inheritParams proBatch
#' @param sample_annotation optional sample annotation data frame.
#' @param batch_col optional batch column to report per-batch summaries.
#' @param n_pcs number of principal components to use.
#' @param center logical; center variables before PCA.
#' @param scale. logical; scale variables before PCA.
#' @param robust logical; use robust covariance via \code{MASS::cov.rob}.
#' @param cutoff probability for the chi-square cutoff used to flag outliers.
#'
#' @return A data frame with \code{sample_id}, \code{batch} (if provided),
#'   \code{mdist}, and \code{is_outlier}. Attributes include \code{cutoff},
#'   \code{df}, and \code{mdist_cutoff}. The object has class
#'   \code{"pb_outliers"}.
#' @export
#'
#' @examples
#' data("example_ecoli_data", package = "proBatch")
#' meta <- example_ecoli_data$all_metadata
#' dm <- as.matrix(example_ecoli_data$all_protein_groups)[1:100, ]
#' outliers <- detect_outlier_samples(
#'     dm,
#'     sample_annotation = meta,
#'     batch_col = "Lab",
#'     n_pcs = 5
#' )
detect_outlier_samples <- function(data_matrix,
                                   sample_annotation = NULL,
                                   batch_col = NULL,
                                   n_pcs = 10,
                                   center = TRUE,
                                   scale. = TRUE,
                                   robust = TRUE,
                                   cutoff = 0.99) {
    if (is.null(data_matrix)) {
        stop("data_matrix must be provided.")
    }
    data_matrix <- .pb_matrix_for_pca(data_matrix)

    sample_ids <- colnames(data_matrix)
    if (is.null(sample_ids)) {
        stop("data_matrix must have column names for sample IDs.")
    }

    if (!is.numeric(n_pcs) || length(n_pcs) != 1 || n_pcs < 1) {
        stop("`n_pcs` must be a positive integer.")
    }
    if (!is.numeric(cutoff) || length(cutoff) != 1 || cutoff <= 0 || cutoff >= 1) {
        stop("`cutoff` must be between 0 and 1.")
    }

    scores_info <- .pb_pca_scores(
        data_matrix = data_matrix,
        n_pcs = n_pcs,
        center = center,
        scale. = scale.
    )
    scores <- scores_info$scores
    n_pcs <- scores_info$n_pcs

    cov_center <- colMeans(scores)
    cov_matrix <- stats::cov(scores)

    if (isTRUE(robust)) {
        if (requireNamespace("MASS", quietly = TRUE)) {
            cov_fit <- tryCatch(
                MASS::cov.rob(scores),
                error = function(e) NULL
            )
            if (!is.null(cov_fit)) {
                cov_center <- cov_fit$center
                cov_matrix <- cov_fit$cov
            } else {
                warning("Robust covariance estimation failed; falling back to cov().")
            }
        } else {
            warning("Package 'MASS' not available; falling back to cov().")
        }
    }

    mdist <- .pb_safe_mahalanobis(scores, cov_center, cov_matrix)
    mdist_cutoff <- stats::qchisq(cutoff, df = n_pcs)
    is_outlier <- mdist > mdist_cutoff

    out <- data.frame(
        sample_id = sample_ids,
        mdist = mdist,
        is_outlier = is_outlier,
        stringsAsFactors = FALSE
    )

    if (!is.null(batch_col) && !is.null(sample_annotation)) {
        annotation <- .pb_align_annotation_to_matrix(
            sample_annotation = sample_annotation,
            sample_ids = sample_ids
        )
        if (!(batch_col %in% names(annotation))) {
            stop(sprintf("batch_col '%s' not found in sample_annotation.", batch_col))
        }
        out$batch <- annotation[[batch_col]]
        out <- out[, c("sample_id", "batch", "mdist", "is_outlier"), drop = FALSE]
        batch_summary <- .pb_outlier_batch_summary(out)
        attr(out, "batch_summary") <- batch_summary
    }

    attr(out, "cutoff") <- cutoff
    attr(out, "df") <- n_pcs
    attr(out, "mdist_cutoff") <- mdist_cutoff
    class(out) <- c("pb_outliers", "data.frame")
    out
}

#' Detect sub-batches within batches
#'
#' @inheritParams proBatch
#' @param n_pcs number of principal components to use.
#' @param method clustering method; \code{"hclust"} or \code{"kmeans"}.
#' @param k_max maximum number of clusters to evaluate per batch.
#'
#' @return A list with per-sample assignments, per-batch summaries, and
#'   suggestions. The object has class \code{"pb_subbatches"}.
#' @export
#'
#' @examples
#' data("example_ecoli_data", package = "proBatch")
#' meta <- example_ecoli_data$all_metadata
#' dm <- as.matrix(example_ecoli_data$all_protein_groups)[1:200, ]
#' subbatches <- subbatch_detection(
#'     dm,
#'     sample_annotation = meta,
#'     batch_col = "Lab",
#'     n_pcs = 5,
#'     method = "kmeans",
#'     k_max = 4
#' )
subbatch_detection <- function(data_matrix,
                               sample_annotation,
                               batch_col,
                               n_pcs = 10,
                               method = c("hclust", "kmeans"),
                               k_max = 6) {
    if (is.null(sample_annotation) || !is.data.frame(sample_annotation)) {
        stop("sample_annotation must be a data.frame.")
    }
    if (is.null(batch_col) || !nzchar(batch_col)) {
        stop("`batch_col` must be provided.")
    }
    if (!is.numeric(k_max) || length(k_max) != 1 || k_max < 2) {
        stop("`k_max` must be an integer >= 2.")
    }
    method <- match.arg(method)

    data_matrix <- .pb_matrix_for_pca(data_matrix)
    sample_ids <- colnames(data_matrix)
    if (is.null(sample_ids)) {
        stop("data_matrix must have column names for sample IDs.")
    }

    annotation <- .pb_align_annotation_to_matrix(
        sample_annotation = sample_annotation,
        sample_ids = sample_ids
    )
    if (!(batch_col %in% names(annotation))) {
        stop(sprintf("batch_col '%s' not found in sample_annotation.", batch_col))
    }
    batch_values <- annotation[[batch_col]]

    keep <- !is.na(batch_values)
    if (!all(keep)) {
        warning("Samples with NA batch labels are removed before subbatch detection.")
        data_matrix <- data_matrix[, keep, drop = FALSE]
        sample_ids <- sample_ids[keep]
        batch_values <- batch_values[keep]
    }

    scores_info <- .pb_pca_scores(
        data_matrix = data_matrix,
        n_pcs = n_pcs,
        center = TRUE,
        scale. = TRUE
    )
    scores <- scores_info$scores

    assignments <- data.frame(
        sample_id = sample_ids,
        batch = batch_values,
        subbatch = NA_character_,
        stringsAsFactors = FALSE
    )

    summary_rows <- list()
    suggestions <- character()

    batches <- unique(batch_values)
    for (batch in batches) {
        idx <- which(batch_values == batch)
        n_batch <- length(idx)
        batch_scores <- scores[idx, , drop = FALSE]

        chosen_k <- 1
        cluster <- rep(1L, n_batch)
        within_disp <- 0

        if (n_batch >= 3) {
            max_k <- min(k_max, n_batch - 1)
            if (max_k >= 2) {
                ch_scores <- rep(NA_real_, max_k)
                clusters <- vector("list", max_k)

                for (k in 2:max_k) {
                    if (method == "kmeans") {
                        km <- stats::kmeans(batch_scores, centers = k, nstart = 10)
                        cluster_k <- km$cluster
                    } else {
                        hc <- stats::hclust(stats::dist(batch_scores), method = "ward.D2")
                        cluster_k <- stats::cutree(hc, k = k)
                    }
                    clusters[[k]] <- cluster_k
                    ch_scores[k] <- .pb_ch_index(batch_scores, cluster_k)
                }

                valid_idx <- which(is.finite(ch_scores))
                valid_idx <- valid_idx[valid_idx >= 2]
                if (length(valid_idx)) {
                    chosen_k <- valid_idx[which.max(ch_scores[valid_idx])]
                    cluster <- clusters[[chosen_k]]
                    within_disp <- .pb_within_dispersion(batch_scores, cluster)
                }
            }
        }

        sizes <- table(cluster)
        subbatch_labels <- paste0(batch, "_sub", as.integer(cluster))

        assignments$subbatch[idx] <- subbatch_labels

        summary_rows[[length(summary_rows) + 1]] <- data.frame(
            batch = batch,
            n_samples = n_batch,
            k = chosen_k,
            cluster_sizes = paste(names(sizes), as.integer(sizes), sep = ":", collapse = ", "),
            within_dispersion = within_disp,
            stringsAsFactors = FALSE
        )

        if (chosen_k > 1 && min(as.integer(sizes)) >= 2) {
            suggestions <- c(suggestions, sprintf(
                "Consider splitting batch %s into %d sub-batches.",
                batch,
                chosen_k
            ))
        }
    }

    summary_df <- if (length(summary_rows)) {
        do.call(rbind, summary_rows)
    } else {
        data.frame(
            batch = character(),
            n_samples = integer(),
            k = integer(),
            cluster_sizes = character(),
            within_dispersion = numeric(),
            stringsAsFactors = FALSE
        )
    }

    if (!length(suggestions)) {
        suggestions <- "No strong sub-batch structure detected."
    }

    out <- list(
        assignments = assignments,
        summary = summary_df,
        suggestions = suggestions
    )
    class(out) <- "pb_subbatches"
    out
}

.pb_design_collect_cols <- function(batch_col, condition_col, covariates) {
    cols <- c(batch_col, condition_col, covariates)
    cols <- cols[!is.na(cols) & nzchar(cols)]
    unique(cols)
}

.pb_has_default_rownames <- function(df) {
    if (is.null(df) || !is.data.frame(df)) {
        return(FALSE)
    }
    rn <- attr(df, "row.names")
    n <- nrow(df)
    is.integer(rn) && length(rn) == 2L &&
        is.na(rn[1L]) && identical(rn[2L], -n)
}

.pb_is_nested <- function(child, parent) {
    pairs <- split(parent, child)
    all(vapply(pairs, function(x) length(unique(x)) == 1, logical(1)))
}

.pb_align_annotation_to_matrix <- function(sample_annotation, sample_ids) {
    if (is.null(sample_annotation)) {
        return(NULL)
    }

    sample_annotation <- as.data.frame(sample_annotation)

    candidate_cols <- c("sample_id", "FullRunName", "Run")
    id_col <- NULL

    if (!is.null(rownames(sample_annotation)) &&
        !anyNA(rownames(sample_annotation)) &&
        setequal(sample_ids, rownames(sample_annotation))) {
        id_col <- NULL
    } else {
        matches <- names(sample_annotation)[vapply(sample_annotation, function(x) {
            vals <- as.character(x)
            !anyNA(vals) && length(unique(vals)) == length(vals) &&
                setequal(sample_ids, vals)
        }, logical(1))]

        if (length(matches)) {
            preferred <- intersect(candidate_cols, matches)
            if (length(preferred)) {
                id_col <- preferred[1]
            } else {
                id_col <- matches[1]
            }
            if (length(matches) > 1) {
                warning(sprintf(
                    "Multiple columns match sample IDs; using '%s'.",
                    id_col
                ))
            }
        } else if (!is.null(rownames(sample_annotation)) &&
            !anyNA(rownames(sample_annotation)) &&
            all(sample_ids %in% rownames(sample_annotation))) {
            id_col <- NULL
        } else {
            stop("Unable to align sample_annotation to data_matrix; provide rownames or a sample ID column matching data_matrix colnames.")
        }
    }

    .align_sample_annotation(sample_annotation, sample_ids, sample_id_col = id_col)
}

.pb_matrix_for_pca <- function(data_matrix) {
    mat <- as.matrix(data_matrix)
    if (!is.numeric(mat)) {
        stop("data_matrix must be numeric.")
    }
    if (anyNA(mat)) {
        mat <- mat[stats::complete.cases(mat), , drop = FALSE]
    }
    if (!nrow(mat) || !ncol(mat)) {
        stop("data_matrix has no complete features remaining for PCA.")
    }
    mat
}

.pb_pca_scores <- function(data_matrix, n_pcs, center, scale.) {
    max_pcs <- min(ncol(data_matrix) - 1, nrow(data_matrix))
    if (max_pcs < 1) {
        stop("Not enough samples/features for PCA.")
    }
    if (n_pcs > max_pcs) {
        warning(sprintf("Reducing n_pcs to %d to fit available data.", max_pcs))
        n_pcs <- max_pcs
    }
    pca <- stats::prcomp(t(data_matrix), center = center, scale. = scale.)
    scores <- pca$x
    scores <- scores[, seq_len(n_pcs), drop = FALSE]
    list(scores = scores, n_pcs = n_pcs)
}

.pb_safe_mahalanobis <- function(scores, center, cov_matrix) {
    mdist <- tryCatch(
        stats::mahalanobis(scores, center, cov_matrix),
        error = function(e) NA_real_
    )
    if (anyNA(mdist)) {
        if (requireNamespace("MASS", quietly = TRUE)) {
            inv_cov <- MASS::ginv(cov_matrix)
            mdist <- stats::mahalanobis(scores, center, inv_cov, inverted = TRUE)
        } else {
            stop("Mahalanobis distance failed; install MASS or provide a full-rank covariance matrix.")
        }
    }
    mdist
}

.pb_outlier_batch_summary <- function(outlier_df) {
    batches <- unique(outlier_df$batch)
    summary <- lapply(batches, function(batch) {
        idx <- outlier_df$batch == batch
        data.frame(
            batch = batch,
            n_samples = sum(idx),
            n_outliers = sum(outlier_df$is_outlier[idx]),
            pct_outliers = 100 * mean(outlier_df$is_outlier[idx]),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, summary)
}

.pb_ch_index <- function(scores, cluster) {
    k <- length(unique(cluster))
    n <- nrow(scores)
    if (k <= 1 || k >= n) {
        return(NA_real_)
    }
    overall <- colMeans(scores)
    within <- 0
    between <- 0
    for (cl in unique(cluster)) {
        idx <- cluster == cl
        center <- colMeans(scores[idx, , drop = FALSE])
        diff <- scores[idx, , drop = FALSE] -
            matrix(center, nrow = sum(idx), ncol = ncol(scores), byrow = TRUE)
        within <- within + sum(diff^2)
        between <- between + sum(idx) * sum((center - overall)^2)
    }
    if (within == 0) {
        return(NA_real_)
    }
    (between / (k - 1)) / (within / (n - k))
}

.pb_within_dispersion <- function(scores, cluster) {
    within <- 0
    for (cl in unique(cluster)) {
        idx <- cluster == cl
        center <- colMeans(scores[idx, , drop = FALSE])
        diff <- scores[idx, , drop = FALSE] -
            matrix(center, nrow = sum(idx), ncol = ncol(scores), byrow = TRUE)
        within <- within + sum(diff^2)
    }
    within
}

.pb_toposort <- function(adjacency) {
    if (is.null(adjacency) || !nrow(adjacency)) {
        return(NULL)
    }
    adjacency <- adjacency != 0
    nodes <- rownames(adjacency)
    incoming <- colSums(adjacency)
    order <- character()

    while (length(nodes)) {
        roots <- nodes[incoming[nodes] == 0]
        if (!length(roots)) {
            return(NULL)
        }
        order <- c(order, roots)
        for (node in roots) {
            children <- names(which(adjacency[node, ]))
            incoming[children] <- incoming[children] - 1
        }
        nodes <- setdiff(nodes, roots)
    }
    order
}

.pb_nested_suggestion <- function(adjacency, equivalent = NULL) {
    if (is.null(adjacency) || !nrow(adjacency)) {
        if (length(equivalent)) {
            eq_text <- vapply(equivalent, function(x) paste(x, collapse = " <-> "), character(1))
            return(sprintf(
                "Equivalent partitions detected: %s",
                paste(eq_text, collapse = "; ")
            ))
        }
        return("No nesting relations detected.")
    }
    incoming <- colSums(adjacency)
    outgoing <- rowSums(adjacency)
    higher <- names(which(outgoing > 0 & incoming == 0))
    lower <- names(which(incoming > 0 & outgoing == 0))
    intermediate <- names(which(incoming > 0 & outgoing > 0))
    standalone <- names(which(incoming == 0 & outgoing == 0))

    suggestions <- character()
    if (length(higher)) {
        suggestions <- c(suggestions, sprintf(
            "Higher-level candidates (often random effects): %s",
            paste(higher, collapse = ", ")
        ))
    }
    if (length(lower)) {
        suggestions <- c(suggestions, sprintf(
            "Lower-level candidates (often fixed/technical): %s",
            paste(lower, collapse = ", ")
        ))
    }
    if (length(intermediate)) {
        suggestions <- c(suggestions, sprintf(
            "Intermediate nesting levels: %s",
            paste(intermediate, collapse = ", ")
        ))
    }
    if (length(standalone)) {
        suggestions <- c(suggestions, sprintf(
            "No nesting detected for: %s",
            paste(standalone, collapse = ", ")
        ))
    }
    if (length(equivalent)) {
        eq_text <- vapply(equivalent, function(x) paste(x, collapse = " <-> "), character(1))
        suggestions <- c(suggestions, sprintf(
            "Equivalent partitions detected: %s",
            paste(eq_text, collapse = "; ")
        ))
    }
    if (!length(suggestions)) {
        suggestions <- "No nesting relations detected."
    }
    suggestions
}
