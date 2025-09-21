#' Plot missing-value heatmap(s)
#'
#' Functions for visualising the missingness pattern of assay intensities as a binary heatmap.
#' The `ProBatchFeatures` method supports drawing multiple assays at once by
#' arranging the resulting heatmaps into a user-controlled grid layout.
#'
#' @param x A data container. For the `ProBatchFeatures` method this must be a
#'   `ProBatchFeatures` object. The default method accepts any matrix-like input
#'   (including `SummarizedExperiment`).
#' @param pbf_name Character scalar or vector with assay names to plot. When
#'   `NULL`, the most recent assay returned by [pb_current_assay()] is used.
#'   Only used by the `ProBatchFeatures` method.
#' @param sample_annotation Optional data frame with sample-level metadata. Row
#'   names (or the column specified via `sample_id_col`) must match the column
#'   names of the intensity matrix. When `x` is a `ProBatchFeatures` object the
#'   sample annotation defaults to `as.data.frame(colData(x))`.
#' @param sample_id_col Optional column in `sample_annotation` providing unique
#'   sample identifiers. Use this when the data frame lacks row names matching
#'   the assay column names.
#' @param color_by Optional column name in `sample_annotation` used to annotate
#'   heatmap columns. Use `NULL` (default) or the string "No" to omit the
#'   annotation bar.
#' @param label_by Optional column name (or character vector) used for column
#'   labels. Use `NULL` for default assay column names or the string "No" to
#'   suppress column labels entirely.
#' @param cluster_samples,cluster_features Logical flags controlling whether the
#'   heatmap columns/rows are clustered.
#' @param show_row_dend,show_column_dend Logical, whether dendrograms should be
#'   drawn for the clustered rows/columns.
#' @param missing_color,valid_color Colours used for missing (`0`) and observed
#'   (`1`) values respectively.
#' @param col_vector Optional vector of colours that will be recycled to colour
#'   the unique values of `color_by`.
#' @param drop_complete Logical, drop features without any missing values before
#'   plotting. Defaults to `TRUE` to focus on missingness patterns.
#' @param nrow,ncol Integers controlling the layout when multiple assays are
#'   plotted. If both are `NULL`, a roughly square layout is chosen
#'   automatically.
#' @param draw Logical, draw the heatmap(s). Set to `FALSE` to obtain the grob(s)
#'   without plotting. For multiple assays, the arranged grob is returned
#'   invisibly when `draw = TRUE`.
#' @param ... Additional parameters forwarded to [pheatmap::pheatmap()].
#'
#' @return For a single assay the returned value is the `pheatmap` object. When
#'   multiple assays are requested a list is returned invisibly with elements
#'   `grob` (the arranged heatmaps) and `heatmaps` (individual `pheatmap`
#'   objects). Assays without missing values are skipped with a warning.
#' @export
plot_NA_heatmap <- function(x, ...) UseMethod("plot_NA_heatmap")

#' @rdname plot_NA_heatmap
#' @method plot_NA_heatmap default
#' @export
plot_NA_heatmap.default <- function(
    x,
    sample_annotation = NULL,
    sample_id_col = NULL,
    color_by = NULL,
    label_by = NULL,
    cluster_samples = TRUE,
    cluster_features = TRUE,
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    missing_color = "black",
    valid_color = "grey90",
    col_vector = NULL,
    drop_complete = TRUE,
    draw = TRUE,
    main = NULL,
    ...) {
    data_matrix <- x
    plot_params <- list(...)

    if (methods::is(data_matrix, "SummarizedExperiment")) {
        if (is.null(sample_annotation)) {
            sample_annotation <- as.data.frame(SummarizedExperiment::colData(data_matrix))
        }
        data_matrix <- SummarizedExperiment::assay(data_matrix)
    }

    if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
    }

    if (!nrow(data_matrix) || !ncol(data_matrix)) {
        warning("Input matrix has zero rows or columns; nothing to plot.")
        return(NULL)
    }

    binary_matrix <- .pb_binary_missing_matrix(data_matrix, drop_complete = drop_complete)
    if (nrow(binary_matrix) == 0L) {
        msg <- if (is.null(main)) {
            "No rows with missing values remain after filtering."
        } else {
            sprintf("Assay '%s' has no rows with missing values after filtering.", main)
        }
        warning(msg)
        return(NULL)
    }

    sample_order <- colnames(binary_matrix)
    ann_info <- .pb_prepare_sample_annotation(
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        color_by = color_by,
        label_by = label_by,
        sample_order = sample_order,
        col_vector = col_vector
    )

    main <- if (is.null(main)) NA else main
    treeheight_row <- if (cluster_features && show_row_dend) 50 else 0
    treeheight_col <- if (cluster_samples && show_column_dend) 50 else 0

    annotation_col <- ann_info$annotation_col
    annotation_colors <- ann_info$annotation_colors
    labels_col <- ann_info$labels_col
    show_column_names <- ann_info$show_column_names

    if (!"fontsize" %in% names(plot_params)) plot_params$fontsize <- 10
    if (!"fontsize_row" %in% names(plot_params)) plot_params$fontsize_row <- 0.6 * plot_params$fontsize
    if (!"fontsize_col" %in% names(plot_params)) plot_params$fontsize_col <- 0.6 * plot_params$fontsize

    res <- do.call(pheatmap::pheatmap, c(
        list(
            mat = binary_matrix,
            cluster_rows = cluster_features,
            cluster_cols = cluster_samples,
            show_colnames = show_column_names,
            labels_col = labels_col,
            color = c(missing_color, valid_color),
            breaks = c(-0.5, 0.5, 1.5),
            legend_breaks = c(0, 1),
            legend_labels = c("Missing", "Valid"),
            annotation_col = annotation_col,
            annotation_colors = annotation_colors,
            treeheight_row = treeheight_row,
            treeheight_col = treeheight_col,
            silent = !draw,
            main = main
        ),
        plot_params
    ))
    if (draw && isTRUE(res$silent)) {
        grid::grid.newpage()
        grid::grid.draw(res$gtable)
    }

    res
}

#' @rdname plot_NA_heatmap
#' @method plot_NA_heatmap ProBatchFeatures
#' @export
plot_NA_heatmap.ProBatchFeatures <- function(
    x,
    pbf_name = NULL,
    color_by = NULL,
    label_by = NULL,
    sample_id_col = NULL,
    cluster_samples = TRUE,
    cluster_features = TRUE,
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    missing_color = "black",
    valid_color = "grey90",
    col_vector = NULL,
    drop_complete = TRUE,
    nrow = NULL,
    ncol = NULL,
    draw = TRUE,
    ...) {
    object <- x
    assays <- if (is.null(pbf_name)) pb_current_assay(object) else pbf_name

    if (!length(assays)) {
        stop("Provide at least one `pbf_name` to plot.")
    }

    sample_annotation <- as.data.frame(SummarizedExperiment::colData(object))
    if (!is.null(sample_id_col)) {
        if (!sample_id_col %in% names(sample_annotation)) {
            stop("Column '", sample_id_col, "' not found in `colData(object)`.")
        }
        sample_annotation[[sample_id_col]] <- as.character(sample_annotation[[sample_id_col]])
        rownames(sample_annotation) <- sample_annotation[[sample_id_col]]
    } else if (is.null(rownames(sample_annotation)) || any(is.na(rownames(sample_annotation)))) {
        reference_matrix <- pb_assay_matrix(object, assay = assays[[1]])
        rownames(sample_annotation) <- colnames(reference_matrix)
    }

    heatmaps <- vector("list", length(assays))
    names(heatmaps) <- assays
    produced <- logical(length(assays))

    for (idx in seq_along(assays)) {
        assay_nm <- assays[[idx]]
        data_matrix <- pb_assay_matrix(object, assay = assay_nm)

        # if number of rows or columns is bigger than 5000, select a random subset of 5000
        if (nrow(data_matrix) > 5000) {
            set.seed(123)
            row_idx <- sort(sample(seq_len(nrow(data_matrix)), 5000))
            data_matrix <- data_matrix[row_idx, , drop = FALSE]
            warning("Assay '", assay_nm, "' has more than 5000 rows; plotting a random subset of 5000 rows.")
        }
        if (ncol(data_matrix) > 5000) {
            set.seed(123)
            col_idx <- sort(sample(seq_len(ncol(data_matrix)), 5000))
            data_matrix <- data_matrix[, col_idx, drop = FALSE]
            warning("Assay '", assay_nm, "' has more than 5000 columns; plotting a random subset of 5000 columns.")
        }



        res <- plot_NA_heatmap.default(
            data_matrix,
            sample_annotation = sample_annotation,
            sample_id_col = NULL,
            color_by = color_by,
            label_by = label_by,
            cluster_samples = cluster_samples,
            cluster_features = cluster_features,
            show_row_dend = show_row_dend,
            show_column_dend = show_column_dend,
            missing_color = missing_color,
            valid_color = valid_color,
            col_vector = col_vector,
            drop_complete = drop_complete,
            draw = draw && length(assays) == 1L,
            main = assay_nm,
            ...
        )
        if (is.null(res)) {
            warning("Skipping assay '", assay_nm, "' because it has no rows with missing values after filtering.")
        } else {
            heatmaps[[idx]] <- res
            produced[[idx]] <- TRUE
        }
    }

    heatmaps <- heatmaps[produced]
    if (!length(heatmaps)) {
        return(invisible(NULL))
    }

    if (length(heatmaps) == 1L) {
        return(heatmaps[[1L]])
    }

    if (!requireNamespace("gridExtra", quietly = TRUE)) {
        stop("Install the `gridExtra` package to arrange multiple heatmaps: install.packages(\"gridExtra\").")
    }

    layout <- .pb_missing_layout(length(heatmaps), nrow = nrow, ncol = ncol)
    grob_list <- lapply(heatmaps, function(ht) ht$gtable)
    arranged <- do.call(
        gridExtra::arrangeGrob,
        c(list(grobs = grob_list, nrow = layout$nrow, ncol = layout$ncol), list())
    )
    if (draw) {
        grid::grid.newpage()
        grid::grid.draw(arranged)
    }

    invisible(list(grob = arranged, heatmaps = heatmaps))
}

#' Plot intensity density by missingness
#'
#' Compare the distribution of average intensities between features with and
#' without missing observations.
#'
#' @inheritParams plot_NA_heatmap
#' @param missing_label,valid_label Labels used to distinguish rows with and
#'   without missing values.
#' @param palette Named vector of colours mapped to `missing_label` and
#'   `valid_label`.
#' @param facet_scales Scaling behaviour passed to [ggplot2::facet_wrap()] when
#'   multiple assays are plotted.
#'
#' @return A `ggplot` object.
#' @export
plot_NA_density <- function(x, ...) UseMethod("plot_NA_density")

#' @rdname plot_NA_density
#' @method plot_NA_density default
#' @export
plot_NA_density.default <- function(
    x,
    missing_label = "Missing Value",
    valid_label = "Valid Value",
    palette = c(`Missing Value` = "#A92C23", `Valid Value` = "#345995"),
    ...) {
    data_matrix <- x
    if (methods::is(data_matrix, "SummarizedExperiment")) {
        data_matrix <- SummarizedExperiment::assay(data_matrix)
    }
    if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
    }
    if (!nrow(data_matrix) || !ncol(data_matrix)) {
        warning("Input matrix has zero rows or columns; nothing to plot.")
        return(ggplot2::ggplot())
    }

    has_observation <- rowSums(!is.na(data_matrix)) > 0
    if (!any(has_observation)) {
        warning("All rows contain only missing values; nothing to plot.")
        return(ggplot2::ggplot())
    }

    row_means <- rowMeans(data_matrix[has_observation, , drop = FALSE], na.rm = TRUE)
    missing_flag <- apply(data_matrix[has_observation, , drop = FALSE], 1, function(v) any(is.na(v)))
    df <- data.frame(
        mean = row_means,
        Type = ifelse(missing_flag, missing_label, valid_label),
        stringsAsFactors = FALSE
    )
    df <- df[is.finite(df$mean), , drop = FALSE]
    if (!nrow(df)) {
        warning("No finite mean intensities available for plotting.")
        return(ggplot2::ggplot())
    }

    palette <- .pb_match_palette(palette, c(missing_label, valid_label))

    ggplot2::ggplot(df, ggplot2::aes(x = mean, colour = Type)) +
        ggplot2::geom_density(na.rm = TRUE) +
        ggplot2::labs(x = "Intensity", y = "Density", colour = "Value Type") +
        ggplot2::scale_colour_manual(values = palette, breaks = c(missing_label, valid_label))
}

#' @rdname plot_NA_density
#' @method plot_NA_density ProBatchFeatures
#' @export
plot_NA_density.ProBatchFeatures <- function(
    x,
    pbf_name = NULL,
    missing_label = "Missing Value",
    valid_label = "Valid Value",
    palette = c(`Missing Value` = "#A92C23", `Valid Value` = "#345995"),
    nrow = NULL,
    ncol = NULL,
    facet_scales = "free_y",
    ...) {
    object <- x
    assays <- if (is.null(pbf_name)) pb_current_assay(object) else pbf_name
    if (!length(assays)) {
        stop("Provide at least one `pbf_name` to plot.")
    }

    df_list <- lapply(assays, function(assay_nm) {
        mat <- pb_assay_matrix(object, assay = assay_nm)
        .pb_missing_density_df(mat, assay_nm, missing_label, valid_label)
    })
    keep <- vapply(df_list, function(df) !is.null(df) && nrow(df) > 0, logical(1))
    if (!any(keep)) {
        warning("No finite mean intensities available across the requested assays.")
        return(ggplot2::ggplot())
    }

    combined <- do.call(rbind, df_list[keep])
    combined$pbf_name <- factor(combined$pbf_name, levels = assays[keep])

    palette <- .pb_match_palette(palette, c(missing_label, valid_label))
    p <- ggplot2::ggplot(combined, ggplot2::aes(x = mean, colour = Type)) +
        ggplot2::geom_density(na.rm = TRUE) +
        ggplot2::labs(x = "Intensity", y = "Density", colour = "Value Type") +
        ggplot2::scale_colour_manual(values = palette, breaks = c(missing_label, valid_label))

    if (length(unique(combined$pbf_name)) > 1L) {
        layout <- .pb_missing_layout(length(unique(combined$pbf_name)), nrow = nrow, ncol = ncol)
        p <- p + ggplot2::facet_wrap(~pbf_name, nrow = layout$nrow, ncol = layout$ncol, scales = facet_scales)
    }

    p
}

#' Plot missing-value frequency distribution
#'
#' Display how many features are observed in a given number of samples.
#'
#' @inheritParams plot_NA_heatmap
#' @param fill Colour used for the columns in the frequency plot.
#' @param facet_scales Scaling behaviour for facets when plotting multiple
#'   assays.
#'
#' @return A `ggplot` object showing the frequency distribution.
#' @export
plot_NA_frequency <- function(x, ...) UseMethod("plot_NA_frequency")

#' @rdname plot_NA_frequency
#' @method plot_NA_frequency default
#' @export
plot_NA_frequency.default <- function(
    x,
    show_percent = FALSE,
    fill = "#345995",
    ...) {
    data_matrix <- x
    if (methods::is(data_matrix, "SummarizedExperiment")) {
        data_matrix <- SummarizedExperiment::assay(data_matrix)
    }
    if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
    }
    if (!nrow(data_matrix) || !ncol(data_matrix)) {
        warning("Input matrix has zero rows or columns; nothing to plot.")
        return(ggplot2::ggplot())
    }

    valid_counts <- rowSums(!is.na(data_matrix))
    freq_df <- as.data.frame(table(valid_counts), stringsAsFactors = FALSE)
    if (!nrow(freq_df)) {
        warning("No frequency data available.")
        return(ggplot2::ggplot())
    }
    freq_df$valid_counts <- as.integer(as.character(freq_df$valid_counts))

    if (show_percent) {
        total <- sum(freq_df$Freq)
        if (total == 0) {
            warning("Total count is zero; cannot compute percentages.")
            return(ggplot2::ggplot())
        }
        freq_df$Percent <- as.numeric(freq_df$Freq) / total * 100
    }

    ggplot2::ggplot(
        freq_df, ggplot2::aes_string(
            x = "valid_counts",
            y = if (show_percent) "Percent" else "Freq"
        )
    ) +
        ggplot2::geom_col(fill = fill) +
        ggplot2::labs(
            x = "Identified in Number of Samples",
            y = if (show_percent) "Percent of Features" else "Number of Features"
        )
}

#' @rdname plot_NA_frequency
#' @method plot_NA_frequency ProBatchFeatures
#' @export
plot_NA_frequency.ProBatchFeatures <- function(
    x,
    pbf_name = NULL,
    fill = "#345995",
    nrow = NULL,
    ncol = NULL,
    facet_scales = "free_y",
    show_percent = FALSE,
    ...) {
    object <- x
    assays <- if (is.null(pbf_name)) pb_current_assay(object) else pbf_name
    if (!length(assays)) {
        stop("Provide at least one `pbf_name` to plot.")
    }

    df_list <- lapply(assays, function(assay_nm) {
        mat <- pb_assay_matrix(object, assay = assay_nm)
        .pb_missing_frequency_df(mat, assay_nm)
    })
    keep <- vapply(df_list, function(df) !is.null(df) && nrow(df) > 0, logical(1))
    if (!any(keep)) {
        warning("No frequency data available for the requested assays.")
        return(ggplot2::ggplot())
    }

    if (show_percent) {
        df_list <- lapply(df_list, function(df) {
            if (is.null(df) || nrow(df) == 0) {
                return(NULL)
            }
            total <- sum(df$count)
            if (total == 0) {
                return(NULL)
            }
            df$Percent <- as.numeric(df$count) / total * 100
            df
        })
        # Recompute which assays still have valid data after percent conversion
        keep <- vapply(df_list, function(df) !is.null(df) && nrow(df) > 0, logical(1))
        if (!any(keep)) {
            warning("No frequency data available for the requested assays.")
            return(ggplot2::ggplot())
        }
    }

    combined <- do.call(rbind, df_list[keep])
    combined$pbf_name <- factor(combined$pbf_name, levels = assays[keep])

    p <- ggplot2::ggplot(combined, ggplot2::aes(
        x = valid_counts,
        y = if (show_percent) Percent else count
    )) +
        ggplot2::geom_col(fill = fill) +
        ggplot2::labs(
            x = "Identified in Number of Samples",
            y = if (show_percent) "Percent" else "Number of Features"
        )

    if (length(unique(combined$pbf_name)) > 1L) {
        layout <- .pb_missing_layout(length(unique(combined$pbf_name)), nrow = nrow, ncol = ncol)
        p <- p + ggplot2::facet_wrap(~pbf_name, nrow = layout$nrow, ncol = layout$ncol, scales = facet_scales)
    }

    p
}

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

.pb_binary_missing_matrix <- function(data_matrix, drop_complete = TRUE) {
    binary <- ifelse(is.na(data_matrix), 0, 1)
    if (drop_complete) {
        keep <- rowSums(binary) < ncol(binary)
        binary <- binary[keep, , drop = FALSE]
    }
    binary
}

.pb_prepare_sample_annotation <- function(sample_annotation,
                                          sample_id_col,
                                          color_by,
                                          label_by,
                                          sample_order,
                                          col_vector) {
    if (is.null(sample_annotation)) {
        sample_annotation <- data.frame(row.names = sample_order)
    } else {
        sample_annotation <- as.data.frame(sample_annotation, stringsAsFactors = FALSE)
        if (!is.null(sample_id_col)) {
            if (!sample_id_col %in% names(sample_annotation)) {
                stop("Column '", sample_id_col, "' not found in `sample_annotation`.")
            }
            sample_annotation[[sample_id_col]] <- as.character(sample_annotation[[sample_id_col]])
            rownames(sample_annotation) <- sample_annotation[[sample_id_col]]
        }
    }

    if (is.null(rownames(sample_annotation))) {
        stop("`sample_annotation` must have row names that match the sample order or provide `sample_id_col`.")
    }

    missing_samples <- setdiff(sample_order, rownames(sample_annotation))
    if (length(missing_samples)) {
        stop("`sample_annotation` is missing entries for: ", paste(missing_samples, collapse = ", "))
    }

    sample_annotation <- sample_annotation[sample_order, , drop = FALSE]

    annotation_col <- NA
    annotation_colors <- NA
    labels_col <- NULL
    show_column_names <- TRUE

    if (!is.null(color_by) && !(identical(tolower(color_by), "no") && length(color_by) == 1L)) {
        if (!color_by %in% colnames(sample_annotation)) {
            stop("Column '", color_by, "' not found in `sample_annotation`.")
        }
        annotation_col <- data.frame(sample_annotation[[color_by]], row.names = sample_order, stringsAsFactors = FALSE)
        colnames(annotation_col) <- color_by
        levels_values <- unique(annotation_col[[color_by]])
        annotation_colors <- list()
        annotation_colors[[color_by]] <- .pb_build_annotation_colors(levels_values, col_vector)
    }

    if (!is.null(label_by)) {
        if (isFALSE(label_by) || (length(label_by) == 1L && identical(tolower(label_by), "no"))) {
            show_column_names <- FALSE
        } else if (length(label_by) == 1L) {
            if (!label_by %in% colnames(sample_annotation)) {
                stop("Column '", label_by, "' not found in `sample_annotation`.")
            }
            labels_col <- as.character(sample_annotation[[label_by]])
        } else if (length(label_by) == length(sample_order)) {
            labels_col <- as.character(label_by)
        } else {
            stop("`label_by` must be a column name or a vector with length equal to the number of samples.")
        }

        if (!is.null(labels_col) && anyDuplicated(labels_col)) {
            warning("Duplicate column labels detected; consider providing unique labels via `label_by`.")
        }
    }

    list(
        annotation_col = annotation_col,
        annotation_colors = annotation_colors,
        labels_col = labels_col,
        show_column_names = show_column_names
    )
}

.pb_build_annotation_colors <- function(values, col_vector = NULL) {
    values <- unique(as.character(values))
    n <- length(values)
    if (n == 0L) {
        return(character())
    }
    if (is.null(col_vector)) {
        qual_info <- RColorBrewer::brewer.pal.info
        qual_info <- qual_info[qual_info$category == "qual", , drop = FALSE]
        palette <- unlist(mapply(RColorBrewer::brewer.pal, qual_info$maxcolors, rownames(qual_info)))
        palette <- rev(palette)
        if (length(palette) < n) {
            palette <- rep_len(palette, n)
        }
    } else {
        palette <- rep_len(col_vector, n)
    }
    stats::setNames(palette[seq_len(n)], values)
}

.pb_missing_layout <- function(k, nrow = NULL, ncol = NULL) {
    if (k <= 0L) {
        stop("Layout requires at least one panel.")
    }
    if (is.null(nrow) && is.null(ncol)) {
        ncol <- ceiling(sqrt(k))
        nrow <- ceiling(k / ncol)
    } else if (is.null(nrow)) {
        nrow <- ceiling(k / ncol)
    } else if (is.null(ncol)) {
        ncol <- ceiling(k / nrow)
    }
    list(nrow = max(1L, as.integer(nrow)), ncol = max(1L, as.integer(ncol)))
}

.pb_missing_density_df <- function(data_matrix, assay_nm, missing_label, valid_label) {
    if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
    }
    if (!nrow(data_matrix) || !ncol(data_matrix)) {
        return(NULL)
    }
    has_observation <- rowSums(!is.na(data_matrix)) > 0
    if (!any(has_observation)) {
        return(NULL)
    }
    row_means <- rowMeans(data_matrix[has_observation, , drop = FALSE], na.rm = TRUE)
    missing_flag <- apply(data_matrix[has_observation, , drop = FALSE], 1, function(v) any(is.na(v)))
    df <- data.frame(
        mean = row_means,
        Type = ifelse(missing_flag, missing_label, valid_label),
        pbf_name = assay_nm,
        stringsAsFactors = FALSE
    )
    df[is.finite(df$mean), , drop = FALSE]
}

.pb_missing_frequency_df <- function(data_matrix, assay_nm) {
    if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
    }
    if (!nrow(data_matrix) || !ncol(data_matrix)) {
        return(NULL)
    }
    valid_counts <- rowSums(!is.na(data_matrix))
    freq_df <- as.data.frame(table(valid_counts), stringsAsFactors = FALSE)
    if (!nrow(freq_df)) {
        return(NULL)
    }
    freq_df$valid_counts <- as.integer(as.character(freq_df$valid_counts))
    data.frame(
        valid_counts = freq_df$valid_counts,
        count = freq_df$Freq,
        pbf_name = assay_nm,
        stringsAsFactors = FALSE
    )
}

.pb_match_palette <- function(palette, values) {
    if (is.null(palette)) {
        return(NULL)
    }
    if (is.null(names(palette))) {
        names(palette) <- values[seq_along(palette)]
    }
    missing <- setdiff(values, names(palette))
    if (length(missing)) {
        palette <- c(palette, stats::setNames(rep_len(palette, length(missing)), missing))
    }
    palette[values]
}
