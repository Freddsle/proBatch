#' cluster the data matrix to visually inspect which confounder dominates
#'
#' @inheritParams proBatch
#' @param distance distance metric used for clustering
#' @param agglomeration agglomeration methods as used by \code{hclust}
#' @param label_samples if \code{TRUE} sample IDs (column names of
#' \code{data_matrix}) will be printed
#' @param label_font size of the font. Is active if \code{label_samples} is
#' \code{TRUE}, ignored otherwise
#' @param fill_the_missing numeric value determining how  missing values
#' should be substituted. If \code{NULL}, features with missing values are
#' excluded.
#' @param ... other parameters of \code{plotDendroAndColors} from
#' \code{WGCNA} package
#'
#' @name plot_hierarchical_clustering
#'
#' @return No return
#' @examples
#' # Load necessary datasets
#' data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#'
#' selected_batches <- example_sample_annotation$MS_batch %in%
#'     c("Batch_1", "Batch_2")
#' selected_samples <- example_sample_annotation$FullRunName[selected_batches]
#' test_matrix <- example_proteome_matrix[, selected_samples]
#'
#' hierarchical_clustering_plot <- plot_hierarchical_clustering(
#'     example_proteome_matrix, example_sample_annotation,
#'     factors_to_plot = c("MS_batch", "Diet", "DateTime"),
#'     color_list = NULL,
#'     distance = "euclidean", agglomeration = "complete",
#'     label_samples = FALSE
#' )
#'
#' # with defined color scheme:
#' color_list <- sample_annotation_to_colors(example_sample_annotation,
#'     factor_columns = c("MS_batch", "Strain", "Diet", "digestion_batch"),
#'     numeric_columns = c("DateTime", "order")
#' )
#' hierarchical_clustering_plot <- plot_hierarchical_clustering(
#'     example_proteome_matrix, example_sample_annotation,
#'     factors_to_plot = c("MS_batch", "Strain", "DateTime", "digestion_batch"),
#'     color_list = color_list,
#'     distance = "euclidean", agglomeration = "complete",
#'     label_samples = FALSE
#' )
#'
#' @export
#'
#' @seealso \code{\link[stats]{hclust}},
#'   \code{\link{sample_annotation_to_colors}},
#'   \code{\link[WGCNA]{plotDendroAndColors}}
plot_hierarchical_clustering.default <- function(data_matrix, sample_annotation,
                                                 sample_id_col = "FullRunName",
                                                 color_list = NULL,
                                                 factors_to_plot = NULL,
                                                 fill_the_missing = 0,
                                                 distance = "euclidean",
                                                 agglomeration = "complete",
                                                 label_samples = TRUE, label_font = .2,
                                                 filename = NULL,
                                                 width = 38, height = 25,
                                                 units = c("cm", "in", "mm"),
                                                 plot_title = NULL,
                                                 ...) {
    df_long <- matrix_to_long(data_matrix, sample_id_col = sample_id_col)
    df_long <- check_sample_consistency(sample_annotation, sample_id_col, df_long,
        merge = FALSE
    )
    data_matrix <- long_to_matrix(df_long, sample_id_col = sample_id_col)
    rm(df_long)

    warning_message <- "Hierarchical clustering cannot operate with missing values
                      in the matrix"
    data_matrix <- handle_missing_values(
        data_matrix, warning_message,
        fill_the_missing
    )

    dist_matrix <- dist(t(as.matrix(data_matrix)), method = distance)
    hierarchical_clust <- hclust(dist_matrix, method = agglomeration)
    if (label_samples) {
        cex.dendroLabels <- label_font
        if (ncol(data_matrix) > 80) {
            warning("Too many samples, adjust the font with `label_font` argument or
              remove labels by setting `label_samples = FALSE` in
              function call")
        }
    } else {
        cex.dendroLabels <- 0.9
    }

    factors_without_colors <- setdiff(factors_to_plot, names(color_list))
    if (length(factors_without_colors) > 0) {
        warning("color_list for samples annotation not defined, inferring
             automatically. Numeric/factor columns are guessed, for more
            controlled color mapping use sample_annotation_to_colors()")
        color_list_new <- sample_annotation_to_colors(sample_annotation,
            sample_id_col = sample_id_col,
            factor_columns = factors_without_colors,
            numeric_columns = NULL
        )
        color_list <- c(color_list, color_list_new)
    }


    if (length(setdiff(names(color_list), factors_to_plot)) > 0 &&
        !is.null(factors_to_plot)) {
        color_list <- color_list[factors_to_plot]
    }

    # transform color list to color_df
    sample_annotation <- sample_annotation %>%
        filter(!!sym(sample_id_col) %in% colnames(data_matrix))
    color_df <- color_list_to_df(color_list, sample_annotation, sample_id_col)

    if (is.null(filename)) {
        plotDendroAndColors(hierarchical_clust, color_df,
            rowTextAlignment = "left",
            main = plot_title,
            hang = -0.1, addGuide = TRUE, ,
            dendroLabels = if (label_samples) NULL else FALSE,
            cex.dendroLabels = cex.dendroLabels,
            ...
        )
    } else {
        units_adjusted <- adjust_units(units, width, height)
        units <- units_adjusted$unit
        width <- units_adjusted$width
        height <- units_adjusted$height

        if (is.na(width)) {
            width <- 7
        }
        if (is.na(height)) {
            height <- 7
        }
        if (file_ext(filename) == "pdf") {
            pdf(file = filename, width = width, height = height, title = plot_title)
        } else if (file_ext(filename) == "png") {
            png(
                filename = filename, width = width, height = height, units = units,
                res = 300
            )
        } else {
            stop("currently only pdf and png extensions for filename are implemented")
        }

        plotDendroAndColors(hierarchical_clust, color_df,
            rowTextAlignment = "left",
            main = plot_title,
            hang = -0.1, addGuide = TRUE,
            dendroLabels = if (label_samples) NULL else FALSE,
            cex.dendroLabels = cex.dendroLabels,
            ...
        )

        dev.off()
    }
}

#' @rdname plot_hierarchical_clustering
#' @method plot_hierarchical_clustering ProBatchFeatures
#' @export
plot_hierarchical_clustering.ProBatchFeatures <- function(x, pbf_name = NULL,
                                                          sample_annotation = NULL,
                                                          sample_id_col = "FullRunName",
                                                          plot_title = NULL,
                                                          ...) {
    object <- x
    assay_name <- if (is.null(pbf_name)) pb_current_assay(object) else pbf_name
    data_matrix <- pb_assay_matrix(object, pbf_name)
    if (is.null(sample_annotation)) {
        sample_annotation <- as.data.frame(colData(object))
    }

    plot_title <- if (is.null(plot_title)) assay_name else plot_title

    plot_hierarchical_clustering.default(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        plot_title = plot_title,
        ...
    )
}

#' @export
plot_hierarchical_clustering <- function(x, ...) UseMethod("plot_hierarchical_clustering")

#' Plot the heatmap of samples (cols) vs features (rows)
#'
#' @inheritParams proBatch
#' @param factors_to_plot vector of technical and biological factors to be
#' plotted in this diagnostic plot (assumed to be present in
#' \code{sample_annotation})
#' @param fill_the_missing numeric value that the missing values are
#'   substituted with, or \code{NULL} if features with missing values are to be
#'   excluded.
#' @param color_for_missing special color to make missing values.
#' Usually black or white, depending on \code{heatmap_color}
#' @param heatmap_color vector of colors used in heatmap (typicall a gradient)
#' @param cluster_rows boolean value determining if rows
#' should be clustered
#' @param cluster_cols boolean value determining if columns
#' should be clustered
#' @param factors_of_feature_ann vector of factors that characterize features,
#' as listed in \code{peptide_annotation}
#' @param color_list_features list, as returned by
#' \code{sample_annotation_to_colors},
#' but mapping \code{peptide_annotation} where each item contains a color vector
#' for each factor to be mapped to the color.
#' @param ... other parameters of \code{link[pheatmap]{pheatmap}}
#'
#' @return object returned by \code{link[pheatmap]{pheatmap}}
#' @export
#'
#' @examples
#' # Load necessary datasets
#' data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#'
#' log_transformed_matrix <- log_transform_dm(example_proteome_matrix)
#' heatmap_plot <- plot_heatmap_diagnostic(log_transformed_matrix,
#'     example_sample_annotation,
#'     factors_to_plot = c("MS_batch", "digestion_batch", "Diet", "DateTime"),
#'     cluster_cols = TRUE, cluster_rows = FALSE,
#'     show_rownames = FALSE, show_colnames = FALSE
#' )
#'
#' color_list <- sample_annotation_to_colors(example_sample_annotation,
#'     factor_columns = c(
#'         "MS_batch", "EarTag", "Strain",
#'         "Diet", "digestion_batch", "Sex"
#'     ),
#'     numeric_columns = c("DateTime", "order")
#' )
#'
#' log_transformed_matrix <- log_transform_dm(example_proteome_matrix)
#' heatmap_plot <- plot_heatmap_diagnostic(log_transformed_matrix,
#'     example_sample_annotation,
#'     factors_to_plot = c("MS_batch", "digestion_batch", "Diet", "DateTime"),
#'     cluster_cols = TRUE, cluster_rows = FALSE,
#'     color_list = color_list,
#'     show_rownames = FALSE, show_colnames = FALSE
#' )
#'
#' @seealso \code{\link{sample_annotation_to_colors}},
#' \code{\link[pheatmap]{pheatmap}}
#'
#' @name plot_heatmap_diagnostic
plot_heatmap_diagnostic.default <- function(data_matrix, sample_annotation = NULL,
                                            sample_id_col = "FullRunName",
                                            factors_to_plot = NULL,
                                            fill_the_missing = -1,
                                            color_for_missing = "black",
                                            heatmap_color = colorRampPalette(
                                                rev(brewer.pal(
                                                    n = 7,
                                                    name = "RdYlBu"
                                                ))
                                            )(100),
                                            cluster_rows = TRUE, cluster_cols = FALSE,
                                            color_list = NULL,
                                            peptide_annotation = NULL,
                                            feature_id_col = NULL,
                                            factors_of_feature_ann = NULL,
                                            color_list_features = NULL,
                                            filename = NULL, width = 7, height = 7,
                                            units = c("cm", "in", "mm"),
                                            plot_title = NULL,
                                            ...) {
    df_long <- matrix_to_long(data_matrix, sample_id_col = sample_id_col)
    df_long <- check_sample_consistency(
        sample_annotation, sample_id_col, df_long,
        merge = FALSE
    )
    data_matrix <- long_to_matrix(df_long, sample_id_col = sample_id_col)
    rm(df_long)

    # infer the color scheme for sample annotation (cols)
    if (is.null(color_list) && !is.null(sample_annotation)) {
        warning("color_list for samples (cols) not defined, inferring automatically.
            Numeric/factor columns are guessed, for more controlled color
            mapping use sample_annotation_to_colors()")
        color_list <- sample_annotation_to_colors(
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            factor_columns = factors_to_plot,
            numeric_columns = NULL,
            guess_factors = TRUE
        )
    }

    if (is.null(factors_of_feature_ann) && !is.null(peptide_annotation)) {
        # in case c("KEGG_pathway","evolutionary_distance") are present in the annotation, use them
        factors_of_feature_ann <- intersect(
            c("KEGG_pathway", "evolutionary_distance"),
            names(peptide_annotation)
        )
    }

    # infer the color scheme for feature annotation (rows)
    if (is.null(color_list_features) && !is.null(peptide_annotation)) {
        warning("color_list_features for features (rows) not defined, inferring
            automatically. Numeric/factor columns are guessed, for more
            controlled color mapping use sample_annotation_to_colors()")

        if (is.null(factors_of_feature_ann)) {
            factors_of_feature_ann <- names(peptide_annotation)[sapply(
                peptide_annotation,
                function(x) is.factor(x) || is.character(x)
            )]
        }
        if (is.null(feature_id_col)) {
            stop("feature_id_col must be specified when peptide_annotation is provided")
        }


        color_list_features <- sample_annotation_to_colors(peptide_annotation,
            sample_id_col = feature_id_col,
            factor_columns = factors_of_feature_ann,
            numeric_columns = NULL,
            guess_factors = TRUE
        )
    }

    p <- plot_heatmap_generic(
        data_matrix,
        column_annotation_df = sample_annotation,
        row_annotation_df = peptide_annotation,
        fill_the_missing = fill_the_missing,
        col_ann_id_col = sample_id_col,
        row_ann_id_col = feature_id_col,
        columns_for_cols = factors_to_plot,
        columns_for_rows = factors_of_feature_ann,
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        annotation_color_cols = color_list,
        annotation_color_rows = color_list_features,
        heatmap_color = heatmap_color,
        color_for_missing = color_for_missing,
        filename = filename, width = width, height = height,
        units = units,
        plot_title = plot_title,
        ...
    )
    return(p)
}

#' @rdname plot_heatmap_diagnostic
#' @method plot_heatmap_diagnostic ProBatchFeatures
#' @export
plot_heatmap_diagnostic.ProBatchFeatures <- function(x, pbf_name = NULL,
                                                     sample_annotation = NULL,
                                                     sample_id_col = "FullRunName",
                                                     peptide_annotation = NULL,
                                                     feature_id_col = "peptide_group_label",
                                                     plot_title = NULL,
                                                     return_gridExtra = FALSE,
                                                     plot_ncol = NULL,
                                                     ...) {
    object <- x
    assays <- .pb_assays_to_plot(object, pbf_name)
    dots <- list(...)

    filename_list <- NULL
    if ("filename" %in% names(dots)) {
        filename_list <- .pb_split_arg_by_assay(dots$filename, assays)
        dots$filename <- NULL
    }
    if (length(assays) > 1L && !"silent" %in% names(dots)) {
        dots$silent <- TRUE
    }

    default_sample_annotation <- as.data.frame(colData(object))
    sample_ann_list <- .pb_split_arg_by_assay(sample_annotation, assays)
    peptide_ann_list <- .pb_split_arg_by_assay(peptide_annotation, assays)
    titles <- .pb_resolve_titles(assays, plot_title, default_fun = function(x) x)

    plot_list <- vector("list", length(assays))
    names(plot_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        data_matrix <- pb_assay_matrix(object, assay_nm)
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- default_sample_annotation
        }
        peptide_ann <- peptide_ann_list[[i]]
        if (is.null(peptide_ann) && assay_nm %in% names(object)) {
            peptide_ann <- as.data.frame(rowData(object[[assay_nm]]))
            peptide_ann[[feature_id_col]] <- rownames(peptide_ann)
        }

        call_args <- dots
        if (!is.null(filename_list)) {
            fn <- filename_list[[i]]
            if (!is.null(fn)) {
                call_args$filename <- fn
            }
        }

        call_args <- c(list(
            data_matrix = data_matrix,
            sample_annotation = sample_ann,
            sample_id_col = sample_id_col,
            peptide_annotation = peptide_ann,
            feature_id_col = feature_id_col,
            plot_title = titles[i]
        ), call_args)

        plot_list[[i]] <- do.call(plot_heatmap_diagnostic.default, call_args)
    }

    .pb_arrange_plot_list(plot_list, convert_fun = function(x) x$gtable, plot_ncol = plot_ncol, return_gridExtra = return_gridExtra)
}

#' @export
plot_heatmap_diagnostic <- function(x, ...) UseMethod("plot_heatmap_diagnostic")

#' Plot the heatmap
#'
#' @inheritParams proBatch
#' @param data_matrix the matrix of data to be plotted
#' @param column_annotation_df data frame annotating columns of
#' \code{data_matrix}
#' @param row_annotation_df data frame annotating rows of \code{data_matrix}
#' @param col_ann_id_col column of \code{column_annotation_df} whose values are
#' unique identifiers of columns in \code{data_matrix}
#' @param row_ann_id_col column of \code{row_annotation_df} whose values are
#' unique identifiers of rows in \code{data_matrix}
#' @param columns_for_cols vector of factors (columns) of
#' \code{column_annotation_df}
#' that will be mapped to color annotation of heatmap columns
#' @param columns_for_rows vector of factors (columns) of
#' \code{row_annotation_df}
#' that will be mapped to color annotation of heatmap rows
#' @param cluster_rows boolean: whether the rows should be clustered
#' @param cluster_cols boolean: whether the rows should be clustered
#' @param annotation_color_cols list of color vectors for column annotation,
#' for each factor to be plotted; for factor-like variables a named vector
#' (names should correspond to the levels of factors). Advisable to supply here
#' color list returned by \code{sample_annotation_to_colors}
#' @param annotation_color_rows list of color vectors for row annotation,
#' for each factor to be plotted; for factor-like variables a named vector
#' (names should correspond to the levels of factors). Advisable to supply here
#' color list returned by \code{sample_annotation_to_colors}
#' @param fill_the_missing numeric value that the missing values are
#'   substituted with, or \code{NULL} if features with missing values are to be
#'   excluded.
#' @param color_for_missing special color to make missing values.
#' Usually black or white, depending on \code{heatmap_color}
#' @param heatmap_color vector of colors used in heatmap (typicall a gradient)
#' @param ... other parameters of \code{link[pheatmap]{pheatmap}}
#'
#' @name plot_heatmap_generic
#'
#' @return pheatmap-type object
#' @export
#'
#' @examples
#' data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#' p <- plot_heatmap_generic(log_transform_dm(example_proteome_matrix),
#'     column_annotation_df = example_sample_annotation,
#'     columns_for_cols = c("MS_batch", "digestion_batch", "Diet", "DateTime"),
#'     plot_title = "test_heatmap",
#'     show_rownames = FALSE, show_colnames = FALSE
#' )
#'
plot_heatmap_generic.default <- function(data_matrix,
                                         column_annotation_df = NULL,
                                         row_annotation_df = NULL,
                                         col_ann_id_col = NULL,
                                         row_ann_id_col = NULL,
                                         columns_for_cols = c(
                                             "MS_batch", "Diet",
                                             "DateTime", "order"
                                         ),
                                         columns_for_rows = c(
                                             "KEGG_pathway",
                                             "WGCNA_module",
                                             "evolutionary_distance"
                                         ),
                                         cluster_rows = FALSE, cluster_cols = TRUE,
                                         annotation_color_cols = NULL,
                                         annotation_color_rows = NULL,
                                         fill_the_missing = -1,
                                         color_for_missing = "black",
                                         heatmap_color = colorRampPalette(
                                             rev(brewer.pal(
                                                 n = 7,
                                                 name = "RdYlBu"
                                             ))
                                         )(100),
                                         filename = NULL, width = 7, height = 7,
                                         units = c("cm", "in", "mm"),
                                         plot_title = NULL,
                                         ...) {
    # deal with the missing values
    warning_message <- "Heatmap cannot operate with missing values in the matrix"
    data_matrix <- handle_missing_values(
        data_matrix, warning_message,
        fill_the_missing
    )
    if (is.null(fill_the_missing) & any(is.na(data_matrix)) &
        (!cluster_rows | !cluster_cols)) {
        message("With NAs removed, clustering of heatmap will work,
              specify: cluster_rows = T, cluster_cols = T")
    }

    if (!is.null(fill_the_missing)) {
        heatmap_color <- c(color_for_missing, heatmap_color)
    }

    if (is.null(col_ann_id_col)) {
        col_ann_id_col <- "FullRunName"
        message(sprintf("Column %s is not in the data, using default", col_ann_id_col))
    }
    if (is.null(row_ann_id_col)) {
        row_ann_id_col <- "peptide_group_label"
        message(sprintf("Column %s is not in the data, using default", row_ann_id_col))
    }

    # if columns_for_cols is NULL, add default columns
    if (is.null(columns_for_cols)) {
        columns_for_cols <- intersect(
            c("MS_batch", "Diet", "DateTime", "order"),
            names(column_annotation_df)
        )
    }
    # if columns_for_rows is NULL, add default columns
    if (is.null(columns_for_rows)) {
        message("columns_for_rows is NULL, adding default columns if present")
        columns_for_rows <- intersect(
            c("KEGG_pathway", "WGCNA_module", "evolutionary_distance"),
            names(row_annotation_df)
        )
        if (length(columns_for_rows) < 1L) {
            # warning("No default columns for row annotation found in row_annotation_df")
            row_annotation_df <- NULL
            columns_for_rows <- NULL
            annotation_color_rows <- NULL
        }
    }

    annotation_col <- NA
    annotation_row <- NA
    if (!is.null(column_annotation_df)) {
        if (!is.null(columns_for_cols)) {
            annotation_col <- column_annotation_df %>%
                select(all_of(c(col_ann_id_col, columns_for_cols)))
        } else {
            annotation_col <- column_annotation_df
        }
        annotation_col <- annotation_col %>%
            mutate_if(is.POSIXct, as.numeric) %>%
            remove_rownames() %>%
            column_to_rownames(var = col_ann_id_col)
        if (is.data.frame(annotation_col) && ncol(annotation_col) == 0) {
            annotation_col <- NULL
        }
    }

    if (!is.null(row_annotation_df)) {
        if (!is.null(columns_for_rows)) {
            annotation_row <- row_annotation_df %>%
                select(all_of(c(row_ann_id_col, columns_for_rows)))
        } else {
            annotation_row <- row_annotation_df
        }

        annotation_row <- annotation_row %>%
            mutate_if(is.POSIXct, as.numeric) %>%
            remove_rownames() %>%
            column_to_rownames(var = row_ann_id_col)
        if (is.data.frame(annotation_row) && ncol(annotation_row) == 0) {
            annotation_row <- NULL
        }
    }

    if (is.null(column_annotation_df) || is.null(row_annotation_df)) {
        warning("annotation_row and / or annotation_col are not specified for heatmap
            (annotation of rows/cols such as sample annotation will not be plotted)")
    }

    if (!identical(annotation_col, NA) & (is.data.frame(annotation_col) | is.matrix(annotation_col))) {
        if (!setequal(rownames(annotation_col), colnames(data_matrix))) {
            warning("coloring by column annotation will not work: annotation rownames do not match data matrix column names")
        }
    }

    if (!identical(annotation_row, NA) & (is.data.frame(annotation_row) | is.matrix(annotation_row))) {
        if (!setequal(rownames(annotation_row), rownames(data_matrix))) {
            warning("coloring by row annotation will not work: annotation rownames do not match data matrix column names")
        }
    }

    if (is.null(plot_title)) {
        plot_title <- NA
    }
    if (is.null(filename)) {
        filename <- NA
    }
    units_adjusted <- adjust_units(units, width, height)
    units <- units_adjusted$unit
    width <- units_adjusted$width
    height <- units_adjusted$height

    if (is.list(annotation_color_cols) && !is.null(annotation_col)) {
        keep_cols <- intersect(names(annotation_color_cols), colnames(annotation_col))
        annotation_color_cols <- annotation_color_cols[keep_cols]
    } else {
        annotation_color_cols <- list()
    }
    if (is.list(annotation_color_rows) && !is.null(annotation_row)) {
        keep_rows <- intersect(names(annotation_color_rows), colnames(annotation_row))
        annotation_color_rows <- annotation_color_rows[keep_rows]
    } else {
        annotation_color_rows <- list()
    }

    annotation_color_list <- c(annotation_color_cols, annotation_color_rows)
    if (!length(annotation_color_list)) {
        annotation_color_list <- NA
    }

    p <- pheatmap(
        data_matrix,
        cluster_rows = cluster_rows, cluster_cols = cluster_cols,
        color = heatmap_color,
        annotation_col = annotation_col,
        annotation_row = annotation_row,
        annotation_colors = annotation_color_list,
        filename = filename, width = width, height = height,
        main = plot_title, ...
    )
    return(p)
}

#' @rdname plot_heatmap_generic
#' @method plot_heatmap_generic ProBatchFeatures
#' @export
plot_heatmap_generic.ProBatchFeatures <- function(x, pbf_name = NULL,
                                                  column_annotation_df = NULL,
                                                  row_annotation_df = NULL,
                                                  col_ann_id_col = NULL,
                                                  row_ann_id_col = NULL,
                                                  plot_title = NULL,
                                                  return_gridExtra = FALSE,
                                                  plot_ncol = NULL,
                                                  ...) {
    object <- x
    assays <- .pb_assays_to_plot(object, pbf_name)
    dots <- list(...)

    filename_list <- NULL
    if ("filename" %in% names(dots)) {
        filename_list <- .pb_split_arg_by_assay(dots$filename, assays)
        dots$filename <- NULL
    }
    if (length(assays) > 1L && !"silent" %in% names(dots)) {
        dots$silent <- TRUE
    }

    default_col_ann <- as.data.frame(colData(object))
    col_ann_list <- .pb_split_arg_by_assay(column_annotation_df, assays)
    row_ann_list <- .pb_split_arg_by_assay(row_annotation_df, assays)
    titles <- .pb_resolve_titles(assays, plot_title, default_fun = function(x) x)

    plot_list <- vector("list", length(assays))
    names(plot_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        data_matrix <- pb_assay_matrix(object, assay_nm)
        col_ann <- col_ann_list[[i]]
        if (is.null(col_ann)) {
            col_ann <- default_col_ann
        }
        row_ann <- row_ann_list[[i]]

        if (is.null(row_ann) && assay_nm %in% names(object) && (!is.null(rowData(object[[assay_nm]])))) {
            row_ann <- as.data.frame(rowData(object[[assay_nm]]))
            row_ann[[row_ann_id_col]] <- rownames(row_ann)
        }

        call_args <- dots
        if (!is.null(filename_list)) {
            fn <- filename_list[[i]]
            if (!is.null(fn)) {
                call_args$filename <- fn
            }
        }

        call_args <- c(list(
            data_matrix = data_matrix,
            column_annotation_df = col_ann,
            row_annotation_df = row_ann,
            col_ann_id_col = col_ann_id_col,
            row_ann_id_col = row_ann_id_col,
            plot_title = titles[i]
        ), call_args)
        plot_list[[i]] <- do.call(plot_heatmap_generic.default, call_args)
    }

    .pb_arrange_plot_list(plot_list, convert_fun = function(x) x$gtable, plot_ncol = plot_ncol, return_gridExtra = return_gridExtra)
}

#' @export
plot_heatmap_generic <- function(x, ...) UseMethod("plot_heatmap_generic")

#' Calculate variance distribution by variable
#'
#' @inheritParams proBatch
#' @param factors_for_PVCA vector of factors from \code{sample_annotation}, that
#'   are used in PVCA analysis
#' @param pca_threshold the percentile value of the minimum amount of the
#'   variabilities that the selected principal components need to explain
#' @param variance_threshold the percentile value of weight each of the factors
#'   needs to explain (the rest will be lumped together)
#' @param fill_the_missing numeric value determining how  missing values
#' should be substituted. If \code{NULL}, features with missing values are
#' excluded.
#'
#' @name calculate_PVCA
#' @return data frame of weights of Principal Variance Components
#' @export
#'
#' @examples
#' data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#' matrix_test <- example_proteome_matrix[1:150, ]
#' pvca_df <- calculate_PVCA(matrix_test, example_sample_annotation,
#'     factors_for_PVCA = c("MS_batch", "digestion_batch", "Diet", "Sex", "Strain"),
#'     pca_threshold = .6, variance_threshold = .01, fill_the_missing = -1
#' )
calculate_PVCA.default <- function(data_matrix, sample_annotation,
                                   feature_id_col = "peptide_group_label",
                                   sample_id_col = "FullRunName",
                                   factors_for_PVCA = c(
                                       "MS_batch", "digestion_batch",
                                       "Diet", "Sex", "Strain"
                                   ),
                                   pca_threshold = .6, variance_threshold = .01,
                                   fill_the_missing = -1) {
    df_long <- matrix_to_long(data_matrix, sample_id_col = sample_id_col)
    df_long <- check_sample_consistency(sample_annotation, sample_id_col, df_long,
        batch_col = NULL, order_col = NULL,
        facet_col = NULL, merge = FALSE
    )
    data_matrix <- long_to_matrix(df_long, sample_id_col = sample_id_col)

    # if factors_for_PVCA is NULL, add default columns
    if (is.null(factors_for_PVCA)) {
        factors_for_PVCA <- intersect(
            c("MS_batch", "digestion_batch", "Diet", "Sex", "Strain"),
            names(sample_annotation)
        )
    }

    sample_annotation <- sample_annotation %>%
        select(all_of(c(sample_id_col, factors_for_PVCA))) %>%
        mutate_if(is.POSIXct, as.numeric) %>%
        as.data.frame() %>%
        column_to_rownames(var = sample_id_col)

    data_matrix <- check_feature_id_col_in_dm(feature_id_col, data_matrix)

    warning_message <- "PVCA cannot operate with missing values in the matrix"
    data_matrix <- handle_missing_values(
        data_matrix, warning_message,
        fill_the_missing
    )

    covrts.annodf <- Biobase::AnnotatedDataFrame(data = sample_annotation)
    data_matrix <- data_matrix[, rownames(sample_annotation)]
    expr_set <- Biobase::ExpressionSet(
        assayData = data_matrix,
        phenoData = covrts.annodf
    )
    pvcaAssess <- pvcaBatchAssess(expr_set, factors_for_PVCA,
        threshold = pca_threshold
    )
    pvcaAssess_df <- data.frame(
        weights = as.vector(pvcaAssess$dat),
        label = pvcaAssess$label,
        stringsAsFactors = FALSE
    )

    label_of_small <- sprintf("Below %1.0f%%", 100 * variance_threshold)
    if (sum(pvcaAssess_df$weights < variance_threshold) > 1) {
        small_weights <- pvcaAssess_df$weights < variance_threshold
        pvca_res_small <- sum(pvcaAssess_df$weights[small_weights])
        pvca_res <- pvcaAssess_df[pvcaAssess_df$weights >= variance_threshold, ]
        pvca_res_add <- data.frame(weights = pvca_res_small, label = label_of_small)
        pvca_res <- rbind(pvca_res, pvca_res_add)
    } else {
        pvca_res <- pvcaAssess_df
    }
    return(pvca_res)
}

#' @rdname calculate_PVCA
#' @method calculate_PVCA ProBatchFeatures
#' @export
calculate_PVCA.ProBatchFeatures <- function(x, pbf_name = NULL,
                                            sample_annotation = NULL,
                                            feature_id_col = "peptide_group_label",
                                            sample_id_col = "FullRunName",
                                            ...) {
    object <- x
    data_matrix <- pb_assay_matrix(object, pbf_name)

    if (is.null(sample_annotation)) {
        sample_annotation <- as.data.frame(colData(object))
    }

    calculate_PVCA.default(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        ...
    )
}

#' @export
calculate_PVCA <- function(x, ...) UseMethod("calculate_PVCA")

#' Plot variance distribution by variable
#'
#' @inheritParams proBatch
#' @param technical_factors vector \code{sample_annotation} column names that
#' are technical covariates
#' @param biological_factors vector \code{sample_annotation} column names, that
#'   are biologically meaningful covariates
#' @param colors_for_bars four-item color vector, specifying colors for the
#'   following categories: c('residual', 'biological', 'biol:techn',
#'   'technical')
#' @param pca_threshold the percentile value of the minimum amount of the
#'   variabilities that the selected principal components need to explain
#' @param variance_threshold the percentile value of weight each of the
#'  covariates needs to explain (the rest will be lumped together)
#' @param fill_the_missing numeric value determining how  missing values
#' should be substituted. If \code{NULL}, features with missing values are
#' excluded.
#' If \code{NULL}, features with missing values are excluded.
#' @param base_size base size of the text in the plot
#'
#' @name plot_PVCA
#' @return \code{ggplot} object with the plot
#' @export
#'
#' @examples
#' data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#' matrix_test <- example_proteome_matrix[1:150, ]
#' pvca_plot <- plot_PVCA(matrix_test, example_sample_annotation,
#'     technical_factors = c("MS_batch", "digestion_batch"),
#'     biological_factors = c("Diet", "Sex", "Strain")
#' )
#'
#' \donttest{
#' pvca_plot <- plot_PVCA(matrix_test, example_sample_annotation,
#'     technical_factors = c("MS_batch", "digestion_batch"),
#'     biological_factors = c("Diet", "Sex", "Strain"),
#'     filename = "test_PVCA.png", width = 28, height = 22, units = "cm"
#' )
#' }
#'
#' @seealso \code{\link{sample_annotation_to_colors}},
#' \code{\link[ggplot2]{ggplot}}
plot_PVCA.default <- function(data_matrix, sample_annotation,
                              feature_id_col = "peptide_group_label",
                              sample_id_col = "FullRunName",
                              technical_factors = c("MS_batch", "instrument"),
                              biological_factors = c("cell_line", "drug_dose"),
                              fill_the_missing = -1,
                              pca_threshold = .6, variance_threshold = .01,
                              colors_for_bars = NULL,
                              filename = NULL, width = NA, height = NA,
                              units = c("cm", "in", "mm"),
                              plot_title = NULL,
                              theme = "classic",
                              base_size = 15) {
    pvca_res <- prepare_PVCA_df(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        technical_factors = technical_factors,
        biological_factors = biological_factors,
        fill_the_missing = fill_the_missing,
        pca_threshold = pca_threshold,
        variance_threshold = variance_threshold
    )

    gg <- plot_PVCA.df(
        pvca_res = pvca_res, colors_for_bars = colors_for_bars,
        filename = filename, width = width, height = height, units = units,
        plot_title = plot_title,
        theme = theme, base_size = base_size
    )
    return(gg)
}

#' @rdname plot_PVCA
#' @method plot_PVCA ProBatchFeatures
#' @export
plot_PVCA.ProBatchFeatures <- function(x, pbf_name = NULL,
                                       sample_annotation = NULL,
                                       feature_id_col = "peptide_group_label",
                                       sample_id_col = "FullRunName",
                                       plot_title = NULL,
                                       return_gridExtra = FALSE,
                                       plot_ncol = NULL,
                                       ...) {
    object <- x
    assays <- .pb_assays_to_plot(object, pbf_name)
    dots <- list(...)

    filename_list <- NULL
    if ("filename" %in% names(dots)) {
        filename_list <- .pb_split_arg_by_assay(dots$filename, assays)
        dots$filename <- NULL
    }

    if (is.null(sample_annotation)) {
        sample_annotation <- as.data.frame(colData(object))
        rownames(sample_annotation) <- NULL
    }
    sample_ann_list <- .pb_split_arg_by_assay(sample_annotation, assays)
    titles <- .pb_resolve_titles(assays, plot_title, default_fun = function(x) x)

    plot_list <- vector("list", length(assays))
    names(plot_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        data_matrix <- pb_assay_matrix(object, assay_nm)
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- as.data.frame(colData(object))
            rownames(sample_ann) <- NULL
        }

        call_args <- dots
        if (!is.null(filename_list)) {
            fn <- filename_list[[i]]
            if (!is.null(fn)) {
                call_args$filename <- fn
            }
        }

        call_args <- c(list(
            data_matrix = data_matrix,
            sample_annotation = sample_ann,
            feature_id_col = feature_id_col,
            sample_id_col = sample_id_col,
            plot_title = titles[i]
        ), call_args)

        plot_list[[i]] <- do.call(plot_PVCA.default, call_args)
    }

    .pb_arrange_plot_list(plot_list, convert_fun = ggplot2::ggplotGrob, plot_ncol = plot_ncol, return_gridExtra = return_gridExtra)
}

#' @export
plot_PVCA <- function(x, ...) UseMethod("plot_PVCA")

#' prepare the weights of Principal Variance Components
#'
#' @inheritParams proBatch
#' @param technical_factors vector \code{sample_annotation} column names that
#' are technical covariates
#' @param biological_factors vector \code{sample_annotation} column names, that
#'   are biologically meaningful covariates
#' @param pca_threshold the percentile value of the minimum amount of the
#'   variabilities that the selected principal components need to explain
#' @param variance_threshold the percentile value of weight each of the
#'  covariates needs to explain (the rest will be lumped together)
#' @param fill_the_missing numeric value determining how  missing values
#' should be substituted. If \code{NULL}, features with missing values are
#' excluded.
#' If \code{NULL}, features with missing values are excluded.
#'
#' @return data frame with weights and factors, combined in a way ready for plotting
#' @export
#'
#' @examples
#' data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#' matrix_test <- example_proteome_matrix[1:150, ]
#' pvca_df_res <- prepare_PVCA_df(matrix_test, example_sample_annotation,
#'     technical_factors = c("MS_batch", "digestion_batch"),
#'     biological_factors = c("Diet", "Sex", "Strain"),
#'     pca_threshold = .6, variance_threshold = .01, fill_the_missing = -1
#' )
#' @name prepare_PVCA_df
prepare_PVCA_df.default <- function(data_matrix, sample_annotation,
                                    feature_id_col = "peptide_group_label",
                                    sample_id_col = "FullRunName",
                                    technical_factors = c("MS_batch", "instrument"),
                                    biological_factors = c("cell_line", "drug_dose"),
                                    fill_the_missing = -1,
                                    pca_threshold = .6, variance_threshold = .01) {
    factors_for_PVCA <- c(technical_factors, biological_factors)

    pvca_res <- calculate_PVCA(
        data_matrix,
        sample_annotation,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        factors_for_PVCA = factors_for_PVCA,
        pca_threshold = pca_threshold,
        variance_threshold = variance_threshold,
        fill_the_missing = fill_the_missing
    )

    tech_interactions <- expand.grid(technical_factors, technical_factors) %>%
        mutate(tech_interactions = paste(Var1, Var2, sep = ":")) %>%
        pull(tech_interactions)
    biol_interactions <- expand.grid(biological_factors, biological_factors) %>%
        mutate(biol_interactions = paste(Var1, Var2, sep = ":")) %>%
        pull(biol_interactions)

    label_of_small <- sprintf("Below %1.0f%%", 100 * variance_threshold)
    technical_factors <- c(technical_factors, tech_interactions)
    biological_factors <- c(biological_factors, biol_interactions)
    pvca_res <- pvca_res %>% mutate(category = ifelse(label %in% technical_factors,
        "technical",
        ifelse(label %in% biological_factors,
            "biological",
            ifelse(label %in% c(label_of_small, "resid"),
                "residual", "biol:techn"
            )
        )
    ))

    pvca_res <- pvca_res %>%
        arrange(desc(weights)) %>%
        arrange(label == label_of_small) %>%
        arrange(label == "resid")
    return(pvca_res)
}

#' @rdname prepare_PVCA_df
#' @method prepare_PVCA_df ProBatchFeatures
#' @export
prepare_PVCA_df.ProBatchFeatures <- function(x, pbf_name = NULL,
                                             sample_annotation = NULL,
                                             feature_id_col = "peptide_group_label",
                                             sample_id_col = "FullRunName",
                                             ...) {
    object <- x
    data_matrix <- pb_assay_matrix(object, pbf_name)

    if (is.null(sample_annotation)) {
        sample_annotation <- as.data.frame(colData(object))
    }

    prepare_PVCA_df.default(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        ...
    )
}

#' @export
prepare_PVCA_df <- function(x, ...) UseMethod("prepare_PVCA_df")

#' plot PVCA, when the analysis is completed
#'
#' @inheritParams proBatch
#' @param pvca_res data frame of weights of Principal Variance Components, result
#' of \code{calculate_PVCA}
#' @param colors_for_bars four-item color vector, specifying colors for the
#'   following categories: c('residual', 'biological', 'biol:techn',
#'   'technical')
#' @param base_size base size of the text in the plot
#'
#' @return \code{ggplot} object with bars as weights, colored by bio/tech factors
#' @export
#'
#' @examples
#' data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#' matrix_test <- example_proteome_matrix[1:150, ]
#' pvca_df_res <- prepare_PVCA_df(matrix_test, example_sample_annotation,
#'     technical_factors = c("MS_batch", "digestion_batch"),
#'     biological_factors = c("Diet", "Sex", "Strain"),
#'     pca_threshold = .6, variance_threshold = .01, fill_the_missing = -1
#' )
#' colors_for_bars <- c("grey", "green", "blue", "red")
#' names(colors_for_bars) <- c("residual", "biological", "biol:techn", "technical")
#'
#' pvca_plot <- plot_PVCA.df(pvca_df_res, colors_for_bars)
#' @name plot_PVCA.df
plot_PVCA.df.default <- function(pvca_res,
                                 colors_for_bars = NULL,
                                 filename = NULL, width = NA, height = NA,
                                 units = c("cm", "in", "mm"),
                                 plot_title = NULL,
                                 theme = "classic",
                                 base_size = 15) {
    pvca_res <- pvca_res %>%
        mutate(label = factor(label, levels = label))

    gg <- ggplot(pvca_res, aes(x = label, y = weights, fill = category)) +
        geom_bar(stat = "identity", color = "black") +
        ylab("Weighted average proportion variance")

    if (is.null(colors_for_bars)) {
        colors_for_bars <- c("grey", wes_palettes$Rushmore[3:5])
        names(colors_for_bars) <- c(
            "residual", "biological",
            "biol:techn", "technical"
        )
    } else {
        if (length(colors_for_bars) != 4) {
            color_names <- paste(c(
                "residual", "biological", "biol:techn",
                "technical"
            ), collapse = " ")
            warning(sprintf("four colors for: %s were expected", color_names))
        }
    }
    gg <- gg + scale_fill_manual(values = colors_for_bars)

    if (!is.null(plot_title)) {
        gg <- gg + ggtitle(plot_title)
    }

    # Change the theme
    if (!is.null(theme) && theme == "classic") {
        gg <- gg + theme_classic(base_size = base_size)
    } else {
        message("plotting with default ggplot theme, only theme = 'classic'
            implemented")
    }

    gg <- gg +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
            plot.title = element_text(size = round(base_size * 1.2, 0))
        ) +
        xlab(NULL) +
        guides(fill = guide_legend(override.aes = list(color = NA), title = NULL))

    save_ggplot(filename, units, width, height, gg)

    return(gg)
}

#' @rdname plot_PVCA.df
#' @method plot_PVCA.df ProBatchFeatures
#' @export
plot_PVCA.df.ProBatchFeatures <- function(x, pbf_name = NULL,
                                          sample_annotation = NULL,
                                          feature_id_col = "peptide_group_label",
                                          sample_id_col = "FullRunName",
                                          colors_for_bars = NULL,
                                          filename = NULL, width = NA, height = NA,
                                          units = c("cm", "in", "mm"),
                                          plot_title = NULL,
                                          theme = "classic",
                                          base_size = 20,
                                          return_gridExtra = FALSE,
                                          plot_ncol = NULL,
                                          ...) {
    object <- x
    assays <- .pb_assays_to_plot(object, pbf_name)
    prepare_dots <- list(...)

    if (is.null(sample_annotation)) {
        sample_annotation <- as.data.frame(colData(object))
    }
    sample_ann_list <- .pb_split_arg_by_assay(sample_annotation, assays)
    filename_list <- .pb_split_arg_by_assay(filename, assays)
    titles <- .pb_resolve_titles(assays, plot_title, default_fun = function(x) x)

    plot_list <- vector("list", length(assays))
    names(plot_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        data_matrix <- pb_assay_matrix(object, assay_nm)
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- as.data.frame(colData(object))
        }

        prepare_args <- c(list(
            data_matrix = data_matrix,
            sample_annotation = sample_ann,
            feature_id_col = feature_id_col,
            sample_id_col = sample_id_col
        ), prepare_dots)
        pvca_res <- do.call(prepare_PVCA_df.default, prepare_args)

        plot_args <- list(
            pvca_res = pvca_res,
            colors_for_bars = colors_for_bars,
            filename = filename_list[[i]],
            width = width,
            height = height,
            units = units,
            plot_title = titles[i],
            theme = theme,
            base_size = base_size
        )

        plot_list[[i]] <- do.call(plot_PVCA.df.default, plot_args)
    }

    .pb_arrange_plot_list(plot_list, convert_fun = ggplot2::ggplotGrob, plot_ncol = plot_ncol, return_gridExtra = return_gridExtra)
}

#' @export
plot_PVCA.df <- function(x, ...) UseMethod("plot_PVCA.df")

#' plot PCA plot
#'
#' @inheritParams proBatch
#' @param color_by column name (as in \code{sample_annotation}) to color by
#' @param PC_to_plot principal component numbers for x and y axis
#' @param fill_the_missing numeric value determining how  missing values
#' should be substituted. If \code{NULL}, features with missing values are
#' excluded.
#' If \code{NULL}, features with missing values are excluded.
#' @param base_size base size of the text in the plot
#'
#' @return ggplot scatterplot colored by factor levels of column specified in
#'   \code{factor_to_color}
#' @name plot_PCA
#' @export
#'
#' @examples
#' data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#' pca_plot <- plot_PCA(example_proteome_matrix, example_sample_annotation,
#'     color_by = "MS_batch", plot_title = "PCA colored by MS batch"
#' )
#' pca_plot <- plot_PCA(example_proteome_matrix, example_sample_annotation,
#'     color_by = "DateTime", plot_title = "PCA colored by DateTime"
#' )
#'
#' color_list <- sample_annotation_to_colors(example_sample_annotation,
#'     factor_columns = c("MS_batch", "digestion_batch"),
#'     numeric_columns = c("DateTime", "order")
#' )
#' pca_plot <- plot_PCA(example_proteome_matrix, example_sample_annotation,
#'     color_by = "DateTime", color_scheme = color_list[["DateTime"]]
#' )
#'
#' \donttest{
#' pca_plot <- plot_PCA(example_proteome_matrix, example_sample_annotation,
#'     color_by = "DateTime", plot_title = "PCA colored by DateTime",
#'     filename = "test_PCA.png", width = 14, height = 9, units = "cm"
#' )
#' }
#'
#' @seealso \code{\link[ggfortify]{autoplot.pca_common}},
#' \code{\link[ggplot2]{ggplot}}
plot_PCA.default <- function(data_matrix, sample_annotation,
                             feature_id_col = "peptide_group_label",
                             sample_id_col = "FullRunName",
                             color_by = "MS_batch",
                             shape_by = NULL,
                             PC_to_plot = c(1, 2), fill_the_missing = -1,
                             color_scheme = "brewer",
                             filename = NULL, width = NA, height = NA,
                             units = c("cm", "in", "mm"),
                             plot_title = NULL,
                             theme = "classic",
                             base_size = 10, point_size = 3, point_alpha = 0.8) {
    df_long <- matrix_to_long(data_matrix, sample_id_col = sample_id_col)
    df_long <- check_sample_consistency(sample_annotation, sample_id_col, df_long,
        batch_col = color_by, order_col = NULL,
        facet_col = NULL, merge = FALSE
    )
    data_matrix <- long_to_matrix(df_long, sample_id_col = sample_id_col)
    data_matrix <- check_feature_id_col_in_dm(feature_id_col, data_matrix)

    # if any missing values, print a warning and handle them
    if (any(is.na(data_matrix))) {
        warning_message <- "PCA cannot operate with missing values in the matrix"
        # if fill the missing is FALSE, remove rows with NAs
        if (isFALSE(fill_the_missing)) {
            data_matrix <- data_matrix[complete.cases(data_matrix), ]
        }
        data_matrix <- handle_missing_values(
            data_matrix, warning_message,
            fill_the_missing
        )
    }

    if (!is.null(shape_by)) {
        if (length(shape_by) > 1) {
            warning("Shaping by the first column specified")
            shape_by <- shape_by[1]
        }
        if (!shape_by %in% colnames(sample_annotation)) {
            stop(sprintf(
                "Shaping column '%s' not found in sample_annotation",
                shape_by
            ))
        }

        shape_column <- sample_annotation[[shape_by]]
        if (!is.factor(shape_column) && !is.character(shape_column)) {
            sample_annotation[[shape_by]] <- as.factor(shape_column)
        }
    }

    pr_comp_res <- prcomp(t(data_matrix))
    gg <- autoplot(pr_comp_res,
        data = sample_annotation,
        x = PC_to_plot[1], y = PC_to_plot[2],
        size = point_size, alpha = point_alpha
    )

    if (!is.null(shape_by)) {
        gg <- gg + aes(shape = !!rlang::sym(shape_by)) +
            labs(shape = shape_by)
    }

    # add colors
    if (length(color_by) > 1) {
        warning("Coloring by the first column specified")
        color_by <- color_by[1]
    } # TODO: create the gridExtra graph with multiple panels, colored by factors

    if (color_by %in% names(color_scheme)) {
        color_scheme <- color_scheme[[color_by]]
    }

    gg <- color_by_factor(
        color_by_batch = TRUE,
        batch_col = color_by, gg = gg,
        color_scheme = color_scheme,
        sample_annotation = sample_annotation,
        fill_or_color = "color"
    )

    # Add the title
    if (!is.null(plot_title)) {
        gg <- gg + ggtitle(plot_title) +
            theme(plot.title = element_text(face = "bold", hjust = .5))
    }

    # Change the theme
    if (!is.null(theme) && theme == "classic") {
        gg <- gg + theme_classic(base_size = base_size)
    } else {
        message("plotting with default ggplot theme, only theme = 'classic'
            implemented")
    }

    save_ggplot(filename, units, width, height, gg)

    return(gg)
}

#' @rdname plot_PCA
#' @method plot_PCA ProBatchFeatures
#' @export
plot_PCA.ProBatchFeatures <- function(x, pbf_name = NULL,
                                      sample_annotation = NULL,
                                      sample_id_col = "FullRunName",
                                      plot_title = NULL,
                                      return_gridExtra = FALSE,
                                      plot_ncol = NULL,
                                      ...) {
    object <- x
    assays <- .pb_assays_to_plot(object, pbf_name)
    dots <- list(...)

    filename_list <- NULL
    if ("filename" %in% names(dots)) {
        filename_list <- .pb_split_arg_by_assay(dots$filename, assays)
        dots$filename <- NULL
    }

    if (is.null(sample_annotation)) {
        sample_annotation <- as.data.frame(colData(object))
    }
    sample_ann_list <- .pb_split_arg_by_assay(sample_annotation, assays)
    titles <- .pb_resolve_titles(assays, plot_title)

    plot_list <- vector("list", length(assays))
    names(plot_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        data_matrix <- pb_assay_matrix(object, assay_nm)
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- as.data.frame(colData(object))
        }

        call_args <- dots
        if (!is.null(filename_list)) {
            fn <- filename_list[[i]]
            if (!is.null(fn)) {
                call_args$filename <- fn
            }
        }

        call_args <- c(list(
            data_matrix = data_matrix,
            sample_annotation = sample_ann,
            sample_id_col = sample_id_col,
            plot_title = titles[i]
        ), call_args)

        plot_list[[i]] <- do.call(plot_PCA.default, call_args)
    }

    .pb_arrange_plot_list(plot_list, convert_fun = ggplot2::ggplotGrob, plot_ncol = plot_ncol, return_gridExtra = return_gridExtra)
}

#' @export
plot_PCA <- function(x, ...) UseMethod("plot_PCA")
