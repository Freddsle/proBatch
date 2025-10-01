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
    alignment <- .pb_align_matrix_and_annotation(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col
    )
    data_matrix <- alignment$data_matrix
    sample_annotation <- alignment$sample_annotation

    data_matrix <- .pb_handle_missing_wrapper(
        data_matrix = data_matrix,
        warning_message = "Hierarchical clustering cannot operate with missing values in the matrix",
        fill_the_missing = fill_the_missing
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

    color_list <- .pb_resolve_color_list(
        color_list = color_list,
        annotation_df = sample_annotation,
        id_col = sample_id_col,
        columns = factors_to_plot,
        warn_message = "color_list for samples annotation not defined, inferring automatically. Numeric/factor columns are guessed, for more controlled color mapping use sample_annotation_to_colors()"
    )

    color_df <- NULL
    if (!is.null(sample_annotation)) {
        color_source <- color_list
        if (is.null(color_source)) {
            color_source <- list()
        }
        color_df <- color_list_to_df(color_source, sample_annotation, sample_id_col)
    }

    device <- .pb_open_graphics_device(
        filename = filename,
        width = width,
        height = height,
        units = units,
        plot_title = plot_title
    )
    if (isTRUE(device$opened)) {
        on.exit(device$close(), add = TRUE)
    }

    plotDendroAndColors(hierarchical_clust, color_df,
        rowTextAlignment = "left",
        main = plot_title,
        hang = -0.1, addGuide = TRUE,
        dendroLabels = if (label_samples) NULL else FALSE,
        cex.dendroLabels = cex.dendroLabels,
        ...
    )
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
    prep <- .pb_prepare_multi_assay(
        object = object,
        pbf_name = pbf_name,
        dots = list(...),
        plot_title = plot_title,
        default_title_fun = function(x) x
    )
    assays <- prep$assays
    dots <- prep$dots
    filename_list <- prep$filename_list
    split_arg <- prep$split_arg
    titles <- prep$titles

    default_sample_annotation <- as.data.frame(colData(object))
    sample_ann_list <- split_arg(sample_annotation)

    plot_list <- vector("list", length(assays))
    names(plot_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        data_matrix <- pb_assay_matrix(object, assay_nm)
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- default_sample_annotation
        }

        call_args <- .pb_per_assay_dots(dots, filename_list, i)
        call_args <- c(list(
            data_matrix = data_matrix,
            sample_annotation = sample_ann,
            sample_id_col = sample_id_col,
            plot_title = titles[i]
        ), call_args)

        plot_list[[i]] <- do.call(plot_hierarchical_clustering.default, call_args)
    }

    if (length(plot_list) == 1L) {
        return(plot_list[[1L]])
    }

    invisible(plot_list)
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
    alignment <- .pb_align_matrix_and_annotation(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col
    )
    data_matrix <- alignment$data_matrix
    sample_annotation <- alignment$sample_annotation

    # infer the color scheme for sample annotation (cols)
    color_list <- .pb_resolve_color_list(
        color_list = color_list,
        annotation_df = sample_annotation,
        id_col = sample_id_col,
        columns = factors_to_plot,
        warn_message = "color_list for samples (cols) not defined, inferring automatically. Numeric/factor columns are guessed, for more controlled color mapping use sample_annotation_to_colors()"
    )

    if (is.null(factors_of_feature_ann) && !is.null(peptide_annotation)) {
        # in case c("KEGG_pathway","evolutionary_distance") are present in the annotation, use them
        factors_of_feature_ann <- intersect(
            c("KEGG_pathway", "evolutionary_distance"),
            names(peptide_annotation)
        )
    }

    # infer the color scheme for feature annotation (rows)
    if (!is.null(peptide_annotation)) {
        if (is.null(factors_of_feature_ann)) {
            factors_of_feature_ann <- names(peptide_annotation)[vapply(
                peptide_annotation,
                function(x) is.factor(x) || is.character(x),
                logical(1)
            )]
        }
        if (is.null(feature_id_col)) {
            stop("feature_id_col must be specified when peptide_annotation is provided")
        }
    }

    color_list_features <- .pb_resolve_color_list(
        color_list = color_list_features,
        annotation_df = peptide_annotation,
        id_col = feature_id_col,
        columns = factors_of_feature_ann,
        warn_message = "color_list_features for features (rows) not defined, inferring automatically. Numeric/factor columns are guessed, for more controlled color mapping use sample_annotation_to_colors()"
    )

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
    prep <- .pb_prepare_multi_assay(
        object = object,
        pbf_name = pbf_name,
        dots = list(...),
        plot_title = plot_title,
        default_title_fun = function(x) x,
        set_silent = TRUE
    )
    assays <- prep$assays
    dots <- prep$dots
    filename_list <- prep$filename_list
    split_arg <- prep$split_arg
    titles <- prep$titles

    default_sample_annotation <- as.data.frame(colData(object))
    sample_ann_list <- split_arg(sample_annotation)
    peptide_ann_list <- split_arg(peptide_annotation)

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

        call_args <- .pb_per_assay_dots(dots, filename_list, i)

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

    if (is.null(columns_for_rows)) {
        message("columns_for_rows is NULL, adding default columns if present")
    }

    ann_info <- .pb_prepare_pheatmap_annotations(
        data_matrix = data_matrix,
        column_annotation_df = column_annotation_df,
        row_annotation_df = row_annotation_df,
        col_ann_id_col = col_ann_id_col,
        row_ann_id_col = row_ann_id_col,
        columns_for_cols = columns_for_cols,
        columns_for_rows = columns_for_rows,
        annotation_color_cols = annotation_color_cols,
        annotation_color_rows = annotation_color_rows
    )
    annotation_col <- ann_info$annotation_col
    annotation_row <- ann_info$annotation_row
    annotation_color_list <- ann_info$annotation_color_list

    if (is.null(annotation_col) || is.null(annotation_row)) {
        warning("annotation_row and / or annotation_col are not specified for heatmap
            (annotation of rows/cols such as sample annotation will not be plotted)")
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
    prep <- .pb_prepare_multi_assay(
        object = object,
        pbf_name = pbf_name,
        dots = list(...),
        plot_title = plot_title,
        default_title_fun = function(x) x,
        set_silent = TRUE
    )
    assays <- prep$assays
    dots <- prep$dots
    filename_list <- prep$filename_list
    split_arg <- prep$split_arg
    titles <- prep$titles

    default_col_ann <- as.data.frame(colData(object))
    col_ann_list <- split_arg(column_annotation_df)
    row_ann_list <- split_arg(row_annotation_df)

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
            id_col <- if (is.null(row_ann_id_col)) "peptide_group_label" else row_ann_id_col
            row_ann[[id_col]] <- rownames(row_ann)
        }

        call_args <- .pb_per_assay_dots(dots, filename_list, i)

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

    # if factors_for_PVCA is NULL, add default columns
    if (is.null(factors_for_PVCA)) {
        factors_for_PVCA <- intersect(
            c("MS_batch", "digestion_batch", "Diet", "Sex", "Strain"),
            names(sample_annotation)
        )
    }

    sample_annotation <- sample_annotation %>%
        select(all_of(c(sample_id_col, factors_for_PVCA))) %>%
        mutate(across(where(is.POSIXct), as.numeric)) %>%
        as.data.frame() %>%
        remove_rownames() %>%
        column_to_rownames(var = sample_id_col)
    data_matrix <- check_feature_id_col_in_dm(feature_id_col, data_matrix)

    data_matrix <- .pb_handle_missing_wrapper(
        data_matrix = data_matrix,
        warning_message = "PVCA cannot operate with missing values in the matrix",
        fill_the_missing = fill_the_missing
    )

    covrts.annodf <- AnnotatedDataFrame(data = sample_annotation)
    data_matrix <- data_matrix[, rownames(sample_annotation)]
    expr_set <- ExpressionSet(
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
    prep <- .pb_prepare_multi_assay(
        object = object,
        pbf_name = pbf_name,
        dots = list(...),
        plot_title = NULL,
        default_title_fun = function(x) x
    )
    assays <- prep$assays
    dots <- prep$dots
    split_arg <- prep$split_arg

    default_sample_annotation <- as.data.frame(colData(object))
    sample_ann_list <- split_arg(sample_annotation)

    pvca_list <- vector("list", length(assays))
    names(pvca_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        data_matrix <- pb_assay_matrix(object, assay_nm)
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- default_sample_annotation
        }

        call_args <- c(list(
            data_matrix = data_matrix,
            sample_annotation = sample_ann,
            feature_id_col = feature_id_col,
            sample_id_col = sample_id_col
        ), dots)

        pvca_list[[i]] <- do.call(calculate_PVCA.default, call_args)
    }

    if (length(pvca_list) == 1L) {
        return(pvca_list[[1L]])
    }

    pvca_list
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
#' pvca_file <- tempfile("pvca", fileext = ".png")
#' pvca_plot <- plot_PVCA(matrix_test, example_sample_annotation,
#'     technical_factors = c("MS_batch", "digestion_batch"),
#'     biological_factors = c("Diet", "Sex", "Strain"),
#'     filename = pvca_file, width = 28, height = 22, units = "cm"
#' )
#' unlink(pvca_file)
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
    prep <- .pb_prepare_multi_assay(
        object = object,
        pbf_name = pbf_name,
        dots = list(...),
        plot_title = plot_title,
        default_title_fun = function(x) x
    )
    assays <- prep$assays
    dots <- prep$dots
    filename_list <- prep$filename_list
    split_arg <- prep$split_arg
    titles <- prep$titles

    if (is.null(sample_annotation)) {
        sample_annotation <- as.data.frame(colData(object))
        rownames(sample_annotation) <- NULL
    }
    sample_ann_list <- split_arg(sample_annotation)

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

        call_args <- .pb_per_assay_dots(dots, filename_list, i)

        call_args <- c(list(
            data_matrix = data_matrix,
            sample_annotation = sample_ann,
            feature_id_col = feature_id_col,
            sample_id_col = sample_id_col,
            plot_title = titles[i]
        ), call_args)

        plot_list[[i]] <- do.call(plot_PVCA.default, call_args)
    }

    .pb_arrange_plot_list(plot_list, convert_fun = ggplotGrob, plot_ncol = plot_ncol, return_gridExtra = return_gridExtra)
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
    prep <- .pb_prepare_multi_assay(
        object = object,
        pbf_name = pbf_name,
        dots = list(...),
        plot_title = NULL,
        default_title_fun = function(x) x
    )
    assays <- prep$assays
    dots <- prep$dots
    split_arg <- prep$split_arg

    default_sample_annotation <- as.data.frame(colData(object))
    sample_ann_list <- split_arg(sample_annotation)

    pvca_df_list <- vector("list", length(assays))
    names(pvca_df_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        data_matrix <- pb_assay_matrix(object, assay_nm)
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- default_sample_annotation
        }

        call_args <- c(list(
            data_matrix = data_matrix,
            sample_annotation = sample_ann,
            feature_id_col = feature_id_col,
            sample_id_col = sample_id_col
        ), dots)

        pvca_df_list[[i]] <- do.call(prepare_PVCA_df.default, call_args)
    }

    if (length(pvca_df_list) == 1L) {
        return(pvca_df_list[[1L]])
    }

    pvca_df_list
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
    prepare_dots <- list(...)
    prep <- .pb_prepare_multi_assay(
        object = object,
        pbf_name = pbf_name,
        dots = list(),
        plot_title = plot_title,
        default_title_fun = function(x) x
    )
    assays <- prep$assays
    titles <- prep$titles
    filename_list <- prep$split_arg(filename)

    if (is.null(sample_annotation)) {
        sample_annotation <- as.data.frame(colData(object))
    }
    sample_ann_list <- prep$split_arg(sample_annotation)

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

    .pb_arrange_plot_list(plot_list, convert_fun = ggplotGrob, plot_ncol = plot_ncol, return_gridExtra = return_gridExtra)
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
#' pca_file <- tempfile("pca_plot", fileext = ".png")
#' pca_plot <- plot_PCA(example_proteome_matrix, example_sample_annotation,
#'     color_by = "DateTime", plot_title = "PCA colored by DateTime",
#'     filename = pca_file, width = 14, height = 9, units = "cm"
#' )
#' unlink(pca_file)
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
    prep <- .pb_prepare_embedding_inputs(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        feature_id_col = feature_id_col,
        color_by = color_by,
        fill_the_missing = fill_the_missing,
        warning_message = "PCA cannot operate with missing values in the matrix"
    )
    data_matrix <- prep$data_matrix
    sample_annotation <- prep$sample_annotation

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
        gg <- gg + aes(shape = !!sym(shape_by)) +
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
    prep <- .pb_prepare_multi_assay(
        object = object,
        pbf_name = pbf_name,
        dots = list(...),
        plot_title = plot_title
    )
    assays <- prep$assays
    dots <- prep$dots
    filename_list <- prep$filename_list
    split_arg <- prep$split_arg
    titles <- prep$titles

    default_sample_annotation <- as.data.frame(colData(object))
    sample_ann_list <- split_arg(sample_annotation)

    plot_list <- vector("list", length(assays))
    names(plot_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        data_matrix <- pb_assay_matrix(object, assay_nm)
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- default_sample_annotation
        }

        call_args <- .pb_per_assay_dots(dots, filename_list, i)
        call_args <- c(list(
            data_matrix = data_matrix,
            sample_annotation = sample_ann,
            sample_id_col = sample_id_col,
            plot_title = titles[i]
        ), call_args)

        plot_list[[i]] <- do.call(plot_PCA.default, call_args)
    }

    .pb_arrange_plot_list(plot_list, convert_fun = ggplotGrob, plot_ncol = plot_ncol, return_gridExtra = return_gridExtra)
}

#' @export
plot_PCA <- function(x, ...) UseMethod("plot_PCA")

#' Plot t-SNE embedding of samples
#'
#' Generate a t-SNE visualization for a data matrix or `ProBatchFeatures`
#' object. A ggplot scatter plot is produced by default, with optional
#' Plotly-based rendering when requested.
#'
#' @inheritParams plot_PCA
#' @param perplexity positive numeric controlling the effective number of
#'   neighbours in t-SNE.
#' @param initial_dims number of principal components retained to initialise the
#'   t-SNE optimisation.
#' @param max_iter maximum number of optimisation iterations performed by t-SNE.
#' @param tsne_dims integer embedding dimensionality (first two dimensions are
#'   visualised).
#' @param random_seed optional integer passed to `set.seed()` for reproducible
#'   t-SNE initialisation.
#' @param return_gridExtra logical; when `TRUE` and plotting multiple assays,
#'   return a list containing the arranged grob along with the individual
#'   ggplot objects.
#' @param plot_ncol optional integer controlling the number of columns used when
#'   arranging multiple ggplots.
#' @param return_subplots logical; when `TRUE` and `use_plotlyrender = TRUE`,
#'   combine multiple assays into a single subplot layout produced by
#'   `plotly::subplot()`.
#' @param subplot_ncol optional integer specifying the number of subplot columns
#'   when `return_subplots = TRUE`.
#' @param share_axes logical indicating whether subplot axes should be shared in
#'   the Plotly output.
#' @param ... additional arguments forwarded to `Rtsne::Rtsne()` (for the default
#'   method) or to the respective default method when called on
#'   `ProBatchFeatures`.
#'
#' @return A `ggplot` object displaying the t-SNE embedding by default, or a
#'   `plotly` object when `use_plotlyrender = TRUE`.
#' @name plot_TSNE
#'
#' @examples
#' \dontrun{
#' data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#' tsne_plot <- plot_TSNE(example_proteome_matrix, example_sample_annotation,
#'     color_by = "MS_batch", perplexity = 5, max_iter = 500
#' )
#' }
#' @method plot_TSNE default
#' @export
plot_TSNE.default <- function(data_matrix, sample_annotation,
                              feature_id_col = "peptide_group_label",
                              sample_id_col = "FullRunName",
                              color_by = "MS_batch",
                              shape_by = NULL,
                              tsne_dims = 2,
                              perplexity = 30,
                              initial_dims = 50,
                              max_iter = 1000,
                              fill_the_missing = -1,
                              color_scheme = "brewer",
                              plot_title = NULL,
                              point_size = 8,
                              point_alpha = 0.85,
                              random_seed = NULL,
                              ...) {
    if (!requireNamespace("Rtsne", quietly = TRUE)) {
        stop("Package 'Rtsne' is required for plot_TSNE(); install it with install.packages('Rtsne').", call. = FALSE)
    }

    if (tsne_dims < 2) {
        stop("tsne_dims must be >= 2 to create a 2D scatter plot.")
    }

    dots <- list(...)
    use_plotlyrender <- isTRUE(dots$use_plotlyrender)
    dots$use_plotlyrender <- NULL

    sample_annotation <- as.data.frame(sample_annotation)

    prep <- .pb_prepare_embedding_inputs(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        feature_id_col = feature_id_col,
        color_by = color_by,
        fill_the_missing = fill_the_missing,
        warning_message = "t-SNE cannot operate with missing values in the matrix",
        allow_partial_annotation = FALSE,
        check_args = list(batch_col = color_by),
        drop_on_false = TRUE
    )
    data_matrix <- prep$data_matrix
    sample_annotation <- prep$sample_annotation
    sample_ids <- prep$sample_ids

    if (!is.null(shape_by)) {
        if (length(shape_by) > 1) {
            warning("Shaping by the first column specified")
            shape_by <- shape_by[1]
        }
        if (!shape_by %in% colnames(sample_annotation)) {
            stop(sprintf("Shaping column '%s' not found in sample_annotation", shape_by))
        }

        shape_column <- sample_annotation[[shape_by]]
        if (!is.factor(shape_column) && !is.character(shape_column)) {
            sample_annotation[[shape_by]] <- as.factor(shape_column)
        }
    }

    if (length(color_by) > 1) {
        warning("Coloring by the first column specified")
        color_by <- color_by[1]
    }

    n_samples <- length(sample_ids)
    if (n_samples < 2) {
        stop("At least two samples are required to compute t-SNE.")
    }

    max_perplexity <- max(1, floor((n_samples - 1) / 3))
    if (perplexity > max_perplexity) {
        warning(sprintf(
            "perplexity %.2f exceeds the recommended maximum %.2f for %d samples; adjusting to %.2f",
            perplexity, max_perplexity, n_samples, max_perplexity
        ))
        perplexity <- max_perplexity
    }

    if (!is.null(random_seed)) {
        set.seed(random_seed)
    }

    tsne_input <- t(as.matrix(data_matrix))
    tsne_args <- c(list(
        X = tsne_input,
        dims = tsne_dims,
        initial_dims = initial_dims,
        perplexity = perplexity,
        max_iter = max_iter,
        check_duplicates = FALSE
    ), dots)
    tsne_res <- do.call(Rtsne::Rtsne, tsne_args)

    tsne_matrix <- as.matrix(tsne_res$Y)
    if (ncol(tsne_matrix) < 2) {
        stop("The computed t-SNE embedding has fewer than two dimensions and cannot be plotted.")
    }

    axis_labels <- list(
        title = "t-SNE embedding",
        x = "t-SNE 1",
        y = "t-SNE 2"
    )

    if (isTRUE(use_plotlyrender)) {
        if (!requireNamespace("plotly", quietly = TRUE)) {
            stop("Package 'plotly' is required when use_plotlyrender = TRUE; install it with install.packages('plotly').", call. = FALSE)
        }

        plot <- .pb_create_embedding_plotly(
            embedding_matrix = tsne_matrix,
            sample_ids = sample_ids,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            color_by = color_by,
            shape_by = shape_by,
            color_scheme = color_scheme,
            point_size = point_size,
            point_alpha = point_alpha,
            plot_title = plot_title,
            axis_labels = axis_labels
        )

        return(plot)
    }

    plot <- .pb_create_embedding_ggplot(
        embedding_matrix = tsne_matrix,
        sample_ids = sample_ids,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        color_by = color_by,
        shape_by = shape_by,
        color_scheme = color_scheme,
        point_size = point_size,
        point_alpha = point_alpha,
        plot_title = plot_title,
        axis_labels = axis_labels
    )

    return(plot)
}

#' @rdname plot_TSNE
#' @method plot_TSNE ProBatchFeatures
#' @export
plot_TSNE.ProBatchFeatures <- function(x, pbf_name = NULL,
                                       sample_annotation = NULL,
                                       sample_id_col = "FullRunName",
                                       plot_title = NULL,
                                       return_gridExtra = FALSE,
                                       plot_ncol = NULL,
                                       return_subplots = FALSE,
                                       subplot_ncol = NULL,
                                       share_axes = TRUE,
                                       ...) {
    object <- x
    prep <- .pb_prepare_multi_assay(
        object = object,
        pbf_name = pbf_name,
        dots = list(...),
        plot_title = plot_title
    )
    assays <- prep$assays
    dots <- prep$dots
    use_plotlyrender <- isTRUE(dots$use_plotlyrender)
    filename_list <- prep$filename_list
    split_arg <- prep$split_arg
    titles <- prep$titles

    default_sample_annotation <- as.data.frame(colData(object))
    sample_ann_list <- split_arg(sample_annotation)

    plot_list <- vector("list", length(assays))
    names(plot_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        data_matrix <- pb_assay_matrix(object, assay_nm)
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- default_sample_annotation
        }

        call_args <- .pb_per_assay_dots(dots, filename_list, i)
        call_args <- c(list(
            data_matrix = data_matrix,
            sample_annotation = sample_ann,
            sample_id_col = sample_id_col,
            plot_title = titles[i]
        ), call_args)

        plot_list[[i]] <- do.call(plot_TSNE.default, call_args)
    }

    if (length(plot_list) == 1L) {
        return(plot_list[[1L]])
    }

    if (isTRUE(use_plotlyrender)) {
        if (isTRUE(return_gridExtra)) {
            warning("return_gridExtra is ignored when use_plotlyrender = TRUE; returning list of plotly objects instead.")
        }
        if (isTRUE(return_subplots)) {
            if (!requireNamespace("plotly", quietly = TRUE)) {
                stop("Package 'plotly' is required to build subplots; install it with install.packages('plotly').", call. = FALSE)
            }
            n_plots <- length(plot_list)
            ncol <- if (is.null(subplot_ncol)) ceiling(sqrt(n_plots)) else subplot_ncol
            nrow <- ceiling(n_plots / ncol)
            subplot_args <- c(plot_list, list(
                nrows = nrow,
                shareX = share_axes,
                shareY = share_axes,
                titleX = TRUE,
                titleY = TRUE
            ))
            return(do.call(plotly::subplot, subplot_args))
        }

        return(plot_list)
    }

    if (isTRUE(return_subplots)) {
        warning("return_subplots = TRUE is only supported when use_plotlyrender = TRUE; arranging ggplot outputs instead.")
    }

    .pb_arrange_plot_list(plot_list, convert_fun = ggplotGrob, plot_ncol = plot_ncol, return_gridExtra = return_gridExtra)
}

#' @rdname plot_TSNE
#' @export
plot_TSNE <- function(x, ...) UseMethod("plot_TSNE")

#' Plot UMAP embedding of samples
#'
#' Produce a UMAP visualization for a data matrix or `ProBatchFeatures`
#' object. A ggplot representation is returned by default, with optional
#' Plotly-based rendering when requested.
#'
#' @inheritParams plot_PCA
#' @param n_neighbors size of the local neighbourhood used by UMAP.
#' @param min_dist minimum distance between embedded points.
#' @param metric distance metric used by UMAP.
#' @param n_components dimensionality of the embedding (first two components are
#'   shown).
#' @param random_state optional integer passed to `config$random_state` for
#'   reproducibility.
#' @param spread optional numeric controlling how clustered points are in the
#'   embedding (forwarded to the UMAP configuration when supplied).
#' @param learning_rate optional numeric learning rate for the UMAP optimiser
#'   (forwarded to the configuration when supplied).
#' @param return_gridExtra logical; when `TRUE` and plotting multiple assays,
#'   return a list containing the arranged grob along with the individual
#'   ggplot objects.
#' @param plot_ncol optional integer controlling the number of columns used when
#'   arranging multiple ggplots.
#' @param return_subplots logical; when `TRUE` and `use_plotlyrender = TRUE`,
#'   combine multiple assays into a single subplot layout produced by
#'   `plotly::subplot()`.
#' @param subplot_ncol optional integer specifying the number of subplot columns
#'   when `return_subplots = TRUE`.
#' @param share_axes logical indicating whether subplot axes should be shared in
#'   the Plotly output.
#' @param ... additional arguments forwarded to `umap::umap()` (default method)
#'   or the respective default method when called on `ProBatchFeatures`.
#'
#' @return A `ggplot` object displaying the UMAP embedding by default, or a
#'   `plotly` object when `use_plotlyrender = TRUE`.
#' @name plot_UMAP
#'
#' @examples
#' \dontrun{
#' data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#' umap_plot <- plot_UMAP(example_proteome_matrix, example_sample_annotation,
#'     color_by = "MS_batch", n_neighbors = 10
#' )
#' }
#' @method plot_UMAP default
#' @export
plot_UMAP.default <- function(data_matrix, sample_annotation,
                              feature_id_col = "peptide_group_label",
                              sample_id_col = "FullRunName",
                              color_by = "MS_batch",
                              shape_by = NULL,
                              n_neighbors = 15,
                              min_dist = 0.1,
                              metric = "euclidean",
                              n_components = 2,
                              fill_the_missing = -1,
                              color_scheme = "brewer",
                              plot_title = NULL,
                              point_size = 8,
                              point_alpha = 0.85,
                              random_state = NULL,
                              spread = NULL,
                              learning_rate = NULL,
                              ...) {
    if (!requireNamespace("umap", quietly = TRUE)) {
        stop("Package 'umap' is required for plot_UMAP(); install it with install.packages('umap').", call. = FALSE)
    }

    if (n_components < 2) {
        stop("n_components must be >= 2 to create a 2D scatter plot.")
    }

    dots <- list(...)
    use_plotlyrender <- isTRUE(dots$use_plotlyrender)
    dots$use_plotlyrender <- NULL

    sample_annotation <- as.data.frame(sample_annotation)

    prep <- .pb_prepare_embedding_inputs(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        feature_id_col = feature_id_col,
        color_by = color_by,
        fill_the_missing = fill_the_missing,
        warning_message = "UMAP cannot operate with missing values in the matrix",
        allow_partial_annotation = FALSE,
        check_args = list(batch_col = color_by),
        drop_on_false = TRUE
    )
    data_matrix <- prep$data_matrix
    sample_annotation <- prep$sample_annotation
    sample_ids <- prep$sample_ids

    if (!is.null(shape_by)) {
        if (length(shape_by) > 1) {
            warning("Shaping by the first column specified")
            shape_by <- shape_by[1]
        }
        if (!shape_by %in% colnames(sample_annotation)) {
            stop(sprintf("Shaping column '%s' not found in sample_annotation", shape_by))
        }

        shape_column <- sample_annotation[[shape_by]]
        if (!is.factor(shape_column) && !is.character(shape_column)) {
            sample_annotation[[shape_by]] <- as.factor(shape_column)
        }
    }

    if (length(color_by) > 1) {
        warning("Coloring by the first column specified")
        color_by <- color_by[1]
    }

    config <- umap::umap.defaults
    config$n_neighbors <- n_neighbors
    config$min_dist <- min_dist
    config$metric <- metric
    config$n_components <- n_components
    if (!is.null(random_state)) {
        config$random_state <- random_state
    }
    if (!is.null(spread)) {
        config$spread <- spread
    }
    if (!is.null(learning_rate)) {
        config$learning_rate <- learning_rate
    }

    umap_input <- t(as.matrix(data_matrix))
    umap_args <- c(list(
        d = umap_input,
        config = config
    ), dots)
    umap_res <- do.call(umap::umap, umap_args)
    umap_matrix <- as.matrix(umap_res$layout)

    if (ncol(umap_matrix) < 2) {
        stop("The computed UMAP embedding has fewer than two dimensions and cannot be plotted.")
    }

    axis_labels <- list(
        title = "UMAP embedding",
        x = "UMAP 1",
        y = "UMAP 2"
    )

    if (isTRUE(use_plotlyrender)) {
        if (!requireNamespace("plotly", quietly = TRUE)) {
            stop("Package 'plotly' is required when use_plotlyrender = TRUE; install it with install.packages('plotly').", call. = FALSE)
        }

        plot <- .pb_create_embedding_plotly(
            embedding_matrix = umap_matrix,
            sample_ids = sample_ids,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            color_by = color_by,
            shape_by = shape_by,
            color_scheme = color_scheme,
            point_size = point_size,
            point_alpha = point_alpha,
            plot_title = plot_title,
            axis_labels = axis_labels
        )

        return(plot)
    }

    plot <- .pb_create_embedding_ggplot(
        embedding_matrix = umap_matrix,
        sample_ids = sample_ids,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        color_by = color_by,
        shape_by = shape_by,
        color_scheme = color_scheme,
        point_size = point_size,
        point_alpha = point_alpha,
        plot_title = plot_title,
        axis_labels = axis_labels
    )

    return(plot)
}

#' @rdname plot_UMAP
#' @method plot_UMAP ProBatchFeatures
#' @export
plot_UMAP.ProBatchFeatures <- function(x, pbf_name = NULL,
                                       sample_annotation = NULL,
                                       sample_id_col = "FullRunName",
                                       plot_title = NULL,
                                       return_gridExtra = FALSE,
                                       plot_ncol = NULL,
                                       return_subplots = FALSE,
                                       subplot_ncol = NULL,
                                       share_axes = TRUE,
                                       ...) {
    object <- x
    prep <- .pb_prepare_multi_assay(
        object = object,
        pbf_name = pbf_name,
        dots = list(...),
        plot_title = plot_title
    )
    assays <- prep$assays
    dots <- prep$dots
    use_plotlyrender <- isTRUE(dots$use_plotlyrender)
    filename_list <- prep$filename_list
    split_arg <- prep$split_arg
    titles <- prep$titles

    default_sample_annotation <- as.data.frame(colData(object))
    sample_ann_list <- split_arg(sample_annotation)

    plot_list <- vector("list", length(assays))
    names(plot_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        data_matrix <- pb_assay_matrix(object, assay_nm)
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- default_sample_annotation
        }

        call_args <- .pb_per_assay_dots(dots, filename_list, i)
        call_args <- c(list(
            data_matrix = data_matrix,
            sample_annotation = sample_ann,
            sample_id_col = sample_id_col,
            plot_title = titles[i]
        ), call_args)

        plot_list[[i]] <- do.call(plot_UMAP.default, call_args)
    }

    if (length(plot_list) == 1L) {
        return(plot_list[[1L]])
    }

    if (isTRUE(use_plotlyrender)) {
        if (isTRUE(return_gridExtra)) {
            warning("return_gridExtra is ignored when use_plotlyrender = TRUE; returning list of plotly objects instead.")
        }
        if (isTRUE(return_subplots)) {
            if (!requireNamespace("plotly", quietly = TRUE)) {
                stop("Package 'plotly' is required to build subplots; install it with install.packages('plotly').", call. = FALSE)
            }
            n_plots <- length(plot_list)
            ncol <- if (is.null(subplot_ncol)) ceiling(sqrt(n_plots)) else subplot_ncol
            nrow <- ceiling(n_plots / ncol)
            subplot_args <- c(plot_list, list(
                nrows = nrow,
                shareX = share_axes,
                shareY = share_axes,
                titleX = TRUE,
                titleY = TRUE
            ))
            return(do.call(plotly::subplot, subplot_args))
        }

        return(plot_list)
    }

    if (isTRUE(return_subplots)) {
        warning("return_subplots = TRUE is only supported when use_plotlyrender = TRUE; arranging ggplot outputs instead.")
    }

    .pb_arrange_plot_list(plot_list, convert_fun = ggplotGrob, plot_ncol = plot_ncol, return_gridExtra = return_gridExtra)
}

#' @rdname plot_UMAP
#' @export
plot_UMAP <- function(x, ...) UseMethod("plot_UMAP")

# Internal helpers for embedding plots ---------------------------------------

.pb_create_embedding_ggplot <- function(embedding_matrix, sample_ids, sample_annotation,
                                        sample_id_col, color_by, shape_by,
                                        color_scheme, point_size, point_alpha,
                                        plot_title, axis_labels) {
    plot_df <- data.frame(
        sample_id = sample_ids,
        Dim1 = embedding_matrix[, 1],
        Dim2 = embedding_matrix[, 2],
        stringsAsFactors = FALSE
    )

    plot_df <- cbind(plot_df, sample_annotation)

    point_aes <- aes(color = !!sym(color_by))
    if (!is.null(shape_by)) {
        point_aes <- aes(color = !!sym(color_by), shape = !!sym(shape_by))
    }

    gg <- ggplot(plot_df, aes(x = Dim1, y = Dim2)) +
        geom_point(mapping = point_aes, size = point_size, alpha = point_alpha) +
        labs(x = axis_labels$x, y = axis_labels$y, color = color_by)

    if (!is.null(shape_by)) {
        gg <- gg + labs(shape = shape_by)
    }

    title_to_use <- if (is.null(plot_title)) axis_labels$title else plot_title
    if (!is.null(title_to_use) && nzchar(title_to_use)) {
        gg <- gg + ggtitle(title_to_use) +
            theme(plot.title = element_text(face = "bold", hjust = 0.5))
    }

    scheme_to_use <- color_scheme
    if (is.list(scheme_to_use) && !is.null(names(scheme_to_use)) && color_by %in% names(scheme_to_use)) {
        scheme_to_use <- scheme_to_use[[color_by]]
    }

    gg <- color_by_factor(
        color_by_batch = TRUE,
        batch_col = color_by,
        gg = gg,
        color_scheme = scheme_to_use,
        sample_annotation = sample_annotation,
        fill_or_color = "color"
    )

    gg
}

.pb_prepare_annotation_for_samples <- function(sample_annotation, sample_id_col, sample_ids, allow_partial = FALSE) {
    if (is.null(sample_annotation)) {
        stop("sample_annotation must be provided to colour or shape embeddings.")
    }
    if (!sample_id_col %in% names(sample_annotation)) {
        stop(sprintf("Sample ID column '%s' not found in sample_annotation.", sample_id_col))
    }

    sample_annotation <- sample_annotation[!duplicated(sample_annotation[[sample_id_col]]), , drop = FALSE]
    match_idx <- match(sample_ids, sample_annotation[[sample_id_col]])
    if (any(is.na(match_idx))) {
        missing_ids <- sample_ids[is.na(match_idx)]
        if (!isTRUE(allow_partial)) {
            stop(sprintf(
                "Sample annotation is missing entries for the following samples: %s",
                paste(missing_ids, collapse = ", ")
            ))
        }
        warning(sprintf(
            "Sample annotation is missing entries for the following samples: %s",
            paste(missing_ids, collapse = ", ")
        ))
        keep <- !is.na(match_idx)
        match_idx <- match_idx[keep]
        sample_annotation <- sample_annotation[match_idx, , drop = FALSE]
        return(sample_annotation)
    }

    sample_annotation[match_idx, , drop = FALSE]
}

.pb_create_embedding_plotly <- function(embedding_matrix, sample_ids, sample_annotation,
                                        sample_id_col, color_by, shape_by,
                                        color_scheme, point_size, point_alpha,
                                        plot_title, axis_labels) {
    plot_df <- data.frame(
        sample_id = sample_ids,
        Dim1 = embedding_matrix[, 1],
        Dim2 = embedding_matrix[, 2],
        stringsAsFactors = FALSE
    )

    plot_df <- cbind(plot_df, sample_annotation)

    color_info <- .pb_resolve_plotly_color_mapping(plot_df, color_by, color_scheme)
    plot_df$.color_value <- color_info$aes_column

    shape_info <- .pb_resolve_plotly_symbol_mapping(plot_df, shape_by)
    if (!is.null(shape_info)) {
        plot_df$.shape_value <- shape_info$aes_column
    }

    hover_columns <- c(sample_id_col, color_info$legend_title)
    if (!is.null(shape_by)) {
        hover_columns <- c(hover_columns, shape_by)
    }
    hover_columns <- unique(hover_columns)
    plot_df$.hover_text <- .pb_build_embedding_hover_text(plot_df, hover_columns)

    marker <- list(size = point_size, opacity = point_alpha)
    plot_args <- list(
        data = plot_df,
        x = ~Dim1,
        y = ~Dim2,
        type = "scatter",
        mode = "markers",
        text = ~.hover_text,
        hoverinfo = "text",
        color = ~.color_value,
        colors = color_info$palette
    )

    if (!is.null(shape_info)) {
        plot_args$symbol <- ~.shape_value
        plot_args$symbols <- shape_info$symbols
    }

    if (color_info$type == "numeric") {
        marker$colorbar <- list(title = color_info$legend_title)
    }

    plot_args$marker <- marker

    plt <- do.call(plotly::plot_ly, plot_args)

    layout_args <- list(
        title = list(text = if (is.null(plot_title)) axis_labels$title else plot_title),
        xaxis = list(title = axis_labels$x),
        yaxis = list(title = axis_labels$y)
    )
    if (color_info$type == "discrete") {
        layout_args$legend <- list(title = list(text = color_info$legend_title))
    }

    plt <- do.call(plotly::layout, c(list(p = plt), layout_args))

    return(plt)
}

.pb_resolve_plotly_color_mapping <- function(plot_df, color_by, color_scheme) {
    if (is.null(color_by)) {
        stop("Coloring column not defined, please define the color column!")
    }
    if (length(color_by) > 1) {
        warning("Coloring by the first column specified")
        color_by <- color_by[1]
    }
    if (!color_by %in% names(plot_df)) {
        stop(sprintf("Coloring column '%s' not found in the data used for plotting.", color_by))
    }

    column <- plot_df[[color_by]]

    palette <- color_scheme
    if (is.list(color_scheme) && !is.null(names(color_scheme)) && color_by %in% names(color_scheme)) {
        palette <- color_scheme[[color_by]]
    }

    is_factor_col <- is_batch_factor(column, if (is.list(color_scheme) && color_by %in% names(color_scheme)) palette else NULL)
    is_numeric_col <- (!is_factor_col) && (is.numeric(column) || inherits(column, "POSIXct") || inherits(column, "POSIXt") || inherits(column, "Date"))

    if (is_numeric_col) {
        numeric_values <- as.numeric(column)
        if (length(palette) == 1 && identical(palette, "brewer")) {
            palette <- RColorBrewer::brewer.pal(11, "PiYG")
        }
        if (is.null(palette)) {
            palette <- RColorBrewer::brewer.pal(11, "PiYG")
        }
        if (length(palette) < 2) {
            palette <- grDevices::hcl.colors(11, "Inferno")
        }
        return(list(
            aes_column = numeric_values,
            palette = palette,
            legend_title = color_by,
            type = "numeric"
        ))
    }

    factor_values <- as.factor(column)
    n_levels <- length(levels(factor_values))
    if (n_levels == 0) {
        stop(sprintf("Coloring column '%s' has no values after factoring.", color_by))
    }

    if (length(palette) == 1 && identical(palette, "brewer")) {
        if (n_levels <= 9) {
            palette <- RColorBrewer::brewer.pal(max(3, n_levels), "Set1")[seq_len(n_levels)]
        } else if (n_levels <= 12) {
            palette <- RColorBrewer::brewer.pal(n_levels, "Set3")[seq_len(n_levels)]
        } else {
            warning("brewer palettes have maximally 12 colors, generating palette with grDevices::hcl.colors")
            palette <- grDevices::hcl.colors(n_levels, "Set 3")
        }
    } else if (is.null(palette)) {
        palette <- grDevices::hcl.colors(n_levels, "Set 3")
    }

    if (!is.null(names(palette))) {
        palette <- palette[levels(factor_values)]
    }
    if (any(is.na(palette))) {
        fallback <- grDevices::hcl.colors(n_levels, "Set 3")
        palette[is.na(palette)] <- fallback[is.na(palette)]
    }
    if (length(palette) < n_levels) {
        warning("Color scheme has fewer entries than factor levels; colors will be recycled.")
        palette <- rep(palette, length.out = n_levels)
    }

    names(palette) <- levels(factor_values)

    list(
        aes_column = factor_values,
        palette = palette,
        legend_title = color_by,
        type = "discrete"
    )
}

.pb_resolve_plotly_symbol_mapping <- function(plot_df, shape_by) {
    if (is.null(shape_by)) {
        return(NULL)
    }
    if (length(shape_by) > 1) {
        warning("Shaping by the first column specified")
        shape_by <- shape_by[1]
    }
    if (!shape_by %in% names(plot_df)) {
        stop(sprintf("Shaping column '%s' not found in the data used for plotting.", shape_by))
    }

    column <- plot_df[[shape_by]]
    if (!is.factor(column) && !is.character(column)) {
        column <- as.factor(column)
    } else {
        column <- as.factor(column)
    }

    symbol_pool <- c(
        "circle", "square", "diamond", "cross", "x",
        "triangle-up", "triangle-down", "triangle-left", "triangle-right",
        "star", "hexagon", "hexagon2", "hourglass"
    )
    n_levels <- length(levels(column))
    if (n_levels > length(symbol_pool)) {
        warning("Not enough unique plotly symbols; symbols will be recycled.")
    }
    symbols <- rep(symbol_pool, length.out = n_levels)
    names(symbols) <- levels(column)

    list(
        aes_column = column,
        symbols = symbols
    )
}

.pb_build_embedding_hover_text <- function(plot_df, columns) {
    columns <- unique(columns)
    columns <- columns[columns %in% names(plot_df)]
    if (!length(columns)) {
        return(rep("", nrow(plot_df)))
    }
    apply(plot_df[, columns, drop = FALSE], 1, function(row) {
        pieces <- sprintf("%s: %s", columns, as.character(row))
        paste(pieces, collapse = "<br>")
    })
}
