#' @title Plot per-sample mean or boxplots for initial assessment
#' @description Plot per-sample mean or boxplots (showing median and quantiles). In ordered samples,
#' e.g. consecutive MS runs, order-associated effects are visualised.
#' @inheritParams proBatch
#' @details functions for quick visual assessment of trends associated, overall
#'   or specific covariate-associated (see \code{batch_col} and \code{facet_col})
#' @inheritParams proBatch
#' @param color_scheme named vector, names corresponding to unique batch values of
#'  \code{batch_col} in \code{sample_annotation}. Best created with \link{sample_annotation_to_colors}
#' @param vline_color color of vertical lines, typically denoting
#'  different MS batches in ordered runs; should be \code{NULL} for experiments without intrinsic order
#' @param ylimits range of y-axis to compare two plots side by side, if required.
#' @param outliers keep (default) or remove the boxplot outliers
#' @param theme_name Name of the ggplot theme to apply to the resulting plot.
#' @param x Input object supplied to the generics (matrix, long data frame, or `ProBatchFeatures`).
#' @param pbf_name Assay name(s) used when `x` is a `ProBatchFeatures`.
#' @param plot_ncol Number of columns when arranging multiple assay plots.
#' @param return_gridExtra Logical; return arranged grobs instead of a plot list.
#' @param ... Additional arguments forwarded between methods.
#'
#' @return ggplot2 class object. Thus, all aesthetics can be overridden
#'
#' @seealso \code{\link[ggplot2]{ggplot}}, \link{date_to_sample_order}
#' @name plot_sample_mean_or_boxplot
#'
#' @examples
#' data(list = c(
#'     "example_proteome", "example_sample_annotation",
#'     "example_proteome_matrix"
#' ), package = "proBatch")
#'
#' demo_ids <- colnames(example_proteome_matrix)[1:6]
#' demo_matrix <- example_proteome_matrix[, demo_ids]
#' demo_annotation <- example_sample_annotation[
#'     example_sample_annotation$FullRunName %in% demo_ids,
#' ]
#'
#' plot_sample_mean(
#'     demo_matrix,
#'     demo_annotation,
#'     order_col = "order",
#'     batch_col = "MS_batch"
#' )
#'
#' demo_proteome <- example_proteome[
#'     example_proteome$FullRunName %in% demo_ids,
#' ]
#' plot_boxplot(
#'     demo_proteome,
#'     sample_annotation = demo_annotation,
#'     batch_col = "MS_batch"
#' )
#'
.pb_apply_initial_assessment_theme <- function(gg, theme_name, base_size) {
    theme_name <- match.arg(theme_name, choices = c("classic", "minimal", "bw", "light", "dark"))
    if (identical(theme_name, "classic")) {
        return(gg + theme_classic(base_size = base_size))
    }
    if (identical(theme_name, "minimal")) {
        return(gg + theme_minimal(base_size = base_size))
    }
    if (identical(theme_name, "bw")) {
        return(gg + theme_bw(base_size = base_size))
    }
    if (identical(theme_name, "light")) {
        return(gg + theme_light(base_size = base_size))
    }
    gg + theme_dark(base_size = base_size)
}

.pb_finalize_initial_assessment_plot <- function(gg, order_values, ylimits = NULL,
                                                 rotate_x = FALSE,
                                                 color_by_batch = FALSE,
                                                 is_factor = FALSE) {
    if (!is.null(ylimits)) {
        gg <- gg + coord_cartesian(ylim = ylimits)
    }
    if (isTRUE(rotate_x)) {
        gg <- gg + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    }
    if (length(unique(order_values)) > 30 && color_by_batch && is_factor) {
        gg <- gg + theme(legend.position = "top")
    }
    gg
}

#' @rdname plot_sample_mean_or_boxplot
#' @method plot_sample_mean default
#' @export
plot_sample_mean.default <- function(x, sample_annotation,
                                     sample_id_col = "FullRunName",
                                     batch_col = "MS_batch",
                                     color_by_batch = FALSE,
                                     color_scheme = "brewer",
                                     order_col = "order",
                                     vline_color = "grey",
                                     facet_col = NULL,
                                     filename = NULL, width = NA, height = NA,
                                     units = c("cm", "in", "mm"),
                                     plot_title = NULL,
                                     theme_name = c("classic", "minimal", "bw", "light", "dark"),
                                     base_size = 20,
                                     ylimits = NULL,
                                     pbf_name = NULL,
                                     ...) {
    data_matrix <- x
    # stop early if sample_annotation missing
    if (is.null(sample_annotation)) {
        stop("`sample_annotation` must be provided.")
    }
    # require data_matrix column names
    if (is.null(colnames(data_matrix))) {
        stop("`data_matrix` must have column names representing samples.")
    }

    # Create a data frame with sample averages
    sample_average <- colMeans(data_matrix, na.rm = TRUE)
    df_ave <- data.frame(
        Mean_Intensity = sample_average,
        temp_id = colnames(data_matrix),
        stringsAsFactors = FALSE
    )
    names(df_ave)[names(df_ave) == "temp_id"] <- sample_id_col

    # Check the consistency of sample ann. sample IDs and measur. table sample IDs
    df_ave <- check_sample_consistency(
        sample_annotation, sample_id_col, df_ave,
        batch_col, order_col, facet_col
    )
    message("Sample ID is kept only if they have a match in annotation data frame,
            otherwise the sample ID is removed from the plot")

    # initialize is_factor to avoid undefined variable
    is_factor <- FALSE
    # Ensure that batch-coloring-related arguments are defined properly
    # validate batch_col existence, disable coloring if NULL
    if (!is.null(batch_col)) {
        if (!(batch_col %in% names(df_ave))) {
            message("batches cannot be colored as the batch column or sample ID column
                    is not defined, check sample_annotation and data matrix")
            stop("Batch column '", batch_col, "' not found in data.")
        }
        # if coloring by batch is true and color scheme is df and contains column with name of batch_col, keep only this column
        if (color_by_batch && (batch_col %in% names(color_scheme))) {
            color_scheme <- color_scheme[[batch_col]]
        }
    } else if (color_by_batch) {
        message("batches cannot be colored as the batch column is defined as NULL,
                continuing without colors")
        warning("`batch_col` is NULL; disabling `color_by_batch`")
        color_by_batch <- FALSE
    }

    # For order definition and subsequent faceting, facet column has to be in the df
    if (!is.null(facet_col) && !(facet_col %in% names(df_ave))) {
        message(sprintf(
            '"%s" is specified as column for faceting, but is not present
                    in the data, check sample annotation data frame',
            facet_col
        ))
        stop(sprintf("Faceting column '%s' not found in data.", facet_col))
    }

    # Defining sample order for plotting
    sample_order <- define_sample_order(
        order_col, sample_annotation, facet_col,
        batch_col, df_ave,
        sample_id_col, color_by_batch
    )
    order_col <- sample_order$order_col
    df_ave <- sample_order$df_long

    # Convert order to factor if not numeric
    if (!is.numeric(df_ave[[order_col]]) && is.character(df_ave[[order_col]])) {
        df_ave[[order_col]] <- factor(df_ave[[order_col]], levels = unique(df_ave[[order_col]]))
    }

    # Create plot
    # Main plotting of intensity means:
    gg <- ggplot(df_ave, aes(x = !!sym(order_col), y = .data$Mean_Intensity)) +
        geom_point()

    # add colors
    gg <- color_by_factor(
        color_by_batch = color_by_batch,
        batch_col = batch_col, gg = gg,
        color_scheme = color_scheme,
        sample_annotation = df_ave,
        fill_or_color = "color"
    )

    axis_x_label <- if (!is.null(order_col)) {
        order_col
    } else {
        sample_id_col
    }
    gg <- gg + labs(x = axis_x_label, y = "Mean_Intensity")

    # add vertical lines, if required (for order-related effects)
    if (!is.null(batch_col)) {
        batch_vector <- df_ave[[batch_col]]
        is_factor <- is_batch_factor(batch_vector, color_scheme)
    }
    if (!is.null(batch_col) && is_factor && !is.null(vline_color)) {
        gg <- add_vertical_batch_borders(
            order_col, sample_id_col, batch_col,
            vline_color,
            facet_col, df_ave, gg
        )
    }

    # Faceting - plot each "facet factor" in it's own subplot
    if (!is.null(facet_col)) {
        gg <- gg + facet_wrap(as.formula(paste("~", facet_col)),
            dir = "v", scales = "free_x"
        )
    }

    # Add the title and theme
    if (!is.null(plot_title)) {
        gg <- gg + ggtitle(plot_title) +
            theme(plot.title = element_text(face = "bold", hjust = 0.5))
    }
    gg <- .pb_apply_initial_assessment_theme(gg = gg, theme_name = theme_name, base_size = base_size)
    gg <- .pb_finalize_initial_assessment_plot(
        gg = gg,
        order_values = df_ave[[order_col]],
        ylimits = ylimits,
        rotate_x = !is.numeric(df_ave[[order_col]]),
        color_by_batch = color_by_batch,
        is_factor = is_factor
    )

    units <- match.arg(units)
    save_ggplot(filename, units, width, height, gg)
    gg
}


#' @rdname plot_sample_mean_or_boxplot
#' @method plot_boxplot default
#' @export
plot_boxplot.default <- function(x, sample_annotation,
                                 sample_id_col = "FullRunName",
                                 measure_col = "Intensity",
                                 batch_col = "MS_batch",
                                 color_by_batch = TRUE,
                                 color_scheme = "brewer",
                                 order_col = "order",
                                 facet_col = NULL,
                                 filename = NULL, width = NA, height = NA,
                                 units = c("cm", "in", "mm"),
                                 plot_title = NULL,
                                 theme_name = c("classic", "minimal", "bw", "light", "dark"),
                                 base_size = 20,
                                 ylimits = NULL, outliers = TRUE,
                                 pbf_name = NULL,
                                 ...) {
    df_long <- x
    # Validate inputs
    if (is.null(sample_annotation)) {
        stop("`sample_annotation` must be provided.")
    }
    if (!(measure_col %in% names(df_long))) {
        stop("`measure_col` '", measure_col, "' not found in data.")
    }

    # Check the consistency of sample ann. sample IDs and measur. table sample IDs
    df_long <- check_sample_consistency(
        sample_annotation, sample_id_col, df_long,
        batch_col, order_col, facet_col
    )

    # Ensure that batch-coloring-related arguments are defined properly
    is_factor <- FALSE
    if (!is.null(batch_col)) {
        if (!(batch_col %in% names(df_long))) {
            message("batches cannot be colored as the batch column or sample ID column
                    is not defined, check sample_annotation and data matrix")
            stop("Batch column '", batch_col, "' not found in data.")
        }
        # if coloring by batch is true and color scheme is df and contains column with name of batch_col, keep only this column
        if (color_by_batch && (batch_col %in% names(color_scheme))) {
            color_scheme <- color_scheme[[batch_col]]
        }
    } else if (color_by_batch) {
        message("batches cannot be colored as the batch column is defined as NULL,
            continuing without colors")
        warning("`batch_col` is NULL; disabling `color_by_batch`")
        color_by_batch <- FALSE
    }

    # For order definition and subsequent faceting, facet column has to be in the df
    if (!is.null(facet_col) && !(facet_col %in% names(df_long))) {
        message(sprintf(
            '"%s" is specified as column for faceting, but is not present
                    in the data, check sample annotation data frame',
            facet_col
        ))
        stop(sprintf("Faceting column '%s' not found in data.", facet_col))
    }

    order_col_name <- order_col
    # Defining sample order for plotting (even if order_col NULL,
    # it will re-arrange df_long levels as required for plotting)
    # Defining sample order for plotting
    sample_order <- define_sample_order(
        order_col, sample_annotation, facet_col,
        batch_col, df_long,
        sample_id_col, color_by_batch
    )
    order_col <- sample_order$order_col
    df_long <- sample_order$df_long

    # Convert order to factor if appropriate
    if (!is.numeric(df_long[[order_col]]) && is.character(df_long[[order_col]])) {
        df_long[[order_col]] <- factor(df_long[[order_col]],
            levels = unique(df_long[[order_col]])
        )
    }

    # Main plotting of intensity distribution boxplots
    gg <- ggplot(df_long, aes(
        x = !!sym(order_col), y = !!sym(measure_col),
        group = !!sym(order_col)
    ))
    if (outliers) {
        gg <- gg + geom_boxplot(outlier.size = 0.15)
    } else {
        warning("`outliers = FALSE`: outliers will be removed (outlier.shape = NA).")
        gg <- gg + geom_boxplot(outlier.shape = NA)
    }

    # Define the color scheme, Apply fill colors
    gg <- color_by_factor(
        color_by_batch = color_by_batch,
        batch_col = batch_col,
        gg = gg,
        color_scheme = color_scheme,
        sample_annotation = df_long,
        fill_or_color = "fill"
    )

    axis_x_label <- if (!is.null(order_col_name)) {
        order_col_name
    } else if (!is.null(order_col)) {
        order_col
    } else {
        sample_id_col
    }
    gg <- gg + labs(x = axis_x_label, y = measure_col)

    # Plot each "facet factor" in it's own subplot
    if (!is.null(facet_col)) {
        gg <- gg + facet_wrap(as.formula(paste("~", facet_col)),
            dir = "v", scales = "free_x"
        )
    }

    # Add the title
    if (!is.null(plot_title)) {
        gg <- gg + ggtitle(plot_title) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    } else if (!is.null(pbf_name)) {
        gg <- gg + ggtitle(pbf_name) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    }

    if (!is.null(batch_col)) {
        batch_vector <- df_long[[batch_col]]
        is_factor <- is_batch_factor(batch_vector, color_scheme)
    }
    gg <- .pb_apply_initial_assessment_theme(gg = gg, theme_name = theme_name, base_size = base_size)
    gg <- .pb_finalize_initial_assessment_plot(
        gg = gg,
        order_values = df_long[[order_col]],
        ylimits = ylimits,
        rotate_x = !is.numeric(df_long[[order_col]]),
        color_by_batch = color_by_batch,
        is_factor = is_factor
    )

    units <- match.arg(units)
    save_ggplot(filename, units, width, height, gg)
    gg
}


#' @rdname plot_sample_mean_or_boxplot
#' @method plot_sample_mean ProBatchFeatures
#' @export
plot_sample_mean.ProBatchFeatures <- function(x, pbf_name = NULL, plot_title = NULL, ...) {
    object <- x # Use 'x' as per convention

    assay_name <- .pb_resolve_assay_for_input(object = object, pbf_name = pbf_name)
    data_matrix <- pb_assay_matrix(object, assay = assay_name)
    sample_annotation <- .pb_default_sample_annotation(
        object = object,
        sample_annotation = NULL,
        sample_id_col = "FullRunName",
        sample_ids = colnames(data_matrix)
    )

    plot_title <- if (is.null(plot_title)) assay_name else plot_title

    # Call the default method with the extracted data
    plot_sample_mean.default(
        x = data_matrix,
        sample_annotation = sample_annotation,
        plot_title = plot_title,
        ...
    )
}

#' @rdname plot_sample_mean_or_boxplot
#' @method plot_boxplot ProBatchFeatures
#' @export
plot_boxplot.ProBatchFeatures <- function(x, pbf_name = NULL, sample_id_col = NULL, plot_title = NULL, plot_ncol = NULL, return_gridExtra = FALSE, ...) {
    object <- x # Use 'x' as per convention

    if (is.null(sample_id_col)) {
        stop("`sample_id_col` must be provided.")
    }
    assays <- .pb_assays_to_plot(object, pbf_name)
    dots <- list(...)

    filename_list <- NULL
    if ("filename" %in% names(dots)) {
        filename_list <- .pb_split_arg_by_assay(dots$filename, assays)
        dots$filename <- NULL
    }

    sample_annotation <- as.data.frame(colData(object))
    titles <- .pb_resolve_titles(assays, plot_title, default_fun = function(x) x)

    plot_list <- vector("list", length(assays))
    names(plot_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        message("Extracting data from assay: ", assay_nm)
        df_long <- pb_as_long(
            object,
            feature_id_col = "Feature",
            sample_id_col = sample_id_col,
            measure_col = "Intensity",
            pbf_name = assay_nm
        )
        df_long <- df_long[!is.na(df_long$Intensity), ]

        # Drop sample_annotation columns from df_long if they exist to avoid duplication, except sample_id_col
        overlap_cols <- setdiff(intersect(names(sample_annotation), names(df_long)), sample_id_col)
        if (length(overlap_cols) > 0) {
            df_long <- df_long[, !names(df_long) %in% overlap_cols, drop = FALSE]
        }

        call_args <- dots
        if (!is.null(filename_list)) {
            fn <- filename_list[[i]]
            if (!is.null(fn)) {
                call_args$filename <- fn
            }
        }

        call_args <- c(list(
            x = df_long,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            measure_col = "Intensity",
            plot_title = titles[i],
            pbf_name = assay_nm
        ), call_args)

        plot_list[[i]] <- do.call(plot_boxplot.default, call_args)
    }

    .pb_arrange_plot_list(plot_list, plot_ncol = plot_ncol, convert_fun = ggplotGrob, return_gridExtra = return_gridExtra)
}

#' Generic function for plotting per-sample mean or boxplots for initial assessment
#' @rdname plot_sample_mean_or_boxplot
#' @export
plot_sample_mean <- function(x, ...) UseMethod("plot_sample_mean")

#' Generic function for plotting per-sample mean or boxplots for initial assessment
#' @rdname plot_sample_mean_or_boxplot
#' @export
plot_boxplot <- function(x, ...) UseMethod("plot_boxplot")
