.pb_assays_to_plot <- function(object, pbf_name) {
    stopifnot(is(object, "ProBatchFeatures"))
    if (is.null(pbf_name) || length(pbf_name) == 0L) {
        assays <- names(object)
    } else {
        assays <- unique(pbf_name)
    }
    if (!length(assays)) {
        stop("Provide at least one `pbf_name` to plot.")
    }
    assays
}

.pb_split_arg_by_assay <- function(arg, assays) {
    n <- length(assays)
    res <- vector("list", n)
    if (is.null(arg)) {
        return(res)
    }
    has_names <- !is.null(names(arg))
    arg_len <- length(arg)
    list_like <- is.list(arg) && !is.data.frame(arg)
    for (i in seq_len(n)) {
        assay <- assays[[i]]
        val <- NULL
        if (is.data.frame(arg) || is.matrix(arg) || inherits(arg, "DataFrame")) {
            val <- arg
        } else if (list_like) {
            if (has_names && assay %in% names(arg)) {
                val <- arg[[assay]]
            } else if (arg_len >= i) {
                val <- arg[[i]]
            } else if (arg_len >= 1L) {
                val <- arg[[arg_len]]
            }
        } else {
            if (has_names && assay %in% names(arg)) {
                val <- arg[[assay]]
            } else if (arg_len >= i) {
                val <- arg[[i]]
            } else if (arg_len >= 1L) {
                val <- arg[[arg_len]]
            }
        }
        res[[i]] <- val
    }
    res
}

.pb_default_title <- function(assay) {
    parts <- strsplit(assay, "::", fixed = TRUE)[[1]]
    if (length(parts) == 2L) {
        paste(parts, collapse = ", ")
    } else {
        assay
    }
}

.pb_resolve_titles <- function(assays, plot_title, default_fun = .pb_default_title) {
    title_list <- .pb_split_arg_by_assay(plot_title, assays)
    vapply(seq_along(assays), function(i) {
        val <- title_list[[i]]
        if (is.null(val) || length(val) == 0L) {
            return(default_fun(assays[[i]]))
        }
        val_chr <- as.character(val)
        if (!length(val_chr) || !nzchar(val_chr[1])) {
            default_fun(assays[[i]])
        } else {
            val_chr[1]
        }
    }, character(1), USE.NAMES = FALSE)
}

.pb_prepare_multi_assay <- function(object, pbf_name, dots, plot_title,
                                    default_title_fun = .pb_default_title,
                                    set_silent = FALSE) {
    assays <- .pb_assays_to_plot(object, pbf_name)

    filename_list <- NULL
    if ("filename" %in% names(dots)) {
        filename_list <- .pb_split_arg_by_assay(dots$filename, assays)
        dots$filename <- NULL
    }

    if (isTRUE(set_silent) && length(assays) > 1L && !"silent" %in% names(dots)) {
        dots$silent <- TRUE
    }

    list(
        assays = assays,
        dots = dots,
        titles = .pb_resolve_titles(assays, plot_title, default_fun = default_title_fun),
        filename_list = filename_list,
        split_arg = function(arg) .pb_split_arg_by_assay(arg, assays)
    )
}

.pb_per_assay_dots <- function(dots, filename_list, index) {
    call_args <- dots
    if (!is.null(filename_list)) {
        fn <- filename_list[[index]]
        if (!is.null(fn)) {
            call_args$filename <- fn
        }
    }
    call_args
}

.pb_arrange_plot_list <- function(plot_list, convert_fun = NULL, draw = TRUE, plot_ncol = NULL, return_gridExtra = FALSE) {
    if (!length(plot_list)) {
        return(invisible(NULL))
    }
    keep <- !vapply(plot_list, is.null, logical(1))
    plot_list <- plot_list[keep]
    if (!length(plot_list)) {
        return(invisible(NULL))
    }
    if (length(plot_list) == 1L) {
        return(plot_list[[1L]])
    }
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
        stop("Install the `gridExtra` package to arrange multiple plots: install.packages(\"gridExtra\").")
    }
    grobs <- if (is.null(convert_fun)) {
        plot_list
    } else {
        lapply(plot_list, convert_fun)
    }
    n <- length(grobs)
    ncol <- if (is.null(plot_ncol)) ceiling(sqrt(n)) else plot_ncol
    nrow <- ceiling(n / ncol)
    arranged <- do.call(
        gridExtra::arrangeGrob,
        list(grobs = grobs, nrow = nrow, ncol = ncol)
    )
    if (isTRUE(draw)) {
        grid::grid.newpage()
        grid::grid.draw(arranged)
    }

    if (isTRUE(return_gridExtra)) {
        return(invisible(list(grob = arranged, plots = plot_list)))
    } else {
        if (requireNamespace("ggplotify", quietly = TRUE)) {
            out <- ggplotify::as.ggplot(arranged) # true ggplot object
            return(invisible(out))
        } else if (requireNamespace("cowplot", quietly = TRUE)) {
            out <- cowplot::ggdraw(arranged) # ggplot object wrapper
            return(invisible(out))
        } else {
            # Last-resort: return the TableGrob (user can grid.draw() it or ggsave() still works)
            message("Returning a grid TableGrob object instead. Use grid.draw() it or ggsave() to plot or save. \nInstall the `ggplotify` or `cowplot` package to get a ggplot object instead.")
            return(invisible(arranged))
        }
    }
}

.pb_align_matrix_and_annotation <- function(data_matrix,
                                            sample_annotation,
                                            sample_id_col,
                                            check_args = list(),
                                            allow_partial_annotation = TRUE) {
    df_long <- matrix_to_long(data_matrix, sample_id_col = sample_id_col)
    check_call <- modifyList(list(
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        df_long = df_long,
        merge = FALSE
    ), check_args)
    df_long <- do.call(check_sample_consistency, check_call)
    data_matrix <- long_to_matrix(df_long, sample_id_col = sample_id_col)

    sample_ids <- colnames(data_matrix)
    annotation <- sample_annotation
    if (!is.null(annotation)) {
        annotation <- .pb_prepare_annotation_for_samples(
            annotation,
            sample_id_col = sample_id_col,
            sample_ids = sample_ids,
            allow_partial = allow_partial_annotation
        )
    }

    list(
        data_matrix = data_matrix,
        sample_annotation = annotation,
        sample_ids = sample_ids
    )
}

.pb_handle_missing_wrapper <- function(data_matrix, warning_message, fill_the_missing, drop_on_false = FALSE) {
    if (!anyNA(data_matrix)) {
        return(data_matrix)
    }

    if (isFALSE(fill_the_missing) && isTRUE(drop_on_false)) {
        data_matrix <- data_matrix[complete.cases(data_matrix), , drop = FALSE]
    }

    handle_missing_values(
        data_matrix = data_matrix,
        warning_message = warning_message,
        fill_the_missing = fill_the_missing
    )
}

.pb_resolve_color_list <- function(color_list,
                                   annotation_df,
                                   id_col,
                                   columns,
                                   warn_message = NULL,
                                   numeric_columns = NULL,
                                   guess_factors = TRUE) {
    if (is.null(columns) || !length(columns)) {
        return(color_list)
    }

    existing <- color_list
    if (is.null(existing)) {
        existing <- list()
    }

    missing_cols <- setdiff(columns, names(existing))
    if (length(missing_cols) && !is.null(annotation_df)) {
        if (!is.null(warn_message)) {
            warning(warn_message)
        }
        inferred <- sample_annotation_to_colors(
            sample_annotation = annotation_df,
            sample_id_col = id_col,
            factor_columns = missing_cols,
            numeric_columns = numeric_columns,
            guess_factors = guess_factors
        )
        existing <- c(existing, inferred)
    }

    keep <- intersect(names(existing), columns)
    if (!length(keep)) {
        return(NULL)
    }
    existing[keep]
}

.pb_prepare_annotation_df <- function(annotation_df,
                                      id_col,
                                      columns,
                                      auto_columns = NULL) {
    if (is.null(annotation_df)) {
        return(list(df = NULL, columns = NULL))
    }
    if (is.null(id_col) || !id_col %in% names(annotation_df)) {
        return(list(df = NULL, columns = NULL))
    }

    target_columns <- columns
    if (is.null(target_columns) && !is.null(auto_columns)) {
        target_columns <- intersect(auto_columns, names(annotation_df))
        if (!length(target_columns)) {
            target_columns <- NULL
        }
    }

    ann <- annotation_df
    if (!is.null(target_columns)) {
        ann <- ann %>%
            select(all_of(c(id_col, target_columns)))
    }

    ann <- ann %>%
        mutate_if(is.POSIXct, as.numeric) %>%
        remove_rownames() %>%
        column_to_rownames(var = id_col)

    if (is.data.frame(ann) && ncol(ann) == 0) {
        ann <- NULL
    }

    list(df = ann, columns = target_columns)
}

.pb_filter_annotation_colors <- function(color_list, annotation_df) {
    if (!is.list(color_list) || is.null(annotation_df)) {
        return(list())
    }
    keep <- intersect(names(color_list), colnames(annotation_df))
    if (!length(keep)) {
        return(list())
    }
    color_list[keep]
}

.pb_prepare_pheatmap_annotations <- function(data_matrix,
                                             column_annotation_df,
                                             row_annotation_df,
                                             col_ann_id_col,
                                             row_ann_id_col,
                                             columns_for_cols,
                                             columns_for_rows,
                                             annotation_color_cols,
                                             annotation_color_rows) {
    col_info <- .pb_prepare_annotation_df(
        column_annotation_df,
        id_col = col_ann_id_col,
        columns = columns_for_cols,
        auto_columns = intersect(
            c("MS_batch", "Diet", "DateTime", "order"),
            names(column_annotation_df)
        )
    )
    annotation_col <- col_info$df
    columns_for_cols <- col_info$columns

    row_info <- .pb_prepare_annotation_df(
        row_annotation_df,
        id_col = row_ann_id_col,
        columns = columns_for_rows,
        auto_columns = intersect(
            c("KEGG_pathway", "WGCNA_module", "evolutionary_distance"),
            names(row_annotation_df)
        )
    )
    annotation_row <- row_info$df
    columns_for_rows <- row_info$columns

    if (!is.null(annotation_col) && !setequal(rownames(annotation_col), colnames(data_matrix))) {
        warning("coloring by column annotation will not work: annotation rownames do not match data matrix column names")
    }

    if (!is.null(annotation_row) && !setequal(rownames(annotation_row), rownames(data_matrix))) {
        warning("coloring by row annotation will not work: annotation rownames do not match data matrix column names")
    }

    annotation_color_cols <- .pb_filter_annotation_colors(annotation_color_cols, annotation_col)
    annotation_color_rows <- .pb_filter_annotation_colors(annotation_color_rows, annotation_row)

    annotation_color_list <- c(annotation_color_cols, annotation_color_rows)
    if (!length(annotation_color_list)) {
        annotation_color_list <- NA
    }

    list(
        annotation_col = annotation_col,
        annotation_row = annotation_row,
        annotation_color_list = annotation_color_list,
        columns_for_cols = columns_for_cols,
        columns_for_rows = columns_for_rows
    )
}

.pb_open_graphics_device <- function(filename, width, height, units, plot_title = NULL, png_res = 300) {
    if (is.null(filename)) {
        return(list(opened = FALSE, close = function() invisible(NULL)))
    }

    units_adjusted <- adjust_units(units, width, height)
    width <- units_adjusted$width
    height <- units_adjusted$height
    unit <- units_adjusted$unit

    if (is.na(width)) {
        width <- 7
    }
    if (is.na(height)) {
        height <- 7
    }

    ext <- tolower(tools::file_ext(filename))
    if (!nzchar(ext)) {
        stop("filename must have an extension")
    }

    if (ext == "pdf") {
        grDevices::pdf(file = filename, width = width, height = height, title = plot_title)
    } else if (ext == "png") {
        grDevices::png(filename = filename, width = width, height = height, units = unit, res = png_res)
    } else {
        stop("currently only pdf and png extensions for filename are implemented")
    }

    list(
        opened = TRUE,
        close = function() grDevices::dev.off()
    )
}

.pb_prepare_embedding_inputs <- function(data_matrix,
                                         sample_annotation,
                                         sample_id_col,
                                         feature_id_col,
                                         color_by = NULL,
                                         fill_the_missing = -1,
                                         warning_message = NULL,
                                         allow_partial_annotation = FALSE,
                                         check_args = list(),
                                         drop_on_false = TRUE) {
    alignment <- .pb_align_matrix_and_annotation(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        check_args = check_args,
        allow_partial_annotation = allow_partial_annotation
    )

    data_matrix <- alignment$data_matrix
    sample_annotation <- alignment$sample_annotation

    data_matrix <- check_feature_id_col_in_dm(feature_id_col, data_matrix)

    if (!is.null(warning_message)) {
        data_matrix <- .pb_handle_missing_wrapper(
            data_matrix = data_matrix,
            warning_message = warning_message,
            fill_the_missing = fill_the_missing,
            drop_on_false = drop_on_false
        )
    }

    list(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_ids = alignment$sample_ids
    )
}
