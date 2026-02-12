.pb_assays_to_plot <- function(object, pbf_name) {
    .pb_resolve_assays_for_input(
        object = object,
        pbf_name = pbf_name,
        default = "all",
        deduplicate = TRUE
    )
}

.pb_split_arg_by_assay <- function(arg, assays) {
    n <- length(assays)
    res <- vector("list", n)

    if (!n || is.null(arg)) {
        return(res)
    }

    if (is.data.frame(arg) || is.matrix(arg) || inherits(arg, "DataFrame")) {
        for (i in seq_len(n)) {
            res[[i]] <- arg
        }
        return(res)
    }

    list_like <- is.list(arg) && !is.data.frame(arg)

    if (!list_like && is.atomic(arg)) {
        # Treat atomic vectors (e.g., character vectors of column names)
        # as shared arguments unless explicitly named.
        if (length(arg) > 1L && is.null(names(arg))) {
            for (i in seq_len(n)) {
                res[[i]] <- arg
            }
            return(res)
        }
        values <- as.list(arg)
    } else {
        values <- if (list_like) arg else as.list(arg)
    }

    arg_len <- length(values)
    names_present <- !is.null(names(values))

    for (i in seq_len(n)) {
        assay <- assays[[i]]
        val <- NULL

        if (names_present && assay %in% names(values)) {
            val <- values[[assay]]
        } else if (arg_len >= i) {
            val <- values[[i]]
        } else if (arg_len >= 1L) {
            val <- values[[arg_len]]
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

.pb_break_long_assay_title <- function(title, max_length = 35L) {
    if (is.na(title) || !nzchar(title) || nchar(title) <= max_length) {
        return(title)
    }

    hyphen_pos <- gregexpr("-", title, fixed = TRUE)[[1]]
    hyphen_pos <- hyphen_pos[hyphen_pos > 1 & hyphen_pos < nchar(title)]
    if (!length(hyphen_pos)) {
        return(title)
    }

    mid_point <- nchar(title) / 2
    idx <- hyphen_pos[which.min(abs(hyphen_pos - mid_point))]
    left <- substr(title, 1, idx)
    right <- substr(title, idx + 1, nchar(title))
    right <- trimws(right)

    paste0(left, "\n\t", right)
}

.pb_intragroup_assay_label_parts <- function(title) {
    if (is.null(title) || !nzchar(title) || !grepl("::", title, fixed = TRUE)) {
        return(list(level = NULL, steps = character()))
    }
    parts <- strsplit(title, "::", fixed = TRUE)[[1]]
    if (length(parts) < 2L) {
        return(list(level = NULL, steps = character()))
    }
    level <- trimws(parts[[1]])
    remainder <- trimws(paste(parts[-1], collapse = "::"))
    if (!nzchar(remainder)) {
        return(list(level = level, steps = character()))
    }
    steps <- strsplit(remainder, "_", fixed = TRUE)[[1]]
    steps <- trimws(steps)
    steps <- steps[nzchar(steps)]
    if (!length(steps)) {
        return(list(level = level, steps = character()))
    }
    connectors <- tolower(steps) == "on"
    steps <- steps[!connectors]
    list(level = level, steps = steps)
}

.pb_intragroup_assay_level <- function(title) {
    parts <- .pb_intragroup_assay_label_parts(title)
    if (is.null(parts$level) || !nzchar(parts$level)) {
        return("")
    }
    parts$level
}

.pb_intragroup_split_long_underscore <- function(line, max_len = 10L) {
    if (is.null(line) || !nzchar(line) || nchar(line) <= max_len || !grepl("_", line, fixed = TRUE)) {
        return(line)
    }
    parts <- strsplit(line, "_", fixed = TRUE)[[1]]
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    if (!length(parts)) {
        return(line)
    }
    parts
}

.pb_intragroup_format_assay_label <- function(title, include_level = TRUE) {
    if (is.null(title) || !nzchar(title)) {
        return(title)
    }
    parts <- .pb_intragroup_assay_label_parts(title)
    lines <- character()
    if (include_level && !is.null(parts$level) && nzchar(parts$level)) {
        lines <- c(lines, paste0(parts$level, " ::"))
    }
    if (length(parts$steps)) {
        lines <- c(lines, parts$steps)
    }
    if (!length(lines) && nchar(title) > 10L && grepl("_", title, fixed = TRUE)) {
        splits <- strsplit(title, "_", fixed = TRUE)[[1]]
        splits <- trimws(splits)
        splits <- splits[nzchar(splits)]
        if (length(splits)) {
            lines <- splits
        }
    }
    if (!length(lines)) {
        return(.pb_break_long_assay_title(title))
    }
    split_lines <- unlist(lapply(lines, .pb_intragroup_split_long_underscore), use.names = FALSE)
    paste(split_lines, collapse = "\n")
}

.pb_refactor_assay_titles <- function(titles, use_shared_title = TRUE, max_length = 35L) {
    if (!length(titles)) {
        return(list(titles = titles, shared_title = NULL))
    }

    titles <- gsub("_on_", "-", titles, fixed = TRUE)
    titles <- trimws(titles)
    original_titles <- titles

    shared_title <- NULL

    if (use_shared_title && length(titles) > 1L) {
        split_titles <- strsplit(titles, "-", fixed = TRUE)
        split_titles <- lapply(split_titles, trimws)
        min_len <- min(lengths(split_titles))

        shared_parts <- character()
        if (min_len > 0) {
            for (i in seq_len(min_len)) {
                parts_i <- vapply(split_titles, function(x) x[[i]], character(1))
                parts_i <- trimws(parts_i)
                unique_parts <- unique(parts_i[nzchar(parts_i)])
                if (length(unique_parts) == 1L) {
                    shared_parts <- c(shared_parts, unique_parts)
                } else {
                    break
                }
            }
        }

        if (length(shared_parts)) {
            shared_title <- paste(shared_parts, collapse = " - ")
            remove_len <- length(shared_parts)
            titles <- vapply(seq_along(split_titles), function(idx) {
                parts <- split_titles[[idx]]
                if (length(parts) > remove_len) {
                    remainder <- parts[-seq_len(remove_len)]
                    remainder <- remainder[nzchar(remainder)]
                    if (length(remainder)) {
                        return(paste(remainder, collapse = " - "))
                    }
                }
                original_titles[[idx]]
            }, character(1), USE.NAMES = FALSE)
        }
    }

    titles <- trimws(titles)
    titles <- vapply(titles, .pb_break_long_assay_title, character(1), max_length = max_length)

    list(titles = titles, shared_title = shared_title)
}

.pb_plot_title_values <- function(title) {
    if (is.null(title)) {
        return(character())
    }
    values <- unlist(title, recursive = TRUE, use.names = FALSE)
    if (!length(values)) {
        return(character())
    }

    values <- as.character(values)
    values <- values[!is.na(values)]
    values <- trimws(values)
    values[nzchar(values)]
}

.pb_collapsed_plot_title <- function(title) {
    values <- .pb_plot_title_values(title)
    if (!length(values)) {
        return(NULL)
    }
    paste(values, collapse = "\n")
}

.pb_first_plot_title <- function(title) {
    values <- .pb_plot_title_values(title)
    if (!length(values)) {
        return(NULL)
    }
    values[[1]]
}

.pb_prepare_multi_assay <- function(object, pbf_name, dots, plot_title,
                                    default_title_fun = .pb_default_title,
                                    set_silent = FALSE,
                                    refactor_titles = TRUE) {
    assays <- .pb_assays_to_plot(object, pbf_name)

    filename_list <- NULL
    if ("filename" %in% names(dots)) {
        filename_list <- .pb_split_arg_by_assay(dots$filename, assays)
        dots$filename <- NULL
    }

    if (isTRUE(set_silent) && length(assays) > 1L && !"silent" %in% names(dots)) {
        dots$silent <- TRUE
    }

    resolved_titles <- .pb_resolve_titles(assays, plot_title, default_fun = default_title_fun)
    if (isTRUE(refactor_titles)) {
        title_info <- .pb_refactor_assay_titles(
            titles = resolved_titles,
            use_shared_title = is.null(plot_title)
        )
    } else {
        title_info <- list(titles = resolved_titles, shared_title = NULL)
    }

    list(
        assays = assays,
        dots = dots,
        titles = title_info$titles,
        shared_title = title_info$shared_title,
        filename_list = filename_list,
        split_arg = function(arg) .pb_split_arg_by_assay(arg, assays)
    )
}

.pb_attach_shared_title <- function(plot_list, shared_title) {
    if (!is.null(shared_title) && nzchar(shared_title)) {
        attr(plot_list, "pb_shared_title") <- shared_title
    } else if (!is.null(attr(plot_list, "pb_shared_title", exact = TRUE))) {
        attr(plot_list, "pb_shared_title") <- NULL
    }
    plot_list
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

.pb_prepare_shape_column <- function(shape_by, sample_annotation, data_label = "sample_annotation") {
    if (is.null(shape_by) || !length(shape_by)) {
        return(list(shape_by = NULL, sample_annotation = sample_annotation))
    }

    if (length(shape_by) > 1) {
        warning("Shaping by the first column specified")
        shape_by <- shape_by[1]
    }

    if (!shape_by %in% colnames(sample_annotation)) {
        stop(sprintf("Shaping column '%s' not found in %s", shape_by, data_label))
    }

    shape_column <- sample_annotation[[shape_by]]
    if (!is.factor(shape_column) && !is.character(shape_column)) {
        sample_annotation[[shape_by]] <- as.factor(shape_column)
    }

    list(shape_by = shape_by, sample_annotation = sample_annotation)
}

.pb_pop_use_plotlyrender <- function(dots) {
    if (is.null(dots)) {
        dots <- list()
    }

    use_plotlyrender <- isTRUE(dots$use_plotlyrender)
    dots$use_plotlyrender <- NULL

    list(use_plotlyrender = use_plotlyrender, dots = dots)
}

.pb_finalize_embedding_collection <- function(plot_list,
                                              use_plotlyrender,
                                              return_gridExtra,
                                              plot_ncol,
                                              return_subplots,
                                              subplot_ncol,
                                              share_axes) {
    shared_title <- attr(plot_list, "pb_shared_title", exact = TRUE)
    if (!length(plot_list)) {
        return(invisible(NULL))
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
            subplot_obj <- do.call(plotly::subplot, subplot_args)
            if (!is.null(shared_title) && nzchar(shared_title)) {
                subplot_obj <- plotly::layout(
                    subplot_obj,
                    title = list(text = shared_title)
                )
            }
            return(subplot_obj)
        }

        return(.pb_attach_shared_title(plot_list, shared_title))
    }

    if (isTRUE(return_subplots)) {
        warning("return_subplots = TRUE is only supported when use_plotlyrender = TRUE; arranging ggplot outputs instead.")
    }

    .pb_arrange_plot_list(
        plot_list,
        convert_fun = ggplotGrob,
        plot_ncol = plot_ncol,
        return_gridExtra = return_gridExtra
    )
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
        lapply(plot_list, function(x) {
            if (grid::is.grob(x) || inherits(x, "gtable")) {
                x
            } else {
                convert_fun(x)
            }
        })
    }
    n <- length(grobs)
    ncol <- if (is.null(plot_ncol)) ceiling(sqrt(n)) else plot_ncol
    nrow <- ceiling(n / ncol)
    shared_title <- attr(plot_list, "pb_shared_title", exact = TRUE)
    top_arg <- list()
    if (!is.null(shared_title) && nzchar(shared_title)) {
        top_arg <- list(top = grid::textGrob(
            shared_title,
            gp = grid::gpar(fontface = "bold")
        ))
    }
    arranged <- do.call(
        gridExtra::arrangeGrob,
        c(list(grobs = grobs, nrow = nrow, ncol = ncol), top_arg)
    )
    if (isTRUE(draw)) {
        grid::grid.newpage()
        grid::grid.draw(arranged)
    }

    if (isTRUE(return_gridExtra)) {
        return(invisible(list(grob = arranged, plots = plot_list)))
    } else {
        return(.pb_convert_arranged_grob(arranged, invisible_out = TRUE))
    }
}

.pb_convert_arranged_grob <- function(arranged, invisible_out = FALSE) {
    if (requireNamespace("ggplotify", quietly = TRUE)) {
        out <- ggplotify::as.ggplot(arranged)
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
        out <- cowplot::ggdraw(arranged)
    } else {
        message("Returning a grid TableGrob object instead. Use grid.draw() it or ggsave() to plot or save. \nInstall the `ggplotify` or `cowplot` package to get a ggplot object instead.")
        out <- arranged
    }

    if (isTRUE(invisible_out)) {
        return(invisible(out))
    }
    out
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
                                      auto_columns = NULL,
                                      target_ids = NULL) {
    if (is.null(annotation_df)) {
        return(list(df = NULL, columns = NULL, missing_ids = character()))
    }

    annotation_df <- as.data.frame(annotation_df, stringsAsFactors = FALSE)
    if (is.null(id_col) || !id_col %in% names(annotation_df)) {
        rn <- rownames(annotation_df)
        if (!is.null(rn) && !anyNA(rn) && length(rn) == nrow(annotation_df)) {
            id_col <- if (!is.null(id_col) && nzchar(id_col)) {
                id_col
            } else {
                ".pb_annotation_id"
            }
            annotation_df[[id_col]] <- rn
        } else {
            return(list(
                df = NULL,
                columns = NULL,
                missing_ids = if (is.null(target_ids)) character() else as.character(target_ids)
            ))
        }
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

    ann_ids <- as.character(ann[[id_col]])
    keep <- !is.na(ann_ids) & nzchar(ann_ids) & !duplicated(ann_ids)
    ann <- ann[keep, , drop = FALSE]
    ann_ids <- as.character(ann[[id_col]])

    missing_ids <- character()
    if (!is.null(target_ids)) {
        target_ids <- as.character(target_ids)
        idx <- match(target_ids, ann_ids)
        missing_ids <- target_ids[is.na(idx)]
        ann <- ann[idx[!is.na(idx)], , drop = FALSE]
    }

    ann <- ann %>%
        mutate(across(where(is.POSIXct), as.numeric)) %>%
        remove_rownames() %>%
        column_to_rownames(var = id_col)

    if (is.data.frame(ann)) {
        ann[] <- lapply(ann, function(col) {
            if (is.factor(col)) {
                droplevels(col)
            } else {
                col
            }
        })
    }

    if (is.data.frame(ann) && (ncol(ann) == 0 || nrow(ann) == 0)) {
        ann <- NULL
    }

    list(df = ann, columns = target_columns, missing_ids = unique(missing_ids))
}

.pb_filter_annotation_colors <- function(color_list, annotation_df) {
    if (!is.list(color_list) || is.null(annotation_df)) {
        return(list())
    }
    keep <- intersect(names(color_list), colnames(annotation_df))
    if (!length(keep)) {
        return(list())
    }
    aligned <- vector("list", length(keep))
    names(aligned) <- keep

    for (col in keep) {
        colors <- color_list[[col]]
        values <- annotation_df[[col]]
        if (is.null(values)) {
            next
        }
        if (is.character(values)) {
            values <- factor(values)
        }
        if (is.factor(values) && !is.null(names(colors))) {
            levels_needed <- levels(droplevels(values))
            missing <- setdiff(levels_needed, names(colors))
            if (length(missing)) {
                extra_cols <- grDevices::hcl.colors(length(missing), "Set 3")
                names(extra_cols) <- missing
                colors <- c(colors, extra_cols)
            }
            colors <- colors[levels_needed]
        }
        aligned[[col]] <- colors
    }
    aligned
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
        ),
        target_ids = colnames(data_matrix)
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
        ),
        target_ids = rownames(data_matrix)
    )
    annotation_row <- row_info$df
    columns_for_rows <- row_info$columns

    if (!is.null(annotation_col) && length(col_info$missing_ids)) {
        warning(
            "column annotation is missing entries for ",
            length(col_info$missing_ids),
            " matrix columns"
        )
    }

    if (!is.null(annotation_row) && length(row_info$missing_ids)) {
        warning(
            "row annotation is missing entries for ",
            length(row_info$missing_ids),
            " matrix rows"
        )
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
