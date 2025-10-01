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
        arrangeGrob,
        list(grobs = grobs, nrow = nrow, ncol = ncol)
    )
    if (isTRUE(draw)) {
        grid.newpage()
        grid.draw(arranged)
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
