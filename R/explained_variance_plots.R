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
#' @param pbf_name Assay name(s) used when `data_matrix` is a `ProBatchFeatures`.
#' @param ... Additional arguments forwarded between methods.
#'
#' @name calculate_PVCA
#' @return data frame of weights of Principal Variance Components
#' @export
#'
#' @examples
#' data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#' matrix_test <- na.omit(example_proteome_matrix)[1:50, ]
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
    if (!is.null(factors_for_PVCA)) {
        factors_for_PVCA <- unique(factors_for_PVCA)
        factors_for_PVCA <- factors_for_PVCA[factors_for_PVCA %in% names(sample_annotation)]
        if (length(factors_for_PVCA)) {
            level_counts <- vapply(sample_annotation[factors_for_PVCA], function(x) {
                length(unique(x[!is.na(x)]))
            }, integer(1))
            factors_for_PVCA <- factors_for_PVCA[level_counts > 1L]
        }
        if (!length(factors_for_PVCA)) {
            stop("No PVCA factors with more than one sampled level.")
        }
    }
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
    pvcaAssess <- tryCatch(
        pvcaBatchAssess(expr_set, factors_for_PVCA, threshold = pca_threshold),
        error = function(e) {
            msg <- conditionMessage(e)
            if (!grepl("minDataPointsPerStratum", msg, fixed = TRUE) &&
                !grepl("too small for reliable estimation", msg, fixed = TRUE)) {
                stop(e)
            }
            # Pad rows for tiny matrices to satisfy vsn2's minimum data requirement.
            pad_target <- max(nrow(data_matrix) * 2L, 50L)
            for (i in seq_len(3L)) {
                idx <- rep(seq_len(nrow(data_matrix)), length.out = pad_target)
                padded_matrix <- data_matrix[idx, , drop = FALSE]
                if (!is.null(rownames(padded_matrix))) {
                    rownames(padded_matrix) <- make.unique(rownames(padded_matrix))
                }
                expr_set_pad <- ExpressionSet(
                    assayData = padded_matrix,
                    phenoData = covrts.annodf
                )
                res <- tryCatch(
                    pvcaBatchAssess(expr_set_pad, factors_for_PVCA, threshold = pca_threshold),
                    error = function(err) err
                )
                if (!inherits(res, "error")) {
                    return(res)
                }
                msg <- conditionMessage(res)
                if (!grepl("minDataPointsPerStratum", msg, fixed = TRUE) &&
                    !grepl("too small for reliable estimation", msg, fixed = TRUE)) {
                    stop(res)
                }
                pad_target <- pad_target * 2L
            }
            stop(e)
        }
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
calculate_PVCA.ProBatchFeatures <- function(data_matrix, pbf_name = NULL,
                                            sample_annotation = NULL,
                                            feature_id_col = "peptide_group_label",
                                            sample_id_col = "FullRunName",
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
calculate_PVCA <- function(data_matrix, ...) UseMethod("calculate_PVCA")

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
#' @param data_matrix Input object: matrix-like data or a `ProBatchFeatures` instance.
#' @param pbf_name Assay name(s) used when `data_matrix` is a `ProBatchFeatures`.
#' @param return_gridExtra Logical; return arranged grobs instead of a plot list.
#' @param plot_ncol Number of columns when arranging multiple assay plots.
#' @param stacked_bar logical; when `TRUE` and multiple `pbf_name` values are
#'   provided, combines all assays in a single stacked bar chart (supported for
#'   `ProBatchFeatures` inputs only).
#' @param stacked_plot_title optional character vector that annotates the stacked
#'   PVCA plot when `stacked_bar = TRUE`; elements are joined with newlines.
#' @param sort_stacked optional factor name used to order stacked bars when
#'   `stacked_bar = TRUE`; assays are ordered by the explained variance of this
#'   factor in descending order.
#' @param ... Additional arguments passed to lower-level methods.
#' @param base_size base size of the text in the plot
#' @param add_values logical; when `TRUE`, annotates each bar with its rounded weight.
#' @param path_to_save_results optional path to save the PVCA results data frame as CSV.
#'
#' @name plot_PVCA
#' @return \code{ggplot} object with the plot
#' @export
#'
#' @examples
#' data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#' matrix_test <- na.omit(example_proteome_matrix)[1:50, ]
#'
#' pvca_file <- tempfile("pvca", fileext = ".png")
#' pvca_plot <- plot_PVCA(
#'     matrix_test,
#'     example_sample_annotation,
#'     technical_factors = c("MS_batch", "digestion_batch"),
#'     biological_factors = c("Diet", "Sex", "Strain"),
#'     filename = pvca_file, # save to file, can be NULL
#'     width = 12, height = 8, units = "cm"
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
                              base_size = 15,
                              path_to_save_results = NULL,
                              ...) {
    dots <- list(...)
    add_values <- FALSE
    if ("add_values" %in% names(dots)) {
        add_values <- isTRUE(dots$add_values)
        dots$add_values <- NULL
    }

    pvca_res <- prepare_PVCA_df(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        technical_factors = technical_factors,
        biological_factors = biological_factors,
        fill_the_missing = fill_the_missing,
        pca_threshold = pca_threshold,
        variance_threshold = variance_threshold,
        path_to_save_results = path_to_save_results
    )

    plot_args <- list(
        df = pvca_res,
        colors_for_bars = colors_for_bars,
        filename = filename,
        width = width,
        height = height,
        units = units,
        plot_title = plot_title,
        theme = theme,
        base_size = base_size,
        add_values = add_values
    )

    gg <- do.call(plot_PVCA.df, plot_args)
    if (!is.null(gg$data) && "category" %in% names(gg$data)) {
        gg$data$category <- as.character(gg$data$category)
    }
    return(gg)
}

#' @rdname plot_PVCA
#' @method plot_PVCA ProBatchFeatures
#' @export
plot_PVCA.ProBatchFeatures <- function(data_matrix, pbf_name = NULL,
                                       sample_annotation = NULL,
                                       feature_id_col = "peptide_group_label",
                                       sample_id_col = "FullRunName",
                                       plot_title = NULL,
                                       return_gridExtra = FALSE,
                                       plot_ncol = NULL,
                                       stacked_bar = FALSE,
                                       stacked_plot_title = "Plot of weighted average proportion variance vs effects in PVCA",
                                       sort_stacked = NULL,
                                       category_order = NULL,
                                       path_to_save_results = NULL,
                                       ...) {
    object <- data_matrix
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
    shared_title <- prep$shared_title

    if (is.null(sample_annotation)) {
        sample_annotation <- as.data.frame(colData(object), stringsAsFactors = FALSE)
        rownames(sample_annotation) <- NULL
    }
    sample_ann_list <- split_arg(sample_annotation)
    default_sample_annotation <- as.data.frame(colData(object), stringsAsFactors = FALSE)
    rownames(default_sample_annotation) <- NULL

    stacked_requested <- isTRUE(stacked_bar)
    use_stacked <- stacked_requested && length(assays) >= 2L

    if (use_stacked) {
        prepare_dots <- dots
        drop_args <- c(
            "colors_for_bars", "filename", "width", "height", "units",
            "plot_title", "theme", "base_size", "add_values",
            "stacked_bar", "sort_stacked", "return_gridExtra", "plot_ncol",
            "path_to_save_results"
        )
        for (nm in drop_args) {
            prepare_dots[[nm]] <- NULL
        }

        pvca_df_list <- vector("list", length(assays))
        names(pvca_df_list) <- assays

        for (i in seq_along(assays)) {
            assay_nm <- assays[[i]]
            data_matrix <- pb_assay_matrix(object, assay_nm)
            sample_ann <- sample_ann_list[[i]]
            if (is.null(sample_ann)) {
                sample_ann <- default_sample_annotation
            }
            sample_ann <- as.data.frame(sample_ann, stringsAsFactors = FALSE)

            path_to_save_results_assay <- NULL
            if (!is.null(path_to_save_results)) {
                path_to_save_results_assay <- file.path(path_to_save_results, assay_nm)
            }

            call_args <- c(list(
                data_matrix = data_matrix,
                sample_annotation = sample_ann,
                feature_id_col = feature_id_col,
                sample_id_col = sample_id_col,
                path_to_save_results = path_to_save_results_assay
            ), prepare_dots)

            pvca_df_list[[i]] <- do.call(prepare_PVCA_df.default, call_args)
        }

        stacked_filename <- NULL
        if (!is.null(filename_list)) {
            idx <- which(!vapply(filename_list, is.null, logical(1)))
            if (length(idx)) {
                stacked_filename <- filename_list[[idx[1]]]
            }
        }

        width_arg <- if ("width" %in% names(dots)) dots$width else NA
        height_arg <- if ("height" %in% names(dots)) dots$height else NA
        units_arg <- if ("units" %in% names(dots)) dots$units else c("cm", "in", "mm")
        theme_arg <- if ("theme" %in% names(dots)) dots$theme else "classic"
        base_size_arg <- if ("base_size" %in% names(dots)) dots$base_size else 15
        colors_arg <- dots$colors_for_bars

        stacked_title <- .pb_collapsed_plot_title(stacked_plot_title)
        if (is.null(stacked_title)) {
            stacked_title <- .pb_first_plot_title(plot_title)
            if (is.null(stacked_title)) {
                stacked_title <- shared_title
            }
        }

        gg <- .pb_plot_pvca_stacked_bar(
            pvca_df_list = pvca_df_list,
            assays = assays,
            colors_for_bars = colors_arg,
            plot_title = stacked_title,
            theme = theme_arg,
            base_size = base_size_arg,
            filename = stacked_filename,
            width = width_arg,
            height = height_arg,
            units = units_arg,
            sort_label = sort_stacked,
            category_order = category_order
        )
        return(gg)
    }

    plot_list <- vector("list", length(assays))
    names(plot_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        data_matrix <- pb_assay_matrix(object, assay_nm)
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- default_sample_annotation
        }
        sample_ann <- as.data.frame(sample_ann, stringsAsFactors = FALSE)

        path_to_save_results_assay <- NULL
        if (!is.null(path_to_save_results)) {
            path_to_save_results_assay <- file.path(path_to_save_results, assay_nm)
        }

        call_args <- .pb_per_assay_dots(dots, filename_list, i)

        call_args <- c(list(
            data_matrix = data_matrix,
            sample_annotation = sample_ann,
            feature_id_col = feature_id_col,
            sample_id_col = sample_id_col,
            plot_title = titles[i],
            path_to_save_results = path_to_save_results_assay
        ), call_args)

        plot_list[[i]] <- do.call(plot_PVCA.default, call_args)
    }

    plot_list <- .pb_attach_shared_title(plot_list, shared_title)

    .pb_arrange_plot_list(plot_list, convert_fun = ggplotGrob, plot_ncol = plot_ncol, return_gridExtra = return_gridExtra)
}

.pb_plot_pvca_stacked_bar <- function(pvca_df_list,
                                      assays,
                                      colors_for_bars = NULL,
                                      plot_title = NULL,
                                      theme = "classic",
                                      base_size = 15,
                                      filename = NULL,
                                      width = NA,
                                      height = NA,
                                      units = c("cm", "in", "mm"),
                                      sort_label = NULL,
                                      category_order = c("biological", "biol:techn", "residual", "technical")) {
    if (!length(pvca_df_list)) {
        return(invisible(NULL))
    }

    combined_list <- list()
    label_levels <- character()

    for (i in seq_along(pvca_df_list)) {
        df <- pvca_df_list[[i]]
        if (is.null(df) || !nrow(df)) {
            next
        }
        assay_nm <- if (!is.null(names(pvca_df_list)) && nzchar(names(pvca_df_list)[i])) {
            names(pvca_df_list)[i]
        } else if (length(assays) >= i) {
            assays[[i]]
        } else {
            paste0("assay_", i)
        }
        df <- df %>%
            mutate(label = as.character(label), assay = assay_nm)

        levs <- levels(df$label)
        if (is.null(levs)) {
            levs <- unique(df$label)
        }
        label_levels <- c(label_levels, levs)
        combined_list[[length(combined_list) + 1L]] <- df
    }

    stacked_df <- bind_rows(combined_list)
    if (!nrow(stacked_df)) {
        return(invisible(NULL))
    }

    if (length(label_levels)) {
        label_levels <- unique(label_levels)
        label_levels <- c(label_levels, setdiff(unique(stacked_df$label), label_levels))
    } else {
        label_levels <- unique(stacked_df$label)
    }

    stacked_df <- stacked_df %>%
        mutate(label = factor(label, levels = label_levels))

    if (is.null(category_order)) {
        category_order <- unique(stacked_df$category)
        if ("residual" %in% category_order) {
            non_res <- setdiff(category_order, "residual")
            # if there are at least two non-residual categories, perform the rotation,
            # otherwise just put 'residual' at the end (handles length 0 or 1 safely)
            if (length(non_res) >= 2) {
                category_order <- c(non_res[-1], "residual", non_res[1])
            } else {
                category_order <- c(non_res, "residual")
            }
        }
    } else {
        category_order <- category_order[category_order %in% unique(stacked_df$category)]
        # if none of the supplied categories match, fall back to observed categories
        if (!length(category_order)) {
            category_order <- unique(stacked_df$category)
            if ("residual" %in% category_order) {
                category_order <- c(setdiff(category_order, "residual"), "residual")
            }
        }
    }
    stacked_df <- stacked_df %>%
        mutate(category = factor(category, levels = rev(category_order)))

    present_assays <- unique(stacked_df$assay)
    assay_levels <- assays[assays %in% present_assays]
    assay_levels <- unique(c(assay_levels, setdiff(present_assays, assay_levels)))

    if (!is.null(sort_label)) {
        sort_label_chr <- as.character(sort_label)[1]
        sort_weights <- stacked_df %>%
            filter(label == sort_label_chr) %>%
            group_by(assay) %>%
            summarise(weight = sum(weights, na.rm = TRUE), .groups = "drop") %>%
            arrange(desc(weight))
        if (nrow(sort_weights)) {
            assay_levels <- c(sort_weights$assay, setdiff(assay_levels, sort_weights$assay))
        }
    }

    display_levels <- rev(assay_levels)
    stacked_df <- stacked_df %>%
        mutate(assay = factor(assay, levels = display_levels))

    colors_vec <- colors_for_bars
    if (is.list(colors_vec) && !is.data.frame(colors_vec)) {
        colors_vec <- colors_vec[[1]]
    }
    if (is.null(colors_vec)) {
        colors_vec <- c("grey", wes_palettes$Rushmore[3:5])
        names(colors_vec) <- c(
            "residual", "biological",
            "biol:techn", "technical"
        )
    } else if (length(colors_vec) != 4) {
        color_names <- paste(c(
            "residual", "biological", "biol:techn",
            "technical"
        ), collapse = " ")
        warning(sprintf("four colors for: %s were expected", color_names))
    }

    totals <- stacked_df %>%
        group_by(assay) %>%
        summarise(total = sum(weights, na.rm = TRUE), .groups = "drop")
    max_weight <- max(totals$total, na.rm = TRUE)
    if (!is.finite(max_weight)) {
        max_weight <- 0
    }

    gg <- ggplot(stacked_df, aes(x = weights, y = assay, fill = category)) +
        geom_col(color = "black") +
        xlab("Weighted average proportion variance") +
        ylab(NULL)

    if (max_weight > 0) {
        gg <- gg + expand_limits(x = if (max_weight * 1.05 <= 1) max_weight * 1.05 else 1.05)
    }

    gg <- gg + scale_fill_manual(values = colors_vec, breaks = category_order, limits = category_order)

    if (!is.null(plot_title)) {
        gg <- gg + ggtitle(plot_title)
    }

    if (!is.null(theme) && theme == "classic") {
        gg <- gg + theme_classic(base_size = base_size)
    } else {
        message("plotting with default ggplot theme, only theme = 'classic'
            implemented")
    }

    gg <- gg +
        theme(
            axis.title.y = element_blank(),
            axis.text.y = element_text(hjust = 1),
            plot.title = element_text(size = round(base_size * 1.2, 0))
        ) +
        guides(fill = guide_legend(override.aes = list(color = NA), title = NULL))

    save_ggplot(filename, units, width, height, gg)

    gg
}

#' Plot stacked PVCA results from saved CSV files
#'
#' When PVCA analysis has been run previously and the aggregated weights were
#' written to disk (one `PVCA_results_aggregated.csv` per assay), this helper
#' rebuilds the stacked bar summary without rerunning the analysis.
#'
#' @param pvca_dir Directory that contains per-assay
#'   `PVCA_results_aggregated.csv` files (the function searches recursively).
#' @param sort_stacked Optional factor label used to order the rows of the
#'   stacked plot (matching `sort_stacked` in [plot_PVCA.ProBatchFeatures]).
#' @param colors_for_bars Character vector of four colors (residual, biological,
#'   biol:techn, technical) passed to the stacked bar helper.
#' @param plot_title Optional plot title.
#' @param stacked_plot_title Optional character vector used to annotate the stacked
#'   plot title; elements are joined with newline characters.
#' @param theme Plot theme; only `"classic"` is currently implemented.
#' @param base_size Base font size passed to `theme_classic()`.
#' @param filename Optional path to save the stacked plot.
#' @param width Plot width forwarded to [save_ggplot()].
#' @param height Plot height forwarded to [save_ggplot()].
#' @param units Plot units forwarded to [save_ggplot()].
#' @param category_order Optional character vector specifying the order of categories in the stacked plot.
#'
#' @return A `ggplot` object showing the stacked PVCA weights.
#' @export
plot_PVCA_stacked_from_saved <- function(pvca_dir,
                                         sort_stacked = NULL,
                                         colors_for_bars = NULL,
                                         plot_title = NULL,
                                         stacked_plot_title = "Plot of weighted average proportion variance vs effects in PVCA",
                                         theme = "classic",
                                         base_size = 15,
                                         filename = NULL,
                                         width = NA,
                                         height = NA,
                                         units = c("cm", "in", "mm"),
                                         category_order = NULL) {
    if (missing(pvca_dir) || is.null(pvca_dir) || !nzchar(pvca_dir)) {
        stop("`pvca_dir` must point to a directory with saved PVCA results.")
    }
    if (!dir.exists(pvca_dir)) {
        stop(sprintf("PVCA directory '%s' does not exist.", pvca_dir))
    }

    csv_paths <- dir(
        path = pvca_dir,
        pattern = "^PVCA_results_aggregated\\.csv$",
        recursive = TRUE,
        full.names = TRUE
    )
    if (!length(csv_paths)) {
        stop(sprintf("No PVCA_results_aggregated.csv files found under '%s'.", pvca_dir))
    }

    pvca_dfs <- lapply(csv_paths, function(path) {
        utils::read.csv(path, stringsAsFactors = FALSE)
    })
    has_rows <- vapply(pvca_dfs, function(df) {
        is.data.frame(df) && nrow(df) > 0
    }, logical(1))
    if (!any(has_rows)) {
        stop("Found PVCA files, but none contain any rows of data.")
    }

    csv_paths <- csv_paths[has_rows]
    pvca_dfs <- pvca_dfs[has_rows]
    assay_names <- basename(dirname(csv_paths))
    assay_names <- make.unique(assay_names, sep = "_")
    names(pvca_dfs) <- assay_names

    stacked_title <- .pb_collapsed_plot_title(stacked_plot_title)
    if (is.null(stacked_title)) {
        stacked_title <- .pb_collapsed_plot_title(plot_title)
    }
    units <- match.arg(units)
    gg <- .pb_plot_pvca_stacked_bar(
        pvca_df_list = pvca_dfs,
        assays = assay_names,
        colors_for_bars = colors_for_bars,
        plot_title = stacked_title,
        theme = theme,
        base_size = base_size,
        filename = filename,
        width = width,
        height = height,
        units = units,
        sort_label = sort_stacked,
        category_order = category_order
    )

    gg
}

#' @export
plot_PVCA <- function(data_matrix, ...) UseMethod("plot_PVCA")

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
#' @param pbf_name Assay name(s) used when `data_matrix` is a `ProBatchFeatures`.
#' @param path_to_save_results optional path to save the PVCA results data frame as CSV.
#' @param ... Additional arguments forwarded between methods.
#'
#' @return data frame with weights and factors, combined in a way ready for plotting
#' @export
#'
#' @examples
#' data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#' matrix_test <- na.omit(example_proteome_matrix)[1:50, ]
#'
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
                                    pca_threshold = .6,
                                    variance_threshold = .01,
                                    path_to_save_results = NULL,
                                    ...) {
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

    pvca_res <- pvca_res %>%
        mutate(label = as.character(label))

    tech_interactions <- expand.grid(technical_factors, technical_factors) %>%
        mutate(tech_interactions = paste(Var1, Var2, sep = ":")) %>%
        pull(tech_interactions)
    biol_interactions <- expand.grid(biological_factors, biological_factors) %>%
        mutate(biol_interactions = paste(Var1, Var2, sep = ":")) %>%
        pull(biol_interactions)

    label_of_small <- sprintf("Below %1.0f%%", 100 * variance_threshold)
    technical_factors <- c(technical_factors, tech_interactions)
    biological_factors <- c(biological_factors, biol_interactions)
    residual_labels <- c("Residuals", "Residual", "resid", "residual")
    if (!is.null(label_of_small)) {
        residual_labels <- c(residual_labels, label_of_small)
    }
    residual_labels <- unique(tolower(residual_labels))

    pvca_res <- pvca_res %>%
        mutate(category = case_when(
            label %in% technical_factors ~ "technical",
            label %in% biological_factors ~ "biological",
            tolower(label) %in% residual_labels ~ "residual",
            TRUE ~ "biol:techn"
        ))

    pvca_res <- pvca_res %>%
        arrange(desc(weights)) %>%
        arrange(label == label_of_small) %>%
        arrange(label == "resid")

    # if path_to_save_results is provided, save the pvca_res data frame there
    if (!is.null(path_to_save_results)) {
        if (!dir.exists(path_to_save_results)) {
            dir.create(path_to_save_results, recursive = TRUE)
        }
        pvca_res_file <- file.path(path_to_save_results, "PVCA_results_aggregated.csv")
        write.csv(pvca_res, pvca_res_file, row.names = FALSE)
    }

    return(pvca_res)
}

#' @rdname prepare_PVCA_df
#' @method prepare_PVCA_df ProBatchFeatures
#' @export
prepare_PVCA_df.ProBatchFeatures <- function(data_matrix, pbf_name = NULL,
                                             sample_annotation = NULL,
                                             feature_id_col = "peptide_group_label",
                                             sample_id_col = "FullRunName",
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
    dots <- prep$dots
    split_arg <- prep$split_arg

    default_sample_annotation <- as.data.frame(colData(object), stringsAsFactors = FALSE)
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
        sample_ann <- as.data.frame(sample_ann, stringsAsFactors = FALSE)

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
prepare_PVCA_df <- function(data_matrix, ...) UseMethod("prepare_PVCA_df")

#' plot PVCA, when the analysis is completed
#'
#' @inheritParams proBatch
#' @param df Data frame of PVCA weights, typically the result of `calculate_PVCA()`.
#' @param colors_for_bars four-item color vector, specifying colors for the
#'   following categories: c('residual', 'biological', 'biol:techn',
#'   'technical')
#' @param base_size base size of the text in the plot
#' @param pbf_name Assay name(s) used when `df` is a `ProBatchFeatures`.
#' @param return_gridExtra Logical; return arranged grobs instead of a plot list.
#' @param plot_ncol Number of columns when arranging multiple assay plots.
#' @param stacked_bar logical; when `TRUE` and multiple `pbf_name` entries are
#'   supplied for a `ProBatchFeatures` object, a single stacked bar chart is
#'   produced instead of subplots.
#' @param stacked_plot_title optional character vector providing lines for the
#'   stacked plot title; entries are joined with newlines before rendering.
#' @param sort_stacked optional factor name; when `stacked_bar = TRUE`, assays
#'   are ordered by the explained variance of this factor (descending).
#' @param ... Additional arguments forwarded to `prepare_PVCA_df()`.
#' @param add_values logical; when `TRUE`, annotates each bar with its rounded weight.
#'
#' @return \code{ggplot} object with bars as weights, colored by bio/tech factors
#' @export
#'
#' @examples
#' data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#' matrix_test <- na.omit(example_proteome_matrix)[1:50, ]
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
plot_PVCA.df.default <- function(df,
                                 colors_for_bars = NULL,
                                 filename = NULL, width = NA, height = NA,
                                 units = c("cm", "in", "mm"),
                                 plot_title = NULL,
                                 theme = "classic",
                                 base_size = 15, add_values = FALSE, ...) {
    pvca_res <- df
    pvca_res <- pvca_res %>%
        mutate(label = as.character(label))
    label_levels <- unique(pvca_res$label)
    pvca_res <- pvca_res %>%
        mutate(label = factor(label, levels = label_levels))
    label_category_map <- pvca_res %>%
        distinct(label, category) %>%
        arrange(label)
    category_levels <- unique(label_category_map$category)
    if (!length(category_levels)) {
        category_levels <- unique(pvca_res$category)
    }
    pvca_res <- pvca_res %>%
        mutate(category = factor(category, levels = category_levels))
    max_weight <- max(pvca_res$weights, na.rm = TRUE)
    if (!is.finite(max_weight)) {
        max_weight <- 0
    }

    gg <- ggplot(pvca_res, aes(x = label, y = weights, fill = category)) +
        geom_bar(stat = "identity", color = "black") +
        ylab("Weighted average proportion variance")
    if (max_weight > 0) {
        gg <- gg + expand_limits(y = if (max_weight * 1.05 <= 1) max_weight * 1.05 else 1.05)
    }

    default_cat_names <- c("residual", "biological", "biol:techn", "technical")
    if (is.null(colors_for_bars)) {
        colors_for_bars <- c("grey", wes_palettes$Rushmore[3:5])
        names(colors_for_bars) <- default_cat_names
    } else if (length(colors_for_bars) != 4) {
        color_names <- paste(default_cat_names, collapse = " ")
        warning(sprintf("four colors for: %s were expected", color_names))
    } else {
        color_names <- names(colors_for_bars)
        if (is.null(color_names) || any(!nzchar(color_names))) {
            names(colors_for_bars) <- default_cat_names
        }
    }
    color_breaks <- levels(pvca_res$category)
    gg <- gg + scale_fill_manual(values = colors_for_bars, breaks = color_breaks, limits = color_breaks)

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

    if (isTRUE(add_values)) {
        gg <- gg +
            geom_text(aes(label = sprintf("%.2f", weights)),
                vjust = -0.3,
                size = base_size / 3
            )
    }

    save_ggplot(filename, units, width, height, gg)

    return(gg)
}

#' @rdname plot_PVCA.df
#' @method plot_PVCA.df ProBatchFeatures
#' @export
plot_PVCA.df.ProBatchFeatures <- function(df, pbf_name = NULL,
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
                                          stacked_bar = FALSE,
                                          stacked_plot_title = "Plot of weighted average proportion variance vs effects in PVCA",
                                          sort_stacked = NULL,
                                          category_order = NULL,
                                          ...) {
    object <- df
    dots <- list(...)
    if (!"filename" %in% names(dots)) {
        dots <- c(list(filename = filename), dots)
    }

    prep <- .pb_prepare_multi_assay(
        object = object,
        pbf_name = pbf_name,
        dots = dots,
        plot_title = plot_title,
        default_title_fun = function(x) x
    )
    assays <- prep$assays
    prepare_dots <- prep$dots
    filename_list <- prep$filename_list
    split_arg <- prep$split_arg
    titles <- prep$titles
    shared_title <- prep$shared_title
    add_values_arg <- NULL
    if ("add_values" %in% names(prepare_dots)) {
        add_values_arg <- prepare_dots$add_values
        prepare_dots$add_values <- NULL
    }
    add_values_list <- split_arg(add_values_arg)

    default_sample_annotation <- as.data.frame(colData(object), stringsAsFactors = FALSE)
    rownames(default_sample_annotation) <- NULL

    if (is.null(sample_annotation)) {
        sample_annotation <- default_sample_annotation
    }
    sample_ann_list <- split_arg(sample_annotation)

    if (is.list(colors_for_bars) && !is.data.frame(colors_for_bars)) {
        colors_for_bars_list <- split_arg(colors_for_bars)
    } else {
        colors_for_bars_list <- rep(list(colors_for_bars), length(assays))
    }

    pvca_res_list <- vector("list", length(assays))
    names(pvca_res_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        data_matrix <- pb_assay_matrix(object, assay_nm)
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- default_sample_annotation
        }
        sample_ann <- as.data.frame(sample_ann, stringsAsFactors = FALSE)

        prepare_args <- c(list(
            data_matrix = data_matrix,
            sample_annotation = sample_ann,
            feature_id_col = feature_id_col,
            sample_id_col = sample_id_col
        ), prepare_dots)
        pvca_res <- do.call(prepare_PVCA_df.default, prepare_args)
        pvca_res_list[[i]] <- pvca_res
    }

    stacked_requested <- isTRUE(stacked_bar)
    use_stacked <- stacked_requested && length(assays) >= 2L
    if (use_stacked) {
        stacked_filename <- filename
        if (!is.null(filename_list)) {
            idx <- which(!vapply(filename_list, is.null, logical(1)))
            if (length(idx)) {
                stacked_filename <- filename_list[[idx[1]]]
            }
        }
        stacked_title <- .pb_collapsed_plot_title(stacked_plot_title)
        if (is.null(stacked_title)) {
            stacked_title <- .pb_first_plot_title(plot_title)
            if (is.null(stacked_title)) {
                stacked_title <- shared_title
            }
        }
        stacked_colors <- NULL
        color_idx <- which(!vapply(colors_for_bars_list, is.null, logical(1)))
        if (length(color_idx)) {
            stacked_colors <- colors_for_bars_list[[color_idx[1]]]
        }

        gg <- .pb_plot_pvca_stacked_bar(
            pvca_df_list = pvca_res_list,
            assays = assays,
            colors_for_bars = stacked_colors,
            plot_title = stacked_title,
            theme = theme,
            base_size = base_size,
            filename = stacked_filename,
            width = width,
            height = height,
            units = units,
            sort_label = sort_stacked,
            category_order = category_order
        )
        return(gg)
    }

    plot_list <- vector("list", length(assays))
    names(plot_list) <- assays

    for (i in seq_along(assays)) {
        pvca_res <- pvca_res_list[[i]]
        fn <- filename
        if (!is.null(filename_list)) {
            fn <- filename_list[[i]]
        }

        plot_args <- list(
            df = pvca_res,
            colors_for_bars = colors_for_bars_list[[i]],
            filename = fn,
            width = width,
            height = height,
            units = units,
            plot_title = titles[i],
            theme = theme,
            base_size = base_size
        )
        add_values_val <- add_values_list[[i]]
        if (!is.null(add_values_val)) {
            plot_args$add_values <- isTRUE(add_values_val)
        }

        plot_list[[i]] <- do.call(plot_PVCA.df.default, plot_args)
    }

    plot_list <- .pb_attach_shared_title(plot_list, shared_title)

    .pb_arrange_plot_list(plot_list, convert_fun = ggplotGrob, plot_ncol = plot_ncol, return_gridExtra = return_gridExtra)
}

#' @export
plot_PVCA.df <- function(df, ...) UseMethod("plot_PVCA.df")

#' Calculate per-feature variance partition contributions
#'
#' @inheritParams proBatch
#' @param model_formula model formula passed to
#'   `variancePartition::fitExtractVarPartModel()`.
#' @param model_variables vector of \code{sample_annotation} column names used to
#'   build a formula when `model_formula` is `NULL`.
#' @param fill_the_missing numeric value determining how missing values should be
#'   substituted. If \code{NULL}, features with missing values are excluded.
#' @param ... Additional arguments forwarded to
#'   `variancePartition::fitExtractVarPartModel()` (e.g., `BPPARAM`).
#'
#' @return data frame with columns \code{feature_id}, \code{label} and
#'   \code{variance_explained}.
#' @name calculate_variance_partition
#' @export
#' @importFrom variancePartition fitExtractVarPartModel
#'
#' @examples
#' if (requireNamespace("variancePartition", quietly = TRUE)) {
#'     data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#'     matrix_test <- na.omit(example_proteome_matrix)[1:30, ]
#'     vp_res <- calculate_variance_partition(
#'         matrix_test,
#'         example_sample_annotation,
#'         model_formula = ~ Diet + Sex
#'     )
#'     # Additional arguments (e.g., BPPARAM = BiocParallel::SnowParam(4)) can be
#'     # supplied via `...` and will be forwarded to variancePartition.
#' }
calculate_variance_partition.default <- function(data_matrix, sample_annotation,
                                                 feature_id_col = "peptide_group_label",
                                                 sample_id_col = "FullRunName",
                                                 model_formula = NULL,
                                                 model_variables = NULL,
                                                 fill_the_missing = -1,
                                                 ...) {
    fit_args <- list(...)

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

    if (is.null(model_formula)) {
        if (is.null(model_variables)) {
            model_variables <- c(
                "MS_batch", "digestion_batch", "Diet",
                "Sex", "Strain"
            )
        }
        model_variables <- model_variables[
            model_variables %in% names(sample_annotation)
        ]
        if (length(model_variables)) {
            model_levels <- vapply(sample_annotation[model_variables], function(x) {
                length(unique(x[!is.na(x)]))
            }, integer(1))
            model_variables <- model_variables[model_levels > 1]
        }
        if (!length(model_variables)) {
            stop("No valid variables supplied via 'model_formula' or 'model_variables'.")
        }
        model_formula <- reformulate(model_variables)
    } else {
        model_formula <- as.formula(model_formula)
        vars_in_model <- all.vars(model_formula)
        missing_vars <- setdiff(vars_in_model, names(sample_annotation))
        if (length(missing_vars) > 0) {
            stop(sprintf(
                "Variables missing from sample_annotation: %s",
                paste(missing_vars, collapse = ", ")
            ))
        }
    }

    vars_in_model <- all.vars(model_formula)
    sample_annotation <- sample_annotation %>%
        select(all_of(unique(c(sample_id_col, vars_in_model)))) %>%
        mutate(across(where(is.POSIXct), as.numeric)) %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        remove_rownames() %>%
        column_to_rownames(var = sample_id_col)

    data_matrix <- check_feature_id_col_in_dm(feature_id_col, data_matrix)
    data_matrix <- .pb_handle_missing_wrapper(
        data_matrix = data_matrix,
        warning_message = "variancePartition cannot operate with missing values in the matrix",
        fill_the_missing = fill_the_missing
    )
    data_matrix <- as.matrix(data_matrix)
    data_matrix <- data_matrix[, rownames(sample_annotation), drop = FALSE]

    fit_call <- c(list(
        exprObj = data_matrix,
        formula = model_formula,
        data = sample_annotation
    ), fit_args)

    var_part <- do.call(variancePartition::fitExtractVarPartModel, fit_call)

    var_part_df <- as.data.frame(var_part, stringsAsFactors = FALSE)
    if (!nrow(var_part_df)) {
        return(tibble(
            feature_id = character(),
            label = character(),
            variance_explained = numeric()
        ))
    }
    label_levels <- colnames(var_part_df)
    var_part_df <- var_part_df %>%
        rownames_to_column(var = "feature_id") %>%
        pivot_longer(
            cols = -feature_id,
            names_to = "label",
            values_to = "variance_explained"
        ) %>%
        mutate(label = factor(label, levels = label_levels))

    return(var_part_df)
}

#' @rdname calculate_variance_partition
#' @method calculate_variance_partition ProBatchFeatures
#' @export
calculate_variance_partition.ProBatchFeatures <- function(data_matrix, pbf_name = NULL,
                                                          sample_annotation = NULL,
                                                          feature_id_col = "peptide_group_label",
                                                          sample_id_col = "FullRunName",
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
    dots <- prep$dots
    split_arg <- prep$split_arg

    default_sample_annotation <- as.data.frame(colData(object))
    sample_ann_list <- split_arg(sample_annotation)

    vp_list <- vector("list", length(assays))
    names(vp_list) <- assays

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

        vp_list[[i]] <- do.call(calculate_variance_partition.default, call_args)
    }

    if (length(vp_list) == 1L) {
        return(vp_list[[1L]])
    }

    vp_list
}

#' @export
calculate_variance_partition <- function(data_matrix, ...) UseMethod("calculate_variance_partition")

#' Prepare variance partition results for plotting
#'
#' @inheritParams proBatch
#' @inheritParams calculate_variance_partition
#' @param technical_factors vector with technically driven covariates.
#' @param biological_factors vector with biologically driven covariates.
#' @param variance_threshold numeric threshold; values below the threshold are
#'   summed per feature and reported as a single group.
#' @param path_to_save_results optional path to save the variance partition
#'   results as CSV.
#' @param ... Additional arguments forwarded to `calculate_variance_partition()`
#'   and ultimately to `variancePartition::fitExtractVarPartModel()` (e.g.,
#'   `BPPARAM`).
#'
#' @return data frame ready for plotting variance partition distributions.
#' @export
#'
#' @examples
#' if (requireNamespace("variancePartition", quietly = TRUE)) {
#'     data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#'     matrix_test <- na.omit(example_proteome_matrix)[1:30, ]
#'     vp_df <- prepare_variance_partition_df(
#'         matrix_test,
#'         example_sample_annotation,
#'         technical_factors = "MS_batch",
#'         biological_factors = c("Diet", "Sex"),
#'         variance_threshold = 0.05
#'     )
#' }
#' @name prepare_variance_partition_df
prepare_variance_partition_df.default <- function(data_matrix, sample_annotation,
                                                  feature_id_col = "peptide_group_label",
                                                  sample_id_col = "FullRunName",
                                                  technical_factors = c("MS_batch", "instrument"),
                                                  biological_factors = c("cell_line", "drug_dose"),
                                                  fill_the_missing = -1,
                                                  model_formula = NULL,
                                                  model_variables = NULL,
                                                  variance_threshold = NULL,
                                                  path_to_save_results = NULL,
                                                  ...) {
    if (is.null(model_formula)) {
        if (is.null(model_variables)) {
            model_variables <- c(technical_factors, biological_factors)
        }
    } else {
        model_variables <- NULL
    }

    vp_res <- calculate_variance_partition(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        model_formula = model_formula,
        model_variables = model_variables,
        fill_the_missing = fill_the_missing,
        ...
    )

    if (!nrow(vp_res)) {
        return(vp_res)
    }

    label_levels <- levels(vp_res$label)
    if (is.null(label_levels)) {
        label_levels <- unique(vp_res$label)
    }
    vp_res <- vp_res %>% mutate(label = as.character(label))

    label_of_small <- NULL
    if (!is.null(variance_threshold)) {
        label_of_small <- sprintf("Below %1.0f%%", 100 * variance_threshold)
        small_res <- vp_res %>%
            filter(variance_explained < variance_threshold) %>%
            group_by(feature_id) %>%
            summarise(variance_explained = sum(variance_explained, na.rm = TRUE), .groups = "drop") %>%
            filter(variance_explained > 0) %>%
            mutate(label = label_of_small)

        vp_res <- vp_res %>%
            filter(variance_explained >= variance_threshold)

        if (nrow(small_res) > 0) {
            vp_res <- dplyr::bind_rows(vp_res, small_res)
        }
    }

    if (!is.null(label_of_small) && !label_of_small %in% label_levels) {
        label_levels <- c(label_levels, label_of_small)
    }

    residual_labels <- c("Residuals", "Residual", "resid", "residual")
    if (!is.null(label_of_small)) {
        residual_labels <- c(residual_labels, label_of_small)
    }
    residual_labels <- unique(tolower(residual_labels))

    vp_res <- vp_res %>%
        mutate(category = ifelse(label %in% technical_factors,
            "technical",
            ifelse(label %in% biological_factors,
                "biological",
                ifelse(
                    tolower(label) %in% residual_labels |
                        grepl("^below\\s+\\d+%$", tolower(label)),
                    "residual", "biol:techn"
                )
            )
        ))

    present_labels <- unique(vp_res$label)
    final_levels <- c(label_levels[label_levels %in% present_labels], setdiff(present_labels, label_levels))
    vp_res <- vp_res %>%
        mutate(label = factor(label, levels = final_levels))

    if (!is.null(path_to_save_results)) {
        if (!dir.exists(path_to_save_results)) {
            dir.create(path_to_save_results, recursive = TRUE)
        }
        vp_res_file <- file.path(path_to_save_results, "variance_partition_results_long.csv")
        write.csv(vp_res, vp_res_file, row.names = FALSE)
    }

    vp_res
}

#' @rdname prepare_variance_partition_df
#' @method prepare_variance_partition_df ProBatchFeatures
#' @export
prepare_variance_partition_df.ProBatchFeatures <- function(data_matrix, pbf_name = NULL,
                                                           sample_annotation = NULL,
                                                           feature_id_col = "peptide_group_label",
                                                           sample_id_col = "FullRunName",
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
    dots <- prep$dots
    split_arg <- prep$split_arg

    default_sample_annotation <- as.data.frame(colData(object), stringsAsFactors = FALSE)
    sample_ann_list <- split_arg(sample_annotation)

    vp_df_list <- vector("list", length(assays))
    names(vp_df_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        data_matrix <- pb_assay_matrix(object, assay_nm)
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- default_sample_annotation
        }
        sample_ann <- as.data.frame(sample_ann, stringsAsFactors = FALSE)

        call_args <- c(list(
            data_matrix = data_matrix,
            sample_annotation = sample_ann,
            feature_id_col = feature_id_col,
            sample_id_col = sample_id_col
        ), dots)

        vp_df_list[[i]] <- do.call(prepare_variance_partition_df.default, call_args)
    }

    if (length(vp_df_list) == 1L) {
        return(vp_df_list[[1L]])
    }

    vp_df_list
}

#' @export
prepare_variance_partition_df <- function(data_matrix, ...) UseMethod("prepare_variance_partition_df")

#' Plot variance partition distributions
#'
#' @inheritParams proBatch
#' @inheritParams prepare_variance_partition_df
#' @inheritParams prepare_variance_partition_df
#' @param colors_for_boxes four-item color vector, specifying colors for the
#'   following categories: c('residual', 'biological', 'biol:techn', 'technical').
#' @param plot_title optional plot title.
#' @param theme ggplot2 theme to apply. Currently, only `"classic"` is supported.
#' @param base_size base size of the text in the plot.
#' @param add_medians logical; when `TRUE`, annotates each box with its rounded
#'   median value.
#' @param median_position numeric adjustment applied to the median labels when
#'   they are nudged upwards.
#' @param show_legend logical; controls legend display.
#' @param y_limits numeric vector of length two specifying y-axis limits.
#' @param return_gridExtra logical; when `TRUE`, returns arranged grobs instead of
#'   a list of plots (for multi-assay objects).
#' @param plot_ncol number of columns when arranging multiple assay plots.
#' @param organize_pbfs logical; when `TRUE` (multi-assay inputs only), derives and
#'   applies shared y-axis limits across all subplots.
#' @param path_to_save_results optional path to save the prepared results as CSV.
#' @param ... Additional arguments forwarded to `prepare_variance_partition_df()`
#'   (and further to `variancePartition::fitExtractVarPartModel()`).
#'
#' @name plot_variance_partition
#' @return ggplot object (or list of ggplot objects for multi-assay inputs).
#' @export
#'
#' @examples
#' if (requireNamespace("variancePartition", quietly = TRUE)) {
#'     data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")
#'     matrix_test <- na.omit(example_proteome_matrix)[1:30, ]
#'     plot_variance_partition(
#'         matrix_test,
#'         example_sample_annotation,
#'         technical_factors = "MS_batch",
#'         biological_factors = c("Diet", "Sex"),
#'         variance_threshold = 0.05
#'     )
#' }
plot_variance_partition.default <- function(data_matrix, sample_annotation,
                                            feature_id_col = "peptide_group_label",
                                            sample_id_col = "FullRunName",
                                            technical_factors = c("MS_batch", "instrument"),
                                            biological_factors = c("cell_line", "drug_dose"),
                                            fill_the_missing = -1,
                                            model_formula = NULL,
                                            model_variables = NULL,
                                            variance_threshold = NULL,
                                            colors_for_boxes = NULL,
                                            filename = NULL, width = NA, height = NA,
                                            units = c("cm", "in", "mm"),
                                            plot_title = NULL,
                                            theme = "classic",
                                            base_size = 15,
                                            add_medians = FALSE,
                                            median_position = 0.05,
                                            show_legend = TRUE,
                                            y_limits = c(0, 1),
                                            path_to_save_results = NULL,
                                            ...) {
    vp_res <- prepare_variance_partition_df(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        technical_factors = technical_factors,
        biological_factors = biological_factors,
        fill_the_missing = fill_the_missing,
        model_formula = model_formula,
        model_variables = model_variables,
        variance_threshold = variance_threshold,
        path_to_save_results = path_to_save_results,
        ...
    )

    plot_args <- list(
        df = vp_res,
        colors_for_boxes = colors_for_boxes,
        filename = filename,
        width = width,
        height = height,
        units = units,
        plot_title = plot_title,
        theme = theme,
        base_size = base_size,
        add_medians = add_medians,
        median_position = median_position,
        show_legend = show_legend,
        y_limits = y_limits
    )

    gg <- do.call(plot_variance_partition.df, plot_args)

    return(gg)
}

#' @rdname plot_variance_partition
#' @method plot_variance_partition ProBatchFeatures
#' @export
plot_variance_partition.ProBatchFeatures <- function(data_matrix, pbf_name = NULL,
                                                     sample_annotation = NULL,
                                                     feature_id_col = "peptide_group_label",
                                                     sample_id_col = "FullRunName",
                                                     colors_for_boxes = NULL,
                                                     filename = NULL, width = NA, height = NA,
                                                     units = c("cm", "in", "mm"),
                                                     plot_title = NULL,
                                                     theme = "classic",
                                                     base_size = 20,
                                                     return_gridExtra = FALSE,
                                                     plot_ncol = NULL,
                                                     organize_pbfs = FALSE,
                                                     path_to_save_results = NULL,
                                                     ...) {
    object <- data_matrix
    dots <- list(...)
    if (!"filename" %in% names(dots)) {
        dots <- c(list(filename = filename), dots)
    }

    prep <- .pb_prepare_multi_assay(
        object = object,
        pbf_name = pbf_name,
        dots = dots,
        plot_title = plot_title,
        default_title_fun = function(x) x
    )
    assays <- prep$assays
    prepare_dots <- prep$dots
    filename_list <- prep$filename_list
    split_arg <- prep$split_arg
    titles <- prep$titles
    shared_title <- prep$shared_title
    add_medians_arg <- NULL
    if ("add_medians" %in% names(prepare_dots)) {
        add_medians_arg <- prepare_dots$add_medians
        prepare_dots$add_medians <- NULL
    }
    add_medians_list <- split_arg(add_medians_arg)

    median_position_arg <- NULL
    if ("median_position" %in% names(prepare_dots)) {
        median_position_arg <- prepare_dots$median_position
        prepare_dots$median_position <- NULL
    }
    median_position_list <- split_arg(median_position_arg)

    show_legend_arg <- NULL
    if ("show_legend" %in% names(prepare_dots)) {
        show_legend_arg <- prepare_dots$show_legend
        prepare_dots$show_legend <- NULL
    }
    show_legend_list <- split_arg(show_legend_arg)

    y_limits_arg <- NULL
    y_limits_provided <- "y_limits" %in% names(prepare_dots)
    if (y_limits_provided) {
        y_limits_arg <- prepare_dots$y_limits
        prepare_dots$y_limits <- NULL
    }
    y_limits_list <- split_arg(y_limits_arg)

    default_sample_annotation <- as.data.frame(colData(object), stringsAsFactors = FALSE)
    rownames(default_sample_annotation) <- NULL

    if (is.null(sample_annotation)) {
        sample_annotation <- default_sample_annotation
    }
    sample_ann_list <- split_arg(sample_annotation)

    if (is.list(colors_for_boxes) && !is.data.frame(colors_for_boxes)) {
        colors_for_boxes_list <- split_arg(colors_for_boxes)
    } else {
        colors_for_boxes_list <- rep(list(colors_for_boxes), length(assays))
    }

    plot_list <- vector("list", length(assays))
    names(plot_list) <- assays

    vp_df_list <- vector("list", length(assays))
    names(vp_df_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        data_matrix <- pb_assay_matrix(object, assay_nm)
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- default_sample_annotation
        }
        sample_ann <- as.data.frame(sample_ann, stringsAsFactors = FALSE)

        path_to_save_results_assay <- NULL
        if (!is.null(path_to_save_results)) {
            path_to_save_results_assay <- file.path(path_to_save_results, assay_nm)
        }

        prepare_args <- c(list(
            data_matrix = data_matrix,
            sample_annotation = sample_ann,
            feature_id_col = feature_id_col,
            sample_id_col = sample_id_col,
            path_to_save_results = path_to_save_results_assay
        ), prepare_dots)
        vp_df_list[[i]] <- do.call(prepare_variance_partition_df.default, prepare_args)
    }

    shared_y_limits <- NULL
    if (length(assays) > 1L && isTRUE(organize_pbfs)) {
        shared_y_limits <- .pb_shared_variance_limits(vp_df_list)
    }

    for (i in seq_along(assays)) {
        vp_res <- vp_df_list[[i]]
        fn <- filename
        if (!is.null(filename_list)) {
            fn <- filename_list[[i]]
        }

        plot_args <- list(
            df = vp_res,
            colors_for_boxes = colors_for_boxes_list[[i]],
            filename = fn,
            width = width,
            height = height,
            units = units,
            plot_title = titles[i],
            theme = theme,
            base_size = base_size
        )
        add_medians_val <- add_medians_list[[i]]
        if (!is.null(add_medians_val)) {
            plot_args$add_medians <- isTRUE(add_medians_val)
        }
        median_position_val <- median_position_list[[i]]
        if (!is.null(median_position_val)) {
            plot_args$median_position <- median_position_val
        }
        show_legend_val <- show_legend_list[[i]]
        if (!is.null(show_legend_val)) {
            plot_args$show_legend <- isTRUE(show_legend_val)
        }
        y_limits_val <- y_limits_list[[i]]
        if (y_limits_provided) {
            plot_args["y_limits"] <- list(y_limits_val)
        } else if (!is.null(shared_y_limits)) {
            plot_args["y_limits"] <- list(shared_y_limits)
        }

        plot_list[[i]] <- do.call(plot_variance_partition.df.default, plot_args)
    }

    plot_list <- .pb_attach_shared_title(plot_list, shared_title)

    .pb_arrange_plot_list(plot_list, convert_fun = ggplotGrob, plot_ncol = plot_ncol, return_gridExtra = return_gridExtra)
}

#' @export
plot_variance_partition <- function(data_matrix, ...) UseMethod("plot_variance_partition")

#' Plot variance partition results stored as a data frame
#'
#' @inheritParams plot_variance_partition
#' @param df data frame produced by `prepare_variance_partition_df()` or a
#'   compatible structure.
#' @param ... Additional arguments ignored.
#'
#' @name plot_variance_partition.df
#' @return ggplot object
#' @export
plot_variance_partition.df.default <- function(df,
                                               colors_for_boxes = NULL,
                                               filename = NULL, width = NA, height = NA,
                                               units = c("cm", "in", "mm"),
                                               plot_title = NULL,
                                               theme = "classic",
                                               base_size = 15,
                                               add_medians = FALSE,
                                               median_position = 0.05,
                                               show_legend = TRUE,
                                               y_limits = c(0, 1),
                                               ...) {
    units <- match.arg(units)
    vp_res <- df
    if (!"variance_explained" %in% names(vp_res)) {
        stop("Input data frame must contain a 'variance_explained' column.")
    }
    label_levels <- levels(vp_res$label)
    if (is.null(label_levels)) {
        label_levels <- unique(vp_res$label)
    }
    vp_res <- vp_res %>%
        mutate(label = factor(label, levels = label_levels))

    label_category_map <- vp_res %>%
        distinct(label, category) %>%
        arrange(label)
    category_levels <- unique(label_category_map$category)
    if (!length(category_levels)) {
        category_levels <- unique(vp_res$category)
    }
    vp_res <- vp_res %>%
        mutate(category = factor(category, levels = category_levels))

    gg <- ggplot(vp_res, aes(x = label, y = variance_explained, fill = category)) +
        geom_boxplot()

    default_cat_names <- c("residual", "biological", "biol:techn", "technical")
    if (is.null(colors_for_boxes)) {
        colors_for_boxes <- c("grey", wes_palettes$Rushmore[3:5])
        names(colors_for_boxes) <- default_cat_names
    } else if (length(colors_for_boxes) != 4) {
        color_names <- paste(default_cat_names, collapse = " ")
        warning(sprintf("four colors for: %s were expected", color_names))
    } else {
        color_names <- names(colors_for_boxes)
        if (is.null(color_names) || any(!nzchar(color_names))) {
            names(colors_for_boxes) <- default_cat_names
        }
    }
    color_breaks <- levels(vp_res$category)
    gg <- gg + scale_fill_manual(values = colors_for_boxes, breaks = color_breaks, limits = color_breaks)

    if (!is.null(plot_title)) {
        gg <- gg + ggtitle(plot_title)
    }

    if (!is.null(y_limits)) {
        gg <- gg + coord_cartesian(ylim = y_limits)
        gg$coordinates$ylim <- y_limits
    }

    if (!is.null(theme) && theme == "classic") {
        gg <- gg + theme_classic(base_size = base_size)
    } else {
        message("plotting with default ggplot theme, only theme = 'classic' implemented")
    }

    legend_position <- "right"
    if (!isTRUE(show_legend)) {
        legend_position <- "none"
    }

    gg <- gg +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
            plot.title = element_text(size = round(base_size * 1.2, 0)),
            legend.position = legend_position
        ) +
        ylab("Variance explained") +
        xlab(NULL)

    if (isTRUE(add_medians)) {
        medians <- vp_res %>%
            group_by(label) %>%
            summarise(median_value = median(variance_explained, na.rm = TRUE), .groups = "drop")
        upper_limit <- 1
        if (!is.null(y_limits) && length(y_limits) == 2 && all(is.finite(y_limits))) {
            upper_limit <- max(y_limits)
        }
        gg <- gg +
            geom_text(
                data = medians,
                aes(
                    x = label,
                    y = median_value,
                    label = sprintf("%.2f", median_value)
                ),
                inherit.aes = FALSE,
                nudge_y = ifelse(medians$median_value > (upper_limit - median_position),
                    -median_position,
                    median_position
                ),
                size = base_size / 3,
                color = "black"
            )
    }

    save_ggplot(filename, units, width, height, gg)

    gg
}

#' @rdname plot_variance_partition.df
#' @method plot_variance_partition.df ProBatchFeatures
#' @export
plot_variance_partition.df.ProBatchFeatures <- function(df, pbf_name = NULL,
                                                        sample_annotation = NULL,
                                                        feature_id_col = "peptide_group_label",
                                                        sample_id_col = "FullRunName",
                                                        colors_for_boxes = NULL,
                                                        filename = NULL, width = NA, height = NA,
                                                        units = c("cm", "in", "mm"),
                                                        plot_title = NULL,
                                                        theme = "classic",
                                                        base_size = 20,
                                                        return_gridExtra = FALSE,
                                                        plot_ncol = NULL,
                                                        organize_pbfs = FALSE,
                                                        ...) {
    object <- df
    dots <- list(...)
    if (!"filename" %in% names(dots)) {
        dots <- c(list(filename = filename), dots)
    }

    prep <- .pb_prepare_multi_assay(
        object = object,
        pbf_name = pbf_name,
        dots = dots,
        plot_title = plot_title,
        default_title_fun = function(x) x
    )
    assays <- prep$assays
    prepare_dots <- prep$dots
    filename_list <- prep$filename_list
    split_arg <- prep$split_arg
    titles <- prep$titles
    shared_title <- prep$shared_title
    add_medians_arg <- NULL
    if ("add_medians" %in% names(prepare_dots)) {
        add_medians_arg <- prepare_dots$add_medians
        prepare_dots$add_medians <- NULL
    }
    add_medians_list <- split_arg(add_medians_arg)

    median_position_arg <- NULL
    if ("median_position" %in% names(prepare_dots)) {
        median_position_arg <- prepare_dots$median_position
        prepare_dots$median_position <- NULL
    }
    median_position_list <- split_arg(median_position_arg)

    show_legend_arg <- NULL
    if ("show_legend" %in% names(prepare_dots)) {
        show_legend_arg <- prepare_dots$show_legend
        prepare_dots$show_legend <- NULL
    }
    show_legend_list <- split_arg(show_legend_arg)

    y_limits_arg <- NULL
    y_limits_provided <- "y_limits" %in% names(prepare_dots)
    if (y_limits_provided) {
        y_limits_arg <- prepare_dots$y_limits
        prepare_dots$y_limits <- NULL
    }
    y_limits_list <- split_arg(y_limits_arg)

    default_sample_annotation <- as.data.frame(colData(object), stringsAsFactors = FALSE)
    rownames(default_sample_annotation) <- NULL

    if (is.null(sample_annotation)) {
        sample_annotation <- default_sample_annotation
    }
    sample_ann_list <- split_arg(sample_annotation)

    if (is.list(colors_for_boxes) && !is.data.frame(colors_for_boxes)) {
        colors_for_boxes_list <- split_arg(colors_for_boxes)
    } else {
        colors_for_boxes_list <- rep(list(colors_for_boxes), length(assays))
    }

    plot_list <- vector("list", length(assays))
    names(plot_list) <- assays

    vp_df_list <- vector("list", length(assays))
    names(vp_df_list) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        data_matrix <- pb_assay_matrix(object, assay_nm)
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- default_sample_annotation
        }
        sample_ann <- as.data.frame(sample_ann, stringsAsFactors = FALSE)

        prepare_args <- c(list(
            data_matrix = data_matrix,
            sample_annotation = sample_ann,
            feature_id_col = feature_id_col,
            sample_id_col = sample_id_col
        ), prepare_dots)
        vp_df_list[[i]] <- do.call(prepare_variance_partition_df.default, prepare_args)
    }

    shared_y_limits <- NULL
    if (length(assays) > 1L && isTRUE(organize_pbfs)) {
        shared_y_limits <- .pb_shared_variance_limits(vp_df_list)
    }

    for (i in seq_along(assays)) {
        vp_res <- vp_df_list[[i]]
        fn <- filename
        if (!is.null(filename_list)) {
            fn <- filename_list[[i]]
        }

        plot_args <- list(
            df = vp_res,
            colors_for_boxes = colors_for_boxes_list[[i]],
            filename = fn,
            width = width,
            height = height,
            units = units,
            plot_title = titles[i],
            theme = theme,
            base_size = base_size
        )
        add_medians_val <- add_medians_list[[i]]
        if (!is.null(add_medians_val)) {
            plot_args$add_medians <- isTRUE(add_medians_val)
        }
        median_position_val <- median_position_list[[i]]
        if (!is.null(median_position_val)) {
            plot_args$median_position <- median_position_val
        }
        show_legend_val <- show_legend_list[[i]]
        if (!is.null(show_legend_val)) {
            plot_args$show_legend <- isTRUE(show_legend_val)
        }
        y_limits_val <- y_limits_list[[i]]
        if (y_limits_provided) {
            plot_args["y_limits"] <- list(y_limits_val)
        } else if (!is.null(shared_y_limits)) {
            plot_args["y_limits"] <- list(shared_y_limits)
        }

        plot_list[[i]] <- do.call(plot_variance_partition.df.default, plot_args)
    }

    plot_list <- .pb_attach_shared_title(plot_list, shared_title)

    .pb_arrange_plot_list(plot_list, convert_fun = ggplotGrob, plot_ncol = plot_ncol, return_gridExtra = return_gridExtra)
}

#' @export
plot_variance_partition.df <- function(df, ...) UseMethod("plot_variance_partition.df")

.pb_shared_variance_limits <- function(vp_list) {
    if (!length(vp_list)) {
        return(NULL)
    }
    values <- unlist(lapply(vp_list, function(df) {
        if (is.null(df) || !"variance_explained" %in% names(df)) {
            return(numeric())
        }
        as.numeric(df$variance_explained)
    }), use.names = FALSE)
    values <- values[is.finite(values)]
    if (!length(values)) {
        return(NULL)
    }
    range(values)
}
