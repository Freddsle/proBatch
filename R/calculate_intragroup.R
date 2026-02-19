## Intragroup diagnostics backed by PRONE

#' Plot intragroup variation diagnostics
#'
#' proBatch provides wrappers around the intragroup evaluation plots implemented
#' in the [PRONE](https://github.com/daisybio/PRONE) package so that they can be
#' applied directly to matrices or `ProBatchFeatures` pipelines. The PRONE
#' functions are used as-is for the individual metrics, while proBatch aligns the
#' supplied sample annotation, collects the per-metric data, arranges the plots,
#' and optionally stores the tabular results.
#'
#' @inheritParams proBatch
#' @param data_matrix Numeric matrix (features × samples) when calling the
#'   `default` method.
#' @param group_col Sample annotation column(s) that define intragroup
#'   membership. Provide a single column name or a character vector; every
#'   requested column must exist in `sample_annotation` and contain no missing
#'   values.
#' @param metrics Character vector choosing one or more of
#'   `c("correlation", "PCV", "PMAD", "PEV")`.
#' @param correlation_method Correlation method passed to
#'   [PRONE::plot_intragroup_correlation()].
#' @param pcv_diff,pmad_diff,pev_diff Logical scalars forwarded to the
#'   corresponding PRONE plotting routines to request `diff = TRUE`.
#' @param assay_label Optional label used when constructing the temporary
#'   `SummarizedExperiment` that is passed to PRONE. When `NULL`, the function
#'   reuses the assay name (for `ProBatchFeatures`) or a generic label.
#' @param metric_plot_ncol Number of columns used when arranging the per-metric
#'   plots for a single assay.
#' @param path_to_save_results Optional directory where the per-metric data are
#'   written as CSV files. When operating on a `ProBatchFeatures` instance, set
#'   this to a parent directory; one sub-directory per assay is created.
#' @param fill_the_missing Missing-value policy applied to `data_matrix` before
#'   calling PRONE. Defaults to `FALSE` (keep missing values). See
#'   [handle_missing_values()] for the supported options.
#' @param filename Optional output path for the combined plot of all metrics
#'   (single assay).
#' @param base_size Base font size used when annotating metric titles.
#' @param pbf_name Names of assays to process when `data_matrix` is a
#'   `ProBatchFeatures` object (see [pb_assay_matrix()]).
#' @param plot_ncol Number of columns when arranging multiple assays (only for
#'   the `ProBatchFeatures` method).
#' @param return_gridExtra Logical; request the `gridExtra` object instead of a
#'   ggplot-ready object for multi-assay layouts.
#' @param ... Additional arguments forwarded to the `default` method (for
#'   example `metrics`, `correlation_method`, `pcv_diff`, and sizing arguments).
#'
#' @return A `ggplot` object showing the requested intragroup metric(s). The
#'   object carries two attributes: `pb_intragroup_metrics` (a list of exported
#'   metric tables) and `pb_intragroup_saved_files` (paths to saved CSV files
#'   when `path_to_save_results` is provided).
#'
#' @importFrom PRONE plot_intragroup_correlation plot_intragroup_PCV plot_intragroup_PMAD plot_intragroup_PEV
#' @name plot_intragroup_variation
#' @export
plot_intragroup_variation.default <- function(data_matrix,
                                              sample_annotation,
                                              group_col,
                                              feature_id_col = "peptide_group_label",
                                              sample_id_col = "FullRunName",
                                              fill_the_missing = FALSE,
                                              metrics = "correlation",
                                              correlation_method = c("pearson", "spearman", "kendall"),
                                              pcv_diff = TRUE,
                                              pmad_diff = FALSE,
                                              pev_diff = FALSE,
                                              assay_label = NULL,
                                              metric_plot_ncol = 2,
                                              filename = NULL,
                                              width = NA,
                                              height = NA,
                                              units = c("cm", "in", "mm"),
                                              path_to_save_results = NULL,
                                              plot_title = NULL,
                                              base_size = 12) {
    .pb_requireNamespace("PRONE")

    if (missing(group_col) || is.null(group_col)) {
        stop("`group_col` must contain at least one column name.")
    }
    if (is.null(sample_annotation)) {
        stop("`sample_annotation` is required for intragroup diagnostics.")
    }
    group_cols <- .pb_intragroup_normalize_group_cols(group_col)
    units <- match.arg(units)

    metric_info <- .pb_intragroup_metric_setup(
        metrics = metrics,
        correlation_method = correlation_method,
        pcv_diff = pcv_diff,
        pmad_diff = pmad_diff,
        pev_diff = pev_diff
    )
    metric_names <- metric_info$names
    if ("correlation" %in% metric_names && isTRUE(pcv_diff) && !("PCV" %in% metric_names)) {
        suppressWarnings(
            suppressMessages(
                try(
                    .pb_intragroup_process_matrix(
                        data_matrix = data_matrix,
                        sample_annotation = sample_annotation,
                        sample_id_col = sample_id_col,
                        feature_id_col = feature_id_col,
                        fill_the_missing = fill_the_missing,
                        group_cols = group_cols,
                        assay_label = assay_label,
                        assay_display = plot_title,
                        plot_title = plot_title,
                        call_info = metric_info$calls[["PCV"]],
                        metric_name = "PCV",
                        path_to_save_results = NULL
                    ),
                    silent = TRUE
                )
            )
        )
    }

    metric_plots <- list()
    metrics_export <- list()
    saved_paths <- list()
    for (metric_nm in metric_names) {
        collect <- .pb_intragroup_process_matrix(
            data_matrix = data_matrix,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            feature_id_col = feature_id_col,
            fill_the_missing = fill_the_missing,
            group_cols = group_cols,
            assay_label = assay_label,
            assay_display = plot_title,
            plot_title = plot_title,
            call_info = metric_info$calls[[metric_nm]],
            metric_name = metric_nm,
            path_to_save_results = path_to_save_results
        )

        metric_plots[[metric_nm]] <- .pb_build_intragroup_plot(
            data = collect$data,
            metric_label = metric_info$titles[[metric_nm]],
            group_levels = group_cols,
            base_size = base_size,
            facet_ncol = metric_plot_ncol
        )
        metrics_export[[metric_nm]] <- collect$metrics[[metric_nm]]
        if (!is.null(collect$saved_files) && metric_nm %in% names(collect$saved_files)) {
            saved_paths[[metric_nm]] <- collect$saved_files[[metric_nm]]
        }
    }

    intragroup_plot <- .pb_arrange_plot_list(
        plot_list = metric_plots,
        convert_fun = ggplot2::ggplotGrob,
        draw = FALSE,
        plot_ncol = metric_plot_ncol,
        return_gridExtra = FALSE
    )
    if (!length(saved_paths)) {
        saved_paths <- NULL
    }

    attr(intragroup_plot, "pb_intragroup_metrics") <- metrics_export
    attr(intragroup_plot, "pb_intragroup_saved_files") <- saved_paths

    save_ggplot(filename, units, width, height, intragroup_plot)

    intragroup_plot
}

#' @rdname plot_intragroup_variation
#' @method plot_intragroup_variation ProBatchFeatures
#' @export
plot_intragroup_variation.ProBatchFeatures <- function(data_matrix,
                                                       pbf_name = NULL,
                                                       sample_annotation = NULL,
                                                       group_col,
                                                       feature_id_col = "peptide_group_label",
                                                       sample_id_col = "FullRunName",
                                                       fill_the_missing = FALSE,
                                                       plot_title = NULL,
                                                       return_gridExtra = FALSE,
                                                       plot_ncol = NULL,
                                                       path_to_save_results = NULL,
                                                       metric_plot_ncol = 2,
                                                       ...) {
    .pb_requireNamespace("PRONE")

    if (missing(group_col) || is.null(group_col)) {
        stop("`group_col` must contain at least one column name.")
    }

    group_cols <- .pb_intragroup_normalize_group_cols(group_col)
    if (isTRUE(return_gridExtra)) {
        warning("return_gridExtra is ignored; returning a ggplot object instead.")
    }

    dots_input <- list(...)
    filename <- .pb_intragroup_pop_arg(dots_input, "filename", default = NULL)
    width <- .pb_intragroup_pop_arg(dots_input, "width", default = NA)
    height <- .pb_intragroup_pop_arg(dots_input, "height", default = NA)
    units <- .pb_intragroup_pop_units(dots_input, default = "cm")

    object <- data_matrix
    prep <- .pb_prepare_multi_assay(
        object = object,
        pbf_name = pbf_name,
        dots = dots_input,
        plot_title = plot_title,
        default_title_fun = function(x) x
    )
    assays <- prep$assays
    dots <- prep$dots
    split_arg <- prep$split_arg
    titles <- prep$titles
    shared_title <- prep$shared_title

    base_size <- .pb_intragroup_pop_arg(dots, "base_size", default = 12)
    metrics_arg <- .pb_intragroup_pop_arg(dots, "metrics", default = "correlation")
    correlation_method_arg <- .pb_intragroup_pop_arg(dots, "correlation_method", default = c("pearson", "spearman", "kendall"))
    pcv_diff_arg <- isTRUE(.pb_intragroup_pop_arg(dots, "pcv_diff", default = TRUE))
    pmad_diff_arg <- isTRUE(.pb_intragroup_pop_arg(dots, "pmad_diff", default = FALSE))
    pev_diff_arg <- isTRUE(.pb_intragroup_pop_arg(dots, "pev_diff", default = FALSE))
    if (length(dots)) {
        warning(sprintf("Ignoring unsupported arguments: %s", paste(names(dots), collapse = ", ")))
    }

    metric_info <- .pb_intragroup_metric_setup(
        metrics = metrics_arg,
        correlation_method = correlation_method_arg,
        pcv_diff = pcv_diff_arg,
        pmad_diff = pmad_diff_arg,
        pev_diff = pev_diff_arg
    )
    metric_names <- metric_info$names

    default_sample_annotation <- .pb_default_sample_annotation(
        object = object,
        sample_id_col = sample_id_col,
        drop_rownames = TRUE
    )
    sample_ann_list <- split_arg(sample_annotation)

    metrics_export <- vector("list", length(assays))
    names(metrics_export) <- assays
    saved_paths <- vector("list", length(assays))
    names(saved_paths) <- assays
    metric_plots <- list()

    facet_cols <- plot_ncol
    if (is.null(facet_cols)) {
        facet_cols <- metric_plot_ncol
    }

    for (metric_nm in metric_names) {
        data_pieces <- vector("list", length(assays))

        for (i in seq_along(assays)) {
            assay_nm <- assays[[i]]
            dm <- pb_assay_matrix(object, assay_nm)
            sample_ann <- sample_ann_list[[i]]
            if (is.null(sample_ann)) {
                sample_ann <- default_sample_annotation
            }
            sample_ann <- as.data.frame(sample_ann, stringsAsFactors = FALSE)

            display_label <- titles[[i]]
            if (is.null(display_label) || !length(display_label)) {
                display_label <- assay_nm
            } else {
                display_label <- as.character(display_label)[1]
                if (!nzchar(display_label)) {
                    display_label <- assay_nm
                }
            }

            assay_path <- NULL
            if (!is.null(path_to_save_results)) {
                assay_path <- file.path(path_to_save_results, assay_nm)
            }

            collect <- .pb_intragroup_process_matrix(
                data_matrix = dm,
                sample_annotation = sample_ann,
                sample_id_col = sample_id_col,
                feature_id_col = feature_id_col,
                fill_the_missing = fill_the_missing,
                group_cols = group_cols,
                assay_label = assay_nm,
                assay_display = display_label,
                plot_title = titles[[i]],
                call_info = metric_info$calls[[metric_nm]],
                metric_name = metric_nm,
                path_to_save_results = assay_path
            )

            data_pieces[[i]] <- collect$data
            metrics_export[[i]][[metric_nm]] <- collect$metrics[[metric_nm]]
            if (!is.null(collect$saved_files) && metric_nm %in% names(collect$saved_files)) {
                saved_paths[[i]][[metric_nm]] <- collect$saved_files[[metric_nm]]
            }
        }

        combined_df <- dplyr::bind_rows(data_pieces)
        if (!nrow(combined_df)) {
            stop("No data available to plot intragroup variation.")
        }

        metric_plots[[metric_nm]] <- .pb_build_intragroup_plot(
            data = combined_df,
            metric_label = metric_info$titles[[metric_nm]],
            group_levels = group_cols,
            base_size = base_size,
            facet_ncol = facet_cols,
            shared_title = shared_title
        )
    }

    combined_plot <- .pb_arrange_plot_list(
        plot_list = metric_plots,
        convert_fun = ggplot2::ggplotGrob,
        draw = FALSE,
        plot_ncol = metric_plot_ncol,
        return_gridExtra = FALSE
    )
    saved_paths <- lapply(saved_paths, function(x) {
        if (length(x)) {
            x
        } else {
            NULL
        }
    })

    attr(combined_plot, "pb_intragroup_metrics") <- metrics_export
    attr(combined_plot, "pb_intragroup_saved_files") <- saved_paths

    save_ggplot(filename, units, width, height, combined_plot)

    combined_plot
}

#' @export
plot_intragroup_variation <- function(data_matrix, ...) UseMethod("plot_intragroup_variation")

# Internal helpers -----------------------------------------------------------

.pb_prepare_intragroup_inputs <- function(data_matrix,
                                          sample_annotation,
                                          sample_id_col,
                                          feature_id_col,
                                          group_col,
                                          fill_the_missing = FALSE) {
    alignment <- .pb_align_matrix_and_annotation(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        allow_partial_annotation = FALSE
    )
    annotation <- as.data.frame(alignment$sample_annotation, stringsAsFactors = FALSE)
    missing_cols <- setdiff(group_col, colnames(annotation))
    if (length(missing_cols)) {
        stop(sprintf(
            "group_col value(s) not found in sample_annotation: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }
    for (col in group_col) {
        if (anyNA(annotation[[col]])) {
            stop(sprintf("group_col '%s' contains missing values for some samples.", col))
        }
    }
    aligned_matrix <- check_feature_id_col_in_dm(feature_id_col, alignment$data_matrix)
    if (!isFALSE(fill_the_missing)) {
        aligned_matrix <- .pb_handle_missing_wrapper(
            data_matrix = aligned_matrix,
            warning_message = "Intragroup diagnostics cannot operate with missing values; adjust `fill_the_missing` to control preprocessing.",
            fill_the_missing = fill_the_missing
        )
    }
    list(
        data_matrix = aligned_matrix,
        sample_annotation = annotation,
        sample_ids = alignment$sample_ids
    )
}

.pb_intragroup_assay_label <- function(assay_label, plot_title) {
    label <- assay_label
    if (is.null(label) || !nzchar(label)) {
        label <- plot_title
    }
    if (is.null(label) || !nzchar(label)) {
        label <- "assay"
    }
    as.character(label)[1]
}

.pb_intragroup_matrix_to_se <- function(data_matrix,
                                        sample_annotation,
                                        sample_id_col,
                                        feature_id_col,
                                        assay_label) {
    data_matrix <- as.matrix(data_matrix)
    sample_ids <- colnames(data_matrix)
    annot_ids <- as.character(sample_annotation[[sample_id_col]])
    if (anyNA(annot_ids)) {
        stop("Sample annotation contains missing IDs.")
    }
    if (!all(sample_ids %in% annot_ids)) {
        missing_ids <- setdiff(sample_ids, annot_ids)
        stop(sprintf(
            "Sample annotation is missing %s: %s",
            ifelse(length(missing_ids) == 1, "ID", "IDs"),
            paste(missing_ids, collapse = ", ")
        ))
    }
    sample_annotation <- sample_annotation[match(sample_ids, annot_ids), , drop = FALSE]
    rownames(sample_annotation) <- sample_ids
    col_data <- S4Vectors::DataFrame(sample_annotation)

    feature_ids <- rownames(data_matrix)
    if (is.null(feature_ids)) {
        feature_ids <- sprintf("feature_%s", seq_len(nrow(data_matrix)))
        rownames(data_matrix) <- feature_ids
    }
    row_data <- S4Vectors::DataFrame(feature_id = rownames(data_matrix))

    SummarizedExperiment::SummarizedExperiment(
        assays = S4Vectors::SimpleList(stats::setNames(list(data_matrix), assay_label)),
        rowData = row_data,
        colData = col_data
    )
}

.pb_save_intragroup_metrics <- function(metrics_data, base_path, assay_label) {
    if (is.null(base_path) || !nzchar(base_path) || !length(metrics_data)) {
        return(NULL)
    }
    if (!dir.exists(base_path)) {
        dir.create(base_path, recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(base_path)) {
        warning(sprintf("Unable to create directory '%s'; skipping metric export.", base_path))
        return(NULL)
    }
    saved <- list()
    safe_label <- .pb_sanitize_filename_component(assay_label)
    for (nm in names(metrics_data)) {
        metric_df <- metrics_data[[nm]]
        if (is.null(metric_df)) {
            next
        }
        file_name <- sprintf("intragroup_%s_%s.csv", tolower(nm), safe_label)
        file_path <- file.path(base_path, file_name)
        utils::write.csv(metric_df, file = file_path, row.names = FALSE)
        saved[[nm]] <- file_path
    }
    saved
}

.pb_sanitize_filename_component <- function(value) {
    if (is.null(value) || !nzchar(value)) {
        return("assay")
    }
    clean <- gsub("[^a-zA-Z0-9]+", "_", value)
    clean <- gsub("_+", "_", clean)
    trimws(gsub("^_|_$", "", clean))
}

.pb_intragroup_normalize_group_cols <- function(group_col) {
    vals <- unique(as.character(group_col))
    vals <- vals[nzchar(vals)]
    if (!length(vals)) {
        stop("`group_col` must contain at least one non-empty column name.")
    }
    vals
}

.pb_intragroup_metric_setup <- function(metrics,
                                        correlation_method,
                                        pcv_diff,
                                        pmad_diff,
                                        pev_diff) {
    allowed_metrics <- c("correlation", "PCV", "PMAD", "PEV")
    metric_choice <- match.arg(metrics, allowed_metrics, several.ok = TRUE)
    metric_choice <- unique(metric_choice)
    correlation_method <- match.arg(correlation_method, c("pearson", "spearman", "kendall"))

    metric_titles <- c(
        correlation = "Intragroup correlation",
        PCV = "Intragroup PCV",
        PMAD = "Intragroup PMAD",
        PEV = "Intragroup PEV"
    )
    metric_calls <- list(
        correlation = list(fun = PRONE::plot_intragroup_correlation, args = list(method = correlation_method)),
        PCV = list(fun = PRONE::plot_intragroup_PCV, args = list(diff = pcv_diff)),
        PMAD = list(fun = PRONE::plot_intragroup_PMAD, args = list(diff = pmad_diff)),
        PEV = list(fun = PRONE::plot_intragroup_PEV, args = list(diff = pev_diff))
    )

    list(
        names = metric_choice,
        titles = metric_titles[metric_choice],
        calls = metric_calls
    )
}

.pb_intragroup_process_matrix <- function(data_matrix,
                                          sample_annotation,
                                          sample_id_col,
                                          feature_id_col,
                                          fill_the_missing,
                                          group_cols,
                                          assay_label,
                                          assay_display,
                                          plot_title,
                                          call_info,
                                          metric_name,
                                          path_to_save_results = NULL) {
    prep <- .pb_prepare_intragroup_inputs(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        feature_id_col = feature_id_col,
        group_col = group_cols,
        fill_the_missing = fill_the_missing
    )
    aligned_matrix <- prep$data_matrix
    aligned_annotation <- prep$sample_annotation

    se_label <- .pb_intragroup_assay_label(assay_label, plot_title)
    display_label <- assay_display
    if (is.null(display_label) || !nzchar(display_label)) {
        display_label <- se_label
    }
    se <- .pb_intragroup_matrix_to_se(
        data_matrix = aligned_matrix,
        sample_annotation = aligned_annotation,
        sample_id_col = sample_id_col,
        feature_id_col = feature_id_col,
        assay_label = se_label
    )

    call_info_local <- call_info
    group_data <- list()
    export_data <- list()
    for (grp in group_cols) {
        call_result <- .pb_intragroup_call_metric(
            call_info = call_info_local,
            se = se,
            assay_label = se_label,
            group_col = grp,
            metric_name = metric_name
        )
        gg <- call_result$plot
        if (isTRUE(call_result$fallback_to_no_diff)) {
            call_info_local$args$diff <- FALSE
        }
        extracted <- .pb_intragroup_extract_plot_data(gg)
        if (is.null(extracted)) {
            next
        }

        values <- extracted$value
        if (length(values)) {
            group_data[[length(group_data) + 1L]] <- data.frame(
                value = values,
                group_column = grp,
                assay_display = display_label,
                stringsAsFactors = FALSE
            )
        }

        export_df <- as.data.frame(extracted$raw, stringsAsFactors = FALSE)
        export_df$group_column <- grp
        export_data[[length(export_data) + 1L]] <- export_df
    }

    if (!length(group_data)) {
        stop("PRONE did not return any data to plot.")
    }

    combined <- dplyr::bind_rows(group_data)
    export_df <- if (length(export_data)) {
        dplyr::bind_rows(export_data)
    } else {
        data.frame(stringsAsFactors = FALSE)
    }
    metrics <- list()
    metrics[[metric_name]] <- export_df

    saved_files <- .pb_save_intragroup_metrics(
        metrics_data = metrics,
        base_path = path_to_save_results,
        assay_label = se_label
    )

    list(
        data = combined,
        metrics = metrics,
        saved_files = saved_files
    )
}

.pb_intragroup_call_metric <- function(call_info,
                                       se,
                                       assay_label,
                                       group_col,
                                       metric_name) {
    call_args <- c(list(se = se, ain = assay_label, condition = group_col), call_info$args)
    result <- tryCatch(
        do.call(call_info$fun, call_args),
        error = function(e) e
    )
    if (!inherits(result, "error")) {
        return(list(plot = result, fallback_to_no_diff = FALSE))
    }

    error_msg <- conditionMessage(result)
    allow_retry <- isTRUE(call_info$args$diff) &&
        grepl("Difference not applicable", error_msg, fixed = TRUE)

    if (!allow_retry) {
        stop(result)
    }

    call_args$diff <- FALSE
    warning(sprintf(
        "Metric '%s' with diff = TRUE failed for assay '%s' and grouping '%s'; retrying with diff = FALSE.",
        metric_name,
        assay_label,
        group_col
    ), call. = FALSE)

    list(
        plot = do.call(call_info$fun, call_args),
        fallback_to_no_diff = TRUE
    )
}

.pb_intragroup_extract_plot_data <- function(gg) {
    if (is.null(gg)) {
        return(NULL)
    }
    data <- gg$data
    if (is.null(data) || !nrow(data)) {
        built <- ggplot2::ggplot_build(gg)
        data <- built$plot$data
    }
    if (is.null(data) || !nrow(data)) {
        return(NULL)
    }
    mapping <- gg$mapping
    value_col <- NULL
    if (!is.null(mapping$y) && rlang::is_symbol(mapping$y)) {
        value_col <- rlang::as_name(mapping$y)
    }
    if (is.null(value_col) || !value_col %in% names(data)) {
        numeric_cols <- names(Filter(is.numeric, data))
        if (!length(numeric_cols)) {
            return(NULL)
        }
        value_col <- numeric_cols[[1]]
    }
    list(
        value = as.numeric(data[[value_col]]),
        raw = data
    )
}

.pb_build_intragroup_plot <- function(data,
                                      metric_label,
                                      group_levels,
                                      base_size,
                                      facet_ncol = 2,
                                      shared_title = NULL) {
    if (is.null(data) || !nrow(data)) {
        stop("No intragroup data available to build the plot.")
    }
    df <- data
    df$group_column <- factor(df$group_column, levels = group_levels)
    assay_levels <- unique(df$assay_display)
    df$assay_display <- factor(df$assay_display, levels = assay_levels)

    n_assays <- length(levels(df$assay_display))
    color_values <- .pb_intragroup_group_colors(levels(df$group_column))

    if (n_assays == 1L) {
        gg <- ggplot2::ggplot(df, ggplot2::aes(x = group_column, y = value, fill = group_column))
    } else {
        gg <- ggplot2::ggplot(df, ggplot2::aes(x = assay_display, y = value, fill = group_column))
    }

    gg <- gg + ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.9, color = "#333333", size = 0.3)
    gg <- gg + ggplot2::scale_fill_manual(values = color_values, drop = FALSE, name = "Grouping factor")
    gg <- gg + ggplot2::theme_classic(base_size = base_size) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
            strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
            strip.text = ggplot2::element_text(face = "bold")
        )
    gg <- gg + ggplot2::ylab(metric_label)

    if (n_assays == 1L) {
        assay_name <- levels(df$assay_display)[1]
        assay_line <- .pb_intragroup_format_assay_label(assay_name)
        x_axis_labels <- levels(df$group_column)
        gg <- gg + ggplot2::labs(
            title = sprintf("%s\nAssay: %s", metric_label, assay_line),
            subtitle = NULL
        )
        gg <- gg + ggplot2::xlab("Grouping factor")
    } else {
        gg <- gg + ggplot2::facet_wrap(~group_column, ncol = facet_ncol)
        subtitle_parts <- character()
        if (!is.null(shared_title) && nzchar(shared_title)) {
            subtitle_parts <- c(subtitle_parts, shared_title)
        }
        level_label <- .pb_intragroup_common_level(assay_levels)
        include_level_in_labels <- is.null(level_label)
        x_axis_labels <- vapply(
            assay_levels,
            .pb_intragroup_format_assay_label,
            character(1),
            include_level = include_level_in_labels
        )
        x_axis_label_map <- stats::setNames(x_axis_labels, assay_levels)
        gg <- gg + ggplot2::scale_x_discrete(labels = function(x) {
            unname(x_axis_label_map[as.character(x)])
        })
        if (!is.null(level_label) && nzchar(level_label) && !level_label %in% subtitle_parts) {
            subtitle_parts <- c(subtitle_parts, level_label)
        }
        subtitle <- if (length(subtitle_parts)) paste(subtitle_parts, collapse = "\n") else NULL
        gg <- gg + ggplot2::labs(title = metric_label, subtitle = subtitle)
        gg <- gg + ggplot2::xlab("Assay")
    }

    if (.pb_intragroup_should_rotate_x_labels(x_axis_labels, max_chars = 15L)) {
        gg <- gg + ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
        )
    }

    gg
}

.pb_intragroup_common_level <- function(assay_labels) {
    if (!length(assay_labels)) {
        return(NULL)
    }
    levels <- vapply(
        assay_labels,
        .pb_intragroup_assay_level,
        character(1),
        USE.NAMES = FALSE
    )
    levels <- unique(levels[nzchar(levels)])
    if (length(levels) == 1L) {
        levels[[1]]
    } else {
        NULL
    }
}

.pb_intragroup_should_rotate_x_labels <- function(labels, max_chars = 15L) {
    if (is.null(labels) || !length(labels)) {
        return(FALSE)
    }
    labels <- as.character(labels)
    labels <- labels[!is.na(labels) & nzchar(labels)]
    if (!length(labels)) {
        return(FALSE)
    }

    max_chars_per_label <- vapply(labels, function(label) {
        lines <- strsplit(label, "\n", fixed = TRUE)[[1]]
        lines <- trimws(lines)
        lines <- lines[nzchar(lines)]
        if (!length(lines)) {
            return(0L)
        }
        max(nchar(lines, type = "chars"))
    }, integer(1))

    any(max_chars_per_label > max_chars)
}

.pb_intragroup_group_colors <- function(group_levels) {
    n <- length(group_levels)
    if (!n) {
        return(character())
    }
    cols <- scales::hue_pal()(max(n, 3))
    cols <- cols[seq_len(n)]
    stats::setNames(cols, group_levels)
}

.pb_intragroup_pop_arg <- function(dots, name, default = NULL) {
    if (!is.list(dots)) {
        stop("dots must be a list")
    }
    dots_name <- deparse(substitute(dots))
    value <- default
    if (length(names(dots)) && name %in% names(dots)) {
        if (!is.null(dots[[name]])) {
            value <- dots[[name]]
        }
        dots[[name]] <- NULL
    }
    assign(dots_name, dots, parent.frame())
    value
}

.pb_intragroup_pop_units <- function(dots, default = "cm") {
    opts <- c("cm", "in", "mm")
    dots_name <- deparse(substitute(dots))
    value <- default
    if (!is.null(dots$units)) {
        value <- dots$units
        dots$units <- NULL
    }
    assign(dots_name, dots, parent.frame())
    match.arg(as.character(value), opts)
}
