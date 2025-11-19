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
#' @export
plot_intragroup_variation.default <- function(data_matrix,
                                              sample_annotation,
                                              group_col,
                                              feature_id_col = "peptide_group_label",
                                              sample_id_col = "FullRunName",
                                              metrics = c("correlation", "PCV", "PMAD", "PEV"),
                                              correlation_method = c("pearson", "spearman", "kendall"),
                                              pcv_diff = FALSE,
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

    collect <- .pb_intragroup_process_matrix(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        feature_id_col = feature_id_col,
        group_cols = group_cols,
        assay_label = assay_label,
        assay_display = plot_title,
        plot_title = plot_title,
        call_info = metric_info$call,
        metric_name = metric_info$name,
        path_to_save_results = path_to_save_results
    )

    intragroup_plot <- .pb_build_intragroup_plot(
        data = collect$data,
        metric_label = metric_info$title,
        group_levels = group_cols,
        base_size = base_size,
        facet_ncol = metric_plot_ncol
    )

    attr(intragroup_plot, "pb_intragroup_metrics") <- collect$metrics
    attr(intragroup_plot, "pb_intragroup_saved_files") <- collect$saved_files

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
    metrics_arg <- .pb_intragroup_pop_arg(dots, "metrics", default = c("correlation", "PCV", "PMAD", "PEV"))
    correlation_method_arg <- .pb_intragroup_pop_arg(dots, "correlation_method", default = c("pearson", "spearman", "kendall"))
    pcv_diff_arg <- isTRUE(.pb_intragroup_pop_arg(dots, "pcv_diff", default = FALSE))
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

    default_sample_annotation <- as.data.frame(colData(object), stringsAsFactors = FALSE)
    rownames(default_sample_annotation) <- NULL
    sample_ann_list <- split_arg(sample_annotation)

    data_pieces <- vector("list", length(assays))
    metrics_export <- vector("list", length(assays))
    names(metrics_export) <- assays
    saved_paths <- vector("list", length(assays))
    names(saved_paths) <- assays

    for (i in seq_along(assays)) {
        assay_nm <- assays[[i]]
        dm <- pb_assay_matrix(object, assay_nm)
        sample_ann <- sample_ann_list[[i]]
        if (is.null(sample_ann)) {
            sample_ann <- default_sample_annotation
        }
        sample_ann <- as.data.frame(sample_ann, stringsAsFactors = FALSE)

        assay_path <- NULL
        if (!is.null(path_to_save_results)) {
            assay_path <- file.path(path_to_save_results, assay_nm)
        }

        collect <- .pb_intragroup_process_matrix(
            data_matrix = dm,
            sample_annotation = sample_ann,
            sample_id_col = sample_id_col,
            feature_id_col = feature_id_col,
            group_cols = group_cols,
            assay_label = assay_nm,
            assay_display = titles[[i]],
            plot_title = titles[[i]],
            call_info = metric_info$call,
            metric_name = metric_info$name,
            path_to_save_results = assay_path
        )

        data_pieces[[i]] <- collect$data
        metrics_export[[i]] <- collect$metrics
        saved_paths[[i]] <- collect$saved_files
    }

    combined_df <- dplyr::bind_rows(data_pieces)
    if (!nrow(combined_df)) {
        stop("No data available to plot intragroup variation.")
    }

    facet_cols <- plot_ncol
    if (is.null(facet_cols)) {
        facet_cols <- metric_plot_ncol
    }

    combined_plot <- .pb_build_intragroup_plot(
        data = combined_df,
        metric_label = metric_info$title,
        group_levels = group_cols,
        base_size = base_size,
        facet_ncol = facet_cols,
        shared_title = shared_title
    )

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
                                          group_col) {
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
    if (length(metric_choice) != 1L) {
        stop("Provide exactly one metric per call when plotting intragroup variation.")
    }
    metric_choice <- metric_choice[[1]]
    correlation_method <- match.arg(correlation_method)

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
        name = metric_choice,
        title = metric_titles[[metric_choice]],
        call = metric_calls[[metric_choice]]
    )
}

.pb_intragroup_process_matrix <- function(data_matrix,
                                          sample_annotation,
                                          sample_id_col,
                                          feature_id_col,
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
        group_col = group_cols
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

    group_data <- list()
    export_data <- list()
    for (grp in group_cols) {
        call_args <- c(list(se = se, ain = se_label, condition = grp), call_info$args)
        gg <- do.call(call_info$fun, call_args)
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
    if (!is.null(mapping$y)) {
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
        assay_line <- .pb_break_long_assay_title(assay_name)
        gg <- gg + ggplot2::labs(
            title = sprintf("%s\nAssay: %s", metric_label, assay_line),
            subtitle = NULL
        )
        gg <- gg + ggplot2::xlab("Grouping factor")
    } else {
        gg <- gg + ggplot2::facet_wrap(~group_column, ncol = facet_ncol)
        gg <- gg + ggplot2::scale_x_discrete(labels = function(x) vapply(x, .pb_break_long_assay_title, character(1)))
        subtitle <- NULL
        if (!is.null(shared_title) && nzchar(shared_title)) {
            subtitle <- shared_title
        }
        gg <- gg + ggplot2::labs(title = metric_label, subtitle = subtitle)
        gg <- gg + ggplot2::xlab("Assay")
    }

    gg
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
    value <- default
    if (!is.null(dots[[name]])) {
        value <- dots[[name]]
        dots[[name]] <- NULL
    }
    assign(deparse(substitute(dots)), dots, parent.frame())
    value
}

.pb_intragroup_pop_units <- function(dots, default = "cm") {
    opts <- c("cm", "in", "mm")
    value <- default
    if (!is.null(dots$units)) {
        value <- dots$units
        dots$units <- NULL
    }
    assign(deparse(substitute(dots)), dots, parent.frame())
    match.arg(as.character(value), opts)
}
