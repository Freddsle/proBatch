#' Visualise correlation matrix
#'
#' Plot correlation of selected  samples or peptides
#' @description recommended for heatmap-type visualisation of correlation matrix
#' with <100 items. With >50 samples and ~10 replicate pairs distribution plots
#' may be more informative.
#'
#' @inheritParams proBatch
#' @param corr_matrix square correlation matrix
#' @param cluster_rows boolean values determining if rows should be clustered or \code{hclust} object
#' @param cluster_cols boolean values determining if columns should be clustered or \code{hclust} object
#' @param heatmap_color vector of colors used in heatmap.
#' @param annotation data frame with \code{peptide_annotation} for protein
#' correlation heatmap or \code{sample_annotation} for sample correlation heatmap
#' @param annotation_id_col \code{feature_id_col} for protein correlation heatmap
#' or \code{sample_id_col} for sample correlation heatmap
#' @param ... parameters for the \code{\link[pheatmap]{pheatmap}} visualisation,
#'  for details see examples and help to corresponding functions
#'
#' @return \code{pheatmap} object
#'
#' @export
#'
#' @seealso \code{\link[pheatmap]{pheatmap}},
#' \code{\link{plot_sample_corr_distribution}},
#' \code{\link{plot_peptide_corr_distribution}}
#'
#' @examples
#' data("example_proteome_matrix", package = "proBatch")
#' peptides <- c("10231_QDVDVWLWQQEGSSK_2", "10768_RLESELDGLR_2")
#' data_matrix_sub <- example_proteome_matrix[peptides, ]
#' corr_matrix <- cor(t(data_matrix_sub), use = "complete.obs")
#' corr_matrix_plot <- plot_corr_matrix(corr_matrix)
#'
plot_corr_matrix <- function(corr_matrix,
                             annotation = NULL,
                             annotation_id_col = "FullRunName",
                             factors_to_plot = NULL,
                             cluster_rows = FALSE, cluster_cols = FALSE,
                             heatmap_color = colorRampPalette(
                                 rev(brewer.pal(n = 7, name = "RdYlBu"))
                             )(100),
                             color_list = NULL,
                             filename = NULL, width = 7, height = 7,
                             units = c("cm", "in", "mm"),
                             plot_title = NULL, ...) {
    # infer the color scheme for annotation (cols & rows)
    color_list <- .pb_resolve_color_list(
        color_list = color_list,
        annotation_df = annotation,
        id_col = annotation_id_col,
        columns = factors_to_plot,
        warn_message = "color_list for annotation (cols & rows) not defined, inferring automatically. Numeric/factor columns are guessed, for more controlled color mapping use sample_annotation_to_colors()"
    )

    if (cluster_rows != cluster_cols) {
        warning("different arguments for clustering of rows and columns, this will make
            correlation matrix heatmap asymmetrical!")
    }
    p <- plot_heatmap_generic(
        corr_matrix,
        column_annotation_df = annotation,
        row_annotation_df = annotation,
        fill_the_missing = NULL,
        col_ann_id_col = annotation_id_col,
        row_ann_id_col = annotation_id_col,
        columns_for_cols = factors_to_plot,
        columns_for_rows = factors_to_plot,
        cluster_rows = cluster_rows, cluster_cols = cluster_cols,
        annotation_color_cols = color_list,
        annotation_color_rows = color_list,
        heatmap_color = heatmap_color,
        filename = filename,
        width = width, height = height,
        units = units,
        plot_title = plot_title,
        ...
    )
    return(p)
}

.pb_corr_resolve_sample_input <- function(data_matrix,
                                          sample_annotation,
                                          sample_id_col,
                                          pbf_name = NULL,
                                          sample_annotation_missing = FALSE,
                                          require_annotation = FALSE) {
    if (is(data_matrix, "ProBatchFeatures")) {
        object <- data_matrix
        assay_name <- .pb_resolve_assay_for_input(
            object = object,
            pbf_name = pbf_name
        )
        data_matrix <- pb_assay_matrix(object, assay = assay_name)
        sample_annotation <- .pb_default_sample_annotation(
            object = object,
            sample_annotation = if (sample_annotation_missing) NULL else sample_annotation,
            sample_id_col = sample_id_col,
            sample_ids = colnames(data_matrix)
        )
    } else if (sample_annotation_missing && require_annotation) {
        stop("`sample_annotation` must be provided.")
    }
    list(data_matrix = data_matrix, sample_annotation = sample_annotation)
}

.pb_corr_resolve_feature_input <- function(data_matrix,
                                           peptide_annotation,
                                           feature_id_col,
                                           pbf_name = NULL,
                                           peptide_annotation_missing = FALSE,
                                           require_annotation = FALSE) {
    if (is(data_matrix, "ProBatchFeatures")) {
        object <- data_matrix
        assay_name <- .pb_resolve_assay_for_input(
            object = object,
            pbf_name = pbf_name
        )
        data_matrix <- pb_assay_matrix(object, assay = assay_name)
        peptide_annotation <- .pb_default_feature_annotation(
            object = object,
            assay_name = assay_name,
            feature_annotation = if (peptide_annotation_missing) NULL else peptide_annotation,
            feature_id_col = feature_id_col
        )
    } else if (peptide_annotation_missing && require_annotation) {
        stop("`peptide_annotation` must be provided.")
    }
    list(data_matrix = data_matrix, peptide_annotation = peptide_annotation)
}

.pb_finalize_corr_distribution_plot <- function(gg,
                                                corr_distribution,
                                                plot_title,
                                                theme,
                                                base_size,
                                                filename,
                                                units,
                                                width,
                                                height) {
    if (!is.null(plot_title)) {
        gg <- gg + ggtitle(plot_title)
    }
    if (("Step" %in% names(corr_distribution)) &&
        length(unique(corr_distribution$Step)) > 1) {
        if (length(unique(corr_distribution$Step)) <= 4) {
            gg <- gg + facet_grid(. ~ Step)
        } else {
            gg <- gg + facet_grid(Step ~ .)
        }
    }
    if (!is.null(theme) && theme == "classic") {
        gg <- gg + theme_classic(base_size = base_size)
    } else {
        message("plotting with default ggplot theme, only theme = 'classic' implemented")
    }
    gg <- gg + theme(plot.title = element_text(hjust = .5, face = "bold"))
    save_ggplot(filename, units, width, height, gg)
    gg
}

.pb_corr_distribution_from_input <- function(data_matrix,
                                             builder,
                                             step_as_factor = FALSE) {
    if (!is.list(data_matrix)) {
        return(builder(data_matrix))
    }
    corr_distribution <- lapply(seq_len(length(data_matrix)), function(i) {
        out <- builder(data_matrix[[i]])
        out$Step <- names(data_matrix)[i]
        out
    })
    corr_distribution <- do.call(rbind, corr_distribution)
    if (isTRUE(step_as_factor)) {
        corr_distribution <- corr_distribution %>%
            mutate(Step = factor(Step, levels = names(data_matrix)))
    }
    corr_distribution
}

#' Peptide correlation matrix (heatmap)
#'
#' Plots correlation plot of peptides from a single protein
#'
#' @inheritParams proBatch
#' @param protein_name the name of the protein
#' @param cluster_rows boolean values determining if rows should be clustered or \code{hclust} object
#' @param cluster_cols boolean values determining if columns should be clustered or \code{hclust} object
#' @param heatmap_color vector of colors used in heatmap.
#' @param ... parameters for the corrplot visualisation
#'
#' @return \code{pheatmap} object
#'
#' @export
#' @examples
#' data(list = c("example_peptide_annotation", "example_proteome_matrix"), package = "proBatch")
#' protein_corrplot_plot <- plot_protein_corrplot(example_proteome_matrix,
#'     protein_name = "Haao", peptide_annotation = example_peptide_annotation,
#'     protein_col = "Gene"
#' )
#'
plot_protein_corrplot <- function(data_matrix,
                                  protein_name,
                                  peptide_annotation = NULL,
                                  protein_col = "ProteinName",
                                  feature_id_col = "peptide_group_label",
                                  factors_to_plot = c("ProteinName"),
                                  cluster_rows = FALSE, cluster_cols = FALSE,
                                  heatmap_color = colorRampPalette(
                                      rev(brewer.pal(n = 7, name = "RdYlBu"))
                                  )(100),
                                  color_list = NULL,
                                  filename = NULL,
                                  width = NA, height = NA,
                                  units = c("cm", "in", "mm"),
                                  plot_title = NULL, ...) {
    resolved <- .pb_corr_resolve_feature_input(
        data_matrix = data_matrix,
        peptide_annotation = peptide_annotation,
        feature_id_col = feature_id_col
    )
    data_matrix <- resolved$data_matrix
    peptide_annotation <- resolved$peptide_annotation

    if (is.null(peptide_annotation)) {
        stop("`peptide_annotation` must be provided for protein correlation plots.")
    }
    if (!feature_id_col %in% names(peptide_annotation)) {
        stop(sprintf("Feature ID column '%s' was not found in `peptide_annotation`.", feature_id_col))
    }
    if (!protein_col %in% names(peptide_annotation)) {
        stop(sprintf("Protein column '%s' was not found in `peptide_annotation`.", protein_col))
    }

    peptides <- peptide_annotation %>%
        filter(!!(sym(feature_id_col)) %in% rownames(data_matrix)) %>%
        filter(!!(sym(protein_col)) %in% protein_name) %>%
        pull(!!sym(feature_id_col)) %>%
        as.character()

    peptides <- unique(peptides[!is.na(peptides) & nzchar(peptides)])
    if (!length(peptides)) {
        stop(
            "No peptides from the selected protein(s) were found in `data_matrix`.",
            call. = FALSE
        )
    }

    data_matrix_sub <- data_matrix[peptides, , drop = FALSE]
    corr_matrix <- cor(t(data_matrix_sub), use = "pairwise.complete.obs")
    if (is.null(dim(corr_matrix))) {
        corr_matrix <- matrix(
            corr_matrix,
            nrow = 1L,
            ncol = 1L,
            dimnames = list(peptides[[1]], peptides[[1]])
        )
    }

    peptide_annotation <- peptide_annotation %>%
        filter(!!(sym(protein_col)) %in% protein_name) %>%
        filter(!!(sym(feature_id_col)) %in% rownames(corr_matrix)) %>%
        arrange(!!sym(protein_col))

    ordered_peptides <- as.character(peptide_annotation[[feature_id_col]])
    if (!length(ordered_peptides)) {
        stop(
            "No annotated peptides from the selected protein(s) could be aligned to the correlation matrix.",
            call. = FALSE
        )
    }

    corr_matrix <- corr_matrix[
        ordered_peptides,
        ordered_peptides,
        drop = FALSE
    ]

    if (nrow(corr_matrix) < 2L || ncol(corr_matrix) < 2L) {
        if (isTRUE(cluster_rows) || isTRUE(cluster_cols)) {
            message("Only one peptide available; disabling clustering for correlation heatmap.")
        }
        cluster_rows <- FALSE
        cluster_cols <- FALSE
    }

    if (is.null(plot_title) && length(protein_name) == 1) {
        plot_title <- sprintf("Correlation matrix of peptides from %s", protein_name)
    } else if (is.null(plot_title) && length(protein_name) > 1) {
        plot_title <- sprintf("Peptide correlation matrix of %s proteins", paste(protein_name, collapse = ", "))
    }

    plot_corr_matrix(corr_matrix,
        annotation = peptide_annotation,
        annotation_id_col = feature_id_col,
        factors_to_plot = factors_to_plot,
        cluster_rows = cluster_rows, cluster_cols = cluster_cols,
        heatmap_color = heatmap_color,
        color_list = color_list,
        plot_title = plot_title,
        filename = filename, width = width,
        height = height, units = units, ...
    )
}

#' Sample correlation matrix (heatmap)
#'
#' Plot correlation of selected samples
#'
#' @inheritParams proBatch
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames, or a
#'   `ProBatchFeatures` object. When `data_matrix` is a
#'   `ProBatchFeatures` object, `pbf_name` is used (or all assays when
#'   `pbf_name = NULL`).
#' @param sample_annotation data frame with sample-level metadata. When
#'   `data_matrix` is a matrix, this argument is required. When
#'   `data_matrix` is a `ProBatchFeatures` object and
#'   `sample_annotation` is not provided, `as.data.frame(colData(data_matrix))`
#'   is used.
#' @param pbf_name Assay name(s) used when `data_matrix` is a
#'   `ProBatchFeatures` object. If `NULL`, all assays are plotted.
#' @param samples_to_plot string vector of samples in
#' \code{data_matrix} to be used in the plot
#' @param cluster_rows boolean values determining if rows should be clustered or \code{hclust} object
#' @param cluster_cols boolean values determining if columns should be clustered or \code{hclust} object
#' @param heatmap_color vector of colors used in heatmap.
#' @param ... parameters for the \code{\link[pheatmap]{pheatmap}} visualisation, for details see
#'   examples and help to corresponding functions
#'
#' @return \code{pheatmap} object for a single assay, or an arranged plot object
#'   when multiple assays are plotted.
#'
#' @export
#'
#' @examples
#' data(list = c("example_sample_annotation", "example_proteome_matrix"), package = "proBatch")
#' specified_samples <- example_sample_annotation$FullRunName[
#'     which(example_sample_annotation$order %in% 110:115)
#' ]
#'
#' sample_corr_heatmap <- plot_sample_corr_heatmap(example_proteome_matrix,
#'     samples_to_plot = specified_samples,
#'     factors_to_plot = c("MS_batch", "Diet", "DateTime", "digestion_batch"),
#'     cluster_rows = FALSE, cluster_cols = FALSE,
#'     annotation_names_col = TRUE, annotation_legend = FALSE,
#'     show_colnames = FALSE
#' )
#'
#' @seealso \code{\link[pheatmap]{pheatmap}}
#'

plot_sample_corr_heatmap <- function(data_matrix, samples_to_plot = NULL,
                                     sample_annotation = NULL,
                                     sample_id_col = "FullRunName",
                                     factors_to_plot = NULL,
                                     cluster_rows = FALSE, cluster_cols = FALSE,
                                     heatmap_color = colorRampPalette(
                                         rev(brewer.pal(n = 7, name = "RdYlBu"))
                                     )(100),
                                     color_list = NULL,
                                     filename = NULL,
                                     width = NA, height = NA,
                                     units = c("cm", "in", "mm"),
                                     plot_title = sprintf(
                                         "Correlation matrix of%s samples",
                                         ifelse(is.null(samples_to_plot), "", " selected")
                                     ),
                                     pbf_name = NULL, ...) {
    plot_title_missing <- missing(plot_title)

    if (is(data_matrix, "ProBatchFeatures")) {
        object <- data_matrix
        prep <- .pb_prepare_multi_assay(
            object = object,
            pbf_name = pbf_name,
            dots = c(list(filename = filename), list(...)),
            plot_title = if (isTRUE(plot_title_missing)) NULL else plot_title,
            default_title_fun = function(x) x,
            set_silent = TRUE
        )
        assays <- prep$assays
        dots <- prep$dots
        filename_list <- prep$filename_list
        split_arg <- prep$split_arg
        titles <- prep$titles
        shared_title <- prep$shared_title

        default_sample_annotation <- .pb_default_sample_annotation(
            object = object,
            sample_id_col = sample_id_col
        )
        sample_ann_list <- split_arg(sample_annotation)

        plot_list <- vector("list", length(assays))
        names(plot_list) <- assays

        for (i in seq_along(assays)) {
            assay_nm <- assays[[i]]
            assay_matrix <- pb_assay_matrix(object, assay = assay_nm)
            sample_ann <- sample_ann_list[[i]]
            if (is.null(sample_ann)) {
                sample_ann <- default_sample_annotation
            }

            call_args <- .pb_per_assay_dots(dots, filename_list, i)
            call_args <- c(list(
                data_matrix = assay_matrix,
                samples_to_plot = samples_to_plot,
                sample_annotation = sample_ann,
                sample_id_col = sample_id_col,
                factors_to_plot = factors_to_plot,
                cluster_rows = cluster_rows,
                cluster_cols = cluster_cols,
                heatmap_color = heatmap_color,
                color_list = color_list,
                width = width,
                height = height,
                units = units,
                plot_title = titles[i]
            ), call_args)

            plot_list[[i]] <- do.call(.pb_plot_sample_corr_heatmap_single, call_args)
        }

        plot_list <- .pb_attach_shared_title(plot_list, shared_title)

        return(.pb_arrange_plot_list(
            plot_list = plot_list,
            convert_fun = function(x) x$gtable
        ))
    }

    resolved <- .pb_corr_resolve_sample_input(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        pbf_name = pbf_name
    )

    .pb_plot_sample_corr_heatmap_single(
        data_matrix = resolved$data_matrix,
        samples_to_plot = samples_to_plot,
        sample_annotation = resolved$sample_annotation,
        sample_id_col = sample_id_col,
        factors_to_plot = factors_to_plot,
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        heatmap_color = heatmap_color,
        color_list = color_list,
        plot_title = plot_title,
        filename = filename,
        width = width,
        height = height,
        units = units,
        ...
    )
}

.pb_plot_sample_corr_heatmap_single <- function(data_matrix,
                                                samples_to_plot,
                                                sample_annotation,
                                                sample_id_col,
                                                factors_to_plot,
                                                cluster_rows,
                                                cluster_cols,
                                                heatmap_color,
                                                color_list = NULL,
                                                filename = NULL,
                                                width = NA,
                                                height = NA,
                                                units = c("cm", "in", "mm"),
                                                plot_title = NULL,
                                                ...) {
    if (!is.null(samples_to_plot)) {
        if (!all(samples_to_plot %in% colnames(data_matrix))) {
            missing_samples <- setdiff(samples_to_plot, colnames(data_matrix))
            stop(sprintf(
                "The following samples are not in data matrix and can not
                    be used in sample correlation plotting %s",
                paste(missing_samples, collapse = ";\n")
            ))
        }
        corr_matrix <- cor(data_matrix[, samples_to_plot], use = "complete.obs")
    } else {
        corr_matrix <- cor(data_matrix, use = "complete.obs")
    }

    if (!is.null(sample_annotation)) {
        if (!all(samples_to_plot %in% sample_annotation[[sample_id_col]])) {
            warning("some of the samples are not in annotation, this may lead to problems in color annotation")
        }
    }

    plot_corr_matrix(
        corr_matrix,
        annotation = sample_annotation,
        annotation_id_col = sample_id_col,
        factors_to_plot = factors_to_plot,
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        heatmap_color = heatmap_color,
        color_list = color_list,
        plot_title = plot_title,
        filename = filename,
        width = width,
        height = height,
        units = units,
        ...
    )
}

get_sample_corr_df <- function(cor_proteome, sample_annotation,
                               sample_id_col = "FullRunName",
                               biospecimen_id_col = "EarTag",
                               batch_col = "MS_batch") {
    comb_to_keep <- data.frame(t(combn(colnames(cor_proteome), 2)))
    names(comb_to_keep) <- paste(sample_id_col, seq_len(2), sep = "_")

    spec_cols <- c(biospecimen_id_col, batch_col)

    if (!all(spec_cols %in% names(sample_annotation))) {
        missing_cols <- setdiff(spec_cols, names(sample_annotation))
        stop(sprintf(
            "Columns %s are not in sample_annotation. Please provide valid biospecimen_id_col and batch_col.",
            paste(missing_cols, collapse = ", ")
        ))
    }

    first_sample_col <- paste(sample_id_col, "1", sep = "_")
    second_sample_col <- paste(sample_id_col, "2", sep = "_")

    corr_distribution <- cor_proteome %>%
        as.data.frame() %>%
        rownames_to_column(var = first_sample_col) %>%
        pivot_longer(
            cols = -all_of(first_sample_col),
            names_to = second_sample_col,
            values_to = "correlation",
            values_drop_na = FALSE
        ) %>%
        merge(comb_to_keep) %>%
        merge(sample_annotation %>% select(all_of(c(sample_id_col, spec_cols))),
            by.x = paste(sample_id_col, "1", sep = "_"),
            by.y = sample_id_col, all.x = TRUE
        ) %>%
        setnames(
            old = spec_cols,
            new = paste(spec_cols, 1, sep = "")
        ) %>%
        merge(sample_annotation %>% select(all_of(c(sample_id_col, spec_cols))),
            by.x = paste(sample_id_col, "2", sep = "_"),
            by.y = sample_id_col, all.x = TRUE
        ) %>%
        setnames(
            old = spec_cols,
            new = paste(spec_cols, 2, sep = "")
        ) %>%
        mutate(replicate = (!!sym(paste(biospecimen_id_col, "1", sep = "")) ==
            !!sym(paste(biospecimen_id_col, "2", sep = "")))) %>%
        mutate(
            batch_the_same = (!!sym(paste(batch_col, "1", sep = "")) ==
                !!sym(paste(batch_col, "2", sep = ""))),
            batches = paste(!!sym(paste(batch_col, "1", sep = "")),
                !!sym(paste(batch_col, "2", sep = "")),
                sep = ":"
            )
        ) %>%
        mutate(batch_replicate = ifelse(replicate,
            ifelse(batch_the_same,
                "same_batch\nsame_biospecimen",
                "same_biospecimen\ndiff_batch"
            ),
            ifelse(batch_the_same,
                "same_batch\ndiff_biospecimen",
                "diff_batch\ndiff_biospecimen"
            )
        ))
    return(corr_distribution)
}


#' Calculates correlation for all pairs of the samples in data matrix, labels
#' as replicated/same_batch/unrelated in output columns (see "Value").
#'
#' @inheritParams proBatch
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames, or a
#'   `ProBatchFeatures` object. When `data_matrix` is a
#'   `ProBatchFeatures` object, `pbf_name` is used (or the current assay when
#'   `pbf_name = NULL`).
#' @param sample_annotation data frame with sample-level metadata. When
#'   `data_matrix` is a matrix, this argument is required. When
#'   `data_matrix` is a `ProBatchFeatures` object and
#'   `sample_annotation` is not provided, `as.data.frame(colData(data_matrix))`
#'   is used.
#' @param pbf_name Assay name used when `data_matrix` is a `ProBatchFeatures`
#'   object. If `NULL`, [pb_current_assay()] is used.
#' @param repeated_samples vector of sample IDs to evaluate, if \code{NULL},
#' all samples are taken into account for plotting
#' @param biospecimen_id_col column in \code{sample_annotation}
#' that defines a unique bio ID, which is usually a
#' combination of conditions or groups.
#'  Tip: if such ID is absent, but can be defined from several columns,
#'  create new \code{biospecimen_id} column
#'
#' @return dataframe with the following columns, that
#' are suggested to use for plotting in
#' \code{\link{plot_sample_corr_distribution}} as \code{plot_param}:
#' \enumerate{
#' \item \code{replicate}
#' \item \code{batch_the_same}
#' \item \code{batch_replicate}
#' \item \code{batches}
#' }
#' other columns are: \enumerate{
#' \item \code{sample_id_1} & \code{sample_id_2}, both
#' generated from \code{sample_id_col} variable
#' \item \code{correlation} - correlation of two corresponding samples
#' \item \code{batch_1} & \code{batch_2} or analogous,
#' created the same as \code{sample_id_1}
#' }
#'
#' @examples
#' data(list = c("example_sample_annotation", "example_proteome_matrix"), package = "proBatch")
#' corr_distribution <- calculate_sample_corr_distr(
#'     data_matrix = example_proteome_matrix,
#'     sample_annotation = example_sample_annotation,
#'     batch_col = "MS_batch", biospecimen_id_col = "EarTag"
#' )
#'
#' @export
#'
calculate_sample_corr_distr <- function(data_matrix, sample_annotation,
                                        repeated_samples = NULL,
                                        biospecimen_id_col = "EarTag",
                                        sample_id_col = "FullRunName",
                                        batch_col = "MS_batch",
                                        pbf_name = NULL) {
    sample_annotation_missing <- missing(sample_annotation)

    resolved <- .pb_corr_resolve_sample_input(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        pbf_name = pbf_name,
        sample_annotation_missing = sample_annotation_missing,
        require_annotation = TRUE
    )
    data_matrix <- resolved$data_matrix
    sample_annotation <- resolved$sample_annotation

    df_long <- matrix_to_long(data_matrix, sample_id_col = sample_id_col)
    df_long <- check_sample_consistency(
        sample_annotation,
        sample_id_col,
        df_long,
        batch_col,
        order_col = NULL,
        facet_col = NULL,
        merge = FALSE
    )
    data_matrix <- long_to_matrix(df_long, sample_id_col = sample_id_col)

    if (!is.null(repeated_samples)) {
        message("calculating correlation of repeated samples only")
        corr_matrix <- cor(data_matrix[, repeated_samples],
            use = "pairwise.complete.obs"
        )
    } else {
        corr_matrix <- cor(data_matrix, use = "pairwise.complete.obs")
    }
    corr_distribution <- get_sample_corr_df(
        cor_proteome = corr_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        biospecimen_id_col = biospecimen_id_col,
        batch_col = batch_col
    )
    return(corr_distribution)
}

#' @name plot_sample_corr_distribution
#' @rdname plot_sample_corr_distribution
#' @title Create violin plot of sample correlation distribution
#'
#' @description Useful to visualize within batch vs within replicate
#' vs non-related sample correlation
#'
#' @inheritParams proBatch
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames, or a
#'   `ProBatchFeatures` object. When `data_matrix` is a
#'   `ProBatchFeatures` object, `pbf_name` is used (or the current assay when
#'   `pbf_name = NULL`).
#' @param sample_annotation data frame with sample-level metadata. When
#'   `data_matrix` is a matrix, this argument is required. When
#'   `data_matrix` is a `ProBatchFeatures` object and
#'   `sample_annotation` is not provided, `as.data.frame(colData(data_matrix))`
#'   is used.
#' @param pbf_name Assay name used when `data_matrix` is a `ProBatchFeatures`
#'   object. If `NULL`, [pb_current_assay()] is used.
#' @param repeated_samples if \code{NULL}, correlation of all samples is plotted
#' @param biospecimen_id_col column in \code{sample_annotation}
#' that captures the biological sample,
#' that (possibly) was profiled several times as technical replicates.
#' Tip: if such ID is absent, but can be defined from several columns,
#' create new \code{biospecimen_id} column
#' @param plot_param columns, defined in correlation_df, which is output of
#' \code{calculate_sample_corr_distr}, specifically,  \enumerate{
#' \item \code{replicate}
#' \item \code{batch_the_same}
#' \item \code{batch_replicate}
#' \item \code{batches}
#' }
#' @param corr_distribution data frame with correlation distribution,
#' as returned by \code{calculate_sample_corr_distr}
#'
#' @return \code{ggplot} type object with violin plot
#' for each \code{plot_param}
#'
#'

#'
#' @seealso \code{\link{calculate_sample_corr_distr}},
#' \code{\link[ggplot2]{ggplot}}
NULL

#' @rdname plot_sample_corr_distribution
#'
#' @examples
#' data(list = c("example_sample_annotation", "example_proteome_matrix"), package = "proBatch")
#' sample_corr_distribution_plot <- plot_sample_corr_distribution(
#'     example_proteome_matrix,
#'     example_sample_annotation,
#'     batch_col = "MS_batch",
#'     biospecimen_id_col = "EarTag",
#'     plot_param = "batch_replicate"
#' )
#'
#' @export
#'
plot_sample_corr_distribution <- function(data_matrix, sample_annotation,
                                          repeated_samples = NULL,
                                          sample_id_col = "FullRunName",
                                          batch_col = "MS_batch",
                                          biospecimen_id_col = "EarTag",
                                          filename = NULL, width = NA, height = NA,
                                          units = c("cm", "in", "mm"),
                                          plot_title = "Sample correlation distribution",
                                          plot_param = "batch_replicate",
                                          theme = "classic",
                                          pbf_name = NULL) {
    sample_annotation_missing <- missing(sample_annotation)

    resolved <- .pb_corr_resolve_sample_input(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        pbf_name = pbf_name,
        sample_annotation_missing = sample_annotation_missing,
        require_annotation = TRUE
    )
    data_matrix <- resolved$data_matrix
    sample_annotation <- resolved$sample_annotation

    corr_distribution <- .pb_corr_distribution_from_input(
        data_matrix = data_matrix,
        builder = function(dm) {
            calculate_sample_corr_distr(
                data_matrix = dm,
                repeated_samples = repeated_samples,
                sample_annotation = sample_annotation,
                sample_id_col = sample_id_col,
                biospecimen_id_col = biospecimen_id_col,
                batch_col = batch_col
            )
        }
    )
    gg <- plot_sample_corr_distribution.corrDF(
        corr_distribution = corr_distribution,
        filename = filename,
        width = width, height = height,
        units = units,
        plot_title = plot_title,
        plot_param = plot_param,
        theme = theme
    )
    return(gg)
}

#' @rdname plot_sample_corr_distribution
#'
#' @examples
#' data(list = c("example_sample_annotation", "example_proteome_matrix"), package = "proBatch")
#' corr_distribution <- calculate_sample_corr_distr(
#'     data_matrix = example_proteome_matrix,
#'     sample_annotation = example_sample_annotation,
#'     batch_col = "MS_batch", biospecimen_id_col = "EarTag"
#' )
#' sample_corr_distribution_plot <- plot_sample_corr_distribution.corrDF(corr_distribution,
#'     plot_param = "batch_replicate"
#' )
#'
#' @export
#'
plot_sample_corr_distribution.corrDF <- function(corr_distribution,
                                                 filename = NULL, width = NA, height = NA,
                                                 units = c("cm", "in", "mm"),
                                                 plot_title = "Sample correlation distribution",
                                                 plot_param = "batch_replicate",
                                                 theme = "classic", base_size = 20) {
    gg <- ggplot(corr_distribution, aes(x = !!sym(plot_param), y = correlation)) +
        geom_violin(scale = "width") +
        geom_boxplot(width = .1) +
        theme(axis.title.x = element_blank())

    if (plot_param == "batches") {
        gg <- gg + theme(axis.text.x = element_text(angle = 90))
    }
    if (plot_param == "batch_replicate") {
        gg <- gg + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }

    .pb_finalize_corr_distribution_plot(
        gg = gg,
        corr_distribution = corr_distribution,
        plot_title = plot_title,
        theme = theme,
        base_size = base_size,
        filename = filename,
        units = units,
        width = width,
        height = height
    )
}

get_peptide_corr_df <- function(peptide_cor, peptide_annotation,
                                protein_col = "ProteinName",
                                feature_id_col = "peptide_group_label") {
    comb_to_keep <- data.frame(t(combn(colnames(peptide_cor), 2)))
    names(comb_to_keep) <- paste(feature_id_col, seq_len(2), sep = "_")

    first_feature_col <- paste(feature_id_col, "1", sep = "_")
    second_feature_col <- paste(feature_id_col, "2", sep = "_")

    corr_distribution <- peptide_cor %>%
        as.data.frame() %>%
        rownames_to_column(var = first_feature_col) %>%
        pivot_longer(
            cols = -all_of(first_feature_col),
            names_to = second_feature_col,
            values_to = "correlation",
            values_drop_na = FALSE
        ) %>%
        filter(!is.na(correlation)) %>%
        merge(comb_to_keep) %>%
        merge(
            peptide_annotation %>%
                select(all_of(c(feature_id_col, protein_col))),
            by.x = paste(feature_id_col, "1", sep = "_"),
            by.y = feature_id_col, all.x = TRUE
        ) %>%
        setnames(
            old = protein_col,
            new = paste(protein_col, 1, sep = "")
        ) %>%
        merge(
            peptide_annotation %>%
                select(all_of(c(feature_id_col, protein_col))),
            by.x = paste(feature_id_col, "2", sep = "_"),
            by.y = feature_id_col, all.x = TRUE
        ) %>%
        setnames(
            old = protein_col,
            new = paste(protein_col, 2, sep = "")
        ) %>%
        mutate(same_protein = (!!sym(paste(protein_col, "1", sep = "")) ==
            !!sym(paste(protein_col, "2", sep = "")))) %>%
        mutate(same_protein = ifelse(same_protein,
            "same protein", "different proteins"
        ))

    return(corr_distribution)
}


#' Calculate peptide correlation between and within peptides of one protein
#'
#' @inheritParams proBatch
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames, or a
#'   `ProBatchFeatures` object. When `data_matrix` is a
#'   `ProBatchFeatures` object, `pbf_name` is used (or the current assay when
#'   `pbf_name = NULL`).
#' @param peptide_annotation long format data frame with peptide ID and their
#'   corresponding protein and/or gene annotations. When `data_matrix` is a
#'   matrix, this argument is required. When `data_matrix` is a
#'   `ProBatchFeatures` object and `peptide_annotation` is not provided,
#'   rowData from the selected assay is used.
#' @param pbf_name Assay name used when `data_matrix` is a `ProBatchFeatures`
#'   object. If `NULL`, [pb_current_assay()] is used.
#'
#' @return dataframe with peptide correlation coefficients
#' that are suggested to use for plotting in
#' \code{\link{plot_peptide_corr_distribution}} as \code{plot_param}:
#'
#' @examples
#' data(list = c("example_peptide_annotation", "example_proteome_matrix"), package = "proBatch")
#' selected_genes <- c("BOVINE_A1ag", "BOVINE_FetuinB", "Cyfip1")
#' gene_filter <- example_peptide_annotation$Gene %in% selected_genes
#' peptides_ann <- example_peptide_annotation$peptide_group_label
#' selected_peptides <- peptides_ann[gene_filter]
#' matrix_test <- example_proteome_matrix[selected_peptides, ]
#' pep_annotation_sel <- example_peptide_annotation[gene_filter, ]
#' corr_distribution <- calculate_peptide_corr_distr(matrix_test,
#'     pep_annotation_sel,
#'     protein_col = "Gene"
#' )
#'
#' @export
#'
calculate_peptide_corr_distr <- function(data_matrix, peptide_annotation,
                                         protein_col = "ProteinName",
                                         feature_id_col = "peptide_group_label",
                                         pbf_name = NULL) {
    peptide_annotation_missing <- missing(peptide_annotation)

    resolved <- .pb_corr_resolve_feature_input(
        data_matrix = data_matrix,
        peptide_annotation = peptide_annotation,
        feature_id_col = feature_id_col,
        pbf_name = pbf_name,
        peptide_annotation_missing = peptide_annotation_missing,
        require_annotation = TRUE
    )
    data_matrix <- resolved$data_matrix
    peptide_annotation <- resolved$peptide_annotation

    if (is.null(peptide_annotation)) {
        stop("`peptide_annotation` must be provided.")
    }
    if (!feature_id_col %in% names(peptide_annotation)) {
        stop(sprintf("Feature ID column '%s' was not found in `peptide_annotation`.", feature_id_col))
    }
    if (!protein_col %in% names(peptide_annotation)) {
        stop(sprintf("Protein column '%s' was not found in `peptide_annotation`.", protein_col))
    }

    corr_matrix <- cor(t(data_matrix), use = "pairwise.complete.obs")
    corr_distribution <- get_peptide_corr_df(
        peptide_cor = corr_matrix,
        peptide_annotation = peptide_annotation,
        protein_col = protein_col,
        feature_id_col = feature_id_col
    )
    return(corr_distribution)
}

#' @name plot_peptide_corr_distribution
#' @rdname plot_peptide_corr_distribution
#' @title Create violin plot of peptide correlation distribution
#'
#' @description Plot distribution of peptide correlations within one
#' protein and between proteins
#'
#' @inheritParams proBatch
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames, or a
#'   `ProBatchFeatures` object. When `data_matrix` is a
#'   `ProBatchFeatures` object, `pbf_name` is used (or the current assay when
#'   `pbf_name = NULL`).
#' @param peptide_annotation long format data frame with peptide ID and their
#'   corresponding protein and/or gene annotations. When `data_matrix` is a
#'   matrix, this argument is required. When `data_matrix` is a
#'   `ProBatchFeatures` object and `peptide_annotation` is not provided,
#'   rowData from the selected assay is used.
#' @param pbf_name Assay name used when `data_matrix` is a `ProBatchFeatures`
#'   object. If `NULL`, [pb_current_assay()] is used.
#' @param corr_distribution data frame with peptide correlation distribution
#'
#' @return \code{ggplot} object (violin plot of peptide correlation)
#'
#' @seealso \code{\link{calculate_peptide_corr_distr}}, \code{\link[ggplot2]{ggplot}}
NULL

#' @rdname plot_peptide_corr_distribution
#'
#' @examples
#' data(list = c("example_peptide_annotation", "example_proteome_matrix"), package = "proBatch")
#' peptide_corr_distribution <- plot_peptide_corr_distribution(
#'     example_proteome_matrix,
#'     example_peptide_annotation,
#'     protein_col = "Gene"
#' )
#'
#' @export
#'
plot_peptide_corr_distribution <- function(data_matrix, peptide_annotation,
                                           protein_col = "ProteinName",
                                           feature_id_col = "peptide_group_label",
                                           filename = NULL, width = NA, height = NA,
                                           units = c("cm", "in", "mm"),
                                           plot_title = "Distribution of peptide correlation",
                                           theme = "classic",
                                           pbf_name = NULL) {
    peptide_annotation_missing <- missing(peptide_annotation)

    resolved <- .pb_corr_resolve_feature_input(
        data_matrix = data_matrix,
        peptide_annotation = peptide_annotation,
        feature_id_col = feature_id_col,
        pbf_name = pbf_name,
        peptide_annotation_missing = peptide_annotation_missing,
        require_annotation = TRUE
    )
    data_matrix <- resolved$data_matrix
    peptide_annotation <- resolved$peptide_annotation

    corr_distribution <- .pb_corr_distribution_from_input(
        data_matrix = data_matrix,
        builder = function(dm) {
            calculate_peptide_corr_distr(
                dm,
                peptide_annotation,
                protein_col,
                feature_id_col
            )
        },
        step_as_factor = TRUE
    )
    p <- plot_peptide_corr_distribution.corrDF(
        corr_distribution = corr_distribution,
        theme = theme,
        filename = filename,
        width = width, height = height,
        units = units,
        plot_title = plot_title
    )
    return(p)
}


#' @rdname plot_peptide_corr_distribution
#'
#' @examples
#' data(list = c("example_peptide_annotation", "example_proteome_matrix"), package = "proBatch")
#' corr_distribution <- calculate_peptide_corr_distr(
#'     example_proteome_matrix,
#'     example_peptide_annotation,
#'     protein_col = "Gene"
#' )
#' peptide_corr_distribution <- plot_peptide_corr_distribution.corrDF(corr_distribution)
#'
#' @export
#'
plot_peptide_corr_distribution.corrDF <- function(corr_distribution,
                                                  filename = NULL, width = NA, height = NA,
                                                  units = c("cm", "in", "mm"),
                                                  plot_title = "Correlation of peptides",
                                                  theme = "classic",
                                                  base_size = 20) {
    median_same_prot <- corr_distribution %>%
        filter(same_protein == "same protein") %>%
        summarize(median = median(correlation, na.rm = TRUE)) %>%
        pull(median)
    gg <- ggplot(corr_distribution, aes(
        x = same_protein,
        y = correlation
    )) +
        geom_violin(scale = "width") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
        geom_hline(
            yintercept = median_same_prot, linetype = "dotted",
            color = "tomato1"
        ) +
        geom_boxplot(width = .1) +
        xlab(NULL)

    .pb_finalize_corr_distribution_plot(
        gg = gg,
        corr_distribution = corr_distribution,
        plot_title = plot_title,
        theme = theme,
        base_size = base_size,
        filename = filename,
        units = units,
        width = width,
        height = height
    )
}
