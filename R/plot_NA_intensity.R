#' Plot missingness versus mean intensity per feature
#'
#' For each feature (row), compute the mean observed intensity and the
#' proportion of missing values, optionally stratified by a sample-level
#' grouping variable.
#' A scatter plot with a spline-smoothed trend and per-group Spearman
#' correlations is returned.
#'
#' @param x A data container.
#'   For the `ProBatchFeatures` method this must be a `ProBatchFeatures`
#'   object.
#'   The default method accepts any matrix-like input (including
#'   `SummarizedExperiment`).
#' @param sample_annotation Optional data frame with sample-level metadata.
#'   Row names (or the column specified via `sample_id_col`) must match the
#'   column names of the intensity matrix.
#'   When `x` is a `SummarizedExperiment`, `colData(x)` is used by default.
#' @param sample_id_col Optional column in `sample_annotation` providing
#'   unique sample identifiers.
#' @param color_by Column name in `sample_annotation` used to stratify
#'   features into sample groups.
#'   When `NULL`, all samples are treated as a single group.
#' @param color_scheme Colour mapping for groups.
#'   Accepts `"brewer"` (default), a named vector of colours, or a named list
#'   such as returned by [sample_annotation_to_colors()].
#' @param col_vector Optional vector of colours recycled across groups defined
#'   by `color_by`.
#' @param spline_df Degrees of freedom for the natural spline fitted via
#'   [stats::glm()] with a quasi-binomial family.
#'   Set to `0` to suppress the trend line.
#' @param point_alpha,point_size Transparency and size of individual points.
#' @param pbf_name Character scalar or vector with assay names to plot.
#'   When `NULL`, the most recent assay returned by [pb_current_assay()] is
#'   used.
#'   Only used by the `ProBatchFeatures` method.
#' @param nrow,ncol Integers controlling the layout when multiple assays are
#'   plotted.
#' @param facet_scales Scaling behaviour passed to [ggplot2::facet_wrap()]
#'   when multiple assays are plotted.
#' @param ... Currently unused.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' plot_NA_intensity(pbf, color_by = "Group")
#' }
#' @export
plot_NA_intensity <- function(x, ...) UseMethod("plot_NA_intensity")

#' @rdname plot_NA_intensity
#' @method plot_NA_intensity default
#' @export
plot_NA_intensity.default <- function(
  x,
  sample_annotation = NULL,
  sample_id_col = NULL,
  color_by = NULL,
  color_scheme = "brewer",
  col_vector = NULL,
  spline_df = 3L,
  point_alpha = 0.15,
  point_size = 0.6,
  ...
) {
    data_matrix <- x

    if (is(data_matrix, "SummarizedExperiment")) {
        if (is.null(sample_annotation)) {
            sample_annotation <- as.data.frame(colData(data_matrix))
        }
        data_matrix <- assay(data_matrix)
    }

    data_matrix <- .pb_coerce_missing_plot_matrix(data_matrix)
    if (is.null(data_matrix)) {
        return(ggplot())
    }

    stats_df <- .pb_NA_intensity_stats(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        color_by = color_by
    )
    if (is.null(stats_df) || !nrow(stats_df)) {
        warning("No finite mean intensities available for plotting.")
        return(ggplot())
    }

    .pb_plot_NA_intensity(
        stats_df = stats_df,
        color_by = color_by,
        color_scheme = .pb_resolve_missing_density_color_scheme(
            color_scheme = color_scheme,
            color_by     = color_by,
            col_vector   = col_vector
        ),
        spline_df = spline_df,
        point_alpha = point_alpha,
        point_size = point_size
    )
}

#' @rdname plot_NA_intensity
#' @method plot_NA_intensity ProBatchFeatures
#' @export
plot_NA_intensity.ProBatchFeatures <- function(
  x,
  pbf_name = NULL,
  color_by = NULL,
  sample_id_col = NULL,
  color_scheme = "brewer",
  col_vector = NULL,
  spline_df = 3L,
  point_alpha = 0.15,
  point_size = 0.6,
  nrow = NULL,
  ncol = NULL,
  facet_scales = "free_y",
  ...
) {
    object <- x
    assays <- .pb_missing_requested_assays(object, pbf_name)
    sample_annotation <- .pb_prepare_missing_heatmap_sample_annotation(
        object        = object,
        assays        = assays,
        sample_id_col = sample_id_col
    )

    collected <- .pb_collect_missing_assay_data(
        object, assays,
        function(mat, assay_nm) {
            df <- .pb_NA_intensity_stats(
                data_matrix       = mat,
                sample_annotation = sample_annotation,
                sample_id_col     = NULL,
                color_by          = color_by
            )
            if (!is.null(df) && nrow(df)) {
                df$pbf_name <- assay_nm
            }
            df
        }
    )
    keep <- collected$keep
    if (!any(keep)) {
        warning("No finite mean intensities across the requested assays.")
        return(ggplot())
    }

    combined <- do.call(rbind, collected$df_list[keep])
    combined$pbf_name <- factor(combined$pbf_name, levels = assays[keep])

    .pb_plot_NA_intensity(
        stats_df = combined,
        color_by = color_by,
        color_scheme = .pb_resolve_missing_density_color_scheme(
            color_scheme = color_scheme,
            color_by     = color_by,
            col_vector   = col_vector
        ),
        spline_df = spline_df,
        point_alpha = point_alpha,
        point_size = point_size,
        nrow = nrow,
        ncol = ncol,
        facet_scales = facet_scales
    )
}


# ---- internal helpers ------------------------------------------------

#' Compute per-feature mean intensity and proportion missing, optionally per
#' sample group.
#'
#' @return A data frame with columns `mean_intensity`, `prop_missing`,
#'   `n_samples`, and optionally `.group`.
#' @noRd
.pb_NA_intensity_stats <- function(data_matrix,
                                   sample_annotation = NULL,
                                   sample_id_col = NULL,
                                   color_by = NULL) {
    if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
    }
    if (!nrow(data_matrix) || !ncol(data_matrix)) {
        return(NULL)
    }

    if (!.pb_missing_grouping_disabled(color_by)) {
        group_info <- .pb_prepare_grouped_missing_annotation(
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            color_by = color_by,
            sample_order = colnames(data_matrix),
            empty_message = paste0(
                "Stratified missingness-intensity plots require one or ",
                "more metadata columns in `color_by`."
            )
        )
        df_list <- lapply(seq_along(group_info$group_indices), function(idx) {
            cols <- group_info$group_indices[[idx]]
            sub_df <- .pb_NA_intensity_stats_ungrouped(
                data_matrix[, cols, drop = FALSE]
            )
            if (!is.null(sub_df) && nrow(sub_df)) {
                sub_df$.group <- group_info$annotation$.group_key[[idx]]
            }
            sub_df
        })
        df_list <- Filter(function(d) !is.null(d) && nrow(d) > 0, df_list)
        if (!length(df_list)) {
            return(NULL)
        }
        return(do.call(rbind, df_list))
    }

    .pb_NA_intensity_stats_ungrouped(data_matrix)
}

#' @noRd
.pb_NA_intensity_stats_ungrouped <- function(data_matrix) {
    n_samples <- as.integer(ncol(data_matrix))
    n_valid <- rowSums(!is.na(data_matrix))
    has_obs <- n_valid > 0L

    if (!any(has_obs)) {
        return(NULL)
    }

    means <- rowMeans(data_matrix[has_obs, , drop = FALSE], na.rm = TRUE)
    props <- 1 - n_valid[has_obs] / n_samples

    df <- data.frame(
        mean_intensity = means,
        prop_missing = props,
        n_samples = n_samples,
        stringsAsFactors = FALSE
    )
    df[is.finite(df$mean_intensity), , drop = FALSE]
}


#' Build the ggplot for missingness vs intensity.
#' @noRd
.pb_plot_NA_intensity <- function(stats_df,
                                  color_by,
                                  color_scheme,
                                  spline_df = 3L,
                                  point_alpha = 0.15,
                                  point_size = 0.6,
                                  nrow = NULL,
                                  ncol = NULL,
                                  facet_scales = "free_y") {
    has_group <- ".group" %in% names(stats_df)

    if (has_group) {
        n_groups <- length(unique(stats_df$.group))
        p <- ggplot(
            stats_df,
            aes(
                x     = .data$mean_intensity,
                y     = .data$prop_missing,
                color = .data$.group
            )
        )
        p <- color_discrete(
            color_scheme = color_scheme,
            batch_col = ".group",
            n_batches = n_groups,
            fill_or_color = "color",
            gg = p
        )
    } else {
        p <- ggplot(
            stats_df,
            aes(x = .data$mean_intensity, y = .data$prop_missing)
        )
    }

    p <- p +
        geom_point(alpha = point_alpha, size = point_size) +
        scale_y_continuous(
            limits = c(0, 1),
            labels = scales::label_percent(accuracy = 1),
            expand = expansion(mult = c(0, 0.05))
        ) +
        labs(
            x     = "Mean log2 intensity (non-missing values)",
            y     = "% missing values",
            title = "Missingness vs intensity per feature",
            color = if (has_group) .pb_missing_density_group_title(color_by)
        )

    if (spline_df > 0L) {
        smooth_formula <- y ~ splines::ns(x, df = spline_df)
        p <- p + geom_smooth(
            aes(weight = .data$n_samples),
            method = "glm",
            method.args = list(family = quasibinomial()),
            formula = smooth_formula,
            se = FALSE,
            linewidth = 1.2
        )
    }

    ## Spearman correlation annotation
    cor_df <- .pb_NA_intensity_cor_labels(stats_df, has_group)
    if (!is.null(cor_df) && nrow(cor_df)) {
        if (has_group) {
            p <- p + geom_text(
                data = cor_df,
                aes(
                    x     = .data$x,
                    y     = .data$y,
                    label = .data$label,
                    color = .data$.group
                ),
                hjust = 1.1,
                fontface = "bold",
                size = 4,
                inherit.aes = FALSE,
                show.legend = FALSE
            )
        } else {
            p <- p + geom_text(
                data = cor_df,
                aes(x = .data$x, y = .data$y, label = .data$label),
                hjust = 1.1,
                fontface = "bold",
                size = 4,
                inherit.aes = FALSE,
                show.legend = FALSE
            )
        }
    }

    if ("pbf_name" %in% names(stats_df) &&
        length(unique(stats_df$pbf_name)) > 1L) {
        layout <- .pb_missing_layout(
            length(unique(stats_df$pbf_name)),
            nrow = nrow, ncol = ncol
        )
        p <- p + facet_wrap(
            ~pbf_name,
            nrow   = layout$nrow,
            ncol   = layout$ncol,
            scales = facet_scales
        )
    }

    p + .pb_missing_density_theme() +
        theme(legend.position = "right")
}


#' Compute per-group Spearman correlations for the annotation layer.
#' @noRd
.pb_NA_intensity_cor_labels <- function(stats_df, has_group) {
    if (has_group) {
        groups <- unique(stats_df$.group)
        rhos <- vapply(groups, function(g) {
            sub <- stats_df[stats_df$.group == g, , drop = FALSE]
            if (nrow(sub) < 3L) {
                return(NA_real_)
            }
            cor(sub$mean_intensity, sub$prop_missing,
                method = "spearman", use = "complete.obs"
            )
        }, numeric(1))
        data.frame(
            .group = groups,
            label = paste0("\u03C1 = ", round(rhos, 2)),
            x = Inf,
            y = seq(0.95, by = -0.06, length.out = length(groups)),
            stringsAsFactors = FALSE
        )
    } else {
        if (nrow(stats_df) < 3L) {
            return(NULL)
        }
        rho <- cor(stats_df$mean_intensity, stats_df$prop_missing,
            method = "spearman", use = "complete.obs"
        )
        data.frame(
            label = paste0("\u03C1 = ", round(rho, 2)),
            x = Inf,
            y = 0.95,
            stringsAsFactors = FALSE
        )
    }
}
