# Helper: compute CV (%) over arbitrary groups
compute_cv <- function(data, measure_col, group_vars, cv_name) {
    data %>%
        group_by(across(all_of(group_vars))) %>%
        mutate(!!sym(cv_name) := 100 * sd(.data[[measure_col]], na.rm = TRUE) /
            mean(.data[[measure_col]], na.rm = TRUE)) %>%
        ungroup()
}

#' Calculate CV distribution for each feature
#'
#' @inheritParams proBatch
#' @inheritParams transform_raw_data
#' @param biospecimen_id_col column in \code{sample_annotation}
#' that defines a unique bio ID, which is usually a
#' combination of conditions or groups.
#'  Tip: if such ID is absent, but can be defined from several columns,
#'  create new \code{biospecimen_id} column
#' @param unlog (logical) whether to reverse log transformation of the original data
#'
#' @return data frame with Total CV for each feature & (optionally) per-batch CV
#' @export
#'
#' @examples
#' data(list = c("example_sample_annotation", "example_proteome"), package = "proBatch")
#' CV_df <- calculate_feature_CV(example_proteome,
#'     sample_annotation = example_sample_annotation,
#'     measure_col = "Intensity",
#'     batch_col = "MS_batch"
#' )
calculate_feature_CV <- function(df_long, sample_annotation = NULL,
                                 feature_id_col = "peptide_group_label",
                                 sample_id_col = "FullRunName",
                                 measure_col = "Intensity", batch_col = NULL,
                                 biospecimen_id_col = NULL,
                                 unlog = TRUE, log_base = 2, offset = 0) {
    # Handle sample_annotation and check for sample_id_col
    if (!is.null(sample_annotation) && is.null(sample_id_col)) {
        message("sample_id_col is not specified, using FullRunName as default")
        sample_id_col <- "FullRunName"
    }
    if (!is.null(sample_annotation)) {
        df_long <- check_sample_consistency(sample_annotation, sample_id_col, df_long)
    }
    # Biospecimen ID fallback
    if (is.null(biospecimen_id_col)) {
        warning("considering all samples as replicates!")
        biospecimen_id_col <- "biospecimen_id"
        df_long[[biospecimen_id_col]] <- "replication"
    } else {
        if (!(biospecimen_id_col %in% names(df_long))) {
            stop("biospecimen ID, indicating replicates, is not in the data (df_long or sample_annotation)")
        }
    }

    # Optional unlog
    if (unlog) {
        message("reversing log-transformation for CV calculation!")
        df_long <- unlog_df(df_long, log_base = log_base, offset = offset, measure_col = measure_col)
    }

    # Filter out features with ≤ 2 total measurements
    if (!is.null(batch_col)) {
        df_long <- df_long %>%
            group_by(!!sym(feature_id_col), !!sym(batch_col), !!sym(biospecimen_id_col)) %>%
            mutate(n_total = sum(!is.na(!!sym(measure_col)))) %>%
            ungroup()
    } else {
        df_long <- df_long %>%
            group_by(!!sym(feature_id_col), !!sym(biospecimen_id_col)) %>%
            mutate(n_total = sum(!is.na(!!sym(measure_col)))) %>%
            ungroup()
    }
    if (any(df_long$n_total <= 2)) {
        # how many features have 2 or less measurements? unique peptides
        n_peptides <- df_long %>%
            filter(n_total <= 2) %>%
            distinct(!!sym(feature_id_col)) %>%
            nrow()
        message(paste0("Cannot calculate CV for ", n_peptides, " peptides with 2 or less measurements, removing those peptides"))
        df_long <- df_long %>%
            filter(n_total > 2)
    }
    df_long$n_total <- NULL

    # Compute CV
    # Detect presence of Step and build grouping vectors
    has_step <- "Step" %in% names(df_long)
    base_group <- feature_id_col
    step_group <- if (has_step) "Step" else NULL

    perbatch_groups <- c(base_group, batch_col, step_group) %>% compact()
    total_groups <- c(base_group, step_group) %>% compact()

    # Compute per‐batch CV (if batch_col given)
    if (!is.null(batch_col)) {
        df_long <- compute_cv(df_long, measure_col, perbatch_groups, "CV_perBatch")
    } else {
        warning("`batch_col` not specified – skipping per‐batch CV, only total CV will be calculated.")
    }
    # Compute total CV
    df_long <- compute_cv(df_long, measure_col, total_groups, "CV_total")

    # Final select + distinct
    select_cols <- c(feature_id_col, if (has_step) "Step", "CV_total", if (!is.null(batch_col)) "CV_perBatch")
    CV_df <- df_long %>%
        select(all_of(select_cols)) %>%
        distinct()

    return(CV_df)
}

#' Plot the distribution (boxplots) of per-batch per-step CV of features
#'
#' @inheritParams proBatch
#'
#' @param CV_df data frame with Total CV for each feature & (optionally) per-batch CV
#' @param log_y_scale (logical) whether to display the CV on log-scale
#'
#' @return ggplot object
#' @export
plot_CV_distr.df <- function(CV_df,
                             plot_title = NULL,
                             filename = NULL, theme = "classic", log_y_scale = TRUE) {
    if ("Step" %in% names(CV_df)) {
        gg <- ggplot(CV_df, aes(x = !!sym("Step"), y = !!sym("CV_total"))) +
            geom_boxplot()
    } else {
        gg <- ggplot(CV_df, aes(y = CV_total)) +
            geom_boxplot()
    }
    if (!is.null(plot_title)) {
        gg <- gg + ggtitle(plot_title)
    }
    if (theme == "classic") {
        gg <- gg + theme_classic()
    }
    if (log_y_scale) {
        gg <- gg + scale_y_log10()
    }
    if (!is.null(filename)) {
        ggsave(filename = filename, plot = gg)
    }
    return(gg)
}

#' Plot CV distribution to compare various steps of the analysis
#'
#' @inheritParams proBatch
#' @inheritParams transform_raw_data
#' @param df_long as in \code{df_long} for the rest of the package, but, when it
#' has entries for intensity, represented in \code{measure_col} for several steps,
#' e.g. raw, normalized, batch corrected data, as seen in column \code{Step}, then
#' multi-step CV comparison can be carried out.
#' @param biospecimen_id_col column in \code{sample_annotation}
#' that defines a unique bio ID, which is usually a
#' combination of conditions or groups.
#'  Tip: if such ID is absent, but can be defined from several columns,
#'  create new \code{biospecimen_id} column
#' @param unlog (logical) whether to reverse log transformation of the original data
#'
#' @return \code{ggplot} object with the boxplot of CVs on one or several steps
#' @export
#'
#' @examples
#' data(list = c("example_sample_annotation", "example_proteome"), package = "proBatch")
#' CV_plot <- plot_CV_distr(example_proteome,
#'     sample_annotation = example_sample_annotation,
#'     measure_col = "Intensity", batch_col = "MS_batch",
#'     plot_title = NULL, filename = NULL, theme = "classic"
#' )
plot_CV_distr <- function(df_long, sample_annotation = NULL,
                          feature_id_col = "peptide_group_label",
                          sample_id_col = "FullRunName",
                          measure_col = "Intensity",
                          biospecimen_id_col = "EarTag",
                          batch_col = NULL,
                          unlog = TRUE,
                          log_base = 2,
                          offset = 1,
                          plot_title = NULL,
                          filename = NULL, theme = "classic") {
    CV_df <- calculate_feature_CV(
        df_long = df_long,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col,
        batch_col = batch_col,
        biospecimen_id_col = biospecimen_id_col,
        unlog = unlog,
        log_base = log_base,
        offset = offset
    )
    # keep only finite CV values - check, message, and filter
    if (any(!is.finite(CV_df$CV_total))) {
        message(
            "Some CV values are not finite, filtering them out - number of such features: ",
            sum(!is.finite(CV_df$CV_total))
        )
        CV_df <- CV_df %>%
            filter(is.finite(CV_total))
    }
    # Check if CV_df is empty
    if (nrow(CV_df) == 0) {
        stop("No features with sufficient measurements to calculate CV.")
    }

    gg <- plot_CV_distr.df(
        CV_df,
        plot_title = plot_title, filename = filename,
        theme = theme
    )
    return(gg)
}
