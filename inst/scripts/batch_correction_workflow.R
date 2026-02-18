#!/usr/bin/env Rscript
#' Automated Batch Effect Correction and Comprehensive Diagnostics Workflow
#' 
#' This script accepts a ProBatchFeatures object, performs batch effect correction,
#' and generates comprehensive diagnostic plots and metrics.
#' 
#' @description
#' The workflow includes:
#' - Baseline diagnostics on raw data
#' - Batch effect correction (multiple methods supported)
#' - Post-correction diagnostics
#' - Comparative analysis (raw vs corrected)
#' - Systematic output saving with meaningful structure
#' 
#' All calculations use try() for robustness - failures are logged but don't halt the workflow.

# ==============================================================================
# CONFIGURATION PARAMETERS - EDIT THIS SECTION
# ==============================================================================

# --- Input data ---
# Provide either:
# Option 1: pbf_object directly (if running interactively)
pbf_object <- NULL  # ProBatchFeatures object (set to NULL if loading from file)

# Option 2: path to saved RDS file containing ProBatchFeatures object
pbf_file_path <- NULL  # Example: "/path/to/pbf_object.rds"

# Assay prefix in pbf_object to use (if pbf_object is provided)
assay_prefix <- "PGs"

# --- Annotation data ---
# If not already in pbf_object, provide paths to annotation tables
sample_annotation_path <- NULL  # Path to sample metadata CSV/RDS
peptide_annotation_path <- NULL  # Path to feature/peptide annotation CSV/RDS

# --- Column identifiers ---
sample_id_col <- "FullRunName"          # Sample ID column name
feature_id_col <- "peptide_group_label" # Feature ID column name
measure_col <- "Intensity"              # Intensity column name for long format

# --- Technical factors (batch variables) ---
batch_col <- "MS_batch"                 # Primary batch column
technical_factors <- c("MS_batch", "digestion_batch")  # All technical factors for PVCA

# --- Biological factors (condition/covariates) ---
condition_col <- "Diet"                 # Primary biological condition
biological_factors <- c("Diet", "Sex", "Strain")  # All biological factors for PVCA
biospecimen_id_col <- "EarTag"         # Biological replicate identifier (for CV, correlation)

# --- Covariates for correction ---
correction_covariates <- c("Diet", "Sex")  # Protected biological factors during correction

# --- Other metadata columns ---
order_col <- "order"                    # Run order column (for drift diagnostics)
datetime_col <- "DateTime"              # Date/time column (optional)

# --- Factors for visualization ---
factors_to_plot <- c("MS_batch", "Diet", "Sex", "Strain")  # For heatmaps, clustering
factor_columns_colors <- c("MS_batch", "Diet", "Sex", "Strain", "digestion_batch")
numeric_columns_colors <- c("DateTime", "order")

# --- Feature-level selectors (for example plots) ---
# Provide specific feature/protein names, or set to NULL for auto-selection
example_protein_name <- NULL      # Example: "Haao" or NULL for auto
example_spike_in_name <- NULL     # Example: "BOVINE_A1ag" or NULL for auto
example_irt_pattern <- "iRT"      # Pattern to identify iRT peptides
example_feature_id <- NULL        # Specific feature ID or NULL for auto

# --- Batch correction method ---
correction_method <- "limmaRBE"   # Options: "limmaRBE", "combat", "centerMean", "centerMedian", etc.
# See ?pb_transform and ?correct_batch_effects for available methods

# --- Output configuration ---
output_base_dir <- file.path(tempdir(), "batch_correction_results")  # Base output directory
results_prefix <- "BEC"           # Prefix for result files
save_plots <- TRUE                # Save plots to files
save_metrics <- TRUE              # Save metric tables to CSV
save_objects <- TRUE              # Save R objects (pbf, matrices) to RDS

# --- Plot settings ---
plot_format <- "png"              # Options: "png", "pdf", "svg"
plot_width <- 10                  # inches
plot_height <- 6                  # inches
plot_dpi <- 300                   # resolution for raster formats

# --- Computational settings ---
set_seed <- 42                    # Random seed for reproducibility
n_pcs_embeddings <- 10            # Number of PCs for t-SNE/UMAP
perplexity_tsne <- 10             # t-SNE perplexity
n_neighbors_umap <- 10            # UMAP n_neighbors
fill_missing <- -1                # Value for missing data in PVCA/variance partition (-1 or 0)
variance_threshold <- 0.05        # PVCA/variance partition threshold

# --- Sample subsetting (for speed in large datasets) ---
samples_for_corr_heatmap <- NULL  # Vector of sample IDs or NULL for all (or first N)
n_samples_corr_heatmap <- 20      # If samples_for_corr_heatmap is NULL, use first N

# --- Caching ---
cache_expensive <- TRUE           # Cache expensive operations (embeddings, variance partition)

# ==============================================================================
# SETUP AND VALIDATION
# ==============================================================================

# Load required packages
required_packages <- c("proBatch", "dplyr", "tibble", "ggplot2", "gridExtra")
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is required but not installed.", pkg))
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Set random seed
if (!is.null(set_seed)) {
    set.seed(set_seed)
}


# Logging function
log_msg <- function(msg, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
}

# Create output directory structure
if (save_plots || save_metrics || save_objects) {
    output_dirs <- list(
        base = output_base_dir,
        raw = file.path(output_base_dir, "01_raw_diagnostics"),
        corrected = file.path(output_base_dir, "02_corrected_diagnostics"),
        comparison = file.path(output_base_dir, "03_raw_vs_corrected"),
        plots_raw = file.path(output_base_dir, "01_raw_diagnostics", "plots"),
        plots_corr = file.path(output_base_dir, "02_corrected_diagnostics", "plots"),
        plots_comp = file.path(output_base_dir, "03_raw_vs_corrected", "plots"),
        metrics_raw = file.path(output_base_dir, "01_raw_diagnostics", "metrics"),
        metrics_corr = file.path(output_base_dir, "02_corrected_diagnostics", "metrics"),
        metrics_comp = file.path(output_base_dir, "03_raw_vs_corrected", "metrics"),
        objects = file.path(output_base_dir, "04_objects"),
        pvca_raw = file.path(output_base_dir, "01_raw_diagnostics", "pvca_results"),
        pvca_corr = file.path(output_base_dir, "02_corrected_diagnostics", "pvca_results")
    )
    
    for (dir_path in output_dirs) {
        if (!dir.exists(dir_path)) {
            dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
        }
    }

    log_msg(sprintf("Output directories created under: %s", output_base_dir))
}

log_msg("Starting batch correction workflow")

# ==============================================================================
# LOAD DATA
# ==============================================================================

log_msg("Loading input data")

# Load pbf_object
if (is.null(pbf_object)) {
    if (!is.null(pbf_file_path) && file.exists(pbf_file_path)) {
        log_msg(sprintf("Loading ProBatchFeatures object from: %s", pbf_file_path))
        pbf_object <- readRDS(pbf_file_path)
    } else {
        stop("No valid pbf_object or pbf_file_path provided.")
    }
}

# Validate pbf_object
if (!inherits(pbf_object, "ProBatchFeatures")) {
    stop("Input object must be of class 'ProBatchFeatures'.")
}

# Extract or load annotations
sample_annotation <- colData(pbf_object) %>% as.data.frame() %>%
    dplyr::select(
        all_of(sample_id_col), 
        all_of(technical_factors), 
        all_of(biological_factors)
    )
rownames(sample_annotation) <- sample_annotation[[sample_id_col]]

if (is.null(sample_annotation) && !is.null(sample_annotation_path)) {
    log_msg(sprintf("Loading sample annotation from: %s", sample_annotation_path))
    if (grepl("\\.rds$", sample_annotation_path, ignore.case = TRUE)) {
        sample_annotation <- readRDS(sample_annotation_path)
    } else {
        sample_annotation <- read.csv(sample_annotation_path, stringsAsFactors = FALSE)
    }
}

peptide_annotation <- NULL
if (!is.null(peptide_annotation_path) && file.exists(peptide_annotation_path)) {
    log_msg(sprintf("Loading feature annotation from: %s", peptide_annotation_path))
    if (grepl("\\.rds$", peptide_annotation_path, ignore.case = TRUE)) {
        peptide_annotation <- readRDS(peptide_annotation_path)
    } else {
        peptide_annotation <- read.csv(peptide_annotation_path, stringsAsFactors = FALSE)
    }
}

# Ensure pbf_object is log-transformed for diagnostics
current_assay <- proBatch::pb_current_assay(pbf_object)
if (!grepl("log", current_assay, ignore.case = TRUE)) {
    log_msg("Applying log2 transformation")
    pbf_object <- proBatch::pb_transform(
        pbf_object,
        from = current_assay,
        steps = "log2",
        store_fast_steps = TRUE
    )
} else {
    log_msg("Data already log-transformed")
}

######################################################################################
######################################################################################
######################################################################################

# Extract matrices and long format
data_matrix_raw_log <- proBatch::pb_assay_matrix(
    pbf_object, 
    assay = paste0(assay_prefix, "::log2_on_raw")
)
df_long_raw <- proBatch::pb_as_long(
    pbf_object,
    feature_id_col = feature_id_col,
    sample_id_col = sample_id_col,
    measure_col = measure_col,
    pbf_name = paste0(assay_prefix, "::raw")
)

# Auto-select example features if not specified
if (is.null(example_feature_id)) {
    example_feature_id <- rownames(data_matrix_raw_log)[1]
}

if (!is.null(peptide_annotation)) {
    if (is.null(example_protein_name) && "Gene" %in% names(peptide_annotation)) {
        example_protein_name <- unique(peptide_annotation$Gene)[1]
    }
    if (is.null(example_spike_in_name) && "Gene" %in% names(peptide_annotation)) {
        spike_candidates <- grep("BOVINE|SPIKE", peptide_annotation$Gene, 
                                value = TRUE, ignore.case = TRUE)
        if (length(spike_candidates) > 0) {
            example_spike_in_name <- spike_candidates[1]
        }
    }
}

# Prepare sample subset for correlation heatmaps
if (is.null(samples_for_corr_heatmap)) {
    samples_for_corr_heatmap <- head(colnames(data_matrix_raw_log), n_samples_corr_heatmap)
}


######################################################################################
######################################################################################
######################################################################################

# Generate color scheme
color_list <- proBatch::sample_annotation_to_colors(
    pbf_log,
    factor_columns = factor_columns_colors,
    numeric_columns = numeric_columns_colors
)

log_msg(sprintf("Data dimensions: %d features x %d samples", 
               nrow(data_matrix_raw_log), ncol(data_matrix_raw_log)))

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Safe function wrapper with try()
safe_try <- function(expr, error_msg = "Error occurred", return_on_error = NULL) {
    result <- try(expr, silent = TRUE)
    if (inherits(result, "try-error")) {
        log_msg(sprintf("%s: %s", error_msg, as.character(result)), level = "WARN")
        return(return_on_error)
    }
    return(result)
}

# Save plot helper
save_plot_helper <- function(plot_obj, filename, subdir = "plots", width = plot_width, 
                            height = plot_height, format = plot_format, dpi = plot_dpi) {
    if (!save_plots) return(invisible(NULL))
    
    filepath <- file.path(subdir, sprintf("%s.%s", filename, format))
    
    safe_try({
        if (format == "pdf") {
            pdf(filepath, width = width, height = height)
        } else if (format == "png") {
            png(filepath, width = width, height = height, units = "in", res = dpi)
        } else if (format == "svg") {
            svg(filepath, width = width, height = height)
        }
        print(plot_obj)
        dev.off()
        log_msg(sprintf("Saved plot: %s", filepath))
    }, error_msg = sprintf("Failed to save plot: %s", filepath))
}

# Save metric helper
save_metric_helper <- function(data_obj, filename, subdir = "metrics") {
    if (!save_metrics) return(invisible(NULL))
    
    filepath <- file.path(subdir, sprintf("%s.csv", filename))
    
    safe_try({
        if (is.data.frame(data_obj) || is.matrix(data_obj)) {
            write.csv(data_obj, filepath, row.names = TRUE)
            log_msg(sprintf("Saved metric: %s", filepath))
        } else {
            safe_try({
                metric_flatten <- pb_flatten_design_check(data_obj)
                write.csv(metric_flatten, filepath, row.names = FALSE)
                log_msg(sprintf("Saved metric (flattened): %s", filepath))
            }, error_msg = sprintf("Failed to save flattened metric: %s", filepath))
        }
    }, error_msg = sprintf("Failed to save metric: %s", filepath))
}


# Safe classification metrics
safe_classification_metrics <- function(
    x, sample_annotation, 
    known_col = batch_col,
    fill_the_missing = fill_missing) {
    if (is.null(sample_annotation) || !all(known_col %in% names(sample_annotation))) {
        return(NULL)
    }
    for (col in known_col) {
        n_levels <- length(unique(sample_annotation[[col]][!is.na(sample_annotation[[col]])]))
        if (n_levels < 2L) {
            log_msg(sprintf("Not enough levels in column '%s' for classification metrics. Skipping.", col), level = "WARN")
            known_col <- setdiff(known_col, col)
        }
    }
    if (length(known_col) == 0) {
        log_msg("No valid known columns for classification metrics. Skipping.", level = "WARN")
        return(NULL)
    }

    safe_try(
        proBatch::calculate_classification_metrics(
            x, known_col = known_col,
            fill_the_missing = fill_the_missing
        ),
        error_msg = "Classification metrics calculation failed",
        return_on_error = NULL
    )
}

# Flatten design check output for saving
pb_flatten_design_check <- function(x, sep = ", ") {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  empty <- data.frame(section=character(), key=character(), value=character(), stringsAsFactors=FALSE)
  if (is.null(x)) return(empty)

  if (!is.list(x)) return(data.frame(section="value", key="", value=paste(x, collapse=sep), stringsAsFactors=FALSE))

  summ <- x$summary %||% x[setdiff(names(x), c("errors","warnings","summary"))]
  summ <- if (length(summ)) rapply(summ, function(v) if (length(v)) paste(v, collapse=sep) else NA_character_, how="replace") else list()
  flat <- if (length(summ)) unlist(summ, recursive=TRUE, use.names=TRUE) else character(0)

  nm <- names(flat); if (is.null(nm)) nm <- rep("", length(flat))
  parts <- strsplit(nm, "\\.", perl = TRUE)
  sec <- vapply(parts, function(p) if (length(p) && nzchar(p[1])) p[1] else "summary", character(1))
  key <- vapply(parts, function(p) if (length(p) > 1) paste(p[-1], collapse=".") else "", character(1))

  errs  <- as.character(unlist(x$errors   %||% character(0), use.names=FALSE))
  warns <- as.character(unlist(x$warnings %||% character(0), use.names=FALSE))

  out <- do.call(rbind, Filter(Negate(is.null), list(
    if (length(errs))  data.frame(section="errors",   key=as.character(seq_along(errs)),  value=errs,  stringsAsFactors=FALSE) else NULL,
    if (length(warns)) data.frame(section="warnings", key=as.character(seq_along(warns)), value=warns, stringsAsFactors=FALSE) else NULL,
    if (length(flat))  data.frame(section=sec, key=key, value=as.character(unname(flat)), stringsAsFactors=FALSE) else NULL
  )))
  if (is.null(out)) empty else out
}



# ==============================================================================
# BASELINE DIAGNOSTICS (RAW DATA)
# ==============================================================================

log_msg("=" %+% "=" %+% "=" %+% " BASELINE DIAGNOSTICS (RAW DATA) " %+% "=" %+% "=" %+% "=")

# --- Metadata and design checks ---
log_msg("Running metadata and design diagnostics")

meta_summary_raw <- safe_try(
    proBatch::metadata_column_summary(pbf_object),
    error_msg = "metadata_column_summary failed"
)
save_metric_helper(meta_summary_raw, "meta_summary_raw", output_dirs$metrics_raw)

dup_cols_raw <- safe_try(
    proBatch::find_duplicated_columns(pbf_object),
    error_msg = "find_duplicated_columns failed"
)
if (!is.null(dup_cols_raw) && NROW(dup_cols_raw) > 0) { 
    save_metric_helper(dup_cols_raw, "duplicated_columns_raw", output_dirs$metrics_raw)
} else {
    log_msg("No duplicated columns found in raw data")
}

design_check_raw <- safe_try(
    proBatch::validate_batch_design(
        sample_annotation,
        batch_col = batch_col,
        condition_col = condition_col,
        covariates = setdiff(biological_factors, condition_col),
        sample_id_col = sample_id_col,
        strict = FALSE
    ),
    error_msg = "validate_batch_design failed"
)
save_metric_helper(design_check_raw, "design_validation_raw", output_dirs$metrics_raw)

design_summary_raw <- safe_try(
    proBatch::summarize_design(
        sample_annotation,
        batch_col = batch_col,
        condition_col = condition_col,
        covariates = setdiff(biological_factors, condition_col)
    ),
    error_msg = "summarize_design failed"
)
save_metric_helper(design_summary_raw, "design_summary_raw", output_dirs$metrics_raw)

nested_raw <- safe_try(
    proBatch::detect_nested_batches(
        sample_annotation,
        batch_cols = c(batch_col, technical_factors)
    ),
    error_msg = "detect_nested_batches failed"
)
save_metric_helper(nested_raw, "nested_batches_raw", output_dirs$metrics_raw)

# --- Missingness diagnostics ---
log_msg("Generating missingness diagnostics")

p_na_heatmap_raw <- safe_try({
    proBatch::plot_NA_heatmap(
        pbf_object,
        color_by = batch_col,
        label_by = sample_id_col,
        cluster_samples = FALSE,
        cluster_features = TRUE,
        show_rownames = FALSE
    )
}, error_msg = "plot_NA_heatmap (raw) failed")
save_plot_helper(p_na_heatmap_raw, "NA_heatmap_raw", output_dirs$plots_raw)

p_na_density_raw <- safe_try(
    proBatch::plot_NA_density(
        pbf_object, 
        pbf_name = paste0(assay_prefix, "::log2_on_raw") 
    ) + ggtitle("Raw: NA density"),
    error_msg = "plot_NA_density (raw) failed"
)
save_plot_helper(p_na_density_raw, "NA_density_raw", output_dirs$plots_raw)

p_na_freq_raw <- safe_try(
    proBatch::plot_NA_frequency(
        pbf_object, 
        pbf_name = paste0(assay_prefix, "::raw") 
    ) + ggtitle("Raw: NA frequency"),
    error_msg = "plot_NA_frequency (raw) failed"
)
save_plot_helper(p_na_freq_raw, "NA_frequency_raw", output_dirs$plots_raw)

# --- Sample distributions ---
log_msg("Generating sample distribution diagnostics")

p_sample_mean_raw <- safe_try({
    proBatch::plot_sample_mean(
        pbf_object,
        order_col = order_col,
        batch_col = batch_col,
        color_by_batch = TRUE,
        color_scheme = color_list,
        base_size = 15
    )
}, error_msg = "plot_sample_mean (raw) failed")
save_plot_helper(p_sample_mean_raw, "sample_mean_raw", output_dirs$plots_raw)

p_boxplot_raw <- safe_try({
    proBatch::plot_boxplot(
        pbf_object,
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        color_by_batch = TRUE,
        color_scheme = color_list,
        base_size = 10
    )
}, error_msg = "plot_boxplot (raw) failed")
save_plot_helper(p_boxplot_raw, "boxplot_raw", output_dirs$plots_raw)

# --- Heatmaps and clustering ---
log_msg("Generating heatmaps and clustering diagnostics")

p_heatmap_diag_raw <- safe_try({
    proBatch::plot_heatmap_diagnostic(
        pbf_object,
        factors_to_plot = factors_to_plot,
        color_list = color_list,
        show_rownames = FALSE,
        show_colnames = FALSE,
        fill_the_missing = fill_missing
    )
}, error_msg = "plot_heatmap_diagnostic (raw) failed")
save_plot_helper(p_heatmap_diag_raw, "heatmap_diagnostic_raw", output_dirs$plots_raw)

p_hierarchical_raw <- safe_try({
    proBatch::plot_hierarchical_clustering(
        pbf_object,
        factors_to_plot = factors_to_plot,
        color_list = color_list,
        label_samples = FALSE,
        fill_the_missing = fill_missing
    )
}, error_msg = "plot_hierarchical_clustering (raw) failed")
save_plot_helper(p_hierarchical_raw, "hierarchical_clustering_raw", output_dirs$plots_raw)

# --- Embeddings ---
log_msg("Generating embeddings (PCA, t-SNE, UMAP)")

pca_raw <- safe_try({
    proBatch::plot_PCA(
        pbf_object,
        color_by = batch_col,
        shape_by = condition_col,
        fill_the_missing = fill_missing
    )
}, error_msg = "plot_PCA (raw) failed")
save_plot_helper(pca_raw, "PCA_raw", output_dirs$plots_raw)

tsne_raw <- safe_try({
    proBatch::plot_TSNE(
        pbf_object,
        color_by = batch_col,
        shape_by = condition_col,
        perplexity = perplexity_tsne,
        max_iter = 300,
        initial_dims = n_pcs_embeddings,
        fill_the_missing = fill_missing,
        point_size = 4
    )
}, error_msg = "plot_TSNE (raw) failed")
save_plot_helper(tsne_raw, "tSNE_raw", output_dirs$plots_raw)

umap_raw <- safe_try({
    proBatch::plot_UMAP(
        pbf_object,
        color_by = batch_col,
        shape_by = condition_col,
        n_neighbors = n_neighbors_umap,
        min_dist = 0.3,
        random_state = set_seed,
        fill_the_missing = fill_missing,
        point_size = 4
    )
}, error_msg = "plot_UMAP (raw) failed")
save_plot_helper(umap_raw, "UMAP_raw", output_dirs$plots_raw)

# --- PVCA ---
log_msg("Calculating PVCA (variance decomposition)")

pvca_raw <- safe_try({
    proBatch::calculate_PVCA(
        pbf_object,
        factors_for_PVCA = c(technical_factors, biological_factors),
        fill_the_missing = fill_missing
    )
}, error_msg = "calculate_PVCA (raw) failed")
save_metric_helper(pvca_raw, "pvca_raw", output_dirs$metrics_raw)

p_pvca_raw <- safe_try({
    proBatch::plot_PVCA(
        pbf_object,
        technical_factors = technical_factors,
        biological_factors = biological_factors,
        fill_the_missing = fill_missing,
        variance_threshold = (ifelse(
            length(c(technical_factors, biological_factors)) > 5, 0.01, 0)),
        base_size = 10
    )
}, error_msg = "plot_PVCA (raw) failed")
save_plot_helper(p_pvca_raw, "PVCA_raw", output_dirs$plots_raw)

pvca_df_raw <- safe_try({
    proBatch::prepare_PVCA_df(
        pbf_object,
        technical_factors = technical_factors,
        biological_factors = biological_factors,
        fill_the_missing = fill_missing,
        variance_threshold = variance_threshold,
        path_to_save_results = output_dirs$pvca_raw
    )
}, error_msg = "prepare_PVCA_df (raw) failed")
save_metric_helper(pvca_df_raw, "pvca_prepared_raw", output_dirs$metrics_raw)

# --- Variance partition ---
log_msg("Calculating variance partition")

vp_raw <- safe_try({
    proBatch::calculate_variance_partition(
        pbf_object,
        model_variables = c(batch_col, biological_factors),
        fill_the_missing = fill_missing
    )
}, error_msg = "calculate_variance_partition (raw) failed")
save_metric_helper(vp_raw, "variance_partition_raw", output_dirs$metrics_raw)

p_vp_raw <- safe_try({
    proBatch::plot_variance_partition(
        pbf_object,
        technical_factors = batch_col,
        biological_factors = biological_factors,
        variance_threshold = variance_threshold,
        summary_stat = "boxplot",
        fill_the_missing = fill_missing,
        base_size = 10
    )
}, error_msg = "plot_variance_partition (raw) failed")
save_plot_helper(p_vp_raw, "variance_partition_raw", output_dirs$plots_raw)

# --- Correlation analyses ---
log_msg("Calculating correlation diagnostics")

sample_corr_df_raw <- safe_try({
    proBatch::calculate_sample_corr_distr(
        pbf_object,
        sample_annotation,
        batch_col = batch_col,
        biospecimen_id_col = biospecimen_id_col,
        pbf_name = paste0(assay_prefix, "::log2_on_raw")
    )
}, error_msg = "calculate_sample_corr_distr (raw) failed")
save_metric_helper(sample_corr_df_raw, "sample_correlations_raw", output_dirs$metrics_raw)

p_sample_corr_raw <- safe_try({
    proBatch::plot_sample_corr_distribution(
        pbf_object,
        sample_annotation,
        batch_col = batch_col,
        biospecimen_id_col = biospecimen_id_col,
        pbf_name = paste0(assay_prefix, "::log2_on_raw"),
        plot_param = "batch_replicate",
        plot_title = paste0(
            "Sample correlation distribution (raw),\n",
            biospecimen_id_col, " and ", batch_col, " columns used for grouping"
        )
    )
}, error_msg = "plot_sample_corr_distribution (raw) failed")
save_plot_helper(p_sample_corr_raw, "sample_corr_distribution_raw", output_dirs$plots_raw)

if (!is.null(peptide_annotation) && "Gene" %in% names(peptide_annotation)) {
    peptide_corr_df_raw <- safe_try({
        proBatch::calculate_peptide_corr_distr(
            pbf_object,
            pbf_name = paste0(assay_prefix, "::log2_on_raw"),
            peptide_annotation,
            protein_col = "Gene"
        )
    }, error_msg = "calculate_peptide_corr_distr (raw) failed")
    save_metric_helper(peptide_corr_df_raw, "peptide_correlations_raw", output_dirs$metrics_raw)
    
    p_peptide_corr_raw <- safe_try({
        proBatch::plot_peptide_corr_distribution(
            pbf_object,
            pbf_name = paste0(assay_prefix, "::log2_on_raw"),
            peptide_annotation,
            protein_col = "Gene"
        )
    }, error_msg = "plot_peptide_corr_distribution (raw) failed")
    save_plot_helper(p_peptide_corr_raw, "peptide_corr_distribution_raw", output_dirs$plots_raw)
}

# Correlation heatmaps
p_corr_heatmap_raw <- safe_try({
    proBatch::plot_sample_corr_heatmap(
        pbf_object,
        sample_annotation = sample_annotation,
        samples_to_plot = samples_for_corr_heatmap,
        factors_to_plot = c(batch_col, condition_col),
        cluster_rows = TRUE,
        cluster_cols = TRUE
    )
}, error_msg = "plot_sample_corr_heatmap (raw) failed")
save_plot_helper(p_corr_heatmap_raw, "sample_corr_heatmap_raw", output_dirs$plots_raw)

# --- Feature-level diagnostics ---
log_msg("Generating feature-level diagnostics")

p_single_feature_raw <- safe_try({
    proBatch::plot_single_feature(
        feature_name = example_feature_id,
        df_long = pbf_object,
        pbf_name = paste0(assay_prefix, "::log2_on_raw"),
        sample_annotation = sample_annotation,
        batch_col = batch_col,
        color_by_batch = TRUE,
        base_size = 13,
        plot_title = paste0("Single feature plot (raw), feature: ", example_feature_id)
    ) 
}, error_msg = "plot_single_feature (raw) failed")
save_plot_helper(p_single_feature_raw, "single_feature_raw", output_dirs$plots_raw)

if (!is.null(peptide_annotation) && !is.null(example_protein_name)) {
    p_protein_raw <- safe_try({
        proBatch::plot_peptides_of_one_protein(
            protein_name = example_protein_name,
            peptide_annotation = peptide_annotation,
            protein_col = "Gene",
            df_long = pbf_object,
            pbf_name = paste0(assay_prefix, "::log2_on_raw"),
            sample_annotation = sample_annotation,
            batch_col = batch_col,
            color_by_batch = TRUE,
            base_size = 13,
            plot_title = paste0("Peptides of one protein (raw, log2), protein: ", example_protein_name)
        )
    }, error_msg = "plot_peptides_of_one_protein (raw) failed")
    save_plot_helper(p_protein_raw, "protein_peptides_raw", output_dirs$plots_raw)
}

if (!is.null(peptide_annotation) && !is.null(example_spike_in_name)) {
    p_spike_raw <- safe_try({
        proBatch::plot_spike_in(
            spike_ins = example_spike_in_name,
            peptide_annotation = peptide_annotation,
            protein_col = "Gene",
            df_long = pbf_object,
            pbf_name = paste0(assay_prefix, "::log2_on_raw"),
            sample_annotation = sample_annotation,
            batch_col = batch_col,
            color_by_batch = TRUE,
            base_size = 13
        )
    }, error_msg = "plot_spike_in (raw) failed")
    save_plot_helper(p_spike_raw, "spike_in_raw", output_dirs$plots_raw)
}

if (!is.null(peptide_annotation) && !is.null(example_irt_pattern)) {
    p_irt_raw <- safe_try({
        proBatch::plot_iRT(
            irt_pattern = example_irt_pattern,
            peptide_annotation = peptide_annotation,
            protein_col = "ProteinName",
            df_long = pbf_object,
            pbf_name = paste0(assay_prefix, "::log2_on_raw"),
            sample_annotation = sample_annotation,
            batch_col = batch_col,
            color_by_batch = TRUE,
            base_size = 13
        )
    }, error_msg = "plot_iRT (raw) failed")
    save_plot_helper(p_irt_raw, "iRT_raw", output_dirs$plots_raw)
}


# --- QC summaries ---
log_msg("Calculating QC summary metrics")

cv_df_raw <- safe_try({
    proBatch::calculate_feature_CV(
        pbf_object,
        pbf_name = paste0(assay_prefix, "::raw"),
        sample_annotation = sample_annotation,
        batch_col = batch_col,
        biospecimen_id_col = biospecimen_id_col,
        unlog = FALSE
    )
}, error_msg = "calculate_feature_CV (raw) failed")
save_metric_helper(cv_df_raw, "feature_CV_raw", output_dirs$metrics_raw)

p_cv_raw <- safe_try({
    proBatch::plot_CV_distr(
        pbf_object,
        pbf_name = paste0(assay_prefix, "::raw"),
        sample_annotation = sample_annotation,
        batch_col = batch_col,
        biospecimen_id_col = biospecimen_id_col,
        unlog = FALSE,
        plot_title = "Feature CV distribution (raw)"
    ) + ggtitle(
            paste0("Feature CV distribution (raw), ", batch_col, " used for grouping")
        )
}, error_msg = "plot_CV_distr (raw) failed")
save_plot_helper(p_cv_raw, "CV_distribution_raw", output_dirs$plots_raw)

class_metrics_raw <- safe_classification_metrics(
    pbf_object,
    sample_annotation,
    fill_the_missing = fill_missing,
    known_col = c(batch_col, condition_col),
)
save_metric_helper(class_metrics_raw, "classification_metrics_raw", output_dirs$metrics_raw)

intragroup_raw <- safe_try({
    proBatch::plot_intragroup_variation(
        pbf_object,
        group_col = c(condition_col, batch_col),
        fill_the_missing = fill_missing,
        metrics = c("correlation", "PCV")
    )
}, error_msg = "plot_intragroup_variation (raw) failed")
save_plot_helper(intragroup_raw, "intragroup_variation_raw", output_dirs$plots_raw)

outliers_raw <- safe_try({
    proBatch::detect_outlier_samples(
        pbf_object,
        pbf_name = paste0(assay_prefix, "::log2_on_raw"),
        sample_annotation = sample_annotation,
        batch_col = batch_col,
        n_pcs = 5
    )
}, error_msg = "detect_outlier_samples (raw) failed")
save_metric_helper(outliers_raw, "outlier_samples_raw", output_dirs$metrics_raw)

subbatches_raw <- safe_try({
    proBatch::subbatch_detection(
        pbf_object,
        pbf_name = paste0(assay_prefix, "::log2_on_raw"),
        sample_annotation = sample_annotation,
        batch_col = batch_col,
        n_pcs = 5,
        method = "kmeans",
        k_max = 4
    )
}, error_msg = "subbatch_detection (raw) failed")
if (!is.null(subbatches_raw) && "summary" %in% names(subbatches_raw)) {
    save_metric_helper(subbatches_raw$summary, "subbatch_detection_raw", output_dirs$metrics_raw)
}

log_msg("Baseline diagnostics complete")

# ==============================================================================
# BATCH EFFECT CORRECTION
# ==============================================================================


log_msg("=" %+% "=" %+% "=" %+% " BATCH EFFECT CORRECTION " %+% "=" %+% "=" %+% "=")
log_msg(sprintf("Applying correction method: %s", correction_method))

# Validate batch column has at least 2 levels
batch_levels <- unique(sample_annotation[[batch_col]][!is.na(sample_annotation[[batch_col]])])
if (length(batch_levels) < 2L) {
    log_msg("WARNING: Batch column has < 2 levels. Skipping correction.", level = "WARN")
    pbf_object <- pbf_object
} else {
    # Validate and filter covariates
    valid_covariates <- correction_covariates[correction_covariates %in% names(sample_annotation)]
    valid_covariates <- valid_covariates[vapply(valid_covariates, function(col) {
        length(unique(sample_annotation[[col]][!is.na(sample_annotation[[col]])])) >= 2
    }, logical(1))]
    
    # Apply correction
    pbf_object <- safe_try({
        proBatch::pb_transform(
            pbf_object,
            from = paste0(assay_prefix, "::log2_on_raw"),
            steps = correction_method,
            params_list = list(list(
                batch_col = batch_col,
                covariates_cols = valid_covariates
            ))
        )
    }, error_msg = sprintf("Batch correction (%s) failed", correction_method),
       return_on_error = pbf_object)
}

#################################################################################
######################################################################################
#######################################################################################
# Extract corrected data
data_matrix_corrected <- proBatch::pb_assay_matrix(pbf_object)
df_long_corrected <- proBatch::pb_as_long(
    pbf_object,
    feature_id_col = feature_id_col,
    sample_id_col = sample_id_col,
    measure_col = measure_col
)

# Save corrected objects
if (save_objects) {
    saveRDS(pbf_object, file.path(output_dirs$objects, "pbf_corrected.rds"))
    saveRDS(data_matrix_corrected, file.path(output_dirs$objects, "data_matrix_corrected.rds"))
    log_msg("Saved corrected ProBatchFeatures object and matrix")
}

#################################################################################
######################################################################################
#######################################################################################


# ==============================================================================
# POST-CORRECTION DIAGNOSTICS
# ==============================================================================

log_msg("=" %+% "=" %+% "=" %+% " POST-CORRECTION DIAGNOSTICS " %+% "=" %+% "=" %+% "=")

# --- Missingness diagnostics ---
log_msg("Generating post-correction missingness diagnostics")

p_na_density_corr <- safe_try(
    proBatch::plot_NA_density(
        pbf_object,
        from = paste0(assay_prefix, "::", correction_method, "_on_log2_on_raw"),
    ) + ggtitle(
            paste0("Corrected: NA density\n", 
            assay_prefix, "::", correction_method, "_on_log2_on_raw")
        ),
    error_msg = "plot_NA_density (corrected) failed"
)
save_plot_helper(p_na_density_corr, "NA_density_corr", output_dirs$plots_corr)


# --- Sample distributions ---
log_msg("Generating post-correction sample distribution diagnostics")

p_sample_mean_corr <- safe_try({
    proBatch::plot_sample_mean(
        pbf_object,
        pbf_name = paste0(assay_prefix, "::", correction_method, "_on_log2_on_raw"),
        order_col = order_col,
        batch_col = batch_col,
        color_by_batch = TRUE,
        color_scheme = color_list,
        base_size = 15
    ) + ggtitle(
        paste0("Sample means after correction\n", 
                assay_prefix, "::", correction_method, "_on_log2_on_raw")
    )
}, error_msg = "plot_sample_mean (corrected) failed")
save_plot_helper(p_sample_mean_corr, "sample_mean_corr", output_dirs$plots_corr)


# will plot all assays available in the pbf object
p_boxplot_corr <- safe_try({
    proBatch::plot_boxplot(
        pbf_object,
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        color_by_batch = TRUE,
        color_scheme = color_list,
        base_size = 10
    ) + ggtitle(paste0("Boxplot of all assays, colored by ", batch_col))
}, error_msg = "plot_boxplot (corrected) failed")
save_plot_helper(p_boxplot_corr, "boxplot_corr", output_dirs$plots_corr)

# --- Heatmaps and clustering ---
log_msg("Generating post-correction heatmaps and clustering diagnostics")

# also plot all assays in pbf_object
p_heatmap_diag_corr <- safe_try({
    proBatch::plot_heatmap_diagnostic(
        pbf_object,
        factors_to_plot = factors_to_plot,
        color_list = color_list,
        show_rownames = FALSE,
        show_colnames = FALSE,
        fill_the_missing = fill_missing
    )
}, error_msg = "plot_heatmap_diagnostic (corrected) failed")
save_plot_helper(p_heatmap_diag_corr, "heatmap_diagnostic_corr", output_dirs$plots_corr)

p_hierarchical_corr <- safe_try({
    proBatch::plot_hierarchical_clustering(
        pbf_object,
        factors_to_plot = factors_to_plot,
        color_list = color_list,
        label_samples = FALSE,
        fill_the_missing = fill_missing
    )
}, error_msg = "plot_hierarchical_clustering (corrected) failed")
save_plot_helper(p_hierarchical_corr, "hierarchical_clustering_corr", output_dirs$plots_corr)

# --- Embeddings ---
log_msg("Generating post-correction embeddings")

pca_corr <- safe_try({
    proBatch::plot_PCA(
        pbf_object,
        color_by = batch_col,
        shape_by = condition_col,
        fill_the_missing = fill_missing
    )
}, error_msg = "plot_PCA (corrected) failed")
save_plot_helper(pca_corr, "PCA_corr", output_dirs$plots_corr)

tsne_corr <- safe_try({
    proBatch::plot_TSNE(
        pbf_object,
        color_by = batch_col,
        shape_by = condition_col,
        perplexity = perplexity_tsne,
        max_iter = 300,
        initial_dims = n_pcs_embeddings,
        fill_the_missing = fill_missing,
        point_size = 4
    )
}, error_msg = "plot_TSNE (corrected) failed")
save_plot_helper(tsne_corr, "tSNE_corr", output_dirs$plots_corr)

umap_corr <- safe_try({
    proBatch::plot_UMAP(
        pbf_object,
        color_by = batch_col,
        shape_by = condition_col,
        n_neighbors = n_neighbors_umap,
        min_dist = 0.3,
        random_state = set_seed,
        fill_the_missing = fill_missing,
        point_size = 4
    )
}, error_msg = "plot_UMAP (corrected) failed")
save_plot_helper(umap_corr, "UMAP_corr", output_dirs$plots_corr)

# --- PVCA ---
log_msg("Calculating post-correction PVCA")

pvca_corr <- safe_try({
    proBatch::calculate_PVCA(
        pbf_object,
        factors_for_PVCA = c(technical_factors, biological_factors),
        fill_the_missing = fill_missing
    )
}, error_msg = "calculate_PVCA (corrected) failed")
save_metric_helper(pvca_corr, "pvca_corr", output_dirs$metrics_corr)

p_pvca_corr <- safe_try({
    proBatch::plot_PVCA(
        pbf_object,
        technical_factors = technical_factors,
        biological_factors = biological_factors,
        fill_the_missing = fill_missing,
        variance_threshold = (ifelse(
            length(c(technical_factors, biological_factors)) > 5, 0.01, 0)),
        base_size = 10
    )
}, error_msg = "plot_PVCA (corrected) failed")
save_plot_helper(p_pvca_corr, "PVCA_corr", output_dirs$plots_corr)

pvca_df_corr <- safe_try({
    proBatch::prepare_PVCA_df(
        pbf_object,
        technical_factors = technical_factors,
        biological_factors = biological_factors,
        fill_the_missing = fill_missing,
        variance_threshold = variance_threshold,
        path_to_save_results = output_dirs$pvca_corr
    )
}, error_msg = "prepare_PVCA_df (corrected) failed")
save_metric_helper(pvca_df_corr, "pvca_prepared_corr", output_dirs$metrics_corr)

# --- Variance partition ---
log_msg("Calculating post-correction variance partition")

vp_corr <- safe_try({
    proBatch::calculate_variance_partition(
        pbf_object,
        model_variables = c(batch_col, biological_factors),
        fill_the_missing = fill_missing
    )
}, error_msg = "calculate_variance_partition (corrected) failed")
save_metric_helper(vp_corr, "variance_partition_corr", output_dirs$metrics_corr)

p_vp_corr <- safe_try({
    proBatch::plot_variance_partition(
        pbf_object,
        technical_factors = batch_col,
        biological_factors = biological_factors,
        variance_threshold = variance_threshold,
        summary_stat = "boxplot",
        fill_the_missing = fill_missing,
        base_size = 10
    )
}, error_msg = "plot_variance_partition (corrected) failed")
save_plot_helper(p_vp_corr, "variance_partition_corr", output_dirs$plots_corr)

# --- Correlation analyses ---
log_msg("Calculating post-correction correlation diagnostics")

sample_corr_df_corr <- safe_try({
    proBatch::calculate_sample_corr_distr(
        pbf_object,
        sample_annotation,
        batch_col = batch_col,
        biospecimen_id_col = biospecimen_id_col,
        pbf_name = paste0(assay_prefix, "::", correction_method, "_on_log2_on_raw")
    )
}, error_msg = "calculate_sample_corr_distr (corrected) failed")
save_metric_helper(sample_corr_df_corr, "sample_correlations_corr", output_dirs$metrics_corr)

p_sample_corr_corr <- safe_try({
    proBatch::plot_sample_corr_distribution(
        pbf_object,
        sample_annotation,
        batch_col = batch_col,
        biospecimen_id_col = biospecimen_id_col,
        pbf_name = paste0(assay_prefix, "::", correction_method, "_on_log2_on_raw"),
        plot_param = "batch_replicate",
        plot_title = paste0(
            "Sample correlation distribution after correction,\n", 
            biospecimen_id_col, " and ", batch_col, " columns used for grouping"
        )
    )
}, error_msg = "plot_sample_corr_distribution (corrected) failed")
save_plot_helper(p_sample_corr_corr, "sample_corr_distribution_corr", output_dirs$plots_corr)

if (!is.null(peptide_annotation) && "Gene" %in% names(peptide_annotation)) {
    peptide_corr_df_corr <- safe_try({
        proBatch::calculate_peptide_corr_distr(
            pbf_object,
            pbf_name = paste0(assay_prefix, "::", correction_method, "_on_log2_on_raw"),
            peptide_annotation,
            protein_col = "Gene"
        )
    }, error_msg = "calculate_peptide_corr_distr (corrected) failed")
    save_metric_helper(peptide_corr_df_corr, "peptide_correlations_corr", output_dirs$metrics_corr)
    
    p_peptide_corr_corr <- safe_try({
        proBatch::plot_peptide_corr_distribution(
            pbf_object,
            pbf_name = paste0(assay_prefix, "::", correction_method, "_on_log2_on_raw"),
            peptide_annotation,
            protein_col = "Gene"
        )
    }, error_msg = "plot_peptide_corr_distribution (corrected) failed")
    save_plot_helper(p_peptide_corr_corr, "peptide_corr_distribution_corr", output_dirs$plots_corr)
}

# Correlation heatmaps
p_corr_heatmap_corr <- safe_try({
    proBatch::plot_sample_corr_heatmap(
        pbf_object,
        sample_annotation = sample_annotation,
        samples_to_plot = samples_for_corr_heatmap,
        factors_to_plot = c(batch_col, condition_col),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        color_list = color_list,
        show_row_dend = FALSE, show_column_dend = FALSE,
        x_axis_label_size = 4, y_axis_label_size = 4,
    )  + ggtitle("Sample correlation heatmap after correction")
}, error_msg = "plot_sample_corr_heatmap (corrected) failed")
save_plot_helper(p_corr_heatmap_corr, "sample_corr_heatmap_corr", output_dirs$plots_corr)

# --- Feature-level diagnostics ---
log_msg("Generating post-correction feature-level diagnostics")
if (!is.null(peptide_annotation) && !is.null(example_spike_in_name)) {
    p_spike_corr <- safe_try({
        proBatch::plot_spike_in(
            spike_ins = example_spike_in_name,
            peptide_annotation = peptide_annotation,
            protein_col = "Gene",
            df_long = pbf_object,
            pbf_name = paste0(assay_prefix, "::", correction_method, "_on_log2_on_raw"),
            sample_annotation = sample_annotation,
            batch_col = batch_col,
            color_by_batch = TRUE,
            base_size = 13
        )
    }, error_msg = "plot_spike_in (corrected) failed")
    save_plot_helper(p_spike_corr, "spike_in_corr", output_dirs$plots_corr)
}

if (!is.null(peptide_annotation) && !is.null(example_irt_pattern)) {
    p_irt_corr <- safe_try({
        proBatch::plot_iRT(
            irt_pattern = example_irt_pattern,
            peptide_annotation = peptide_annotation,
            protein_col = "ProteinName",
            df_long = pbf_object,
            pbf_name = paste0(assay_prefix, "::", correction_method, "_on_log2_on_raw"),
            sample_annotation = sample_annotation,
            batch_col = batch_col,
            color_by_batch = TRUE,
            base_size = 13
        )
    }, error_msg = "plot_iRT (corrected) failed")
    save_plot_helper(p_irt_corr, "iRT_corr", output_dirs$plots_corr)
}


# --- QC summaries ---
log_msg("Calculating post-correction QC summary metrics")

cv_df_corr <- safe_try({
    proBatch::calculate_feature_CV(
        pbf_object,
        pbf_name = paste0(assay_prefix, "::", correction_method, "_on_log2_on_raw"),
        sample_annotation = sample_annotation,
        batch_col = batch_col,
        biospecimen_id_col = biospecimen_id_col,
        unlog = FALSE
    )
}, error_msg = "calculate_feature_CV (corrected) failed")
save_metric_helper(cv_df_corr, "feature_CV_corr", output_dirs$metrics_corr)

p_cv_corr <- safe_try({
    proBatch::plot_CV_distr(
        pbf_object,
        pbf_name = paste0(assay_prefix, "::", correction_method, "_on_log2_on_raw"),
        sample_annotation = sample_annotation,
        batch_col = batch_col,
        biospecimen_id_col = biospecimen_id_col,
        unlog = FALSE
    ) + ggtitle(paste0("Feature CV distribution after correction\n", 
                    assay_prefix, "::", correction_method, "_on_log2_on_raw"))
}, error_msg = "plot_CV_distr (corrected) failed")
save_plot_helper(p_cv_corr, "CV_distribution_corr", output_dirs$plots_corr)

class_metrics_corr <- safe_classification_metrics(
    pbf_object,
    sample_annotation,
    fill_the_missing = fill_missing,
    known_col = c(batch_col, condition_col)
)
save_metric_helper(class_metrics_corr, "classification_metrics_corr", output_dirs$metrics_corr)

intragroup_corr <- safe_try({
    proBatch::plot_intragroup_variation(
        pbf_object,
        group_col = condition_col,
        fill_the_missing = fill_missing,
        metrics = c("correlation", "PCV")
    )
}, error_msg = "plot_intragroup_variation (corrected) failed")
save_plot_helper(intragroup_corr, "intragroup_variation_corr", output_dirs$plots_corr)

subbatches_corr <- safe_try({
    proBatch::subbatch_detection(
        pbf_object,
        pbf_name = paste0(assay_prefix, "::", correction_method, "_on_log2_on_raw"),
        sample_annotation = sample_annotation,
        batch_col = batch_col,
        n_pcs = 5,
        method = "kmeans",
        k_max = 4
    )
}, error_msg = "subbatch_detection (corrected) failed")
if (!is.null(subbatches_corr) && "summary" %in% names(subbatches_corr)) {
    save_metric_helper(subbatches_corr$summary, "subbatch_detection_corr", output_dirs$metrics_corr)
}

log_msg("Post-correction diagnostics complete")

# ==============================================================================
# COMPARATIVE ANALYSIS (RAW vs CORRECTED)
# ==============================================================================

log_msg("=" %+% "=" %+% "=" %+% " COMPARATIVE ANALYSIS " %+% "=" %+% "=" %+% "=")

# --- Summary metrics table ---
log_msg("Generating summary metrics comparison table")

first_num <- function(x) {
    if (is.null(x) || !is.numeric(x) || length(x) == 0) return(NA_real_)
    x[!is.na(x)][1]
}

pvca_share <- function(pvca_df, patterns) {
    if (is.null(pvca_df) || !is.data.frame(pvca_df)) return(NA_real_)
    matched <- pvca_df[grepl(paste(patterns, collapse = "|"), 
                            pvca_df$label, ignore.case = TRUE), ]
    if (nrow(matched) == 0) return(NA_real_)
    sum(matched$Weighted.average.proportion.variance, na.rm = TRUE)
}

median_corr_by_group <- function(corr_df, group_label) {
    if (is.null(corr_df) || !is.data.frame(corr_df)) return(NA_real_)
    subset_df <- corr_df[corr_df$batch_replicate == group_label, ]
    if (nrow(subset_df) == 0) return(NA_real_)
    median(subset_df$correlation, na.rm = TRUE)
}

count_outliers <- function(x) {
    if (is.null(x)) return(NA_integer_)
    if (is.data.frame(x) && "outlier" %in% names(x)) {
        return(sum(x$outlier, na.rm = TRUE))
    }
    if (is.list(x) && "outliers" %in% names(x)) {
        return(length(x$outliers))
    }
    return(NA_integer_)
}

summary_metrics <- tibble::tibble(
    Metric = c(
        "PVCA: Technical %",
        "PVCA: Biological %",
        "Median within-replicate corr",
        "Median within-batch corr",
        "Silhouette (batch)",
        "ASW (batch)",
        "Number of outliers",
        "Median feature CV (within-replicate)"
    ),
    Raw = c(
        pvca_share(pvca_df_raw, technical_factors) * 100,
        pvca_share(pvca_df_raw, biological_factors) * 100,
        median_corr_by_group(sample_corr_df_raw, "within_replicate"),
        median_corr_by_group(sample_corr_df_raw, "within_batch"),
        first_num(class_metrics_raw$silhouette_score),
        first_num(class_metrics_raw$avg_silhouette_width),
        count_outliers(outliers_raw),
        if (!is.null(cv_df_raw)) median(cv_df_raw$CV_replicate, na.rm = TRUE) else NA_real_
    ),
    Corrected = c(
        pvca_share(pvca_df_corr, technical_factors) * 100,
        pvca_share(pvca_df_corr, biological_factors) * 100,
        median_corr_by_group(sample_corr_df_corr, "within_replicate"),
        median_corr_by_group(sample_corr_df_corr, "within_batch"),
        first_num(class_metrics_corr$silhouette_score),
        first_num(class_metrics_corr$avg_silhouette_width),
        count_outliers(outliers_corr),
        if (!is.null(cv_df_corr)) median(cv_df_corr$CV_replicate, na.rm = TRUE) else NA_real_
    )
) %>%
    mutate(Change = Corrected - Raw)

save_metric_helper(summary_metrics, "summary_metrics_comparison", output_dirs$metrics_comp)

# --- Classification metrics comparison ---
if (!is.null(class_metrics_raw) && !is.null(class_metrics_corr)) {
    class_metrics_compare <- bind_rows(
        transform(class_metrics_raw, Step = "raw"),
        transform(class_metrics_corr, Step = "corrected")
    )
    save_metric_helper(class_metrics_compare, "classification_metrics_comparison", 
                      output_dirs$metrics_comp)
}

# --- Split violin comparison ---
log_msg("Generating split violin comparison")

if (!is.null(sample_corr_df_raw) && !is.null(sample_corr_df_corr)) {
    corr_compare <- bind_rows(
        transform(sample_corr_df_raw, Step = "raw"),
        transform(sample_corr_df_corr, Step = "corrected")
    )
    
    p_split_violin <- safe_try({
        proBatch::plot_split_violin_with_boxplot(
            corr_compare,
            y_col = "correlation",
            col_for_color = "Step",
            col_for_box = "batch_replicate",
            plot_title = sprintf("Sample correlation: raw vs %s", correction_method)
        )
    }, error_msg = "Split violin plot failed")
    save_plot_helper(p_split_violin, "split_violin_correlation", output_dirs$plots_comp)
    
    save_metric_helper(corr_compare, "correlation_comparison", output_dirs$metrics_comp)
}

log_msg("Comparative analysis complete")

# ==============================================================================
# GENERATE SUMMARY REPORT
# ==============================================================================

log_msg("Generating summary report")

report_file <- file.path(output_base_dir, "SUMMARY_REPORT.txt")

sink(report_file)
cat("=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "\n")
cat("BATCH EFFECT CORRECTION - SUMMARY REPORT\n")
cat("=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "\n\n")

cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Output directory:", output_base_dir, "\n\n")

cat("--- Configuration ---\n")
cat("Correction method:", correction_method, "\n")
cat("Batch column:", batch_col, "\n")
cat("Condition column:", condition_col, "\n")
cat("Technical factors:", paste(technical_factors, collapse = ", "), "\n")
cat("Biological factors:", paste(biological_factors, collapse = ", "), "\n")
cat("Biospecimen ID:", biospecimen_id_col, "\n\n")

cat("--- Data dimensions ---\n")
cat("Features:", nrow(data_matrix_raw_log), "\n")
cat("Samples:", ncol(data_matrix_raw_log), "\n\n")

cat("--- Summary metrics ---\n")
print(summary_metrics)
cat("\n")

if (!is.null(design_check_raw)) {
    cat("--- Design validation ---\n")
    print(design_check_raw)
    cat("\n")
}

if (!is.null(nested_raw)) {
    cat("--- Nested batch structure ---\n")
    print(nested_raw)
    cat("\n")
}

cat("=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "\n")
cat("END OF REPORT\n")
cat("=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "=" %+% "\n")
sink()

log_msg(sprintf("Summary report saved: %s", report_file))

# ==============================================================================
# WORKFLOW COMPLETE
# ==============================================================================

log_msg("=" %+% "=" %+% "=" %+% " WORKFLOW COMPLETE " %+% "=" %+% "=" %+% "=")
log_msg(sprintf("All results saved to: %s", output_base_dir))
log_msg("Check SUMMARY_REPORT.txt for key findings")

cat("\n")
cat("Directory structure:\n")
cat("  01_raw_diagnostics/\n")
cat("    plots/          - Diagnostic plots for raw data\n")
cat("    metrics/        - Metric tables (CSV) for raw data\n")
cat("    pvca_results/   - PVCA intermediate results\n")
cat("  02_corrected_diagnostics/\n")
cat("    plots/          - Diagnostic plots for corrected data\n")
cat("    metrics/        - Metric tables (CSV) for corrected data\n")
cat("    pvca_results/   - PVCA intermediate results\n")
cat("  03_raw_vs_corrected/\n")
cat("    plots/          - Comparative plots\n")
cat("    metrics/        - Comparative metric tables\n")
cat("  04_objects/\n")
cat("    pbf_corrected.rds           - Corrected ProBatchFeatures object\n")
cat("    data_matrix_corrected.rds   - Corrected data matrix\n")
cat("  SUMMARY_REPORT.txt            - Summary of key metrics\n")
cat("\n")

# Return key objects (for interactive use)
invisible(list(
    pbf_raw = pbf_log,
    pbf_corrected = pbf_object,
    summary_metrics = summary_metrics,
    output_dir = output_base_dir
))
