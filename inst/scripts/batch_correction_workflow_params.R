#!/usr/bin/env Rscript
# User-editable parameters for batch_correction_workflow.R.
# Copy this file, edit values in `workflow_params`, then run it.

workflow_params <- list(
    # --- Input data ---
    pbf_object = NULL,
    pbf_file_path = NULL,
    assay_prefix = "PGs",

    # --- Annotation data ---
    sample_annotation_path = NULL,
    peptide_annotation_path = NULL,

    # --- Column identifiers ---
    sample_id_col = "FullRunName",
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",

    # --- Technical factors (batch variables) ---
    batch_col = "MS_batch",
    technical_factors = c("MS_batch", "digestion_batch"),

    # --- Biological factors ---
    condition_col = "Diet",
    biological_factors = c("Diet", "Sex", "Strain"),
    biospecimen_id_col = "EarTag",

    # --- Covariates for correction ---
    correction_covariates = c("Diet", "Sex"),

    # --- Other metadata columns ---
    order_col = "order",
    datetime_col = "DateTime",

    # --- Factors for visualization ---
    factors_to_plot = c("MS_batch", "Diet", "Sex", "Strain"),
    factor_columns_colors = c("MS_batch", "Diet", "Sex", "Strain", "digestion_batch"),
    numeric_columns_colors = c("DateTime", "order"),

    # --- Feature-level selectors ---
    example_protein_name = NULL,
    example_spike_in_name = NULL,
    example_irt_pattern = "iRT",
    example_feature_id = NULL,

    # --- Batch correction methods ---
    correction_methods = c("limmaRBE"),
    correction_tasks_yaml = "inst/scripts/pb_tasks.yaml",
    correction_task_labels = NULL,
    enable_rowname_repair_retry = FALSE,

    # --- Output configuration ---
    output_base_dir = file.path(tempdir(), "batch_correction_results"),
    results_prefix = "BEC",
    save_plots = TRUE,
    save_metrics = TRUE,
    save_objects = TRUE,
    # Optional root for NormAE stdout/stderr logs.
    # Keep NULL for current behavior (logs in current working directory).
    # To colocate with workflow outputs, set after this list is created:
    # workflow_params$normae_log_base_dir <- file.path(workflow_params$output_base_dir, "02_corrected_diagnostics", "logs", "normae")
    normae_log_base_dir = NULL,
    # Optional root for (s)PLSDA-batch R output/messages.
    # Keep NULL for current behavior (messages in workflow log).
    # To colocate with workflow outputs, set after this list is created:
    # workflow_params$plsdabatch_log_base_dir <- file.path(workflow_params$output_base_dir, "02_corrected_diagnostics", "logs", "plsdabatch")
    plsdabatch_log_base_dir = NULL,

    # --- Plot settings ---
    plot_format = "png",
    plot_width = 10,
    plot_height = 6,
    plot_dpi = 300,

    # --- Computational settings ---
    set_seed = 42,
    n_pcs_embeddings = 10,
    perplexity_tsne = 10,
    n_neighbors_umap = 10,
    fill_missing = -1,
    variance_threshold = 0.05,

    # --- Sample subsetting ---
    samples_for_corr_heatmap = NULL,
    n_samples_corr_heatmap = 20,

    # --- Caching ---
    cache_expensive = TRUE
)

file_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_args) > 0L) {
    dirname(normalizePath(sub("^--file=", "", file_args[1]), mustWork = FALSE))
} else {
    getwd()
}

workflow_candidates <- c(
    file.path(script_dir, "batch_correction_workflow.R"),
    file.path(script_dir, "inst", "scripts", "batch_correction_workflow.R"),
    file.path(getwd(), "inst", "scripts", "batch_correction_workflow.R")
)
workflow_path <- workflow_candidates[file.exists(workflow_candidates)][1]

if (is.na(workflow_path) || !nzchar(workflow_path)) {
    stop("Could not locate 'batch_correction_workflow.R'.")
}

source(workflow_path, local = globalenv())
