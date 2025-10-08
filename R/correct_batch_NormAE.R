#######################################
##  Functions to run (train) NormAE  ##
#######################################

#' @title Batch correction via NormAE (Python CLI)
#' @description
#' Wrapper around the Python package \code{normae} (Normalization Autoencoder)
#' to remove batch effects using the official CLI. Works on wide matrices and
#' long data.frames by converting to/from matrices.
#'
#' @inheritParams correct_with_ComBat
#' @param format One of \code{"wide"} or \code{"long"}.
#' @param inj_order_col Column in \code{sample_annotation} with injection order.
#'   If \code{NULL}, the order of \code{sample_annotation[[sample_id_col]]} is used.
#' @param qc_col_name Optional column in \code{sample_annotation} marking QC samples.
#'   Logical, factor (level "QC"), or character ("QC") values are treated as QC.
#' @param normae_args Named list of CLI parameters forwarded to \code{normae}.
#'   Entries are translated to \code{--name value} pairs. Default: empty list.
#'   Some defaults are set internally if not provided:
#'     --mz_row set to "".
#'     --rt_row set to "".
#'     --early_stop set to False.
#' @param python_env Optional path to a Python binary.
#' @param conda_env Optional conda environment name. Ignored when \code{python_env}
#'   is supplied. If neither is provided, a \code{"normae"} conda env is created.
#' @return Same type as input: numeric \code{matrix} for \code{format="wide"},
#'   or a long \code{data.frame} for \code{format="long"}.
#' @details Requires Python (>=3.10) with the \pkg{normae} package installed.
#'   The first call without \code{python_env}/\code{conda_env} creates a conda
#'   environment named \code{"normae"} and installs the package from GitHub.
#' @references Rong, Zhiwei, et al. "NormAE: deep adversarial learning model to remove batch effects in
#' liquid chromatography mass spectrometry-based metabolomics data." Analytical chemistry 92.7 (2020): 5082-5090.
#' GitHub: https://github.com/luyiyun/NormAE/tree/release
#' @export
correct_with_NormAE <- function(
    x,
    sample_annotation,
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    format = c("wide", "long"),
    inj_order_col = NULL,
    qc_col_name = NULL,
    keep_all = "default",
    python_env = NULL,
    conda_env = NULL,
    normae_args = list()) {
    format <- match.arg(format)

    if (!is.null(normae_args) && !is.list(normae_args)) {
        stop("normae_args must be a named list.")
    }
    if (is.null(normae_args)) {
        normae_args <- list()
    }

    if (identical(format, "wide")) {
        if (!is.matrix(x)) {
            stop("format='wide' requires a numeric matrix.")
        }
        return(.normae_matrix_step(
            data_matrix = x,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            batch_col = batch_col,
            inj_order_col = inj_order_col,
            qc_col_name = qc_col_name,
            python_env = python_env,
            conda_env = conda_env,
            normae_args = normae_args
        ))
    }

    if (!is.data.frame(x)) {
        stop("format='long' requires a data.frame.")
    }
    df_long <- x
    original_cols <- names(df_long)

    df_long <- check_sample_consistency(
        sample_annotation, sample_id_col, df_long,
        batch_col,
        order_col = NULL, facet_col = NULL, merge = FALSE
    )

    data_matrix <- long_to_matrix(
        df_long,
        feature_id_col = feature_id_col,
        measure_col    = measure_col,
        sample_id_col  = sample_id_col
    )

    corrected_matrix <- .normae_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        inj_order_col = inj_order_col,
        qc_col_name = qc_col_name,
        python_env = python_env,
        conda_env = conda_env,
        normae_args = normae_args
    )

    .post_correction_to_long(
        corrected_matrix, df_long,
        feature_id_col, measure_col, sample_id_col,
        original_cols, keep_all
    )
}

.normae_matrix_step <- function(
    data_matrix, sample_annotation,
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    inj_order_col = NULL,
    qc_col_name = NULL,
    python_env = NULL,
    conda_env = NULL,
    normae_args = list()) {
    if (is.null(sample_annotation)) stop("sample_annotation must be provided for NormAE correction.")
    if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
    }
    if (anyNA(data_matrix)) stop("NormAE requires no NAs; impute/filter before calling.")
    if (!is.numeric(data_matrix)) {
        stop("Input must be coercible to a numeric matrix for NormAE correction.")
    }
    storage.mode(data_matrix) <- "double"
    if (any(data_matrix < 0, na.rm = TRUE)) {
        message("NormAE expects non-negative intensities; negative values were detected.")
    }

    sample_ids <- colnames(data_matrix)
    if (is.null(sample_ids) || length(sample_ids) == 0 ||
        any(!sample_ids %in% as.character(sample_annotation[[sample_id_col]]))) {
        stop("data_matrix must have column names (sample IDs) present in sample_annotation.")
    }

    .run_normae_core(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        inj_col_name = inj_col_name,
        qc_col_name = qc_col_name,
        qc_indicator_value = qc_indicator_value,
        python_env = python_env,
        conda_env = conda_env,
        normae_args = normae_args
    )
}

.run_normae_core <- function(
    data_matrix, # numeric matrix (features × samples)
    sample_annotation,
    sample_id_col,
    batch_col,
    inj_col_name = NULL,
    qc_col_name = NULL,
    qc_indicator_value = NULL,
    python_env = NULL,
    conda_env = NULL,
    normae_args = list()) {
    stopifnot(is.matrix(data_matrix))
    if (is.null(sample_annotation)) stop("sample_annotation is required.")

    # save original rownames and colnames, set temporary if missing
    original_rownames <- rownames(data_matrix)
    sample_ids <- colnames(data_matrix)

    # Align SA to matrix columns (samples)
    sample_annotation <- .align_sample_annotation(
        sample_annotation,
        sample_ids    = colnames(data_matrix),
        sample_id_col = sample_id_col
    )
    # Ensure SA contains required columns
    if (!(batch_col %in% names(sample_annotation))) {
        stop("Batch column '", batch_col, "' is not present in sample_annotation.")
    }
    if (!is.null(inj_col_name)) {
        if (!(inj_col_name %in% names(sample_annotation))) {
            stop("Injection order column '", inj_col_name, "' not found in sample_annotation.")
        }
    }

    sample_annotation <- .preprocess_qc_column(
        sample_annotation,
        qc_col_name = qc_col_name,
        sample_id_col = sample_id_col
    )

    # Starting with python env for NormAE
    python_bin <- .normae_prepare_python(python_env, conda_env)

    # Transpose to samples × features
    meta_df <- as.data.frame(t(data_matrix), check.names = FALSE)
    rownames(meta_df) <- sample_ids

    # Create temporary files for meta and sample info
    meta_path <- tempfile("normae_meta_", fileext = ".csv")
    sample_path <- tempfile("normae_samples_", fileext = ".csv")
    output_dir <- tempfile("normae_output_")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(c(meta_path, sample_path), force = TRUE), add = TRUE)
    on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

    write.csv(meta_df, meta_path, quote = TRUE)
    write.csv(sa_subset, sample_path, quote = TRUE)

    cli_args <- c(
        "--meta_csv", normalizePath(meta_path, winslash = "/", mustWork = TRUE),
        "--sample_csv", normalizePath(sample_path, winslash = "/", mustWork = TRUE),
        "--output_dir", normalizePath(output_dir, winslash = "/", mustWork = TRUE),
        "--batch_indicator_col", batch_col,
        "--order_indicator_col", inj_col_name
    )

    if (!is.null(qc_col_name)) {
        qc_value <- if (is.null(qc_indicator_value)) "QC" else as.character(qc_indicator_value)
        cli_args <- c(
            cli_args,
            "--qc_indicator_col", qc_col_name,
            "--qc_indicator_value", qc_value
        )
    } else {
        cli_args <- c(cli_args, "--qc_indicator_col", "")
    }

    extra_args <- .normae_format_cli_args(normae_args)
    if (length(extra_args)) {
        cli_args <- c(cli_args, extra_args)
    }

    cli_output <- system2(
        command = python_bin,
        args = c("-m", "normae.cli", cli_args),
        stdout = TRUE,
        stderr = TRUE
    )
    status <- attr(cli_output, "status")
    if (is.null(status)) {
        status <- 0L
    }
    if (!identical(status, 0L)) {
        stop(
            paste(
                c(
                    "NormAE CLI execution failed.",
                    "Command output:",
                    cli_output
                ),
                collapse = "\n"
            )
        )
    }

    result_path <- file.path(output_dir, "X_clean.csv")
    if (!file.exists(result_path)) {
        stop("NormAE CLI did not produce 'X_clean.csv' in the output directory.")
    }

    result_df <- read.csv(result_path, check.names = FALSE, row.names = 1)
    result_mat <- as.matrix(result_df)
    storage.mode(result_mat) <- "double"

    if (identical(rownames(result_mat), original_rownames) && identical(colnames(result_mat), sample_ids)) {
        aligned <- result_mat
    } else if (identical(rownames(result_mat), sample_ids) && identical(colnames(result_mat), original_rownames)) {
        aligned <- t(result_mat)
    } else {
        missing_rows <- setdiff(original_rownames, rownames(result_mat))
        missing_cols <- setdiff(sample_ids, colnames(result_mat))
        stop(
            "NormAE output does not match expected dimensions. ",
            "Missing features: ", paste(missing_rows, collapse = ", "),
            "; missing samples: ", paste(missing_cols, collapse = ", ")
        )
    }
    out <- aligned[original_rownames, sample_ids, drop = FALSE]
    if (!identical(dim(out), dim(meta_df))) {
        stop("NormAE returned unexpected dimensions.")
    }
    our <- as.matrix(out)
    dimnames(out) <- dimnames(data_matrix)
    storage.mode(out) <- "double"
    out
}

.normae_prepare_python <- function(python_env = NULL, conda_env = NULL) {
    .pb_requireNamespace("reticulate")

    active_env <- NULL
    if (!is.null(python_env)) {
        reticulate::use_python(python_env, required = TRUE)
    } else if (!is.null(conda_env)) {
        reticulate::use_condaenv(conda_env, required = TRUE)
        active_env <- conda_env
    } else {
        active_env <- "normae"
        if (!reticulate::condaenv_exists(active_env)) {
            message("Creating conda environment '", active_env, "' (Python 3.10) for NormAE...")
            reticulate::conda_create(active_env, python_version = "3.10")
        }
        reticulate::use_condaenv(active_env, required = TRUE)
        if (!reticulate::py_module_available("normae")) {
            message("Installing NormAE into conda environment '", active_env, "'...")
            reticulate::conda_install(active_env, "git+https://github.com/luyiyun/NormAE.git", pip = TRUE)
        }
    }

    config <- reticulate::py_config(required = TRUE)

    if (!reticulate::py_module_available("normae")) {
        if (is.null(active_env)) {
            stop("Python module 'normae' is not available in the selected Python environment.")
        }
        stop(
            "Python module 'normae' is not available in the environment '",
            active_env,
            "'. Install it manually or allow automatic installation by omitting python_env/conda_env."
        )
    }

    if (!reticulate::py_module_available("normae.cli")) {
        stop("Python module 'normae.cli' is required but not available.")
    }

    normalizePath(config$python, winslash = "/", mustWork = TRUE)
}

.normae_format_cli_args <- function(normae_args) {
    # ---- constants -----------------------------------------------------------
    # Args handled elsewhere and therefore NOT emitted here
    reserved <- c(
        "meta_csv", "sample_csv", "output_dir",
        "batch_indicator_col", "order_indicator_col",
        "qc_indicator_value", "qc_indicator_col"
    )

    # Defaults to force-emit when "unspecified"
    default_if_unspecified <- list(
        mz_row          = "",
        rt_row          = "",
        early_stop      = FALSE
    )

    # ---- helpers -------------------------------------------------------------
    is_unspecified <- function(x) {
        is.null(x) || length(x) == 0L || all(is.na(x))
    }

    norm_val <- function(x) {
        if (is.null(x)) {
            return(NULL)
        }
        if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
        if (!length(x)) {
            return(NULL)
        }
        x <- x[!is.na(x)]
        if (!length(x)) {
            return(NULL)
        }
        if (is.logical(x)) {
            return(ifelse(x, "True", "False"))
        } else if (is.numeric(x)) {
            return(format(x, scientific = FALSE, trim = TRUE))
        } else {
            return(as.character(x))
        }
    }

    append_arg <- function(acc, name, value_vec) {
        # Append flag + its values (vector allowed)
        c(acc, paste0("--", name), value_vec)
    }

    # ---- validation & normalization of input --------------------------------
    if (is.null(normae_args) || !length(normae_args)) {
        normae_args <- list()
    }
    arg_names <- names(normae_args)
    if (is.null(arg_names) || any(arg_names == "")) {
        stop("All entries in normae_args must be named.")
    }

    # Inject defaults when unspecified
    for (k in names(default_if_unspecified)) {
        if (!k %in% names(normae_args) || is_unspecified(normae_args[[k]])) {
            normae_args[[k]] <- default_if_unspecified[[k]]
        }
    }

    # ---- build CLI args ------------------------------------------------------
    res <- character(0)
    for (nm in names(normae_args)) {
        if (nm %in% reserved) next
        val <- norm_val(normae_args[[nm]])
        if (is.null(val)) next # skip only truly empty after normalization
        res <- append_arg(res, nm, val)
    }
    res
}


.preprocess_qc_column <- function(sample_annotation, qc_col_name, sample_id_col) {
    if (!is.null(qc_col_name)) {
        if (!(qc_col_name %in% names(sample_annotation))) {
            stop("QC column '", qc_col_name, "' not found in sample_annotation.")
        }
        qc_mask <- .normae_qc_mask(sample_annotation[[qc_col_name]])
        if (any(qc_mask)) {
            qc_col_name <- "..normae_qc"
            sample_annotation[[qc_col_name]] <- ifelse(qc_mask, "QC", "Sample")
            qc_indicator_value <- "QC"
        } else {
            message("QC column supplied but no QC samples detected; QC information will be ignored.")
        }
    }
    sample_annotation
}

.normae_qc_mask <- function(values) {
    if (is.logical(values)) {
        return(!is.na(values) & values)
    }
    if (is.factor(values)) {
        values <- as.character(values)
    }
    if (is.character(values)) {
        lower <- tolower(values)
        return(!is.na(lower) & lower %in% c("qc", "quality_control", "qualitycontrol"))
    }
    logical(length(values))
}
