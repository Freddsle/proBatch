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
#' @param format One of `"long"` or `"wide"`.
#' @param inj_order_col Column in \code{sample_annotation} with injection order.
#'   If \code{NULL}, the order of \code{sample_annotation[[sample_id_col]]} is used.
#' @param qc_col_name Optional column in \code{sample_annotation} marking QC samples.
#'   Logical, factor (level "QC"), or character ("QC") values are treated as QC.
#' @param fill_the_missing Missing-value policy before invoking NormAE. If \code{NULL},
#'   missing values are left unchanged (and NormAE will fail if any remain). Set \code{FALSE}
#'   to keep NA entries untouched, or provide a numeric value / removal token (e.g., \code{"remove"})
#'   to impute or drop before correction.
#' @param normae_args Named list of CLI parameters forwarded to \code{normae}.
#'   Entries are translated to \code{--name value} pairs. Default: empty list.
#'   Some defaults are set internally if not provided:
#'     --mz_row set to "".
#'     --rt_row set to "".
#'     --early_stop set to False.
#' @param python_env Optional path to a Python binary.
#' @param conda_env Optional conda environment name. Ignored when \code{python_env}
#'   is supplied. If neither is provided, a \code{"normae"} conda env is created.
#'
#' @return Same type as input: numeric \code{matrix} for \code{format="wide"},
#'   or a long \code{data.frame} for \code{format="long"}.
#'
#' @details Requires Python (>=3.10) with the \pkg{normae} package installed.
#'   The first call without \code{python_env}/\code{conda_env} creates a conda
#'   environment named \code{"normae"} and installs the package from GitHub.
#'
#' @references Rong, Zhiwei, et al. "NormAE: deep adversarial learning model to remove batch effects in
#' liquid chromatography mass spectrometry-based metabolomics data." Analytical chemistry 92.7 (2020): 5082-5090.
#' GitHub: https://github.com/luyiyun/NormAE/tree/release
#'
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
  fill_the_missing = NULL,
  keep_all = "default",
  python_env = NULL,
  conda_env = NULL,
  normae_args = list()
) {
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
            fill_the_missing = fill_the_missing,
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

    handled <- .handle_missing_for_batch_df(
        df_long = df_long,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col,
        fill_the_missing = fill_the_missing,
        warning_message = "NormAE cannot operate with missing values in the matrix",
        qual_col = NULL,
        qual_value = NULL
    )
    df_long <- handled$df_long
    sample_annotation <- handled$sample_annotation

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
        fill_the_missing = fill_the_missing,
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
  fill_the_missing = NULL,
  python_env = NULL,
  conda_env = NULL,
  normae_args = list()
) {
    if (is.null(sample_annotation)) {
        stop("sample_annotation must be provided for NormAE correction.")
    }

    .run_matrix_method(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        fill_the_missing = fill_the_missing,
        missing_warning = "NormAE cannot operate with missing values in the matrix",
        method_fun = function(data_matrix, sample_annotation) {
            if (anyNA(data_matrix)) {
                stop("NormAE requires no NAs; impute/filter before calling.")
            }
            if (!is.numeric(data_matrix)) {
                stop("Input must be coercible to a numeric matrix for NormAE correction.")
            }
            if (any(data_matrix < 0, na.rm = TRUE)) {
                message("NormAE expects non-negative intensities; negative values were detected.")
            }

            .run_normae_core(
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
        }
    )
}

.run_normae_core <- function(
  data_matrix, # numeric matrix (features × samples)
  sample_annotation,
  sample_id_col,
  batch_col,
  inj_order_col = NULL,
  qc_col_name = NULL,
  python_env = NULL,
  conda_env = NULL,
  normae_args = list()
) {
    stopifnot(is.matrix(data_matrix))
    if (is.null(sample_annotation)) stop("sample_annotation is required.")

    # --- SA coverage + ordering (minimal, deterministic) ---------------------
    sample_ids <- colnames(data_matrix)
    if (is.null(sample_ids) || !length(sample_ids)) {
        stop("data_matrix must have column names (sample IDs).")
    }
    if (!(sample_id_col %in% names(sample_annotation))) {
        stop("Column '", sample_id_col, "' not in sample_annotation.")
    }
    if (!(batch_col %in% names(sample_annotation))) {
        stop("Column '", batch_col, "' not in sample_annotation.")
    }

    # Keep only rows for present samples and order to match data_matrix
    sa <- sample_annotation[sample_annotation[[sample_id_col]] %in% sample_ids, , drop = FALSE]
    sa[[sample_id_col]] <- as.character(sa[[sample_id_col]])
    sa <- sa[match(sample_ids, sa[[sample_id_col]]), , drop = FALSE]
    if (anyNA(sa[[sample_id_col]])) {
        stop("Some matrix samples are absent in sample_annotation.")
    }

    # --- Injection order (if missing, synthesize 1..n) -----------------------
    inj_res <- .normae_prepare_injection_order(sa, sample_id_col, inj_order_col)
    sa <- inj_res$sample_annotation
    inj_col_name <- inj_res$inj_col_name

    # --- QC (optional, tolerant) ---------------------------------------------
    qc_res <- .preprocess_qc_column(sa, qc_col_name, sample_id_col)
    sa <- qc_res$sample_annotation
    qc_col_name <- qc_res$qc_col_name
    qc_indicator_value <- qc_res$qc_indicator_value

    # --- Python resolver (reticulate first, basilisk as fallback) ------------
    python_bin <- .normae_prepare_python(python_env, conda_env)

    # --- I/O in the *working directory* --------------------------------------
    # Prepare input CSV files for NormAE CLI
    meta_df <- as.data.frame(t(data_matrix), check.names = FALSE)
    rownames(meta_df) <- sample_ids
    run_tag <- paste0("normae_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", as.integer(runif(1, 1, 1e9)))
    run_dir <- file.path(getwd(), run_tag)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

    meta_path <- file.path(run_dir, "meta.csv")
    sample_path <- file.path(run_dir, "samples.csv")
    out_dir <- file.path(run_dir, "output")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    write.csv(meta_df, meta_path, quote = TRUE)
    cols_needed <- unique(c(
        sample_id_col, batch_col,
        if (!is.null(inj_col_name)) inj_col_name,
        if (!is.null(qc_col_name)) qc_col_name
    ))
    sa_write <- sa[, cols_needed, drop = FALSE]
    # Convert factors to characters to avoid writing factor levels
    for (nm in names(sa_write)) {
        if (is.factor(sa_write[[nm]])) {
            sa_write[[nm]] <- as.character(sa_write[[nm]])
        }
    }
    write.csv(sa_write, sample_path, quote = TRUE, row.names = FALSE)

    # --- CLI args -----------------------------------------------------------
    cli_args <- c(
        "--meta_csv", normalizePath(meta_path, winslash = "/", mustWork = TRUE),
        "--sample_csv", normalizePath(sample_path, winslash = "/", mustWork = TRUE),
        "--output_dir", normalizePath(out_dir, winslash = "/", mustWork = TRUE),
        "--batch_indicator_col", batch_col
    )
    if (!is.null(inj_col_name) && nzchar(inj_col_name)) {
        cli_args <- c(cli_args, "--order_indicator_col", inj_col_name)
    }
    if (!is.null(qc_col_name) && nzchar(qc_col_name)) {
        cli_args <- c(
            cli_args, "--qc_indicator_col", qc_col_name,
            "--qc_indicator_value", if (is.null(qc_indicator_value)) "QC" else qc_indicator_value
        )
    }
    # Append any additional NormAE CLI arguments
    cli_args <- unname(c(cli_args, .normae_format_cli_args(normae_args)))

    # --- Run NormAE CLI via system call -------------------------------------
    out <- system2(python_bin, c("-m", "normae.cli", cli_args), stdout = TRUE, stderr = TRUE)
    status <- attr(out, "status")
    if (is.null(status)) status <- 0L
    if (!identical(status, 0L)) {
        stop("NormAE CLI failed.\n", paste(out, collapse = "\n"), call. = FALSE)
    }

    # --- Read result and reconcile dimensions -------------------------------
    res_csv <- file.path(out_dir, "X_clean.csv")
    if (!file.exists(res_csv)) {
        stop("NormAE did not produce 'X_clean.csv' in: ", out_dir, call. = FALSE)
    }
    res_df <- read.csv(res_csv, check.names = FALSE, row.names = 1)
    res_mat <- as.matrix(res_df)
    storage.mode(res_mat) <- "double"

    # Check orientation and align output to input matrix dimensions
    aligned <- if (identical(dim(res_mat), dim(data_matrix)) &&
        identical(colnames(res_mat), colnames(data_matrix))) {
        res_mat
    } else if (identical(rownames(res_mat), colnames(data_matrix))) {
        t(res_mat)
    } else {
        stop("Unexpected NormAE output shape; cannot align to input.")
    }
    dimnames(aligned) <- dimnames(data_matrix)

    # Clean up temporary directory
    unlink(run_dir, recursive = TRUE)
    aligned
}

.normae_prepare_python <- function(python_env = NULL, conda_env = NULL) {
    # Ensure at least one Python interface package is available (reticulate or basilisk)
    if (!.pb_requireNamespace("reticulate")) {
        if (!.pb_requireNamespace("basilisk")) {
            stop("NormAE correction requires either the 'reticulate' package or 'basilisk' package.")
        }
        # Fallback to basilisk: create/use a Python environment with NormAE installed
        env <- basilisk::BasiliskEnvironment(
            envname = "normae", pkgname = "proBatch",
            packages = c("python=3.10"),
            pip = c("git+https://github.com/luyiyun/NormAE.git")
        )
        env_path <- basilisk::obtainEnvironmentPath(env)
        py <- basilisk::getPythonBinary(env_path)
        return(py)
    }

    # If we reach here, reticulate is available. Use reticulate's facilities:
    has_cli <- function(py) {
        if (!nzchar(py) || !file.exists(py) || dir.exists(py)) {
            return(FALSE)
        }
        out <- tryCatch(
            system2(py, c("-c", "import importlib.util as u; print('OK' if u.find_spec('normae.cli') else 'NO')"),
                stdout = TRUE, stderr = TRUE
            ),
            error = function(e) character()
        )
        any(grepl("^OK$", out))
    }

    ensure_normae <- function(py) {
        if (has_cli(py)) {
            return(py)
        }
        # Upgrade pip and install NormAE from GitHub
        system2(py, c("-m", "pip", "install", "--upgrade", "pip"), stdout = TRUE, stderr = TRUE)
        cmd <- c("-m", "pip", "install", "git+https://github.com/luyiyun/NormAE.git")
        out <- system2(py, cmd, stdout = TRUE, stderr = TRUE)
        if (!has_cli(py)) {
            stop("Failed to install 'normae' into Python at: ", py,
                "\nOutput:\n", paste(out, collapse = "\n"),
                call. = FALSE
            )
        }
        py
    }

    # 1) If a specific Python binary path is provided, use it
    if (!is.null(python_env)) {
        py <- normalizePath(python_env, winslash = "/", mustWork = FALSE)
        if (!file.exists(py)) {
            stop("Provided python_env does not exist: ", python_env, call. = FALSE)
        }
        return(ensure_normae(py))
    }

    # 2) If a specific conda environment name is provided, use or create it
    if (!is.null(conda_env)) {
        py <- tryCatch(reticulate::conda_python(conda_env, conda = "auto"),
            error = function(e) NULL
        )
        if (is.null(py)) {
            # Create env if missing
            py <- reticulate::conda_create(envname = conda_env, python_version = "3.10", conda = "auto")
        }
        return(ensure_normae(py))
    }

    # 3) If system default Python already has NormAE CLI, use it
    path_py <- Sys.which("python")
    if (nzchar(path_py) && has_cli(path_py)) {
        return(normalizePath(path_py, winslash = "/", mustWork = FALSE))
    }

    # 4) Otherwise, create or use a conda environment named "normae" via reticulate
    if (is.na(reticulate::conda_binary(conda = "auto"))) {
        # Install Miniconda if no conda is available
        reticulate::install_miniconda()
    }
    py <- reticulate::conda_create(envname = "normae", python_version = "3.10", conda = "auto")
    py <- ensure_normae(py)
    if (has_cli(py)) {
        return(py)
    }

    # 5) Last resort: install Python and use virtualenv (for non-Windows platforms)
    reticulate::install_python(version = "3.10")
    if (!"normae" %in% reticulate::virtualenv_list()) {
        reticulate::virtualenv_create("normae")
    }
    vpy <- reticulate::virtualenv_python("normae")
    ensure_normae(vpy)
}

.normae_format_cli_args <- function(normae_args) {
    if (is.null(normae_args)) normae_args <- list()
    # Set internal defaults if not overridden by user
    if (!("mz_row" %in% names(normae_args))) normae_args$mz_row <- ""
    if (!("rt_row" %in% names(normae_args))) normae_args$rt_row <- ""
    if (!("early_stop" %in% names(normae_args))) normae_args$early_stop <- FALSE

    reserved <- c(
        "meta_csv", "sample_csv", "output_dir",
        "batch_indicator_col", "order_indicator_col",
        "qc_indicator_col", "qc_indicator_value"
    )

    args_list <- character(0)
    for (nm in setdiff(names(normae_args), reserved)) {
        val <- normae_args[[nm]]
        if (is.null(val)) next
        # Flatten simple lists into a vector of values
        if (is.list(val)) val <- unlist(val, recursive = TRUE, use.names = FALSE)
        if (!length(val)) next
        # Convert logicals to Python booleans as strings
        val <- vapply(val, function(v) {
            if (isTRUE(v)) "True" else if (identical(v, FALSE)) "False" else as.character(v)
        }, character(1))
        args_list <- c(args_list, paste0("--", nm), val)
    }
    unname(args_list)
}

.normae_prepare_injection_order <- function(sample_annotation, sample_id_col, inj_order_col) {
    if (!(sample_id_col %in% names(sample_annotation))) {
        stop("Sample ID column '", sample_id_col, "' not found in sample_annotation.")
    }

    if (is.null(inj_order_col) || identical(inj_order_col, "")) {
        inj_col_name <- "..normae_order"
        sample_annotation[[inj_col_name]] <- seq_len(nrow(sample_annotation))
        return(list(
            sample_annotation = sample_annotation,
            inj_col_name = inj_col_name
        ))
    }

    if (!(inj_order_col %in% names(sample_annotation))) {
        stop("Injection order column '", inj_order_col, "' not found in sample_annotation.")
    }

    inj_vals <- sample_annotation[[inj_order_col]]
    if (anyNA(inj_vals)) {
        stop("Injection order column '", inj_order_col, "' contains missing values.")
    }
    if (!is.numeric(inj_vals)) {
        suppressWarnings(num_vals <- as.numeric(as.character(inj_vals)))
        if (anyNA(num_vals)) {
            stop("Injection order column '", inj_order_col, "' must be numeric or coercible to numeric.")
        }
        inj_vals <- num_vals
    }
    # Convert to integer if all values are essentially integers
    if (all(inj_vals == round(inj_vals))) {
        inj_vals <- as.integer(round(inj_vals))
    }
    # Update and return
    sample_annotation[[inj_order_col]] <- inj_vals
    list(
        sample_annotation = sample_annotation,
        inj_col_name = inj_order_col
    )
}

.preprocess_qc_column <- function(sample_annotation, qc_col_name, sample_id_col) {
    qc_indicator_value <- NULL
    if (!(sample_id_col %in% names(sample_annotation))) {
        stop("Sample ID column '", sample_id_col, "' not found in sample_annotation.")
    }
    if (!is.null(qc_col_name)) {
        if (!(qc_col_name %in% names(sample_annotation))) {
            stop("QC column '", qc_col_name, "' not found in sample_annotation.")
        }
        qc_mask <- .normae_qc_mask(sample_annotation[[qc_col_name]])
        if (any(qc_mask)) {
            new_col <- "..normae_qc"
            sample_annotation[[new_col]] <- ifelse(qc_mask, "QC", "Subject")
            qc_col_name <- new_col
            qc_indicator_value <- "QC"
        } else {
            message("QC column supplied but no QC samples detected; QC information will be ignored.")
            qc_col_name <- NULL
        }
    }
    list(
        sample_annotation = sample_annotation,
        qc_col_name = qc_col_name,
        qc_indicator_value = qc_indicator_value
    )
}

.normae_qc_mask <- function(values) {
    if (is.logical(values)) {
        return(!is.na(values) & values)
    }
    if (is.factor(values)) values <- as.character(values)
    if (is.character(values)) {
        lower <- tolower(values)
        return(!is.na(lower) & lower %in% c("qc", "quality_control", "qualitycontrol"))
    }
    # Default: no values qualify as QC (non-matching type)
    logical(length(values))
}
