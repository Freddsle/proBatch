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
#'   If \code{NULL}, the column is omitted.
#'   If \code{"sequential"}, the order of \code{sample_annotation[[sample_id_col]]} is used.
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
#' @param conda_env Optional conda/mamba environment name. Ignored when \code{python_env}
#'   is supplied. If neither is provided, a \code{"normae"} conda env is created - NOT TESTED.
#' @param conda_env_path Optional path to an existing conda/mamba environment. When provided,
#'   the environment is treated as ready with NormAE installed and used as-is.
#' @param python_env Optional path to a Python binary.
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
  conda_env_path = NULL,
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
            conda_env_path = conda_env_path,
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
        conda_env_path = conda_env_path,
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
  conda_env_path = NULL,
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
                conda_env_path = conda_env_path,
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
  conda_env_path = NULL,
  normae_args = list()
) {
    stopifnot(is.matrix(data_matrix))
    if (is.null(sample_annotation)) stop("sample_annotation is required.")

    # --- SA coverage + ordering ------------------------------------------------
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
    sa <- sample_annotation[sample_annotation[[sample_id_col]] %in% sample_ids, , drop = FALSE]
    sa[[sample_id_col]] <- as.character(sa[[sample_id_col]])
    sa <- sa[match(sample_ids, sa[[sample_id_col]]), , drop = FALSE]
    if (anyNA(sa[[sample_id_col]])) stop("Some matrix samples are absent in sample_annotation.")

    # --- Injection order + QC ---------------------------------------------------
    if (is.null(inj_order_col)) {
        inj_col_name <- "''"
    } else {
        inj_res <- .normae_prepare_injection_order(sa, sample_id_col, inj_order_col)
        sa <- inj_res$sample_annotation
        inj_col_name <- inj_res$inj_col_name
    }

    qc_res <- .preprocess_qc_column(sa, qc_col_name, sample_id_col)
    sa <- qc_res$sample_annotation
    qc_col_name <- qc_res$qc_col_name
    qc_indicator_value <- qc_res$qc_indicator_value

    # --- Python resolver --------------------------------------------------------
    python_bin <- .normae_prepare_python(python_env, conda_env, conda_env_path)

    # --- I/O: create isolated run dir ------------------------------------------
    meta_df <- as.data.frame(t(data_matrix), check.names = FALSE)
    rownames(meta_df) <- sample_ids
    run_tag <- paste0("normae_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", as.integer(runif(1, 1, 1e9)))
    run_dir <- file.path(getwd(), run_tag)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

    meta_path <- file.path(run_dir, "meta.csv")
    sample_path <- file.path(run_dir, "samples.csv")
    out_dir <- file.path(run_dir, "output")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    write.csv(meta_df, meta_path, quote = TRUE, row.names = TRUE) # explicit
    cols_needed <- unique(c(
        sample_id_col, batch_col,
        if (!is.null(inj_col_name) && inj_col_name != "''") inj_col_name,
        if (!is.null(qc_col_name)) qc_col_name
    ))
    sa_write <- sa[, cols_needed, drop = FALSE]
    for (nm in names(sa_write)) if (is.factor(sa_write[[nm]])) sa_write[[nm]] <- as.character(sa_write[[nm]])
    write.csv(sa_write, sample_path, quote = TRUE, row.names = FALSE)

    # --- CLI args ---------------------------------------------------------------
    cli_args <- c(
        "--meta_csv", normalizePath(meta_path, winslash = "/", mustWork = TRUE),
        "--sample_csv", normalizePath(sample_path, winslash = "/", mustWork = TRUE),
        "--output_dir", normalizePath(out_dir, winslash = "/", mustWork = TRUE),
        "--batch_indicator_col", batch_col
    )
    if (!is.null(inj_col_name) && nzchar(inj_col_name)) {
        cli_args <- c(cli_args, "--order_indicator_col", inj_col_name)
    }
    if (!is.null(qc_col_name) && nzchar(trimws(qc_col_name))) {
        cli_args <- c(
            cli_args,
            "--qc_indicator_col", trimws(qc_col_name),
            "--qc_indicator_value", if (is.null(qc_indicator_value)) "QC" else qc_indicator_value
        )
    } else {
        # Explicitly pass empty QC column to the CLI when no QC column is provided
        cli_args <- c(cli_args, "--qc_indicator_col", "''")
    }
    cli_args <- unname(c(cli_args, .normae_format_cli_args(normae_args)))

    env <- .normae_overlay_env()

    # ---- Prefer the overlay script if present (guaranteed path execution) -----
    overlay_dir <- .normae_overlay_dir()
    if (!is.null(overlay_dir)) {
        overlay_main <- file.path(overlay_dir, "__main__.py")
        cmd <- c(overlay_main, cli_args)
        cat("Running command:", shQuote(python_bin), paste(cmd, collapse = " "), "\n")
        stat <- system2(python_bin, cmd, stdout = "", stderr = "", env = env)
        status <- if (length(stat)) as.integer(stat[[1L]]) else 0L
        if (!identical(status, 0L)) {
            # capture overlay error for diagnostics, then fallback to package CLIs
            out_overlay <- tryCatch(
                system2(python_bin, c(overlay_main, cli_args), stdout = TRUE, stderr = TRUE, env = env),
                error = function(e) paste("Failed to capture overlay output:", conditionMessage(e))
            )
            message("Overlay run failed; falling back to installed 'normae' module.\n", paste(out_overlay, collapse = "\n"))
        } else {
            # success path; proceed to read X_clean.csv
            res_csv <- file.path(out_dir, "X_clean.csv")
            if (!file.exists(res_csv)) {
                stop("NormAE overlay finished but did not produce 'X_clean.csv' in: ", out_dir, call. = FALSE)
            }
            res_df <- read.csv(res_csv, check.names = FALSE, row.names = 1)
            res_mat <- as.matrix(res_df)
            storage.mode(res_mat) <- "double"
            aligned <- if (identical(dim(res_mat), dim(data_matrix)) &&
                identical(colnames(res_mat), colnames(data_matrix))) {
                res_mat
            } else if (identical(rownames(res_mat), colnames(data_matrix))) {
                t(res_mat)
            } else {
                stop("Unexpected NormAE overlay output shape; cannot align to input.")
            }
            dimnames(aligned) <- dimnames(data_matrix)
            return(aligned)
        }
    }

    # ---- Fallbacks: installed CLI (may plot PCA; last resort) ------------------
    run_try <- function(mod) {
        cmd <- c("-m", mod, cli_args)
        cat("Running command:", shQuote(python_bin), paste(cmd, collapse = " "), "\n")
        stat <- system2(python_bin, cmd, stdout = "", stderr = "", env = env)
        status <- if (length(stat)) as.integer(stat[[1L]]) else 0L
        if (is.na(status)) NA_integer_ else status
    }
    status <- run_try("normae")
    if (!identical(status, 0L)) {
        out <- tryCatch(
            system2(python_bin, c("-m", "normae", cli_args), stdout = TRUE, stderr = TRUE, env = env),
            error = function(e) paste("Failed to capture output:", conditionMessage(e))
        )
        status2 <- run_try("normae.cli")
        if (!identical(status2, 0L)) {
            out2 <- tryCatch(
                system2(python_bin, c("-m", "normae.cli", cli_args), stdout = TRUE, stderr = TRUE, env = env),
                error = function(e) paste("Failed to capture output:", conditionMessage(e))
            )
            unlink(run_dir, recursive = TRUE)
            stop("NormAE CLI failed.\n---- Attempt 1: python -m normae ----\n",
                paste(out, collapse = "\n"),
                "\n---- Attempt 2: python -m normae.cli ----\n",
                paste(out2, collapse = "\n"),
                call. = FALSE
            )
        }
    }

    # --- Read result & align ----------------------------------------------------
    res_csv <- file.path(out_dir, "X_clean.csv")
    if (!file.exists(res_csv)) {
        # unlink(run_dir, recursive = TRUE)
        stop("NormAE did not produce 'X_clean.csv' in: ", out_dir, call. = FALSE)
    }
    res_df <- read.csv(res_csv, check.names = FALSE, row.names = 1)
    res_mat <- as.matrix(res_df)
    storage.mode(res_mat) <- "double"

    aligned <- if (identical(dim(res_mat), dim(data_matrix)) &&
        identical(colnames(res_mat), colnames(data_matrix))) {
        res_mat
    } else if (identical(rownames(res_mat), colnames(data_matrix))) {
        t(res_mat)
    } else {
        # unlink(run_dir, recursive = TRUE)
        stop("Unexpected NormAE output shape; cannot align to input.")
    }
    dimnames(aligned) <- dimnames(data_matrix)
    # unlink(run_dir, recursive = TRUE)
    aligned
}

.normae_cli_requirement <- function() {
    "normae @ git+https://github.com/luyiyun/NormAE.git@release"
}

.normae_python_basilisk <- function() {
    if (!.pb_requireNamespace("basilisk")) {
        stop("NormAE requires either 'reticulate' or 'basilisk'; neither is available.")
    }
    # Bioconductor recommends disabling strict version checks while discovering versions.
    basilisk::setBasiliskCheckVersions(FALSE)
    env <- basilisk::BasiliskEnvironment(
        envname = "normae_basilisk",
        pkgname = "proBatch",
        packages = c("python=3.10"),
        pip = c(.normae_cli_requirement())
    )
    env_path <- basilisk::obtainEnvironmentPath(env)
    py <- basilisk::getPythonBinary(env_path)
    py
}

.normae_prepare_python <- function(python_env = NULL, conda_env = NULL, conda_env_path = NULL) {
    basilisk_fallback <- function() .normae_python_basilisk()
    pip_spec <- .normae_cli_requirement()

    has_cli <- function(py) {
        if (is.null(py) || !nzchar(py) || !file.exists(py)) {
            return(FALSE)
        }
        probe <- tempfile("normae-cli-probe-", fileext = ".py")
        on.exit(unlink(probe), add = TRUE)
        # Look for either `normae.__main__` (python -m normae) or a CLI module.
        writeLines(
            c(
                "import importlib.util, sys",
                "mods = ('normae.__main__','normae.cli')",
                "ok = any(importlib.util.find_spec(m) is not None for m in mods)",
                "sys.exit(0 if ok else 1)"
            ),
            probe
        )
        status <- suppressWarnings(tryCatch(system2(py, probe, stdout = FALSE, stderr = FALSE),
            warning = function(w) NA_integer_, error = function(e) NA_integer_
        ))
        if (length(status) != 1L || is.na(status)) {
            return(FALSE)
        }
        identical(as.integer(status), 0L)
    }

    ensure_normae <- function(py) {
        if (has_cli(py)) {
            return(py)
        }
        # Try to get pip up-to-date and install NormAE CLI
        try(system2(py, c("-m", "pip", "install", "--upgrade", "pip"), stdout = TRUE, stderr = TRUE), silent = TRUE)
        out <- tryCatch(
            system2(py, c("-m", "pip", "install", pip_spec), stdout = TRUE, stderr = TRUE),
            error = function(e) character()
        )
        if (has_cli(py)) {
            return(py)
        }
        stop("Failed to install 'normae' into Python (", py, ").\n", paste(out, collapse = "\n"), call. = FALSE)
    }

    # 1) Explicit Python binary
    if (!is.null(python_env)) {
        py <- normalizePath(python_env, winslash = "/", mustWork = FALSE)
        if (!file.exists(py)) stop("Provided python_env does not exist: ", python_env)
        return(ensure_normae(py))
    }

    # 2) Explicit conda env path (assume preinstalled)
    if (!is.null(conda_env_path)) {
        env_path <- normalizePath(conda_env_path, winslash = "/", mustWork = FALSE)
        if (!nzchar(env_path) || (!dir.exists(env_path) && !file.exists(env_path))) {
            stop("Provided conda_env_path does not exist: ", conda_env_path)
        }
        candidates <- if (dir.exists(env_path)) {
            c(
                file.path(env_path, "bin", "python"),
                file.path(env_path, "bin", "python3"),
                file.path(env_path, "bin", "python.exe"),
                file.path(env_path, "python"),
                file.path(env_path, "python3"),
                file.path(env_path, "python.exe"),
                file.path(env_path, "Scripts", "python.exe"),
                file.path(env_path, "Scripts", "python")
            )
        } else {
            env_path
        }
        candidates <- candidates[file.exists(candidates)]
        if (!length(candidates)) stop("Could not locate a Python binary under conda_env_path: ", conda_env_path)
        py <- normalizePath(candidates[[1]], winslash = "/", mustWork = FALSE)
        if (!has_cli(py)) stop("Provided conda_env_path does not contain a NormAE CLI installation.")
        return(py)
    }

    # 3) If reticulate is unavailable -> basilisk immediately
    if (!.pb_requireNamespace("reticulate")) {
        if (!.pb_requireNamespace("basilisk")) {
            stop("NormAE correction requires either 'reticulate' or 'basilisk'.")
        }
        return(basilisk_fallback())
    }

    # 4) Reticulate path (single conda binary pinned across calls to avoid env lookup mismatch)
    conda_bin <- tryCatch(reticulate::conda_binary("auto"), error = function(e) NULL)
    conda_arg <- function() if (!is.null(conda_bin) && !is.na(conda_bin) && nzchar(conda_bin)) conda_bin else "auto"

    py_final <- tryCatch(
        {
            # Named conda env
            if (!is.null(conda_env)) {
                py <- tryCatch(reticulate::conda_python(conda_env, conda = conda_arg()), error = function(e) NULL)
                if (is.null(py)) {
                    reticulate::conda_create(envname = conda_env, python_version = "3.10", conda = conda_arg())
                    py <- tryCatch(reticulate::conda_python(conda_env, conda = conda_arg()), error = function(e) NULL)
                }
                if (is.null(py)) stop("Unable to resolve Python binary for conda environment '", conda_env, "'.")
                return(ensure_normae(py))
            }

            # System python already OK?
            path_py <- Sys.which("python")
            if (nzchar(path_py) && has_cli(path_py)) {
                return(normalizePath(path_py, winslash = "/", mustWork = FALSE))
            }

            # Default 'normae' env via conda
            if (is.null(conda_bin) || is.na(conda_bin)) {
                reticulate::install_miniconda()
                conda_bin <- reticulate::conda_binary("auto")
            }
            py <- tryCatch(reticulate::conda_python("normae", conda = conda_arg()), error = function(e) NULL)
            if (is.null(py)) {
                reticulate::conda_create(envname = "normae", python_version = "3.10", conda = conda_arg())
                py <- tryCatch(reticulate::conda_python("normae", conda = conda_arg()), error = function(e) NULL)
            }
            if (is.null(py)) stop("Unable to resolve Python binary for default 'normae' conda environment.")
            ensure_normae(py)
        },
        error = function(e) {
            if (.pb_requireNamespace("basilisk")) {
                warning("Reticulate setup for NormAE failed (", conditionMessage(e), "); using Basilisk-managed environment.")
                return(basilisk_fallback())
            }
            stop(e)
        }
    )

    if (!has_cli(py_final)) {
        return(basilisk_fallback())
    }
    py_final
}

.normae_format_cli_args <- function(normae_args) {
    if (is.null(normae_args)) normae_args <- list()
    # Set internal defaults if not overridden by user
    if (!("mz_row" %in% names(normae_args))) normae_args$mz_row <- "''"
    if (!("rt_row" %in% names(normae_args))) normae_args$rt_row <- "''"

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

    if (inj_order_col == "sequential") {
        inj_col_name <- "normae_order"
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


.normae_overlay_dir <- function() {
    # 1) explicit option
    opt <- getOption("proBatch.normae_overlay_dir", NULL)
    if (!is.null(opt) && dir.exists(opt)) {
        return(normalizePath(opt, winslash = "/", mustWork = FALSE))
    }
    # 2) env var
    env <- Sys.getenv("PROBATCH_NORMAE_OVERLAY", "")
    if (nzchar(env) && dir.exists(env)) {
        return(normalizePath(env, winslash = "/", mustWork = FALSE))
    }

    # 3) installed package layout: inst/ is flattened; files land under <lib>/proBatch/overrides/normae
    root <- tryCatch(find.package("proBatch"), error = function(e) "")
    if (nzchar(root)) {
        inst_overrides <- file.path(root, "overrides", "normae")
        if (dir.exists(inst_overrides) && file.exists(file.path(inst_overrides, "__main__.py"))) {
            return(normalizePath(inst_overrides, winslash = "/", mustWork = FALSE))
        }
        # 4) development tree layout: keep files under inst/
        dev_overrides <- file.path(root, "inst", "overrides", "normae")
        if (dir.exists(dev_overrides) && file.exists(file.path(dev_overrides, "__main__.py"))) {
            return(normalizePath(dev_overrides, winslash = "/", mustWork = FALSE))
        }
    }
}


.normae_overlay_env <- function(base_env = character(0)) {
    od <- .normae_overlay_dir()
    if (is.null(od)) {
        return(base_env)
    }
    main_py <- file.path(od, "__main__.py")
    if (!file.exists(main_py)) {
        stop("Configured NormAE overlay dir lacks '__main__.py': ", od)
    }
    pkg_root <- normalizePath(dirname(od), winslash = "/", mustWork = FALSE)
    old_pp <- Sys.getenv("PYTHONPATH", "")
    new_pp <- if (nzchar(old_pp)) paste(pkg_root, old_pp, sep = .Platform$path.sep) else pkg_root
    c(base_env, paste0("PYTHONPATH=", new_pp))
}
