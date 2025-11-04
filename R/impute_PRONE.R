#' Impute missing intensities with PRONE (optional)
#'
#' Adapter around `PRONE::impute_se()` that lets proBatch users run PRONE's
#' mixed kNN/MNAR imputation on:
#'   - long data.frames (one feature per row and sample),
#'   - wide data matrices (features in rows, samples in columns).
#'
#' For `ProBatchFeatures` pipelines, use [pb_transform()] with the registered
#' `"PRONEImpute"` step.
#'
#' @param x Object to impute; one of:
#'   \itemize{
#'     \item long data.frame (proBatch long format, see [proBatch()]);
#'     \item numeric matrix/data.frame with features in rows and samples in columns.
#'   }
#' @param sample_annotation Optional sample annotation table; must contain
#'   `sample_id_col` if provided.
#' @param sample_id_col Column in `sample_annotation` that matches matrix column
#'   names or long-format sample IDs. Default `"FullRunName"`.
#' @param feature_id_col Column in long-format data holding peptide/protein IDs.
#'   Default `"feature_id"`.
#' @param measure_col Column with intensities in long-format data.
#' @param condition_col Optional column name in `sample_annotation` that is
#'   passed as the `condition` argument to `PRONE::impute_se()` (controls group
#'   separation before imputation). If `NULL`, PRONE imputes globally.
#'
#' @return Object of the same class as `x`, with imputed values.
#' @rdname imputePRONE
#' @export
imputePRONE <- function(x,
                        sample_annotation = NULL,
                        sample_id_col = "FullRunName",
                        feature_id_col = "feature_id",
                        measure_col = "Intensity",
                        condition_col = NULL) {
    .pb_requireNamespace("PRONE")

    if (is.matrix(x)) {
        return(imputePRONE_dm(
            x = x,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            condition_col = condition_col,
            assay_in = "tmp_raw"
        ))
    }

    if (is.data.frame(x)) {
        cols <- names(x)
        looks_long <- all(c(feature_id_col, sample_id_col, measure_col) %in% cols)
        if (looks_long) {
            return(imputePRONE_df(
                x = x,
                sample_annotation = sample_annotation,
                sample_id_col = sample_id_col,
                feature_id_col = feature_id_col,
                measure_col = measure_col,
                condition_col = condition_col,
                assay_in = "tmp_raw"
            ))
        }
        return(imputePRONE_dm(
            x = x,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            condition_col = condition_col,
            assay_in = "tmp_raw"
        ))
    }

    if (inherits(x, "ProBatchFeatures") || inherits(x, "QFeatures")) {
        stop(
            "imputePRONE(): ProBatchFeatures/QFeatures inputs are not supported. ",
            "Use pb_transform(..., steps = 'PRONEImpute') to apply PRONE within a pipeline."
        )
    }

    stop("imputePRONE(): unsupported object of class ", paste(class(x), collapse = "/"))
}

# ---- long-format front-end -------------------------------------------------

#' @rdname imputePRONE
#' @export
imputePRONE_df <- function(x,
                           sample_annotation = NULL,
                           sample_id_col = "FullRunName",
                           feature_id_col = "feature_id",
                           measure_col = "Intensity",
                           condition_col = NULL,
                           assay_in = "raw") {
    .pb_requireNamespace("PRONE")

    data_matrix <- long_to_matrix(
        df_long = x,
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col
    )

    imputed_matrix <- .prone_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        condition_col = condition_col,
        assay_in = assay_in
    )

    matrix_to_long(
        data_matrix = imputed_matrix,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col
    )
}

# ---- wide-format front-end -------------------------------------------------

#' @rdname imputePRONE
#' @export
imputePRONE_dm <- function(x,
                           sample_annotation = NULL,
                           sample_id_col = "FullRunName",
                           condition_col = NULL,
                           assay_in = "raw") {
    .pb_requireNamespace("PRONE")

    data_matrix <- if (is.matrix(x)) x else as.matrix(x)

    .prone_matrix_step(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        condition_col = condition_col,
        assay_in = assay_in
    )
}

# ---- shared helper --------------------------------------------------------

.prone_matrix_step <- function(data_matrix,
                               sample_annotation = NULL,
                               sample_id_col = "FullRunName",
                               condition_col = NULL,
                               assay_in = "raw") {
    .pb_requireNamespace("PRONE")

    if (!is.matrix(data_matrix)) {
        data_matrix <- as.matrix(data_matrix)
    }
    storage.mode(data_matrix) <- "double"

    if (is.null(colnames(data_matrix))) {
        stop("PRONE imputation requires matrix column names (sample identifiers).", call. = FALSE)
    }
    sample_ids <- colnames(data_matrix)
    sample_id_col_local <- sample_id_col
    placeholder_col <- NULL
    sample_annotation_local <- sample_annotation

    if (is.null(sample_annotation_local)) {
        placeholder_col <- ".pb_prone_sample_id"
        sample_annotation_local <- data.frame(
            .pb_prone_sample_id = sample_ids,
            stringsAsFactors = FALSE
        )
        sample_id_col_local <- placeholder_col
    }

    .run_matrix_method(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation_local,
        sample_id_col = sample_id_col_local,
        fill_the_missing = NULL,
        missing_warning = "PRONE imputation cannot operate without sample identifiers.",
        method_fun = function(data_matrix, sample_annotation) {
            sample_ids_local <- colnames(data_matrix)
            sample_df <- as.data.frame(sample_annotation)

            if (!is.null(placeholder_col) && placeholder_col %in% names(sample_df)) {
                sample_df[[placeholder_col]] <- NULL
            }

            rownames(sample_df) <- sample_ids_local

            condition_arg <- NULL
            if (!is.null(condition_col)) {
                if (is.character(condition_col) && length(condition_col) == 1L) {
                    if (!(condition_col %in% colnames(sample_df))) {
                        stop(
                            "PRONE imputation: condition column '", condition_col,
                            "' not found in sample annotation.",
                            call. = FALSE
                        )
                    }
                    condition_arg <- condition_col
                } else {
                    cond_values <- condition_col
                    if (!is.null(names(cond_values))) {
                        idx <- match(sample_ids_local, names(cond_values))
                        if (anyNA(idx)) {
                            stop(
                                "PRONE imputation: condition vector is missing values for some samples.",
                                call. = FALSE
                            )
                        }
                        cond_values <- cond_values[idx]
                    } else {
                        if (length(cond_values) != length(sample_ids_local)) {
                            stop(
                                "PRONE imputation: unnamed condition vector must match the number of samples.",
                                call. = FALSE
                            )
                        }
                    }

                    cond_values <- unname(as.vector(cond_values))

                    existing_names <- colnames(sample_df)
                    base_name <- ".pb_prone_condition"
                    new_name <- base_name
                    counter <- 1L
                    while (!is.null(existing_names) && new_name %in% existing_names) {
                        counter <- counter + 1L
                        new_name <- paste0(base_name, "_", counter)
                    }

                    sample_df[[new_name]] <- cond_values
                    condition_arg <- new_name
                }
            }

            col_data <- S4Vectors::DataFrame(sample_df)
            rownames(col_data) <- sample_ids_local
            col_data$Column <- rownames(col_data)

            se <- SummarizedExperiment::SummarizedExperiment(
                assays = setNames(list(data_matrix), assay_in),
                colData = col_data
            )

            prone_impute <- .pb_prone_impute_fun()
            # TODO: switch back to PRONE::impute_se once the upstream bug is fixed.
            se_imp <- prone_impute(
                se,
                ain = assay_in,
                condition = condition_arg
            )

            SummarizedExperiment::assay(se_imp, assay_in)
        }
    )
}

.pb_prone_impute_fun <- local({
    override_fun <- NULL
    function() {
        opt_fun <- getOption("proBatch.prone_impute_se", NULL)
        if (is.function(opt_fun)) {
            return(opt_fun)
        }
        if (is.null(override_fun)) {
            override_fun <<- .pb_load_prone_impute()
        }
        override_fun
    }
})

.pb_load_prone_impute <- function() {
    override_path <- .pb_prone_override_path()
    if (!is.null(override_path)) {
        override_env <- new.env(parent = environment())
        tryCatch(
            {
                sys.source(override_path, envir = override_env)
                fun <- override_env$impute_se
                if (is.function(fun)) {
                    return(fun)
                }
            },
            error = function(e) NULL
        )
    }
    PRONE::impute_se
}

.pb_prone_override_path <- function() {
    inst_path <- system.file("overrides/PRONE/Imputation.R", package = "proBatch")
    if (nzchar(inst_path) && file.exists(inst_path)) {
        return(inst_path)
    }
    pkg_root <- tryCatch(find.package("proBatch"), error = function(e) "")
    if (nzchar(pkg_root)) {
        dev_path <- file.path(pkg_root, "inst", "overrides", "PRONE", "Imputation.R")
        if (file.exists(dev_path)) {
            return(dev_path)
        }
    }
    NULL
}
