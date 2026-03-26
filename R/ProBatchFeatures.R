#' ProBatchFeatures: QFeatures subclass with operation log, levels/pipelines, and lazy storage
#'
#' Assay naming convention:
#'   "<level>::<pipeline>"  e.g., "peptide::raw", "protein::median_on_log"
#' Pipelines are strings produced by get_chain(as_string=TRUE), e.g., "combat_on_medianNorm_on_log".
#'
#' Ephemeral "fast" steps are computed but not stored by default (store_fast_steps = FALSE).
#' Use pb_eval() to compute and return data after a step/pipeline without storing.
#' Use pb_transform() to build pipelines and optionally materialize the final assay.
#'
#' @slot chain character() ordered list of steps (e.g., c("log","medianNorm","combat")).
#' @slot oplog S4Vectors::DataFrame with columns:
#' - step (character), fun (character), from (character), to (character),
#' params (list), timestamp (POSIXct), pkg (character)
#'
#' @example inst/examples/ProBatchFeatures-basic.R
#'
#' @import QFeatures
#' @import SummarizedExperiment
#' @import methods
#' @importClassesFrom S4Vectors DataFrame
#' @exportClass ProBatchFeatures
setClass(
    "ProBatchFeatures",
    contains = "QFeatures",
    slots = list(
        chain = "character", # ordered processing steps (last is most recent)
        oplog = "DataFrame" # structured operation log
    )
)

setValidity("ProBatchFeatures", function(object) {
    # Required columns in oplog
    req <- c("step", "fun", "from", "to", "params", "timestamp", "pkg")
    if (!inherits(object@oplog, "DataFrame")) {
        return("oplog must be an S4Vectors::DataFrame")
    }
    if (!all(req %in% colnames(object@oplog))) {
        return("oplog lacks required columns")
    }
    # Light sanity checks
    if (nrow(object@oplog)) {
        if (!is.list(object@oplog$params)) {
            return("oplog$params must be a list column")
        }
        if (!inherits(object@oplog$timestamp, "POSIXct")) {
            return("oplog$timestamp must be POSIXct")
        }
    }
    if (!is.character(object@chain)) {
        return("chain must be character()")
    }
    TRUE
})

# ---------------------------
# Small utilities
# ---------------------------

.pb_now <- function() as.POSIXct(Sys.time(), tz = "UTC")

.pb_make_pipeline_name <- function(steps) {
    if (!length(steps)) {
        return("raw")
    }
    paste(rev(steps), collapse = "_on_")
}

.pb_step_label <- function(step, params) {
    if (!identical(step, "BERT")) {
        return(step)
    }
    method <- NULL
    if (!is.null(params)) {
        if (is.list(params) && length(params)) {
            method <- params[["bert_method"]] %||% params[["method"]]
        }
    }
    method <- tolower(as.character(method %||% "combat"))
    mapping <- c(
        combat = "BERTc",
        limma = "BERTl",
        ref = "BERTr"
    )
    mapped <- mapping[[method]]
    if (!is.null(mapped)) {
        return(mapped)
    }
    if (!nzchar(method) || identical(method, "bert")) {
        return("BERT")
    }
    paste0("BERT", method)
}

.pb_assay_name <- function(level, pipeline) {
    level <- if (is.null(level) || !nzchar(level)) "feature" else level
    pipeline <- if (is.null(pipeline) || !nzchar(pipeline)) "raw" else pipeline
    paste0(level, "::", pipeline)
}

.pb_is_fast_step <- function(step, fast_steps) {
    isTRUE(step %in% fast_steps)
}

# registry for step functions (can be extended at runtime)
.pb_step_registry <- local({
    reg <- new.env(parent = emptyenv())
    # sensible defaults; replace with proBatch functions if you prefer
    reg$log2 <- function(m, pseudo = 1) log2(m + pseudo)
    reg$log <- function(m, base = exp(1), pseudo = 1) log(m + pseudo, base = base)
    reg$medianNorm <- function(m,
                               sample_annotation = NULL,
                               sample_id_col = "FullRunName",
                               group_col = NULL,
                               inside_batch = FALSE,
                               fill_the_missing = NULL) {
        handle_flag <- !is.null(fill_the_missing) || identical(fill_the_missing, FALSE)
        if (handle_flag && anyNA(m)) {
            m <- handle_missing_values(
                data_matrix = m,
                warning_message = "Median normalization: applying requested missing value handling before centering.",
                fill_the_missing = fill_the_missing
            )
            if (!nrow(m) || !ncol(m)) {
                stop("No data remaining after handling missing values for median normalization")
            }
        }
        normalize_sample_medians_dm(
            data_matrix = m,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            group_col = group_col,
            inside_batch = inside_batch
        )
    }
    reg
})

#' Allow to register/override steps at runtime (e.g., map "combat" -> proBatch::combat_dm)
#' @param name character(1) step name
#' @param fun function implementing the step
#' @return NULL (invisible)
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
pb_register_step <- function(name, fun) {
    stopifnot(is.character(name), length(name) == 1, is.function(fun))
    .pb_step_registry[[name]] <- fun
    invisible(TRUE)
}

.pb_get_step_fun <- function(fun_or_name) {
    if (is.function(fun_or_name)) {
        return(fun_or_name)
    }
    if (is.character(fun_or_name) && exists(fun_or_name, envir = .pb_step_registry, inherits = FALSE)) {
        return(get(fun_or_name, envir = .pb_step_registry, inherits = FALSE))
    }
    # fallbacks: search proBatch namespace first, then match.fun
    if (is.character(fun_or_name) && requireNamespace("proBatch", quietly = TRUE)) {
        if (exists(fun_or_name, envir = asNamespace("proBatch"), inherits = FALSE)) {
            return(get(fun_or_name, envir = asNamespace("proBatch"), inherits = FALSE))
        }
    }
    match.fun(fun_or_name)
}

# Choose materialization backend
.pb_materialize_matrix <- function(m, backend = c("auto", "memory", "hdf5"), hdf5_path = NULL) {
    backend <- match.arg(backend)
    if (backend == "memory") {
        return(m)
    }
    if (backend == "hdf5" || backend == "auto") {
        if (!requireNamespace("HDF5Array", quietly = TRUE)) {
            # fallback to memory if HDF5Array not present
            return(m)
        }
        # heuristic: if auto, write to HDF5 for large matrices
        use_hdf5 <- (backend == "hdf5") ||
            (backend == "auto" && (prod(dim(m)) > 5e7L)) # ~50M cells heuristic
        if (!use_hdf5) {
            return(m)
        }
        # writeHDF5Array returns a DelayedMatrix backed by HDF5
        if (is.null(hdf5_path)) {
            return(HDF5Array::writeHDF5Array(m))
        } else {
            return(HDF5Array::writeHDF5Array(m, filepath = hdf5_path))
        }
    }
    m
}

# Harmonize/validate colData before addAssay
.pb_harmonize_colData <- function(object, se, from_assay) {
    obj_cd <- DataFrame(colData(object))
    se_cd <- DataFrame(colData(se))

    samp_obj <- rownames(obj_cd)
    samp_se <- rownames(se_cd)

    if (!setequal(samp_obj, samp_se)) {
        stop(
            "Sample sets differ between object and incoming assay. ",
            "Missing in object: ", paste(setdiff(samp_se, samp_obj), collapse = ", "),
            "; missing in assay: ", paste(setdiff(samp_obj, samp_se), collapse = ", ")
        )
    }

    # Align order
    obj_cd <- obj_cd[samp_se, , drop = FALSE]

    # Compare overlapping columns value-wise
    common_cols <- intersect(colnames(obj_cd), colnames(se_cd))

    # helper equality on vectors (handles POSIXct / factor vs character)
    .vec_equal <- function(a, b) {
        if (inherits(a, "POSIXct") && inherits(b, "POSIXct")) {
            return(identical(as.numeric(a), as.numeric(b))) # ignore tz/attr differences
        }
        if (is.factor(a)) a <- as.character(a)
        if (is.factor(b)) b <- as.character(b)
        # exact value equality including NAs; avoid attributes noise
        return(isTRUE(all.equal(a, b, check.attributes = FALSE)))
    }

    bad <- vapply(common_cols, function(cc) !.vec_equal(obj_cd[[cc]], se_cd[[cc]]), logical(1))
    if (any(bad)) {
        type_in_object <- vapply(
            obj_cd[common_cols[bad]],
            function(x) paste(class(x), collapse = "/"),
            character(1)
        )
        type_in_assay <- vapply(
            se_cd[common_cols[bad]],
            function(x) paste(class(x), collapse = "/"),
            character(1)
        )
        stop(
            "Conflicting colData values in columns: ",
            paste(common_cols[bad], collapse = ", "),
            " (assay '", from_assay, "' vs incoming). ",
            "\nType in object: ", paste(type_in_object, collapse = ", "),
            "\nType in assay: ", paste(type_in_assay, collapse = ", "),
            "\nCommon colData columns overlap: ", paste(common_cols, collapse = ", "),
            "\nEnsure identical values or rename/remove conflicting columns."
        )
    }

    # Merge: keep object columns + append any NEW columns from se
    new_cols <- setdiff(colnames(se_cd), colnames(obj_cd))
    merged_cd <- if (length(new_cols)) cbind(obj_cd, se_cd[, new_cols, drop = FALSE]) else obj_cd

    colData(se) <- merged_cd
    se
}

# ---------------------------
# Constructors
# ---------------------------

#' Construct a ProBatchFeatures object from a wide matrix + sample annotation.
#' @param data_matrix numeric matrix (features x samples)
#' @param sample_annotation data.frame with sample metadata (rows = samples)
#' @param sample_id_col character(1), column in sample_annotation that matches colnames(data_matrix)
#' If missing, rownames(sample_annotation) are used.
#' @param level character label like "peptide"/"protein" (default "feature").
#' @param name character(1), optional; if missing, name is "<level>::raw".
#' If only a single value is provided to the function, without specifying whether it is a name or level, it will be used as the name value.
#' @return A `ProBatchFeatures` object.
#'
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
ProBatchFeatures <- function(
  data_matrix,
  sample_annotation = NULL,
  sample_id_col = "FullRunName",
  name = NULL,
  level = "feature"
  # TODO: add feature chain_sep - currently always "on" - and support it everywhere
) {
    stopifnot(is.matrix(data_matrix) || is.data.frame(data_matrix))
    data_matrix <- as.matrix(data_matrix)
    if (is.null(colnames(data_matrix))) {
        stop("data_matrix must have column names (sample IDs).")
    }

    # colData alignment
    if (!is.null(sample_annotation)) {
        sa <- as.data.frame(sample_annotation, stringsAsFactors = FALSE)
        if (!is.null(sample_id_col) && nzchar(sample_id_col)) {
            if (!sample_id_col %in% colnames(sa)) {
                stop("sample_id_col '", sample_id_col, "' not found in sample_annotation.")
            }
            rn <- as.character(sa[[sample_id_col]])
            if (anyNA(rn) || any(rn == "")) {
                stop("sample_id_col contains NA/empty values.")
            }
            rownames(sa) <- make.unique(rn)
        } else if (is.null(rownames(sa))) {
            stop("Provide rownames(sample_annotation) or a valid sample_id_col.")
        }
        if (!all(colnames(data_matrix) %in% rownames(sa))) {
            miss <- setdiff(colnames(data_matrix), rownames(sa))
            stop("Sample annotation missing for: ", paste(miss, collapse = ", "))
        }
        sa <- sa[colnames(data_matrix), , drop = FALSE]
        cd <- DataFrame(sa)
    } else {
        message("No sample_annotation provided; using empty colData.")
        cd <- DataFrame(row.names = colnames(data_matrix))
    }
    message("Sample annotation has ", ncol(cd), " columns and ", nrow(cd), " samples.")

    # Create SummarizedExperiment
    se <- SummarizedExperiment(
        assays  = list(intensity = data_matrix),
        colData = cd
    )
    # normalize name and always ensure "<level>::<pipeline>"
    if (is.null(name)) {
        name <- .pb_assay_name(level, "raw")
    } else if (!grepl("::", name, fixed = TRUE)) {
        name <- .pb_assay_name(level, name)
    }
    # Use QFeatures constructor to make a QFeatures, then wrap as ProBatchFeatures
    qf <- QFeatures(setNames(list(se), name),
        colData = cd
    )

    # start with empty structured oplog
    empty_log <- DataFrame(
        step      = character(),
        fun       = character(),
        from      = character(),
        to        = character(),
        params    = I(vector("list", 0L)),
        timestamp = as.POSIXct(character()),
        pkg       = character()
    )
    new(
        "ProBatchFeatures",
        qf,
        chain = character(),
        oplog = empty_log
    )
}

#' Construct from LONG df via proBatch::long_to_matrix
#' @param df_long Data frame in long format with feature/sample/value columns.
#' @param sample_annotation Optional sample metadata aligned to the samples.
#' @param feature_id_col Column containing feature identifiers in `df_long`.
#' @param sample_id_col Column containing sample identifiers in `df_long`.
#' @param measure_col Column with the measured intensity values.
#' @param level Character label describing the biological level of the assay.
#' @param name Optional pipeline name; defaults to `<level>::raw` when missing.
#' @return A `ProBatchFeatures` object constructed from the long-format input.
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
ProBatchFeatures_from_long <- function(
  df_long,
  sample_annotation = NULL,
  feature_id_col = "peptide_group_label",
  sample_id_col = "FullRunName",
  measure_col = "Intensity",
  level = "feature",
  name = NULL
) {
    stopifnot(is.data.frame(df_long))
    # 1) long -> wide using existing proBatch utility
    data_matrix <- long_to_matrix(
        df_long,
        feature_id_col = feature_id_col,
        sample_id_col  = sample_id_col,
        measure_col    = measure_col
    )
    # normalize name here too (same logic as above)
    if (is.null(name)) {
        name <- .pb_assay_name(level, "raw")
    } else if (!grepl("::", name, fixed = TRUE)) {
        name <- .pb_assay_name(level, name)
    }
    # 2) delegate to the wide constructor
    ProBatchFeatures(
        data_matrix       = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col     = sample_id_col,
        level             = level,
        name              = name
    )
}

#' Coerce a QFeatures object into ProBatchFeatures
#'
#' Wraps an existing \code{QFeatures} instance into the \code{ProBatchFeatures} subclass,
#' initializing the operation log and optional assay renaming when a single assay is present.
#'
#' @param object A \code{QFeatures} object to wrap.
#' @param level Character scalar used as the default level when renaming a single assay.
#' @param pipeline Character scalar used as the default pipeline when renaming a single assay.
#' @param sample_id_name Optional character scalar indicating the sample ID column name in colData.
#'
#' @return A \code{ProBatchFeatures} object containing the same assays as \code{object}.
#'
#' @examples
#' if (requireNamespace("QFeatures", quietly = TRUE)) {
#'     data(example_proteome_matrix, package = "proBatch")
#'     data(example_sample_annotation, package = "proBatch")
#'     cd <- S4Vectors::DataFrame(example_sample_annotation)
#'     rownames(cd) <- example_sample_annotation$FullRunName
#'     se <- SummarizedExperiment::SummarizedExperiment(
#'         assays = list(intensity = example_proteome_matrix),
#'         colData = cd
#'     )
#'     qf <- QFeatures::QFeatures(
#'         experiments = list(peptideRaw = se),
#'         colData = cd
#'     )
#'     as_ProBatchFeatures(qf, level = "peptide")
#' }
#'
#' @export
as_ProBatchFeatures <- function(object,
                                level = "feature",
                                pipeline = "raw",
                                sample_id_name = NULL) {
    if (!is(object, "QFeatures")) {
        stop("`object` must be a QFeatures object.", call. = FALSE)
    }

    qf <- object
    if (length(qf) == 0L) {
        stop("Cannot coerce a QFeatures object with no assays.", call. = FALSE)
    }

    nm <- names(qf)

    # Normalize 'level' and 'pipeline' (avoid %||%)
    level <- if (is.null(level) || is.na(level) || !nzchar(level)) "feature" else as.character(level)
    pipeline <- if (is.null(pipeline) || is.na(pipeline) || !nzchar(pipeline)) "raw" else as.character(pipeline)

    if (length(qf) == 1L) {
        current_name <- if (length(nm)) nm[[1L]] else ""
        if (!nzchar(current_name) || !grepl("::", current_name, fixed = TRUE)) {
            names(qf) <- .pb_assay_name(level, pipeline)
        }
    } else if (length(nm) && any(!grepl("::", nm, fixed = TRUE))) {
        warning("Some assay names do not follow the '<level>::<pipeline>' convention; consider renaming manually.")
    }

    empty_log <- DataFrame(
        step      = character(),
        fun       = character(),
        from      = character(),
        to        = character(),
        params    = I(vector("list", 0L)),
        timestamp = as.POSIXct(character(), tz = "UTC"),
        pkg       = character()
    )

    cd <- DataFrame(SummarizedExperiment::colData(qf))

    # Ensure colData rownames exist and are consistent with first assay
    if ((is.null(rownames(cd)) || anyNA(rownames(cd)) || any(!nzchar(rownames(cd)))) && length(qf)) {
        first_se <- qf[[1L]]
        if (is(first_se, "SummarizedExperiment")) {
            cn <- colnames(first_se)
            if (length(cn) != nrow(cd)) {
                stop("colData had invalid/missing rownames and its number of rows (", nrow(cd),
                    ") does not match the number of samples in the first assay (", length(cn), ").",
                    call. = FALSE
                )
            }
            rownames(cd) <- cn
        }
    }

    # Sample ID column handling
    final_sample_id <- NULL
    if (!is.null(sample_id_name) && nzchar(sample_id_name)) {
        final_sample_id <- as.character(sample_id_name)
        if (!final_sample_id %in% colnames(cd)) {
            warning("sample_id_name '", final_sample_id, "' not found in colData; initializing with rownames.")
            if (is.null(rownames(cd))) {
                stop("Cannot initialise sample_id_name without colData rownames.", call. = FALSE)
            }
            cd[[final_sample_id]] <- rownames(cd)
        }
    } else if ("sample_id" %in% colnames(cd)) {
        final_sample_id <- "sample_id"
    } else if (!is.null(rownames(cd))) {
        final_sample_id <- "sample_id"
        message("sample_id_name not provided; creating 'sample_id' column from colData rownames.")
        cd[[final_sample_id]] <- rownames(cd)
    }

    if (!is.null(final_sample_id)) {
        cd[[final_sample_id]] <- as.character(cd[[final_sample_id]])
    }

    colData(qf) <- cd

    nm <- names(qf)
    for (idx in seq_along(qf)) {
        se <- qf[[idx]]
        if (!is(se, "SummarizedExperiment")) {
            next
        }

        assay_name <- if (length(nm) >= idx) nm[[idx]] else paste0("assay", idx)
        se <- .pb_harmonize_colData(qf, se, from_assay = assay_name)

        if (length(assays(se)) == 1L) {
            cur <- assayNames(se)
            if (!length(cur) || is.na(cur) || !nzchar(cur)) {
                assayNames(se) <- "intensity"
            }
        }
        qf[[idx]] <- se
    }

    out <- S4Vectors::new2("ProBatchFeatures", qf, chain = character(), oplog = empty_log, check = TRUE)
    validObject(out)
    out
}


`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------------------
# Accessors
# ---------------------------

#' Access the operation log (structured)
#' @param object A `ProBatchFeatures` object.
#' @return S4Vectors::DataFrame
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
get_operation_log <- function(object) {
    stopifnot(is(object, "ProBatchFeatures"))
    object@oplog
}

#' Retrieve operation chain as vector or single string "combat_on_mediannorm_on_log"
#' @param object A `ProBatchFeatures` object.
#' @param as_string logical(1). if `TRUE` returns the chain as a single string
#'   of the form `"combat_on_mediannorm_on_log"`.
#' @return Character vector or string describing the processing chain.
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
get_chain <- function(object, as_string = FALSE) {
    stopifnot(is(object, "ProBatchFeatures"))
    ch <- object@chain
    if (!as_string) {
        return(ch)
    }
    if (!length(ch)) {
        return("")
    }
    paste(rev(ch), collapse = "_on_")
}

#' Pretty pipeline name derived from the assay
#' @param object ProBatchFeatures
#' @param assay character(1) assay name; defaults to current assay
#' @return character(1) pipeline string like "combat_on_medianNorm_on_log2" or "raw"
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
pb_pipeline_name <- function(object, assay = pb_current_assay(object)) {
    stopifnot(is(object, "ProBatchFeatures"))
    nm <- if (length(assay) == 1) assay else pb_current_assay(object)
    # Prefer the assay's own pipeline encoded as "<level>::<pipeline>"
    if (nm %in% names(object)) {
        parts <- strsplit(nm, "::", fixed = TRUE)[[1]]
        if (length(parts) >= 2 && nzchar(parts[2])) {
            return(parts[2])
        }
    }
    # Fallback for legacy objects: use global chain (may be ambiguous with branches)
    ch <- get_chain(object, as_string = TRUE)
    if (nzchar(ch)) ch else "raw"
}

#' Current (latest) assay name
#' @param object A `ProBatchFeatures` object.
#' @return character(1) assay identifier for the most recently stored assay
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
pb_current_assay <- function(object) {
    stopifnot(is(object, "ProBatchFeatures"))
    tail(names(object), 1)
}

.pb_apply_logged_step <- function(base_matrix, step, fun_name, params) {
    step <- as.character(step)
    fun_name <- as.character(fun_name)
    params <- params %||% list()

    if (step %in% c("log", "log2")) {
        log_base <- if (identical(step, "log2")) {
            params$log_base %||% 2
        } else {
            params$log_base %||% params$base %||% exp(1)
        }
        offset <- params$offset %||% params$pseudo %||% 1
        return(log_transform_dm.default(
            x = base_matrix,
            log_base = log_base,
            offset = offset
        ))
    }
    if (identical(step, "unlog")) {
        log_base <- params$log_base %||% params$base %||% 2
        offset <- params$offset %||% 1
        return(unlog_dm.default(
            x = base_matrix,
            log_base = log_base,
            offset = offset
        ))
    }

    fun_candidate <- NULL
    if (length(fun_name) && !is.na(fun_name) && nzchar(fun_name)) {
        fun_candidate <- tryCatch(.pb_get_step_fun(fun_name), error = function(e) NULL)
    }
    if (is.null(fun_candidate) && length(step) && !is.na(step) && nzchar(step)) {
        fun_candidate <- tryCatch(.pb_get_step_fun(step), error = function(e) NULL)
    }
    if (is.null(fun_candidate)) {
        stop(
            "Unable to reconstruct fast step '", step,
            "' (function '", fun_name, "' not found)."
        )
    }
    do.call(fun_candidate, c(list(base_matrix), params))
}

.pb_resolve_assay_from_log <- function(object, assay, name = "intensity", visited = character()) {
    if (assay %in% names(object)) {
        se <- object[[assay]]
        return(list(
            matrix = assay(se, i = name),
            colData = colData(se)
        ))
    }

    if (assay %in% visited) {
        stop("Detected cyclic dependency while resolving assay '", assay, "'.")
    }

    log <- get_operation_log(object)
    if (!nrow(log)) {
        return(NULL)
    }

    matches <- which(as.character(log$to) == assay)
    if (!length(matches)) {
        return(NULL)
    }

    idx <- matches[length(matches)]
    entry <- log[idx, , drop = FALSE]
    from_assay <- as.character(entry$from[[1]])
    params <- entry$params[[1]] %||% list()
    step <- entry$step[[1]]
    fun_name <- entry$fun[[1]]

    base <- .pb_resolve_assay_from_log(object, from_assay, name, c(visited, assay))
    if (is.null(base)) {
        stop("Unable to resolve base assay '", from_assay, "' for '", assay, "'.")
    }

    matrix <- .pb_apply_logged_step(base$matrix, step, fun_name, params)
    list(matrix = matrix, colData = base$colData)
}

.pb_assay_payload <- function(object, assay_name, name = "intensity") {
    if (assay_name %in% names(object)) {
        se <- object[[assay_name]]
        return(list(
            matrix = assay(se, i = name),
            colData = colData(se),
            stored = TRUE
        ))
    }

    resolved <- .pb_resolve_assay_from_log(object, assay_name, name = name)
    if (is.null(resolved)) {
        return(NULL)
    }

    list(
        matrix = resolved$matrix,
        colData = resolved$colData,
        stored = FALSE
    )
}

.pb_coldata_for_assay <- function(object, assay) {
    if (assay %in% names(object)) {
        se <- object[[assay]]
        return(colData(se))
    }

    payload <- .pb_assay_payload(object, assay_name = assay, name = "intensity")
    if (is.null(payload)) {
        stop("Unable to retrieve colData for assay '", assay, "'.")
    }
    payload$colData
}

.pb_enrich_step_params <- function(object, assay, fun, params) {
    if (is.null(params)) {
        params <- list()
    } else if (!is.list(params)) {
        params <- list(params)
    }
    fun_formals <- tryCatch(names(formals(fun)), error = function(...) character())
    needs_sa <- "sample_annotation" %in% fun_formals
    has_sa <- !is.null(params) && length(params) &&
        "sample_annotation" %in% names(params) &&
        !is.null(params[["sample_annotation"]])

    if (needs_sa && !has_sa) {
        cd <- .pb_coldata_for_assay(object, assay)
        params$sample_annotation <- as.data.frame(cd)
    }

    params
}


#' Convenience accessor for assay matrix by name/index (returns the 'intensity' assay)
#' @param object A `ProBatchFeatures` object.
#' @param assay Assay identifier to extract; defaults to the current assay.
#' @param name Assay entry to read from the underlying `SummarizedExperiment`.
#' @return assay data matrix with features in rows and samples in columns
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
pb_assay_matrix <- function(object, assay = NULL, name = "intensity") {
    if (is.null(assay)) {
        assay <- pb_current_assay(object)
        message("`assay` not provided, using the most recent assay: ", assay)
    } else {
        message("Using assay: ", assay)
    }
    payload <- .pb_assay_payload(object, assay_name = assay, name = name)
    if (is.null(payload)) {
        stop("Assay '", assay, "' not found in object or operation log.")
    }
    if (!isTRUE(payload$stored)) {
        message("Assay '", assay, "' not stored; computed from operation log.")
    }
    payload$matrix
}

#' Get current assay as LONG (via proBatch::matrix_to_long)
#' @param object A `ProBatchFeatures` object.
#' @param feature_id_col Column name used for feature identifiers in the long table.
#' @param sample_id_col Column name used for sample identifiers in the long table.
#' @param measure_col Column name containing measured values in the long table.
#' @param pbf_name Assay name whose intensities should be returned in long form.
#' @return tibble/data.frame containing one row per feature-sample combination
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
pb_as_long <- function(
  object,
  feature_id_col = "feature_label",
  sample_id_col = "FullRunName",
  measure_col = "Intensity",
  pbf_name = pb_current_assay(object)
) {
    payload <- .pb_assay_payload(object, assay_name = pbf_name, name = "intensity")
    if (is.null(payload)) {
        stop("Assay '", pbf_name, "' not found in object or operation log.")
    }

    if (isTRUE(payload$stored)) {
        message("Using stored assay '", pbf_name, "'.")
    }

    matrix_to_long(
        data_matrix       = payload$matrix,
        sample_annotation = as.data.frame(payload$colData),
        feature_id_col    = feature_id_col,
        measure_col       = measure_col,
        sample_id_col     = sample_id_col
    )
}

#' Get an assay matrix (wide)
#' @param object A `ProBatchFeatures` object.
#' @param assay Assay identifier to extract; defaults to the current assay.
#' @param name Assay entry name inside the `SummarizedExperiment` to return.
#' @return numeric matrix (wide) corresponding to the requested assay
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
pb_as_wide <- function(object, assay = pb_current_assay(object), name = "intensity") {
    stopifnot(is(object, "ProBatchFeatures"))
    se <- object[[assay]]
    assay(se, i = name)
}

# ---------------------------
# Internal core: add assay + link + log (atomic)
# ---------------------------

.pb_add_log_entry <- function(object, step, fun, from, to, params, pkg = "proBatch") {
    fun_name <- if (is.character(fun)) fun else deparse(substitute(fun))

    if (nrow(object@oplog)) {
        dup <- object@oplog$step == step &
            object@oplog$fun == fun_name &
            object@oplog$from == as.character(from) &
            object@oplog$to == as.character(to) &
            vapply(object@oplog$params, function(p) identical(p, params), logical(1))
        if (any(dup)) {
            return(object)
        }
    }

    entry <- DataFrame(
        step      = step,
        fun       = fun_name,
        from      = as.character(from),
        to        = as.character(to),
        params    = I(list(params)),
        timestamp = .pb_now(),
        pkg       = pkg
    )
    object@oplog <- rbind(object@oplog, entry)
    object@chain <- c(object@chain, step)
    object
}

.pb_add_assay_with_link <- function(object, se, to, from) {
    original_names <- names(object)
    stopifnot(
        is.character(to), length(to) == 1L,
        is.character(from), length(from) == 1L
    )

    has_from <- from %in% names(object)
    if (!has_from) {
        resolved <- .pb_resolve_assay_from_log(object, from, name = "intensity")
        if (is.null(resolved)) {
            stop("Assay '", from, "' not found in object or operation log; cannot link.")
        }
    }

    # Harmonize colData (throws if incompatible)
    se <- .pb_harmonize_colData(object, se, from_assay = from)

    # Add assay
    object <- addAssay(object, se, name = to)

    # Best-effort 1:1 link only for unique feature rownames.
    # Duplicate rownames can trigger addAssayLinkOneToOne() failures.
    ok_link <- FALSE
    if (has_from) {
        r_to <- rownames(assay(object[[to]], "intensity"))
        r_from <- rownames(assay(object[[from]], "intensity"))
        can_link_1to1 <- !is.null(r_to) && !is.null(r_from) &&
            !anyDuplicated(r_to) && !anyDuplicated(r_from) &&
            setequal(r_to, r_from)
        if (can_link_1to1) {
            object <- addAssayLinkOneToOne(object, from = from, to = to)
            ok_link <- TRUE
        }
    }
    metadata(object)$linked_last <- ok_link

    if (!(from %in% original_names) && from %in% names(object)) {
        object <- object[names(object) != from]
    }

    object
}

# ---------------------------
# The small internal helper: ONE step apply
# ---------------------------

#' Apply a single step to an assay, optionally store result, always log
#' @param object ProBatchFeatures
#' @param from assay name (e.g., "peptide::raw")
#' @param step character step id (e.g., "log2", "medianNorm", "combat")
#' @param fun function or name in registry
#' @param params list of parameters for fun
#' @param store logical: store result as new assay?
#' @param new_level optional level label for the new assay (defaults to level parsed from 'from')
#' @param to_override optional assay name override for the stored result
#' @param backend "memory","hdf5","auto"
#' @param hdf5_path optional filepath for HDF5Array
#' @return list(object=updated, assay=assay_name_or_NULL, matrix=the_result_matrix)
#'
#' Note: This is an internal function; users should typically use pb_transform() or pb_eval(). Do not export.
#' @noRd
.pb_apply_step <- function(
  object, from, step, fun, params = list(),
  store = TRUE, new_level = NULL, to_override = NULL,
  backend = c("auto", "memory", "hdf5"),
  hdf5_path = NULL, .base_m = NULL,
  from_data = from
) {
    backend <- match.arg(backend)
    stopifnot(is(object, "ProBatchFeatures"))

    from_parts <- strsplit(from, "::", fixed = TRUE)[[1]]
    base_level <- if (length(from_parts) >= 1) from_parts[1] else "feature"
    new_level <- new_level %||% base_level
    from_pipeline <- if (length(from_parts) >= 2) from_parts[2] else "raw"
    prev_tokens <- if (identical(from_pipeline, "raw")) "raw" else rev(strsplit(from_pipeline, "_on_", fixed = TRUE)[[1]])
    new_pipeline <- .pb_make_pipeline_name(c(prev_tokens, step))
    to <- to_override %||% .pb_assay_name(new_level, new_pipeline)

    # Avoid duplicate identical entries
    fun_name <- if (is.character(fun)) fun else step
    if (to %in% names(object) && nrow(object@oplog)) {
        dup <- object@oplog$step == step &
            object@oplog$fun == fun_name &
            object@oplog$from == from &
            object@oplog$to == to &
            vapply(object@oplog$params, function(p) identical(p, params), logical(1))
        if (any(dup)) {
            return(list(
                object = object,
                assay = to,
                matrix = pb_assay_matrix(object, to),
                to = to
            ))
        }
    }

    base_m <- if (!is.null(.base_m)) {
        .base_m
    } else {
        suppressMessages(
            pb_assay_matrix(object, assay = from_data)
        )
    }

    f <- .pb_get_step_fun(fun)
    params <- .pb_enrich_step_params(object, from_data, f, params)
    res_m <- do.call(f, c(list(base_m), params))

    saved_assay <- NULL
    if (store) {
        mat <- .pb_materialize_matrix(res_m, backend = backend, hdf5_path = hdf5_path)
        cd_from <- .pb_coldata_for_assay(object, from_data)
        se <- SummarizedExperiment(
            assays  = list(intensity = mat),
            colData = cd_from
        )
        object <- .pb_add_assay_with_link(object, se, to = to, from = from_data)
        saved_assay <- to
    }

    object <- .pb_add_log_entry(
        object,
        step = step,
        fun = fun_name,
        from = from, to = to, params = params
    )

    list(object = object, assay = saved_assay, matrix = res_m, to = to)
}

# ---------------------------
# Sequences / Pipelines
# ---------------------------

#' Compute a pipeline and optionally store only the final result
#' @param object A `ProBatchFeatures` object.
#' @param from Assay name to start the pipeline from.
#' @param steps character vector, e.g. c("log2","medianNorm","combat")
#' @param funs optional same-length vector/list of functions/names (default: steps)
#' @param params_list list of parameter lists (same length as steps)
#' @param level Optional level label to assign to the generated assay(s).
#' @param store_fast_steps logical; if FALSE, fast steps are computed but not stored
#' @param fast_steps which steps count as fast (default: c("log","log2","medianNorm"))
#' @param store_intermediate logical; if TRUE store every step (overrides fast behavior)
#' @param final_name optional final assay name override
#' @param backend "memory","hdf5","auto"
#' @param hdf5_path Optional file path used when `backend = "hdf5"`.
#' @return ProBatchFeatures with the requested pipeline added (as log and/or assay)
#'
#' @example inst/examples/ProBatchFeatures-basic.R
#'
#' @export
pb_transform <- function(
  object, from,
  steps,
  funs = NULL,
  params_list = NULL,
  level = NULL,
  store_fast_steps = FALSE,
  fast_steps = c("log", "log2", "medianNorm"),
  store_intermediate = FALSE,
  final_name = NULL,
  backend = c("auto", "memory", "hdf5"),
  hdf5_path = NULL
) {
    backend <- match.arg(backend)
    stopifnot(is(object, "ProBatchFeatures"))
    if (is.null(funs)) funs <- steps
    if (is.null(params_list)) params_list <- replicate(length(steps), list(), simplify = FALSE)
    stopifnot(length(steps) == length(funs), length(steps) == length(params_list))

    cur_from <- from
    cur_from_data <- from
    base_m <- NULL
    last_assay <- NULL

    for (k in seq_along(steps)) {
        step <- steps[[k]]
        fun <- funs[[k]]
        par <- params_list[[k]]

        step_label <- .pb_step_label(step, par)
        is_fast <- .pb_is_fast_step(step_label, fast_steps)
        is_final <- k == length(steps)
        # Only "log"/"log2" are treated as ephemeral fast steps by default.
        is_ephemeral_fast <- is_fast && step_label %in% c("log", "log2")
        store_this <- if (store_intermediate) {
            TRUE
        } else if (is_ephemeral_fast) {
            store_fast_steps
        } else if (is_final) {
            TRUE
        } else if (is_fast) {
            store_fast_steps
        } else {
            TRUE
        }

        use_final_name <- k == length(steps) && !is.null(final_name)
        out <- .pb_apply_step(
            object = object, from = cur_from,
            step = step_label, fun = fun, params = par,
            store = store_this, new_level = if (k == length(steps)) level else level,
            to_override = if (use_final_name) final_name else NULL,
            backend = backend, hdf5_path = hdf5_path,
            .base_m = base_m,
            from_data = cur_from_data
        )
        object <- out$object
        base_m <- out$matrix
        cur_from <- out$to %||% cur_from
        if (store_this) {
            cur_from_data <- out$assay %||% cur_from_data
            last_assay <- cur_from_data
        } else {
            last_assay <- cur_from_data
        }
    }
    # Rename final assay if requested and it exists
    if (!is.null(final_name) && !is.null(last_assay) && last_assay %in% names(object) &&
        !identical(last_assay, final_name)) {
        names(object)[match(last_assay, names(object))] <- final_name
        last_assay <- final_name
        if (nrow(object@oplog)) {
            object@oplog$to[nrow(object@oplog)] <- final_name
        }
    }
    object
}

#' Evaluate a pipeline and return the matrix, without storing
#' @param object ProBatchFeatures
#' @param from assay name (e.g., "peptide::raw")
#' @param steps character vector, e.g. c("log2","medianNorm","combat")
#' @param funs optional same-length vector/list of functions/names (default: steps
#' for registry lookup)
#' @param params_list list of parameter lists (same length as steps)
#' @return numeric matrix (features x samples)
#' @example inst/examples/ProBatchFeatures-basic.R
#'
#' @export
pb_eval <- function(
  object, from,
  steps,
  funs = NULL,
  params_list = NULL
) {
    stopifnot(is(object, "ProBatchFeatures"))
    if (is.null(funs)) funs <- steps
    if (is.null(params_list)) params_list <- replicate(length(steps), list(), simplify = FALSE)
    stopifnot(length(steps) == length(funs), length(steps) == length(params_list))

    m <- pb_assay_matrix(object, from)
    for (k in seq_along(steps)) {
        f <- .pb_get_step_fun(funs[[k]])
        params <- .pb_enrich_step_params(object, from, f, params_list[[k]])
        m <- do.call(f, c(list(m), params))
    }
    m
}

# ---------------------------
# Linking across levels (use aggregateFeatures for many-to-one)
# ---------------------------

#' Aggregate features (e.g., peptide -> protein) and store as new level
#' @param object ProBatchFeatures
#' @param from assay name (e.g., "peptide::raw")
#' @param feature_var name of a column in rowData(from) holding group labels (e.g. protein IDs)
#' @param fun summarization function (e.g., matrixStats::colMedians), or name
#' @param new_level new level label (e.g., "protein")
#' @param new_pipeline optional pipeline name (default carries over from 'from')
#' @return ProBatchFeatures with an additional aggregated assay appended
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
pb_aggregate_level <- function(
  object, from,
  feature_var,
  fun = colMedians,
  new_level = "protein",
  new_pipeline = NULL
) {
    stopifnot(is(object, "ProBatchFeatures"))
    # Let QFeatures handle both aggregation and linkage book-keeping
    from_parts <- strsplit(from, "::", fixed = TRUE)[[1]]
    from_pipeline <- if (length(from_parts) >= 2) from_parts[2] else "raw"
    pipeline <- new_pipeline %||% from_pipeline
    to <- .pb_assay_name(new_level, pipeline)

    obj <- aggregateFeatures(
        object,
        i = from,
        fcol = feature_var,
        name = to,
        fun = fun
    )
    # Log the step
    obj <- .pb_add_log_entry(obj,
        step = paste0("aggregate_", feature_var),
        fun = if (is.character(fun)) fun else "aggregate",
        from = from, to = to, params = list(fcol = feature_var)
    )
    obj
}

#' Add a new level from an external matrix and link to an existing assay
#' @param object ProBatchFeatures
#' @param from assay name (e.g., "peptide::raw")
#' @param new_matrix numeric matrix (features x samples)
#' @param mapping_df data.frame with mapping from 'from' IDs to 'to' IDs
#' @param from_id column in mapping_df for 'from' IDs (e.g., "Precursor.Id")
#' @param to_id column in mapping_df for 'to' IDs (e.g., "Protein.Ids")
#' @param map_strategy how to resolve multiple to-ids per from-id:
#' "as_is" (error if not 1:1), "first" (take first), "longest" (take longest string)
#' @param link_var rowData variable name to use for linking (e.g., "ProteinID")
#' @param to_level e.g. "protein"
#' @param to_pipeline optional pipeline name (default carries over from 'from')
#' @param name optional final assay name override
#' @param backend "memory","hdf5","auto"
#' @param hdf5_path optional filepath for HDF5Array
#' @return ProBatchFeatures with new assay and link added
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
pb_add_level <- function(
  object,
  from, # e.g. "peptide::raw"
  new_matrix, # numeric matrix (features x samples)
  to_level, # e.g. "protein"
  to_pipeline = NULL, # default = carry pipeline from 'from'
  name = NULL, # override final assay name if desired
  mapping_df = NULL, # data.frame with mapping
  from_id = NULL, # column in mapping_df for 'from' IDs (e.g., "Precursor.Id")
  to_id = NULL, # column in mapping_df for 'to' IDs   (e.g., "Protein.Ids")
  map_strategy = c("as_is", "first", "longest"), # how to resolve multiple to-ids per from-id
  link_var = "ProteinID", # rowData variable name to use for linking
  backend = c("auto", "memory", "hdf5"),
  hdf5_path = NULL
) {
    stopifnot(is(object, "ProBatchFeatures"))
    backend <- match.arg(backend)
    map_strategy <- match.arg(map_strategy)

    # ----- Determine target assay name -----
    from_parts <- strsplit(from, "::", fixed = TRUE)[[1]]
    from_pipeline <- if (length(from_parts) >= 2) from_parts[2] else "raw"
    pipeline <- to_pipeline %||% from_pipeline
    to <- if (!is.null(name)) name else .pb_assay_name(to_level, pipeline)
    # Idempotency: if identical assay recorded, do nothing
    if (to %in% names(object) && from %in% names(object)) {
        message("Assay '", to, "' already exists and is linked; skipping addition.")
        return(object)
    }

    # ----- Build SE for new_matrix and align samples to 'from' assay -----
    m <- as.matrix(new_matrix)
    if (is.null(colnames(m))) stop("new_matrix must have column names (sample IDs).")
    from_cd <- colData(object[[from]])
    if (!setequal(rownames(from_cd), colnames(m))) {
        stop("Samples of new_matrix don't match samples in '", from, "'.")
    }
    from_cd <- from_cd[colnames(m), , drop = FALSE]
    se_new <- SummarizedExperiment(
        assays  = list(intensity = .pb_materialize_matrix(m, backend, hdf5_path)),
        colData = from_cd
    )

    se_new <- .pb_harmonize_colData(object, se_new, from_assay = from)

    # ----- Add assay -----
    object <- addAssay(object, se_new, name = to)

    # ----- Linking logic -----
    # Case A: identical rownames -> 1:1 link
    r_from <- rownames(assay(object[[from]], "intensity"))
    r_to <- rownames(assay(object[[to]], "intensity"))
    if (identical(r_from, r_to)) {
        object <- addAssayLinkOneToOne(object, from = from, to = to)
        object <- .pb_add_log_entry(object,
            step = sprintf("add_level(%s)_1to1", to_level),
            fun = "addAssayLinkOneToOne",
            from = from, to = to,
            params = list()
        )
        return(object)
    }

    # Case B: many-to-one (typical peptide->protein) using addAssayLink()
    # We need a variable present in rowData(from) and rowData(to) to match on.
    # For the child (protein) side we set varTo to its rownames:
    rowData(object[[to]])[[link_var]] <- rownames(object[[to]])

    # For the parent (peptide) side, set varFrom using the mapping_df
    # mapping_df[from_id] -> mapping_df[to_id]
    if (is.null(mapping_df) || is.null(from_id) || is.null(to_id)) {
        stop("Provide mapping_df, from_id and to_id to establish cross-level links.")
    }

    # in mapping_df, ensure that parents are fully covered; children need not all be referenced
    r_from_set <- unique(r_from)
    r_to_set <- unique(r_to)
    # every parent must appear at least once in mapping_df
    miss_from <- setdiff(r_from_set, as.character(mapping_df[[from_id]]))
    if (length(miss_from)) {
        stop(
            "Mapping incomplete: ", length(miss_from),
            " parent IDs from '", from, "' have no mapping. Examples: ",
            paste(head(miss_from, 10), collapse = ", ")
        )
    }
    # keep only rows that refer to existing parents/children
    mapping_df <- mapping_df[
        as.character(mapping_df[[from_id]]) %in% r_from_set &
            as.character(mapping_df[[to_id]]) %in% r_to_set, ,
        drop = FALSE
    ]

    # Build a named vector: from_key -> to_key
    # If multiple protein groups per precursor, resolve by map_strategy
    df <- mapping_df[, c(from_id, to_id)]
    names(df) <- c("from_key", "to_key")

    # Normalize to character
    df$from_key <- as.character(df$from_key)
    df$to_key <- as.character(df$to_key)

    # Group mapping rows by parent (peptide/precursor)
    # group by parent and drop duplicate targets per parent
    by_parent <- split(df$to_key, df$from_key)
    by_parent <- lapply(by_parent, unique)

    # Resolve duplicates per from_key - apply chooser per parent
    picked <- vapply(by_parent, .choose_target, character(1), r_to = r_to, map_strategy = map_strategy)

    parent_ids <- rownames(object[[from]])
    parent_keys <- unname(picked[parent_ids])

    if (anyNA(parent_keys)) {
        bad <- parent_ids[is.na(parent_keys)]
        stop(
            "Linking failed: ", length(bad),
            " parents have no exact match under map_strategy='", map_strategy, "'. Examples: ",
            paste(head(bad, 10), collapse = ", ")
        )
    }

    # Write linking variables:
    rowData(object[[to]])[[link_var]] <- r_to
    rowData(object[[from]])[[link_var]] <- parent_keys

    # Add the link by variable
    object <- addAssayLink(object,
        from = from, to = to,
        varFrom = link_var, varTo = link_var
    )

    # ----- Log -----
    # Log
    object <- .pb_add_log_entry(object,
        step = sprintf("add_level(%s)_byVar", to_level),
        fun = "addAssayLink",
        from = from, to = to,
        params = list(
            varFrom = link_var, varTo = link_var,
            map_strategy = map_strategy
        )
    )
    object
}

# Choose for each parent, following your 3-case policy
.choose_target <- function(cands_raw, r_to, map_strategy) {
    # Helper: exact presence among child rows
    .has_child <- function(x) x %in% r_to

    # Keep only exact candidates that exist as child rows, preserving original order
    exact <- unique(cands_raw[.has_child(cands_raw)])

    if (map_strategy == "as_is") {
        # exactly one, and it must exist as a child row
        uniq <- unique(cands_raw)
        if (length(uniq) != 1L || !.has_child(uniq)) {
            stop("map_strategy='as_is' but multiple or zero targets for a parent feature.")
        }
        return(uniq[[1]])
    }

    if (map_strategy == "first") {
        # strictly first exact candidate; no splitting fallback
        if (length(exact)) {
            return(exact[1])
        }
        return(NA_character_)
    }

    if (map_strategy == "longest") {
        # choose among exact groups only:
        # 1) max number of proteins (# of ';' + 1)
        # 2) tie-break by longer string
        # 3) tie-break by first occurrence in original order
        if (length(exact)) {
            nprot <- vapply(exact, function(s) {
                length(strsplit(s, ";", fixed = TRUE)[[1]])
            }, integer(1))
            ord <- order(-nprot, -nchar(exact), seq_along(exact))
            return(exact[ord][1])
        }
        return(NA_character_)
    }

    stop("Unknown map_strategy: ", map_strategy)
}

# ---------------------------
# Preserve subclass on subsetting
# ---------------------------

# Re-wrap QFeatures results into ProBatchFeatures so slots survive
.as_ProBatchFeatures <- function(out, from) {
    if (is(out, "ProBatchFeatures")) {
        return(out)
    }
    new("ProBatchFeatures", out, chain = from@chain, oplog = from@oplog)
}

#' Subset `ProBatchFeatures` objects without dropping metadata.
#'
#' Ensures the `[` method returns a `ProBatchFeatures` instance so the
#' subclass-specific slots remain available after subsetting.
#'
#' @param x A `ProBatchFeatures` object.
#' @param i Row indices passed to the underlying `QFeatures` subset.
#' @param j Column indices passed to the underlying `QFeatures` subset.
#' @param ... Additional arguments forwarded to the next method.
#' @param drop Logical flag controlling dimension dropping; defaults to `TRUE`.
#'
#' @return A `ProBatchFeatures` object containing the requested subset.
#' @rdname ProBatchFeatures-subset
#' @aliases ProBatchFeatures-subset [,ProBatchFeatures,ANY,ANY,ANY-method
#' @export
setMethod(
    "[",
    c("ProBatchFeatures", "ANY", "ANY", "ANY"),
    function(x, i, j, k, ..., drop = TRUE) {
        out <- callNextMethod(x = x, i = i, j = j, k = k, ..., drop = drop)
        .as_ProBatchFeatures(out, from = x)
    }
)

#' Subset samples in a `ProBatchFeatures` object using sample metadata.
#'
#' This is a thin wrapper around the existing `[` method. It resolves sample
#' metadata from `colData(object)`, builds a column mask, and delegates the
#' actual subsetting to the subclass-preserving method above.
#'
#' @param object A `ProBatchFeatures` object.
#' @param sample_id_col Character scalar naming the sample identifier column in
#'   `colData(object)`. When missing from `colData(object)`, rownames are used.
#' @param subset_by Character scalar naming the metadata column used for
#'   filtering. Defaults to `sample_id_col`.
#' @param subset_values Vector of values to keep in `subset_by`.
#'
#' @return A `ProBatchFeatures` object restricted to the selected samples.
#' @examples
#' data(example_proteome_matrix, package = "proBatch")
#' data(example_sample_annotation, package = "proBatch")
#' sample_annotation <- example_sample_annotation
#' sample_annotation$Group <- rep(c("Pool", "Study"), length.out = nrow(sample_annotation))
#' pbf <- ProBatchFeatures(
#'     data_matrix = example_proteome_matrix,
#'     sample_annotation = sample_annotation,
#'     sample_id_col = "FullRunName"
#' )
#' pb_subset_samples(
#'     pbf,
#'     sample_id_col = "FullRunName",
#'     subset_by = "Group",
#'     subset_values = "Pool"
#' )
#' @export
pb_subset_samples <- function(object,
                              sample_id_col = "FullRunName",
                              subset_by = sample_id_col,
                              subset_values) {
    stopifnot(is(object, "ProBatchFeatures"))

    if (!is.character(sample_id_col) || length(sample_id_col) != 1L ||
        is.na(sample_id_col) || !nzchar(sample_id_col)) {
        stop("`sample_id_col` must be a non-empty character scalar.", call. = FALSE)
    }
    if (!is.character(subset_by) || length(subset_by) != 1L ||
        is.na(subset_by) || !nzchar(subset_by)) {
        stop("`subset_by` must be a non-empty character scalar.", call. = FALSE)
    }
    if (missing(subset_values)) {
        stop("`subset_values` must be provided.", call. = FALSE)
    }

    sample_annotation <- .pb_default_sample_annotation(
        object = object,
        sample_id_col = sample_id_col
    )

    if (!subset_by %in% colnames(sample_annotation)) {
        stop("Column '", subset_by, "' not found in `colData(object)`.", call. = FALSE)
    }

    keep <- sample_annotation[[subset_by]] %in% subset_values
    object[, keep, , drop = FALSE]
}

# ---------------------------
# Show
# ---------------------------
setMethod("show", "ProBatchFeatures", function(object) {
    callNextMethod()
    log <- get_operation_log(object)
    if (!length(object@chain)) {
        cat("  Processing chain: unprocessed data (raw) \n")
    } else {
        cat("  Processing chain:\n")
        ch_lines <- paste(utils::capture.output(print(noquote(object@chain))), collapse = "; ")
        cat("  ", ch_lines, "\n", sep = "")
    }
    if (nrow(log)) {
        # for each level, show recorded chains from operation log, from "raw" to latest
        # all operations stored in oplog should be printed - grouped by level (peptide/protein/...)
        .split_level_pipe <- function(x) {
            x <- as.character(x)
            parts <- strsplit(x, "::", fixed = TRUE)
            lvl <- vapply(parts, function(p) if (length(p) >= 2) p[1] else NA_character_, character(1))
            pipe <- vapply(parts, function(p) if (length(p) >= 2 && nzchar(p[2])) p[2] else p[1], character(1))
            list(level = lvl, pipeline = pipe)
        }
        from_lp <- .split_level_pipe(log$from)
        to_lp <- .split_level_pipe(log$to)
        levels <- unique(stats::na.omit(c(from_lp$level, to_lp$level)))

        n_tokens <- function(s) length(strsplit(s, "_on_", fixed = TRUE)[[1]])

        for (lvl in levels) {
            pipes_from <- from_lp$pipeline[from_lp$level == lvl]
            pipes_to <- to_lp$pipeline[to_lp$level == lvl]
            chains <- unique(c(pipes_from, pipes_to))
            chains <- chains[chains != "raw" & !is.na(chains) & nzchar(chains)]
            chains <- chains[order(vapply(chains, n_tokens, integer(1)), chains)]
            if (length(chains)) {
                cat("   - ", lvl, ": ", paste(chains, collapse = ", "), "\n", sep = "")
            }
        }
        cat("  Steps logged: ", nrow(log), " (see get_operation_log())\n", sep = "")
    }
})


# ---------------------------
# List available steps
# ---------------------------
#' List available pb_transform steps
#'
#' Returns the names currently registered in the internal step registry used by
#' \code{pb_transform()} and \code{pb_eval()} (e.g., "log2", "medianNorm",
#' "combat", "limmaRBE", and "BERT" if the BERT package is installed).
#'
#' @param pattern optional regular expression to filter step names.
#' @param details logical(1); if TRUE, return a S4Vectors::DataFrame with
#'   columns: step, pkg (best-effort), env, n_formals.
#' @return character vector of step names (default) or a DataFrame if
#'   \code{details=TRUE}.
#' @examples
#' pb_list_steps()
#' pb_list_steps(details = TRUE)
#' pb_list_steps("^B") # list steps starting with 'B'
#' @export
pb_list_steps <- function(pattern = NULL, details = FALSE) {
    nm <- sort(ls(envir = .pb_step_registry, all.names = FALSE))
    if (!is.null(pattern)) nm <- nm[grepl(pattern, nm)]
    if (!isTRUE(details)) {
        return(nm)
    }

    # Build a lightweight info table about each registered function
    get_pkg <- function(f) {
        en <- environment(f)
        nm <- environmentName(en)
        if (is.null(nm)) {
            return(NA_character_)
        }
        sub("^namespace:", "", nm)
    }
    funs <- mget(nm, envir = .pb_step_registry, inherits = FALSE)
    pkgs <- vapply(funs, get_pkg, character(1))
    envs <- vapply(funs, function(f) environmentName(environment(f)), character(1))
    DataFrame(
        step = nm,
        pkg = pkgs,
        env = envs,
        n_formals = vapply(funs, function(f) length(formals(f)), integer(1)),
        row.names = NULL
    )
}

#' Is a step available in the registry?
#' @param name character(1) step name (e.g., "combat", "BERT")
#' @return logical(1)
#' @examples
#' pb_has_step("medianNorm")
#' pb_has_step("BERT")
#' @export
pb_has_step <- function(name) {
    is.character(name) && length(name) == 1L &&
        exists(name, envir = .pb_step_registry, inherits = FALSE)
}
