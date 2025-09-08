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
#' @import QFeatures
#' @import SummarizedExperiment
#' @import methods
#' @importClassesFrom S4Vectors DataFrame
#' @exportClass ProBatchFeatures
methods::setClass(
    "ProBatchFeatures",
    contains = "QFeatures",
    slots = list(
        chain = "character", # ordered processing steps (last is most recent)
        oplog = "DataFrame" # structured operation log
    )
)

methods::setValidity("ProBatchFeatures", function(object) {
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
    # keep chain length consistent with oplog rows (optional but helps catch mistakes)
    if (length(object@chain) != nrow(object@oplog)) {
        return("chain length must match number of oplog rows")
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
    reg$medianNorm <- function(m) {
        # center each sample by its median
        med <- apply(m, 2, stats::median, na.rm = TRUE)
        sweep(m, 2, med, FUN = "-")
    }
    reg
})

# allow to register/override steps at runtime (e.g., map "combat" -> proBatch::combat_dm)
# Example: pb_register_step("medianNorm", proBatch::median_normalize_dm)
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
#' @export
ProBatchFeatures <- function(
    data_matrix,
    sample_annotation = NULL,
    sample_id_col = "FullRunName",
    name = NULL,
    level = "feature") {
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
        cd <- S4Vectors::DataFrame(sa)
    } else {
        cd <- S4Vectors::DataFrame(row.names = colnames(data_matrix))
    }

    # Create SummarizedExperiment
    se <- SummarizedExperiment::SummarizedExperiment(
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
    qf <- QFeatures::QFeatures(stats::setNames(list(se), name))

    # start with empty structured oplog
    empty_log <- S4Vectors::DataFrame(
        step      = character(),
        fun       = character(),
        from      = character(),
        to        = character(),
        params    = I(vector("list", 0L)),
        timestamp = as.POSIXct(character()),
        pkg       = character()
    )
    methods::new(
        "ProBatchFeatures",
        qf,
        chain = character(),
        oplog = empty_log
    )
}

#' Construct from LONG df via proBatch::long_to_matrix
#' @export
ProBatchFeatures_from_long <- function(
    df_long,
    sample_annotation = NULL,
    feature_id_col = "peptide_group_label",
    sample_id_col = "FullRunName",
    measure_col = "Intensity",
    level = "feature",
    name = NULL) {
    stopifnot(is.data.frame(df_long))
    # 1) long -> wide using existing proBatch utility
    data_matrix <- proBatch::long_to_matrix(
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

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------------------
# Accessors
# ---------------------------

#' Access the operation log (structured)
#' @param object A `ProBatchFeatures` object.
#' @return S4Vectors::DataFrame
#' @export
get_operation_log <- function(object) {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    object@oplog
}

#' Retrieve operation chain as vector or single string "combat_on_mediannorm_on_log"
#' @param object A `ProBatchFeatures` object.
#' @param as_string logical(1). if `TRUE` returns the chain as a single string
#'   of the form `"combat_on_mediannorm_on_log"`.
#' @return Character vector or string describing the processing chain.
#' @export
get_chain <- function(object, as_string = FALSE) {
    stopifnot(methods::is(object, "ProBatchFeatures"))
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
#' @export
pb_pipeline_name <- function(object, assay = pb_current_assay(object)) {
    stopifnot(methods::is(object, "ProBatchFeatures"))
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
#' @export
pb_current_assay <- function(object) {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    utils::tail(names(object), 1)
}

#' Convenience accessor for assay matrix by name/index (returns the 'intensity' assay)
#' @export
pb_assay_matrix <- function(object, assay = pb_current_assay(object), name = "intensity") {
    se <- object[[assay]]
    SummarizedExperiment::assay(se, i = name)
}

#' Get current assay as LONG (via proBatch::matrix_to_long)
#' @export
pb_as_long <- function(
    object,
    feature_id_col = "peptide_group_label",
    sample_id_col = "FullRunName",
    measure_col = "Intensity") {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    assay <- pb_current_assay(object)
    se <- object[[assay]]
    m <- SummarizedExperiment::assay(se, i = "intensity")
    sa <- as.data.frame(SummarizedExperiment::colData(se))

    proBatch::matrix_to_long(
        data_matrix       = m,
        sample_annotation = sa,
        feature_id_col    = feature_id_col,
        measure_col       = measure_col,
        sample_id_col     = sample_id_col
    )
}

#' Get an assay matrix (wide)
#' @export
pb_as_wide <- function(object, assay = pb_current_assay(object), name = "intensity") {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    se <- object[[assay]]
    SummarizedExperiment::assay(se, i = name)
}

# ---------------------------
# Internal core: add assay + link + log (atomic)
# ---------------------------

.pb_add_log_entry <- function(object, step, fun, from, to, params, pkg = "proBatch") {
    entry <- S4Vectors::DataFrame(
        step      = step,
        fun       = if (is.character(fun)) fun else deparse(substitute(fun)),
        from      = as.character(from),
        to        = as.character(to),
        params    = I(list(params)),
        timestamp = .pb_now(),
        pkg       = pkg
    )
    object@oplog <- S4Vectors::rbind(object@oplog, entry)
    object@chain <- c(object@chain, step)
    object
}

.pb_add_assay_with_link <- function(object, se, to = to, from = from) {
    stopifnot(
        is.character(to), length(to) == 1L,
        is.character(from), length(from) == 1L
    )

    # Check if from exists
    if (!from %in% names(object)) {
        stop("Assay '", from, "' not found in object. Was level/pipeline name misspecified?")
    }

    ## 1) Make sure we have sample IDs
    m <- SummarizedExperiment::assay(se, "intensity")
    if (is.null(colnames(m))) stop("New assay matrix must have column names (sample IDs).")
    samp <- colnames(m)

    ## 2) Give the new assay *empty* colData (only rownames) to avoid conflicts
    SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(row.names = samp)

    ## 3) Add the assay under the fully-qualified name `to`
    object <- QFeatures::addAssay(object, se, name = to)

    ## 4) Try 1:1 link only when rownames match exactly
    ok_link <- FALSE
    if (identical(
        rownames(SummarizedExperiment::assay(object[[to]], "intensity")),
        rownames(SummarizedExperiment::assay(object[[from]], "intensity"))
    )) {
        object <- QFeatures::addAssayLinkOneToOne(object, from = from, to = to)
        ok_link <- TRUE
    }
    S4Vectors::metadata(object)$linked_last <- ok_link
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
#' @param backend "memory","hdf5","auto"
#' @param hdf5_path optional filepath for HDF5Array
#' @return list(object=updated, assay=assay_name_or_NULL, matrix=the_result_matrix)
.pb_apply_step <- function(
    object, from, step, fun, params = list(),
    store = TRUE, new_level = NULL,
    backend = c("auto", "memory", "hdf5"),
    hdf5_path = NULL) {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    backend <- match.arg(backend)

    # fetch base matrix
    base_m <- pb_assay_matrix(object, assay = from)

    # run step (reusing user/proBatch functions when provided)
    f <- .pb_get_step_fun(fun)
    res_m <- do.call(f, c(list(base_m), params))

    # new assay name (level::pipeline)
    from_parts <- strsplit(from, "::", fixed = TRUE)[[1]]
    base_level <- if (length(from_parts) >= 1) from_parts[1] else "feature"
    new_level <- new_level %||% base_level

    # compute new pipeline name by appending this step
    # recover historical steps that led to 'from' (from the name). If 'from' is the latest assay,
    # expect its pipeline name to be suffix after '::'
    from_pipeline <- if (length(from_parts) >= 2) from_parts[2] else "raw"
    prev_steps <- if (identical(from_pipeline, "raw")) character() else rev(strsplit(from_pipeline, "_on_")[[1]])
    new_pipeline <- .pb_make_pipeline_name(c(prev_steps, step))
    to <- .pb_assay_name(new_level, new_pipeline)

    saved_assay <- NULL
    if (store) {
        # materialize (possibly to HDF5) then store
        mat <- .pb_materialize_matrix(res_m, backend = backend, hdf5_path = hdf5_path)
        se <- SummarizedExperiment::SummarizedExperiment(
            assays = list(intensity = mat),
            colData = SummarizedExperiment::colData(object[[from]])
        )
        object <- .pb_add_assay_with_link(object, se, name = to, from = from)
        saved_assay <- to
    }

    # always log
    object <- .pb_add_log_entry(object,
        step = step, fun = if (is.character(fun)) fun else step,
        from = from, to = to, params = params
    )

    list(object = object, assay = saved_assay, matrix = res_m)
}

# ---------------------------
# Sequences / Pipelines
# ---------------------------

#' Compute a pipeline and optionally store only the final result
#' @param steps character vector, e.g. c("log2","medianNorm","combat")
#' @param funs optional same-length vector/list of functions/names (default: steps)
#' @param params_list list of parameter lists (same length as steps)
#' @param store_fast_steps logical; if FALSE, fast steps are computed but not stored
#' @param fast_steps which steps count as fast (default: c("log","log2","medianNorm"))
#' @param store_intermediate logical; if TRUE store every step (overrides fast behavior)
#' @param final_name optional final assay name override
#' @param backend "memory","hdf5","auto"
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
    hdf5_path = NULL) {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    backend <- match.arg(backend)
    if (is.null(funs)) funs <- steps
    if (is.null(params_list)) params_list <- replicate(length(steps), list(), simplify = FALSE)
    stopifnot(length(steps) == length(funs), length(steps) == length(params_list))

    cur_from <- from
    last_assay <- NULL
    last_mat <- NULL

    for (k in seq_along(steps)) {
        step <- steps[[k]]
        fun <- funs[[k]]
        par <- params_list[[k]]
        this_is_fast <- .pb_is_fast_step(step, fast_steps)

        store_this <- if (store_intermediate) TRUE else if (this_is_fast) store_fast_steps else TRUE

        out <- .pb_apply_step(
            object = object, from = cur_from,
            step = step, fun = fun, params = par,
            store = store_this, new_level = if (k == length(steps)) level else level,
            backend = backend, hdf5_path = hdf5_path
        )
        object <- out$object
        last_mat <- out$matrix
        last_assay <- out$assay %||% cur_from # if not stored, no new assay name

        # If stored, move the cursor; if not, keep 'from' but the next step should compute on last_mat.
        # To support ephemeral chaining, when not stored temporarily inject a volatile assay name
        # by reusing 'from' matrix in-memory. Will handle this by computing next step on last_mat directly.
        # Implementation detail: when an intermediate wasn't stored, set cur_from to the same original
        # but also attach the last_mat to a local variable to be used as base – realized below.
        if (!isTRUE(store_this)) {
            # monkey patch: replace base fetch in next iteration
            assign(".pb_ephemeral_base", last_mat, inherits = FALSE)
            .pb_fetch_base <- function(obj, assay) get(".pb_ephemeral_base", inherits = FALSE)
            fetch_fun <- .pb_fetch_base
        } else {
            fetch_fun <- pb_assay_matrix
            cur_from <- out$assay
        }

        # override the default base fetcher for next loop iteration
        assign(".pb_fetch_fun", fetch_fun, inherits = FALSE)

        # redefine pb_assay_matrix() binding locally for next iteration:
        unlockBinding("pb_assay_matrix", as.environment("package:proBatch"))
    }

    # Rename final assay if requested and it exists
    if (!is.null(final_name) && !is.null(last_assay) && last_assay %in% names(object)) {
        names(object)[match(last_assay, names(object))] <- final_name
        last_assay <- final_name
    }

    object
}

#' Evaluate a pipeline and return the matrix, without storing
#' @export
pb_eval <- function(
    object, from,
    steps,
    funs = NULL,
    params_list = NULL) {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    if (is.null(funs)) funs <- steps
    if (is.null(params_list)) params_list <- replicate(length(steps), list(), simplify = FALSE)
    stopifnot(length(steps) == length(funs), length(steps) == length(params_list))

    m <- pb_assay_matrix(object, from)
    for (k in seq_along(steps)) {
        f <- .pb_get_step_fun(funs[[k]])
        m <- do.call(f, c(list(m), params_list[[k]]))
    }
    m
}

# ---------------------------
# Linking across levels (use aggregateFeatures for many-to-one)
# ---------------------------

#' Aggregate features (e.g., peptide -> protein) and store as new level
#' @param feature_var name of a column in rowData(from) holding group labels (e.g. protein IDs)
#' @param fun summarization function (e.g., matrixStats::colMedians), or name
#' @param new_level new level label (e.g., "protein")
#' @param new_pipeline optional pipeline name (default carries over from 'from')
#' @export
pb_aggregate_level <- function(
    object, from,
    feature_var,
    fun = matrixStats::colMedians,
    new_level = "protein",
    new_pipeline = NULL) {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    # Let QFeatures handle both aggregation and linkage book-keeping
    from_parts <- strsplit(from, "::", fixed = TRUE)[[1]]
    from_pipeline <- if (length(from_parts) >= 2) from_parts[2] else "raw"
    pipeline <- new_pipeline %||% from_pipeline
    to <- .pb_assay_name(new_level, pipeline)

    obj <- QFeatures::aggregateFeatures(
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

# ---------------------------
# Preserve subclass on subsetting
# ---------------------------

# Re-wrap QFeatures results into ProBatchFeatures so slots survive
.as_ProBatchFeatures <- function(out, from) {
    if (methods::is(out, "ProBatchFeatures")) {
        return(out)
    }
    methods::new("ProBatchFeatures", out, chain = from@chain, oplog = from@oplog)
}

# Override '[' so subsetting doesn't drop subclass or logs
methods::setMethod(
    "[",
    signature(x = "ProBatchFeatures", i = "ANY", j = "ANY", drop = "ANY"),
    function(x, i, j, ..., drop = TRUE) {
        out <- methods::callNextMethod()
        .as_ProBatchFeatures(out, from = x)
    }
)

# ---------------------------
# Show
# ---------------------------
methods::setMethod("show", "ProBatchFeatures", function(object) {
    methods::callNextMethod()
    ch <- get_chain(object, as_string = TRUE)
    if (nzchar(ch)) {
        cat("  Processing chain: ", ch, "\n", sep = "")
    } else {
        cat("  Processing chain: unprocessed data (raw) \n")
    }
    if (nrow(object@oplog)) {
        cat("  Steps logged: ", nrow(object@oplog), " (see get_operation_log())\n", sep = "")
    }
})
