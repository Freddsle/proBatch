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
#' @return No return value; defines the \code{ProBatchFeatures} S4 class.
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

.pb_assay_name <- function(level, pipeline) {
    level <- if (is.null(level) || !nzchar(level)) "feature" else level
    pipeline <- if (is.null(pipeline) || !nzchar(pipeline)) "raw" else pipeline
    paste0(level, "::", pipeline)
}

.pb_assay_parts <- function(name, strict = FALSE) {
    if (!is.character(name) || length(name) != 1L ||
        is.na(name) || !nzchar(name)) {
        stop("Assay name must be a non-empty character scalar.", call. = FALSE)
    }
    separator_positions <- gregexpr("::", name, fixed = TRUE)[[1L]]
    separator_count <- if (separator_positions[[1L]] < 0L) {
        0L
    } else {
        length(separator_positions)
    }
    if (separator_count == 0L && !strict) {
        return(list(level = "feature", pipeline = name))
    }
    if (separator_count != 1L) {
        stop(
            "Assay name '", name,
            "' must use the '<level>::<pipeline>' form.",
            call. = FALSE
        )
    }
    parts <- strsplit(name, "::", fixed = TRUE)[[1L]]
    if (length(parts) != 2L || any(!nzchar(parts))) {
        stop(
            "Assay name '", name,
            "' must use the '<level>::<pipeline>' form.",
            call. = FALSE
        )
    }
    list(level = parts[[1L]], pipeline = parts[[2L]])
}

.pb_is_fast_step <- function(step, fast_steps) {
    isTRUE(step %in% fast_steps)
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

    # helper equality on vectors (handles POSIXct / factor vs character /
    # factor-with-integer-like-levels vs numeric)
    .vec_equal <- function(a, b) {
        if (inherits(a, "POSIXct") && inherits(b, "POSIXct")) {
            return(identical(as.numeric(a), as.numeric(b))) # ignore tz/attr differences
        }
        .factor_numeric_equal <- function(factor_value, numeric_value) {
            factor_numeric <- suppressWarnings(
                as.numeric(as.character(factor_value))
            )
            factor_missing <- is.na(factor_value)

            if (!identical(is.na(factor_numeric), factor_missing) ||
                !identical(factor_missing, is.na(numeric_value))) {
                return(FALSE)
            }

            observed <- factor_numeric[!factor_missing]
            if (any(!is.finite(observed) | observed != trunc(observed))) {
                return(FALSE)
            }

            isTRUE(all.equal(
                factor_numeric,
                as.numeric(numeric_value),
                check.attributes = FALSE
            ))
        }
        if (is.factor(a) && is.numeric(b)) {
            return(.factor_numeric_equal(a, b))
        }
        if (is.numeric(a) && is.factor(b)) {
            return(.factor_numeric_equal(b, a))
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
#' @md
#' @param data_matrix numeric matrix (features x samples)
#' @param sample_annotation Optional data frame with sample metadata (rows =
#'   samples). When `NULL`, the object uses empty sample metadata.
#' @param sample_id_col Character scalar naming the column in
#'   `sample_annotation` that matches `colnames(data_matrix)`. Set it to `NULL`
#'   or `""` to use row names instead.
#' @param level character label like "peptide"/"protein" (default "feature").
#' @param name Optional non-missing, non-empty character scalar. A pipeline
#'   shorthand is combined with `level`; a full name must use
#'   `<level>::<pipeline>`. If missing, the name is `<level>::raw`.
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
    colnames(data_matrix) <- .pb_validate_identifiers(
        colnames(data_matrix),
        "`data_matrix` sample axis"
    )
    if (!is.null(rownames(data_matrix))) {
        rownames(data_matrix) <- .pb_validate_identifiers(
            rownames(data_matrix),
            "`data_matrix` feature axis"
        )
    }

    # colData alignment
    if (!is.null(sample_annotation)) {
        sa <- as.data.frame(sample_annotation, stringsAsFactors = FALSE)
        if (!is.null(sample_id_col) && nzchar(sample_id_col)) {
            if (!sample_id_col %in% colnames(sa)) {
                stop("sample_id_col '", sample_id_col, "' not found in sample_annotation.")
            }
            rn <- .pb_validate_identifiers(
                sa[[sample_id_col]],
                "sample_id_col"
            )
            rownames(sa) <- rn
        } else if (is.null(rownames(sa))) {
            stop("Provide rownames(sample_annotation) or a valid sample_id_col.")
        } else {
            rownames(sa) <- .pb_validate_identifiers(
                rownames(sa),
                "`sample_annotation` row names"
            )
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
    } else {
        .pb_assay_parts(name)
        if (!grepl("::", name, fixed = TRUE)) {
            name <- .pb_assay_name(level, name)
        }
    }
    .pb_assay_parts(name, strict = TRUE)
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
#' @md
#' @param df_long Data frame in long format with feature/sample/value columns.
#' @param sample_annotation Optional sample metadata aligned to the samples.
#' @param feature_id_col Column containing feature identifiers in `df_long`.
#' @param sample_id_col Column containing sample identifiers in `df_long`.
#' @param measure_col Column with the measured intensity values.
#' @param level Character label describing the biological level of the assay.
#' @param name Optional non-missing, non-empty pipeline or full assay name;
#'   defaults to `<level>::raw` when missing.
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
    # 2) delegate naming and construction to the wide constructor
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
    methods::validObject(out)
    out
}


`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------------------
# Accessors
# ---------------------------

#' Access the operation log (structured)
#' @md
#' @param object A `ProBatchFeatures` object.
#' @return S4Vectors::DataFrame
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
get_operation_log <- function(object) {
    stopifnot(is(object, "ProBatchFeatures"))
    object@oplog
}

.pb_log_edge_equal <- function(log, left, right) {
    stable_fields <- c("from", "to", "step", "fun", "pkg")
    same_fields <- vapply(stable_fields, function(field) {
        identical(
            as.character(log[[field]][[left]]),
            as.character(log[[field]][[right]])
        )
    }, logical(1))
    all(same_fields) &&
        identical(log$params[[left]], log$params[[right]])
}

.pb_log_edge_matches <- function(
  log, index, from, to, step, fun, pkg, params
) {
    identical(as.character(log$from[[index]]), as.character(from)) &&
        identical(as.character(log$to[[index]]), as.character(to)) &&
        identical(as.character(log$step[[index]]), as.character(step)) &&
        identical(as.character(log$fun[[index]]), as.character(fun)) &&
        identical(as.character(log$pkg[[index]]), as.character(pkg)) &&
        identical(log$params[[index]], params)
}

.pb_unique_log_edges <- function(log, indices) {
    keep <- integer()
    for (index in indices) {
        duplicate <- length(keep) &&
            any(vapply(
                keep,
                function(previous) {
                    .pb_log_edge_equal(
                        log,
                        previous,
                        index
                    )
                },
                logical(1)
            ))
        if (!duplicate) {
            keep <- c(keep, index)
        }
    }
    keep
}

.pb_lineage_parent_edge <- function(log, assay) {
    matches <- which(
        as.character(log$to) == assay &
            as.character(log$from) != assay
    )
    if (!length(matches)) {
        return(integer())
    }

    origins <- .pb_unique_log_edges(log, matches)
    if (length(origins) > 1L) {
        stop(
            "Ambiguous operation-log lineage for result '",
            assay, "': multiple non-self origins are recorded.",
            call. = FALSE
        )
    }

    tail(matches, 1L)
}

.pb_assay_data_identical <- function(left, right) {
    identical(as.matrix(left), as.matrix(right))
}

.pb_target_retry_status <- function(
  object, to, from, data, step, fun, params = list(),
  pkg = "proBatch", name = "intensity"
) {
    stopifnot(is(object, "ProBatchFeatures"))
    log <- get_operation_log(object)
    stored <- to %in% names(object)
    target_rows <- if (nrow(log)) {
        which(as.character(log$to) == to)
    } else {
        integer()
    }

    if (!stored && !length(target_rows)) {
        return("available")
    }

    .pb_lineage_parent_edge(log, to)
    exact_rows <- target_rows[vapply(
        target_rows,
        function(index) {
            .pb_log_edge_matches(
                log = log,
                index = index,
                from = from,
                to = to,
                step = step,
                fun = fun,
                pkg = pkg,
                params = params
            )
        },
        logical(1)
    )]
    if (!length(exact_rows)) {
        stop(
            "Result target '", to,
            "' already exists or is reserved with a different parent ",
            "or stable lineage origin.",
            call. = FALSE
        )
    }

    .pb_lineage_from_log(object, to)
    payload <- .pb_assay_payload(object, assay_name = to, name = name)
    if (is.null(payload)) {
        stop(
            "Result target '", to,
            "' is reserved but its existing data cannot be resolved.",
            call. = FALSE
        )
    }
    if (!.pb_assay_data_identical(payload$matrix, data)) {
        stop(
            "Result target '", to,
            "' already exists or is reserved with conflicting data.",
            call. = FALSE
        )
    }

    if (stored) "stored_idempotent" else "virtual_idempotent"
}

.pb_lineage_from_log <- function(object, assay) {
    log <- get_operation_log(object)
    node <- as.character(assay)
    visited <- character()
    reverse_steps <- character()
    lineage_nodes <- node

    while (nrow(log)) {
        if (node %in% visited) {
            stop(
                "Cyclic operation-log lineage detected at assay '",
                node, "'.",
                call. = FALSE
            )
        }
        visited <- c(visited, node)
        matches <- which(as.character(log$to) == node)
        if (!length(matches)) {
            break
        }

        self_matches <- matches[as.character(log$from[matches]) == node]
        self_matches <- .pb_unique_log_edges(log, self_matches)
        if (length(self_matches)) {
            reverse_steps <- c(
                reverse_steps,
                rev(as.character(log$step[self_matches]))
            )
        }

        parent <- .pb_lineage_parent_edge(log, node)
        if (!length(parent)) {
            break
        }
        reverse_steps <- c(reverse_steps, as.character(log$step[[parent]]))
        node <- as.character(log$from[[parent]])
        lineage_nodes <- c(lineage_nodes, node)
    }

    list(
        steps = rev(reverse_steps),
        root = node,
        nodes = lineage_nodes
    )
}

.pb_pipeline_tokens <- function(object, assay) {
    lineage <- .pb_lineage_from_log(object, assay)
    root <- .pb_assay_parts(lineage$root)
    root_tokens <- rev(strsplit(
        root$pipeline,
        "_on_",
        fixed = TRUE
    )[[1L]])
    c(root_tokens, lineage$steps)
}

#' Retrieve operation chain as vector or single string "combat_on_mediannorm_on_log"
#' @md
#' @param object A `ProBatchFeatures` object.
#' @param as_string logical(1). if `TRUE` returns the chain as a single string
#'   of the form `"combat_on_mediannorm_on_log"`.
#' @param assay Optional assay identifier. When supplied, derive that assay's
#'   lineage from operation-log edges. When `NULL`, return the legacy global
#'   `chain` slot.
#' @return Character vector or string describing the processing chain.
#' @details Logged lineage is validated before a chain is returned. Ambiguous
#'   or cyclic lineage causes an error, including when `assay = NULL`.
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
get_chain <- function(object, as_string = FALSE, assay = NULL) {
    stopifnot(is(object, "ProBatchFeatures"))
    ch <- if (is.null(assay)) {
        if (nrow(object@oplog)) {
            logged_targets <- unique(as.character(object@oplog$to))
            for (target in logged_targets) {
                .pb_lineage_from_log(object, target)
            }
        }
        object@chain
    } else {
        assay <- as.character(assay)
        if (length(assay) != 1L || is.na(assay) || !nzchar(assay)) {
            stop("`assay` must be one non-empty assay name.", call. = FALSE)
        }
        known <- assay %in% names(object) ||
            any(as.character(object@oplog$to) == assay) ||
            any(as.character(object@oplog$from) == assay)
        if (!known) {
            stop(
                "Assay '", assay,
                "' not found in object or operation log.",
                call. = FALSE
            )
        }
        .pb_lineage_from_log(object, assay)$steps
    }
    if (!as_string) {
        return(ch)
    }
    if (!length(ch)) {
        return("")
    }
    paste(rev(ch), collapse = "_on_")
}

#' Pretty pipeline name derived from the assay
#' @md
#' @param object ProBatchFeatures
#' @param assay character(1) assay name; defaults to current assay
#' @return character(1) pipeline string like "combat_on_medianNorm_on_log2" or "raw"
#' @details The operation-log lineage must be unambiguous and acyclic.
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
pb_pipeline_name <- function(object, assay = pb_current_assay(object)) {
    stopifnot(is(object, "ProBatchFeatures"))
    nm <- if (length(assay) == 1) assay else pb_current_assay(object)
    .pb_make_pipeline_name(.pb_pipeline_tokens(object, nm))
}

#' Current (latest) assay name
#' @md
#' @param object A `ProBatchFeatures` object.
#' @return character(1) assay identifier for the most recently stored assay
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
pb_current_assay <- function(object) {
    stopifnot(is(object, "ProBatchFeatures"))
    tail(names(object), 1)
}

.pb_apply_logged_step <- function(
  base_matrix, step, fun_name, params, pkg = "proBatch",
  sample_annotation = NULL
) {
    step <- as.character(step)
    fun_name <- as.character(fun_name)
    pkg <- as.character(pkg)
    params <- params %||% list()

    if (identical(pkg, "proBatch") &&
        identical(step, "filterNA") &&
        identical(fun_name, "filterNA")) {
        replay_params <- params
        replay_params$inplace <- NULL
        return(do.call(
            .pb_filterNA_matrix,
            c(list(data_matrix = base_matrix), replay_params)
        ))
    }
    if (identical(pkg, "proBatch") &&
        identical(step, "groupfilterNA") &&
        identical(fun_name, "pb_groupfilterNA")) {
        replay_params <- params
        replay_params$inplace <- NULL
        return(do.call(
            .pb_groupfilterNA_matrix,
            c(
                list(
                    data_matrix = base_matrix,
                    sample_annotation = sample_annotation
                ),
                replay_params
            )
        ))
    }

    registered_token <- NULL
    if (length(fun_name) == 1L && !is.na(fun_name) && nzchar(fun_name) &&
        pb_has_step(fun_name)) {
        record <- .pb_lookup_step_record(fun_name)
        if (identical(record$package, pkg)) {
            registered_token <- fun_name
        }
    }

    if (!is.null(registered_token)) {
        resolution <- .pb_resolve_step(
            registered_token,
            package = pkg,
            require_available = TRUE
        )
        fun_candidate <- resolution$fun
    } else if (!identical(pkg, "proBatch")) {
        .pb_resolve_step(
            fun_name,
            package = pkg,
            require_available = TRUE
        )
    } else {
        fun_candidate <- NULL
        if (length(step) == 1L && !is.na(step) && nzchar(step) &&
            pb_has_step(step)) {
            step_record <- .pb_lookup_step_record(step)
            if (identical(step_record$package, pkg)) {
                resolution <- .pb_resolve_step(
                    step,
                    package = pkg,
                    require_available = TRUE
                )
                fun_candidate <- resolution$fun
            }
        }

        if (is.null(fun_candidate) && step %in% c("log", "log2")) {
            log_base <- if (identical(step, "log2")) {
                params$log_base %||% 2
            } else {
                params$log_base %||% params$base %||% exp(1)
            }
            offset <- params$offset %||% params$pseudo %||% 1
            raw_result <- log_transform_dm.default(
                x = base_matrix,
                log_base = log_base,
                offset = offset
            )
            return(.pb_step_result_matrix_parts(
                raw_result,
                paste0("Replayed step '", step, "'")
            )$data)
        }
        if (is.null(fun_candidate) && identical(step, "unlog")) {
            log_base <- params$log_base %||% params$base %||% 2
            offset <- params$offset %||% 1
            raw_result <- unlog_dm.default(
                x = base_matrix,
                log_base = log_base,
                offset = offset
            )
            return(.pb_step_result_matrix_parts(
                raw_result,
                paste0("Replayed step '", step, "'")
            )$data)
        }

        if (is.null(fun_candidate) &&
            length(fun_name) == 1L && !is.na(fun_name) && nzchar(fun_name)) {
            fun_candidate <- tryCatch(
                .pb_get_step_fun(fun_name, use_registry = FALSE),
                error = function(error) NULL
            )
        }
        if (is.null(fun_candidate) &&
            length(step) == 1L && !is.na(step) && nzchar(step)) {
            fun_candidate <- tryCatch(
                .pb_get_step_fun(step, use_registry = FALSE),
                error = function(error) NULL
            )
        }
        if (is.null(fun_candidate)) {
            stop(
                "Unable to reconstruct fast step '", step,
                "' (function '", fun_name, "' not found)."
            )
        }
    }

    call_params <- .pb_enrich_step_params(
        fun = fun_candidate,
        params = params,
        sample_annotation = sample_annotation,
        sample_ids = colnames(base_matrix)
    )
    raw_result <- do.call(fun_candidate, c(list(base_matrix), call_params))
    .pb_step_result_matrix_parts(
        raw_result,
        paste0("Replayed step '", step, "'")
    )$data
}

.pb_resolve_assay_from_log <- function(object, assay, name = "intensity", visited = character()) {
    if (!length(visited)) {
        .pb_lineage_from_log(object, assay)
    }

    if (assay %in% visited) {
        stop("Detected cyclic dependency while resolving assay '", assay, "'.")
    }

    if (assay %in% names(object)) {
        se <- object[[assay]]
        return(list(
            matrix = assay(se, i = name),
            colData = colData(se)
        ))
    }

    log <- get_operation_log(object)
    if (!nrow(log)) {
        return(NULL)
    }

    matches <- which(as.character(log$to) == assay)
    if (!length(matches)) {
        return(NULL)
    }

    idx <- .pb_lineage_parent_edge(log, assay)
    if (!length(idx)) {
        stop(
            "Unable to resolve virtual assay '", assay,
            "' from self-referential operation-log entries.",
            call. = FALSE
        )
    }
    entry <- log[idx, , drop = FALSE]
    from_assay <- as.character(entry$from[[1]])
    params <- entry$params[[1]] %||% list()
    step <- entry$step[[1]]
    fun_name <- entry$fun[[1]]
    pkg <- entry$pkg[[1]]

    base <- .pb_resolve_assay_from_log(object, from_assay, name, c(visited, assay))
    if (is.null(base)) {
        stop("Unable to resolve base assay '", from_assay, "' for '", assay, "'.")
    }

    matrix <- .pb_apply_logged_step(
        base$matrix,
        step,
        fun_name,
        params,
        pkg = pkg,
        sample_annotation = base$colData
    )
    self_matches <- matches[as.character(log$from[matches]) == assay]
    self_matches <- .pb_unique_log_edges(log, self_matches)
    for (self_index in self_matches) {
        self_entry <- log[self_index, , drop = FALSE]
        matrix <- .pb_apply_logged_step(
            matrix,
            self_entry$step[[1]],
            self_entry$fun[[1]],
            self_entry$params[[1]] %||% list(),
            pkg = self_entry$pkg[[1]],
            sample_annotation = base$colData
        )
    }
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

.pb_enrich_step_params <- function(
  object = NULL, assay = NULL, fun, params, sample_annotation = NULL,
  sample_ids = NULL
) {
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

    if (needs_sa || has_sa) {
        cd <- if (has_sa) params[["sample_annotation"]] else sample_annotation
        if (is.null(cd) && needs_sa) {
            if (!is.null(object) && !is.null(assay)) {
                cd <- .pb_coldata_for_assay(object, assay)
            } else {
                stop(
                    "Sample annotation is required to invoke this step.",
                    call. = FALSE
                )
            }
        }

        if (!is.null(cd) && !is.null(sample_ids)) {
            cd <- as.data.frame(cd)
            sample_id_col <- params[["sample_id_col"]]
            if (is.null(sample_id_col)) {
                formal_id_col <- tryCatch(
                    formals(fun)[["sample_id_col"]],
                    error = function(...) NULL
                )
                if (is.character(formal_id_col) &&
                    length(formal_id_col) == 1L &&
                    !is.na(formal_id_col) && nzchar(formal_id_col)) {
                    sample_id_col <- formal_id_col
                }
            }
            cd <- .align_sample_annotation(
                cd,
                sample_ids = sample_ids,
                sample_id_col = sample_id_col
            )
            rownames(cd) <- sample_ids
        }
        if (!is.null(cd)) {
            params$sample_annotation <- as.data.frame(cd)
        }
    }

    params
}

.pb_step_invocation <- function(fun_or_name, step) {
    if (is.character(fun_or_name) && length(fun_or_name) == 1L &&
        !is.na(fun_or_name) && nzchar(fun_or_name) &&
        pb_has_step(fun_or_name)) {
        return(.pb_resolve_step(fun_or_name, require_available = TRUE))
    }

    fun <- .pb_get_step_fun(fun_or_name)
    logged_name <- if (is.character(fun_or_name) &&
        length(fun_or_name) == 1L && !is.na(fun_or_name) &&
        nzchar(fun_or_name)) {
        fun_or_name
    } else {
        as.character(step)
    }
    list(
        input_name = logged_name,
        name = logged_name,
        fun = fun,
        package = "proBatch",
        available = TRUE,
        registered = FALSE
    )
}


#' Convenience accessor for stored or virtual assay matrices
#'
#' Stored assays are read directly. A virtual target recorded in the operation
#' log is replayed through its recorded canonical step and provider; replay
#' fails with guidance if that provider is no longer registered or available.
#' Ambiguous or cyclic operation-log lineage also causes an error.
#' @md
#' @param object A `ProBatchFeatures` object.
#' @param assay Assay identifier to extract; defaults to the current assay.
#' @param name Assay entry to read from the underlying `SummarizedExperiment`.
#' @return Assay data matrix with features in rows and samples in columns.
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
#' @md
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
#' @md
#' @param object A `ProBatchFeatures` object.
#' @param assay Assay identifier to extract; defaults to the current assay.
#' @param name Assay entry name inside the `SummarizedExperiment` to return.
#' @return numeric matrix (wide) corresponding to the requested assay
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
pb_as_wide <- function(object, assay = pb_current_assay(object), name = "intensity") {
    stopifnot(is(object, "ProBatchFeatures"))
    pb_assay_matrix(object, assay = assay, name = name)
}

# ---------------------------
# Internal core: add assay + link + log (atomic)
# ---------------------------

.pb_add_log_entry <- function(object, step, fun, from, to, params, pkg = "proBatch") {
    fun_name <- if (is.character(fun)) fun else deparse(substitute(fun))
    from <- as.character(from)
    to <- as.character(to)
    pkg <- as.character(pkg)
    self_edge <- identical(from, to)
    if (self_edge) {
        target_has_origin <- to %in% names(object) ||
            (nrow(object@oplog) &&
                any(
                    as.character(object@oplog$to) == to &
                        as.character(object@oplog$from) != to
                ))
        if (!target_has_origin) {
            stop(
                "A self operation-log edge requires an existing stored ",
                "or virtual result target '", to, "'.",
                call. = FALSE
            )
        }
    }

    if (nrow(object@oplog)) {
        target_rows <- which(as.character(object@oplog$to) == to)
        exact_rows <- target_rows[vapply(
            target_rows,
            function(index) {
                .pb_log_edge_matches(
                    log = object@oplog,
                    index = index,
                    from = from,
                    to = to,
                    step = step,
                    fun = fun_name,
                    pkg = pkg,
                    params = params
                )
            },
            logical(1)
        )]

        if (!self_edge) {
            parent_rows <- target_rows[
                as.character(object@oplog$from[target_rows]) != to
            ]
            if (length(parent_rows)) {
                all_same <- all(vapply(
                    parent_rows,
                    function(index) {
                        .pb_log_edge_matches(
                            log = object@oplog,
                            index = index,
                            from = from,
                            to = to,
                            step = step,
                            fun = fun_name,
                            pkg = pkg,
                            params = params
                        )
                    },
                    logical(1)
                ))
                if (all_same) {
                    return(object)
                }
                stop(
                    "Operation-log result '", to,
                    "' already has a different non-self parent or stable origin.",
                    call. = FALSE
                )
            }
        } else {
            if (length(exact_rows)) {
                return(object)
            }
        }
    }

    if (!self_edge) {
        ancestry <- .pb_lineage_from_log(object, from)$nodes
        if (to %in% ancestry) {
            stop(
                "Adding operation-log edge '", from, "' -> '", to,
                "' would create cyclic lineage.",
                call. = FALSE
            )
        }
    }

    entry <- DataFrame(
        step      = step,
        fun       = fun_name,
        from      = from,
        to        = to,
        params    = I(list(params)),
        timestamp = .pb_now(),
        pkg       = pkg
    )
    object@oplog <- rbind(object@oplog, entry)
    object@chain <- c(object@chain, step)
    object
}

.pb_add_assay_with_link <- function(
  object, se, to, from, step = NULL, fun = NULL, params = list(),
  pkg = "proBatch", lineage_from = from
) {
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

    target_exists <- to %in% names(object) ||
        (nrow(object@oplog) &&
            any(as.character(object@oplog$to) == to))
    retry_status <- "available"
    if (target_exists) {
        if (is.null(step) || is.null(fun)) {
            stop(
                "Result target '", to,
                "' already exists or is reserved; stable lineage origin ",
                "is required to validate an exact retry.",
                call. = FALSE
            )
        }
        retry_status <- .pb_target_retry_status(
            object = object,
            to = to,
            from = lineage_from,
            data = assay(se, "intensity"),
            step = step,
            fun = fun,
            params = params,
            pkg = pkg
        )
        if (identical(retry_status, "stored_idempotent")) {
            return(object)
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

    invocation <- .pb_step_invocation(fun, step = step)
    fun_name <- invocation$name
    fun_package <- invocation$package
    base_m <- if (!is.null(.base_m)) {
        .base_m
    } else {
        suppressMessages(
            pb_assay_matrix(object, assay = from_data)
        )
    }

    logged_params <- if (is.null(params)) {
        list()
    } else if (is.list(params)) {
        params
    } else {
        list(params)
    }
    call_params <- .pb_enrich_step_params(
        object = object,
        assay = from_data,
        fun = invocation$fun,
        params = logged_params,
        sample_ids = colnames(base_m)
    )
    raw_result <- do.call(invocation$fun, c(list(base_m), call_params))
    result_parts <- .pb_step_result_matrix_parts(
        raw_result,
        paste0("Step '", step, "' result")
    )
    res_m <- .pb_adapter_validate_output(
        result_parts$data,
        input_matrix = base_m,
        missing = "keep",
        allow_unnamed_features = TRUE
    )
    if (!identical(colnames(res_m), colnames(base_m))) {
        stop(
            "Step '", step,
            "' result must preserve every input sample in order. ",
            "Use `pb_subset_samples()` to change the sample set.",
            call. = FALSE
        )
    }
    retry_status <- .pb_target_retry_status(
        object = object,
        to = to,
        from = from,
        data = res_m,
        step = step,
        fun = fun_name,
        params = logged_params,
        pkg = fun_package
    )

    saved_assay <- NULL
    if (store) {
        mat <- .pb_materialize_matrix(res_m, backend = backend, hdf5_path = hdf5_path)
        if (identical(retry_status, "stored_idempotent")) {
            saved_assay <- to
        } else {
            cd_from <- .pb_coldata_for_assay(object, from_data)
            se <- SummarizedExperiment(
                assays  = list(intensity = mat),
                colData = cd_from
            )
            if (isTRUE(result_parts$structured)) {
                se <- .pb_attach_step_artifacts(
                    se,
                    result_parts$artifacts
                )
            }
            object <- .pb_add_assay_with_link(
                object,
                se,
                to = to,
                from = from_data,
                step = step,
                fun = fun_name,
                params = logged_params,
                pkg = fun_package,
                lineage_from = from
            )
            saved_assay <- to
        }
    }

    object <- .pb_add_log_entry(
        object,
        step = step,
        fun = fun_name,
        from = from,
        to = to,
        params = logged_params,
        pkg = fun_package
    )

    list(object = object, assay = saved_assay, matrix = res_m, to = to)
}

# ---------------------------
# Sequences / Pipelines
# ---------------------------

#' Compute a pipeline and optionally store only the final result
#' @md
#' @param object A `ProBatchFeatures` object.
#' @param from Assay name to start the pipeline from.
#' @param steps character vector, e.g. c("log2","medianNorm","combat")
#' @param funs Optional same-length vector or list of functions or registered
#'   canonical names or aliases (default: `steps`). Registered aliases are
#'   resolved to their canonical names for operation logging.
#' @param params_list list of parameter lists (same length as steps)
#' @param level Optional level label to assign to the generated assay(s).
#' @param store_fast_steps Logical; if `FALSE`, intermediate steps listed in
#'   `fast_steps` are computed but not stored. A final `log` or `log2` step is
#'   also virtual unless `final_name` is supplied. Setting
#'   `store_fast_steps = TRUE` materializes it; other final steps are
#'   materialized regardless.
#' @param fast_steps Steps treated as fast (default:
#'   `c("log", "log2", "medianNorm")`).
#' @param store_intermediate logical; if TRUE store every step (overrides fast behavior)
#' @param final_name Optional final assay name override in
#'   `<level>::<pipeline>` form. It must differ from `from`. Stored assays and
#'   virtual operation-log targets share one result namespace, so a conflicting
#'   explicit name is rejected rather than disambiguated unless the request is
#'   an exact idempotent retry. Supplying a name materializes the final result
#'   even for an otherwise ephemeral log step.
#' @param backend "memory","hdf5","auto"
#' @param hdf5_path Optional file path used when `backend = "hdf5"`.
#' @return ProBatchFeatures with the requested pipeline added (as log and/or assay)
#' @details Registered steps record their canonical name and provider package.
#'   When a step returns [pb_step_result()] and its assay is materialized, the
#'   result's artifacts are stored in that assay's metadata under
#'   `pb_step_artifacts`. Step results must retain all samples in input order;
#'   an input-ordered feature subset is allowed.
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
    if (!is.null(final_name)) {
        .pb_assay_parts(final_name, strict = TRUE)
        if (identical(final_name, from)) {
            stop(
                "`final_name` conflicts with the source assay '",
                from, "'.",
                call. = FALSE
            )
        }
    }

    cur_from <- from
    cur_from_data <- from
    base_m <- NULL

    for (k in seq_along(steps)) {
        step <- steps[[k]]
        fun <- funs[[k]]
        par <- params_list[[k]]

        step_label <- step
        is_fast <- .pb_is_fast_step(step_label, fast_steps)
        is_final <- k == length(steps)
        # Only "log"/"log2" are treated as ephemeral fast steps by default.
        is_ephemeral_fast <- is_fast && step_label %in% c("log", "log2")
        store_this <- if (store_intermediate) {
            TRUE
        } else if (is_ephemeral_fast) {
            if (is_final && !is.null(final_name)) TRUE else store_fast_steps
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
            store = store_this, new_level = level,
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
        }
    }
    object
}

#' Evaluate a pipeline and return the matrix, without storing
#' @md
#' @param object ProBatchFeatures
#' @param from assay name (e.g., "peptide::raw")
#' @param steps character vector, e.g. c("log2","medianNorm","combat")
#' @param funs Optional same-length vector or list of functions or registered
#'   canonical names or aliases (default: `steps`).
#' @param params_list list of parameter lists (same length as steps)
#' @return Numeric matrix (features x samples). If a step returns a
#'   [pb_step_result()], only its validated `data` component is returned;
#'   artifacts are not persisted by this non-storing API. Every step must
#'   retain all samples in input order and may return an input-ordered feature
#'   subset.
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
        invocation <- .pb_step_invocation(funs[[k]], step = steps[[k]])
        logged_params <- if (is.null(params_list[[k]])) {
            list()
        } else if (is.list(params_list[[k]])) {
            params_list[[k]]
        } else {
            list(params_list[[k]])
        }
        call_params <- .pb_enrich_step_params(
            object = object,
            assay = from,
            fun = invocation$fun,
            params = logged_params,
            sample_ids = colnames(m)
        )
        raw_result <- do.call(invocation$fun, c(list(m), call_params))
        result_parts <- .pb_step_result_matrix_parts(
            raw_result,
            paste0("Step '", steps[[k]], "' result")
        )
        next_m <- .pb_adapter_validate_output(
            result_parts$data,
            input_matrix = m,
            missing = "keep",
            allow_unnamed_features = TRUE
        )
        if (!identical(colnames(next_m), colnames(m))) {
            stop(
                "Step '", steps[[k]],
                "' result must preserve every input sample in order. ",
                "Use `pb_subset_samples()` to change the sample set.",
                call. = FALSE
            )
        }
        m <- next_m
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
#' @md
#' @param object ProBatchFeatures
#' @param from assay name (e.g., "peptide::raw")
#' @param new_matrix numeric matrix (features x samples)
#' @param mapping_df Data frame mapping `from` IDs to `to` IDs. Required only
#'   when source and target row identifiers are not identical one-to-one.
#' @param from_id Source-ID column in `mapping_df` (for example,
#'   `"Precursor.Id"`); required with `mapping_df`.
#' @param to_id Target-ID column in `mapping_df` (for example,
#'   `"Protein.Ids"`); required with `mapping_df`.
#' @param map_strategy how to resolve multiple to-ids per from-id:
#'   `"as_is"` (error if not one-to-one), `"first"` (take the first exact
#'   target), or `"longest"` (prefer the exact target with the most
#'   semicolon-separated identifiers, then the longest string, then its first
#'   occurrence).
#' @param link_var rowData variable name to use for linking (e.g., "ProteinID")
#' @param to_level e.g. "protein"
#' @param to_pipeline optional pipeline name (default carries over from 'from')
#' @param name Optional final assay name override in
#'   `<level>::<pipeline>` form. Stored assays and virtual operation-log
#'   targets share one result namespace. A collision is accepted only for an
#'   exact idempotent retry with identical data, parent, mapping, and stable
#'   lineage origin.
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
    .pb_assay_parts(from, strict = TRUE)

    # ----- Determine target assay name -----
    from_parts <- strsplit(from, "::", fixed = TRUE)[[1]]
    from_pipeline <- if (length(from_parts) >= 2) from_parts[2] else "raw"
    pipeline <- to_pipeline %||% from_pipeline
    if (!is.null(name)) {
        .pb_assay_parts(name, strict = TRUE)
    }
    to <- if (!is.null(name)) name else .pb_assay_name(to_level, pipeline)
    .pb_assay_parts(to, strict = TRUE)

    # ----- Build SE for new_matrix and align samples to 'from' assay -----
    if (!(from %in% names(object))) {
        stop("Assay '", from, "' is not stored in the object.")
    }
    m <- as.matrix(new_matrix)
    if (!is.numeric(m)) {
        stop("`new_matrix` must be a numeric matrix.", call. = FALSE)
    }
    if (!nrow(m) || !ncol(m)) {
        stop(
            "`new_matrix` must contain at least one feature and one sample.",
            call. = FALSE
        )
    }
    rownames(m) <- .pb_validate_identifiers(
        rownames(m),
        "`new_matrix` feature axis"
    )
    colnames(m) <- .pb_validate_identifiers(
        colnames(m),
        "`new_matrix` sample axis"
    )
    r_from <- .pb_validate_identifiers(
        rownames(assay(object[[from]], "intensity")),
        paste0("Source assay '", from, "' feature axis")
    )
    r_to <- rownames(m)
    one_to_one <- identical(r_from, r_to)
    virtual_target <- !(to %in% names(object)) &&
        nrow(object@oplog) &&
        any(as.character(object@oplog$to) == to)
    if (virtual_target) {
        stop(
            "Result target '", to,
            "' is reserved by the operation log and cannot be used for ",
            "an added level.",
            call. = FALSE
        )
    }

    parent_ids <- NULL
    parent_keys <- NULL
    if (!one_to_one) {
        if (is.null(mapping_df) || is.null(from_id) || is.null(to_id)) {
            stop(
                "Provide mapping_df, from_id and to_id to establish ",
                "cross-level links."
            )
        }
        mapping_df <- as.data.frame(mapping_df, stringsAsFactors = FALSE)
        if (anyDuplicated(names(mapping_df))) {
            stop("`mapping_df` must have unique column names.", call. = FALSE)
        }
        mapping_id_args <- list(from_id = from_id, to_id = to_id)
        invalid_id_arg <- vapply(mapping_id_args, function(value) {
            !is.character(value) ||
                length(value) != 1L ||
                is.na(value) ||
                !nzchar(value)
        }, logical(1))
        if (any(invalid_id_arg)) {
            stop(
                "`", names(mapping_id_args)[which(invalid_id_arg)[1L]],
                "` must be one non-empty mapping column name.",
                call. = FALSE
            )
        }
        missing_mapping_columns <- setdiff(
            unlist(mapping_id_args, use.names = FALSE),
            names(mapping_df)
        )
        if (length(missing_mapping_columns)) {
            stop(
                "`mapping_df` is missing column(s): ",
                paste(missing_mapping_columns, collapse = ", "),
                ".",
                call. = FALSE
            )
        }
        mapping_from <- .pb_validate_identifiers(
            mapping_df[[from_id]],
            paste0("`mapping_df$", from_id, "`"),
            require_unique = FALSE
        )
        mapping_to <- .pb_validate_identifiers(
            mapping_df[[to_id]],
            paste0("`mapping_df$", to_id, "`"),
            require_unique = FALSE
        )

        r_from_set <- unique(r_from)
        r_to_set <- unique(r_to)
        miss_from <- setdiff(
            r_from_set,
            mapping_from
        )
        if (length(miss_from)) {
            stop(
                "Mapping incomplete: ", length(miss_from),
                " parent IDs from '", from, "' have no mapping. Examples: ",
                paste(head(miss_from, 10), collapse = ", ")
            )
        }
        keep_mapping <- mapping_from %in% r_from_set &
            mapping_to %in% r_to_set
        mapping <- data.frame(
            from_key = mapping_from[keep_mapping],
            to_key = mapping_to[keep_mapping],
            stringsAsFactors = FALSE
        )
        by_parent <- split(mapping$to_key, mapping$from_key)
        by_parent <- lapply(by_parent, unique)
        picked <- vapply(
            by_parent,
            .choose_target,
            character(1),
            r_to = r_to,
            map_strategy = map_strategy
        )

        parent_ids <- rownames(object[[from]])
        parent_keys <- unname(picked[parent_ids])
        if (anyNA(parent_keys)) {
            bad <- parent_ids[is.na(parent_keys)]
            stop(
                "Linking failed: ", length(bad),
                " parents have no exact match under map_strategy='",
                map_strategy, "'. Examples: ",
                paste(head(bad, 10), collapse = ", ")
            )
        }
    }

    origin <- if (one_to_one) {
        list(
            step = sprintf("add_level(%s)_1to1", to_level),
            fun = "addAssayLinkOneToOne",
            params = list()
        )
    } else {
        list(
            step = sprintf("add_level(%s)_byVar", to_level),
            fun = "addAssayLink",
            params = list(
                varFrom = link_var,
                varTo = link_var,
                map_strategy = map_strategy,
                from_id = from_id,
                to_id = to_id,
                mapping = stats::setNames(parent_keys, parent_ids)
            )
        )
    }
    retry_status <- .pb_target_retry_status(
        object = object,
        to = to,
        from = from,
        data = m,
        step = origin$step,
        fun = origin$fun,
        params = origin$params,
        pkg = "proBatch"
    )
    if (identical(retry_status, "stored_idempotent")) {
        message("Assay '", to, "' is an exact existing retry; skipping addition.")
        return(object)
    }

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
    if (one_to_one) {
        object <- addAssayLinkOneToOne(object, from = from, to = to)
        object <- .pb_add_log_entry(object,
            step = origin$step,
            fun = origin$fun,
            from = from, to = to,
            params = origin$params
        )
        return(object)
    }

    # Case B: many-to-one (typical peptide->protein) using addAssayLink()
    # We need a variable present in rowData(from) and rowData(to) to match on.
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
        step = origin$step,
        fun = origin$fun,
        from = from, to = to,
        params = origin$params
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
#' @md
#'
#' Ensures the `[` method returns a `ProBatchFeatures` instance so the
#' subclass-specific slots remain available after subsetting.
#'
#' @param x A `ProBatchFeatures` object.
#' @param i Row indices passed to the underlying `QFeatures` subset.
#' @param j Column indices passed to the underlying `QFeatures` subset.
#' @param k Assay indices passed to the underlying `QFeatures` subset.
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
        out <- callNextMethod(
            x = x,
            i = i,
            j = j,
            k = k,
            ...,
            drop = drop
        )
        .as_ProBatchFeatures(out, from = x)
    }
)

#' Subset samples using `ProBatchFeatures` metadata.
#' @md
#'
#' Resolve sample metadata from `colData(object)`, build a column mask, and
#' delegate to the subclass-preserving `[` method so every assay is subset
#' consistently.
#'
#' @param object A `ProBatchFeatures` object.
#' @param sample_id_col Character scalar naming the sample identifier column in
#'   `colData(object)`. When the column is absent, row names are used.
#' @param subset_by Character scalar naming the metadata column used for
#'   filtering. Defaults to `sample_id_col`.
#' @param subset_values Values to retain from `subset_by`.
#'
#' @return A `ProBatchFeatures` object restricted to the selected samples.
#' @examples
#' data(example_proteome_matrix, package = "proBatch")
#' data(example_sample_annotation, package = "proBatch")
#' sample_annotation <- example_sample_annotation
#' sample_annotation$Group <- rep(
#'     c("Pool", "Study"),
#'     length.out = nrow(sample_annotation)
#' )
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
    if (!is(object, "ProBatchFeatures")) {
        stop("`object` must be a `ProBatchFeatures` object.", call. = FALSE)
    }
    if (!is.character(sample_id_col) ||
        length(sample_id_col) != 1L ||
        is.na(sample_id_col) ||
        !nzchar(sample_id_col)) {
        stop(
            "`sample_id_col` must be a non-empty character scalar.",
            call. = FALSE
        )
    }
    if (!is.character(subset_by) ||
        length(subset_by) != 1L ||
        is.na(subset_by) ||
        !nzchar(subset_by)) {
        stop(
            "`subset_by` must be a non-empty character scalar.",
            call. = FALSE
        )
    }
    if (missing(subset_values)) {
        stop("`subset_values` must be provided.", call. = FALSE)
    }

    sample_annotation <- .pb_default_sample_annotation(
        object = object,
        sample_id_col = sample_id_col
    )
    if (!subset_by %in% colnames(sample_annotation)) {
        stop(
            "Column '",
            subset_by,
            "' not found in `colData(object)`.",
            call. = FALSE
        )
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
