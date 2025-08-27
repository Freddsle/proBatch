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
        chain = "character", # ordered processing steps
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
    TRUE
})

# ---------------------------
# Constructors
# ---------------------------

#' Construct a ProBatchFeatures object from a wide matrix + sample annotation
#'
#' @param data_matrix numeric matrix (features x samples) or data.frame
#' @param sample_annotation data.frame with sample metadata (rows = samples)
#' @param sample_id_col character(1), column in sample_annotation that must match colnames(data_matrix).
#' If missing, rownames(sample_annotation) are used.
#' @param name character(1) name of the initial assay added to the object (default: "raw").
#' @return A `ProBatchFeatures` object.
#' @export
ProBatchFeatures <- function(
    data_matrix,
    sample_annotation = NULL,
    sample_id_col = "FullRunName",
    name = "raw") {
    stopifnot(is.matrix(data_matrix) || is.data.frame(data_matrix))
    data_matrix <- as.matrix(data_matrix)

    if (is.null(colnames(data_matrix))) {
        stop("data_matrix must have column names (sample IDs).")
    }

    # Align colData
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

        # Reorder & restrict to present samples; error on mismatch
        if (!all(colnames(data_matrix) %in% rownames(sa))) {
            miss <- setdiff(colnames(data_matrix), rownames(sa))
            stop("Sample annotation missing for: ", paste(miss, collapse = ", "))
        }
        sa <- sa[colnames(data_matrix), , drop = FALSE]
        cd <- S4Vectors::DataFrame(sa)
    } else {
        cd <- S4Vectors::DataFrame(row.names = colnames(data_matrix))
    }

    # Create SummarizedExperiment for the raw data
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = data_matrix),
        colData = cd
    )
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
        chain = character(), # start with empty history
        oplog = empty_log
    )
}

# =========================
# Constructors for LONG data
# =========================

#' Construct a ProBatchFeatures from LONG df (reuses `long_to_matrix`).
#' @param df_long long df with at least feature_id_col, sample_id_col, measure_col
#' @param sample_annotation optional; if provided, must contain sample_id_col
#' @param feature_id_col, sample_id_col, measure_col names in df_long
#' @param name assay name (default "raw")
#' @export
ProBatchFeatures_from_long <- function(
    df_long,
    sample_annotation = NULL,
    feature_id_col = "peptide_group_label",
    sample_id_col = "FullRunName",
    measure_col = "Intensity",
    name = "raw") {
    stopifnot(is.data.frame(df_long))

    # 1) long -> wide using existing proBatch utility
    data_matrix <- proBatch::long_to_matrix(
        df_long,
        feature_id_col = feature_id_col,
        sample_id_col  = sample_id_col,
        measure_col    = measure_col
    )
    # 2) delegate to the wide constructor
    ProBatchFeatures(
        data_matrix        = data_matrix,
        sample_annotation  = sample_annotation,
        sample_id_col      = sample_id_col,
        name               = name
    )
}

# =========================
# Helpers
# =========================

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

#' Pretty pipeline name derived from chain
#' @export
pb_pipeline_name <- function(object) get_chain(object, as_string = TRUE)

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
# =========================
# Format shims
# =========================

#' Get current assay as LONG (via proBatch::matrix_to_long)
#' @export
pb_as_long <- function(
    object,
    feature_id_col = "peptide_group_label",
    sample_id_col = "FullRunName",
    measure_col = "Intensity") {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    m <- pb_assay_matrix(object)
    se <- object[[pb_current_assay(object)]]
    sa <- as.data.frame(SummarizedExperiment::colData(se))

    proBatch::matrix_to_long(
        data_matrix       = m,
        sample_annotation = sa,
        feature_id_col    = feature_id_col,
        measure_col       = measure_col,
        sample_id_col     = sample_id_col
    )
}

#' Get an assay matrix (wide) explicitly
#' @export
pb_as_wide <- function(object) {
    stopifnot(methods::is(object, "ProBatchFeatures"))
    pb_assay_matrix(object)
}

#---------------------------
# Internal helpers (not exported)
#---------------------------

.pb_add_log_entry <- function(object, step, fun, from, to, params, pkg = "proBatch") {
    entry <- S4Vectors::DataFrame(
        step      = step,
        fun       = fun,
        from      = as.character(from),
        to        = as.character(to),
        params    = I(list(params)),
        timestamp = Sys.time(),
        pkg       = pkg
    )
    object@oplog <- S4Vectors::rbind(object@oplog, entry)
    object@chain <- c(object@chain, step)
    object
}

# Add assay and try to link 1:1 to parent.
# Filtering may change row counts; linking will be skipped if not 1:1.
.pb_add_assay_with_link <- function(object, se, name, from) {
    # Ensure new SE colData exactly matches the current global colData for these samples
    # (prevents .updateColDataFromAssays() conflicts, e.g., factor levels / POSIXct attrs)
    samp <- colnames(SummarizedExperiment::assay(se, 1))
    # MultiAssayExperiment exports colData() for QFeatures
    global_cd <- MultiAssayExperiment::colData(object)
    if (!all(samp %in% rownames(global_cd))) {
        stop("Samples in new assay not present in QFeatures colData.")
    }
    SummarizedExperiment::colData(se) <- global_cd[samp, , drop = FALSE]

    # Add assay
    object <- QFeatures::addAssay(object, se, name = name)

    # Try to link 1:1 if rownames match exactly
    ok_link <- FALSE
    from_se <- object[[from]]
    if (identical(
        rownames(SummarizedExperiment::assay(se, 1)),
        rownames(from_se)
    )) {
        had_err <- FALSE
        object <- tryCatch(
            QFeatures::addAssayLinkOneToOne(object, from = from, to = name),
            error = function(e) {
                had_err <<- TRUE
                object
            }
        )
        ok_link <- !had_err
    }

    S4Vectors::metadata(object)$linked_last <- ok_link
    object
}

# ---------------------------
# Show
# ---------------------------
methods::setMethod("show", "ProBatchFeatures", function(object) {
    methods::callNextMethod()
    ch <- get_chain(object, as_string = TRUE)
    if (nzchar(ch)) {
        cat("  Processing chain: ", ch, "\n", sep = "")
    } else {
        cat("  Processing chain: <empty>\n")
    }
    if (nrow(object@oplog)) {
        cat("  Steps logged: ", nrow(object@oplog), " (see get_operation_log())\n", sep = "")
    }
})

#---------------------------
# Show
#---------------------------
methods::setMethod("show", "ProBatchFeatures", function(object) {
    methods::callNextMethod()
    ch <- get_chain(object, as_string = TRUE)
    if (nzchar(ch)) {
        cat("  Processing chain: ", ch, "\n", sep = "")
    } else {
        cat("  Processing chain: <empty>\n")
    }
    if (nrow(object@oplog)) {
        cat("  Steps logged: ", nrow(object@oplog), " (see get_operation_log())\n", sep = "")
    }
})
