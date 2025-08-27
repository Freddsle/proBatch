
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
        chain = "character", # to record processing steps
        oplog = "DataFrame"
    )
)

methods::setValidity("ProBatchFeatures", function(object) {
# Basic checks for oplog columns
req <- c("step","fun","from","to","params","timestamp","pkg")
if (!all(req %in% colnames(object@oplog)))
return("oplog lacks required columns")
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
ProBatchFeatures <- function(data_matrix,
                             sample_annotation = NULL,
                             sample_id_col = "FullRunName",
                             name = "raw") {
    stopifnot(is.matrix(data_matrix) || is.data.frame(data_matrix))
    data_matrix <- as.matrix(data_matrix)

    # Align colData
    if (!is.null(sample_annotation)) {
        sa <- as.data.frame(sample_annotation, stringsAsFactors = FALSE)
        if (is.null(rownames(sa))) {
            if (!sample_id_col %in% colnames(sa)) {
                stop("Provide rownames(sample_annotation) or a valid sample_id_col.")
            }
            rownames(sa) <- sa[[sample_id_col]]
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
    # Use QFeatures constructor to make a QFeatures, then coerce to ProBatchFeatures
    qf <- QFeatures::QFeatures(setNames(list(se), name))

    # start with empty structured oplog
    empty_log <- S4Vectors::DataFrame(
        step = character(),
        fun = character(),
        from = character(),
        to = character(),
        params = I(vector("list", 0L)),
        timestamp = as.POSIXct(character()),
        pkg = character()
    )

    methods::new(
        "ProBatchFeatures",
        qf,
        chain = character(), # start with empty history
        oplog = empty_log
    )
}

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
    tail(QFeatures::names(object), 1)
}

#' Convenience accessor for assay matrix by name/index (returns the 'intensity' assay)
#' @export
pb_assay_matrix <- function(object, assay = pb_current_assay(object)) {
    QFeatures::assay(object, assay)
}


#---------------------------
# Internal helpers (not exported)
#---------------------------

.pb_add_log_entry <- function(object, step, fun, from, to, params, pkg = "proBatch") {
    entry <- S4Vectors::DataFrame(
        step = step,
        fun = fun,
        from = as.character(from),
        to = as.character(to),
        params = I(list(params)),
        timestamp = Sys.time(),
        pkg = pkg
    )
    object@oplog <- S4Vectors::rbind(object@oplog, entry)
    object@chain <- c(object@chain, step)
    object
}

# Add assay and try to link 1:1 to parent (safe for log/norm/ComBat).
# Filtering may change row counts; linking will be skipped if not 1:1.
.pb_add_assay_with_link <- function(object, se, name, from) {
    object <- QFeatures::addAssay(object, se, name = name)

    # Try one-to-one link if same rownames & nrow
    ok_link <- FALSE
    if (all(rownames(SummarizedExperiment::assay(se, 1)) %in%
        rownames(QFeatures::experiment(object, from))) &&
        nrow(SummarizedExperiment::assay(se, 1)) ==
            nrow(QFeatures::experiment(object, from))) {
        # robust: use tryCatch to avoid breaking pipelines on older QFeatures
        object <- tryCatch(
            QFeatures::addAssayLinkOneToOne(object, from = from, to = name),
            error = function(e) object
        )
        ok_link <- TRUE
    }
    attr(object, "linked_last") <- ok_link
    object
}

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
