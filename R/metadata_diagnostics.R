#' Identify duplicated columns in metadata or assays
#'
#' Compare columns cell-by-cell and report the pairs that contain identical
#' values, including matching `NA` placements. The method works on plain data
#' frames as well as on [`ProBatchFeatures`] objects by extracting the relevant
#' annotation tables.
#'
#' @param x A data frame-like or `ProBatchFeatures` object to inspect.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A list of character vectors (each of length two) naming duplicated
#'   column pairs. An empty list is returned when no duplicates are detected.
#'
#' @examples
#' metadata <- data.frame(
#'     sample = LETTERS[1:4],
#'     batch = c("A", "B", "A", "B"),
#'     batch_copy = c("A", "B", "A", "B"),
#'     check.names = FALSE
#' )
#' find_duplicated_columns(metadata)
#'
#' @export
find_duplicated_columns <- function(x, ...) {
    UseMethod("find_duplicated_columns")
}

#' @rdname find_duplicated_columns
#' @method find_duplicated_columns default
#' @export
find_duplicated_columns.default <- function(x, ...) {
    if (!is.data.frame(x)) {
        x <- as.data.frame(x)
    }

    if (ncol(x) < 2L) {
        return(list())
    }

    col_names <- colnames(x)
    if (is.null(col_names)) {
        col_names <- paste0("V", seq_len(ncol(x)))
    }

    groups <- .pb_duplicate_groups(x, col_names)

    if (!length(groups)) {
        return(list())
    }

    duplicates <- vector("list", length = 0L)
    for (cols in groups) {
        combs <- utils::combn(cols, 2L, simplify = FALSE)
        duplicates <- c(duplicates, combs)
    }

    duplicates
}

#' @rdname find_duplicated_columns
#' @param component Character string selecting which part of the object to
#'   inspect. Defaults to sample metadata (`"colData"`); use `"rowData"` to
#'   inspect feature metadata of a specific assay.
#' @param assay Optional single assay name used when `component = "rowData"`. If
#'   omitted, the most recent assay returned by [pb_current_assay()] is used.
#' @param df Optional data frame. When supplied, `component` and `assay` are
#'   ignored and the provided data frame is analysed instead.
#' @method find_duplicated_columns ProBatchFeatures
#' @export
find_duplicated_columns.ProBatchFeatures <- function(
    x,
    component = c("colData", "rowData"),
    assay = NULL,
    df = NULL,
    ...) {
    object <- x

    if (!is.null(df)) {
        if (!is.data.frame(df)) {
            stop("`df` must be a data frame when supplied.")
        }
        return(find_duplicated_columns.default(df, ...))
    }

    component <- match.arg(component)
    if (component == "colData") {
        target <- SummarizedExperiment::colData(object)
    } else {
        chosen_assay <- assay
        if (is.null(chosen_assay)) {
            chosen_assay <- pb_current_assay(object)
        }
        if (!length(chosen_assay)) {
            stop("Provide an assay name via `assay` or ensure the object stores assays.")
        }
        if (length(chosen_assay) != 1L) {
            stop("`assay` must be a single assay name.")
        }
        if (!chosen_assay %in% names(object)) {
            stop("Assay '", chosen_assay, "' not found in object.")
        }
        target <- SummarizedExperiment::rowData(object[[chosen_assay]])
    }

    target_df <- .pb_to_base_df(target)
    find_duplicated_columns.default(target_df, ...)
}

#' Summarise metadata column cardinality and missingness
#'
#' Produce a per-column overview reporting how many distinct non-missing values
#' appear in each column together with the number and percentage of `NA`
#' entries. The summary is available for plain data frames and
#' [`ProBatchFeatures`] objects.
#'
#' @param x A data frame or `ProBatchFeatures` object to inspect.
#' @param sort Logical flag indicating whether to order the resulting table by
#'   descending `pct_NA` and then ascending `n_unique`. Defaults to `TRUE` to
#'   mirror console usage in batch diagnostics.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A data frame with columns `colname`, `n_unique`, `n_NA`, and
#'   `pct_NA`. Percentages are expressed on a 0–100 scale. When `sort = TRUE` the
#'   rows are ordered by decreasing missingness.
#'
#' @examples
#' metadata <- data.frame(
#'     sample = LETTERS[1:4],
#'     batch = c("A", "B", "A", "B"),
#'     batch_copy = c("A", "B", "A", "B"),
#'     check.names = FALSE
#' )
#' metadata_column_summary(metadata)
#'
#' @export
metadata_column_summary <- function(x, ...) {
    UseMethod("metadata_column_summary")
}

#' @rdname metadata_column_summary
#' @method metadata_column_summary default
#' @export
metadata_column_summary.default <- function(x, sort = TRUE, ...) {
    if (!is.data.frame(x)) {
        x <- as.data.frame(x)
    }

    if (!ncol(x)) {
        return(data.frame(
            colname = character(),
            n_unique = integer(),
            n_NA = integer(),
            pct_NA = numeric(),
            stringsAsFactors = FALSE
        ))
    }

    col_names <- colnames(x)
    if (is.null(col_names)) {
        col_names <- paste0("V", seq_len(ncol(x)))
    }

    n_unique <- vapply(x, function(col) length(unique(col[!is.na(col)])), integer(1))
    n_na <- vapply(x, function(col) sum(is.na(col)), integer(1))
    pct_na <- if (nrow(x)) (n_na / nrow(x)) * 100 else rep(NaN, length(n_na))

    summary_df <- data.frame(
        colname = col_names,
        n_unique = n_unique,
        n_NA = n_na,
        pct_NA = pct_na,
        stringsAsFactors = FALSE
    )

    if (sort) {
        ord <- order(-summary_df$pct_NA, summary_df$n_unique, summary_df$colname)
        summary_df <- summary_df[ord, , drop = FALSE]
        rownames(summary_df) <- NULL
    }

    summary_df
}

#' @rdname metadata_column_summary
#' @param component Character string selecting which part of the object to
#'   inspect. Defaults to sample metadata (`"colData"`); use `"rowData"` to
#'   inspect feature metadata of a specific assay.
#' @param assay Optional single assay name used when `component = "rowData"`. If
#'   omitted, the most recent assay returned by [pb_current_assay()] is used.
#' @param df Optional data frame. When supplied, `component` and `assay` are
#'   ignored and the provided data frame is analysed instead.
#' @method metadata_column_summary ProBatchFeatures
#' @export
metadata_column_summary.ProBatchFeatures <- function(
    x,
    component = c("colData", "rowData"),
    assay = NULL,
    df = NULL,
    sort = TRUE,
    ...) {
    object <- x

    if (!is.null(df)) {
        if (!is.data.frame(df)) {
            stop("`df` must be a data frame when supplied.")
        }
        return(metadata_column_summary.default(df, sort = sort, ...))
    }

    component <- match.arg(component)
    if (component == "colData") {
        target <- SummarizedExperiment::colData(object)
    } else {
        chosen_assay <- assay
        if (is.null(chosen_assay)) {
            chosen_assay <- pb_current_assay(object)
        }
        if (!length(chosen_assay)) {
            stop("Provide an assay name via `assay` or ensure the object stores assays.")
        }
        if (length(chosen_assay) != 1L) {
            stop("`assay` must be a single assay name.")
        }
        if (!chosen_assay %in% names(object)) {
            stop("Assay '", chosen_assay, "' not found in object.")
        }
        target <- SummarizedExperiment::rowData(object[[chosen_assay]])
    }

    target_df <- .pb_to_base_df(target)
    metadata_column_summary.data_frame(target_df, sort = sort, ...)
}

#' Filter metadata columns based on duplication and completeness
#'
#' Automate common metadata cleaning steps by dropping duplicated columns while
#' retaining the preferred representative and optionally removing columns with
#' excessive missingness. Works on plain data frames as well as
#' [`ProBatchFeatures`] objects.
#'
#' @param x A data frame or `ProBatchFeatures` object.
#' @param duplicate_keep Strategy for selecting which column to retain within
#'   each duplicated set. Choose between keeping the `"first"`, `"last"`, or the
#'   first column whose name matches `duplicate_pattern` (`"pattern"`).
#' @param duplicate_pattern Character pattern used when
#'   `duplicate_keep = "pattern"`. A regular expression is expected unless
#'   `pattern_fixed = TRUE`.
#' @param pattern_ignore_case Logical flag passed to [base::grepl()] when
#'   matching `duplicate_pattern`. Defaults to `TRUE` as metadata column names are
#'   often case-insensitive.
#' @param pattern_fixed Forwarded to the `fixed` argument of [base::grepl()].
#'   Defaults to `FALSE` to allow regular expressions.
#' @param min_non_na Optional minimum number of non-missing values a column must
#'   contain to be retained. Set to `NULL` to skip this filter.
#' @param max_pct_na Optional maximum percentage (0–100) of missing values a
#'   column may contain to be retained. Set to `NULL` to disable.
#' @param sort Logical flag forwarded to [metadata_column_summary()] while
#'   computing missingness statistics. Has no impact on the returned data frame.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return For data frames, a filtered data frame with attributes
#'   `"dropped_columns"`, `"dropped_duplicates"`, and `"dropped_missing"`
#'   describing the removal decisions. For `ProBatchFeatures` objects, either a
#'   filtered data frame (default) or the modified object when `inplace = TRUE`.
#'
#' @examples
#' metadata <- data.frame(
#'     `factor value intensity` = c(1, 2, 3),
#'     intensity = c(1, 2, 3),
#'     batch = c("A", NA, "B"),
#'     stringsAsFactors = FALSE,
#'     check.names = FALSE
#' )
#' filtered <- filter_metadata_columns(
#'     metadata,
#'     duplicate_keep = "pattern",
#'     duplicate_pattern = "^factor value",
#'     min_non_na = 2
#' )
#' names(filtered)
#'
#' @export
filter_metadata_columns <- function(x, ...) {
    UseMethod("filter_metadata_columns")
}

#' @rdname filter_metadata_columns
#' @method filter_metadata_columns default
#' @export
filter_metadata_columns.default <- function(
    x,
    duplicate_keep = c("first", "last", "pattern"),
    duplicate_pattern = NULL,
    pattern_ignore_case = TRUE,
    pattern_fixed = FALSE,
    min_non_na = NULL,
    max_pct_na = NULL,
    sort = FALSE,
    ...) {
    if (!is.data.frame(x)) {
        x <- as.data.frame(x)
    }

    duplicate_keep <- match.arg(duplicate_keep)
    if (duplicate_keep == "pattern" && (is.null(duplicate_pattern) || !nzchar(duplicate_pattern))) {
        stop("Provide `duplicate_pattern` when `duplicate_keep = 'pattern'`.")
    }

    if (!is.null(min_non_na)) {
        if (!is.numeric(min_non_na) || length(min_non_na) != 1L || is.na(min_non_na)) {
            stop("`min_non_na` must be a single non-missing numeric value.")
        }
        if (min_non_na < 0) {
            stop("`min_non_na` must be >= 0.")
        }
    }

    if (!is.null(max_pct_na)) {
        if (!is.numeric(max_pct_na) || length(max_pct_na) != 1L || is.na(max_pct_na)) {
            stop("`max_pct_na` must be a single non-missing numeric value.")
        }
        if (max_pct_na < 0 || max_pct_na > 100) {
            stop("`max_pct_na` must be between 0 and 100.")
        }
    }

    original_colnames <- colnames(x)
    if (is.null(original_colnames)) {
        col_names <- paste0("V", seq_len(ncol(x)))
    } else {
        col_names <- original_colnames
    }

    groups <- if (ncol(x) >= 2L) .pb_duplicate_groups(x, col_names) else list()
    dropped_duplicates <- character(0)

    if (length(groups)) {
        for (cols in groups) {
            keep <- switch(duplicate_keep,
                first = cols[[1]],
                last = cols[[length(cols)]],
                pattern = {
                    hits <- cols[grepl(duplicate_pattern, cols, ignore.case = pattern_ignore_case, fixed = pattern_fixed)]
                    if (length(hits)) {
                        hits[[1]]
                    } else {
                        cols[[1]]
                    }
                }
            )
            dropped_duplicates <- c(dropped_duplicates, setdiff(cols, keep))
        }
        dropped_duplicates <- unique(dropped_duplicates)
    }
    
    dropped_missing <- character(0)
    if (!is.null(min_non_na) || !is.null(max_pct_na)) {
        missing_summary <- metadata_column_summary.default(x, sort = sort)
        missing_summary$n_non_na <- nrow(x) - missing_summary$n_NA
        if (!is.null(min_non_na)) {
            dropped_missing <- union(dropped_missing, missing_summary$colname[missing_summary$n_non_na < min_non_na])
        }
        if (!is.null(max_pct_na)) {
            dropped_missing <- union(dropped_missing, missing_summary$colname[missing_summary$pct_NA > max_pct_na])
        }
    }
    
    if (length(dropped_duplicates)) {
        message(sprintf("Dropping %s duplicated columns.", length(dropped_duplicates)))
    }
    if (length(dropped_missing)) {
        message(sprintf("Dropping %s columns based on missingness thresholds.", length(dropped_missing)))
    }

    drop_candidates <- union(dropped_duplicates, dropped_missing)

    if (!length(drop_candidates)) {
        attr(x, "dropped_columns") <- character(0)
        attr(x, "dropped_duplicates") <- dropped_duplicates
        attr(x, "dropped_missing") <- dropped_missing
        return(x)
    }

    keep_mask <- !(col_names %in% drop_candidates)
    keep_cols <- col_names[keep_mask]
    if (!length(keep_cols)) {
        stop("All columns would be removed; adjust filtering thresholds or duplicate strategy.")
    }

    filtered <- if (inherits(x, "data.table")) {
        x[, ..keep_cols]
    } else {
        x[, keep_cols, drop = FALSE]
    }

    attr(filtered, "dropped_columns") <- col_names[!keep_mask]
    attr(filtered, "dropped_duplicates") <- dropped_duplicates
    attr(filtered, "dropped_missing") <- dropped_missing

    filtered
}

#' @rdname filter_metadata_columns
#' @param component Character string selecting which part of the object to
#'   inspect. Defaults to sample metadata (`"colData"`); use `"rowData"` to
#'   inspect feature metadata of a specific assay.
#' @param assay Optional single assay name used when `component = "rowData"`. If
#'   omitted, the most recent assay returned by [pb_current_assay()] is used.
#' @param df Optional data frame. When supplied, `component` and `assay` are
#'   ignored and the provided data frame is analysed instead.
#' @param inplace Logical flag; when `TRUE`, the filtered metadata replaces the
#'   corresponding component inside the object. Defaults to `FALSE` to keep the
#'   original object untouched.
#' @method filter_metadata_columns ProBatchFeatures
#' @export
filter_metadata_columns.ProBatchFeatures <- function(
    x,
    component = c("colData", "rowData"),
    assay = NULL,
    df = NULL,
    inplace = FALSE,
    ...) {
    object <- x

    if (!is.null(df)) {
        if (!is.data.frame(df)) {
            stop("`df` must be a data frame when supplied.")
        }
        return(filter_metadata_columns.default(df, ...))
    }

    component <- match.arg(component)
    if (component == "colData") {
        target <- SummarizedExperiment::colData(object)
        base_target <- as.data.frame(target)
        row_ids <- rownames(base_target)
        filtered_df <- filter_metadata_columns.default(base_target, ...)
        kept <- colnames(filtered_df)
        if (inplace) {
            filtered_target <- target[, kept, drop = FALSE]
            SummarizedExperiment::colData(object) <- filtered_target
            info <- list(
                component = "colData",
                assay = NA_character_,
                dropped_columns = attr(filtered_df, "dropped_columns"),
                dropped_duplicates = attr(filtered_df, "dropped_duplicates"),
                dropped_missing = attr(filtered_df, "dropped_missing")
            )
            meta <- metadata(object)
            meta$filter_metadata_columns <- info
            metadata(object) <- meta
            return(object)
        }
        rownames(filtered_df) <- row_ids
        return(filtered_df)
    }

    chosen_assay <- assay
    if (is.null(chosen_assay)) {
        chosen_assay <- pb_current_assay(object)
    }
    if (!length(chosen_assay)) {
        stop("Provide an assay name via `assay` or ensure the object stores assays.")
    }
    if (length(chosen_assay) != 1L) {
        stop("`assay` must be a single assay name.")
    }
    if (!chosen_assay %in% names(object)) {
        stop("Assay '", chosen_assay, "' not found in object.")
    }

    target <- SummarizedExperiment::rowData(object[[chosen_assay]])
    base_target <- as.data.frame(target)
    row_ids <- rownames(base_target)
    filtered_df <- filter_metadata_columns.default(base_target, ...)
    kept <- colnames(filtered_df)

    if (inplace) {
        filtered_target <- target[, kept, drop = FALSE]
        SummarizedExperiment::rowData(object[[chosen_assay]]) <- filtered_target
        info <- list(
            component = "rowData",
            assay = chosen_assay,
            dropped_columns = attr(filtered_df, "dropped_columns"),
            dropped_duplicates = attr(filtered_df, "dropped_duplicates"),
            dropped_missing = attr(filtered_df, "dropped_missing")
        )
        meta <- metadata(object)
        meta$filter_metadata_columns <- info
        metadata(object) <- meta
        return(object)
    }

    rownames(filtered_df) <- row_ids
    filtered_df
}

.pb_normalise_column <- function(column) {
    col <- column
    attributes(col) <- NULL
    if (is.integer(col)) {
        storage.mode(col) <- "double"
    }
    col
}

.pb_column_signature <- function(column) {
    raw_column <- serialize(column, connection = NULL, ascii = FALSE)
    paste(as.integer(raw_column), collapse = ",")
}

.pb_duplicate_groups <- function(x, col_names) {
    if (!length(col_names)) {
        return(list())
    }
    normalised <- lapply(x, .pb_normalise_column)
    names(normalised) <- col_names
    signatures <- vapply(normalised, .pb_column_signature, character(1), USE.NAMES = FALSE)
    groups <- split(col_names, signatures)
    Filter(function(cols) length(cols) > 1L, groups)
}

.pb_to_base_df <- function(x) {
    if (is.null(x)) {
        return(data.frame())
    }
    if (inherits(x, "data.frame")) {
        return(x)
    }
    if (inherits(x, "DataFrame")) {
        return(as.data.frame(x))
    }
    stop("Object cannot be coerced to a data frame.")
}
