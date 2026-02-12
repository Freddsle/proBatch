.pb_resolve_assay_for_input <- function(object,
                                        pbf_name = NULL,
                                        inform_if_default = FALSE) {
    stopifnot(is(object, "ProBatchFeatures"))

    if (is.null(pbf_name) || !length(pbf_name)) {
        assay_name <- pb_current_assay(object)
        if (isTRUE(inform_if_default)) {
            message("`pbf_name` not provided, using the most recent assay: ", assay_name)
        }
        return(assay_name)
    }

    if (length(pbf_name) != 1L) {
        stop("`pbf_name` must contain exactly one non-empty assay name.")
    }

    assay_name <- as.character(pbf_name[[1]])
    if (length(assay_name) != 1L || is.na(assay_name) || !nzchar(assay_name)) {
        stop("`pbf_name` must contain exactly one non-empty assay name.")
    }

    assay_name
}

.pb_resolve_assays_for_input <- function(object,
                                         pbf_name = NULL,
                                         default = c("current", "all"),
                                         deduplicate = FALSE,
                                         inform_if_default = FALSE,
                                         empty_message = "Provide at least one `pbf_name` to plot.") {
    stopifnot(is(object, "ProBatchFeatures"))
    default <- match.arg(default)
    using_default <- is.null(pbf_name) || !length(pbf_name)

    if (using_default) {
        assays <- if (identical(default, "all")) names(object) else pb_current_assay(object)
    } else {
        assays <- as.character(pbf_name)
    }

    assays <- assays[!is.na(assays) & nzchar(assays)]
    if (isTRUE(deduplicate)) {
        assays <- unique(assays)
    }

    if (!length(assays)) {
        stop(empty_message)
    }

    if (using_default && isTRUE(inform_if_default)) {
        if (identical(default, "all")) {
            message("`pbf_name` not provided, using all assays: ", paste(assays, collapse = ", "))
        } else {
            message("`pbf_name` not provided, using the most recent assay: ", assays[[1]])
        }
    }

    assays
}

.pb_default_sample_annotation <- function(object,
                                          sample_annotation = NULL,
                                          sample_id_col = "FullRunName",
                                          sample_ids = NULL,
                                          drop_rownames = FALSE) {
    if (!is.null(sample_annotation)) {
        annotation <- as.data.frame(sample_annotation, stringsAsFactors = FALSE)
    } else {
        annotation <- as.data.frame(colData(object), stringsAsFactors = FALSE)
    }

    if (!sample_id_col %in% names(annotation)) {
        rn <- rownames(annotation)
        if (!is.null(rn) && !anyNA(rn) && length(rn) == nrow(annotation)) {
            annotation[[sample_id_col]] <- rn
        } else if (!is.null(sample_ids) && length(sample_ids) == nrow(annotation)) {
            annotation[[sample_id_col]] <- sample_ids
        }
    }

    if (isTRUE(drop_rownames)) {
        rownames(annotation) <- NULL
    }

    annotation
}

.pb_default_feature_annotation <- function(object,
                                           assay_name,
                                           feature_annotation = NULL,
                                           feature_id_col = "peptide_group_label") {
    if (!is.null(feature_annotation)) {
        annotation <- as.data.frame(feature_annotation, stringsAsFactors = FALSE)
    } else if (assay_name %in% names(object)) {
        annotation <- as.data.frame(rowData(object[[assay_name]]), stringsAsFactors = FALSE)
    } else {
        annotation <- NULL
    }

    if (is.null(annotation)) {
        return(NULL)
    }

    if (!feature_id_col %in% names(annotation)) {
        rn <- rownames(annotation)
        if (!is.null(rn) && !anyNA(rn) && length(rn) == nrow(annotation)) {
            annotation[[feature_id_col]] <- rn
        }
    }

    annotation
}

.pb_pbf_to_long <- function(object,
                            assay_name,
                            feature_id_col = "peptide_group_label",
                            sample_id_col = "FullRunName",
                            measure_col = "Intensity") {
    matrix_to_long(
        data_matrix = pb_assay_matrix(object, assay = assay_name),
        feature_id_col = feature_id_col,
        sample_id_col = sample_id_col,
        measure_col = measure_col
    )
}

.pb_prepare_long_inputs <- function(df_long,
                                    sample_annotation,
                                    sample_id_col,
                                    feature_id_col,
                                    measure_col,
                                    pbf_name = NULL) {
    object <- NULL
    assay_name <- NULL

    if (is(df_long, "ProBatchFeatures")) {
        object <- df_long
        assay_name <- .pb_resolve_assay_for_input(object, pbf_name = pbf_name)
        df_long <- .pb_pbf_to_long(
            object = object,
            assay_name = assay_name,
            feature_id_col = feature_id_col,
            sample_id_col = sample_id_col,
            measure_col = measure_col
        )
        sample_annotation <- .pb_default_sample_annotation(
            object = object,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            sample_ids = unique(df_long[[sample_id_col]])
        )
    }

    list(
        df_long = df_long,
        sample_annotation = sample_annotation,
        object = object,
        assay_name = assay_name
    )
}
