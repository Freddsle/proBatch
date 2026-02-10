.pb_resolve_assay_for_input <- function(object, pbf_name = NULL) {
    stopifnot(is(object, "ProBatchFeatures"))

    if (is.null(pbf_name) || !length(pbf_name)) {
        return(pb_current_assay(object))
    }

    if (length(pbf_name) != 1L || !nzchar(as.character(pbf_name[[1]]))) {
        stop("`pbf_name` must contain exactly one non-empty assay name.")
    }

    as.character(pbf_name[[1]])
}

.pb_default_sample_annotation <- function(object,
                                          sample_annotation = NULL,
                                          sample_id_col = "FullRunName",
                                          sample_ids = NULL) {
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
