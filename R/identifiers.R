.pb_validate_identifiers <- function(
    identifiers,
    context,
    require_unique = TRUE
) {
    if (is.null(identifiers)) {
        stop(context, " must be named.", call. = FALSE)
    }

    identifiers <- as.character(identifiers)
    if (anyNA(identifiers) || any(!nzchar(identifiers))) {
        stop(context, " contains NA or empty identifiers.", call. = FALSE)
    }

    if (require_unique) {
        duplicates <- unique(identifiers[duplicated(identifiers)])
        if (length(duplicates)) {
            stop(
                context,
                " contains duplicate identifiers: ",
                paste(utils::head(duplicates, 10L), collapse = ", "),
                ".",
                call. = FALSE
            )
        }
    }

    identifiers
}

.pb_validate_long_keys <- function(
    value,
    feature_id_col,
    sample_id_col,
    context = "Long input"
) {
    if (!is.data.frame(value)) {
        stop(context, " must be a data frame.", call. = FALSE)
    }
    if (anyDuplicated(names(value))) {
        stop(context, " must have unique column names.", call. = FALSE)
    }

    required <- c(feature_id_col, sample_id_col)
    missing_columns <- setdiff(required, names(value))
    if (length(missing_columns)) {
        stop(
            context,
            " is missing identifier column(s): ",
            paste(missing_columns, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    feature_ids <- .pb_validate_identifiers(
        value[[feature_id_col]],
        paste0(context, " feature identifiers"),
        require_unique = FALSE
    )
    sample_ids <- .pb_validate_identifiers(
        value[[sample_id_col]],
        paste0(context, " sample identifiers"),
        require_unique = FALSE
    )
    keys <- data.frame(
        feature = feature_ids,
        sample = sample_ids,
        stringsAsFactors = FALSE
    )
    duplicate_keys <- duplicated(keys)
    if (any(duplicate_keys)) {
        examples <- unique(keys[duplicate_keys, , drop = FALSE])
        examples <- examples[seq_len(min(10L, nrow(examples))), , drop = FALSE]
        stop(
            context,
            " contains duplicate feature/sample keys: ",
            paste(
                paste0(examples$feature, "/", examples$sample),
                collapse = ", "
            ),
            ". Supply one value per feature/sample pair.",
            call. = FALSE
        )
    }

    invisible(value)
}
