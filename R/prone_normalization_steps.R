# Internal registrar used from .onLoad(): registers PRONE normalization methods
# both as direct names (when no collision with existing proBatch steps) and as
# prefixed aliases (default: "PR<method>").
.pb_prone_default_step_prefix <- function() {
    "PR"
}

.pb_register_prone_normalization_steps <- function(methods = NULL,
                                                   prefix = .pb_prone_default_step_prefix(),
                                                   register_lowercase = TRUE) {
    if (!is.character(prefix) || length(prefix) != 1L || !nzchar(prefix)) {
        stop("`prefix` must be a non-empty character scalar.")
    }
    if (!is.logical(register_lowercase) || length(register_lowercase) != 1L) {
        stop("`register_lowercase` must be TRUE or FALSE.")
    }

    methods <- .pb_prone_normalization_method_names(methods)
    if (!length(methods)) {
        return(invisible(list(
            direct = character(),
            prefixed = character(),
            collisions = character()
        )))
    }

    direct_steps <- make.unique(
        vapply(methods, .pb_prone_norm_step_name, prefix = "", FUN.VALUE = character(1)),
        sep = "."
    )
    prefixed_steps <- make.unique(
        vapply(methods, .pb_prone_norm_step_name, prefix = prefix, FUN.VALUE = character(1)),
        sep = "."
    )

    registered_direct <- character()
    registered_prefixed <- character()
    collisions <- character()

    for (idx in seq_along(methods)) {
        norm_method <- methods[[idx]]
        step_fun <- .pb_make_prone_norm_step_fun(norm_method = norm_method)
        direct_step <- direct_steps[[idx]]
        prefixed_step <- prefixed_steps[[idx]]

        added_direct <- .pb_register_prone_step_name(
            step_name = direct_step,
            step_fun = step_fun,
            register_lowercase = register_lowercase
        )
        if (added_direct) {
            registered_direct <- c(registered_direct, direct_step)
        } else {
            collisions <- c(collisions, direct_step)
        }

        if (!identical(prefixed_step, direct_step)) {
            added_prefixed <- .pb_register_prone_step_name(
                step_name = prefixed_step,
                step_fun = step_fun,
                register_lowercase = register_lowercase
            )
            if (added_prefixed) {
                registered_prefixed <- c(registered_prefixed, prefixed_step)
            }
        }
    }

    invisible(list(
        direct = unname(registered_direct),
        prefixed = unname(registered_prefixed),
        collisions = unname(unique(collisions))
    ))
}

.pb_register_prone_step_name <- function(step_name,
                                         step_fun,
                                         register_lowercase = TRUE) {
    added <- .pb_register_step_if_absent(step_name, step_fun)
    if (added && isTRUE(register_lowercase)) {
        .pb_register_step_if_absent(tolower(step_name), step_fun)
    }
    added
}

.pb_register_step_if_absent <- function(step_name, step_fun) {
    if (!is.character(step_name) || length(step_name) != 1L || !nzchar(step_name)) {
        return(FALSE)
    }
    if (pb_has_step(step_name)) {
        return(FALSE)
    }
    pb_register_step(step_name, step_fun)
    TRUE
}

.pb_prone_norm_step_name <- function(norm_method, prefix = .pb_prone_default_step_prefix()) {
    token <- gsub("[^[:alnum:].]+", "", as.character(norm_method))
    if (!nzchar(token)) {
        token <- "method"
    }
    paste0(prefix, token)
}

.pb_prone_normalization_method_names <- function(methods = NULL) {
    if (is.null(methods)) {
        .pb_requireNamespace("PRONE")
        methods <- PRONE::get_normalization_methods()
    }
    methods <- as.character(methods)
    methods <- trimws(methods)
    methods <- methods[!is.na(methods)]
    methods <- unique(methods[nzchar(methods)])
    methods
}

.pb_make_prone_norm_step_fun <- function(norm_method) {
    force(norm_method)
    function(data_matrix,
             sample_annotation = NULL,
             sample_id_col = "FullRunName",
             on_raw = NULL,
             assay_in = "raw",
             batch = NULL,
             refs = NULL,
             ...) {
        .prone_normalize_matrix_step(
            data_matrix = data_matrix,
            sample_annotation = sample_annotation,
            sample_id_col = sample_id_col,
            on_raw = on_raw,
            norm_method = norm_method,
            assay_in = assay_in,
            batch = batch,
            refs = refs,
            ...
        )
    }
}

.pb_prone_normalize_single_fun <- function() {
    opt_fun <- getOption("proBatch.prone_normalize_se_single", NULL)
    if (is.function(opt_fun)) {
        return(opt_fun)
    }
    .pb_requireNamespace("PRONE")
    PRONE::normalize_se_single
}

.pb_prone_match_supported_args <- function(fun, base_args, dot_args = list()) {
    normalize_formals <- tryCatch(names(formals(fun)), error = function(e) NULL)
    if (is.null(normalize_formals)) {
        return(c(base_args, dot_args))
    }

    has_dots <- "..." %in% normalize_formals
    explicit_formals <- setdiff(normalize_formals, "...")
    if (length(explicit_formals)) {
        base_args <- base_args[intersect(names(base_args), explicit_formals)]
        if (!has_dots && length(dot_args)) {
            dot_args <- dot_args[intersect(names(dot_args), explicit_formals)]
        }
    } else if (!has_dots) {
        base_args <- list()
        dot_args <- list()
    }

    c(base_args, dot_args)
}

.pb_prone_resolve_sample_vector <- function(values,
                                            sample_annotation,
                                            sample_ids,
                                            arg_name) {
    if (is.null(values)) {
        return(NULL)
    }

    values_local <- values
    if (
        is.character(values_local) &&
            length(values_local) == 1L &&
            values_local %in% colnames(sample_annotation)
    ) {
        values_local <- sample_annotation[[values_local]]
    }

    if (!is.null(names(values_local))) {
        idx <- match(sample_ids, names(values_local))
        if (anyNA(idx)) {
            stop(
                "PRONE normalization: `", arg_name,
                "` has names that do not cover all samples."
            )
        }
        values_local <- values_local[idx]
    } else if (length(values_local) != length(sample_ids)) {
        stop(
            "PRONE normalization: `", arg_name,
            "` must have length 1 (column name) or length equal to number of samples."
        )
    }

    unname(values_local)
}

.pb_prone_resolve_refs <- function(refs, sample_annotation, sample_ids) {
    if (is.null(refs)) {
        return(NULL)
    }

    if (is.character(refs) && length(refs) && all(refs %in% sample_ids)) {
        return(unique(refs))
    }

    if (
        is.character(refs) &&
            length(refs) == 1L &&
            refs %in% colnames(sample_annotation)
    ) {
        refs_vec <- sample_annotation[[refs]]
    } else {
        refs_vec <- .pb_prone_resolve_sample_vector(
            values = refs,
            sample_annotation = sample_annotation,
            sample_ids = sample_ids,
            arg_name = "refs"
        )
    }

    if (is.logical(refs_vec) || is.numeric(refs_vec) || is.integer(refs_vec)) {
        idx <- as.logical(refs_vec)
        idx[is.na(idx)] <- FALSE
        return(sample_ids[idx])
    }

    refs_chr <- as.character(refs_vec)
    refs_chr <- refs_chr[!is.na(refs_chr) & nzchar(refs_chr)]
    refs_ids <- intersect(sample_ids, refs_chr)
    if (length(refs_ids)) {
        return(refs_ids)
    }

    ref_labels <- tolower(trimws(as.character(refs_vec)))
    ref_idx <- ref_labels %in% c(
        "ref", "reference", "pooled", "pool", "comref",
        "true", "yes", "1"
    )
    sample_ids[which(ref_idx)]
}

.pb_prone_attach_metadata <- function(se,
                                      sample_annotation,
                                      sample_ids,
                                      batch = NULL,
                                      refs = NULL) {
    batch_vec <- .pb_prone_resolve_sample_vector(
        values = batch,
        sample_annotation = sample_annotation,
        sample_ids = sample_ids,
        arg_name = "batch"
    )
    if (!is.null(batch_vec)) {
        S4Vectors::metadata(se)$batch <- as.vector(batch_vec)
    }

    refs_vec <- .pb_prone_resolve_refs(
        refs = refs,
        sample_annotation = sample_annotation,
        sample_ids = sample_ids
    )
    if (!is.null(refs_vec) && length(refs_vec)) {
        S4Vectors::metadata(se)$refs <- refs_vec
    }

    se
}

.pb_prone_guess_norm_assay <- function(se_norm, ain, norm_method) {
    assay_names <- SummarizedExperiment::assayNames(se_norm)
    if (!length(assay_names)) {
        stop("PRONE normalization returned no assays.")
    }

    candidates <- c(
        norm_method,
        paste0(norm_method, "_on_", ain),
        tolower(norm_method),
        paste0(tolower(norm_method), "_on_", ain)
    )
    matched <- candidates[candidates %in% assay_names]
    if (length(matched)) {
        return(matched[[1]])
    }

    added <- setdiff(assay_names, ain)
    if (length(added) == 1L) {
        return(added[[1]])
    }

    on_raw <- grep(paste0("_on_", ain, "$"), assay_names, value = TRUE)
    if (length(on_raw) == 1L) {
        return(on_raw[[1]])
    }

    if (length(assay_names) == 1L) {
        return(assay_names[[1]])
    }

    stop(
        "Unable to determine normalized assay returned by PRONE for method '",
        norm_method,
        "'. Available assays: ",
        paste(assay_names, collapse = ", ")
    )
}

.prone_normalize_matrix_step <- function(data_matrix,
                                         sample_annotation = NULL,
                                         sample_id_col = "FullRunName",
                                         on_raw = NULL,
                                         norm_method,
                                         assay_in = "raw",
                                         batch = NULL,
                                         refs = NULL,
                                         ...) {
    if (is.null(colnames(data_matrix))) {
        stop("PRONE normalization requires matrix column names (sample identifiers).")
    }

    sample_ids <- colnames(data_matrix)
    annotation_prep <- .pb_prone_prepare_sample_annotation(
        sample_annotation = sample_annotation,
        sample_ids = sample_ids,
        sample_id_col = sample_id_col,
        align = FALSE
    )
    sample_annotation_local <- annotation_prep$sample_annotation
    sample_id_col_local <- annotation_prep$sample_id_col
    placeholder_col <- annotation_prep$placeholder_col

    .run_matrix_method(
        data_matrix = data_matrix,
        sample_annotation = sample_annotation_local,
        sample_id_col = sample_id_col_local,
        fill_the_missing = NULL,
        missing_warning = "PRONE normalization requires sample identifiers.",
        method_fun = function(data_matrix, sample_annotation) {
            sample_ids_local <- colnames(data_matrix)
            sample_df <- as.data.frame(sample_annotation, stringsAsFactors = FALSE)

            if (!is.null(placeholder_col) && placeholder_col %in% names(sample_df)) {
                sample_df[[placeholder_col]] <- NULL
            }

            se <- .pb_prone_make_se(
                data_matrix = data_matrix,
                sample_annotation = sample_df,
                assay_in = assay_in
            )
            se <- .pb_prone_attach_metadata(
                se = se,
                sample_annotation = sample_df,
                sample_ids = sample_ids_local,
                batch = batch,
                refs = refs
            )

            normalize_single <- .pb_prone_normalize_single_fun()
            base_args <- list(
                se = se,
                methods = norm_method,
                method = norm_method,
                norm_method = norm_method,
                on_raw = on_raw,
                ain = assay_in,
                aout = norm_method
            )
            dot_args <- list(...)
            all_args <- .pb_prone_match_supported_args(
                fun = normalize_single,
                base_args = base_args,
                dot_args = dot_args
            )

            se_norm <- do.call(normalize_single, all_args)
            assay_out <- .pb_prone_guess_norm_assay(
                se_norm = se_norm,
                ain = assay_in,
                norm_method = norm_method
            )

            # PRONE often stores assays as data.table. With duplicated feature names
            # in rownames(se_norm), extracting with withDimnames=TRUE can fail because
            # data.frame-like containers require unique row names.
            normalized_mat <- SummarizedExperiment::assay(
                se_norm,
                assay_out,
                withDimnames = FALSE
            )
            if (!is.matrix(normalized_mat)) {
                normalized_mat <- as.matrix(normalized_mat)
            }
            storage.mode(normalized_mat) <- "double"

            .pb_prone_restore_dimnames(
                original_matrix = data_matrix,
                imputed_matrix = normalized_mat
            )
        }
    )
}
