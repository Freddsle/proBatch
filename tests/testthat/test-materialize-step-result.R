pb_test_materialize_fixture <- function() {
    input <- matrix(
        as.numeric(seq_len(9L)),
        nrow = 3L,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:3))
    )
    annotation <- data.frame(
        FullRunName = rev(colnames(input)),
        batch = c("B", "A", "A"),
        row.names = rev(colnames(input)),
        stringsAsFactors = FALSE
    )
    object <- suppressMessages(ProBatchFeatures(
        data_matrix = input,
        sample_annotation = annotation,
        sample_id_col = "FullRunName",
        name = "raw"
    ))

    state <- new.env(parent = emptyenv())
    state$calls <- 0L
    state$structured <- TRUE
    state$artifacts <- list(
        model = list(converged = TRUE, rank = 2L),
        scores = matrix(
            c(0.25, 0.75, 0.5),
            ncol = 1L,
            dimnames = list(colnames(input), "score")
        )
    )
    state$last <- NULL

    provider <- function(
        data_matrix,
        sample_annotation = NULL,
        multiplier = 2
    ) {
        state$calls <- state$calls + 1L
        if (
            !is.null(sample_annotation) &&
                !identical(rownames(sample_annotation), colnames(data_matrix))
        ) {
            stop("sample annotation was not aligned")
        }
        transformed <- data_matrix * multiplier
        result <- if (isTRUE(state$structured)) {
            pb_step_result(transformed, artifacts = state$artifacts)
        } else {
            transformed
        }
        state$last <- result
        result
    }

    list(
        input = input,
        annotation = annotation,
        object = object,
        source = "feature::raw",
        state = state,
        provider = provider
    )
}

pb_test_register_materialize_provider <- function(
    fixture,
    canonical,
    aliases,
    package = "testthat"
) {
    pb_register_step(
        name = canonical,
        fun = fixture$provider,
        package = package,
        kind = "transform",
        aliases = aliases,
        label = "Materialization test provider"
    )
}

pb_test_cleanup_materialize_providers <- function() {
    pb_unregister_steps("testthat")
    pb_unregister_steps("stats")
    invisible(TRUE)
}

test_that("precomputed alias results materialize without provider invocation", {
    pb_test_cleanup_materialize_providers()
    on.exit(pb_test_cleanup_materialize_providers(), add = TRUE)

    fixture <- pb_test_materialize_fixture()
    canonical <- "materialize_success_canonical"
    alias <- "materialize_success_alias"
    pb_test_register_materialize_provider(
        fixture,
        canonical = canonical,
        aliases = alias
    )

    raw <- fixture$object[[fixture$source]]
    raw_row_data <- SummarizedExperiment::rowData(raw)
    raw_row_data$sentinel <- paste0("row-", seq_len(nrow(raw)))
    SummarizedExperiment::rowData(raw) <- raw_row_data
    raw_metadata <- S4Vectors::metadata(raw)
    raw_metadata$sentinel <- list(scope = "raw-assay")
    S4Vectors::metadata(raw) <- raw_metadata
    fixture$object[[fixture$source]] <- raw

    object_metadata <- S4Vectors::metadata(fixture$object)
    object_metadata$sentinel <- list(scope = "object", value = 42L)
    S4Vectors::metadata(fixture$object) <- object_metadata
    object <- suppressMessages(pb_transform(
        fixture$object,
        from = fixture$source,
        steps = "log2",
        store_fast_steps = TRUE
    ))
    before <- object
    before_names <- names(before)
    before_log <- get_operation_log(before)
    before_metadata <- S4Vectors::metadata(before)

    params <- list(
        multiplier = 2L,
        option = "value",
        nested = list(flag = TRUE)
    )
    result <- fixture$provider(
        suppressMessages(pb_assay_matrix(object, fixture$source)),
        multiplier = params$multiplier
    )
    expect_identical(fixture$state$calls, 1L)

    updated <- pb_materialize_step_result(
        object = object,
        from = fixture$source,
        step = alias,
        result = result,
        params = params,
        backend = "memory"
    )
    target <- paste0("feature::", alias, "_on_raw")

    expect_identical(fixture$state$calls, 1L)
    expect_identical(object, before)
    expect_identical(names(updated)[seq_along(before_names)], before_names)
    for (assay_name in before_names) {
        expect_identical(updated[[assay_name]], before[[assay_name]])
        expect_identical(
            class(SummarizedExperiment::assay(
                updated[[assay_name]],
                "intensity"
            )),
            class(SummarizedExperiment::assay(
                before[[assay_name]],
                "intensity"
            ))
        )
    }
    expect_identical(
        SummarizedExperiment::colData(updated),
        SummarizedExperiment::colData(before)
    )
    for (metadata_name in names(before_metadata)) {
        expect_identical(
            S4Vectors::metadata(updated)[[metadata_name]],
            before_metadata[[metadata_name]]
        )
    }

    expect_true(target %in% names(updated))
    expect_true(is.matrix(SummarizedExperiment::assay(
        updated[[target]],
        "intensity"
    )))
    expect_identical(
        suppressMessages(pb_assay_matrix(updated, target)),
        result$data
    )
    expect_identical(pb_step_artifacts(updated, target), result$artifacts)
    expect_identical(
        SummarizedExperiment::colData(updated[[target]]),
        SummarizedExperiment::colData(before[[fixture$source]])
    )

    operation_log <- get_operation_log(updated)
    expect_identical(
        operation_log[seq_len(nrow(before_log)), , drop = FALSE],
        before_log
    )
    entry <- operation_log[nrow(operation_log), , drop = FALSE]
    expect_identical(as.character(entry$from[[1L]]), fixture$source)
    expect_identical(as.character(entry$to[[1L]]), target)
    expect_identical(as.character(entry$step[[1L]]), alias)
    expect_identical(as.character(entry$fun[[1L]]), canonical)
    expect_identical(as.character(entry$pkg[[1L]]), "testthat")
    expect_identical(entry$params[[1L]], params)
    expect_false("sample_annotation" %in% names(entry$params[[1L]]))
    parent_rows <- which(
        as.character(operation_log$to) == target &
            as.character(operation_log$from) != target
    )
    expect_length(parent_rows, 1L)
    expect_identical(get_chain(updated, assay = target), alias)

    explicit_target <- "feature::provider_on_raw"
    explicit <- pb_materialize_step_result(
        object = object,
        from = fixture$source,
        step = alias,
        result = result,
        params = params,
        final_name = explicit_target,
        backend = "memory"
    )
    expect_identical(fixture$state$calls, 1L)
    expect_true(explicit_target %in% names(explicit))
    expect_false(paste0(explicit_target, ".1") %in% names(explicit))
    expect_identical(
        suppressMessages(pb_assay_matrix(explicit, explicit_target)),
        result$data
    )
})

test_that("materialization resolves virtual sources and accepts ordered subsets", {
    pb_test_cleanup_materialize_providers()
    on.exit(pb_test_cleanup_materialize_providers(), add = TRUE)

    fixture <- pb_test_materialize_fixture()
    canonical <- "materialize_subset_canonical"
    alias <- "materialize_subset_alias"
    pb_test_register_materialize_provider(
        fixture,
        canonical = canonical,
        aliases = alias
    )

    virtual_source <- suppressMessages(pb_transform(
        fixture$object,
        from = fixture$source,
        steps = "log2"
    ))
    source <- "feature::log2_on_raw"
    expect_false(source %in% names(virtual_source))
    source_matrix <- suppressMessages(pb_assay_matrix(virtual_source, source))
    result <- fixture$provider(source_matrix, multiplier = 3)
    expect_identical(fixture$state$calls, 1L)

    updated <- pb_materialize_step_result(
        object = virtual_source,
        from = source,
        step = alias,
        result = result,
        params = list(multiplier = 3),
        backend = "memory"
    )
    target <- paste0("feature::", alias, "_on_log2_on_raw")
    expect_identical(fixture$state$calls, 1L)
    expect_true(target %in% names(updated))
    expect_identical(
        suppressMessages(pb_assay_matrix(updated, target)),
        result$data
    )

    subset_result <- pb_step_result(
        result$data[c("f1", "f3"), , drop = FALSE],
        artifacts = list(retained = c("f1", "f3"))
    )
    subsetted <- pb_materialize_step_result(
        object = virtual_source,
        from = source,
        step = alias,
        result = subset_result,
        params = list(multiplier = 3),
        final_name = "feature::ordered_subset",
        backend = "memory"
    )
    expect_identical(
        rownames(suppressMessages(pb_assay_matrix(
            subsetted,
            "feature::ordered_subset"
        ))),
        c("f1", "f3")
    )
    expect_identical(
        colnames(suppressMessages(pb_assay_matrix(
            subsetted,
            "feature::ordered_subset"
        ))),
        colnames(source_matrix)
    )
    expect_identical(
        pb_step_artifacts(subsetted, "feature::ordered_subset"),
        subset_result$artifacts
    )

    reordered_subset <- pb_step_result(
        result$data[c("f3", "f1"), , drop = FALSE]
    )
    expect_error(
        pb_materialize_step_result(
            object = virtual_source,
            from = source,
            step = alias,
            result = reordered_subset,
            final_name = "feature::reordered_subset"
        ),
        "preserve input feature and sample order"
    )
})

test_that("materialization rejects malformed results and incompatible axes", {
    pb_test_cleanup_materialize_providers()
    on.exit(pb_test_cleanup_materialize_providers(), add = TRUE)

    fixture <- pb_test_materialize_fixture()
    canonical <- "materialize_validation_canonical"
    alias <- "materialize_validation_alias"
    pb_test_register_materialize_provider(
        fixture,
        canonical = canonical,
        aliases = alias
    )

    materialize_data <- function(data) {
        pb_materialize_step_result(
            object = fixture$object,
            from = fixture$source,
            step = alias,
            result = pb_step_result(data),
            final_name = "feature::validation_target",
            backend = "memory"
        )
    }

    reordered_features <- fixture$input[c("f2", "f1", "f3"), , drop = FALSE]
    expect_error(
        materialize_data(reordered_features),
        "preserve input feature and sample order"
    )
    unknown_features <- fixture$input
    rownames(unknown_features)[[1L]] <- "unknown"
    expect_error(materialize_data(unknown_features), "unknown identifier")
    duplicate_features <- fixture$input
    rownames(duplicate_features)[[2L]] <- rownames(duplicate_features)[[1L]]
    expect_error(materialize_data(duplicate_features), "duplicate identifiers")
    missing_feature_names <- fixture$input
    rownames(missing_feature_names) <- NULL
    expect_error(materialize_data(missing_feature_names), "must be named")
    missing_feature_id <- fixture$input
    rownames(missing_feature_id)[[1L]] <- NA_character_
    expect_error(
        materialize_data(missing_feature_id),
        "NA or empty identifiers"
    )
    empty_feature_id <- fixture$input
    rownames(empty_feature_id)[[1L]] <- ""
    expect_error(materialize_data(empty_feature_id), "NA or empty identifiers")
    expect_error(
        materialize_data(fixture$input[FALSE, , drop = FALSE]),
        "at least one feature and one sample"
    )

    missing_sample <- fixture$input[, -1L, drop = FALSE]
    expect_error(
        materialize_data(missing_sample),
        "preserve every input sample"
    )
    extra_sample <- cbind(fixture$input, extra = c(10, 11, 12))
    expect_error(materialize_data(extra_sample), "unknown identifier")
    duplicate_samples <- fixture$input
    colnames(duplicate_samples)[[2L]] <- colnames(duplicate_samples)[[1L]]
    expect_error(materialize_data(duplicate_samples), "duplicate identifiers")
    reordered_samples <- fixture$input[, c("s2", "s1", "s3"), drop = FALSE]
    expect_error(
        materialize_data(reordered_samples),
        "preserve input feature and sample order"
    )
    missing_sample_names <- fixture$input
    colnames(missing_sample_names) <- NULL
    expect_error(materialize_data(missing_sample_names), "must be named")
    missing_sample_id <- fixture$input
    colnames(missing_sample_id)[[1L]] <- NA_character_
    expect_error(materialize_data(missing_sample_id), "NA or empty identifiers")
    empty_sample_id <- fixture$input
    colnames(empty_sample_id)[[1L]] <- ""
    expect_error(materialize_data(empty_sample_id), "NA or empty identifiers")

    expect_error(
        materialize_data(as.data.frame(fixture$input)),
        "numeric matrix"
    )
    character_matrix <- matrix(
        letters[seq_along(fixture$input)],
        nrow = nrow(fixture$input),
        dimnames = dimnames(fixture$input)
    )
    expect_error(materialize_data(character_matrix), "numeric matrix")
    latent <- matrix(
        seq_len(6L),
        nrow = ncol(fixture$input),
        dimnames = list(colnames(fixture$input), c("component1", "component2"))
    )
    expect_error(materialize_data(latent), "unknown identifier")
    dummy <- matrix(
        0,
        nrow = 1L,
        dimnames = list("dummy_feature", "dummy_sample")
    )
    expect_error(materialize_data(dummy), "unknown identifier")
    expect_error(materialize_data(2L), "numeric matrix")

    base_args <- list(
        object = fixture$object,
        from = fixture$source,
        step = alias,
        final_name = "feature::invalid_result",
        backend = "memory"
    )
    expect_error(
        do.call(
            pb_materialize_step_result,
            c(base_args, list(result = fixture$input))
        ),
        "must be a pb_step_result"
    )
    expect_error(
        do.call(
            pb_materialize_step_result,
            c(base_args, list(result = 2L))
        ),
        "must be a pb_step_result"
    )
    arbitrary_bench_result <- structure(
        list(data = fixture$input, artifacts = list()),
        class = "probatchbench_result"
    )
    expect_error(
        do.call(
            pb_materialize_step_result,
            c(base_args, list(result = arbitrary_bench_result))
        ),
        "must be a pb_step_result"
    )
    malformed_results <- list(
        structure(
            list(data = fixture$input),
            class = "pb_step_result"
        ),
        structure(
            list(data = fixture$input, artifacts = list(), extra = TRUE),
            class = "pb_step_result"
        ),
        structure(
            list(data = fixture$input, artifacts = "not-a-list"),
            class = "pb_step_result"
        )
    )
    for (malformed in malformed_results) {
        expect_error(
            do.call(
                pb_materialize_step_result,
                c(base_args, list(result = malformed))
            ),
            "Invalid `pb_step_result`"
        )
    }
})

test_that("materialization exposes no identity override or validation bypass", {
    pb_test_cleanup_materialize_providers()
    on.exit(pb_test_cleanup_materialize_providers(), add = TRUE)

    fixture <- pb_test_materialize_fixture()
    canonical <- "materialize_boundary_canonical"
    alias <- "materialize_boundary_alias"
    pb_test_register_materialize_provider(
        fixture,
        canonical = canonical,
        aliases = alias
    )
    result <- pb_step_result(fixture$input, artifacts = list())

    expect_identical(
        formals(pb_materialize_step_result),
        as.pairlist(alist(
            object = ,
            from = ,
            step = ,
            result = ,
            params = list(),
            final_name = NULL,
            backend = c("auto", "memory", "hdf5"),
            hdf5_path = NULL
        ))
    )
    expect_identical(
        formals(pb_transform),
        as.pairlist(alist(
            object = ,
            from = ,
            steps = ,
            funs = NULL,
            params_list = NULL,
            level = NULL,
            store_fast_steps = FALSE,
            fast_steps = c("log", "log2", "medianNorm"),
            store_intermediate = FALSE,
            final_name = NULL,
            backend = c("auto", "memory", "hdf5"),
            hdf5_path = NULL
        ))
    )
    expect_identical(
        formals(pb_eval),
        as.pairlist(alist(
            object = ,
            from = ,
            steps = ,
            funs = NULL,
            params_list = NULL
        ))
    )

    expect_error(
        pb_materialize_step_result(
            fixture$input,
            fixture$source,
            alias,
            result
        ),
        "must be a ProBatchFeatures object"
    )
    expect_error(
        pb_materialize_step_result(
            fixture$object,
            "feature::unknown",
            alias,
            result
        ),
        "not found in object or operation log"
    )
    expect_error(
        pb_materialize_step_result(
            fixture$object,
            fixture$source,
            "not_registered",
            result
        ),
        "not registered"
    )
    for (bad_from in list("", NA_character_, c("feature::raw", "other"))) {
        expect_error(
            pb_materialize_step_result(
                fixture$object,
                bad_from,
                alias,
                result
            ),
            "`from`"
        )
    }
    for (bad_step in list("", NA_character_)) {
        expect_error(
            pb_materialize_step_result(
                fixture$object,
                fixture$source,
                bad_step,
                result
            ),
            "`step`"
        )
    }
    expect_error(
        pb_materialize_step_result(
            fixture$object,
            fixture$source,
            fixture$provider,
            result
        ),
        "`step`"
    )
    expect_error(
        pb_materialize_step_result(
            fixture$object,
            fixture$source,
            c(alias, canonical),
            result
        ),
        "`step`"
    )

    for (bad_params in list(NULL, 1, identity)) {
        expect_error(
            pb_materialize_step_result(
                fixture$object,
                fixture$source,
                alias,
                result,
                params = bad_params
            ),
            "`params` must be a list"
        )
    }
    for (bad_name in list("", NA_character_, c("a", "b"), "malformed")) {
        expect_error(
            pb_materialize_step_result(
                fixture$object,
                fixture$source,
                alias,
                result,
                final_name = bad_name
            ),
            "Assay name"
        )
    }
    expect_error(
        pb_materialize_step_result(
            fixture$object,
            fixture$source,
            alias,
            result,
            final_name = fixture$source
        ),
        "conflicts with the source assay"
    )

    call_args <- list(
        object = fixture$object,
        from = fixture$source,
        step = alias,
        result = result,
        final_name = "feature::no_override"
    )
    for (override in list(
        list(fun = fixture$provider),
        list(package = "stats"),
        list(pkg = "stats"),
        list(provider_callback = fixture$provider),
        list(validate = FALSE)
    )) {
        expect_error(
            do.call(pb_materialize_step_result, c(call_args, override)),
            "unused argument"
        )
    }

    unavailable_owner <- "proBatchUnavailableProviderForMaterializationTest"
    unavailable_step <- "materialize_unavailable_canonical"
    on.exit(pb_unregister_steps(unavailable_owner), add = TRUE)
    pb_register_step(
        unavailable_step,
        fixture$provider,
        package = unavailable_owner
    )
    expect_false(pb_has_step(unavailable_step, available = TRUE))
    unnamed_artifacts <- unname(list(TRUE, 2L))
    unavailable <- pb_materialize_step_result(
        fixture$object,
        fixture$source,
        unavailable_step,
        pb_step_result(fixture$input, unnamed_artifacts),
        final_name = "feature::unavailable_provider_result",
        backend = "memory"
    )
    unavailable_log <- get_operation_log(unavailable)
    unavailable_entry <- unavailable_log[nrow(unavailable_log), , drop = FALSE]
    expect_identical(as.character(unavailable_entry$step), unavailable_step)
    expect_identical(as.character(unavailable_entry$fun), unavailable_step)
    expect_identical(
        as.character(unavailable_entry$pkg),
        unavailable_owner
    )
    expect_identical(
        pb_step_artifacts(unavailable, "feature::unavailable_provider_result"),
        unnamed_artifacts
    )
    expect_null(names(pb_step_artifacts(
        unavailable,
        "feature::unavailable_provider_result"
    )))
    expect_identical(fixture$state$calls, 0L)
})

test_that("stored retries compare every materialization identity component", {
    pb_test_cleanup_materialize_providers()
    on.exit(pb_test_cleanup_materialize_providers(), add = TRUE)

    fixture <- pb_test_materialize_fixture()
    canonical <- "materialize_retry_canonical"
    alternate_canonical <- "materialize_retry_other_canonical"
    alias <- "materialize_retry_alias"
    alternate_alias <- "materialize_retry_alternate_alias"
    pb_test_register_materialize_provider(
        fixture,
        canonical = canonical,
        aliases = c(alias, alternate_alias)
    )

    params <- list(multiplier = 2L, option = "strict")
    result <- fixture$provider(fixture$input, multiplier = params$multiplier)
    target <- "feature::retry_target"
    stored <- pb_materialize_step_result(
        fixture$object,
        fixture$source,
        alias,
        result,
        params = params,
        final_name = target,
        backend = "memory"
    )
    stored_log <- get_operation_log(stored)
    expect_identical(fixture$state$calls, 1L)

    exact <- pb_materialize_step_result(
        stored,
        fixture$source,
        alias,
        result,
        params = params,
        final_name = target,
        backend = "memory"
    )
    expect_identical(exact, stored)
    expect_identical(get_operation_log(exact)$timestamp, stored_log$timestamp)
    expect_identical(fixture$state$calls, 1L)

    retry_path <- tempfile(fileext = ".h5")
    on.exit(unlink(retry_path), add = TRUE)
    alternate_backend <- pb_materialize_step_result(
        stored,
        fixture$source,
        alias,
        result,
        params = params,
        final_name = target,
        backend = "hdf5",
        hdf5_path = retry_path
    )
    expect_identical(alternate_backend, stored)
    expect_false(file.exists(retry_path))

    changed_matrix <- result
    changed_matrix$data[[1L]] <- changed_matrix$data[[1L]] + 1
    expect_error(
        pb_materialize_step_result(
            stored,
            fixture$source,
            alias,
            changed_matrix,
            params = params,
            final_name = target
        ),
        "conflicting data"
    )
    changed_artifacts <- result
    changed_artifacts$artifacts$model$rank <- 3L
    expect_error(
        pb_materialize_step_result(
            stored,
            fixture$source,
            alias,
            changed_artifacts,
            params = params,
            final_name = target
        ),
        "conflicting structured artifacts"
    )
    expect_error(
        pb_materialize_step_result(
            stored,
            fixture$source,
            alias,
            result,
            params = list(multiplier = 2L, option = "relaxed"),
            final_name = target
        ),
        "different parent or stable lineage origin"
    )

    with_alternate_parent <- QFeatures::addAssay(
        stored,
        stored[[fixture$source]],
        name = "feature::alternate_parent"
    )
    expect_error(
        pb_materialize_step_result(
            with_alternate_parent,
            "feature::alternate_parent",
            alias,
            result,
            params = params,
            final_name = target
        ),
        "different parent or stable lineage origin"
    )
    expect_error(
        pb_materialize_step_result(
            stored,
            fixture$source,
            alternate_alias,
            result,
            params = params,
            final_name = target
        ),
        "different parent or stable lineage origin"
    )

    occupied <- QFeatures::addAssay(
        fixture$object,
        fixture$object[[fixture$source]],
        name = "feature::occupied_target"
    )
    expect_error(
        pb_materialize_step_result(
            occupied,
            fixture$source,
            alias,
            result,
            params = params,
            final_name = "feature::occupied_target"
        ),
        "already exists or is reserved"
    )
    virtual_reserved <- suppressMessages(pb_transform(
        fixture$object,
        from = fixture$source,
        steps = "log2"
    ))
    expect_error(
        pb_materialize_step_result(
            virtual_reserved,
            fixture$source,
            alias,
            result,
            params = params,
            final_name = "feature::log2_on_raw"
        ),
        "different parent or stable lineage origin"
    )

    pb_unregister_steps("testthat")
    pb_test_register_materialize_provider(
        fixture,
        canonical = alternate_canonical,
        aliases = alias
    )
    expect_error(
        pb_materialize_step_result(
            stored,
            fixture$source,
            alias,
            result,
            params = params,
            final_name = target
        ),
        "different parent or stable lineage origin"
    )

    pb_unregister_steps("testthat")
    pb_test_register_materialize_provider(
        fixture,
        canonical = canonical,
        aliases = alias,
        package = "stats"
    )
    expect_error(
        pb_materialize_step_result(
            stored,
            fixture$source,
            alias,
            result,
            params = params,
            final_name = target
        ),
        "different parent or stable lineage origin"
    )
    expect_identical(fixture$state$calls, 1L)
    expect_identical(pb_step_artifacts(stored, target), result$artifacts)
    expect_identical(get_operation_log(stored), stored_log)
})

test_that("matching virtual targets materialize without replay", {
    pb_test_cleanup_materialize_providers()
    on.exit(pb_test_cleanup_materialize_providers(), add = TRUE)

    fixture <- pb_test_materialize_fixture()
    canonical <- "materialize_virtual_canonical"
    alias <- "materialize_virtual_alias"
    pb_test_register_materialize_provider(
        fixture,
        canonical = canonical,
        aliases = alias
    )
    params <- list(multiplier = 2L)

    virtual <- suppressMessages(pb_transform(
        fixture$object,
        from = fixture$source,
        steps = c(alias, "log2"),
        params_list = list(params, list()),
        fast_steps = c(alias, "log2"),
        store_fast_steps = FALSE
    ))
    operation_log <- get_operation_log(virtual)
    target <- as.character(operation_log$to[[1L]])
    expect_identical(fixture$state$calls, 1L)
    expect_false(target %in% names(virtual))

    expect_error(
        pb_step_artifacts(virtual, target),
        "virtual and has not been materialized"
    )
    expect_identical(fixture$state$calls, 1L)

    virtual_with_self <- proBatch:::.pb_add_log_entry(
        virtual,
        step = "zeroIsNA",
        fun = "zeroIsNA",
        from = target,
        to = target,
        params = list(),
        pkg = "proBatch"
    )
    expect_error(
        pb_materialize_step_result(
            virtual_with_self,
            fixture$source,
            alias,
            fixture$state$last,
            params = params,
            final_name = target,
            backend = "memory"
        ),
        "includes additional in-place operations"
    )
    expect_identical(fixture$state$calls, 1L)

    materialized <- pb_materialize_step_result(
        virtual,
        fixture$source,
        alias,
        fixture$state$last,
        params = params,
        final_name = target,
        backend = "memory"
    )
    expect_identical(fixture$state$calls, 1L)
    expect_true(target %in% names(materialized))
    expect_identical(get_operation_log(materialized), operation_log)
    expect_identical(
        suppressMessages(pb_assay_matrix(materialized, target)),
        fixture$state$last$data
    )
    expect_identical(
        pb_step_artifacts(materialized, target),
        fixture$state$last$artifacts
    )
    expect_identical(get_chain(materialized, assay = target), alias)
})

test_that("structured transform retries are artifact-aware only when needed", {
    pb_test_cleanup_materialize_providers()
    on.exit(pb_test_cleanup_materialize_providers(), add = TRUE)

    fixture <- pb_test_materialize_fixture()
    canonical <- "materialize_transform_retry_canonical"
    alias <- "materialize_transform_retry_alias"
    pb_test_register_materialize_provider(
        fixture,
        canonical = canonical,
        aliases = alias
    )
    params <- list(multiplier = 2L)
    target <- "feature::structured_transform_retry"

    stored <- pb_transform(
        fixture$object,
        from = fixture$source,
        steps = alias,
        params_list = list(params),
        final_name = target,
        backend = "memory"
    )
    stored_log <- get_operation_log(stored)
    exact <- pb_transform(
        stored,
        from = fixture$source,
        steps = alias,
        params_list = list(params),
        final_name = target,
        backend = "memory"
    )
    expect_identical(exact, stored)

    fixture$state$artifacts$model$rank <- 3L
    expect_error(
        pb_transform(
            stored,
            from = fixture$source,
            steps = alias,
            params_list = list(params),
            final_name = target
        ),
        "conflicting structured artifacts"
    )
    fixture$state$artifacts$model$rank <- 2L

    missing_metadata <- stored
    missing_assay <- missing_metadata[[target]]
    missing_assay_metadata <- S4Vectors::metadata(missing_assay)
    missing_assay_metadata[[proBatch:::.pb_step_artifacts_key]] <- NULL
    S4Vectors::metadata(missing_assay) <- missing_assay_metadata
    missing_metadata[[target]] <- missing_assay
    expect_error(
        pb_transform(
            missing_metadata,
            from = fixture$source,
            steps = alias,
            params_list = list(params),
            final_name = target
        ),
        "conflicting structured artifacts"
    )

    malformed <- stored
    malformed_assay <- malformed[[target]]
    malformed_metadata <- S4Vectors::metadata(malformed_assay)
    malformed_metadata[[proBatch:::.pb_step_artifacts_key]] <- "malformed"
    S4Vectors::metadata(malformed_assay) <- malformed_metadata
    malformed[[target]] <- malformed_assay
    expect_error(
        pb_transform(
            malformed,
            from = fixture$source,
            steps = alias,
            params_list = list(params),
            final_name = target
        ),
        "malformed structured artifact metadata"
    )

    fixture$state$structured <- FALSE
    plain_retry <- pb_transform(
        malformed,
        from = fixture$source,
        steps = alias,
        params_list = list(params),
        final_name = target
    )
    expect_identical(plain_retry, malformed)
    expect_identical(get_operation_log(plain_retry), stored_log)
    expect_identical(
        suppressMessages(pb_assay_matrix(plain_retry, target)),
        fixture$input * 2
    )
})

test_that("artifact access is stored-only, opaque, and serializable", {
    pb_test_cleanup_materialize_providers()
    on.exit(pb_test_cleanup_materialize_providers(), add = TRUE)

    fixture <- pb_test_materialize_fixture()
    canonical <- "materialize_artifact_canonical"
    alias <- "materialize_artifact_alias"
    pb_test_register_materialize_provider(
        fixture,
        canonical = canonical,
        aliases = alias
    )
    result <- fixture$provider(fixture$input)
    target <- "feature::artifact_target"
    updated <- pb_materialize_step_result(
        fixture$object,
        fixture$source,
        alias,
        result,
        final_name = target,
        backend = "memory"
    )

    expect_identical(pb_current_assay(updated), target)
    expect_identical(pb_step_artifacts(updated), result$artifacts)
    expect_identical(pb_step_artifacts(updated, target), result$artifacts)
    expect_identical(pb_step_artifacts(updated, fixture$source), list())
    expect_error(pb_step_artifacts(fixture$input), "ProBatchFeatures object")
    expect_error(pb_step_artifacts(updated, "feature::unknown"), "not stored")
    for (bad_assay in list("", NA_character_, c(target, fixture$source))) {
        expect_error(pb_step_artifacts(updated, bad_assay), "`assay` must be")
    }

    malformed <- updated
    assay <- malformed[[target]]
    assay_metadata <- S4Vectors::metadata(assay)
    assay_metadata[[proBatch:::.pb_step_artifacts_key]] <- 1
    S4Vectors::metadata(assay) <- assay_metadata
    malformed[[target]] <- assay
    expect_error(
        pb_step_artifacts(malformed, target),
        "malformed structured artifact metadata"
    )

    rds_path <- tempfile(fileext = ".rds")
    on.exit(unlink(rds_path), add = TRUE)
    saveRDS(updated, rds_path)
    restored <- readRDS(rds_path)
    expect_identical(pb_step_artifacts(restored, target), result$artifacts)
    expect_identical(
        suppressMessages(pb_assay_matrix(restored, target)),
        result$data
    )
    expect_identical(get_operation_log(restored), get_operation_log(updated))
})

test_that("precomputed results support the optional HDF5 backend", {
    skip_if_not_installed("HDF5Array")
    pb_test_cleanup_materialize_providers()
    on.exit(pb_test_cleanup_materialize_providers(), add = TRUE)

    fixture <- pb_test_materialize_fixture()
    canonical <- "materialize_hdf5_canonical"
    alias <- "materialize_hdf5_alias"
    pb_test_register_materialize_provider(
        fixture,
        canonical = canonical,
        aliases = alias
    )
    result <- fixture$provider(fixture$input)
    path <- tempfile(fileext = ".h5")
    on.exit(unlink(path), add = TRUE)
    target <- "feature::hdf5_target"

    updated <- pb_materialize_step_result(
        fixture$object,
        fixture$source,
        alias,
        result,
        final_name = target,
        backend = "hdf5",
        hdf5_path = path
    )
    stored <- SummarizedExperiment::assay(updated[[target]], "intensity")
    expect_false(is.matrix(stored))
    expect_true(
        methods::is(stored, "HDF5Matrix") ||
            methods::is(stored, "DelayedMatrix")
    )
    expect_identical(as.matrix(stored), result$data)
    expect_identical(pb_step_artifacts(updated, target), result$artifacts)
    expect_identical(fixture$state$calls, 1L)

    source_path <- tempfile(fileext = ".h5")
    on.exit(unlink(source_path), add = TRUE)
    source_matrix <- HDF5Array::writeHDF5Array(
        fixture$input,
        filepath = source_path
    )
    source_col_data <- S4Vectors::DataFrame(
        sample_id = colnames(fixture$input),
        row.names = colnames(fixture$input)
    )
    source_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(intensity = source_matrix),
        colData = source_col_data
    )
    source_qf <- QFeatures::QFeatures(
        experiments = list(raw = source_se),
        colData = source_col_data
    )
    hdf5_source <- suppressMessages(as_ProBatchFeatures(source_qf))
    source_name <- "feature::raw"
    source_before <- hdf5_source[[source_name]]
    source_storage_before <- SummarizedExperiment::assay(
        source_before,
        "intensity"
    )
    source_result <- fixture$provider(as.matrix(source_storage_before))
    preserved <- pb_materialize_step_result(
        hdf5_source,
        source_name,
        alias,
        source_result,
        final_name = "feature::memory_target",
        backend = "memory"
    )
    source_storage_after <- SummarizedExperiment::assay(
        preserved[[source_name]],
        "intensity"
    )
    expect_identical(preserved[[source_name]], source_before)
    expect_identical(class(source_storage_after), class(source_storage_before))
    expect_false(is.matrix(source_storage_after))
    expect_identical(as.matrix(source_storage_after), fixture$input)
    expect_identical(fixture$state$calls, 2L)
})
