local_mocked_prone_normalize <- function(mock_fun) {
    caller_env <- parent.frame()
    prev_opt <- getOption("proBatch.prone_normalize_se_single", NULL)
    withr::defer(
        options(proBatch.prone_normalize_se_single = prev_opt),
        envir = caller_env
    )
    options(proBatch.prone_normalize_se_single = mock_fun)
}

make_prone_norm_test_pbf <- function() {
    dm <- matrix(
        c(
            10, 11, 13, 14, 16, 18,
            20, 19, 21, 22, 24, 23,
            30, 28, 27, 29, 31, 33,
            40, 41, 39, 38, 37, 36,
            50, 49, 48, 47, 46, 45,
            60, 59, 61, 62, 64, 63
        ),
        nrow = 6,
        byrow = TRUE,
        dimnames = list(
            paste0("feat", 1:6),
            paste0("sample", 1:6)
        )
    )

    sa <- data.frame(
        FullRunName = colnames(dm),
        MS_batch = rep(c("B1", "B2"), each = 3),
        stringsAsFactors = FALSE
    )

    ProBatchFeatures(
        data_matrix = dm,
        sample_annotation = sa,
        sample_id_col = "FullRunName",
        name = "raw"
    )
}

make_prone_norm_compare_test_pbf <- function() {
    n_features <- 80L
    n_samples <- 12L
    feature_index <- seq_len(n_features)
    sample_index <- seq_len(n_samples)

    dm <- outer(feature_index, sample_index, function(i, j) {
        base <- 10000 + (i * 73) + (j * 191)
        interaction <- ((i %% 11) * (j %% 5) * 7) + ((i %% 3) * 17)
        as.numeric(base + interaction)
    })
    rownames(dm) <- paste0("feat", feature_index)
    colnames(dm) <- paste0("sample", sample_index)

    sa <- data.frame(
        FullRunName = colnames(dm),
        MS_batch = rep(c("B1", "B2"), each = n_samples / 2),
        stringsAsFactors = FALSE
    )

    ProBatchFeatures(
        data_matrix = dm,
        sample_annotation = sa,
        sample_id_col = "FullRunName",
        name = "raw"
    )
}

.make_direct_prone_se <- function(data_matrix, sample_annotation) {
    sample_df <- as.data.frame(sample_annotation, stringsAsFactors = FALSE)
    if (!("FullRunName" %in% colnames(sample_df))) {
        sample_df$FullRunName <- colnames(data_matrix)
    }
    rownames(sample_df) <- sample_df$FullRunName
    col_data <- S4Vectors::DataFrame(sample_df)
    rownames(col_data) <- sample_df$FullRunName

    feature_ids <- rownames(data_matrix)
    row_data <- S4Vectors::DataFrame(
        "Protein.IDs" = feature_ids,
        "IDs" = feature_ids,
        check.names = FALSE
    )

    SummarizedExperiment::SummarizedExperiment(
        assays = list(raw = data_matrix),
        colData = col_data,
        rowData = row_data
    )
}

.to_numeric_matrix <- function(x) {
    rn <- rownames(x)
    cn <- colnames(x)
    mat <- if (is.matrix(x)) x else as.matrix(x)
    if (!is.null(rn) && length(rn) == nrow(mat)) {
        rownames(mat) <- rn
    }
    if (!is.null(cn) && length(cn) == ncol(mat)) {
        colnames(mat) <- cn
    }
    storage.mode(mat) <- "double"
    mat
}

.run_direct_prone_normalization <- function(se,
                                            method,
                                            on_raw = NULL,
                                            assay_in = "raw",
                                            ...) {
    normalize_fun <- .pb_prone_normalize_single_fun()
    base_args <- .pb_prone_normalization_base_args(
        se = se,
        norm_method = method,
        on_raw = on_raw,
        assay_in = assay_in
    )
    all_args <- .pb_prone_match_supported_args(
        fun = normalize_fun,
        base_args = base_args,
        dot_args = list(...)
    )
    do.call(normalize_fun, all_args)
}

.extract_direct_prone_matrix <- function(se_norm,
                                         input_matrix,
                                         method,
                                         assay_in = "raw") {
    assay_out <- .pb_prone_guess_norm_assay(
        se_norm = se_norm,
        ain = assay_in,
        norm_method = method
    )
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
        original_matrix = input_matrix,
        imputed_matrix = normalized_mat
    )
}

.resolve_prone_step_name <- function(method) {
    prefixed <- .pb_prone_norm_step_name(method)
    candidates <- c(
        prefixed,
        tolower(prefixed)
    )
    steps <- pb_list_steps()
    matched <- candidates[candidates %in% steps]
    if (!length(matched)) {
        stop("PRONE step for method '", method, "' is not registered.")
    }
    matched[[1]]
}

test_that("Direct PRONE normalization matches proBatch PRONE steps for complex methods", {
    skip_if_not_installed("PRONE")

    pbf <- make_prone_norm_compare_test_pbf()
    methods <- c("VSN", "NormicsVSN", "NormicsMedian")
    .pb_register_prone_normalization_steps(
        methods = methods,
        prefix = "PR",
        register_lowercase = TRUE
    )

    raw_mat <- pb_assay_matrix(pbf, "feature::raw")
    sample_ann <- as.data.frame(SummarizedExperiment::colData(pbf), stringsAsFactors = FALSE)

    method_params <- list(
        VSN = list(VSN_quantile = 0.9),
        NormicsVSN = list(
            reduce_correlation_by = 1,
            NormicsVSN_quantile = 0.8,
            top_x = 60
        ),
        NormicsMedian = list(
            reduce_correlation_by = 1,
            NormicsVSN_quantile = 0.8,
            top_x = 60
        )
    )

    for (method in methods) {
        se_direct <- .make_direct_prone_se(
            data_matrix = raw_mat,
            sample_annotation = sample_ann
        )
        direct_args <- c(
            list(
                se = se_direct,
                method = method,
                on_raw = NULL,
                assay_in = "raw"
            ),
            method_params[[method]]
        )
        se_direct_norm <- do.call(.run_direct_prone_normalization, direct_args)
        direct_mat <- .extract_direct_prone_matrix(
            se_norm = se_direct_norm,
            input_matrix = raw_mat,
            method = method,
            assay_in = "raw"
        )

        step_name <- .resolve_prone_step_name(method)
        pb_args <- c(
            list(sample_id_col = "FullRunName"),
            method_params[[method]]
        )
        pb_mat <- .to_numeric_matrix(pb_eval(
            object = pbf,
            from = "feature::raw",
            steps = step_name,
            params_list = list(pb_args)
        ))

        expect_identical(dim(pb_mat), dim(direct_mat), info = method)
        expect_identical(rownames(pb_mat), rownames(direct_mat), info = method)
        expect_identical(colnames(pb_mat), colnames(direct_mat), info = method)
        expect_equal(pb_mat, direct_mat, tolerance = 1e-8, info = method)
    }
})

test_that("proBatch limmaRBE matches PRONE limBE after harmonizing input scale", {
    skip_if_not_installed("PRONE")

    pbf <- make_prone_norm_compare_test_pbf()
    .pb_register_prone_normalization_steps(
        methods = "limBE",
        prefix = "PR",
        register_lowercase = TRUE
    )
    prone_step <- .resolve_prone_step_name("limBE")

    prone_limbe <- suppressWarnings(.to_numeric_matrix(pb_eval(
        object = pbf,
        from = "feature::raw",
        steps = c("log2", prone_step),
        params_list = list(
            list(),
            list(
                sample_id_col = "FullRunName",
                batch = "MS_batch",
                # limBE expects log2 input by default; after a preceding log2 step,
                # point PRONE to the "log2" assay to avoid an extra log2 transform.
                assay_in = "log2",
                on_raw = FALSE
            )
        )
    )))

    probatch_limmarbe <- suppressWarnings(.to_numeric_matrix(pb_eval(
        object = pbf,
        from = "feature::raw",
        steps = c("log2", "limmaRBE"),
        params_list = list(
            list(),
            list(
                sample_id_col = "FullRunName",
                batch_col = "MS_batch"
            )
        )
    )))

    expect_identical(dim(probatch_limmarbe), dim(prone_limbe))
    expect_identical(rownames(probatch_limmarbe), rownames(prone_limbe))
    expect_identical(colnames(probatch_limmarbe), colnames(prone_limbe))

    max_abs_diff <- max(abs(probatch_limmarbe - prone_limbe), na.rm = TRUE)
    expect_lte(max_abs_diff, 1e-8)
    expect_equal(probatch_limmarbe, prone_limbe, tolerance = 1e-8)
})

test_that("Direct PRONE sequential normalization matches proBatch PRONE step chains", {
    skip_if_not_installed("PRONE")

    pbf <- make_prone_norm_compare_test_pbf()
    chain_methods <- c("NormicsMedian", "VSN")
    .pb_register_prone_normalization_steps(
        methods = chain_methods,
        prefix = "PR",
        register_lowercase = TRUE
    )

    raw_mat <- pb_assay_matrix(pbf, "feature::raw")
    sample_ann <- as.data.frame(SummarizedExperiment::colData(pbf), stringsAsFactors = FALSE)

    method_params <- list(
        NormicsMedian = list(
            reduce_correlation_by = 1,
            NormicsVSN_quantile = 0.8,
            top_x = 60
        ),
        VSN = list(VSN_quantile = 0.9)
    )

    direct_mat <- .to_numeric_matrix(raw_mat)
    for (method in chain_methods) {
        se_direct <- .make_direct_prone_se(
            data_matrix = direct_mat,
            sample_annotation = sample_ann
        )
        direct_args <- c(
            list(
                se = se_direct,
                method = method,
                on_raw = NULL,
                assay_in = "raw"
            ),
            method_params[[method]]
        )
        se_direct_norm <- do.call(.run_direct_prone_normalization, direct_args)
        direct_mat <- .extract_direct_prone_matrix(
            se_norm = se_direct_norm,
            input_matrix = direct_mat,
            method = method,
            assay_in = "raw"
        )
    }

    pb_steps <- vapply(chain_methods, .resolve_prone_step_name, FUN.VALUE = character(1))
    pb_params <- lapply(
        chain_methods,
        function(method) {
            c(
                list(sample_id_col = "FullRunName"),
                method_params[[method]]
            )
        }
    )
    pb_chain_mat <- .to_numeric_matrix(pb_eval(
        object = pbf,
        from = "feature::raw",
        steps = pb_steps,
        params_list = pb_params
    ))

    expect_identical(dim(pb_chain_mat), dim(direct_mat))
    expect_identical(rownames(pb_chain_mat), rownames(direct_mat))
    expect_identical(colnames(pb_chain_mat), colnames(direct_mat))
    expect_equal(pb_chain_mat, direct_mat, tolerance = 1e-8)
})

test_that("internal PRONE registration keeps collisions on proBatch and adds PR* aliases", {
    unique_method <- "pbUniqueNormMethod"
    info <- .pb_register_prone_normalization_steps(
        methods = c(unique_method, "medianNorm"),
        prefix = "PR",
        register_lowercase = TRUE
    )

    expect_type(info, "list")
    expect_true(unique_method %in% c(info$direct, info$collisions))
    expect_true("medianNorm" %in% info$collisions)
    expect_true("PRmedianNorm" %in% pb_list_steps())
    expect_true(paste0("PR", unique_method) %in% pb_list_steps())
    expect_true(tolower(paste0("PR", unique_method)) %in% pb_list_steps())
    expect_true(tolower(unique_method) %in% pb_list_steps())
    expect_true("medianNorm" %in% pb_list_steps())
})

test_that("pb_transform chains mocked PRONE normalization with combat", {
    pbf <- make_prone_norm_test_pbf()
    reg_info <- .pb_register_prone_normalization_steps(
        methods = "mockMean",
        prefix = "PR",
        register_lowercase = FALSE
    )
    step_name <- if (length(reg_info$prefixed)) reg_info$prefixed[[1]] else "PRmockMean"

    local_mocked_prone_normalize(function(se, ain, norm_method, ...) {
        mat <- SummarizedExperiment::assay(se, ain)
        SummarizedExperiment::assay(
            se,
            paste0(norm_method, "_on_", ain),
            withDimnames = FALSE
        ) <- mat + 1
        se
    })

    out <- suppressWarnings(
        pb_transform(
            object = pbf,
            from = "feature::raw",
            steps = c(step_name, "combat"),
            params_list = list(
                list(sample_id_col = "FullRunName"),
                list(
                    batch_col = "MS_batch",
                    sample_id_col = "FullRunName",
                    par.prior = TRUE
                )
            )
        )
    )

    norm_assay <- paste0("feature::", step_name, "_on_raw")
    combat_assay <- paste0("feature::combat_on_", step_name, "_on_raw")
    expect_true(norm_assay %in% names(out))
    expect_true(combat_assay %in% names(out))

    raw_mat <- pb_assay_matrix(pbf, "feature::raw")
    norm_mat <- pb_assay_matrix(out, norm_assay)
    combat_mat <- pb_assay_matrix(out, combat_assay)

    expect_equal(norm_mat, raw_mat + 1, ignore_attr = TRUE)
    expect_equal(dim(combat_mat), dim(raw_mat))
})

test_that("PRONE normalization methods run systematically with BEC steps in pb_eval", {
    pbf <- make_prone_norm_test_pbf()
    reg_info <- .pb_register_prone_normalization_steps(
        methods = c("mockA", "mockB"),
        prefix = "PR",
        register_lowercase = FALSE
    )
    step_names <- if (length(reg_info$prefixed)) reg_info$prefixed else c("PRmockA", "PRmockB")
    shift_map <- c(mockA = 0.5, mockB = 1.5)

    local_mocked_prone_normalize(function(se, ain, norm_method, ...) {
        mat <- SummarizedExperiment::assay(se, ain)
        shift <- unname(shift_map[[norm_method]])
        SummarizedExperiment::assay(
            se,
            paste0(norm_method, "_on_", ain),
            withDimnames = FALSE
        ) <- mat + shift
        se
    })

    param_by_bec <- list(
        combat = list(batch_col = "MS_batch", sample_id_col = "FullRunName", par.prior = TRUE),
        limmaRBE = list(batch_col = "MS_batch", sample_id_col = "FullRunName")
    )
    available_steps <- pb_list_steps()
    bec_steps <- names(param_by_bec)[names(param_by_bec) %in% available_steps]
    expect_gt(length(bec_steps), 0L)

    for (prone_step in step_names) {
        for (bec_step in bec_steps) {
            result <- suppressWarnings(
                pb_eval(
                    object = pbf,
                    from = "feature::raw",
                    steps = c(prone_step, bec_step),
                    params_list = list(
                        list(sample_id_col = "FullRunName"),
                        param_by_bec[[bec_step]]
                    )
                )
            )
            expect_true(is.matrix(result))
            expect_equal(dim(result), dim(pb_assay_matrix(pbf, "feature::raw")))
            expect_identical(colnames(result), colnames(pb_assay_matrix(pbf, "feature::raw")))
            expect_identical(rownames(result), rownames(pb_assay_matrix(pbf, "feature::raw")))
        }
    }
})

test_that("PRONE normalization forwards on_raw explicitly", {
    pbf <- make_prone_norm_test_pbf()
    reg_info <- .pb_register_prone_normalization_steps(
        methods = "mockOnRaw",
        prefix = "PR",
        register_lowercase = FALSE
    )
    step_name <- if (length(reg_info$prefixed)) reg_info$prefixed[[1]] else "PRmockOnRaw"

    captured <- new.env(parent = emptyenv())
    local_mocked_prone_normalize(function(se, methods, on_raw, ...) {
        captured$methods <- methods
        captured$on_raw <- on_raw
        mat <- SummarizedExperiment::assay(se, 1)
        SummarizedExperiment::assay(
            se,
            paste0(methods, "_on_", SummarizedExperiment::assayNames(se)[1]),
            withDimnames = FALSE
        ) <- mat
        se
    })

    out <- pb_eval(
        object = pbf,
        from = "feature::raw",
        steps = step_name,
        params_list = list(list(sample_id_col = "FullRunName", on_raw = FALSE))
    )

    expect_true(is.matrix(out))
    expect_identical(captured$methods, "mockOnRaw")
    expect_identical(captured$on_raw, FALSE)
})

test_that("PRONE normalization provides rowData expected by Normics internals", {
    pbf <- make_prone_norm_test_pbf()
    reg_info <- .pb_register_prone_normalization_steps(
        methods = "mockNormics",
        prefix = "PR",
        register_lowercase = FALSE
    )
    step_name <- if (length(reg_info$prefixed)) reg_info$prefixed[[1]] else "PRmockNormics"

    local_mocked_prone_normalize(function(se, methods, reduce_correlation_by = 1, ...) {
        dt <- data.table::as.data.table(SummarizedExperiment::assay(se, 1))
        rowdata <- data.table::as.data.table(SummarizedExperiment::rowData(se))

        expect_true("Protein.IDs" %in% colnames(rowdata))

        dt_reduced <- dt[seq(1, nrow(dt), by = reduce_correlation_by), ]
        rowdata_reduced <- rowdata[seq(1, nrow(rowdata), by = reduce_correlation_by), ]
        expect_identical(nrow(rowdata_reduced), nrow(dt_reduced))

        mat <- SummarizedExperiment::assay(se, 1)
        SummarizedExperiment::assay(
            se,
            paste0(methods, "_on_", SummarizedExperiment::assayNames(se)[1]),
            withDimnames = FALSE
        ) <- mat
        se
    })

    out <- pb_eval(
        object = pbf,
        from = "feature::raw",
        steps = step_name,
        params_list = list(list(sample_id_col = "FullRunName", reduce_correlation_by = 1))
    )

    expect_true(is.matrix(out))
})

test_that("PRONE normalization forwards batch and refs metadata", {
    pbf <- make_prone_norm_test_pbf()
    reg_info <- .pb_register_prone_normalization_steps(
        methods = "mockMetadata",
        prefix = "PR",
        register_lowercase = FALSE
    )
    step_name <- if (length(reg_info$prefixed)) reg_info$prefixed[[1]] else "PRmockMetadata"
    sample_ids <- colnames(pb_assay_matrix(pbf, "feature::raw"))

    captured <- new.env(parent = emptyenv())
    local_mocked_prone_normalize(function(se, methods, ...) {
        captured$methods <- methods
        captured$batch <- S4Vectors::metadata(se)$batch
        captured$batch_values <- as.character(SummarizedExperiment::colData(se)[[captured$batch]])
        captured$refs <- S4Vectors::metadata(se)$refs
        mat <- SummarizedExperiment::assay(se, 1)
        SummarizedExperiment::assay(
            se,
            paste0(methods, "_on_", SummarizedExperiment::assayNames(se)[1]),
            withDimnames = FALSE
        ) <- mat
        se
    })

    out <- pb_eval(
        object = pbf,
        from = "feature::raw",
        steps = step_name,
        params_list = list(list(
            sample_id_col = "FullRunName",
            batch = "MS_batch",
            refs = sample_ids[1:2]
        ))
    )

    expect_true(is.matrix(out))
    expect_identical(captured$methods, "mockMetadata")
    expect_identical(captured$batch, "MS_batch")
    expect_identical(captured$batch_values, c("B1", "B1", "B1", "B2", "B2", "B2"))
    expect_identical(captured$refs, sample_ids[1:2])
})

test_that("PRONE normalization forwards condition metadata as a column key", {
    pbf <- make_prone_norm_test_pbf()
    reg_info <- .pb_register_prone_normalization_steps(
        methods = "mockCondition",
        prefix = "PR",
        register_lowercase = FALSE
    )
    step_name <- if (length(reg_info$prefixed)) reg_info$prefixed[[1]] else "PRmockCondition"
    condition_vec <- rep(c("ctrl", "case"), each = 3)

    captured <- new.env(parent = emptyenv())
    local_mocked_prone_normalize(function(se, methods, ...) {
        captured$methods <- methods
        captured$condition <- S4Vectors::metadata(se)$condition
        captured$condition_values <- as.character(
            SummarizedExperiment::colData(se)[[captured$condition]]
        )
        mat <- SummarizedExperiment::assay(se, 1)
        SummarizedExperiment::assay(
            se,
            paste0(methods, "_on_", SummarizedExperiment::assayNames(se)[1]),
            withDimnames = FALSE
        ) <- mat
        se
    })

    out <- pb_eval(
        object = pbf,
        from = "feature::raw",
        steps = step_name,
        params_list = list(list(
            sample_id_col = "FullRunName",
            condition = condition_vec
        ))
    )

    expect_true(is.matrix(out))
    expect_identical(captured$methods, "mockCondition")
    expect_true(startsWith(captured$condition, ".pb_prone_condition"))
    expect_identical(captured$condition_values, condition_vec)
})

test_that("PRONE EigenMS step errors clearly when condition is missing", {
    pbf <- make_prone_norm_test_pbf()
    reg_info <- .pb_register_prone_normalization_steps(
        methods = "EigenMS",
        prefix = "PR",
        register_lowercase = FALSE
    )
    step_name <- if (length(reg_info$prefixed)) reg_info$prefixed[[1]] else "PREigenMS"

    local_mocked_prone_normalize(function(se, methods, ain, ...) {
        mat <- SummarizedExperiment::assay(se, ain)
        SummarizedExperiment::assay(
            se,
            paste0(methods, "_on_", ain),
            withDimnames = FALSE
        ) <- mat
        se
    })

    expect_error(
        pb_eval(
            object = pbf,
            from = "feature::raw",
            steps = step_name,
            params_list = list(list(sample_id_col = "FullRunName"))
        ),
        "requires `condition`"
    )
})

test_that("PRONE normalization supports backends that use `ains`", {
    pbf <- make_prone_norm_test_pbf()
    reg_info <- .pb_register_prone_normalization_steps(
        methods = "mockAins",
        prefix = "PR",
        register_lowercase = FALSE
    )
    step_name <- if (length(reg_info$prefixed)) reg_info$prefixed[[1]] else "PRmockAins"

    captured <- new.env(parent = emptyenv())
    local_mocked_prone_normalize(function(se, methods, ains, ...) {
        captured$methods <- methods
        captured$ains <- ains
        mat <- SummarizedExperiment::assay(se, ains)
        SummarizedExperiment::assay(
            se,
            paste0(methods, "_on_", ains),
            withDimnames = FALSE
        ) <- mat
        se
    })

    out <- pb_eval(
        object = pbf,
        from = "feature::raw",
        steps = step_name,
        params_list = list(list(sample_id_col = "FullRunName"))
    )

    expect_true(is.matrix(out))
    expect_identical(captured$methods, "mockAins")
    expect_identical(captured$ains, "raw")
})

test_that("PRONE normalization keeps base args when backend exposes only dots", {
    pbf <- make_prone_norm_test_pbf()
    reg_info <- .pb_register_prone_normalization_steps(
        methods = "mockDotsOnly",
        prefix = "PR",
        register_lowercase = FALSE
    )
    step_name <- if (length(reg_info$prefixed)) reg_info$prefixed[[1]] else "PRmockDotsOnly"

    captured <- new.env(parent = emptyenv())
    local_mocked_prone_normalize(function(...) {
        args <- list(...)
        captured$arg_names <- names(args)

        se <- args$se
        methods <- args$methods
        ain <- args$ain

        mat <- SummarizedExperiment::assay(se, ain)
        SummarizedExperiment::assay(
            se,
            paste0(methods, "_on_", ain),
            withDimnames = FALSE
        ) <- mat
        se
    })

    out <- pb_eval(
        object = pbf,
        from = "feature::raw",
        steps = step_name,
        params_list = list(list(sample_id_col = "FullRunName"))
    )

    expect_true(is.matrix(out))
    expect_true("se" %in% captured$arg_names)
    expect_true("methods" %in% captured$arg_names)
    expect_true("ain" %in% captured$arg_names)
    expect_true("ains" %in% captured$arg_names)
})

make_prone_norm_assay_se <- function(assay_names) {
    assays <- lapply(assay_names, function(name) {
        matrix(
            c(1, 2, 3, 4),
            nrow = 2,
            dimnames = list(c("f1", "f2"), c("s1", "s2"))
        )
    })
    names(assays) <- assay_names
    SummarizedExperiment::SummarizedExperiment(assays = assays)
}

test_that(".pb_prone_guess_norm_assay resolves known fallback cases", {
    se_by_added_assay <- make_prone_norm_assay_se(c("raw", "something_else"))
    expect_identical(
        .pb_prone_guess_norm_assay(se_by_added_assay, ain = "raw", norm_method = "mock"),
        "something_else"
    )

    se_by_on_raw <- make_prone_norm_assay_se(c("raw", "norm_on_raw", "extra"))
    expect_identical(
        .pb_prone_guess_norm_assay(se_by_on_raw, ain = "raw", norm_method = "mock"),
        "norm_on_raw"
    )

    se_single <- make_prone_norm_assay_se("only_assay")
    expect_identical(
        .pb_prone_guess_norm_assay(se_single, ain = "raw", norm_method = "mock"),
        "only_assay"
    )
})

test_that(".pb_prone_guess_norm_assay errors on ambiguous assay outputs", {
    se_ambiguous <- make_prone_norm_assay_se(c("raw", "assayA", "assayB"))

    expect_error(
        .pb_prone_guess_norm_assay(se_ambiguous, ain = "raw", norm_method = "mock"),
        "Unable to determine normalized assay returned by PRONE"
    )
})

# -----------------------------------------------------------------------
# Verify PRONE normalization does not introduce all-NA rows/columns
# -----------------------------------------------------------------------

test_that("PRONE normalization step does not introduce all-NA rows when input has none", {
    skip_if_not_installed("PRONE")

    # Matrix with partial NAs but no all-NA rows or columns
    dm <- matrix(
        c(
            10, NA, 30,
            40, 50, NA,
            70, 80, 90
        ),
        nrow = 3, byrow = TRUE,
        dimnames = list(paste0("feat", 1:3), paste0("s", 1:3))
    )
    sa <- data.frame(
        FullRunName = paste0("s", 1:3),
        stringsAsFactors = FALSE
    )

    # Confirm no all-NA rows/cols initially
    expect_false(any(apply(dm, 1, function(r) all(is.na(r)))))
    expect_false(any(apply(dm, 2, function(c) all(is.na(c)))))

    # Mock PRONE normalize to return a "normalized" assay
    local_mocked_prone_normalize(function(se, methods, method, norm_method,
                                          on_raw, ain, ains,
                                          aout, combination_pattern, ...) {
        nm <- methods %||% method %||% norm_method %||% "MockNorm"
        mat <- SummarizedExperiment::assay(se, ain)
        mat <- mat * 2 # simple transformation, preserves NAs
        SummarizedExperiment::assay(se, nm) <- mat
        se
    })

    res <- proBatch:::.prone_normalize_matrix_step(
        data_matrix = dm,
        sample_annotation = sa,
        sample_id_col = "FullRunName",
        norm_method = "MockNorm",
        assay_in = "raw"
    )

    # Output must not have new all-NA rows or columns
    expect_false(any(apply(res, 1, function(r) all(is.na(r)))),
        info = "PRONE normalization should not introduce all-NA rows"
    )
    expect_false(any(apply(res, 2, function(c) all(is.na(c)))),
        info = "PRONE normalization should not introduce all-NA columns"
    )
    expect_identical(dim(res), dim(dm))
    expect_identical(dimnames(res), dimnames(dm))
})
