test_that("pb_transform applies combat step via registry", {
    pb_test_load_example_data()

    batch1_ids <- head(example_sample_annotation$FullRunName[example_sample_annotation$MS_batch == "Batch_1"], 4)
    batch2_ids <- head(example_sample_annotation$FullRunName[example_sample_annotation$MS_batch == "Batch_2"], 4)
    sample_ids <- unique(c(batch1_ids, batch2_ids))
    sample_ids <- sample_ids[!is.na(sample_ids)]

    expect_gt(length(sample_ids), 4L)

    sample_ann <- example_sample_annotation[
        match(sample_ids, example_sample_annotation$FullRunName), ,
        drop = FALSE
    ]

    feature_ids <- head(rownames(example_proteome_matrix), 6)
    sub_matrix <- example_proteome_matrix[feature_ids, sample_ids, drop = FALSE]

    pbf <- ProBatchFeatures(
        data_matrix = sub_matrix,
        sample_annotation = sample_ann,
        sample_id_col = "FullRunName",
        name = "raw"
    )

    params <- list(
        sample_annotation = as.data.frame(colData(pbf)),
        batch_col = "MS_batch",
        sample_id_col = "FullRunName",
        par.prior = TRUE
    )

    pbf2 <- pb_transform(
        pbf,
        from = "feature::raw",
        steps = "combat",
        params_list = list(params)
    )

    expect_true("feature::combat_on_raw" %in% names(pbf2))

    log <- get_operation_log(pbf2)
    combat_idx <- which(as.character(log$step) == "combat")
    expect_gt(length(combat_idx), 0L)
    combat_entry <- log[combat_idx[length(combat_idx)], , drop = FALSE]
    expect_identical(as.character(combat_entry$from), "feature::raw")
    expect_identical(as.character(combat_entry$to), "feature::combat_on_raw")
    expect_equal(combat_entry$params[[1]]$batch_col, "MS_batch")

    m_raw <- pb_assay_matrix(pbf, "feature::raw")
    m_combat <- pb_assay_matrix(pbf2, "feature::combat_on_raw")
    expect_equal(dim(m_combat), dim(m_raw))
    expect_gt(sum(abs(m_combat - m_raw)), 0)
})

test_that("pb_transform applies limmaRBE step via registry", {
    pb_test_load_example_data()

    batch_ids <- example_sample_annotation$FullRunName[example_sample_annotation$MS_batch %in% c("Batch_1", "Batch_2")]
    sample_ids <- unique(batch_ids)
    sample_ids <- sample_ids[!is.na(sample_ids)]

    sample_ann <- example_sample_annotation[
        match(sample_ids, example_sample_annotation$FullRunName), ,
        drop = FALSE
    ]

    sub_matrix <- example_proteome_matrix[, sample_ids, drop = FALSE]
    # drop rows with NAs
    sub_matrix <- sub_matrix[complete.cases(sub_matrix), , drop = FALSE]

    pbf <- ProBatchFeatures(
        data_matrix = sub_matrix,
        sample_annotation = sample_ann,
        sample_id_col = "FullRunName",
        name = "raw"
    )

    params <- list(
        sample_annotation = as.data.frame(colData(pbf)),
        batch_col = "MS_batch",
        covariates_cols = c("Diet", "Sex")
    )

    pbf2 <- pb_transform(
        pbf,
        from = "feature::raw",
        steps = "limmaRBE",
        params_list = list(params)
    )

    expect_true("feature::limmaRBE_on_raw" %in% names(pbf2))

    log <- get_operation_log(pbf2)
    limma_idx <- which(as.character(log$step) == "limmaRBE")
    expect_gt(length(limma_idx), 0L)
    limma_entry <- log[limma_idx[length(limma_idx)], , drop = FALSE]
    expect_identical(as.character(limma_entry$from), "feature::raw")
    expect_identical(as.character(limma_entry$to), "feature::limmaRBE_on_raw")
    expect_equal(limma_entry$params[[1]]$covariates_cols, c("Diet", "Sex"))

    m_raw <- pb_assay_matrix(pbf, "feature::raw")
    m_limma <- pb_assay_matrix(pbf2, "feature::limmaRBE_on_raw")
    expect_equal(dim(m_limma), dim(m_raw))
    expect_gt(sum(abs(m_limma - m_raw), na.rm = TRUE), 0)
})

test_that("medianNorm step handles fill_the_missing via pb_transform", {
    pb_test_load_example_data()

    sub_matrix <- example_proteome_matrix[1:5, 1:5, drop = FALSE]
    sub_matrix[1, 2] <- NA_real_

    sample_ann <- example_sample_annotation[
        match(colnames(sub_matrix), example_sample_annotation$FullRunName), ,
        drop = FALSE
    ]

    pbf <- ProBatchFeatures(
        data_matrix = sub_matrix,
        sample_annotation = sample_ann,
        sample_id_col = "FullRunName",
        name = "raw"
    )

    params <- list(
        sample_annotation = as.data.frame(colData(pbf)),
        sample_id_col = "FullRunName",
        fill_the_missing = 0
    )

    suppressWarnings(
        pbf2 <- pb_transform(
            pbf,
            from = "feature::raw",
            steps = "medianNorm",
            params_list = list(params),
            store_fast_steps = TRUE
        )
    )

    assay_name <- "feature::medianNorm_on_raw"
    expect_true(assay_name %in% names(pbf2))

    m_norm <- pb_assay_matrix(pbf2, assay_name)
    expect_false(anyNA(m_norm))

    log <- get_operation_log(pbf2)
    median_idx <- which(as.character(log$step) == "medianNorm")
    expect_gt(length(median_idx), 0L)
    median_entry <- log[median_idx[length(median_idx)], , drop = FALSE]
    expect_identical(median_entry$params[[1]]$fill_the_missing, 0)
})

test_that("loessLimmaRBE step is registered", {
    testthat::skip_if_not_installed("proBatch")
    expect_true(pb_has_step("loessLimmaRBE"))
})

test_that("pb_transform applies loessLimmaRBE step via registry", {
    testthat::skip_if_not_installed("proBatch")
    pb_test_load_example_data()

    batch1_ids <- example_sample_annotation$FullRunName[example_sample_annotation$MS_batch == "Batch_1"]
    batch2_ids <- example_sample_annotation$FullRunName[example_sample_annotation$MS_batch == "Batch_2"]
    batch1_ids <- unique(batch1_ids[!is.na(batch1_ids)])
    batch2_ids <- unique(batch2_ids[!is.na(batch2_ids)])
    keep_n <- min(length(batch1_ids), length(batch2_ids), 12L)
    testthat::skip_if(keep_n < 6L, "Need at least 6 samples per batch for loessLimmaRBE test")

    sample_ids <- c(head(batch1_ids, keep_n), head(batch2_ids, keep_n))
    sample_ids <- unique(sample_ids)

    sample_ann <- example_sample_annotation[
        match(sample_ids, example_sample_annotation$FullRunName), ,
        drop = FALSE
    ]

    candidate_matrix <- example_proteome_matrix[, sample_ids, drop = FALSE]
    complete_feature_ids <- rownames(candidate_matrix)[complete.cases(candidate_matrix)]
    testthat::skip_if(
        length(complete_feature_ids) < 8L,
        "Need at least 8 complete features for loessLimmaRBE test"
    )
    sub_matrix <- candidate_matrix[head(complete_feature_ids, 8), , drop = FALSE]
    testthat::skip_if(nrow(sub_matrix) < 2L, "Need at least two complete features for loessLimmaRBE test")

    candidate_covariates <- c("Diet", "Sex")
    valid_covariates <- candidate_covariates[vapply(
        candidate_covariates,
        function(column_name) {
            column_name %in% names(sample_ann) &&
                length(unique(sample_ann[[column_name]][!is.na(sample_ann[[column_name]])])) >= 2L
        },
        logical(1)
    )]
    if (!length(valid_covariates)) {
        valid_covariates <- NULL
    }

    pbf <- ProBatchFeatures(
        data_matrix = sub_matrix,
        sample_annotation = sample_ann,
        sample_id_col = "FullRunName",
        name = "raw"
    )

    params <- list(
        batch_col = "MS_batch",
        sample_id_col = "FullRunName",
        order_col = "order",
        covariates_cols = valid_covariates,
        fill_the_missing = FALSE,
        min_measurements = 6,
        no_fit_imputed = FALSE,
        span = 1
    )

    pbf2 <- suppressWarnings(
        pb_transform(
            pbf,
            from = "feature::raw",
            steps = "loessLimmaRBE",
            params_list = list(params)
        )
    )

    expect_true("feature::loessLimmaRBE_on_raw" %in% names(pbf2))

    log <- get_operation_log(pbf2)
    loess_idx <- which(as.character(log$step) == "loessLimmaRBE")
    expect_gt(length(loess_idx), 0L)
    loess_entry <- log[loess_idx[length(loess_idx)], , drop = FALSE]
    expect_identical(as.character(loess_entry$from), "feature::raw")
    expect_identical(as.character(loess_entry$to), "feature::loessLimmaRBE_on_raw")
    expect_identical(loess_entry$params[[1]]$order_col, "order")

    m_step <- pb_assay_matrix(pbf2, "feature::loessLimmaRBE_on_raw")
    m_expected <- suppressWarnings(
        correct_batch_effects(
            x = pb_assay_matrix(pbf, "feature::raw"),
            sample_annotation = as.data.frame(colData(pbf)),
            format = "wide",
            continuous_func = "loess_regression",
            discrete_func = "removeBatchEffect",
            batch_col = "MS_batch",
            sample_id_col = "FullRunName",
            order_col = "order",
            covariates_cols = valid_covariates,
            fill_the_missing = FALSE,
            min_measurements = 6,
            no_fit_imputed = FALSE,
            span = 1
        )
    )

    expect_equal(dim(m_step), dim(m_expected))
    expect_equal(m_step, m_expected, tolerance = 1e-8)
})
