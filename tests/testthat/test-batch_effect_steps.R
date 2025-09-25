test_that("pb_transform applies combat step via registry", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

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
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

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
