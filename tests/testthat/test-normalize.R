test_that("quantile_normalize", {
    data(example_proteome_matrix, package = "proBatch")

    matrix_test <- example_proteome_matrix[1:10, ]
    quant_normalized_matrix <- quantile_normalize_dm(matrix_test)

    quant_matrix <- normalize.quantiles(matrix_test)

    expect_equal(quant_matrix[1:10], quant_normalized_matrix[1:10], ignore_attr = TRUE)
})

test_that("normalize_sample_medians", {
    data(example_proteome, package = "proBatch")

    samples <- c(
        "I170925_BXD66_HF_ET1560_SW_Run001",
        "I170925_BXD101_CD_ET1424_SW_Run002",
        "I170925_BXD63_HF_ET1715_SW_Run003",
        "I170925_BXD61_HF_ET1721_SW_Run004"
    )
    sampled_rows <- example_proteome$FullRunName %in% samples
    sample_df <- example_proteome[sampled_rows, ]
    median_centered_data <- normalize_sample_medians_df(sample_df,
        keep_all = "all"
    )

    expect_equal(median_centered_data$median_global[1], 30277.6, tolerance = 1e-2, ignore_attr = TRUE)
    expect_equal(median_centered_data$median_run[2], 29252.14, tolerance = 1e-2, ignore_attr = TRUE)

    expect_equal(median_centered_data$Intensity[1], 2728887.42, tolerance = 1e-2, ignore_attr = TRUE)
    expect_equal(median_centered_data$Intensity[549], 35191.95, tolerance = 1e-2, ignore_attr = TRUE)
    expect_equal(median_centered_data$preNorm_Intensity[549], 25479.8, ignore_attr = TRUE)
})

test_that("normalize_sample_medians inside batch using metadata", {
    mat <- matrix(
        c(
            10, 20, 30, 40,
            15, 25, 35, 45,
            NA, NA, 5, 10,
            NA, NA, NA, NA
        ),
        nrow = 4,
        byrow = TRUE
    )
    rownames(mat) <- paste0("F", seq_len(nrow(mat)))
    colnames(mat) <- paste0("S", seq_len(ncol(mat)))
    sample_annot <- data.frame(
        FullRunName = colnames(mat),
        Batch = c("A", "A", "B", "B"),
        stringsAsFactors = FALSE
    )

    expected <- matrix(
        c(
            15, 15, 32.5, 32.5,
            20, 20, 37.5, 37.5,
            NA, NA, 7.5, 2.5,
            NA, NA, NA, NA
        ),
        nrow = 4,
        byrow = TRUE,
        dimnames = dimnames(mat)
    )

    inside_batch_norm <- normalize_sample_medians_dm(
        mat,
        sample_annotation = sample_annot,
        group_col = "Batch",
        inside_batch = TRUE
    )
    expect_equal(inside_batch_norm, expected, tolerance = 1e-8, ignore_attr = TRUE)

    df_long <- matrix_to_long(mat)
    inside_batch_df <- normalize_sample_medians_df(
        df_long,
        sample_annotation = sample_annot,
        group_col = "Batch",
        inside_batch = TRUE
    )
    reconstructed <- long_to_matrix(
        inside_batch_df,
        feature_id_col = "peptide_group_label",
        sample_id_col = "FullRunName",
        measure_col = "Intensity"
    )
    expect_equal(reconstructed, expected, tolerance = 1e-8, ignore_attr = TRUE)

    stats_df <- unique(inside_batch_df[, c("FullRunName", "median_run", "median_global")])
    stats_df <- stats_df[order(stats_df$FullRunName), ]
    expect_equal(stats_df$median_run, c(12.5, 22.5, 30, 40))
    expect_equal(stats_df$median_global, c(17.5, 17.5, 32.5, 32.5))

    pbf <- ProBatchFeatures(
        data_matrix = mat,
        sample_annotation = sample_annot
    )
    pbf_norm <- pb_eval(
        pbf,
        from = "feature::raw",
        steps = "medianNorm",
        params_list = list(list(group_col = "Batch", inside_batch = TRUE))
    )
    expect_equal(pbf_norm, expected, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("normalize_data_dm forwards offset to log transform", {
    data(example_proteome_matrix, package = "proBatch")
    mat <- example_proteome_matrix[1:5, 1:4]
    res_quant <- normalize_data_dm(mat, normalize_func = "quantile", log_base = 2, offset = 0.5)
    expected_quant <- quantile_normalize_dm(log_transform_dm(mat, log_base = 2, offset = 0.5))
    expect_equal(res_quant, expected_quant, ignore_attr = TRUE)

    res_median <- normalize_data_dm(mat, normalize_func = "medianCentering", log_base = 2, offset = 0.5)
    expected_median <- normalize_sample_medians_dm(log_transform_dm(mat, log_base = 2, offset = 0.5))
    expect_equal(res_median, expected_median, ignore_attr = TRUE)
})
