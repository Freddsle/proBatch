test_that("validate_batch_design flags duplicates, missingness, singleton batches, and nesting", {
    df <- data.frame(
        sample_id = c("s1", "s2", "s3", "s4", "s5", "s5"),
        batch = c("b1", "b1", "b2", "b2", "b2", "b3"),
        condition = c("A", "A", "A", "A", "A", "B"),
        cov1 = c(1, 2, NA, 4, 5, 6),
        stringsAsFactors = FALSE
    )

    res <- validate_batch_design(
        df,
        batch_col = "batch",
        condition_col = "condition",
        covariates = "cov1",
        sample_id_col = "sample_id",
        strict = FALSE
    )

    expect_s3_class(res, "pb_design_check")
    expect_true(any(grepl("Duplicated", res$errors)))
    expect_true(any(grepl("missing values", res$errors, ignore.case = TRUE)))
    expect_true(any(grepl("Single-sample", res$warnings)))
    expect_true(any(grepl("nested", res$warnings, ignore.case = TRUE)))
})

test_that("summarize_design detects rank deficiency with confounded factors", {
    df <- data.frame(
        batch = factor(rep(c("b1", "b2"), each = 3)),
        condition = factor(rep(c("b1", "b2"), each = 3))
    )

    summary <- summarize_design(
        df,
        batch_col = "batch",
        condition_col = "condition"
    )

    expect_s3_class(summary, "pb_design_summary")
    expect_true(summary$deficiency > 0)
    expect_true(length(summary$aliased_terms) > 0)
})

test_that("summarize_design drops non-informative one-level terms", {
    df <- data.frame(
        batch = factor(rep("b1", 6)),
        condition = factor(rep(c("A", "B"), each = 3)),
        covariate = factor(rep("only_level", 6)),
        stringsAsFactors = FALSE
    )

    summary <- summarize_design(
        df,
        batch_col = "batch",
        condition_col = "condition",
        covariates = "covariate"
    )

    expect_s3_class(summary, "pb_design_summary")
    expect_equal(summary$design_matrix_dim[["n_terms"]], 2)
    expect_true(any(grepl("Dropped non-informative design terms", summary$notes)))
})

test_that("detect_nested_batches identifies nesting in synthetic metadata", {
    df <- data.frame(
        site = c("S1", "S1", "S2", "S2", "S3", "S3"),
        run = c("R1", "R1", "R1", "R1", "R2", "R2"),
        plate = c("P1", "P2", "P1", "P2", "P1", "P2"),
        stringsAsFactors = FALSE
    )

    res <- detect_nested_batches(df, batch_cols = c("site", "run", "plate"))

    expect_s3_class(res, "pb_nested_batches")
    expect_true(res$adjacency["run", "site"])
    expect_false(res$adjacency["site", "run"])
    expect_true(is.null(res$ordering) || match("run", res$ordering) < match("site", res$ordering))
})

test_that("detect_outlier_samples flags an injected outlier", {
    mat <- matrix(0, nrow = 2, ncol = 6)
    colnames(mat) <- paste0("s", 1:6)
    mat[, 6] <- 100

    out <- detect_outlier_samples(mat, n_pcs = 1, robust = FALSE, cutoff = 0.95)

    expect_s3_class(out, "pb_outliers")
    expect_true(out$is_outlier[out$sample_id == "s6"])
    expect_identical(attr(out, "outlier_mode"), "classical")
})

test_that("detect_outlier_samples falls back when robust fit flags whole batches", {
    skip_if_not_installed("MASS")

    set.seed(7)
    n_features <- 400
    batch_sizes <- c(40, 40, 10, 10, 10)
    batch_labels <- paste0("B", seq_along(batch_sizes))
    shift_means <- c(0, 0.2, 1.6, 2.0, 2.4)

    mat_list <- lapply(seq_along(batch_labels), function(i) {
        batch_noise <- matrix(
            rnorm(n_features * batch_sizes[i], sd = 0.9),
            nrow = n_features
        )
        batch_shift <- rnorm(n_features, mean = shift_means[i], sd = 0.2)
        batch_noise + batch_shift
    })
    mat <- do.call(cbind, mat_list)
    colnames(mat) <- paste0("s", seq_len(ncol(mat)))

    meta <- data.frame(
        batch = rep(batch_labels, times = batch_sizes),
        stringsAsFactors = FALSE
    )
    rownames(meta) <- colnames(mat)

    expect_warning(
        out <- detect_outlier_samples(
            mat,
            sample_annotation = meta,
            batch_col = "batch",
            n_pcs = 5,
            robust = TRUE,
            cutoff = 0.99
        ),
        "falling back to classical covariance"
    )

    batch_summary <- attr(out, "batch_summary")
    expect_true(max(batch_summary$pct_outliers, na.rm = TRUE) < 100)
    expect_identical(attr(out, "outlier_mode"), "classical_fallback")
})

test_that("subbatch_detection splits a batch with two clusters", {
    set.seed(42)
    n_features <- 40
    cluster1 <- matrix(rnorm(n_features * 10, mean = 0), nrow = n_features)
    cluster2 <- matrix(rnorm(n_features * 10, mean = 5), nrow = n_features)
    mat <- cbind(cluster1, cluster2)
    colnames(mat) <- paste0("s", 1:20)

    meta <- data.frame(
        sample_id = colnames(mat),
        batch = rep("B1", ncol(mat)),
        stringsAsFactors = FALSE
    )

    res <- subbatch_detection(
        mat,
        sample_annotation = meta,
        batch_col = "batch",
        n_pcs = 5,
        method = "kmeans",
        k_max = 3
    )

    expect_s3_class(res, "pb_subbatches")
    b1_row <- res$summary[res$summary$batch == "B1", , drop = FALSE]
    expect_true(nrow(b1_row) == 1)
    expect_true(b1_row$k >= 2)
    expect_true(length(unique(res$assignments$subbatch)) >= 2)
})

test_that("design diagnostics accept ProBatchFeatures inputs", {
    mat_outlier <- matrix(0, nrow = 2, ncol = 6)
    colnames(mat_outlier) <- paste0("s", 1:6)
    mat_outlier[, 6] <- 100
    meta_outlier <- data.frame(
        FullRunName = colnames(mat_outlier),
        batch = "B1",
        stringsAsFactors = FALSE
    )

    pbf_outlier <- suppressMessages(ProBatchFeatures(
        data_matrix = mat_outlier,
        sample_annotation = meta_outlier,
        sample_id_col = "FullRunName",
        name = "feature::raw"
    ))

    outliers <- detect_outlier_samples(
        pbf_outlier,
        batch_col = "batch",
        n_pcs = 1,
        robust = FALSE,
        cutoff = 0.95
    )
    expect_s3_class(outliers, "pb_outliers")
    expect_true(outliers$is_outlier[outliers$sample_id == "s6"])

    set.seed(42)
    n_features <- 40
    cluster1 <- matrix(rnorm(n_features * 10, mean = 0), nrow = n_features)
    cluster2 <- matrix(rnorm(n_features * 10, mean = 5), nrow = n_features)
    mat_subbatch <- cbind(cluster1, cluster2)
    colnames(mat_subbatch) <- paste0("s", 1:20)
    meta_subbatch <- data.frame(
        FullRunName = colnames(mat_subbatch),
        batch = rep("B1", ncol(mat_subbatch)),
        stringsAsFactors = FALSE
    )

    pbf_subbatch <- suppressMessages(ProBatchFeatures(
        data_matrix = mat_subbatch,
        sample_annotation = meta_subbatch,
        sample_id_col = "FullRunName",
        name = "feature::raw"
    ))

    subbatches <- subbatch_detection(
        pbf_subbatch,
        batch_col = "batch",
        n_pcs = 5,
        method = "kmeans",
        k_max = 3
    )
    expect_s3_class(subbatches, "pb_subbatches")
    b1_row <- subbatches$summary[subbatches$summary$batch == "B1", , drop = FALSE]
    expect_true(nrow(b1_row) == 1)
    expect_true(b1_row$k >= 2)
})

test_that("design diagnostics honor explicit pbf_name assay selection", {
    set.seed(1)
    n_features <- 120
    n_per_batch <- 20

    batch1 <- matrix(rexp(n_features * n_per_batch, rate = 1 / 100), nrow = n_features)
    batch2 <- matrix(rexp(n_features * n_per_batch, rate = 1 / 100), nrow = n_features)
    shifted_idx <- sample(seq_len(n_features), size = round(0.35 * n_features))
    batch2[shifted_idx, ] <- batch2[shifted_idx, ] * runif(length(shifted_idx), 1.8, 5.0)

    mat_raw <- cbind(batch1, batch2)
    colnames(mat_raw) <- paste0("s", seq_len(ncol(mat_raw)))
    mat_log <- log2(mat_raw + 1)

    meta <- data.frame(
        FullRunName = colnames(mat_raw),
        batch = rep(c("B1", "B2"), each = n_per_batch),
        stringsAsFactors = FALSE
    )

    pbf <- suppressMessages(ProBatchFeatures(
        data_matrix = mat_raw,
        sample_annotation = meta,
        sample_id_col = "FullRunName",
        name = "feature::raw"
    ))
    pbf <- pb_transform(pbf, from = "feature::raw", steps = "log2")

    raw_assay <- "feature::raw"
    log_assay <- as.character(utils::tail(get_operation_log(pbf)$to, 1))
    expect_false(identical(raw_assay, log_assay))

    out_raw_mat <- detect_outlier_samples(
        mat_raw,
        sample_annotation = meta,
        batch_col = "batch",
        n_pcs = 5,
        robust = FALSE,
        cutoff = 0.99
    )
    out_log_mat <- detect_outlier_samples(
        mat_log,
        sample_annotation = meta,
        batch_col = "batch",
        n_pcs = 5,
        robust = FALSE,
        cutoff = 0.99
    )
    out_raw_pbf <- detect_outlier_samples(
        pbf,
        pbf_name = raw_assay,
        sample_annotation = meta,
        batch_col = "batch",
        n_pcs = 5,
        robust = FALSE,
        cutoff = 0.99
    )
    out_log_pbf <- detect_outlier_samples(
        pbf,
        pbf_name = log_assay,
        sample_annotation = meta,
        batch_col = "batch",
        n_pcs = 5,
        robust = FALSE,
        cutoff = 0.99
    )

    expect_equal(out_raw_pbf$is_outlier, out_raw_mat$is_outlier)
    expect_equal(out_log_pbf$is_outlier, out_log_mat$is_outlier)
    expect_false(identical(out_raw_pbf$is_outlier, out_log_pbf$is_outlier))

    set.seed(123)
    sub_raw_mat <- subbatch_detection(
        mat_raw,
        sample_annotation = meta,
        batch_col = "batch",
        n_pcs = 5,
        method = "kmeans",
        k_max = 4
    )
    set.seed(123)
    sub_log_mat <- subbatch_detection(
        mat_log,
        sample_annotation = meta,
        batch_col = "batch",
        n_pcs = 5,
        method = "kmeans",
        k_max = 4
    )
    set.seed(123)
    sub_raw_pbf <- subbatch_detection(
        pbf,
        pbf_name = raw_assay,
        sample_annotation = meta,
        batch_col = "batch",
        n_pcs = 5,
        method = "kmeans",
        k_max = 4
    )
    set.seed(123)
    sub_log_pbf <- subbatch_detection(
        pbf,
        pbf_name = log_assay,
        sample_annotation = meta,
        batch_col = "batch",
        n_pcs = 5,
        method = "kmeans",
        k_max = 4
    )

    expect_equal(sub_raw_pbf$summary$k, sub_raw_mat$summary$k)
    expect_equal(sub_log_pbf$summary$k, sub_log_mat$summary$k)
    expect_false(identical(sub_raw_pbf$summary$k, sub_log_pbf$summary$k))
})
