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
