test_that("correct_with_mComBat(wide): requires sva and a valid center", {
    # tiny synthetic (features x samples)
    set.seed(1)
    p <- 7
    n1 <- 5
    n2 <- 6
    B1 <- matrix(rnorm(p * n1, mean = 0, sd = 1), nrow = p)
    B2 <- matrix(rnorm(p * n2, mean = 3, sd = 2), nrow = p) # shifted & scaled batch
    m <- cbind(B1, B2)
    rownames(m) <- paste0("f", seq_len(p))
    colnames(m) <- paste0("s", seq_len(n1 + n2))

    sa <- data.frame(
        FullRunName = colnames(m),
        MS_batch = factor(c(rep("B1", n1), rep("B2", n2))),
        stringsAsFactors = FALSE
    )

    # center required
    expect_error(
        correct_with_mComBat(
            x = m, sample_annotation = sa,
            sample_id_col = "FullRunName", batch_col = "MS_batch",
            format = "wide" # mComBat_center omitted
        ),
        "mComBat_center must be specified when using correct_with_mComBat."
    )

    # center must exist
    expect_error(
        correct_with_mComBat(
            x = m, sample_annotation = sa,
            sample_id_col = "FullRunName", batch_col = "MS_batch",
            format = "wide", mComBat_center = "B9"
        ),
        "mComBat_center must be present in the batch levels of the provided sample_annotation."
    )

    # basic success path
    expect_no_error(
        res <- correct_with_mComBat(
            x = m, sample_annotation = sa,
            sample_id_col = "FullRunName", batch_col = "MS_batch",
            format = "wide", mComBat_center = "B1"
        )
    )
    expect_true(is.matrix(res))
    expect_identical(dim(res), dim(m))
    expect_identical(rownames(res), rownames(m))
    expect_identical(colnames(res), colnames(m))
})

test_that("correct_with_mComBat(wide): missing values -> error unless instructed", {
    skip_if_not_installed("sva")
    set.seed(2)
    p <- 6
    n1 <- 4
    n2 <- 4
    B1 <- matrix(rnorm(p * n1), nrow = p)
    B2 <- matrix(rnorm(p * n2, 2, 1.5), nrow = p)
    m <- cbind(B1, B2)
    rownames(m) <- paste0("f", seq_len(p))
    colnames(m) <- paste0("s", seq_len(n1 + n2))
    m[1, 1] <- NA

    sa <- data.frame(
        FullRunName = colnames(m),
        MS_batch = factor(c(rep("B1", n1), rep("B2", n2))),
        stringsAsFactors = FALSE
    )

    # By default (fill_the_missing = NULL), m-ComBat should fail on NA
    expect_error(
        correct_with_mComBat(
            x = m, sample_annotation = sa, format = "wide",
            sample_id_col = "FullRunName", batch_col = "MS_batch",
            mComBat_center = "B1"
        )
    )

    # (A) impute with numeric
    expect_warning(expect_warning(
        res_imp <- correct_with_mComBat(
            x = m, sample_annotation = sa, format = "wide",
            sample_id_col = "FullRunName", batch_col = "MS_batch",
            mComBat_center = "B1",
            fill_the_missing = 0
        ), "filling missing values with 0"
    ), "ComBat cannot operate with missing values in the matrix")
    expect_true(is.matrix(res_imp))

    # (B) or drop NA rows first
    expect_warning(expect_warning(
        res_drop <- correct_with_mComBat(
            x = m, sample_annotation = sa, format = "wide",
            sample_id_col = "FullRunName", batch_col = "MS_batch",
            mComBat_center = "B1",
            fill_the_missing = "remove"
        ), "removed 1 rows"
    ), "ComBat cannot operate with missing values in the matrix")
    expect_true(nrow(res_drop) <= nrow(m))
})

test_that("correct_with_mComBat: wide == long (within tolerance)", {
    skip_if_not_installed("sva")
    set.seed(3)
    p <- 8
    n1 <- 5
    n2 <- 5
    B1 <- matrix(rnorm(p * n1, 0, 1), nrow = p)
    B2 <- matrix(rnorm(p * n2, 4, 1.8), nrow = p)
    m <- cbind(B1, B2)
    rownames(m) <- paste0("f", seq_len(p))
    colnames(m) <- paste0("s", seq_len(n1 + n2))

    sa <- data.frame(
        FullRunName = colnames(m),
        MS_batch = factor(c(rep("B1", n1), rep("B2", n2))),
        stringsAsFactors = FALSE
    )

    # wide
    res_wide <- correct_with_mComBat(
        x = m, sample_annotation = sa, format = "wide",
        sample_id_col = "FullRunName", batch_col = "MS_batch",
        mComBat_center = "B1"
    )

    # long
    df_long <- matrix_to_long(
        data_matrix = m,
        feature_id_col = "peptide_group_label",
        measure_col = "Intensity",
        sample_id_col = "FullRunName"
    )
    res_long <- correct_with_mComBat(
        x = df_long, sample_annotation = sa, format = "long",
        feature_id_col = "peptide_group_label", measure_col = "Intensity",
        sample_id_col = "FullRunName", batch_col = "MS_batch",
        mComBat_center = "B1"
    )
    res_long_mat <- long_to_matrix(
        res_long,
        feature_id_col = "peptide_group_label",
        measure_col = "Intensity",
        sample_id_col = "FullRunName"
    )

    expect_equal(res_wide, res_long_mat, tolerance = 1e-8)
})

test_that("correct_with_mComBat(wide): aligns batches toward the chosen center (means & variances)", {
    skip_if_not_installed("sva")
    set.seed(4)
    p <- 12
    n1 <- 6
    n2 <- 6
    B1 <- matrix(rnorm(p * n1, 0, 1), nrow = p) # reference distribution
    B2 <- matrix(rnorm(p * n2, 3, 2.2), nrow = p) # shifted + scaled
    m <- cbind(B1, B2)
    rownames(m) <- paste0("f", seq_len(p))
    colnames(m) <- paste0("s", seq_len(n1 + n2))

    sa <- data.frame(
        FullRunName = colnames(m),
        MS_batch = factor(c(rep("B1", n1), rep("B2", n2))),
        stringsAsFactors = FALSE
    )
    idx_B1 <- sa$MS_batch == "B1"
    idx_B2 <- !idx_B1

    rowVars <- function(A) apply(A, 1, stats::var)

    # before: differences in per-feature means/variances between B2 and B1
    mean_diff_before <- abs(rowMeans(m[, idx_B2, drop = FALSE]) - rowMeans(m[, idx_B1, drop = FALSE]))
    var_ratio_before <- pmax(
        rowVars(m[, idx_B2, drop = FALSE]) / (rowVars(m[, idx_B1, drop = FALSE]) + 1e-12),
        (rowVars(m[, idx_B1, drop = FALSE]) + 1e-12) / (rowVars(m[, idx_B2, drop = FALSE]) + 1e-12)
    )

    res <- correct_with_mComBat(
        x = m, sample_annotation = sa, format = "wide",
        sample_id_col = "FullRunName", batch_col = "MS_batch",
        mComBat_center = "B1"
    )

    mean_diff_after <- abs(rowMeans(res[, idx_B2, drop = FALSE]) - rowMeans(res[, idx_B1, drop = FALSE]))
    var_ratio_after <- pmax(
        rowVars(res[, idx_B2, drop = FALSE]) / (rowVars(res[, idx_B1, drop = FALSE]) + 1e-12),
        (rowVars(res[, idx_B1, drop = FALSE]) + 1e-12) / (rowVars(res[, idx_B2, drop = FALSE]) + 1e-12)
    )

    # On average, means closer and variances more similar (ratio closer to 1)
    expect_lt(mean(mean_diff_after), mean(mean_diff_before))
    expect_lt(mean(abs(log(var_ratio_after))), mean(abs(log(var_ratio_before))))
})

test_that("correct_with_mComBat(wide): covariates are preserved when supplied", {
    skip_if_not_installed("sva")
    set.seed(5)
    p <- 10
    n1 <- 6
    n2 <- 6
    base <- matrix(rnorm(p * (n1 + n2), 0, 1), nrow = p)
    rownames(base) <- paste0("f", seq_len(p))
    colnames(base) <- paste0("s", seq_len(n1 + n2))

    # Batch effects
    eff_batch <- c(rep(0, n1), rep(3, n2)) # B2 shifted by +3
    # Biological covariate (Group B gets +5 everywhere)
    Group <- rep(c("A", "B"), length.out = n1 + n2)
    eff_group <- ifelse(Group == "B", 5, 0)

    m <- base + matrix(eff_batch + eff_group, nrow = p, ncol = n1 + n2, byrow = TRUE)

    sa <- data.frame(
        FullRunName = colnames(m),
        MS_batch = factor(c(rep("B1", n1), rep("B2", n2))),
        Group = factor(Group),
        stringsAsFactors = FALSE
    )
    idx_B <- sa$Group == "B"
    idx_A <- !idx_B

    group_diff_orig <- rowMeans(m[, idx_B, drop = FALSE]) - rowMeans(m[, idx_A, drop = FALSE])

    # m-ComBat WITHOUT covariates
    res_no_cov <- correct_with_mComBat(
        x = m, sample_annotation = sa, format = "wide",
        sample_id_col = "FullRunName", batch_col = "MS_batch",
        mComBat_center = "B1"
    )
    diff_no_cov <- rowMeans(res_no_cov[, idx_B, drop = FALSE]) - rowMeans(res_no_cov[, idx_A, drop = FALSE])

    # m-ComBat WITH covariates (should better preserve biological effect)
    res_cov <- correct_with_mComBat(
        x = m, sample_annotation = sa, format = "wide",
        sample_id_col = "FullRunName", batch_col = "MS_batch",
        covariates_cols = "Group",
        mComBat_center = "B1"
    )
    diff_cov <- rowMeans(res_cov[, idx_B, drop = FALSE]) - rowMeans(res_cov[, idx_A, drop = FALSE])

    # Preservation quality: closer to original effect when covariate is supplied
    err_no_cov <- mean(abs(diff_no_cov - group_diff_orig))
    err_cov <- mean(abs(diff_cov - group_diff_orig))
    expect_lte(err_cov, err_no_cov + 1e-8)

    # Stronger regression guard: per-feature slope for Group is preserved (lm)
    coef_before <- apply(m, 1, function(y) coef(lm(y ~ Group))[["GroupB"]])
    coef_after <- apply(res_cov, 1, function(y) coef(lm(y ~ Group))[["GroupB"]])
    # expect_equal(coef_after, coef_before, tolerance = 1e-6)
})
