# -------------------------
# Optional integration tests (mocked)
# -------------------------
expect_matrix_like <- function(actual, template, storage_mode = "double") {
    expect_true(is.matrix(actual))
    expect_identical(dim(actual), dim(template))
    expect_identical(dimnames(actual), dimnames(template))
    if (!is.null(storage_mode)) {
        expect_identical(storage.mode(actual), storage_mode)
    }
}
local_fake_ruviiic_matrix_step <- function(fake_core) {
    caller_env <- parent.frame()
    testthat::local_mocked_bindings(
        .pb_requireNamespace = function(pkg) invisible(TRUE),
        .ruviiic_matrix_step = fake_core,
        .package = "proBatch",
        .env = caller_env
    )
}
test_that("correct_with_RUVIII_C(wide) preserves dimnames and returns numeric matrix (mocked)", {
    m <- matrix(
        c(
            1, 4, 6, 8,
            2, 5, 7, 9
        ),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3", "s4"))
    )
    sa <- data.frame(
        FullRunName = c("s1", "s2", "s3", "s4"),
        replicate_id = c("r1", "r1", "r2", "r2"),
        stringsAsFactors = FALSE
    )
    fake_core <- function(data_matrix, ...) {
        storage.mode(data_matrix) <- "double"
        data_matrix + 5
    }
    local_fake_ruviiic_matrix_step(fake_core)
    out <- proBatch::correct_with_RUVIII_C(
        x = m,
        sample_annotation = sa,
        feature_id_col = "peptide_group_label",
        measure_col = "Intensity",
        sample_id_col = "FullRunName",
        replicate_col = "replicate_id",
        negative_control_features = "f1",
        k = 1L,
        format = "wide"
    )
    expect_matrix_like(out, m)
    expect_equal(out, m + 5)
})
test_that("correct_with_RUVIII_C(long) returns corrected data.frame when with_extra = FALSE (mocked)", {
    df <- data.frame(
        peptide_group_label = rep(c("p1", "p2"), each = 2),
        FullRunName = rep(c("s1", "s2"), times = 2),
        Intensity = as.numeric(1:4),
        stringsAsFactors = FALSE
    )
    sa <- data.frame(
        FullRunName = c("s1", "s2"),
        replicate_id = c("r1", "r2"),
        stringsAsFactors = FALSE
    )
    fake_core <- function(data_matrix, ...) {
        storage.mode(data_matrix) <- "double"
        data_matrix + 3
    }
    local_fake_ruviiic_matrix_step(fake_core)
    res <- proBatch::correct_with_RUVIII_C(
        x = df,
        sample_annotation = sa,
        feature_id_col = "peptide_group_label",
        measure_col = "Intensity",
        sample_id_col = "FullRunName",
        replicate_col = "replicate_id",
        negative_control_features = c("p1"),
        k = 1,
        format = "long"
    )
    expect_s3_class(res, "data.frame")
    expect_true("preBatchCorr_Intensity" %in% names(res))
    key <- interaction(res$peptide_group_label, res$FullRunName)
    expected_idx <- match(
        interaction(df$peptide_group_label, df$FullRunName),
        key
    )
    expect_equal(res$Intensity[expected_idx], df$Intensity + 3)
    expect_equal(res$preBatchCorr_Intensity[expected_idx], df$Intensity)
})
test_that("correct_with_RUVIII_C(long) adds preBatchCorr_* and returns list when with_extra = TRUE (mocked)", {
    df <- data.frame(
        peptide_group_label = rep(c("p1", "p2"), each = 4),
        FullRunName = rep(c("s1", "s2", "s3", "s4"), times = 2),
        Intensity = as.numeric(1:8),
        stringsAsFactors = FALSE
    )
    sa <- data.frame(
        FullRunName = c("s1", "s2", "s3", "s4"),
        replicate_id = c("r1", "r1", "r2", "r2"),
        stringsAsFactors = FALSE
    )
    fake_core <- function(data_matrix, ...) {
        list(
            newY = data_matrix + 2,
            residualDimensions = setNames(c(2L, 3L), rownames(data_matrix))
        )
    }
    local_fake_ruviiic_matrix_step(fake_core)
    res <- proBatch::correct_with_RUVIII_C(
        x = df,
        sample_annotation = sa,
        feature_id_col = "peptide_group_label",
        measure_col = "Intensity",
        sample_id_col = "FullRunName",
        replicate_col = "replicate_id",
        negative_control_features = c("p1"),
        k = 1,
        format = "long",
        with_extra = TRUE,
        keep_all = "minimal"
    )
    expect_true(is.list(res))
    expect_true("corrected_long" %in% names(res))
    expect_matrix_like(res$newY, proBatch::long_to_matrix(df, sample_id_col = "FullRunName"))
    expect_s3_class(res$corrected_long, "data.frame")
    minimal_cols <- c("peptide_group_label", "FullRunName", "Intensity", "preBatchCorr_Intensity")
    expect_true(all(minimal_cols %in% names(res$corrected_long)))
    # ensure correction applied (+2)
    idx <- match(
        interaction(df$peptide_group_label, df$FullRunName),
        interaction(res$corrected_long$peptide_group_label, res$corrected_long$FullRunName)
    )
    expect_equal(res$corrected_long$Intensity[idx], df$Intensity + 2)
    expect_equal(res$corrected_long$preBatchCorr_Intensity[idx], df$Intensity)
})
test_that(".ruviiic_matrix_step builds design matrix and transposes output (mocked core)", {
    captured <- new.env(parent = emptyenv())
    fake_core <- function(...) {
        args <- list(...)
        captured$args <- args
        Y <- args$Y
        matrix(
            seq_len(nrow(Y) * ncol(Y)),
            nrow = nrow(Y),
            ncol = ncol(Y),
            dimnames = dimnames(Y)
        )
    }
    testthat::local_mocked_bindings(
        .run_RUVIIIC_core = fake_core,
        .package = "proBatch"
    )
    m <- matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    sa <- data.frame(
        FullRunName = c("s2", "s1"),
        replicate_id = c("rB", "rA"),
        stringsAsFactors = FALSE
    )
    out <- proBatch:::`.ruviiic_matrix_step`(
        data_matrix = m,
        sample_annotation = sa,
        sample_id_col = "FullRunName",
        replicate_col = "replicate_id",
        negative_control_features = "f1",
        k = 1,
        version = "CPP"
    )
    expect_matrix_like(out, m)
    expect_equal(rownames(captured$args$Y), c("s1", "s2"))
    expect_equal(colnames(captured$args$Y), c("f1", "f2"))
    expect_equal(rownames(captured$args$M), c("s1", "s2"))
    expect_true(all(levels(factor(sa$replicate_id)) %in% colnames(captured$args$M)))
    expect_equal(captured$args$controls, "f1")
    expect_equal(captured$args$to_correct, rownames(m))
})
test_that(".ruviiic_matrix_step rejects missing replicate column", {
    m <- matrix(
        1:4,
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    sa <- data.frame(FullRunName = c("s1", "s2"), stringsAsFactors = FALSE)
    expect_error(
        proBatch:::`.ruviiic_matrix_step`(
            data_matrix = m,
            sample_annotation = sa,
            sample_id_col = "FullRunName",
            replicate_col = "replicate_id",
            negative_control_features = "f1",
            k = 1
        ),
        "replicate_col 'replicate_id' is not present",
        fixed = TRUE
    )
})
test_that(".ruviiic_matrix_step rejects replicate identifiers with missing values", {
    m <- matrix(
        1:4,
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    sa <- data.frame(
        FullRunName = c("s1", "s2"),
        replicate_id = c("r1", NA),
        stringsAsFactors = FALSE
    )
    expect_error(
        proBatch:::`.ruviiic_matrix_step`(
            data_matrix = m,
            sample_annotation = sa,
            sample_id_col = "FullRunName",
            replicate_col = "replicate_id",
            negative_control_features = "f1",
            k = 1
        ),
        "replicate_col contains missing values; provide replicate identifiers for all samples.",
        fixed = TRUE
    )
})
test_that(".ruviiic_matrix_step validates to_correct and controls against available features", {
    m <- matrix(
        1:4,
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    sa <- data.frame(
        FullRunName = c("s1", "s2"),
        replicate_id = c("r1", "r1"),
        stringsAsFactors = FALSE
    )
    expect_error(
        proBatch:::`.ruviiic_matrix_step`(
            data_matrix = m,
            sample_annotation = sa,
            sample_id_col = "FullRunName",
            replicate_col = "replicate_id",
            negative_control_features = "missing",
            k = 1
        ),
        "negative_control_features missing from data_matrix: missing",
        fixed = TRUE
    )
    expect_error(
        proBatch:::`.ruviiic_matrix_step`(
            data_matrix = m,
            sample_annotation = sa,
            sample_id_col = "FullRunName",
            replicate_col = "replicate_id",
            negative_control_features = "f1",
            to_correct = "ghost",
            k = 1
        ),
        "to_correct contains features absent from data_matrix: ghost",
        fixed = TRUE
    )
})
test_that(".check_ruviiic_inputs validates required arguments", {
    expect_error(
        proBatch:::`.check_ruviiic_inputs`(
            negative_control_features = "f1",
            k = 1
        ),
        "replicate_col must be provided and non-empty.",
        fixed = TRUE
    )
    expect_error(
        proBatch:::`.check_ruviiic_inputs`(
            replicate_col = "replicate_id",
            k = 1
        ),
        "negative_control_features must be a non-empty character vector.",
        fixed = TRUE
    )
    expect_error(
        proBatch:::`.check_ruviiic_inputs`(
            replicate_col = "replicate_id",
            negative_control_features = "f1",
            k = 0
        ),
        "k must be >= 1 for RUV-III-C correction.",
        fixed = TRUE
    )
})
