test_that("correct_with_BERT(wide): preserves dimnames and returns numeric matrix (mocked)", {
    testthat::skip_if_not_installed("BERT")
    m <- matrix(
        c(
            1, NA, 3,
            4, 5, 6
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
    )

    sa <- data.frame(
        FullRunName = c("s2", "s1", "s3"), # shuffled to test alignment
        MS_batch = c("B1", "B1", "B2"),
        stringsAsFactors = FALSE
    )

    # Fake core that adds +10 to all numeric cells and returns same shape
    fake_core <- function(data_matrix, ...) {
        storage.mode(data_matrix) <- "double"
        data_matrix + 10
    }

    local_mocked_bindings(
        .pb_requireNamespace = function(pkg) invisible(TRUE),
        .bert_matrix_step    = fake_core,
        .package             = "proBatch"
    )

    out <- proBatch::correct_with_BERT(
        x = m, sample_annotation = sa,
        sample_id_col = "FullRunName",
        batch_col = "MS_batch",
        format = "wide"
    )

    expect_true(is.matrix(out))
    expect_identical(dim(out), dim(m))
    expect_identical(dimnames(out), dimnames(m))
    expect_identical(storage.mode(out), "double")
    expect_equal(out, m + 10, tolerance = 1e-12)
})

test_that("correct_with_BERT(long): adds preBatchCorr_* and respects keep_all = 'minimal' (mocked)", {
    testthat::skip_if_not_installed("BERT")
    # Toy long with two features × three samples (complete)
    df <- data.frame(
        peptide_group_label = rep(c("p1", "p2"), each = 3),
        FullRunName         = rep(c("s1", "s2", "s3"), times = 2),
        Intensity           = as.numeric(1:6),
        stringsAsFactors    = FALSE
    )

    sa <- data.frame(
        FullRunName = c("s1", "s2", "s3"),
        MS_batch = c("B1", "B1", "B2"),
        stringsAsFactors = FALSE
    )

    # Core returns +1 everywhere
    fake_core <- function(data_matrix, ...) {
        storage.mode(data_matrix) <- "double"
        data_matrix + 1
    }

    local_mocked_bindings(
        .pb_requireNamespace = function(pkg) invisible(TRUE),
        .bert_matrix_step    = fake_core,
        .package             = "proBatch"
    )

    out <- correct_with_BERT(
        x = df, sample_annotation = sa,
        feature_id_col = "peptide_group_label",
        sample_id_col = "FullRunName",
        measure_col = "Intensity",
        batch_col = "MS_batch",
        format = "long",
        keep_all = "minimal"
    )
    expect_true(is.data.frame(out))

    # Minimal set of columns (order not enforced)
    expect_true(all(c(
        "peptide_group_label", "FullRunName", "Intensity",
        "preBatchCorr_Intensity"
    ) %in% names(out)))

    # preBatchCorr_* must match the original measure by (feature, sample)
    idx <- match(
        interaction(df$peptide_group_label, df$FullRunName),
        interaction(out$peptide_group_label, out$FullRunName)
    )
    expect_equal(out$preBatchCorr_Intensity[idx], df$Intensity)

    # Corrected measure is +1 vs original
    expect_equal(out$Intensity[idx], df$Intensity + 1)
})

test_that(".bert_matrix_step: rejects non-numeric and <2 features (mocked/no BERT)", {
    testthat::skip_if_not_installed("BERT")
    m_bad <- matrix(
        c("a", "b", "c", "d"),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    sa <- data.frame(FullRunName = c("s1", "s2"), MS_batch = c("B1", "B2"))

    expect_error(
        proBatch:::`.bert_matrix_step`(
            data_matrix = m_bad, sample_annotation = sa,
            sample_id_col = "FullRunName", batch_col = "MS_batch"
        ),
        "Input must be coercible to a numeric matrix for BERT correction",
        fixed = TRUE
    )

    m_tiny <- matrix(1:2, nrow = 1, dimnames = list("f1", c("s1", "s2")))
    expect_error(
        proBatch:::`.bert_matrix_step`(
            data_matrix = m_tiny, sample_annotation = sa,
            sample_id_col = "FullRunName", batch_col = "MS_batch"
        ),
        " at least two not-NA features.",
        ignore.case = TRUE
    )
})

test_that(".run_BERT_core: fails if batch column missing; aligns SA to matrix columns", {
    testthat::skip_if_not_installed("BERT")
    m <- matrix(
        c(
            1.0, 2.0, 3.0, 4.0,
            4.0, 5.0, 6.0, 7.0
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3", "s4"))
    )
    # Missing batch column -> error
    sa_bad <- data.frame(FullRunName = c("s3", "s2", "s1", "s4"))

    expect_error(
        proBatch:::`.run_BERT_core`(
            data_matrix = m, sample_annotation = sa_bad,
            sample_id_col = "FullRunName", batch_col = "MS_batch",
            bert_method = "limma"
        ),
        "Batch column is not present",
        fixed = TRUE
    )

    # Alignment check with shuffled SA
    sa <- data.frame(
        FullRunName = c("s3", "s2", "s1", "s4"), # shuffled on purpose
        MS_batch = c(1, 2, 1, 2),
        stringsAsFactors = FALSE
    )
    out <- proBatch:::`.run_BERT_core`(
        data_matrix = m, sample_annotation = sa,
        sample_id_col = "FullRunName", batch_col = "MS_batch",
        bert_method = "limma"
    )
    expect_true(is.matrix(out))
    expect_identical(rownames(out), rownames(m))
    expect_identical(colnames(out), colnames(m))
})

test_that("correct_with_BERT(long/wide): format argument checking and plumbing (mocked)", {
    testthat::skip_if_not_installed("BERT")
    sa <- data.frame(FullRunName = c("s1", "s2"), MS_batch = c("B1", "B2"))

    # format='wide' but x not a matrix -> error
    expect_error(
        proBatch::correct_with_BERT(
            x = data.frame(x = 1:2), sample_annotation = sa,
            format = "wide"
        ),
        "requires a numeric matrix",
        ignore.case = TRUE
    )

    # format='long' but x not a data.frame -> error
    expect_error(
        proBatch::correct_with_BERT(
            x = matrix(1, 1, 1), sample_annotation = sa,
            format = "long"
        ),
        "requires a data.frame",
        ignore.case = TRUE
    )
})

# -------------------------
# Optional integration test (only if BERT is available)
# -------------------------
test_that("Integration: .bert_matrix_step runs with BERT on tiny data (method='limma')", {
    testthat::skip_if_not_installed("BERT")

    m <- matrix(
        c(
            1, NA, 3, 4,
            4, 5, 6, 7
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3", "s4"))
    )
    sa <- data.frame(
        FullRunName = c("s1", "s2", "s3", "s4"),
        MS_batch = c("B1", "B1", "B2", "B2"),
        Condition = c("A", "A", "B", "B"), # covariate
        stringsAsFactors = FALSE
    )

    out <- proBatch:::`.bert_matrix_step`(
        data_matrix       = m,
        sample_annotation = sa,
        sample_id_col     = "FullRunName",
        batch_col         = "MS_batch",
        covariates_cols   = "Condition",
        bert_method       = "limma", # avoids ComBat priors & works with NA per BERT
        combatmode        = 1
    )

    expect_true(is.matrix(out))
    expect_identical(dim(out), dim(m))
    expect_identical(dimnames(out), dimnames(m))
    expect_equal(storage.mode(out), "double")
})
