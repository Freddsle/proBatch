# -------------------------
# Optional integration test (only if BERT is available)
# -------------------------

expect_matrix_like <- function(actual, template, storage_mode = "double") {
    expect_true(is.matrix(actual))
    expect_identical(dim(actual), dim(template))
    expect_identical(dimnames(actual), dimnames(template))
    if (!is.null(storage_mode)) {
        expect_identical(storage.mode(actual), storage_mode)
    }
}

local_fake_bert_core <- function(fake_core) {
    caller_env <- parent.frame()
    testthat::local_mocked_bindings(
        .pb_requireNamespace = function(pkg) invisible(TRUE),
        .bert_matrix_step    = fake_core,
        .package             = "proBatch",
        .env                 = caller_env
    )
}
test_that("correct_with_BERT(wide): preserves dimnames and returns numeric matrix (mocked)", {
    testthat::skip_if_not_installed("BERT")
    m <- matrix(
        c(
            1, NA, 3, 7,
            4, 5, 6, 8
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3", "s4"))
    )

    sa <- data.frame(
        FullRunName = c("s2", "s1", "s3", "s4"), # shuffled to test alignment
        MS_batch = c("B1", "B1", "B2", "B2"),
        stringsAsFactors = FALSE
    )

    # Fake core that adds +10 to all numeric cells and returns same shape
    fake_core <- function(data_matrix, ...) {
        storage.mode(data_matrix) <- "double"
        data_matrix + 10
    }

    local_fake_bert_core(fake_core)

    out <- proBatch::correct_with_BERT(
        x = m, sample_annotation = sa,
        sample_id_col = "FullRunName",
        batch_col = "MS_batch",
        format = "wide"
    )

    expect_matrix_like(out, m)
    expect_equal(out, m + 10, tolerance = 1e-12)
})

test_that(".run_BERT_core builds BERT input consistent with manual pipeline", {
    testthat::skip_if_not_installed("BERT")

    m <- matrix(
        c(
            1, 2, 3, 4,
            5, 6, 7, 8
        ),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s2", "s1", "s3", "s4"))
    )

    sa <- data.frame(
        ric_id = c("s1", "s2", "s3", "s4"),
        batch = c("B2", "B1", NA, "B10"),
        label = c("L1", "L2", "L1", NA),
        stringsAsFactors = FALSE
    )

    captured <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        BERT = function(data, method, combatmode = NULL, batchname, samplename,
                        covariatename, referencename, verify = NULL, ...) {
            captured$data <- data
            captured$args <- list(
                method = method,
                combatmode = combatmode,
                batchname = batchname,
                samplename = samplename,
                covariatename = covariatename,
                referencename = referencename,
                verify = verify
            )
            captured$dots <- list(...)
            data
        },
        .package = "BERT"
    )

    out <- proBatch:::.run_BERT_core(
        data_matrix = m,
        sample_annotation = sa,
        sample_id_col = "ric_id",
        batch_col = "batch",
        bert_method = "ComBat",
        labelname = "label",
        referencename = "REF"
    )

    sa_aligned <- proBatch:::.align_sample_annotation(
        sa,
        sample_ids = colnames(m),
        sample_id_col = "ric_id"
    )
    sa_aligned$Sample <- colnames(m)
    sa_aligned$batch <- as.numeric(as.factor(sa_aligned$batch))
    sa_aligned$label <- as.numeric(as.factor(sa_aligned$label))

    expected <- cbind(
        sa_aligned[, c("Sample", "batch", "label"), drop = FALSE],
        as.data.frame(t(m), check.names = TRUE)
    )
    rownames(expected) <- seq_len(nrow(expected))

    expect_identical(captured$data, expected)
    expect_false("ric_id" %in% names(captured$data))
    expect_identical(captured$args$batchname, "batch")
    expect_identical(captured$args$samplename, "Sample")
    expect_identical(captured$args$referencename, "REF")
    expect_identical(captured$dots$labelname, "label")
    expect_null(captured$dots$referencename)
    expect_identical(out, m)
})


# preBatchCorr_* must match the original measure by (feature, sample)


# Missing batch column -> error


# Alignment check with shuffled SA

# shuffled on purpose


# format='wide' but x not a matrix -> error


# format='long' but x not a data.frame -> error


# avoids ComBat priors & works with NA per BERT


# -------------------------
# Optional integration test (only if PLSDAbatch is available)
# -------------------------


test_that("correct_with_PLSDA_batch(wide): forwards args to PLSDAbatch (mocked)", {
    testthat::skip_if_not_installed("PLSDAbatch")

    m <- matrix(
        c(
            1, 2, 3, 4,
            5, 6, 7, 8,
            9, 10, 11, 12
        ),
        nrow = 3,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:4))
    )
    sa <- data.frame(
        FullRunName = c("s4", "s1", "s2", "s3"),
        MS_batch = c("B2", "B1", "B1", "B2"),
        Sex = c("F", "F", "F", "M"),
        stringsAsFactors = FALSE
    )

    captured <- new.env(parent = emptyenv())
    captured$calls <- list()

    testthat::local_mocked_bindings(
        .select_ncomp_trt = function(...) 2L,
        .select_keepX_trt = function(X, Y.trt, ncomp_trt, ...) rep(1L, ncomp_trt),
        .package = "proBatch"
    )

    testthat::local_mocked_bindings(
        PLSDA_batch = function(X, Y.trt, Y.bat,
                               ncomp.trt, ncomp.bat,
                               keepX.trt = NULL, keepX.bat = NULL,
                               max.iter = NULL, tol = NULL,
                               near.zero.var = NULL, balance = NULL, ...) {
            captured$calls[[length(captured$calls) + 1L]] <- list(
                X = X,
                Y.trt = Y.trt,
                Y.bat = Y.bat,
                ncomp.trt = ncomp.trt,
                ncomp.bat = ncomp.bat,
                keepX.trt = keepX.trt,
                keepX.bat = keepX.bat,
                max.iter = max.iter,
                tol = tol,
                near.zero.var = near.zero.var,
                balance = balance
            )
            list(
                X.nobatch = X,
                explained_variance.bat = list(
                    X = rep(0.5, ncomp.bat),
                    Y = rep(0.5, ncomp.bat)
                )
            )
        },
        .package = "PLSDAbatch"
    )

    # 1) basic wide run (auto-select ncomp_trt/bat; unbalanced=FALSE expected)
    out <- correct_with_PLSDA_batch(
        x = m, sample_annotation = sa,
        sample_id_col = "FullRunName",
        batch_col = "MS_batch",
        effect_col = "Sex",
        format = "wide"
    )
    expect_matrix_like(out, m)
    last_call <- tail(captured$calls, 1)[[1]]
    expect_false(last_call$balance)

    # 2) effect_col = NULL (should pass Y.trt = NULL into PLSDAbatch)
    out2 <- correct_with_PLSDA_batch(
        x = m, sample_annotation = sa,
        sample_id_col = "FullRunName",
        batch_col = "MS_batch",
        effect_col = NULL,
        format = "wide"
    )
    expect_matrix_like(out2, m)
    last_call <- tail(captured$calls, 1)[[1]]
    expect_null(last_call$Y.trt)

    # 3) sPLSDA path with tuning (keepX_trt inferred via tune.splsda)
    out3 <- correct_with_PLSDA_batch(
        x = m, sample_annotation = sa,
        sample_id_col = "FullRunName",
        batch_col = "MS_batch",
        effect_col = "Sex",
        ncomp_trt = 2L,
        run_splsda = TRUE,
        format = "wide"
    )
    expect_matrix_like(out3, m)
    last_call <- tail(captured$calls, 1)[[1]]
    expect_equal(last_call$keepX.trt, rep(1L, 2L))
})

test_that("correct_with_PLSDA_batch(wide): basic validation errors are informative", {
    testthat::skip_if_not_installed("PLSDAbatch")
    library(PLSDAbatch)
    m <- matrix(c(1, 2, NA, 4),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    sa <- data.frame(
        FullRunName = colnames(m),
        MS_batch = c("B1", "B2"),
        stringsAsFactors = FALSE
    )

    expect_error(
        correct_with_PLSDA_batch(
            x = m, sample_annotation = sa,
            sample_id_col = "FullRunName",
            batch_col = "MS_batch",
            effect_col = "Condition",
            format = "wide"
        ),
        "requires no NAs",
        ignore.case = TRUE
    )

    expect_error(
        correct_with_PLSDA_batch(
            x = as.data.frame(m), sample_annotation = sa,
            sample_id_col = "FullRunName",
            batch_col = "MS_batch",
            effect_col = "Condition",
            format = "wide"
        ),
        "requires a numeric matrix",
        ignore.case = TRUE
    )

    expect_error(
        correct_with_PLSDA_batch(
            x = matrix(c(1, 2, 3, 4), 2,
                dimnames = list(c("f1", "f2"), c("s1", "s2"))
            ),
            sample_annotation = sa[, "FullRunName", drop = FALSE], # drop batch col
            sample_id_col = "FullRunName",
            batch_col = "MS_batch",
            format = "wide"
        ),
        "Batch column 'MS_batch' is not present in sample_annotation",
        fixed = TRUE
    )
})

test_that("correct_with_PLSDA_batch(long): round-trip keeps rows and adds preBatchCorr_*", {
    testthat::skip_if_not_installed("PLSDAbatch")
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    # drop NA rows to avoid triggering validation errors
    m <- example_proteome
    m <- m[!apply(is.na(m), 1, any), ]
    df_long <- m
    sa <- example_sample_annotation

    library(PLSDAbatch)

    out <- correct_with_PLSDA_batch(
        x = df_long,
        sample_annotation = sa,
        sample_id_col = "FullRunName",
        feature_id_col = "peptide_group_label",
        measure_col = "Intensity",
        batch_col = "MS_batch",
        effect_col = "Sex",
        format = "long"
    )

    expect_true(is.data.frame(out))
    expect_true(all(c("peptide_group_label", "FullRunName", "Intensity") %in% names(out)))
    precol <- "preBatchCorr_Intensity"
    expect_true(precol %in% names(out))
    # row count is preserved (features * samples that were present)
    expect_equal(nrow(out), nrow(df_long))

    # minimal columns request
    out_min <- correct_with_PLSDA_batch(
        x = df_long,
        sample_annotation = sa,
        sample_id_col = "FullRunName",
        feature_id_col = "peptide_group_label",
        measure_col = "Intensity",
        batch_col = "MS_batch",
        effect_col = "Sex",
        format = "long",
        keep_all = "minimal"
    )
    expect_true(setequal(
        names(out_min),
        c("FullRunName", "peptide_group_label", "Intensity", precol)
    ))
})

test_that("correct_with_PLSDA_batch(wide, integration) runs with PLSDAbatch installed", {
    testthat::skip_if_not_installed("PLSDAbatch")

    m <- matrix(
        c(
            1, 2, 3, 4, 1,
            4, 5, 6, 7, 4
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("f1", "f2"), paste0("s", 1:5))
    )
    sa <- data.frame(
        FullRunName = colnames(m),
        MS_batch = c("B1", "B1", "B2", "B2", "B1"),
        Condition = c("A", "A", "B", "B", "A"),
        stringsAsFactors = FALSE
    )

    expect_warning(
        out <- correct_with_PLSDA_batch(
            x = m, sample_annotation = sa,
            sample_id_col = "FullRunName",
            batch_col = "MS_batch",
            effect_col = "Condition",
            format = "wide"
        ),
        "fewer than 3 samples",
        ignore.case = TRUE
    )
    expect_matrix_like(out, m)
})

test_that("correct_with_PLSDA_batch(long, integration) runs with PLSDAbatch installed", {
    testthat::skip_if_not_installed("PLSDAbatch")

    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    # drop NA rows to avoid triggering validation errors
    m <- example_proteome
    m <- m[!apply(is.na(m), 1, any), ]
    df_long <- m
    sa <- example_sample_annotation

    out <- correct_with_PLSDA_batch(
        x = df_long, sample_annotation = sa,
        sample_id_col = "FullRunName",
        feature_id_col = "peptide_group_label",
        measure_col = "Intensity",
        batch_col = "MS_batch",
        effect_col = "Sex",
        format = "long"
    )
    expect_true(is.data.frame(out))
    expect_true(all(c("peptide_group_label", "FullRunName", "Intensity", "preBatchCorr_Intensity") %in% names(out)))
})

test_that("ProBatchFeatures smoke-test (skipped if class/methods unavailable)", {
    testthat::skip_if_not_installed("PLSDAbatch")

    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")
    # drop NA rows to avoid triggering validation errors
    m <- example_proteome_matrix
    m <- m[!apply(is.na(m), 1, any), ]

    pbf <- suppressMessages(
        ProBatchFeatures(
            data_matrix = m,
            sample_annotation = example_sample_annotation,
            sample_id_col = "FullRunName",
            name = "protein::raw"
        )
    )

    params <- list(
        batch_col = "MS_batch",
        sample_id_col = "FullRunName",
        effect_col = "Sex"
    )

    pb_corrected <- suppressMessages(pb_transform(
        pbf,
        from = "protein::raw",
        steps = "PLSDAbatch",
        params_list = list(params)
    ))

    expect_true(inherits(pb_corrected, "ProBatchFeatures"))
    expect_true("protein::PLSDAbatch_on_raw" %in% names(pb_corrected))

    corrected <- suppressMessages(pb_assay_matrix(pb_corrected, "protein::PLSDAbatch_on_raw"))
    expect_matrix_like(corrected, m)
})


test_that("correct_with_PLSDAbatch(wide)", {
    testthat::skip_if_not_installed("PLSDAbatch")
    m <- matrix(
        c(
            1, 2, 3, 4, 1,
            4, 5, 6, 7, 4
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(
            c("f1", "f2"),
            c("s1", "s2", "s3", "s4", "s5")
        )
    )

    sa <- data.frame(
        FullRunName = c("s1", "s2", "s3", "s4", "s5"),
        MS_batch = c("B1", "B1", "B2", "B2", "B1"),
        Condition = c("A", "A", "B", "B", "A"),
        stringsAsFactors = FALSE
    )

    expect_warning(
        out <- correct_with_PLSDA_batch(
            x = m, sample_annotation = sa,
            sample_id_col = "FullRunName",
            batch_col = "MS_batch",
            format = "wide",
            effect_col = "Condition",
        ),
        "Some treatment levels have fewer than 3 samples"
    )

    expect_true(is.matrix(out))
    expect_identical(dim(out), dim(m))
    expect_identical(dimnames(out), dimnames(m))
    expect_identical(storage.mode(out), "double")
})
