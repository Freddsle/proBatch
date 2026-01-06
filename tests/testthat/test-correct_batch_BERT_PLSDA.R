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


test_that(, {
    # drop NA rows to avoid triggering validation errors


    # Capture arguments reaching PLSDAbatch::PLSDA_batch


    # samples x features
    # record what the wrapper passed


    # minimal structure used by .select_ncomp_bat()


    # 1) basic wide run (auto-select ncomp_trt/bat; unbalanced=FALSE expected)
    out <- correct_with_PLSDA_batch(,
        batch_col = "MS_batch",
        effect_col = "Sex",
        format = "wide"
    )

    expect_matrix_like(out, m)

    # 2) effect_col = NULL (should pass Y.trt = NULL into PLSDAbatch)
    out2 <- correct_with_PLSDA_batch(
        x = m, sample_annotation = sa,
        sample_id_col = "FullRunName",
        batch_col = "MS_batch",
        effect_col = NULL,
        format = "wide"
    )
    expect_matrix_like(out2, m)

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
