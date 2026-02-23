make_prone_test_pbf <- function() {
    testthat::skip_if_not_installed("QFeatures")
    testthat::skip_if_not_installed("SummarizedExperiment")
    testthat::skip_if_not_installed("S4Vectors")

    dm <- matrix(
        c(
            10, NA, 30,
            40, 50, NA
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(
            c("feat1", "feat2"),
            paste0("sample", 1:3)
        )
    )

    sa <- data.frame(
        FullRunName = paste0("sample", 1:3),
        Condition = c("ctrl", "case", "ctrl"),
        stringsAsFactors = FALSE
    )

    ProBatchFeatures(
        data_matrix = dm,
        sample_annotation = sa,
        sample_id_col = "FullRunName",
        name = "raw"
    )
}

local_mocked_prone <- function(mock_fun) {
    caller_env <- parent.frame()
    prev_opt <- getOption("proBatch.prone_impute_se", NULL)
    withr::defer(
        options(proBatch.prone_impute_se = prev_opt),
        envir = caller_env
    )
    options(proBatch.prone_impute_se = mock_fun)
    testthat::local_mocked_bindings(
        impute_se = mock_fun,
        .package = "PRONE",
        .env = caller_env
    )
}

test_that("pb_transform with PRONEImpute forwards condition column name", {
    skip_if_not_installed("PRONE")

    pbf <- make_prone_test_pbf()
    expected_condition <- as.character(SummarizedExperiment::colData(pbf)$Condition)
    expected_samples <- colnames(pb_assay_matrix(pbf))

    captured <- new.env(parent = emptyenv())
    local_mocked_prone(function(se, ain, condition) {
        captured$condition <- condition
        captured$colData <- as.data.frame(SummarizedExperiment::colData(se))
        captured$assay <- SummarizedExperiment::assay(se, ain)
        se
    })

    out <- pb_transform(
        object = pbf,
        from = "feature::raw",
        steps = "PRONEImpute",
        params_list = list(list(condition_col = "Condition"))
    )

    expect_s4_class(out, "ProBatchFeatures")
    expect_equal(captured$condition, "Condition")
    expect_equal(
        rownames(captured$colData),
        expected_samples
    )
    expect_equal(
        as.character(captured$colData$Condition),
        expected_condition
    )
    expect_identical(
        captured$assay,
        pb_assay_matrix(pbf),
        ignore_attr = TRUE
    )
    expect_false(any(grepl("^\\.pb_prone_condition", names(captured$colData))))
})

test_that("pb_transform with PRONEImpute accepts condition vectors", {
    skip_if_not_installed("PRONE")

    pbf <- make_prone_test_pbf()
    sample_ids <- colnames(pb_assay_matrix(pbf))
    custom_condition <- setNames(c("grpA", "grpB", "grpA"), sample_ids)

    captured <- new.env(parent = emptyenv())
    local_mocked_prone(function(se, ain, condition) {
        captured$condition <- condition
        captured$colData <- as.data.frame(SummarizedExperiment::colData(se))
        captured$assay <- SummarizedExperiment::assay(se, ain)
        se
    })

    out <- pb_transform(
        object = pbf,
        from = "feature::raw",
        steps = "PRONEImpute",
        params_list = list(list(condition_col = custom_condition))
    )

    expect_s4_class(out, "ProBatchFeatures")
    expect_true(grepl("^\\.pb_prone_condition", captured$condition))
    expect_true(captured$condition %in% names(captured$colData))
    expect_equal(
        captured$colData[[captured$condition]],
        unname(custom_condition[rownames(captured$colData)])
    )
    expect_identical(
        captured$assay,
        pb_assay_matrix(pbf),
        ignore_attr = TRUE
    )
})

test_that("imputePRONE_dm returns the imputed assay matrix", {
    skip_if_not_installed("PRONE")
    skip_if_not_installed("SummarizedExperiment")

    dm <- matrix(
        c(1, NA, 3, 4),
        nrow = 2,
        dimnames = list(
            c("featA", "featB"),
            c("sample1", "sample2")
        )
    )

    sa <- data.frame(
        FullRunName = colnames(dm),
        stringsAsFactors = FALSE
    )

    expected <- dm
    expected["featB", "sample1"] <- 2

    local_mocked_prone(function(se, ain, condition) {
        se_out <- se
        SummarizedExperiment::assay(
            se_out,
            paste0(ain, "_imputed"),
            withDimnames = FALSE
        ) <- expected
        se_out
    })

    res <- imputePRONE_dm(
        dm,
        sample_annotation = sa,
        sample_id_col = "FullRunName",
        assay_in = "raw"
    )

    expect_equal(
        res,
        expected,
        ignore_attr = TRUE
    )
    expect_false(anyNA(res))
    expect_identical(storage.mode(res), "double")
})

test_that(".pb_prone_restore_dimnames realigns matrices by feature/sample IDs", {
    original <- matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("feat1", "feat2"), c("sample1", "sample2"))
    )

    scrambled <- original[c("feat2", "feat1"), c("sample2", "sample1"), drop = FALSE]

    restored <- proBatch:::.pb_prone_restore_dimnames(
        original_matrix = original,
        imputed_matrix = scrambled
    )

    expect_identical(rownames(restored), rownames(original))
    expect_identical(colnames(restored), colnames(original))
    expect_equal(restored, original)
})

test_that(".pb_prone_restore_dimnames pads missing features with NA", {
    original <- matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("feat1", "feat2"), c("sample1", "sample2"))
    )

    partial <- original["feat2", , drop = FALSE]

    restored <- proBatch:::.pb_prone_restore_dimnames(
        original_matrix = original,
        imputed_matrix = partial
    )

    expect_identical(rownames(restored), rownames(original))
    expect_identical(colnames(restored), colnames(original))
    expect_true(all(is.na(restored["feat1", ])))
    expect_equal(restored["feat2", ], original["feat2", ])
})

test_that("imputePRONE_dm restores dimnames when PRONE output drops them", {
    skip_if_not_installed("PRONE")
    skip_if_not_installed("SummarizedExperiment")

    dm <- matrix(
        c(1, NA, 3, 4),
        nrow = 2,
        dimnames = list(
            c("featA", "featB"),
            c("sample1", "sample2")
        )
    )

    sa <- data.frame(
        FullRunName = colnames(dm),
        stringsAsFactors = FALSE
    )

    expected <- dm
    expected["featB", "sample1"] <- 2
    no_dimnames <- unname(expected)

    local_mocked_prone(function(se, ain, condition) {
        se_out <- se
        SummarizedExperiment::assay(
            se_out,
            paste0(ain, "_imputed"),
            withDimnames = FALSE
        ) <- no_dimnames
        se_out
    })

    res <- imputePRONE_dm(
        dm,
        sample_annotation = sa,
        sample_id_col = "FullRunName",
        assay_in = "raw"
    )

    expect_identical(rownames(res), rownames(dm))
    expect_identical(colnames(res), colnames(dm))
    expect_equal(res, expected, ignore_attr = TRUE)
})

test_that("imputePRONE_dm pads missing feature rows returned by PRONE", {
    skip_if_not_installed("PRONE")
    skip_if_not_installed("SummarizedExperiment")

    dm <- matrix(
        c(1, NA, 3, 4),
        nrow = 2,
        dimnames = list(
            c("featA", "featB"),
            c("sample1", "sample2")
        )
    )

    sa <- data.frame(
        FullRunName = colnames(dm),
        stringsAsFactors = FALSE
    )

    local_mocked_prone(function(se, ain, condition) {
        # Return a valid SummarizedExperiment with fewer rows (features).
        # This mimics an imputer that dropped some features.
        se_out <- se["featB", , drop = FALSE]
        SummarizedExperiment::assay(
            se_out,
            paste0(ain, "_imputed"),
            withDimnames = FALSE
        ) <- matrix(c(10, 20), nrow = 1)
        se_out
    })

    res <- imputePRONE_dm(
        dm,
        sample_annotation = sa,
        sample_id_col = "FullRunName",
        assay_in = "raw"
    )

    expect_identical(rownames(res), rownames(dm))
    expect_identical(colnames(res), colnames(dm))
    expect_true(all(is.na(res["featA", ])))
    expect_equal(res["featB", ], c(sample1 = 10, sample2 = 20))
})
