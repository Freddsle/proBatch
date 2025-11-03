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
