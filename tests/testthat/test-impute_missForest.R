expect_matrix_like <- function(actual, template, storage_mode = "double") {
    testthat::expect_true(is.matrix(actual))
    testthat::expect_identical(dim(actual), dim(template))
    testthat::expect_identical(dimnames(actual), dimnames(template))
    if (!is.null(storage_mode)) {
        testthat::expect_identical(storage.mode(actual), storage_mode)
    }
}

## helper to mock the internal matrix step, like local_fake_omicsgmf_step() did
local_fake_missforest_step <- function(fake_step) {
    caller_env <- parent.frame()
    testthat::local_mocked_bindings(
        .pb_requireNamespace = function(...) invisible(TRUE),
        .missforest_matrix_step = fake_step,
        .package = "proBatch",
        .env = caller_env
    )
}

## we also need the ProBatchFeatures helper, copied from your omicsGMF test
make_test_pbf <- function() {
    testthat::skip_if_not_installed("QFeatures")
    testthat::skip_if_not_installed("SummarizedExperiment")
    testthat::skip_if_not_installed("S4Vectors")

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            raw = matrix(
                c(
                    1, NA, 3,
                    4, 5, NA
                ),
                nrow = 2,
                byrow = TRUE,
                dimnames = list(
                    c("feat1", "feat2"),
                    c("s1", "s2", "s3")
                )
            )
        ),
        colData = S4Vectors::DataFrame(
            FullRunName = c("s1", "s2", "s3"),
            Batch       = c("A", "A", "B")
        )
    )

    qf <- QFeatures::QFeatures(list(raw = se))

    testthat::skip_if_not("as_ProBatchFeatures" %in% getNamespaceExports("proBatch"))
    proBatch::as_ProBatchFeatures(qf)
}

testthat::test_that("imputeMissForest(wide): forwards params and preserves dimnames (mocked)", {
    m <- matrix(
        as.numeric(1:9),
        nrow = 3,
        dimnames = list(
            paste0("feat", 1:3),
            paste0("sample", 1:3)
        )
    )

    captured <- new.env(parent = emptyenv())
    fake_step <- function(data_matrix, missforest_args = list()) {
        captured$data_matrix <- data_matrix
        captured$missforest_args <- missforest_args
        storage.mode(data_matrix) <- "double"
        data_matrix + 1
    }

    local_fake_missforest_step(fake_step)

    out <- proBatch::imputeMissForest(
        x = m,
        maxiter = 3L,
        ntree = 250L
    )

    expect_matrix_like(out, m)
    testthat::expect_equal(out, m + 1, tolerance = 1e-12)
    testthat::expect_identical(captured$data_matrix, m)
    ## missforest_args should contain exactly what we passed through ...
    testthat::expect_identical(captured$missforest_args$maxiter, 3L)
    testthat::expect_identical(captured$missforest_args$ntree, 250L)
})

testthat::test_that("imputeMissForest(long): roundtrips long -> matrix -> long and keeps IDs (mocked)", {
    df <- data.frame(
        feature_id = rep(c("pep1", "pep2"), each = 3),
        FullRunName = rep(paste0("s", 1:3), times = 2),
        Intensity = c(1, NA, 3, 4, 5, NA),
        stringsAsFactors = FALSE
    )

    sa <- data.frame(
        FullRunName = c("s2", "s1", "s3"),
        Condition = c("B", "A", "C"),
        stringsAsFactors = FALSE
    )

    ## what long_to_matrix() would produce
    expected_matrix <- proBatch::long_to_matrix(
        df_long        = df,
        feature_id_col = "feature_id",
        sample_id_col  = "FullRunName",
        measure_col    = "Intensity"
    )

    captured <- new.env(parent = emptyenv())
    fake_step <- function(data_matrix, missforest_args = list()) {
        captured$data_matrix <- data_matrix
        captured$missforest_args <- missforest_args
        storage.mode(data_matrix) <- "double"
        data_matrix + 10
    }

    local_fake_missforest_step(fake_step)

    out <- proBatch::imputeMissForest(
        x = df,
        sample_annotation = sa,
        sample_id_col = "FullRunName",
        feature_id_col = "feature_id",
        measure_col = "Intensity",
        parallelize = "no"
    )

    testthat::expect_true(is.data.frame(out))
    ## should at least contain these 3 core columns
    testthat::expect_true(all(c("feature_id", "FullRunName", "Intensity") %in% names(out)))
    ## no rows should be lost
    testthat::expect_equal(nrow(out), nrow(df))

    ## match rows by (feature_id, FullRunName)
    key_in <- interaction(df$feature_id, df$FullRunName, drop = TRUE)
    key_out <- interaction(out$feature_id, out$FullRunName, drop = TRUE)

    idx <- match(key_in, key_out)
    testthat::expect_equal(out$Intensity[idx], df$Intensity + 10)

    ## we actually called the matrix step with the matrix we expected
    expect_matrix_like(captured$data_matrix, expected_matrix)
    ## and we forwarded extra args
    testthat::expect_identical(captured$missforest_args$parallelize, "no")
})

testthat::test_that("imputeMissForest(ProBatchFeatures) calls .pb_apply_step with correct params", {
    pbf <- make_test_pbf()

    captured <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        .pb_requireNamespace = function(...) invisible(TRUE),
        .pb_apply_step = function(object, from, step, fun, params) {
            captured$object <- object
            captured$from <- from
            captured$step <- step
            captured$fun <- fun
            captured$params <- params
            ## mimic the real helper: return updated object + assay name
            list(object = object, assay = "missForestImpute_1")
        },
        .package = "proBatch"
    )

    out <- proBatch::imputeMissForest(
        x = pbf,
        pbf_name = NULL,
        maxiter = 2L
    )

    testthat::expect_s4_class(out, "ProBatchFeatures")
    testthat::expect_identical(captured$fun, "missForestImpute")
    ## missForest args forwarded via missforest_args bundle
    testthat::expect_identical(
        captured$params,
        list(missforest_args = list(maxiter = 2L))
    )
})

testthat::test_that("imputeMissForest(ProBatchFeatures) merges missforest_args with dots", {
    pbf <- make_test_pbf()

    captured <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        .pb_requireNamespace = function(...) invisible(TRUE),
        .pb_apply_step = function(object, from, step, fun, params) {
            captured$params <- params
            list(object = object, assay = "missForestImpute_1")
        },
        .package = "proBatch"
    )

    out <- proBatch::imputeMissForest(
        x = pbf,
        missforest_args = list(maxiter = 1L, mtry = 4L),
        ntree = 200L
    )

    testthat::expect_s4_class(out, "ProBatchFeatures")
    args <- captured$params$missforest_args
    testthat::expect_identical(args$maxiter, 1L)
    testthat::expect_identical(args$mtry, 4L)
    testthat::expect_identical(args$ntree, 200L)
})

testthat::test_that(".missforest_matrix_step() actually runs missForest when installed", {
    testthat::skip_if_not_installed("missForest")

    m <- matrix(
        c(
            1, NA, 3,
            4, 5, NA,
            NA, NA, NA
        ), # <- all-NA row to test the remove/reinsert branch
        nrow = 3,
        byrow = TRUE,
        dimnames = list(
            c("f1", "f2", "f_allNA"),
            c("s1", "s2", "s3")
        )
    )

    res <- proBatch:::`.missforest_matrix_step`(m)

    ## original dim + names are back
    expect_matrix_like(res, m)
    ## row that was all NA must stay NA (we decided to keep it NA)
    testthat::expect_true(all(is.na(res["f_allNA", ])))
    ## other rows should have no NA if missForest could impute them
    testthat::expect_false(any(is.na(res[c("f1", "f2"), ])))
    ## diagnostics from missForest should be there
    testthat::expect_true(!is.null(attr(res, "missForest_OOBerror")))
})

testthat::test_that(".missforest_matrix_step() errors when columns have no names", {
    testthat::skip_if_not_installed("missForest")

    m <- matrix(c(1, NA, 3, 4), nrow = 2)
    testthat::expect_error(
        proBatch:::`.missforest_matrix_step`(m),
        "requires matrix column names",
        fixed = TRUE
    )
})
