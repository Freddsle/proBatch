expect_matrix_like <- function(actual, template, storage_mode = "double") {
    expect_true(is.matrix(actual))
    expect_identical(dim(actual), dim(template))
    expect_identical(dimnames(actual), dimnames(template))
    if (!is.null(storage_mode)) {
        expect_identical(storage.mode(actual), storage_mode)
    }
}

local_fake_omicsgmf_correct_step <- function(fake_step) {
    caller_env <- parent.frame()
    testthat::local_mocked_bindings(
        .pb_requireNamespace = function(...) invisible(TRUE),
        .omicsgmf_correct_matrix_step = fake_step,
        .package = "proBatch",
        .env = caller_env
    )
}

test_that("correct_with_omicsGMF(wide): forwards parameters and preserves dimnames (mocked)", {
    m <- matrix(
        as.numeric(1:9),
        nrow = 3,
        dimnames = list(
            paste0("feat", 1:3),
            paste0("sample", 1:3)
        )
    )

    sa <- data.frame(
        FullRunName = c("sample3", "sample1", "sample2"),
        Condition = c("C", "A", "B"),
        stringsAsFactors = FALSE
    )

    captured <- new.env(parent = emptyenv())
    fake_step <- function(data_matrix,
                          sample_annotation,
                          sample_id_col,
                          design_formula,
                          family,
                          ncomponents,
                          gmf_args = list(),
                          impute_args = list()) {
        captured$data_matrix <- data_matrix
        captured$sample_annotation <- sample_annotation
        captured$args <- list(
            sample_id_col  = sample_id_col,
            design_formula = design_formula,
            family         = family,
            ncomponents    = ncomponents,
            gmf_args       = gmf_args,
            impute_args    = impute_args
        )
        storage.mode(data_matrix) <- "double"
        data_matrix + 0.25
    }

    local_fake_omicsgmf_correct_step(fake_step)

    out <- correct_with_omicsGMF(
        x = m,
        sample_annotation = sa,
        sample_id_col = "FullRunName",
        design_formula = ~Condition,
        family = poisson(),
        ncomponents = 3L,
        gmf_args = list(name = "custom_dimred"),
        impute_args = list(name = "custom_imputed"),
        format = "wide"
    )

    expect_matrix_like(out, m)
    expect_equal(out, m + 0.25, tolerance = 1e-12)

    expect_identical(captured$data_matrix, m)
    expect_identical(
        captured$sample_annotation$FullRunName,
        sa$FullRunName
    )
    expect_identical(captured$args$sample_id_col, "FullRunName")
    expect_identical(captured$args$design_formula, stats::as.formula("~Condition"))
    expect_identical(captured$args$family$family, poisson()$family)
    expect_identical(captured$args$ncomponents, 3L)
    expect_identical(captured$args$gmf_args, list(name = "custom_dimred"))
    expect_identical(captured$args$impute_args, list(name = "custom_imputed"))
})

test_that("correct_with_omicsGMF(long): adds preBatchCorr_* and respects keep_all = 'minimal' (mocked)", {
    df <- data.frame(
        peptide_group_label = rep(c("pep1", "pep2"), each = 3),
        FullRunName = rep(paste0("sample", 1:3), times = 2),
        Intensity = as.numeric(1:6),
        stringsAsFactors = FALSE
    )

    sa <- data.frame(
        FullRunName = c("sample2", "sample1", "sample3"),
        Condition = c("B", "A", "C"),
        stringsAsFactors = FALSE
    )

    expected_matrix <- proBatch::long_to_matrix(
        df,
        feature_id_col = "peptide_group_label",
        sample_id_col = "FullRunName",
        measure_col = "Intensity"
    )

    captured <- new.env(parent = emptyenv())
    fake_step <- function(data_matrix,
                          sample_annotation,
                          sample_id_col,
                          design_formula,
                          family,
                          ncomponents,
                          gmf_args = list(),
                          impute_args = list()) {
        captured$data_matrix <- data_matrix
        captured$sample_annotation <- sample_annotation
        captured$args <- list(
            sample_id_col  = sample_id_col,
            design_formula = design_formula,
            family         = family,
            ncomponents    = ncomponents,
            gmf_args       = gmf_args,
            impute_args    = impute_args
        )
        storage.mode(data_matrix) <- "double"
        data_matrix + 1
    }

    local_fake_omicsgmf_correct_step(fake_step)

    out <- correct_with_omicsGMF(
        x = df,
        sample_annotation = sa,
        feature_id_col = "peptide_group_label",
        sample_id_col = "FullRunName",
        measure_col = "Intensity",
        design_formula = ~Condition,
        family = gaussian(),
        ncomponents = 2L,
        gmf_args = list(name = "gmf_name"),
        impute_args = list(name = "imputed_name"),
        keep_all = "minimal",
        format = "long"
    )

    expect_true(is.data.frame(out))
    expect_equal(
        sort(names(out)),
        sort(c("peptide_group_label", "FullRunName", "Intensity", "preBatchCorr_Intensity"))
    )

    idx <- match(
        interaction(df$peptide_group_label, df$FullRunName),
        interaction(out$peptide_group_label, out$FullRunName)
    )
    expect_equal(out$preBatchCorr_Intensity[idx], df$Intensity)
    expect_equal(out$Intensity[idx], df$Intensity + 1)

    expect_matrix_like(captured$data_matrix, expected_matrix)
    expect_equal(captured$data_matrix, expected_matrix, tolerance = 1e-12)
    expect_identical(
        captured$sample_annotation$FullRunName,
        colnames(expected_matrix)
    )
    expect_identical(captured$args$gmf_args, list(name = "gmf_name"))
    expect_identical(captured$args$impute_args, list(name = "imputed_name"))
})

test_that(".omicsgmf_correct_matrix_step reconstructs from GMF components (mocked fit)", {
    testthat::skip_if_not_installed("SingleCellExperiment")
    testthat::skip_if_not_installed("SummarizedExperiment")
    testthat::skip_if_not_installed("S4Vectors")

    m <- matrix(
        c(1, 4, 2, 5, 3, 6),
        nrow = 2,
        dimnames = list(
            c("f1", "f2"),
            c("s1", "s2", "s3")
        )
    )

    sa <- data.frame(
        FullRunName = c("s2", "s1", "s3"),
        stringsAsFactors = FALSE
    )

    gmf_results <- matrix(
        c(
            1, 2,
            3, 4,
            5, 6
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(NULL, c("comp1", "comp2"))
    )
    rotation <- matrix(
        c(
            0.5, 0.1,
            0.2, 0.4
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("f1", "f2"), c("comp1", "comp2"))
    )
    attr(gmf_results, "rotation") <- rotation

    fake_sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(dummy = matrix(0, nrow = nrow(m), ncol = ncol(m))),
        colData = S4Vectors::DataFrame(
            FullRunName = colnames(m)
        )
    )
    SingleCellExperiment::reducedDim(fake_sce, "GMF") <- gmf_results

    captured <- new.env(parent = emptyenv())
    fake_fit <- function(data_matrix,
                         sample_annotation,
                         design_formula,
                         family,
                         ncomponents,
                         gmf_args,
                         impute_args) {
        captured$data_matrix <- data_matrix
        captured$sample_annotation <- sample_annotation
        captured$args <- list(
            design_formula = design_formula,
            family = family,
            ncomponents = ncomponents,
            gmf_args = gmf_args,
            impute_args = impute_args
        )
        list(
            sce = fake_sce,
            dimred_name = "GMF",
            imputed = matrix(NA_real_, nrow = nrow(m), ncol = ncol(m)),
            imputed_assay = "omicsGMF_imputed"
        )
    }

    caller_env <- parent.frame()
    testthat::local_mocked_bindings(
        .pb_requireNamespace = function(...) invisible(TRUE),
        .omicsgmf_fit_and_impute = fake_fit,
        .package = "proBatch",
        .env = caller_env
    )

    out <- proBatch:::`.omicsgmf_correct_matrix_step`(
        data_matrix = m,
        sample_annotation = sa,
        sample_id_col = "FullRunName",
        design_formula = ~1,
        family = gaussian(),
        ncomponents = 2L,
        gmf_args = list(),
        impute_args = list()
    )

    expect_matrix_like(out, m)

    expected <- t(gmf_results %*% t(rotation))
    rownames(expected) <- rownames(m)
    colnames(expected) <- colnames(m)
    expect_equal(out, expected, tolerance = 1e-12)

    expect_matrix_like(captured$data_matrix, m)
    expect_equal(captured$data_matrix, m, tolerance = 1e-12)
    expect_identical(
        captured$sample_annotation$FullRunName,
        colnames(m)
    )
    expect_identical(captured$args$ncomponents, 2L)
    expect_identical(captured$args$gmf_args, list())
    expect_identical(captured$args$impute_args, list())
})

test_that("correct_with_omicsGMF validates format argument before dispatch (mocked)", {
    fake_step <- function(...) stop("should not be called")
    local_fake_omicsgmf_correct_step(fake_step)

    sa <- data.frame(
        FullRunName = c("s1", "s2"),
        stringsAsFactors = FALSE
    )

    expect_error(
        correct_with_omicsGMF(
            x = data.frame(a = 1:2),
            sample_annotation = sa,
            ncomponents = 2L,
            format = "wide"
        ),
        "requires a numeric matrix",
        ignore.case = TRUE
    )

    expect_error(
        correct_with_omicsGMF(
            x = matrix(1:4, nrow = 2),
            sample_annotation = sa,
            ncomponents = 2L,
            format = "long"
        ),
        "requires a data.frame",
        ignore.case = TRUE
    )
})
