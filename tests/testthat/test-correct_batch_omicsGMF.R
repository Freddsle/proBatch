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
                          batch_col = NULL,
                          family,
                          ncomponents,
                          gmf_args = list(),
                          impute_args = list()) {
        captured$data_matrix <- data_matrix
        captured$sample_annotation <- sample_annotation
        captured$args <- list(
            sample_id_col  = sample_id_col,
            design_formula = design_formula,
            batch_col      = batch_col,
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
    expect_null(captured$args$batch_col)
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
                          batch_col = NULL,
                          family,
                          ncomponents,
                          gmf_args = list(),
                          impute_args = list()) {
        captured$data_matrix <- data_matrix
        captured$sample_annotation <- sample_annotation
        captured$args <- list(
            sample_id_col  = sample_id_col,
            design_formula = design_formula,
            batch_col      = batch_col,
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
    expect_null(captured$args$batch_col)
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
    # Strip omicsGMF_* latent attributes that are now attached for downstream
    # consumers; this test only verifies numeric reconstruction.
    out_numeric <- out
    for (a in grep("^omicsGMF_", names(attributes(out_numeric)), value = TRUE)) {
        attr(out_numeric, a) <- NULL
    }
    expect_equal(out_numeric, expected, tolerance = 1e-12)

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

test_that(".omicsgmf_correct_matrix_step preserves non-batch design terms when batch_col is provided", {
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
        MS_batch = c("b2", "b1", "b2"),
        Diet = c("A", "A", "B"),
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

    X <- cbind(
        "(Intercept)" = c(1, 1, 1),
        "DietB" = c(0, 0, 1),
        "MS_batchb2" = c(1, 0, 1),
        "DietB:MS_batchb2" = c(0, 0, 1)
    )
    # omicsGMF convention: Beta has one row per feature and one column per
    # design term, so the modelled fixed-effect mean is X %*% t(Beta).
    Beta <- matrix(
        c(
            1.0, 0.1, -0.3, 0.7,
            0.5, 0.2, 0.4, -0.6
        ),
        nrow = 2,
        byrow = TRUE
    )
    attr(gmf_results, "X") <- X
    attr(gmf_results, "Beta") <- Beta

    fake_sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(dummy = matrix(0, nrow = nrow(m), ncol = ncol(m))),
        colData = S4Vectors::DataFrame(
            FullRunName = colnames(m)
        )
    )
    SingleCellExperiment::reducedDim(fake_sce, "GMF") <- gmf_results

    fake_fit <- function(data_matrix,
                         sample_annotation,
                         design_formula,
                         family,
                         ncomponents,
                         gmf_args,
                         impute_args) {
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
        design_formula = ~ MS_batch + Diet,
        batch_col = "MS_batch",
        family = gaussian(),
        ncomponents = 2L,
        gmf_args = list(),
        impute_args = list()
    )

    expect_matrix_like(out, m)

    # New semantics: subtract only the batch-attributable part of the modelled
    # mean (X[, batch] %*% t(Beta[, batch])) from the observed data.
    batch_idx <- which(colnames(X) %in% c("MS_batchb2", "DietB:MS_batchb2"))
    batch_mean <- X[, batch_idx, drop = FALSE] %*%
        t(Beta[, batch_idx, drop = FALSE])
    expected <- m - t(batch_mean)
    rownames(expected) <- rownames(m)
    colnames(expected) <- colnames(m)

    # Strip omicsGMF_* latent attributes that are now attached for downstream
    # consumers; this test only verifies numeric reconstruction.
    out_numeric <- out
    for (a in grep("^omicsGMF_", names(attributes(out_numeric)), value = TRUE)) {
        attr(out_numeric, a) <- NULL
    }
    expect_equal(out_numeric, expected, tolerance = 1e-12)
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


# -------------------------
# Regression tests: omicsGMFcor fill_the_missing pass-through
# -------------------------

test_that(".omicsgmf_correct_matrix_step accepts and forwards fill_the_missing", {
    m <- matrix(
        as.numeric(1:6),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
    )
    sa <- data.frame(
        FullRunName = c("s1", "s2", "s3"),
        MS_batch = c("B1", "B1", "B2"),
        stringsAsFactors = FALSE
    )

    captured <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        .pb_require_omicsgmf_stack = function(...) invisible(TRUE),
        .run_matrix_method = function(data_matrix, sample_annotation, sample_id_col,
                                      fill_the_missing, missing_warning, method_fun, ...) {
            captured$fill_the_missing <- fill_the_missing
            captured$data_matrix <- data_matrix
            data_matrix
        },
        .package = "proBatch"
    )

    # Default — fill_the_missing = NULL
    proBatch:::.omicsgmf_correct_matrix_step(
        data_matrix = m, sample_annotation = sa,
        sample_id_col = "FullRunName", batch_col = "MS_batch",
        ncomponents = 2L
    )
    expect_null(captured$fill_the_missing)

    # Explicit FALSE — used to raise "unused argument" before the fix
    proBatch:::.omicsgmf_correct_matrix_step(
        data_matrix = m, sample_annotation = sa,
        sample_id_col = "FullRunName", batch_col = "MS_batch",
        ncomponents = 2L,
        fill_the_missing = FALSE
    )
    expect_identical(captured$fill_the_missing, FALSE)

    # "remove" policy is forwarded as-is
    proBatch:::.omicsgmf_correct_matrix_step(
        data_matrix = m, sample_annotation = sa,
        sample_id_col = "FullRunName", batch_col = "MS_batch",
        ncomponents = 2L,
        fill_the_missing = "remove"
    )
    expect_identical(captured$fill_the_missing, "remove")
})

test_that("correct_with_omicsGMF forwards fill_the_missing to the matrix step", {
    m <- matrix(
        as.numeric(1:9),
        nrow = 3,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:3))
    )
    sa <- data.frame(
        FullRunName = paste0("s", 1:3),
        MS_batch = c("B1", "B1", "B2"),
        stringsAsFactors = FALSE
    )

    captured <- new.env(parent = emptyenv())
    fake_step <- function(data_matrix, sample_annotation, sample_id_col,
                          design_formula, batch_col = NULL, family,
                          ncomponents, gmf_args = list(), impute_args = list(),
                          fill_the_missing = NULL, ...) {
        captured$fill_the_missing <- fill_the_missing
        data_matrix
    }
    local_fake_omicsgmf_correct_step(fake_step)

    correct_with_omicsGMF(
        x = m, sample_annotation = sa,
        sample_id_col = "FullRunName",
        design_formula = ~1,
        batch_col = "MS_batch",
        ncomponents = 2L,
        fill_the_missing = FALSE,
        format = "wide"
    )
    expect_identical(captured$fill_the_missing, FALSE)
})

# -------------------------
# Latent attributes are attached to the corrected matrix
# -------------------------

test_that(".omicsgmf_reconstruct_corrected_matrix attaches latent attributes (latent-only path)", {
    gmf_results <- matrix(
        c(1, 0, 0, 1, 0.5, 0.5),
        nrow = 3, ncol = 2,
        dimnames = list(NULL, c("comp1", "comp2"))
    )
    rotation <- matrix(
        c(0.6, 0.4, 0.3, 0.7),
        nrow = 2, ncol = 2,
        dimnames = list(c("f1", "f2"), c("comp1", "comp2"))
    )
    attr(gmf_results, "rotation") <- rotation

    data_matrix <- matrix(
        c(1, 4, 2, 5, 3, 6),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
    )

    out <- proBatch:::`.omicsgmf_reconstruct_corrected_matrix`(
        gmf_results = gmf_results,
        data_matrix = data_matrix,
        batch_col = NULL
    )

    expect_true(is.matrix(out))
    expect_identical(dim(out), dim(data_matrix))

    expect_true(!is.null(attr(out, "omicsGMF_scores")),
        info = "omicsGMF_scores attribute must be present"
    )
    expect_true(!is.null(attr(out, "omicsGMF_loadings")),
        info = "omicsGMF_loadings attribute must be present"
    )
    expect_identical(attr(out, "omicsGMF_dimred_name"), "omicsGMF")

    scores <- attr(out, "omicsGMF_scores")
    loadings <- attr(out, "omicsGMF_loadings")
    expect_identical(dim(scores), dim(gmf_results))
    expect_identical(dim(loadings), dim(rotation))
})

test_that(".omicsgmf_reconstruct_corrected_matrix attaches latent attributes (batch-subtracted path)", {
    n_samples <- 3L
    n_features <- 2L
    n_comp <- 2L

    gmf_results <- matrix(
        c(1, 0, 0, 1, 0.5, 0.5),
        nrow = n_samples, ncol = n_comp,
        dimnames = list(NULL, c("comp1", "comp2"))
    )
    rotation <- matrix(
        c(0.6, 0.4, 0.3, 0.7),
        nrow = n_features, ncol = n_comp,
        dimnames = list(c("f1", "f2"), c("comp1", "comp2"))
    )
    attr(gmf_results, "rotation") <- rotation

    X <- matrix(
        c(1, 1, 1, 0, 1, 0),
        nrow = n_samples, ncol = 2,
        dimnames = list(NULL, c("(Intercept)", "MS_batchB2"))
    )
    Beta <- matrix(
        c(1, 0.5, 0.2, 0.8),
        nrow = n_features, ncol = 2,
        dimnames = list(c("f1", "f2"), c("(Intercept)", "MS_batchB2"))
    )
    attr(gmf_results, "X") <- X
    attr(gmf_results, "Beta") <- Beta

    data_matrix <- matrix(
        c(1, 4, 2, 5, 3, 6),
        nrow = n_features,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
    )

    out <- proBatch:::`.omicsgmf_reconstruct_corrected_matrix`(
        gmf_results = gmf_results,
        data_matrix = data_matrix,
        batch_col = "MS_batch"
    )

    expect_true(is.matrix(out))
    expect_true(!is.null(attr(out, "omicsGMF_scores")))
    expect_true(!is.null(attr(out, "omicsGMF_loadings")))
    expect_identical(dim(attr(out, "omicsGMF_scores")), c(n_samples, n_comp))
    expect_identical(dim(attr(out, "omicsGMF_loadings")), c(n_features, n_comp))
    # design terms also carried through
    expect_identical(attr(out, "omicsGMF_design_X"), X)
    expect_identical(attr(out, "omicsGMF_design_Beta"), Beta)
})

# -------------------------
# impute_and_correct_with_omicsGMF: reuse + fallback semantics
# -------------------------

test_that("impute_and_correct_with_omicsGMF delegates to correct_with_omicsGMF with use_imputed = TRUE", {
    testthat::skip_if_not_installed("SingleCellExperiment")
    testthat::skip_if_not_installed("SummarizedExperiment")
    testthat::skip_if_not_installed("S4Vectors")

    m <- matrix(c(1, NA, 2, 5, 3, 6),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
    )
    sa <- data.frame(
        FullRunName = c("s1", "s2", "s3"),
        MS_batch = c("b1", "b2", "b2"),
        stringsAsFactors = FALSE
    )

    # Imputed matrix the fake fit returns.
    imputed <- matrix(c(1, 10, 2, 5, 3, 6),
        nrow = 2,
        dimnames = dimnames(m)
    )

    gmf_results <- matrix(c(1, 0, 0, 1, 0.5, 0.5),
        nrow = 3, ncol = 2,
        dimnames = list(NULL, c("c1", "c2"))
    )
    rotation <- matrix(c(0.6, 0.4, 0.3, 0.7),
        nrow = 2, ncol = 2,
        dimnames = list(c("f1", "f2"), c("c1", "c2"))
    )
    attr(gmf_results, "rotation") <- rotation
    X <- cbind("(Intercept)" = c(1, 1, 1), "MS_batchb2" = c(0, 1, 1))
    Beta <- matrix(c(1, 0.5, 0.2, 0.8),
        nrow = 2, ncol = 2,
        dimnames = list(c("f1", "f2"), c("(Intercept)", "MS_batchb2"))
    )
    attr(gmf_results, "X") <- X
    attr(gmf_results, "Beta") <- Beta

    fake_sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(dummy = matrix(0, nrow = nrow(m), ncol = ncol(m))),
        colData = S4Vectors::DataFrame(FullRunName = colnames(m))
    )
    SingleCellExperiment::reducedDim(fake_sce, "GMF") <- gmf_results

    fake_fit <- function(data_matrix, sample_annotation, design_formula,
                         family, ncomponents, gmf_args, impute_args) {
        list(
            sce = fake_sce, dimred_name = "GMF",
            imputed = imputed, imputed_assay = "omicsGMF_imputed"
        )
    }

    testthat::local_mocked_bindings(
        .pb_requireNamespace = function(...) invisible(TRUE),
        .omicsgmf_fit_and_impute = fake_fit,
        .package = "proBatch"
    )

    out <- impute_and_correct_with_omicsGMF(
        x = m, sample_annotation = sa, sample_id_col = "FullRunName",
        design_formula = ~MS_batch, batch_col = "MS_batch",
        ncomponents = 2L, format = "wide"
    )

    # Should equal imputed - batch_mean (operating on the imputed matrix, not raw).
    batch_idx <- which(colnames(X) == "MS_batchb2")
    batch_mean <- X[, batch_idx, drop = FALSE] %*% t(Beta[, batch_idx, drop = FALSE])
    expected <- imputed - t(batch_mean)
    out_numeric <- out
    for (a in grep("^omicsGMF_", names(attributes(out_numeric)), value = TRUE)) {
        attr(out_numeric, a) <- NULL
    }
    expect_equal(out_numeric, expected, tolerance = 1e-12)
    expect_false(any(is.na(out_numeric)))
})

test_that(".omicsgmf_reconstruct_corrected_matrix fallback_to_data preserves imputed values when no batch", {
    gmf_results <- matrix(c(1, 0, 0, 1, 0.5, 0.5),
        nrow = 3, ncol = 2,
        dimnames = list(NULL, c("c1", "c2"))
    )
    rotation <- matrix(c(0.6, 0.4, 0.3, 0.7),
        nrow = 2, ncol = 2,
        dimnames = list(c("f1", "f2"), c("c1", "c2"))
    )
    attr(gmf_results, "rotation") <- rotation

    imputed <- matrix(c(1, 10, 2, 5, 3, 6),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
    )

    # fallback_to_data = TRUE -> returns the input (imputed) matrix.
    out_imp <- proBatch:::`.omicsgmf_reconstruct_corrected_matrix`(
        gmf_results = gmf_results, data_matrix = imputed,
        batch_col = NULL, fallback_to_data = TRUE
    )
    out_imp_num <- out_imp
    for (a in grep("^omicsGMF_", names(attributes(out_imp_num)), value = TRUE)) {
        attr(out_imp_num, a) <- NULL
    }
    expect_equal(out_imp_num, imputed, tolerance = 1e-12)

    # Default (FALSE) keeps historical latent-only fallback (not equal to input).
    out_latent <- proBatch:::`.omicsgmf_reconstruct_corrected_matrix`(
        gmf_results = gmf_results, data_matrix = imputed, batch_col = NULL
    )
    out_latent_num <- out_latent
    for (a in grep("^omicsGMF_", names(attributes(out_latent_num)), value = TRUE)) {
        attr(out_latent_num, a) <- NULL
    }
    expected_latent <- t(gmf_results %*% t(rotation))
    # Inner reconstruction does not restore sample (column) names; the caller
    # .omicsgmf_correct_matrix_step does. Compare numerics only.
    expect_equal(unname(out_latent_num), unname(expected_latent), tolerance = 1e-12)
})
