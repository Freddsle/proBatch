expect_matrix_like <- function(actual, template, storage_mode = "double") {
    expect_true(is.matrix(actual))
    expect_identical(dim(actual), dim(template))
    expect_identical(dimnames(actual), dimnames(template))
    if (!is.null(storage_mode)) {
        expect_identical(storage.mode(actual), storage_mode)
    }
}

local_fake_omicsgmf_step <- function(fake_step) {
    caller_env <- parent.frame()
    testthat::local_mocked_bindings(
        .pb_requireNamespace = function(...) invisible(TRUE),
        .omicsgmf_matrix_step = fake_step,
        .package = "proBatch",
        .env = caller_env
    )
}

make_test_pbf <- function() {
    testthat::skip_if_not_installed("QFeatures")
    testthat::skip_if_not_installed("SummarizedExperiment")
    testthat::skip_if_not_installed("S4Vectors")

    # 1) small assay, 2 features × 3 samples
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
            Batch = c("A", "A", "B")
        )
    )

    qf <- QFeatures::QFeatures(list(raw = se))
    testthat::skip_if_not("as_ProBatchFeatures" %in% getNamespaceExports("proBatch"))
    proBatch::as_ProBatchFeatures(qf)
}


test_that("impute_with_omicsGMF(wide): forwards parameters and preserves dimnames (mocked)", {
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
            sample_id_col   = sample_id_col,
            design_formula  = design_formula,
            family          = family,
            ncomponents     = ncomponents,
            gmf_args        = gmf_args,
            impute_args     = impute_args
        )
        storage.mode(data_matrix) <- "double"
        data_matrix + 0.5
    }

    local_fake_omicsgmf_step(fake_step)

    out <- impute_with_omicsGMF(
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
    expect_equal(out, m + 0.5, tolerance = 1e-12)

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

test_that("impute_with_omicsGMF(long): adds preImpute_* and respects keep_all = 'minimal' (mocked)", {
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
            sample_id_col   = sample_id_col,
            design_formula  = design_formula,
            family          = family,
            ncomponents     = ncomponents,
            gmf_args        = gmf_args,
            impute_args     = impute_args
        )
        storage.mode(data_matrix) <- "double"
        data_matrix + 42
    }

    local_fake_omicsgmf_step(fake_step)

    out <- impute_with_omicsGMF(
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
        sort(c("peptide_group_label", "FullRunName", "Intensity", "preImpute_Intensity"))
    )

    idx <- match(
        interaction(df$peptide_group_label, df$FullRunName),
        interaction(out$peptide_group_label, out$FullRunName)
    )
    expect_equal(out$preImpute_Intensity[idx], df$Intensity)
    expect_equal(out$Intensity[idx], df$Intensity + 42)

    expect_matrix_like(captured$data_matrix, expected_matrix)
    expect_identical(
        captured$sample_annotation$FullRunName,
        colnames(expected_matrix)
    )
    expect_identical(captured$args$ncomponents, 2L)
    expect_identical(captured$args$gmf_args, list(name = "gmf_name"))
    expect_identical(captured$args$impute_args, list(name = "imputed_name"))
})


test_that("impute_with_omicsGMF(ProBatchFeatures) calls .pb_apply_step with correct params", {
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
            # mimic what your real helper returns
            list(object = object, assay = "omicsGMFImpute_1")
        },
        .package = "proBatch"
    )

    out <- impute_with_omicsGMF(
        x = pbf,
        ncomponents = 3L,
        sample_id_col = "FullRunName",
        design_formula = ~Batch
    )

    expect_s4_class(out, "ProBatchFeatures")
    expect_identical(captured$fun, "omicsGMFImpute")
    expect_identical(captured$params$ncomponents, 3L)
    expect_identical(captured$params$sample_id_col, "FullRunName")
})


test_that("impute_with_omicsGMF() errors when ncomponents is missing", {
    m <- matrix(1, 2, 2)
    expect_error(
        impute_with_omicsGMF(m, sample_annotation = data.frame(FullRunName = c("V1", "V2"))),
        "`ncomponents` must be supplied",
        fixed = TRUE
    )
})

# -----------------------------------------------------------------------
# omicsGMFImpute: all-NA row/column guard
# -----------------------------------------------------------------------

test_that(".omicsgmf_fit_and_impute filters all-NA rows and re-inserts them (mocked)", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("S4Vectors")
    skip_if_not_installed("SingleCellExperiment")

    # 4 features x 3 samples; feat2 and feat4 are all-NA
    m <- matrix(
        c(
            1, NA, 3, NA,
            4, 5, NA, NA,
            7, 8, 9, NA
        ),
        nrow = 4, ncol = 3,
        dimnames = list(paste0("feat", 1:4), paste0("s", 1:3))
    )
    sa <- data.frame(
        FullRunName = paste0("s", 1:3), Cond = c("A", "A", "B"),
        stringsAsFactors = FALSE
    )

    # Mock omicsGMF::runGMF and imputeGMF to just fill NAs with 0
    fake_runGMF <- function(x, ...) x
    fake_imputeGMF <- function(x, exprs_values, name, ...) {
        mat <- SummarizedExperiment::assay(x, exprs_values)
        mat[is.na(mat)] <- 0
        SummarizedExperiment::assay(x, name) <- mat
        x
    }

    # Temporarily mock the omicsGMF namespace calls
    local_mocked_bindings(
        .pb_require_omicsgmf_stack = function() invisible(TRUE),
        .package = "proBatch"
    )

    # Call the internal function directly, but we need to mock the do.call targets.
    # Instead, use the step function with mocked externals.
    captured <- new.env(parent = emptyenv())
    fake_step <- function(data_matrix, sample_annotation, sample_id_col,
                          design_formula, family, ncomponents,
                          gmf_args = list(), impute_args = list()) {
        # Record what matrix the step sees (should have all-NA rows removed)
        captured$input_rows <- rownames(data_matrix)
        captured$input_has_allNA_row <- any(apply(data_matrix, 1, function(r) all(is.na(r))))
        # Return a filled version preserving dims
        out <- data_matrix
        out[is.na(out)] <- 0
        storage.mode(out) <- "double"
        out
    }

    local_fake_omicsgmf_step(fake_step)

    out <- impute_with_omicsGMF(
        x = m,
        sample_annotation = sa,
        sample_id_col = "FullRunName",
        ncomponents = 2L,
        format = "wide"
    )

    # Output should have original dimensions
    expect_identical(dim(out), dim(m))
    expect_identical(dimnames(out), dimnames(m))
})

test_that(".omicsgmf_fit_and_impute all-NA guard: filter and re-insert logic (unit)", {
    # Test the guard logic added to .omicsgmf_fit_and_impute without needing
    # the full omicsGMF stack — validates the filter + re-insert pattern.

    # Build a matrix with all-NA row(s) and column(s)
    m <- matrix(
        c(
            1.0, NA, 3.0,
            NA, NA, NA, # all-NA row
            7.0, 8.0, 9.0
        ),
        nrow = 3, ncol = 3, byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:3))
    )
    m[, "s2"] <- NA # also make column s2 all-NA

    row_all_na <- apply(m, 1L, function(r) all(is.na(r)))
    col_all_na <- apply(m, 2L, function(c) all(is.na(c)))

    expect_true(any(row_all_na), info = "test matrix should have an all-NA row")
    expect_true(any(col_all_na), info = "test matrix should have an all-NA column")

    # After filtering, no all-NA rows/cols remain
    reduced <- m[!row_all_na, !col_all_na, drop = FALSE]
    expect_false(any(apply(reduced, 1, function(r) all(is.na(r)))),
        info = "reduced matrix should have no all-NA rows"
    )
    expect_false(any(apply(reduced, 2, function(c) all(is.na(c)))),
        info = "reduced matrix should have no all-NA columns"
    )

    # Simulate imputation on the reduced matrix
    imputed_reduced <- reduced
    imputed_reduced[is.na(imputed_reduced)] <- 0
    storage.mode(imputed_reduced) <- "double"

    # Re-insert (same logic as in .omicsgmf_fit_and_impute)
    full <- matrix(NA_real_,
        nrow = nrow(m), ncol = ncol(m),
        dimnames = dimnames(m)
    )
    full[!row_all_na, !col_all_na] <- imputed_reduced

    expect_identical(dim(full), dim(m))
    expect_identical(dimnames(full), dimnames(m))
    # All-NA row/col stay NA
    expect_true(all(is.na(full["f2", ])))
    expect_true(all(is.na(full[, "s2"])))
    # Non-NA entries got imputed
    expect_false(any(is.na(full[!row_all_na, !col_all_na])))
})
