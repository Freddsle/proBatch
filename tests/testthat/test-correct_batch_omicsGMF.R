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
                          max_rank,
                          run_rank,
                          rank_args,
                          gmf_args,
                          impute_args,
                          fill_the_missing) {
        captured$data_matrix <- data_matrix
        captured$sample_annotation <- sample_annotation
        captured$args <- list(
            sample_id_col = sample_id_col,
            design_formula = design_formula,
            family = family,
            ncomponents = ncomponents,
            max_rank = max_rank,
            run_rank = run_rank,
            rank_args = rank_args,
            gmf_args = gmf_args,
            impute_args = impute_args,
            fill_the_missing = fill_the_missing
        )
        storage.mode(data_matrix) <- "double"
        data_matrix + 0.5
    }

    local_fake_omicsgmf_step(fake_step)

    out <- impute_with_omicsGMF(
        x = m,
        sample_annotation = sa,
        sample_id_col = "FullRunName",
        design_formula = ~ Condition,
        family = poisson(),
        ncomponents = 3L,
        max_rank = 8L,
        run_rank = FALSE,
        rank_args = list(seed = 101),
        gmf_args = list(name = "custom_dimred"),
        impute_args = list(name = "custom_imputed"),
        fill_the_missing = FALSE,
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
    expect_identical(captured$args$max_rank, 8L)
    expect_false(captured$args$run_rank)
    expect_identical(captured$args$rank_args, list(seed = 101))
    expect_identical(captured$args$gmf_args, list(name = "custom_dimred"))
    expect_identical(captured$args$impute_args, list(name = "custom_imputed"))
    expect_false(captured$args$fill_the_missing)
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
                          max_rank,
                          run_rank,
                          rank_args,
                          gmf_args,
                          impute_args,
                          fill_the_missing) {
        captured$data_matrix <- data_matrix
        captured$sample_annotation <- sample_annotation
        captured$args <- list(
            sample_id_col = sample_id_col,
            design_formula = design_formula,
            family = family,
            ncomponents = ncomponents,
            max_rank = max_rank,
            run_rank = run_rank,
            rank_args = rank_args,
            gmf_args = gmf_args,
            impute_args = impute_args,
            fill_the_missing = fill_the_missing
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
        design_formula = ~ Condition,
        family = gaussian(),
        max_rank = 5L,
        rank_args = list(alpha = 0.1),
        gmf_args = list(name = "gmf_name"),
        impute_args = list(name = "imputed_name"),
        fill_the_missing = NULL,
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
    expect_identical(captured$args$rank_args, list(alpha = 0.1))
    expect_identical(captured$args$gmf_args, list(name = "gmf_name"))
    expect_identical(captured$args$impute_args, list(name = "imputed_name"))
})

