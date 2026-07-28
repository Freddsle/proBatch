testthat::skip_if_not_installed("QFeatures")

make_test_pbf <- function(mat, sa_extra = NULL) {
    stopifnot(!is.null(colnames(mat)))
    sa <- data.frame(Sample = colnames(mat), stringsAsFactors = FALSE)
    if (!is.null(sa_extra)) {
        stopifnot(is.data.frame(sa_extra))
        if (nrow(sa_extra) != nrow(sa)) {
            stop("`sa_extra` must have one row per sample.")
        }
        sa <- cbind(sa, sa_extra)
    }
    suppressMessages(ProBatchFeatures(
        data_matrix = mat,
        sample_annotation = sa,
        sample_id_col = "Sample",
        name = "raw",
        level = "feature"
    ))
}

test_that("pb_zeroIsNA converts zeros to missing values and logs the step", {
    mat <- matrix(
        c(0, 1, 2, 0, 3, 4),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)
    assay_name <- pb_current_assay(pbf)
    zero_count <- sum(assay(pbf[[assay_name]], "intensity") == 0)
    expect_gt(zero_count, 0L)

    res <- suppressMessages(pb_zeroIsNA(pbf))
    expect_s4_class(res, "ProBatchFeatures")

    updated <- assay(res[[assay_name]], "intensity")
    expect_equal(sum(is.na(updated)), zero_count)

    log <- get_operation_log(res)
    expect_equal(nrow(log), 1L)
    expect_identical(as.character(log$step), "zeroIsNA")
    expect_identical(as.character(log$from), assay_name)
    expect_identical(as.character(log$to), assay_name)
    expect_identical(log$params[[1]], list())
})

test_that("pb_infIsNA replaces infinities and records the operation", {
    mat <- matrix(
        c(Inf, 2, 3, 4, 5, Inf),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)
    assay_name <- pb_current_assay(pbf)
    inf_count <- sum(is.infinite(assay(pbf[[assay_name]], "intensity")))
    expect_gt(inf_count, 0L)

    res <- suppressMessages(pb_infIsNA(pbf))
    expect_s4_class(res, "ProBatchFeatures")

    updated <- assay(res[[assay_name]], "intensity")
    expect_equal(sum(is.na(updated)), inf_count)

    log <- get_operation_log(res)
    expect_equal(nrow(log), 1L)
    expect_identical(as.character(log$step), "infIsNA")
    expect_identical(as.character(log$from), assay_name)
    expect_identical(as.character(log$to), assay_name)
    expect_identical(log$params[[1]], list())
})

test_that("pb_nNA returns per-assay results for multiple assays", {
    mat <- matrix(
        c(1, NA, 3, 4, 5, NA),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)
    assay_name <- pb_current_assay(pbf)

    # Create an additional stored assay with different missing pattern
    se_dup <- pbf[[assay_name]]
    alt_mat <- assay(se_dup, "intensity")
    alt_mat[1, 1] <- NA
    assay(se_dup, "intensity") <- alt_mat
    new_name <- paste0(assay_name, "_alt")
    pbf <- addAssay(pbf, se_dup, name = new_name)

    res <- pb_nNA(pbf, c(assay_name, new_name))
    expect_type(res, "list")
    expect_identical(names(res), c(assay_name, new_name, "nNA"))

    qf1 <- QFeatures(setNames(list(pbf[[assay_name]]), assay_name))
    qf2 <- QFeatures(setNames(list(pbf[[new_name]]), new_name))

    expect_type(res[[1]], "list")
    expect_type(res[[2]], "list")

    expect_identical(
        res[[1]],
        nNA(qf1, i = assay_name)
    )
    expect_identical(
        res[[2]],
        nNA(qf2, i = new_name)
    )
})

test_that("pb_filterNA stores filtered assays when not operating in place", {
    mat <- matrix(
        c(1, NA, 3, 4, 5, 6),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)
    original_names <- names(pbf)
    assay_name <- pb_current_assay(pbf)

    res <- suppressMessages(pb_filterNA(pbf, inplace = FALSE))
    expect_s4_class(res, "ProBatchFeatures")

    expect_setequal(original_names, names(pbf))
    expect_equal(length(names(res)), length(original_names) + 1L)

    new_name <- setdiff(names(res), original_names)
    expect_length(new_name, 1L)
    expect_match(new_name, paste0(assay_name, "_filteredNA"))

    filtered <- assay(res[[new_name]], "intensity")
    expect_true(is.matrix(filtered))
    expect_false(anyNA(filtered))

    log <- get_operation_log(res)
    expect_equal(nrow(log), 1L)
    expect_identical(as.character(log$step), "filterNA")
    expect_identical(as.character(log$from), assay_name)
    expect_identical(as.character(log$to), new_name)
    expect_true(isFALSE(log$params[[1]]$inplace) || is.null(log$params[[1]]$inplace))
})

test_that("pb_filterNA modifies stored assays in place when requested", {
    mat <- matrix(
        c(1, NA, 3, 4, 5, 6),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)
    assay_name <- pb_current_assay(pbf)

    res <- suppressMessages(pb_filterNA(pbf, inplace = TRUE))
    expect_s4_class(res, "ProBatchFeatures")

    expect_identical(names(res), names(pbf))

    filtered <- assay(res[[assay_name]], "intensity")
    expect_true(is.matrix(filtered))
    expect_false(anyNA(filtered))

    log <- get_operation_log(res)
    expect_equal(nrow(log), 1L)
    expect_identical(as.character(log$step), "filterNA")
    expect_identical(as.character(log$from), assay_name)
    expect_identical(as.character(log$to), assay_name)
    expect_identical(log$params[[1]]$inplace, TRUE)
})

test_that("pb_filterNA validates final_name length when creating new assays", {
    mat <- matrix(
        c(1, NA, 3, 4),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:2), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)

    expect_error(
        pb_filterNA(pbf, inplace = FALSE, final_name = c("a", "b")),
        "`final_name` must be length 1 or match `pbf_name`"
    )
})

test_that("pb_groupfilterNA retains union of group-wise valid features in place", {
    mat <- matrix(
        c(
            1, 2, 3, 4,
            NA, NA, 5, 6,
            1, NA, 2, NA
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:4))
    )
    sa_extra <- data.frame(
        Batch = rep(c("B1", "B2"), each = 2),
        stringsAsFactors = FALSE
    )
    pbf <- make_test_pbf(mat, sa_extra = sa_extra)
    assay_name <- pb_current_assay(pbf)

    res <- suppressMessages(pb_groupfilterNA(
        pbf,
        group_cols = "Batch",
        min_valid = 2L,
        inplace = TRUE
    ))

    expect_s4_class(res, "ProBatchFeatures")
    expect_identical(names(res), names(pbf))

    filtered <- assay(res[[assay_name]], "intensity")
    expect_true(is.matrix(filtered))
    expect_identical(rownames(filtered), c("f1", "f2"))

    log <- get_operation_log(res)
    expect_equal(nrow(log), 1L)
    expect_identical(as.character(log$step), "groupfilterNA")
    expect_identical(as.character(log$from), assay_name)
    expect_identical(as.character(log$to), assay_name)
    expect_identical(log$params[[1]]$group_cols, "Batch")
    expect_identical(log$params[[1]]$min_valid, 2L)
    expect_null(log$params[[1]]$pNA)
})

test_that("pb_groupfilterNA stores union of group-wise valid features when not operating in place", {
    mat <- matrix(
        c(
            1, 2, 3, 4,
            NA, NA, 5, 6,
            1, NA, 2, NA
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:4))
    )
    sa_extra <- data.frame(
        Batch = rep(c("B1", "B2"), each = 2),
        stringsAsFactors = FALSE
    )
    pbf <- make_test_pbf(mat, sa_extra = sa_extra)
    assay_name <- pb_current_assay(pbf)
    original_names <- names(pbf)

    res <- suppressMessages(pb_groupfilterNA(
        pbf,
        group_cols = "Batch",
        min_valid = 2L,
        inplace = FALSE,
        final_name = "filtered_group"
    ))

    expect_s4_class(res, "ProBatchFeatures")
    expect_setequal(original_names, names(pbf))
    expect_equal(length(names(res)), length(original_names) + 1L)

    new_name <- setdiff(names(res), original_names)
    expect_identical(new_name, "filtered_group")

    filtered <- assay(res[[new_name]], "intensity")
    expect_identical(rownames(filtered), c("f1", "f2"))

    log <- get_operation_log(res)
    expect_equal(nrow(log), 1L)
    expect_identical(as.character(log$step), "groupfilterNA")
    expect_identical(as.character(log$from), assay_name)
    expect_identical(as.character(log$to), new_name)
    expect_identical(log$params[[1]]$group_cols, "Batch")
    expect_identical(log$params[[1]]$min_valid, 2L)
    expect_null(log$params[[1]]$pNA)
})

test_that("pb_groupfilterNA respects pNA thresholds when supplied", {
    mat <- matrix(
        c(
            1, 2, 3, 4,
            NA, 2, NA, 4,
            NA, NA, 5, 6
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:4))
    )
    sa_extra <- data.frame(
        Batch = rep(c("B1", "B2"), each = 2),
        stringsAsFactors = FALSE
    )
    pbf <- make_test_pbf(mat, sa_extra = sa_extra)
    assay_name <- pb_current_assay(pbf)

    res <- suppressMessages(pb_groupfilterNA(
        pbf,
        group_cols = "Batch",
        min_valid = NULL,
        pNA = 0.25,
        inplace = TRUE
    ))

    filtered <- assay(res[[assay_name]], "intensity")
    expect_identical(rownames(filtered), c("f1", "f3"))

    log <- get_operation_log(res)
    expect_equal(nrow(log), 1L)
    expect_null(log$params[[1]]$min_valid)
    expect_equal(log$params[[1]]$pNA, 0.25)
})

test_that("pb_groupfilterNA combines min_valid and pNA per group", {
    mat <- matrix(
        c(
            1, 2, 3, 4, 5, 6,
            NA, 2, 3, 7, 8, 9,
            NA, NA, 5, 6, NA, NA
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:6))
    )
    sa_extra <- data.frame(
        Batch = c(rep("B1", 3), rep("B2", 3)),
        stringsAsFactors = FALSE
    )
    pbf <- make_test_pbf(mat, sa_extra = sa_extra)
    assay_name <- pb_current_assay(pbf)

    res <- suppressMessages(pb_groupfilterNA(
        pbf,
        group_cols = "Batch",
        min_valid = 2L,
        pNA = 0.25,
        inplace = TRUE
    ))

    filtered <- assay(res[[assay_name]], "intensity")
    expect_identical(rownames(filtered), c("f1", "f2"))

    log <- get_operation_log(res)
    expect_equal(nrow(log), 1L)
    expect_equal(log$params[[1]]$min_valid, 2L)
    expect_equal(log$params[[1]]$pNA, 0.25)
})

test_that("pb_groupfilterNA validates presence of grouping columns", {
    mat <- matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:2), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)

    expect_error(
        pb_groupfilterNA(pbf, group_cols = "Batch"),
        "missing group column(s)",
        fixed = TRUE
    )
})

test_that("pb_groupfilterNA errors when a group has too few samples", {
    mat <- matrix(
        c(
            1, 2, 3,
            4, 5, NA
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:2), paste0("s", 1:3))
    )
    sa_extra <- data.frame(
        Batch = c("B1", "B1", "B2"),
        stringsAsFactors = FALSE
    )
    pbf <- make_test_pbf(mat, sa_extra = sa_extra)

    expect_error(
        pb_groupfilterNA(pbf, group_cols = "Batch", min_valid = 2L),
        "requires at least 2",
        fixed = TRUE
    )
})

test_that("pb_missing helpers error on non-materialised assays with guidance", {
    mat <- matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:2), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)
    assay_name <- pb_current_assay(pbf)

    expect_error(
        pb_zeroIsNA(pbf, pbf_name = "feature::not_there"),
        "are not stored in the object",
        fixed = TRUE
    )

    logged_only <- proBatch:::.pb_add_log_entry(
        pbf,
        step = "log2",
        fun = "log_transform_dm",
        from = assay_name,
        to = paste0(assay_name, "_log2"),
        params = list()
    )

    expect_error(
        pb_zeroIsNA(logged_only, pbf_name = paste0(assay_name, "_log2")),
        "store_fast_steps = TRUE",
        fixed = TRUE
    )
})
