toy_matrix_nai <- matrix(
    # fmt: skip
    c(
        10, NA, 12,
        7, 8, 9,
        NA, 5, 6
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("prot", 1:3), paste0("S", 1:3))
)

toy_matrix_nai2 <- matrix(
    # fmt: skip
    c(
        NA, 2, 3,
        4, NA, 6,
        7, 8, 9
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("prot", 1:3), paste0("S", 1:3))
)

toy_sa_nai <- data.frame(
    FullRunName = colnames(toy_matrix_nai),
    Condition = c("A", "B", "B"),
    stringsAsFactors = FALSE
)
rownames(toy_sa_nai) <- toy_sa_nai$FullRunName

pbf_toy_nai <- ProBatchFeatures(
    data_matrix = toy_matrix_nai,
    sample_annotation = toy_sa_nai,
    sample_id_col = "FullRunName",
    name = "toy"
)

toy_assay_nai <- pb_current_assay(pbf_toy_nai)
toy_assay_nai_alt <- paste0(toy_assay_nai, "_alt")

se_alt_nai <- SummarizedExperiment(
    assays = list(intensity = toy_matrix_nai2),
    colData = colData(pbf_toy_nai[[toy_assay_nai]])
)

pbf_multi_nai <- proBatch:::.pb_add_assay_with_link(
    pbf_toy_nai,
    se = se_alt_nai,
    to = toy_assay_nai_alt,
    from = toy_assay_nai
)


# ------------------------------------------------------------------
# Internal helper: .pb_NA_intensity_stats_ungrouped
# ------------------------------------------------------------------

test_that(".pb_NA_intensity_stats_ungrouped computes mean and prop_missing", {
    res <- proBatch:::.pb_NA_intensity_stats_ungrouped(toy_matrix_nai)

    expect_s3_class(res, "data.frame")
    expect_named(res, c("mean_intensity", "prop_missing", "n_samples"))
    expect_equal(nrow(res), 3L)

    idx1 <- which(rownames(res) == "prot1")
    expect_equal(res$mean_intensity[idx1], 11)
    expect_equal(res$prop_missing[idx1], 1 / 3)

    idx2 <- which(rownames(res) == "prot2")
    expect_equal(res$prop_missing[idx2], 0)

    idx3 <- which(rownames(res) == "prot3")
    expect_equal(res$mean_intensity[idx3], 5.5)
    expect_equal(res$prop_missing[idx3], 1 / 3)
})

test_that(".pb_NA_intensity_stats_ungrouped returns NULL for all-NA matrix", {
    all_na <- matrix(
        NA_real_,
        nrow = 2,
        ncol = 3,
        dimnames = list(c("a", "b"), c("s1", "s2", "s3"))
    )
    expect_null(proBatch:::.pb_NA_intensity_stats_ungrouped(all_na))
})

test_that(".pb_NA_intensity_stats_ungrouped returns NULL for empty matrix", {
    empty <- matrix(numeric(0), nrow = 0, ncol = 3)
    expect_null(proBatch:::.pb_NA_intensity_stats_ungrouped(empty))
})


# ------------------------------------------------------------------
# Internal helper: .pb_NA_intensity_stats
# ------------------------------------------------------------------

test_that(".pb_NA_intensity_stats returns ungrouped result without color_by", {
    res <- proBatch:::.pb_NA_intensity_stats(
        data_matrix = toy_matrix_nai,
        sample_annotation = toy_sa_nai,
        sample_id_col = "FullRunName",
        color_by = NULL
    )

    expect_s3_class(res, "data.frame")
    expect_false(".group" %in% names(res))
    expect_equal(nrow(res), 3L)
})

test_that(".pb_NA_intensity_stats adds groups when color_by is supplied", {
    res <- proBatch:::.pb_NA_intensity_stats(
        data_matrix = toy_matrix_nai,
        sample_annotation = toy_sa_nai,
        sample_id_col = "FullRunName",
        color_by = "Condition"
    )

    expect_s3_class(res, "data.frame")
    expect_true(".group" %in% names(res))
    expect_setequal(unique(res$.group), c("A", "B"))

    group_a <- res[res$.group == "A", , drop = FALSE]
    expect_equal(group_a$n_samples, rep(1L, nrow(group_a)))

    group_b <- res[res$.group == "B", , drop = FALSE]
    expect_equal(group_b$n_samples, rep(2L, nrow(group_b)))
})


# ------------------------------------------------------------------
# Internal helper: .pb_NA_intensity_cor_labels
# ------------------------------------------------------------------

test_that(".pb_NA_intensity_cor_labels computes per-group rho", {
    stats_df <- data.frame(
        mean_intensity = c(1, 2, 3, 4, 5, 6),
        prop_missing = c(0.5, 0.4, 0.3, 0.1, 0.3, 0.5),
        n_samples = 3L,
        .group = rep(c("X", "Y"), each = 3),
        stringsAsFactors = FALSE
    )
    cor_df <- proBatch:::.pb_NA_intensity_cor_labels(
        stats_df,
        has_group = TRUE
    )

    expect_s3_class(cor_df, "data.frame")
    expect_equal(nrow(cor_df), 2L)
    expect_setequal(cor_df$.group, c("X", "Y"))
    expect_true(all(grepl("^rho == ", cor_df$label)))
    expect_false(any(grepl("[^ -~]", cor_df$label)))
    expect_true(all(cor_df$x == Inf))
})

test_that(".pb_NA_intensity_cor_labels returns one row without groups", {
    stats_df <- data.frame(
        mean_intensity = c(1, 2, 3),
        prop_missing = c(0.6, 0.3, 0.1),
        n_samples = 5L,
        stringsAsFactors = FALSE
    )
    cor_df <- proBatch:::.pb_NA_intensity_cor_labels(
        stats_df,
        has_group = FALSE
    )

    expect_s3_class(cor_df, "data.frame")
    expect_equal(nrow(cor_df), 1L)
    expect_true(grepl("^rho == ", cor_df$label))
    expect_false(grepl("[^ -~]", cor_df$label))
})

test_that(".pb_NA_intensity_cor_labels needs three ungrouped rows", {
    stats_df <- data.frame(
        mean_intensity = c(1, 2),
        prop_missing = c(0.5, 0.3),
        n_samples = 3L,
        stringsAsFactors = FALSE
    )
    expect_null(proBatch:::.pb_NA_intensity_cor_labels(
        stats_df,
        has_group = FALSE
    ))
})

test_that(".pb_NA_intensity_cor_labels computes labels per assay", {
    stats_df <- data.frame(
        mean_intensity = rep(c(1, 2, 3), 2),
        prop_missing = c(0.9, 0.5, 0.1, 0.1, 0.5, 0.9),
        n_samples = 3L,
        pbf_name = rep(c("first", "second"), each = 3),
        stringsAsFactors = FALSE
    )

    cor_df <- proBatch:::.pb_NA_intensity_cor_labels(
        stats_df,
        has_group = FALSE
    )

    expect_identical(as.character(cor_df$pbf_name), c("first", "second"))
    expect_identical(cor_df$label, c("rho == -1", "rho == 1"))
})


# ------------------------------------------------------------------
# Public default method
# ------------------------------------------------------------------

test_that("plot_NA_intensity returns a ggplot for a matrix", {
    p <- plot_NA_intensity(toy_matrix_nai)

    expect_s3_class(p, "ggplot")
    expect_equal(p$labels$x, "Mean log2 intensity (non-missing values)")
    expect_equal(p$labels$y, "% missing values")
    expect_true(nrow(p$data) > 0)
})

test_that("plot_NA_intensity renders portable plotmath correlation labels", {
    p <- plot_NA_intensity(toy_matrix_nai, spline_df = 0)
    text_layers <- Filter(
        function(layer) inherits(layer$geom, "GeomText"),
        p$layers
    )

    expect_length(text_layers, 1L)
    expect_true(isTRUE(text_layers[[1]]$geom_params$parse))
    expect_false(any(grepl("[^ -~]", text_layers[[1]]$data$label)))
    expect_no_error(ggplot2::ggplotGrob(p))
})

test_that("plot_NA_intensity returns an empty plot for an all-NA matrix", {
    all_na <- matrix(
        NA_real_,
        nrow = 2,
        ncol = 3,
        dimnames = list(c("a", "b"), c("s1", "s2", "s3"))
    )
    expect_warning(
        p <- plot_NA_intensity(all_na),
        "No finite mean intensities"
    )
    expect_s3_class(p, "ggplot")
})

test_that("plot_NA_intensity warns on zero-dimensional input", {
    empty <- matrix(numeric(0), nrow = 0, ncol = 0)
    expect_warning(
        p <- plot_NA_intensity(empty),
        "zero rows or columns"
    )
    expect_s3_class(p, "ggplot")
})

test_that("plot_NA_intensity colours by group", {
    p <- plot_NA_intensity(
        toy_matrix_nai,
        sample_annotation = toy_sa_nai,
        sample_id_col = "FullRunName",
        color_by = "Condition"
    )

    expect_s3_class(p, "ggplot")
    expect_true(".group" %in% names(p$data))
    expect_setequal(unique(p$data$.group), c("A", "B"))
})

test_that("plot_NA_intensity respects a custom color scheme", {
    p <- plot_NA_intensity(
        toy_matrix_nai,
        sample_annotation = toy_sa_nai,
        sample_id_col = "FullRunName",
        color_by = "Condition",
        color_scheme = c(A = "#111111", B = "#222222"),
        spline_df = 0
    )

    expect_no_warning(built <- ggplot_build(p))
    expect_setequal(unique(built$data[[1]]$colour), c("#111111", "#222222"))
})

test_that("plot_NA_intensity suppresses smoothing with spline_df zero", {
    p <- plot_NA_intensity(toy_matrix_nai, spline_df = 0)

    layer_classes <- vapply(
        p$layers,
        function(layer) class(layer$geom)[1],
        character(1)
    )
    expect_false("GeomSmooth" %in% layer_classes)
})

test_that("plot_NA_intensity includes smoothing with positive spline_df", {
    p <- plot_NA_intensity(toy_matrix_nai, spline_df = 3)

    layer_classes <- vapply(
        p$layers,
        function(layer) class(layer$geom)[1],
        character(1)
    )
    expect_true("GeomSmooth" %in% layer_classes)
})

test_that("plot_NA_intensity accepts a SummarizedExperiment", {
    se <- SummarizedExperiment(
        assays = list(intensity = toy_matrix_nai),
        colData = DataFrame(toy_sa_nai)
    )
    p <- plot_NA_intensity(se, color_by = "Condition")

    expect_s3_class(p, "ggplot")
    expect_true(".group" %in% names(p$data))
})

test_that("plot_NA_intensity passes point aesthetics through", {
    p <- plot_NA_intensity(
        toy_matrix_nai,
        point_alpha = 0.5,
        point_size = 2
    )

    point_layer <- p$layers[[1]]
    expect_equal(point_layer$aes_params$alpha, 0.5)
    expect_equal(point_layer$aes_params$size, 2)
})


# ------------------------------------------------------------------
# ProBatchFeatures method
# ------------------------------------------------------------------

test_that("plot_NA_intensity supports one ProBatchFeatures assay", {
    p <- plot_NA_intensity(pbf_toy_nai)

    expect_s3_class(p, "ggplot")
    expect_true(nrow(p$data) > 0)
})

test_that("plot_NA_intensity supports ProBatchFeatures grouping", {
    p <- plot_NA_intensity(pbf_toy_nai, color_by = "Condition")

    expect_s3_class(p, "ggplot")
    expect_true(".group" %in% names(p$data))
    expect_setequal(unique(p$data$.group), c("A", "B"))
})

test_that("plot_NA_intensity facets multiple ProBatchFeatures assays", {
    p <- plot_NA_intensity(
        pbf_multi_nai,
        pbf_name = c(toy_assay_nai, toy_assay_nai_alt),
        nrow = 1,
        ncol = 2
    )

    expect_s3_class(p, "ggplot")
    expect_setequal(
        unique(as.character(p$data$pbf_name)),
        c(toy_assay_nai, toy_assay_nai_alt)
    )

    text_layers <- Filter(
        function(layer) inherits(layer$geom, "GeomText"),
        p$layers
    )
    expect_length(text_layers, 1L)
    expect_setequal(
        as.character(text_layers[[1]]$data$pbf_name),
        c(toy_assay_nai, toy_assay_nai_alt)
    )
})

test_that("plot_NA_intensity facets grouped ProBatchFeatures assays", {
    p <- plot_NA_intensity(
        pbf_multi_nai,
        pbf_name = c(toy_assay_nai, toy_assay_nai_alt),
        color_by = "Condition",
        nrow = 1,
        ncol = 2
    )

    expect_s3_class(p, "ggplot")
    expect_setequal(
        unique(as.character(p$data$pbf_name)),
        c(toy_assay_nai, toy_assay_nai_alt)
    )
    expect_setequal(unique(p$data$.group), c("A", "B"))
})
