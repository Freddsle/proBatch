test_that("hierarchical_clustering", {
    pb_test_load_example_data()

    matrix_test <- example_proteome_matrix[1:10, ]

    # exclude columns EarTag Strain Sex RunDate RunTime digestion_batch
    example_sample_annotation <- example_sample_annotation %>%
        select("FullRunName", "MS_batch", "Diet")

    color_list <- sample_annotation_to_colors(
        example_sample_annotation,
        sample_id_col = "FullRunName",
        factor_columns = c("MS_batch", "Diet")
    )

    pb_test_expect_warnings(
        hiearchical <- plot_hierarchical_clustering(matrix_test,
            sample_annotation = example_sample_annotation,
            factors_to_plot = c("MS_batch", "Diet"),
            distance = "euclidean", agglomeration = "complete",
            color_list = color_list,
            label_samples = FALSE
        ),
        c(
            "filling missing values with 0",
            "Hierarchical clustering cannot operate with missing values",
            "color list and sample annotation have different factors,"
        ),
        fixed = TRUE
    )

    expect_identical(hiearchical$mar[[1]], 1)
    expect_identical(hiearchical$mar[[2]], 5)
    expect_identical(hiearchical$mar[[3]], 0)
    expect_identical(hiearchical$mar[[4]], 1)
})


test_that("heatmap_plot", {
    pb_test_load_example_data()

    matrix_test <- example_proteome_matrix[1:20, ]
    example_sample_annotation <- example_sample_annotation %>%
        select("FullRunName", "MS_batch", "Sex", "digestion_batch", "Diet")


    color_list <- sample_annotation_to_colors(
        example_sample_annotation,
        sample_id_col = "FullRunName",
        factor_columns = c("MS_batch", "Sex", "digestion_batch", "Diet")
    )

    suppressWarnings(
        heatmap <- plot_heatmap_diagnostic(
            matrix_test,
            sample_annotation = example_sample_annotation,
            factors_to_plot = c("MS_batch", "Sex", "digestion_batch", "Diet"),
            cluster_cols = TRUE,
            show_rownames = TRUE, show_colnames = FALSE,
            color_list = color_list
        )
    )

    expect_equal(heatmap$tree_row$method, "complete")
    expect_equal(heatmap$tree_row$dist.method, "euclidean")

    expect_equal(heatmap$tree_row$labels[1], "10062_NVGVSFYADKPEVTQEQK_3")
    expect_equal(heatmap$tree_row$labels[2], "10063_NVGVSFYADKPEVTQEQKK_3")

    expect_equal(heatmap$gtable$layout$name[1], "col_tree")
    expect_equal(heatmap$gtable$layout$name[3], "matrix")
    expect_equal(heatmap$gtable$layout$name[5], "col_annotation")
    expect_equal(heatmap$gtable$layout$name[8], "legend")

    expect_equal(heatmap$gtable$layout$t, c(2, 4, 4, 4, 3, 3, 3, 3))
})


test_that("plot_PCA warns on missing values and honors x/y PC selection", {
    pb_test_load_example_data()

    pb_test_expect_warnings(
        pca <- plot_PCA(
            example_proteome_matrix, example_sample_annotation,
            color_by = "MS_batch", plot_title = "PCA colored by MS batch",
            fill_the_missing = -1
        ),
        c(
            "filling missing values with -1",
            "PCA cannot operate with missing values in the matrix",
            "color_scheme will be inferred automatically"
        ),
        fixed = TRUE
    )

    expect_equal(pca$labels$y, "PC2 (14.24%)")
    expect_equal(pca$labels$x, "PC1 (69.50%)")
    expect_equal(pca$labels$colour, "MS_batch")

    pca_xy <- suppressWarnings(plot_PCA(
        example_proteome_matrix, example_sample_annotation,
        color_by = "MS_batch",
        fill_the_missing = -1,
        x_nPC = 2,
        y_nPC = 3
    ))

    expect_match(pca_xy$labels$x, "^PC2 ")
    expect_match(pca_xy$labels$y, "^PC3 ")
})

test_that("plot_PCA supports marginal density plots", {
    pb_test_load_example_data()

    pca <- suppressMessages(suppressWarnings(plot_PCA(
        example_proteome_matrix, example_sample_annotation,
        color_by = "MS_batch",
        fill_the_missing = -1,
        marginal_density = TRUE
    )))

    expect_true(inherits(pca, "ggplot") || grid::is.grob(pca))
})

test_that("plot_PCA ProBatchFeatures arranges assays and preserves subset order", {
    skip_if_not_installed("gridExtra")
    pbf <- pb_test_make_pbf(n_rows = 40, n_cols = 6, add_log2 = TRUE)

    res <- suppressWarnings(plot_PCA(
        pbf,
        sample_id_col = "FullRunName",
        return_gridExtra = TRUE
    ))

    expect_type(res, "list")
    expect_named(res$plots, names(pbf))
    expect_equal(length(res$plots), length(names(pbf)))
    expect_true(all(vapply(res$plots, inherits, logical(1), "ggplot")))

    subset_assays <- rev(names(pbf))
    res_subset <- suppressWarnings(plot_PCA(
        pbf,
        pbf_name = subset_assays,
        sample_id_col = "FullRunName",
        return_gridExtra = TRUE
    ))

    expect_type(res_subset, "list")
    expect_equal(names(res_subset$plots), subset_assays)
    expect_true(all(vapply(res_subset$plots, inherits, logical(1), "ggplot")))
})

test_that("plot_heatmap_diagnostic ProBatchFeatures arranges assays and preserves subset order", {
    skip_if_not_installed("gridExtra")
    pbf <- pb_test_make_pbf(n_rows = 30, n_cols = 5, add_log2 = TRUE)

    res <- suppressWarnings(plot_heatmap_diagnostic(
        pbf,
        sample_id_col = "FullRunName",
        factors_to_plot = c("MS_batch"),
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        return_gridExtra = TRUE
    ))

    expect_type(res, "list")
    expect_equal(length(res$plots), length(names(pbf)))
    expect_true(all(vapply(res$plots, function(x) inherits(x, "pheatmap"), logical(1))))

    subset_assays <- names(pbf)[2:1]
    res_subset <- suppressWarnings(plot_heatmap_diagnostic(
        pbf,
        pbf_name = subset_assays,
        sample_id_col = "FullRunName",
        factors_to_plot = c("MS_batch"),
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        return_gridExtra = TRUE
    ))

    expect_type(res_subset, "list")
    expect_equal(names(res_subset$plots), subset_assays)
    expect_true(all(vapply(res_subset$plots, function(x) inherits(x, "pheatmap"), logical(1))))
})

test_that("plot_PCA ProBatchFeatures returns ggplot for single assay", {
    pbf <- pb_test_make_pbf(n_rows = 40, n_cols = 6, add_log2 = TRUE)

    single_assay <- names(pbf)[1]
    res <- suppressWarnings(plot_PCA(
        pbf,
        pbf_name = single_assay,
        sample_id_col = "FullRunName",
        x_nPC = 2,
        y_nPC = 3
    ))

    expect_s3_class(res, "ggplot")
    expect_match(res$labels$x, "^PC2 ")
    expect_match(res$labels$y, "^PC3 ")
})

test_that("plot_heatmap_diagnostic ProBatchFeatures single assay returns pheatmap", {
    pbf <- pb_test_make_pbf(n_rows = 30, n_cols = 5, add_log2 = FALSE)

    res <- suppressWarnings(plot_heatmap_diagnostic(
        pbf,
        pbf_name = names(pbf)[1],
        sample_id_col = "FullRunName",
        factors_to_plot = c("MS_batch"),
        cluster_rows = FALSE,
        cluster_cols = FALSE
    ))

    expect_true(inherits(res, "pheatmap"))
})
