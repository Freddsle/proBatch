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

    expect_warning(
        expect_warning(
            expect_warning(
                hiearchical <- plot_hierarchical_clustering(matrix_test,
                    sample_annotation = example_sample_annotation,
                    factors_to_plot = c("MS_batch", "Diet"),
                    distance = "euclidean", agglomeration = "complete",
                    color_list = color_list,
                    label_samples = FALSE
                ),
                "filling missing values with 0"
            ),
            "Hierarchical clustering cannot operate with missing values"
        ),
        "color list and sample annotation have different factors,"
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

    expect_warning(
        expect_warning(
            expect_warning(
                pca <- plot_PCA(
                    example_proteome_matrix, example_sample_annotation,
                    color_by = "MS_batch", plot_title = "PCA colored by MS batch",
                    fill_the_missing = -1
                ),
                "filling missing values with -1"
            ),
            "PCA cannot operate with missing values in the matrix"
        ), "color_scheme will be inferred automatically"
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

test_that("plot_PCA supports Bayesian PCA on incomplete matrices", {
    skip_if_not_installed("pcaMethods")

    pb_test_load_example_data()

    color_list <- sample_annotation_to_colors(
        example_sample_annotation,
        sample_id_col = "FullRunName",
        factor_columns = "MS_batch"
    )

    warnings_seen <- character()
    pca_bpca <- withCallingHandlers(
        plot_PCA(
            example_proteome_matrix,
            example_sample_annotation,
            color_by = "MS_batch",
            color_scheme = color_list[["MS_batch"]],
            pca_method = "bpca",
            bpca_nPcs = 3,
            x_nPC = 2,
            y_nPC = 3
        ),
        warning = function(w) {
            warnings_seen <<- c(warnings_seen, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )

    expect_s3_class(pca_bpca, "ggplot")
    expect_match(pca_bpca$labels$x, "^PC2")
    expect_match(pca_bpca$labels$y, "^PC3")
    expect_false(any(grepl(
        "PCA cannot operate with missing values in the matrix",
        warnings_seen,
        fixed = TRUE
    )))
    expect_false(any(grepl(
        "filling missing values with -1",
        warnings_seen,
        fixed = TRUE
    )))
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

test_that("plot_TSNE returns ggplot by default", {
    skip_if_not_installed("Rtsne")

    pb_test_load_example_data()

    tsne_plot <- plot_TSNE(
        example_proteome_matrix, example_sample_annotation,
        color_by = "MS_batch",
        fill_the_missing = -1,
        perplexity = 2,
        max_iter = 250,
        random_seed = 123
    )

    expect_s3_class(tsne_plot, "ggplot")
})

test_that("plot_TSNE use_plotlyrender returns plotly object", {
    skip_if_not_installed("plotly")
    skip_if_not_installed("Rtsne")

    pb_test_load_example_data()

    expect_warning(
        expect_warning(
            tsne_plot <- plot_TSNE(
                example_proteome_matrix, example_sample_annotation,
                color_by = "MS_batch",
                fill_the_missing = -1,
                perplexity = 2,
                max_iter = 250,
                random_seed = 123,
                use_plotlyrender = TRUE
            ),
            "filling missing values with -1"
        ),
        "t-SNE cannot operate with missing values in the matrix"
    )

    expect_s3_class(tsne_plot, "plotly")
})

test_that("plot_UMAP returns ggplot by default", {
    skip_if_not_installed("umap")

    pb_test_load_example_data()

    expect_warning(
        expect_warning(
            expect_warning(
                umap_plot <- plot_UMAP(
                    example_proteome_matrix, example_sample_annotation,
                    color_by = "MS_batch",
                    fill_the_missing = -1,
                    n_neighbors = 5,
                    random_state = 42
                ),
                "filling missing values with -1"
            ),
            "UMAP cannot operate with missing values in the matrix"
        ),
        "color_scheme will be inferred automatically"
    )

    expect_s3_class(umap_plot, "ggplot")
})

test_that("plot_UMAP use_plotlyrender returns plotly object", {
    skip_if_not_installed("plotly")
    skip_if_not_installed("umap")

    pb_test_load_example_data()

    expect_warning(
        expect_warning(
            umap_plot <- plot_UMAP(
                example_proteome_matrix, example_sample_annotation,
                color_by = "MS_batch",
                fill_the_missing = -1,
                n_neighbors = 5,
                random_state = 42,
                use_plotlyrender = TRUE
            ),
            "filling missing values with -1"
        ),
        "UMAP cannot operate with missing values in the matrix"
    )

    expect_s3_class(umap_plot, "plotly")
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

test_that("plot_TSNE ProBatchFeatures arranges multiple assays by default", {
    skip_if_not_installed("Rtsne")

    pbf <- pb_test_make_pbf(n_rows = 30, n_cols = 6, add_log2 = TRUE)

    res <- suppressWarnings(plot_TSNE(
        pbf,
        sample_id_col = "FullRunName",
        perplexity = 2,
        max_iter = 250,
        random_seed = 321,
        return_gridExtra = TRUE
    ))

    expect_type(res, "list")
    expect_named(res$plots, names(pbf))
    expect_true(all(vapply(res$plots, inherits, logical(1), "ggplot")))
})

test_that("plot_TSNE ProBatchFeatures use_plotlyrender returns plotly objects", {
    skip_if_not_installed("plotly")
    skip_if_not_installed("Rtsne")

    pbf <- pb_test_make_pbf(n_rows = 30, n_cols = 6, add_log2 = TRUE)

    res <- suppressWarnings(plot_TSNE(
        pbf,
        sample_id_col = "FullRunName",
        perplexity = 2,
        max_iter = 250,
        random_seed = 321,
        use_plotlyrender = TRUE
    ))

    expect_type(res, "list")
    expect_named(res, names(pbf))
    expect_true(all(vapply(res, inherits, logical(1), "plotly")))
})

test_that("plot_UMAP ProBatchFeatures arranges multiple assays by default", {
    skip_if_not_installed("umap")

    pbf <- pb_test_make_pbf(n_rows = 30, n_cols = 6, add_log2 = TRUE)

    res <- suppressWarnings(plot_UMAP(
        pbf,
        sample_id_col = "FullRunName",
        n_neighbors = 5,
        random_state = 99,
        return_gridExtra = TRUE
    ))

    expect_type(res, "list")
    expect_named(res$plots, names(pbf))
    expect_true(all(vapply(res$plots, inherits, logical(1), "ggplot")))
})

test_that("plot_UMAP ProBatchFeatures use_plotlyrender returns plotly objects", {
    skip_if_not_installed("plotly")
    skip_if_not_installed("umap")

    pbf <- pb_test_make_pbf(n_rows = 30, n_cols = 6, add_log2 = TRUE)

    res <- suppressWarnings(plot_UMAP(
        pbf,
        sample_id_col = "FullRunName",
        n_neighbors = 5,
        random_state = 99,
        use_plotlyrender = TRUE
    ))

    expect_type(res, "list")
    expect_named(res, names(pbf))
    expect_true(all(vapply(res, inherits, logical(1), "plotly")))
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
