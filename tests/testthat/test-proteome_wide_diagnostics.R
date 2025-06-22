test_that("hierarchical_clustering", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    matrix_test <- example_proteome_matrix[1:10, ]

    # exclude columns EarTag Strain Sex RunDate RunTime digestion_batch
    example_sample_annotation <- example_sample_annotation %>%
        dplyr::select("FullRunName", "MS_batch", "Diet")

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
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    matrix_test <- example_proteome_matrix[1:20, ]
    example_sample_annotation <- example_sample_annotation %>%
        dplyr::select("FullRunName", "MS_batch", "Sex", "digestion_batch", "Diet")


    color_list <- sample_annotation_to_colors(
        example_sample_annotation,
        sample_id_col = "FullRunName",
        factor_columns = c("MS_batch", "Sex", "digestion_batch", "Diet")
    )

    expect_warning(
        expect_warning(
            heatmap <- plot_heatmap_diagnostic(
                matrix_test,
                sample_annotation = example_sample_annotation,
                factors_to_plot = c("MS_batch", "Sex", "digestion_batch", "Diet"),
                cluster_cols = TRUE,
                show_rownames = TRUE, show_colnames = FALSE,
                color_list = color_list
            ),
            "filling missing values with -1"
        ),
        "Heatmap cannot operate with missing values in the matrix"
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


test_that("pvca_plot", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    matrix_test <- example_proteome_matrix[1:150, ]
    expect_warning(
        expect_warning(
            pvca <- plot_PVCA(matrix_test, example_sample_annotation,
                technical_factors = c("MS_batch", "digestion_batch"),
                biological_factors = c("Diet", "Sex", "Strain")
            ),
            "PVCA cannot operate with missing values in the matrix"
        ),
        "filling missing values with -1"
    )

    expect_equal(pvca$data$weights[1], 0.39166175, tolerance = 3e-2, ignore_attr = TRUE)
    expect_equal(as.character(pvca$data$label[3]), "MS_batch", ignore_attr = TRUE)
    expect_equal(as.character(pvca$data$label[2]), "Sex:Strain", ignore_attr = TRUE)

    expect_equal(pvca$data$category[1], "biological")
    expect_equal(pvca$data$category[3], "technical")
})


test_that("pca_plot", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

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
    expect_equal(pca$labels$x, "PC1 (69.5%)")
    expect_equal(pca$labels$colour, "MS_batch")
})
