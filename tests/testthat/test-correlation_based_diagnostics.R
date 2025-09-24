test_that("corr_matrix_plots", {
    data(example_proteome_matrix, package = "proBatch")
    peptides <- c("10231_QDVDVWLWQQEGSSK_2", "10768_RLESELDGLR_2")

    matrix_test <- example_proteome_matrix[peptides, ]
    corr_matrix <- cor(t(matrix_test), use = "complete.obs")
    expect_warning(
        corr_matrix_pheatmap <- plot_corr_matrix(corr_matrix),
        "annotation_row and / or annotation_col are not specified for heatmap"
    )
    expect_s3_class(corr_matrix_pheatmap, "pheatmap")
    expect_equal(corr_matrix_pheatmap$gtable$layout$name[4], "legend", ignore_attr = TRUE)
})


test_that("protein_corrplot_plots", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_peptide_annotation, package = "proBatch")

    expect_warning(
        color_list <- sample_annotation_to_colors(
            sample_annotation = example_peptide_annotation,
            factor_columns = "Gene"
        ),
        "The following columns will not be mapped to colors: peptide_group_label ProteinName ;"
    )

    corrplot <- plot_protein_corrplot(
        example_proteome_matrix,
        protein_name = "Haao",
        peptide_annotation = example_peptide_annotation,
        protein_col = "Gene",
        cluster_rows = TRUE, cluster_cols = TRUE,
        color_list = color_list
    )

    expect_equal(corrplot$tree_row$method, "complete", ignore_attr = TRUE)
    expect_equal(corrplot$tree_row$dist.method, "euclidean", ignore_attr = TRUE)

    expect_equal(corrplot$tree_row$labels[1], "10231_QDVDVWLWQQEGSSK_2", ignore_attr = TRUE)
    expect_equal(corrplot$tree_row$labels[2], "10768_RLESELDGLR_2", ignore_attr = TRUE)
    expect_equal(corrplot$tree_row$labels[3], "1131_AQGSVALSVTQDPAR_2", ignore_attr = TRUE)
})


test_that("sample_corr_heatmap", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    specified_samples <- example_sample_annotation$FullRunName[
        which(example_sample_annotation$order %in% 110:115)
    ]

    expect_warning(
        sample_heatmap <- plot_sample_corr_heatmap(
            example_proteome_matrix,
            samples_to_plot = specified_samples,
            cluster_rows = TRUE, cluster_cols = TRUE,
            annotation_names_col = TRUE, annotation_legend = FALSE,
            show_colnames = FALSE
        ), 
        "annotation_row and / or annotation_col are not specified for heatmap"
    )

    expect_equal(sample_heatmap$tree_row$method, "complete", ignore_attr = TRUE)
    expect_equal(sample_heatmap$tree_row$dist.method, "euclidean", ignore_attr = TRUE)

    expect_equal(sample_heatmap$tree_col$method, "complete", ignore_attr = TRUE)
    expect_equal(sample_heatmap$tree_col$dist.method, "euclidean", ignore_attr = TRUE)

    expect_equal(sample_heatmap$tree_row$labels[1], "I171013_BXD61_CD_ET2145_Run113", ignore_attr = TRUE)
    expect_equal(sample_heatmap$tree_row$labels[2], "I171013_BXD70_HF_ET1728_Run114", ignore_attr = TRUE)
    expect_equal(sample_heatmap$tree_row$labels[3], "I171016_BXD89_CD_ET2078_Run115", ignore_attr = TRUE)

    expect_equal(sample_heatmap$gtable$layout[[1]], c(1, 2, 4, 4, 4, 3), ignore_attr = TRUE)
    expect_equal(sample_heatmap$gtable$layout$name[[4]], "matrix", ignore_attr = TRUE)
})


test_that("sample_distribution_plot", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    matrix_test <- example_proteome_matrix[1:20, ]
    sample_dist <- plot_sample_corr_distribution(matrix_test,
        example_sample_annotation,
        batch_col = "MS_batch",
        biospecimen_id_col = "EarTag",
        plot_param = "batch_replicate"
    )

    expect_equal(sample_dist$labels$x, "batch_replicate", ignore_attr = TRUE)
    expect_equal(sample_dist$labels$y, "correlation", ignore_attr = TRUE)

    expect_s3_class(sample_dist$plot_env$corr_distribution, "data.frame")
    expect_equal(sample_dist$plot_env$plot_param, "batch_replicate", ignore_attr = TRUE)
    expect_s3_class(sample_dist$plot_env$gg, "ggplot")
})

test_that("calculate_sample_corr_distribution", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    matrix_test <- example_proteome_matrix[1:20, ]
    corr_distribution <- calculate_sample_corr_distr(
        data_matrix = matrix_test,
        repeated_samples = NULL,
        sample_annotation = example_sample_annotation,
        biospecimen_id_col = "EarTag",
        sample_id_col = "FullRunName",
        batch_col = "MS_batch"
    )

    expect_s3_class(corr_distribution, "data.frame")

    sample_cols <- paste("FullRunName", seq_len(2), sep = "_")

    expect_equal("batch_replicate" %in% names(corr_distribution), TRUE, ignore_attr = TRUE)
    expect_equal("correlation" %in% names(corr_distribution), TRUE, ignore_attr = TRUE)
    expect_equal("replicate" %in% names(corr_distribution), TRUE, ignore_attr = TRUE)
    expect_equal(all(sample_cols %in% names(corr_distribution)), TRUE, ignore_attr = TRUE)
})


test_that("peptide_distribution_plots", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_peptide_annotation, package = "proBatch")

    matrix_test <- example_proteome_matrix[1:20, ]
    peptide_dist <- plot_peptide_corr_distribution(
        data_matrix = matrix_test,
        peptide_annotation = example_peptide_annotation,
        protein_col = "Gene"
    )

    expect_equal(peptide_dist$labels$x, NULL, ignore_attr = TRUE)
    expect_equal(peptide_dist$labels$y, "correlation", ignore_attr = TRUE)

    expect_s3_class(peptide_dist$plot_env$corr_distribution, "data.frame")
    expect_equal(peptide_dist$plot_env$median_same_prot, 0.7337642, tolerance = 1e-6)
})
