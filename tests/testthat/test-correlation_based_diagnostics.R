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
        "The following columns will not be mapped to colors: peptide_group_label, ProteinName"
    )
    suppressWarnings(
        corrplot <- plot_protein_corrplot(
            example_proteome_matrix,
            protein_name = "Haao",
            peptide_annotation = example_peptide_annotation,
            protein_col = "Gene",
            cluster_rows = TRUE, cluster_cols = TRUE,
            color_list = color_list
        )
    )

    expect_equal(corrplot$tree_row$method, "complete", ignore_attr = TRUE)
    expect_equal(corrplot$tree_row$dist.method, "euclidean", ignore_attr = TRUE)

    expect_equal(corrplot$tree_row$labels[1], "10231_QDVDVWLWQQEGSSK_2", ignore_attr = TRUE)
    expect_equal(corrplot$tree_row$labels[2], "10768_RLESELDGLR_2", ignore_attr = TRUE)
    expect_equal(corrplot$tree_row$labels[3], "1131_AQGSVALSVTQDPAR_2", ignore_attr = TRUE)
})


test_that("sample_corr_heatmap", {
    pb_test_load_example_data()

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
    pb_test_load_example_data()

    matrix_test <- example_proteome_matrix[1:20, ]
    sample_dist <- plot_sample_corr_distribution(matrix_test,
        example_sample_annotation,
        batch_col = "MS_batch",
        biospecimen_id_col = "EarTag",
        plot_param = "batch_replicate"
    )

    expect_equal(as_label(sample_dist$mapping$x), "batch_replicate", ignore_attr = TRUE)
    expect_equal(as_label(sample_dist$mapping$y), "correlation", ignore_attr = TRUE)

    expect_s3_class(sample_dist$plot_env$corr_distribution, "data.frame")
    expect_equal(sample_dist$plot_env$plot_param, "batch_replicate", ignore_attr = TRUE)
    expect_s3_class(sample_dist$plot_env$gg, "ggplot")
})

test_that("calculate_sample_corr_distribution", {
    pb_test_load_example_data()

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

    expect_equal(as_label(peptide_dist$mapping$x), "same_protein", ignore_attr = TRUE)
    expect_equal(as_label(peptide_dist$mapping$y), "correlation", ignore_attr = TRUE)

    expect_s3_class(peptide_dist$plot_env$corr_distribution, "data.frame")
    expect_equal(peptide_dist$plot_env$median_same_prot, 0.7337642, tolerance = 1e-6)
})

make_corr_pbf_fixture <- function() {
    pb_test_load_example_data()
    data(example_peptide_annotation, package = "proBatch")

    matrix_small <- example_proteome_matrix[1:120, 1:20]
    sample_ann <- example_sample_annotation[
        match(colnames(matrix_small), example_sample_annotation$FullRunName),
    ]
    pbf <- suppressMessages(ProBatchFeatures(
        data_matrix = matrix_small,
        sample_annotation = sample_ann,
        sample_id_col = "FullRunName",
        name = "feature::raw"
    ))
    pbf <- suppressMessages(pb_transform(
        pbf,
        from = "feature::raw",
        steps = "log2",
        store_fast_steps = TRUE
    ))

    assay_name <- pb_current_assay(pbf)
    matrix_log <- suppressMessages(pb_assay_matrix(pbf, assay = assay_name))
    peptide_ann <- example_peptide_annotation[
        match(rownames(matrix_log), example_peptide_annotation$peptide_group_label),
    ]
    peptide_ann <- peptide_ann[!is.na(peptide_ann$peptide_group_label), , drop = FALSE]

    list(
        pbf = pbf,
        matrix_small = matrix_small,
        sample_ann = sample_ann,
        peptide_ann = peptide_ann
    )
}

test_that("sample correlation diagnostics accept PBF with rowname-only annotation", {
    fixture <- make_corr_pbf_fixture()
    pbf <- fixture$pbf
    matrix_small <- fixture$matrix_small
    sample_ann <- fixture$sample_ann

    sample_ann_row <- sample_ann
    rownames(sample_ann_row) <- sample_ann_row$FullRunName
    sample_ann_row$FullRunName <- NULL

    sample_corr_df <- calculate_sample_corr_distr(
        pbf,
        sample_annotation = sample_ann_row,
        batch_col = "MS_batch",
        biospecimen_id_col = "EarTag"
    )
    expect_s3_class(sample_corr_df, "data.frame")
    expect_true("correlation" %in% names(sample_corr_df))

    sample_corr_plot <- plot_sample_corr_distribution(
        pbf,
        sample_annotation = sample_ann_row,
        batch_col = "MS_batch",
        biospecimen_id_col = "EarTag",
        plot_param = "batch_replicate"
    )
    expect_s3_class(sample_corr_plot, "ggplot")

    sample_corr_heatmap <- plot_sample_corr_heatmap(
        pbf,
        sample_annotation = sample_ann_row,
        samples_to_plot = colnames(matrix_small)[1:6],
        sample_id_col = "FullRunName"
    )
    expect_s3_class(sample_corr_heatmap, "pheatmap")
})

test_that("calculate_sample_corr_distr() supports explicit PBF assay selection", {
    fixture <- make_corr_pbf_fixture()
    pbf <- fixture$pbf
    sample_ann <- fixture$sample_ann

    current_assay <- pb_current_assay(pbf)
    raw_assay <- names(pbf)[1]

    corr_default <- calculate_sample_corr_distr(
        pbf,
        sample_annotation = sample_ann,
        batch_col = "MS_batch",
        biospecimen_id_col = "EarTag"
    )
    corr_current <- calculate_sample_corr_distr(
        pbf,
        sample_annotation = sample_ann,
        batch_col = "MS_batch",
        biospecimen_id_col = "EarTag",
        pbf_name = current_assay
    )
    expect_equal(corr_default$correlation, corr_current$correlation)

    corr_raw <- calculate_sample_corr_distr(
        pbf,
        sample_annotation = sample_ann,
        batch_col = "MS_batch",
        biospecimen_id_col = "EarTag",
        pbf_name = raw_assay
    )
    corr_raw_expected <- calculate_sample_corr_distr(
        pb_assay_matrix(pbf, assay = raw_assay),
        sample_annotation = sample_ann,
        batch_col = "MS_batch",
        biospecimen_id_col = "EarTag"
    )
    expect_equal(corr_raw$correlation, corr_raw_expected$correlation)
})

test_that("plot_sample_corr_distribution() supports explicit PBF assay selection", {
    fixture <- make_corr_pbf_fixture()
    pbf <- fixture$pbf
    sample_ann <- fixture$sample_ann

    current_assay <- pb_current_assay(pbf)
    raw_assay <- names(pbf)[1]

    plot_default <- plot_sample_corr_distribution(
        pbf,
        sample_annotation = sample_ann,
        batch_col = "MS_batch",
        biospecimen_id_col = "EarTag",
        plot_param = "batch_replicate"
    )
    plot_current <- plot_sample_corr_distribution(
        pbf,
        sample_annotation = sample_ann,
        batch_col = "MS_batch",
        biospecimen_id_col = "EarTag",
        plot_param = "batch_replicate",
        pbf_name = current_assay
    )
    expect_equal(
        plot_default$plot_env$corr_distribution$correlation,
        plot_current$plot_env$corr_distribution$correlation
    )

    plot_raw <- plot_sample_corr_distribution(
        pbf,
        sample_annotation = sample_ann,
        batch_col = "MS_batch",
        biospecimen_id_col = "EarTag",
        plot_param = "batch_replicate",
        pbf_name = raw_assay
    )
    plot_raw_expected <- plot_sample_corr_distribution(
        pb_assay_matrix(pbf, assay = raw_assay),
        sample_annotation = sample_ann,
        batch_col = "MS_batch",
        biospecimen_id_col = "EarTag",
        plot_param = "batch_replicate"
    )
    expect_equal(
        plot_raw$plot_env$corr_distribution$correlation,
        plot_raw_expected$plot_env$corr_distribution$correlation
    )
})

test_that("peptide correlation diagnostics accept PBF with rowname-only annotation", {
    fixture <- make_corr_pbf_fixture()
    pbf <- fixture$pbf
    peptide_ann <- fixture$peptide_ann
    rownames(peptide_ann) <- peptide_ann$peptide_group_label
    peptide_ann$peptide_group_label <- NULL

    peptide_corr_df <- calculate_peptide_corr_distr(
        pbf,
        peptide_annotation = peptide_ann,
        protein_col = "Gene"
    )
    expect_s3_class(peptide_corr_df, "data.frame")
    expect_true("same_protein" %in% names(peptide_corr_df))

    peptide_corr_plot <- plot_peptide_corr_distribution(
        pbf,
        peptide_annotation = peptide_ann,
        protein_col = "Gene"
    )
    expect_s3_class(peptide_corr_plot, "ggplot")

    gene_counts <- table(peptide_corr_df$Gene1)
    gene_target <- names(gene_counts[gene_counts >= 2])[1]
    if (!length(gene_target) || is.na(gene_target)) {
        skip("Need at least one protein with >=2 peptides for correlation heatmap.")
    }

    protein_corrplot <- suppressWarnings(plot_protein_corrplot(
        pbf,
        protein_name = gene_target,
        peptide_annotation = peptide_ann,
        protein_col = "Gene",
        cluster_rows = TRUE,
        cluster_cols = TRUE
    ))
    expect_s3_class(protein_corrplot, "pheatmap")
})
