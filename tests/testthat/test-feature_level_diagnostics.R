test_that("single_feature_plot", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    expect_warning(
        single_feature <- plot_single_feature(
            feature_name = "46213_NVGVSFYADKPEVTQEQK_2",
            df_long = example_proteome,
            sample_annotation = example_sample_annotation
        ),
        "inferring order-related batch borders for a plot;"
    )

    expect_equal(single_feature$plot_env$feature_name, "46213_NVGVSFYADKPEVTQEQK_2")

    expect_equal(as_label(single_feature$mapping$x), "order")
    expect_equal(as_label(single_feature$mapping$y), "Intensity")

    expect_equal(single_feature$plot_env$batch_col, "MS_batch")
    expect_equal(single_feature$plot_env$color_by_batch, FALSE)
    expect_equal(single_feature$plot_env$order_col, "order")
    expect_equal(single_feature$plot_env$vline_color, "red")
})


test_that("peptides_of_one_protein_plot", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")
    data(example_peptide_annotation, package = "proBatch")

    expect_warning(
        color_scheme <- sample_annotation_to_colors(
            example_sample_annotation,
            sample_id_col = "FullRunName",
            factor_columns = c("MS_batch", "Diet")
        ),
        "The following columns will not be mapped to colors: EarTag, Strain, Sex, RunDate, RunTime, digestion_batch;"
    )

    expect_warning(
        peptides_plot <- plot_peptides_of_one_protein(
            protein_name = "Haao",
            peptide_annotation = example_peptide_annotation,
            protein_col = "Gene", df_long = example_proteome,
            sample_annotation = example_sample_annotation,
            color_by_batch = TRUE,
            order_col = "order", sample_id_col = "FullRunName",
            batch_col = "MS_batch",
            color_scheme = color_scheme
        ),
        "inferring order-related batch borders for a plot;"
    )

    expect_equal(peptides_plot$plot_env$feature_name[1], "10231_QDVDVWLWQQEGSSK_2")
    expect_equal(peptides_plot$plot_env$feature_name[2], "10768_RLESELDGLR_2")

    expect_equal(as_label(peptides_plot$mapping$x), "order")
    expect_equal(as_label(peptides_plot$mapping$y), "Intensity")

    expect_equal(peptides_plot$plot_env$batch_col, "MS_batch")
    expect_equal(peptides_plot$plot_env$color_by_batch, TRUE)
    expect_equal(peptides_plot$plot_env$order_col, "order")
    expect_equal(peptides_plot$plot_env$vline_color, "red")
})


test_that("spike_in_peptides_plot", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")
    data(example_peptide_annotation, package = "proBatch")

    spike_in <- plot_spike_in(
        spike_ins = "BOVINE_A1ag",
        peptide_annotation = example_peptide_annotation,
        protein_col = "Gene", df_long = example_proteome,
        sample_annotation = example_sample_annotation,
        plot_title = "Spike-in BOVINE protein peptides",
        vline_color = NULL
    )

    expect_equal(spike_in$plot_env$feature_name[1], "10062_NVGVSFYADKPEVTQEQK_3")
    expect_equal(spike_in$plot_env$feature_name[2], "10063_NVGVSFYADKPEVTQEQKK_3")

    expect_equal(as_label(spike_in$mapping$x), "order")
    expect_equal(as_label(spike_in$mapping$y), "Intensity")

    expect_equal(spike_in$plot_env$batch_col, "MS_batch")
    expect_equal(spike_in$plot_env$color_by_batch, FALSE)
    expect_equal(spike_in$plot_env$order_col, "order")
    expect_equal(spike_in$plot_env$vline_color, NULL)
})


test_that("iRT_peptides_plot", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")
    data(example_peptide_annotation, package = "proBatch")

    expect_warning(
        iRT <- plot_iRT(
            irt_pattern = "iRT",
            peptide_annotation = example_peptide_annotation,
            protein_col = "Gene",
            df_long = example_proteome,
            sample_annotation = example_sample_annotation
        ),
        "inferring order-related batch borders for a plot;"
    )

    expect_equal(iRT$plot_env$feature_name[1], "1146_ADVTPADFSEWSK_3")
    expect_equal(iRT$plot_env$feature_name[2], "12476_TPVISGGPYEYR_2")

    expect_equal(as_label(iRT$mapping$x), "order")
    expect_equal(as_label(iRT$mapping$y), "Intensity")

    expect_equal(iRT$plot_env$batch_col, "MS_batch")
    expect_equal(iRT$plot_env$color_by_batch, FALSE)
    expect_equal(iRT$plot_env$order_col, "order")
    expect_equal(iRT$plot_env$vline_color, "red")
})


test_that("fitting_trend_plots", {
    data(example_proteome, package = "proBatch")
    pb_test_load_example_data()

    short_df <- example_proteome[example_proteome$peptide_group_label %in%
        unique(example_proteome$peptide_group_label)[1:3], ]

    expect_warning(
        loess_fit <- adjust_batch_trend_df(short_df, example_sample_annotation, span = 0.7),
        "`qual_col` is NULL, setting `no_fit_imputed = FALSE` so imputed flags are ignored."
    )

    expect_warning(
        fit_plot <- plot_with_fitting_curve(
            feature_name = "10062_NVGVSFYADKPEVTQEQK_3",
            fit_df = loess_fit, fit_value_col = "fit",
            df_long = example_proteome,
            sample_annotation = example_sample_annotation
        ),
        "inferring order-related batch borders for a plot;"
    )

    expect_equal(fit_plot$plot_env$feature_name[1], "10062_NVGVSFYADKPEVTQEQK_3")

    expect_equal(as_label(fit_plot$mapping$x), "order")
    expect_equal(as_label(fit_plot$mapping$y), "Intensity")

    expect_equal(fit_plot$plot_env$batch_col, "MS_batch")
    expect_equal(fit_plot$plot_env$color_by_batch, FALSE)
    expect_equal(fit_plot$plot_env$order_col, "order")
    expect_equal(fit_plot$plot_env$vline_color, "grey")
})

make_feature_pbf_fixture <- function() {
    pb_test_load_example_data()
    data(example_peptide_annotation, package = "proBatch")

    matrix_small <- example_proteome_matrix[1:100, 1:16]
    sample_ann <- example_sample_annotation[
        match(colnames(matrix_small), example_sample_annotation$FullRunName),
    ]
    peptide_ann <- example_peptide_annotation[
        match(rownames(matrix_small), example_peptide_annotation$peptide_group_label),
    ]
    peptide_ann <- peptide_ann[!is.na(peptide_ann$peptide_group_label), , drop = FALSE]

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

    list(
        pbf = pbf,
        matrix_small = matrix_small,
        sample_ann = sample_ann,
        peptide_ann = peptide_ann
    )
}

test_that("feature-level diagnostics accept PBF and rowname-based annotation", {
    fixture <- make_feature_pbf_fixture()
    pbf <- fixture$pbf
    matrix_small <- fixture$matrix_small
    sample_ann <- fixture$sample_ann
    peptide_ann <- fixture$peptide_ann

    sample_ann_row <- sample_ann
    rownames(sample_ann_row) <- sample_ann_row$FullRunName
    sample_ann_row$FullRunName <- NULL

    feature_id <- rownames(matrix_small)[1]
    single_feature <- suppressWarnings(plot_single_feature(
        feature_name = feature_id,
        df_long = pbf,
        sample_annotation = sample_ann_row
    ))
    expect_s3_class(single_feature, "ggplot")

    gene_counts <- table(peptide_ann$Gene)
    gene_target <- names(gene_counts[gene_counts >= 2])[1]
    if (!length(gene_target) || is.na(gene_target)) {
        skip("Need at least one protein with >=2 peptides for protein panel test.")
    }

    rownames(peptide_ann) <- peptide_ann$peptide_group_label
    peptide_ann$peptide_group_label <- NULL

    protein_plot <- suppressWarnings(plot_peptides_of_one_protein(
        protein_name = gene_target,
        peptide_annotation = peptide_ann,
        protein_col = "Gene",
        df_long = pbf,
        sample_annotation = sample_ann_row
    ))
    expect_s3_class(protein_plot, "ggplot")
})

test_that("feature-level single/protein plots support explicit PBF assay selection", {
    fixture <- make_feature_pbf_fixture()
    pbf <- fixture$pbf
    sample_ann <- fixture$sample_ann
    peptide_ann <- fixture$peptide_ann

    current_assay <- pb_current_assay(pbf)
    raw_assay <- names(pbf)[1]
    feature_id <- rownames(fixture$matrix_small)[1]

    single_default <- suppressWarnings(plot_single_feature(
        feature_name = feature_id,
        df_long = pbf,
        sample_annotation = sample_ann
    ))
    single_current <- suppressWarnings(plot_single_feature(
        feature_name = feature_id,
        df_long = pbf,
        sample_annotation = sample_ann,
        pbf_name = current_assay
    ))
    expect_equal(single_default$data$Intensity, single_current$data$Intensity)

    raw_long <- matrix_to_long(
        data_matrix = pb_assay_matrix(pbf, assay = raw_assay),
        feature_id_col = "peptide_group_label",
        sample_id_col = "FullRunName",
        measure_col = "Intensity"
    )
    single_raw <- suppressWarnings(plot_single_feature(
        feature_name = feature_id,
        df_long = pbf,
        sample_annotation = sample_ann,
        pbf_name = raw_assay
    ))
    single_raw_expected <- suppressWarnings(plot_single_feature(
        feature_name = feature_id,
        df_long = raw_long,
        sample_annotation = sample_ann
    ))
    expect_equal(single_raw$data$Intensity, single_raw_expected$data$Intensity)

    gene_counts <- table(peptide_ann$Gene)
    gene_target <- names(gene_counts[gene_counts >= 2])[1]
    if (!length(gene_target) || is.na(gene_target)) {
        skip("Need at least one protein with >=2 peptides for assay selection test.")
    }

    protein_default <- suppressWarnings(plot_peptides_of_one_protein(
        protein_name = gene_target,
        peptide_annotation = peptide_ann,
        protein_col = "Gene",
        df_long = pbf,
        sample_annotation = sample_ann
    ))
    protein_current <- suppressWarnings(plot_peptides_of_one_protein(
        protein_name = gene_target,
        peptide_annotation = peptide_ann,
        protein_col = "Gene",
        df_long = pbf,
        sample_annotation = sample_ann,
        pbf_name = current_assay
    ))
    expect_equal(protein_default$data$Intensity, protein_current$data$Intensity)

    protein_raw <- suppressWarnings(plot_peptides_of_one_protein(
        protein_name = gene_target,
        peptide_annotation = peptide_ann,
        protein_col = "Gene",
        df_long = pbf,
        sample_annotation = sample_ann,
        pbf_name = raw_assay
    ))
    protein_raw_expected <- suppressWarnings(plot_peptides_of_one_protein(
        protein_name = gene_target,
        peptide_annotation = peptide_ann,
        protein_col = "Gene",
        df_long = raw_long,
        sample_annotation = sample_ann
    ))
    expect_equal(protein_raw$data$Intensity, protein_raw_expected$data$Intensity)
})

test_that("fitting diagnostics accept ProBatchFeatures inputs", {
    fixture <- make_feature_pbf_fixture()
    pbf <- fixture$pbf
    matrix_small <- fixture$matrix_small
    sample_ann <- fixture$sample_ann

    expect_warning(
        fit_df <- adjust_batch_trend_df(
            pbf,
            sample_annotation = sample_ann,
            span = 0.7
        ),
        "`qual_col` is NULL, setting `no_fit_imputed = FALSE` so imputed flags are ignored."
    )
    expect_s3_class(fit_df, "data.frame")

    feature_id <- rownames(matrix_small)[1]
    fit_df <- fit_df[fit_df$peptide_group_label %in% feature_id, , drop = FALSE]
    fit_plot <- suppressWarnings(plot_with_fitting_curve(
        feature_name = feature_id,
        fit_df = fit_df,
        df_long = pbf,
        sample_annotation = sample_ann
    ))
    expect_s3_class(fit_plot, "ggplot")
})
