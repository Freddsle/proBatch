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
        "EarTag Strain Sex RunDate RunTime digestion_batch "
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
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    short_df <- example_proteome[example_proteome$peptide_group_label %in%
        unique(example_proteome$peptide_group_label)[1:3], ]

    expect_warning(
        loess_fit <- adjust_batch_trend_df(short_df, example_sample_annotation, span = 0.7),
        "imputed value flag column is NULL, changing no_fit_imputed to FALSE"
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
