test_that("sample_mean_plots", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    matrix <- example_proteome_matrix[1:20, ]
    expect_warning(
        expect_warning(
            meanplot <- plot_sample_mean(matrix, example_sample_annotation,
                order_col = "order", batch_col = "MS_batch", color_by_batch = TRUE
            ),
            "inferring order-related batch borders for a plot;"
        ),
        "color_scheme will be inferred automatically"
    )

    expect_equal(meanplot$labels$colour, "MS_batch")
    expect_equal(meanplot$labels$x, "order")
    expect_equal(meanplot$labels$y, "Mean_Intensity")
})


test_that("boxplot_plots", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    proteome <- example_proteome[1:20, ]
    expect_warning(
        expect_warning(
            boxplot <- plot_boxplot(proteome, example_sample_annotation,
                batch_col = "MS_batch"
            ),
            "color_scheme will be inferred automatically"
        ),
        "Mismatch between sample_annotation and df_long sample"
    )

    expect_equal(boxplot$labels$fill, "MS_batch")
    expect_equal(rlang::as_name(boxplot$mapping$group), "order")
    expect_equal(rlang::as_name(boxplot$mapping$x), "order")
    expect_equal(rlang::as_name(boxplot$mapping$y), "Intensity")

    expect_equal(boxplot$plot_env$color_by_batch, TRUE)
    expect_equal(boxplot$plot_env$facet_col, NULL)
})



test_that("mean plot adds vertical lines and y limits", {
    matrix <- matrix(1:9, nrow = 3)
    colnames(matrix) <- paste0("S", 1:3)
    sample_annotation <- data.frame(
        FullRunName = paste0("S", 1:3),
        MS_batch = factor(c("A", "B", "B")),
        order = 1:3
    )

    expect_warning(meanplot <- plot_sample_mean(
        matrix,
        sample_annotation,
        order_col = "order",
        batch_col = "MS_batch",
        vline_color = "blue",
        ylimits = c(0, 10)
    ))

    expect_true(any(sapply(meanplot$layers, function(x) inherits(x$geom, "GeomVline"))))
    expect_equal(meanplot$coordinates$limits$y, c(0, 10))
})


test_that("axis rotation when order is character", {
    matrix <- matrix(1:9, nrow = 3)
    colnames(matrix) <- paste0("S", 1:3)
    sample_annotation <- data.frame(
        FullRunName = paste0("S", 1:3),
        MS_batch = factor(c("A", "A", "B")),
        order = 1:3
    )

    expect_warning(meanplot <- plot_sample_mean(
        matrix,
        sample_annotation,
        order_col = "FullRunName",
        batch_col = "MS_batch"
    ))

    expect_equal(meanplot$theme$axis.text.x$angle, 90)
})


test_that("boxplot without outliers", {
    sample_annotation <- data.frame(
        FullRunName = paste0("S", 1:4),
        MS_batch = factor(c("A", "A", "B", "B")),
        order = 1:4
    )
    df_long <- data.frame(
        FullRunName = sample_annotation$FullRunName,
        Intensity = 1:4,
        peptide_group_label = paste0("p", 1:4)
    )
    df_long <- merge(df_long, sample_annotation, by = "FullRunName")

    color_scheme <- sample_annotation_to_colors(
        sample_annotation,
        factor_columns = c("MS_batch"),
        numeric_columns = c("order")
    )
    expect_warning(boxplot <- plot_boxplot(
        df_long,
        sample_annotation,
        order_col = "order",
        batch_col = "MS_batch",
        outliers = FALSE,
        color_scheme = color_scheme
    ), "outliers will be removed")

    expect_true(is.null(boxplot$layers[[1]]$geom_params$outlier.shape) ||
        is.na(boxplot$layers[[1]]$geom_params$outlier.shape))
})


test_that("plot_sample_mean with ProBatchFeatures", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    pbf <- ProBatchFeatures(
        data_matrix = log2(example_proteome_matrix + 1),
        sample_annotation = example_sample_annotation,
        sample_id_col = "FullRunName",
        name = "raw"
    )
    expect_warning(
        meanplot <- plot_sample_mean(pbf,
            sample_id_col = "FullRunName",
            batch_col = "MS_batch",
            pbf_name = "feature::raw"
        ),
        "inferring order-related batch borders for a plot"
    )
    expect_equal(meanplot$labels$x, "order")
    expect_equal(meanplot$labels$y, "Mean_Intensity")
    expect_equal(meanplot$plot_env$color_by_batch, FALSE)
    expect_equal(meanplot$plot_env$facet_col, NULL)
})

test_that("plot_boxplot ProBatchFeatures handles multiple assays", {
    skip_if_not_installed("gridExtra")
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    matrix_small <- example_proteome_matrix[1:30, 1:5]
    sample_ids <- colnames(matrix_small)
    sample_ann <- example_sample_annotation[match(sample_ids, example_sample_annotation$FullRunName), ]

    pbf <- suppressMessages(ProBatchFeatures(
        data_matrix = matrix_small,
        sample_annotation = sample_ann,
        sample_id_col = "FullRunName",
        name = "feature::raw"
    ))
    pbf <- suppressMessages(pb_transform(pbf,
        from = "feature::raw",
        steps = "log2",
        store_fast_steps = TRUE
    ))

    res <- suppressWarnings(plot_boxplot(
        pbf,
        sample_id_col = "FullRunName",
        return_gridExtra = TRUE
    ))

    expect_type(res, "list")
    expect_equal(length(res$plots), length(names(pbf)))
    expect_true(all(vapply(res$plots, inherits, logical(1), "ggplot")))
})

test_that("plot_boxplot ProBatchFeatures returns ggplot for single assay", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    matrix_small <- example_proteome_matrix[1:30, 1:5]
    sample_ids <- colnames(matrix_small)
    sample_ann <- example_sample_annotation[match(sample_ids, example_sample_annotation$FullRunName), ]

    pbf <- suppressMessages(ProBatchFeatures(
        data_matrix = matrix_small,
        sample_annotation = sample_ann,
        sample_id_col = "FullRunName",
        name = "feature::raw"
    ))

    res <- suppressWarnings(plot_boxplot(
        pbf,
        pbf_name = names(pbf)[1],
        sample_id_col = "FullRunName"
    ))

    expect_s3_class(res, "ggplot")
})

test_that("plot_boxplot ProBatchFeatures respects assay subset order", {
    skip_if_not_installed("gridExtra")
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    matrix_small <- example_proteome_matrix[1:30, 1:5]
    sample_ids <- colnames(matrix_small)
    sample_ann <- example_sample_annotation[match(sample_ids, example_sample_annotation$FullRunName), ]

    pbf <- suppressMessages(ProBatchFeatures(
        data_matrix = matrix_small,
        sample_annotation = sample_ann,
        sample_id_col = "FullRunName",
        name = "feature::raw"
    ))
    pbf <- suppressMessages(pb_transform(pbf,
        from = "feature::raw",
        steps = "log2",
        store_fast_steps = TRUE
    ))

    subset_assays <- names(pbf)[2:1]
    res <- suppressWarnings(plot_boxplot(
        pbf,
        pbf_name = subset_assays,
        sample_id_col = "FullRunName",
        return_gridExtra = TRUE
    ))

    expect_type(res, "list")
    expect_equal(names(res$plots), subset_assays)
    expect_true(all(vapply(res$plots, inherits, logical(1), "ggplot")))
})
