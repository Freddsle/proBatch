test_that("sample_mean_plots", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    matrix <- example_proteome_matrix[1:20, ]
    meanplot <- plot_sample_mean(matrix, example_sample_annotation,
        order_col = "order", batch_col = "MS_batch", color_by_batch = TRUE
    )

    expect_equal(meanplot$labels$colour, "MS_batch")
    expect_equal(meanplot$labels$x, "order")
    expect_equal(meanplot$labels$y, "Mean_Intensity")
})


test_that("boxplot_plots", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    proteome <- example_proteome[1:20, ]
    expect_warning(boxplot <- plot_boxplot(proteome, example_sample_annotation,
        batch_col = "MS_batch"
    ))

    expect_equal(boxplot$labels$fill, "MS_batch")
    expect_equal(boxplot$label$group, "order")
    expect_equal(boxplot$label$x, "order")
    expect_equal(boxplot$label$y, "Intensity")

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

    expect_warning(boxplot <- plot_boxplot(
        df_long,
        sample_annotation,
        order_col = "order",
        batch_col = "MS_batch",
        outliers = FALSE
    ))

    expect_equal(boxplot$layers[[1]]$geom_params$outlier.shape, NA)
})
