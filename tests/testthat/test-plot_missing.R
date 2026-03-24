toy_matrix <- matrix(
    c(
        10, NA, 12,
        7, 8, 9,
        NA, 5, 6
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("prot", 1:3), paste0("S", 1:3))
)

toy_matrix2 <- matrix(
    c(
        NA, 2, 3,
        4, NA, 6,
        7, 8, 9
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("prot", 1:3), paste0("S", 1:3))
)

toy_sa <- data.frame(
    FullRunName = colnames(toy_matrix),
    Condition = c("A", "B", "B"),
    Lab = c("L1", "L1", "L2"),
    Label = paste("Sample", seq_len(ncol(toy_matrix))),
    stringsAsFactors = FALSE
)
rownames(toy_sa) <- toy_sa$FullRunName

pbf_toy <- ProBatchFeatures(
    data_matrix = toy_matrix,
    sample_annotation = toy_sa,
    sample_id_col = "FullRunName",
    name = "toy"
)

toy_assay <- pb_current_assay(pbf_toy)
toy_assay_alt <- paste0(toy_assay, "_alt")

se_alt <- SummarizedExperiment(
    assays = list(intensity = toy_matrix2),
    colData = colData(pbf_toy[[toy_assay]])
)

pbf_multi <- proBatch:::.pb_add_assay_with_link(
    pbf_toy,
    se = se_alt,
    to = toy_assay_alt,
    from = toy_assay
)

test_that(".pb_truncate_heatmap_labels shortens long labels by default", {
    labels <- c(
        "short_label",
        "123456789012345",
        "1234567890123456",
        NA_character_
    )

    expect_equal(
        .pb_truncate_heatmap_labels(labels),
        c("short_label", "123456789012345", "1234567890...", NA_character_)
    )
    expect_null(.pb_truncate_heatmap_labels(NULL))
})

test_that("plot_NA_heatmap.default returns pheatmap", {
    skip_if_not_installed("pheatmap")

    res <- plot_NA_heatmap(
        toy_matrix,
        sample_annotation = toy_sa,
        color_by = "Condition",
        label_by = "Label",
        cluster_samples = FALSE,
        cluster_features = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        drop_complete = FALSE,
        draw = FALSE
    )

    expect_s3_class(res, "pheatmap")
})


test_that("plot_NA_heatmap.default supports multiple color_by columns", {
    skip_if_not_installed("pheatmap")

    res <- plot_NA_heatmap(
        toy_matrix,
        sample_annotation = toy_sa,
        color_by = c("Condition", "Label"),
        cluster_samples = FALSE,
        cluster_features = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        drop_complete = FALSE,
        draw = FALSE
    )

    expect_s3_class(res, "pheatmap")
})

test_that(".pb_group_missing_matrix aggregates observed fractions by group", {
    grouped <- .pb_group_missing_matrix(
        data_matrix = toy_matrix,
        sample_annotation = toy_sa,
        sample_id_col = "FullRunName",
        color_by = "Condition",
        drop_complete = FALSE
    )

    expected <- matrix(
        c(
            1, 1, 0,
            0.5, 1, 1
        ),
        nrow = 3,
        dimnames = list(
            rownames(toy_matrix),
            c("Condition=A", "Condition=B")
        )
    )

    expect_equal(grouped$matrix, expected)
    expect_equal(grouped$annotation$Condition, c("A", "B"))
    expect_equal(grouped$annotation$.group_size, c(1L, 2L))
    expect_equal(
        grouped$annotation$.group_label,
        c("Condition=A (n=1)", "Condition=B (n=2)")
    )
})

test_that(".pb_group_missing_matrix supports multi-column grouping", {
    grouped <- .pb_group_missing_matrix(
        data_matrix = toy_matrix,
        sample_annotation = toy_sa,
        sample_id_col = "FullRunName",
        color_by = c("Condition", "Lab"),
        drop_complete = TRUE
    )

    expected <- matrix(
        c(
            1, 0,
            0, 1,
            1, 1
        ),
        nrow = 2,
        dimnames = list(
            c("prot1", "prot3"),
            c(
                "Condition=A | Lab=L1",
                "Condition=B | Lab=L1",
                "Condition=B | Lab=L2"
            )
        )
    )

    expect_equal(grouped$matrix, expected)
    expect_equal(grouped$annotation$Condition, c("A", "B", "B"))
    expect_equal(grouped$annotation$Lab, c("L1", "L1", "L2"))
    expect_equal(grouped$annotation$.group_size, c(1L, 1L, 1L))
})

test_that("plot_grouped_NA_heatmap.default supports grouped observed fractions", {
    skip_if_not_installed("pheatmap")

    res <- plot_grouped_NA_heatmap(
        toy_matrix,
        sample_annotation = toy_sa,
        color_by = c("Condition", "Lab"),
        cluster_samples = FALSE,
        cluster_features = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        drop_complete = FALSE,
        draw = FALSE
    )

    expect_s3_class(res, "pheatmap")
})

test_that("plot_grouped_NA_heatmap.default requires color_by", {
    expect_error(
        plot_grouped_NA_heatmap(
            toy_matrix,
            sample_annotation = toy_sa,
            draw = FALSE
        ),
        "require one or more metadata columns in `color_by`"
    )
})


test_that("plot_NA_heatmap.ProBatchFeatures arranges multiple assays", {
    skip_if_not_installed("pheatmap")
    skip_if_not_installed("gridExtra")

    res <- plot_NA_heatmap(
        pbf_multi,
        pbf_name = c(toy_assay, toy_assay_alt),
        color_by = "Condition",
        label_by = "FullRunName",
        cluster_samples = FALSE,
        cluster_features = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        drop_complete = FALSE,
        draw = FALSE
    )

    expect_type(res, "list")
    expect_s3_class(res$grob, "gtable")
    expect_length(res$heatmaps, 2L)
    expect_equal(names(res$heatmaps), c(toy_assay, toy_assay_alt))
    expect_true(all(vapply(res$heatmaps, inherits, logical(1), what = "pheatmap")))
})

test_that("plot_NA_heatmap.ProBatchFeatures accepts explicit main without collision", {
    skip_if_not_installed("pheatmap")

    res <- plot_NA_heatmap(
        pbf_toy,
        color_by = "Condition",
        label_by = "FullRunName",
        cluster_samples = FALSE,
        cluster_features = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        drop_complete = FALSE,
        draw = FALSE,
        main = "Custom NA heatmap title"
    )

    expect_s3_class(res, "pheatmap")
})

test_that("plot_grouped_NA_heatmap.ProBatchFeatures supports grouped observed fractions", {
    skip_if_not_installed("pheatmap")
    skip_if_not_installed("gridExtra")

    res <- plot_grouped_NA_heatmap(
        pbf_multi,
        pbf_name = c(toy_assay, toy_assay_alt),
        color_by = c("Condition", "Lab"),
        cluster_samples = FALSE,
        cluster_features = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        drop_complete = FALSE,
        draw = FALSE
    )

    expect_type(res, "list")
    expect_s3_class(res$grob, "gtable")
    expect_length(res$heatmaps, 2L)
    expect_true(all(vapply(res$heatmaps, inherits, logical(1), what = "pheatmap")))
})


test_that("plot_NA_density.default summarises intensity distribution", {
    density_plot <- plot_NA_density(toy_matrix)

    expect_s3_class(density_plot, "ggplot")
    expect_equal(
        sort(unique(density_plot$data$Type)),
        sort(c("Missing Value", "Valid Value"))
    )
    expect_equal(density_plot$labels$x, "Intensity")
    expect_equal(density_plot$labels$y, "Density")
})


test_that("plot_NA_density.ProBatchFeatures facets assays", {
    density_plot <- plot_NA_density(
        pbf_multi,
        pbf_name = c(toy_assay, toy_assay_alt),
        nrow = 1,
        ncol = 2
    )

    expect_s3_class(density_plot, "ggplot")
    expect_setequal(
        unique(as.character(density_plot$data$pbf_name)),
        c(toy_assay, toy_assay_alt)
    )
})


test_that("plot_NA_frequency.default summarises observation counts", {
    frequency_plot <- plot_NA_frequency(toy_matrix)

    expect_s3_class(frequency_plot, "ggplot")
    freq_data <- frequency_plot$data
    freq_data <- freq_data[order(freq_data$valid_counts), , drop = FALSE]
    expected <- data.frame(
        valid_counts = c(2L, 3L),
        Freq = c(2L, 1L)
    )
    expect_equal(freq_data, expected, ignore_attr = TRUE)
})


test_that("plot_NA_frequency.ProBatchFeatures facets assays", {
    frequency_plot <- plot_NA_frequency(
        pbf_multi,
        pbf_name = c(toy_assay, toy_assay_alt),
        nrow = 1,
        ncol = 2
    )

    expect_s3_class(frequency_plot, "ggplot")
    freq_data <- frequency_plot$data
    freq_data$pbf_name <- as.character(freq_data$pbf_name)
    freq_data <- freq_data[order(freq_data$pbf_name, freq_data$valid_counts), , drop = FALSE]
    expected <- data.frame(
        valid_counts = rep(c(2L, 3L), 2L),
        count = rep(c(2L, 1L), 2L),
        pbf_name = rep(c(toy_assay, toy_assay_alt), each = 2L),
        stringsAsFactors = FALSE
    )
    expect_equal(freq_data, expected, ignore_attr = TRUE)
})
