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

test_that(".pb_group_missing_matrix binarizes grouped fractions at 0.5 when requested", {
    grouped <- .pb_group_missing_matrix(
        data_matrix = toy_matrix,
        sample_annotation = toy_sa,
        sample_id_col = "FullRunName",
        color_by = "Condition",
        drop_complete = FALSE,
        force_binarization = TRUE
    )

    expected <- matrix(
        c(
            1, 1, 0,
            1, 1, 1
        ),
        nrow = 3,
        dimnames = list(
            rownames(toy_matrix),
            c("Condition=A", "Condition=B")
        )
    )

    expect_equal(grouped$matrix, expected)
})

test_that(".pb_group_missing_matrix applies numeric binarization before drop_complete", {
    grouped <- .pb_group_missing_matrix(
        data_matrix = toy_matrix,
        sample_annotation = toy_sa,
        sample_id_col = "FullRunName",
        color_by = "Condition",
        drop_complete = TRUE,
        force_binarization = 0.75
    )

    expected <- matrix(
        c(
            1, 0,
            0, 1
        ),
        nrow = 2,
        dimnames = list(
            c("prot1", "prot3"),
            c("Condition=A", "Condition=B")
        )
    )

    expect_equal(grouped$matrix, expected)
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

test_that(".pb_group_missing_matrix preserves one-feature matrix dimensions", {
    skip_if_not_installed("pheatmap")

    one_feature <- toy_matrix[1L, , drop = FALSE]
    grouped <- .pb_group_missing_matrix(
        data_matrix = one_feature,
        sample_annotation = toy_sa,
        sample_id_col = "FullRunName",
        color_by = "Condition",
        drop_complete = FALSE
    )

    expected <- matrix(
        c(1, 0.5),
        nrow = 1L,
        dimnames = list(
            rownames(one_feature),
            c("Condition=A", "Condition=B")
        )
    )
    expect_equal(grouped$matrix, expected)

    result <- plot_grouped_NA_heatmap(
        one_feature,
        sample_annotation = toy_sa,
        color_by = "Condition",
        drop_complete = FALSE,
        draw = FALSE
    )
    expect_s3_class(result, "pheatmap")
})

test_that("plot_grouped_NA_heatmap supports a single group with default clustering", {
    skip_if_not_installed("pheatmap")

    one_group_annotation <- toy_sa
    one_group_annotation$Condition <- "only"

    result <- plot_grouped_NA_heatmap(
        toy_matrix,
        sample_annotation = one_group_annotation,
        color_by = "Condition",
        drop_complete = FALSE,
        draw = FALSE
    )

    expect_s3_class(result, "pheatmap")
})

test_that(".pb_group_missing_matrix keeps delimiter-bearing tuples separate", {
    collision_matrix <- matrix(
        c(1, NA_real_),
        nrow = 1L,
        dimnames = list("feature", c("S1", "S2"))
    )
    collision_annotation <- data.frame(
        First = c("x", "x | y"),
        Second = c("y | z", "z"),
        row.names = colnames(collision_matrix),
        stringsAsFactors = FALSE
    )

    grouped <- .pb_group_missing_matrix(
        data_matrix = collision_matrix,
        sample_annotation = collision_annotation,
        sample_id_col = NULL,
        color_by = c("First", "Second"),
        drop_complete = FALSE
    )

    expect_equal(dim(grouped$matrix), c(1L, 2L))
    expect_equal(unname(grouped$matrix), matrix(c(1, 0), nrow = 1L))
})

test_that(".pb_group_missing_matrix distinguishes missing and literal NA labels", {
    missing_matrix <- matrix(
        c(1, NA_real_),
        nrow = 1L,
        dimnames = list("feature", c("S1", "S2"))
    )
    missing_annotation <- data.frame(
        Group = c(NA_character_, "<NA>"),
        row.names = colnames(missing_matrix),
        stringsAsFactors = FALSE
    )

    grouped <- .pb_group_missing_matrix(
        data_matrix = missing_matrix,
        sample_annotation = missing_annotation,
        sample_id_col = NULL,
        color_by = "Group",
        drop_complete = FALSE
    )

    expect_equal(ncol(grouped$matrix), 2L)
    expect_equal(sort(as.numeric(grouped$matrix)), c(0, 1))
    expect_length(unique(colnames(grouped$matrix)), 2L)
    expect_true(any(startsWith(colnames(grouped$matrix), "Group=<missing>")))
    expect_true(any(startsWith(colnames(grouped$matrix), "Group=<NA>")))
    expect_equal(sum(is.na(grouped$annotation$Group)), 1L)
    expect_equal(sum(grouped$annotation$Group == "<NA>", na.rm = TRUE), 1L)

    prepared <- .pb_prepare_sample_annotation(
        sample_annotation = grouped$annotation,
        sample_id_col = NULL,
        color_by = "Group",
        label_by = ".group_label",
        sample_order = colnames(grouped$matrix),
        col_vector = NULL,
        allow_color_disable = FALSE
    )
    expect_equal(sum(is.na(prepared$annotation_col$Group)), 1L)
    expect_named(prepared$annotation_colors$Group, "<NA>")
})

test_that(".pb_group_missing_matrix protects generated metadata columns", {
    for (reserved_name in c(".group_key", ".group_size", ".group_label")) {
        reserved_annotation <- data.frame(
            value = c("small", "large", "large"),
            row.names = colnames(toy_matrix),
            stringsAsFactors = FALSE
        )
        names(reserved_annotation) <- reserved_name

        expect_error(
            .pb_group_missing_matrix(
                data_matrix = toy_matrix,
                sample_annotation = reserved_annotation,
                sample_id_col = NULL,
                color_by = reserved_name,
                drop_complete = FALSE
            ),
            "reserved names"
        )
    }
})

test_that(".pb_group_missing_matrix accepts a grouping column named No", {
    skip_if_not_installed("pheatmap")

    no_annotation <- data.frame(
        No = c("first", "second", "second"),
        row.names = colnames(toy_matrix),
        stringsAsFactors = FALSE
    )

    grouped <- .pb_group_missing_matrix(
        data_matrix = toy_matrix,
        sample_annotation = no_annotation,
        sample_id_col = NULL,
        color_by = "No",
        drop_complete = FALSE
    )

    expect_equal(ncol(grouped$matrix), 2L)
    result <- plot_grouped_NA_heatmap(
        toy_matrix,
        sample_annotation = no_annotation,
        color_by = "No",
        cluster_samples = FALSE,
        cluster_features = FALSE,
        draw = FALSE
    )
    expect_s3_class(result, "pheatmap")

    prepared <- .pb_prepare_sample_annotation(
        sample_annotation = grouped$annotation,
        sample_id_col = NULL,
        color_by = grouped$color_by,
        label_by = ".group_label",
        sample_order = colnames(grouped$matrix),
        col_vector = NULL,
        allow_color_disable = FALSE
    )
    expect_named(prepared$annotation_col, "No")
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

test_that("plot_grouped_NA_heatmap.default validates force_binarization", {
    expect_error(
        plot_grouped_NA_heatmap(
            toy_matrix,
            sample_annotation = toy_sa,
            color_by = "Condition",
            force_binarization = 1.1,
            draw = FALSE
        ),
        "`force_binarization` must be FALSE, TRUE, or a single numeric value between 0 and 1."
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
    expect_equal(names(res$heatmaps), c(toy_assay, toy_assay_alt))
    expect_true(all(vapply(res$heatmaps, inherits, logical(1), what = "pheatmap")))
})

test_that("multiple requested heatmaps remain arranged with one surviving assay", {
    skip_if_not_installed("pheatmap")
    skip_if_not_installed("gridExtra")

    complete_matrix <- matrix(
        seq_len(length(toy_matrix)),
        nrow = nrow(toy_matrix),
        dimnames = dimnames(toy_matrix)
    )
    complete_assay <- paste0(toy_assay, "_complete")
    complete_se <- SummarizedExperiment(
        assays = list(intensity = complete_matrix),
        colData = colData(pbf_toy[[toy_assay]])
    )
    pbf_with_complete <- proBatch:::.pb_add_assay_with_link(
        pbf_toy,
        se = complete_se,
        to = complete_assay,
        from = toy_assay
    )

    plot_file <- tempfile(fileext = ".pdf")
    grDevices::pdf(plot_file)
    on.exit(
        {
            grDevices::dev.off()
            unlink(plot_file)
        },
        add = TRUE
    )

    expect_warning(
        res <- plot_grouped_NA_heatmap(
            pbf_with_complete,
            pbf_name = c(toy_assay, complete_assay),
            color_by = "Condition",
            cluster_samples = FALSE,
            cluster_features = FALSE,
            drop_complete = TRUE,
            draw = TRUE
        ),
        "no rows with missing values|Skipping assay"
    )

    expect_type(res, "list")
    expect_s3_class(res$grob, "gtable")
    expect_length(res$heatmaps, 1L)
    expect_named(res$heatmaps, toy_assay)
    expect_s3_class(res$heatmaps[[1L]], "pheatmap")
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

test_that("plot_NA_density.default supports grouped densities", {
    density_plot <- plot_NA_density(
        toy_matrix,
        sample_annotation = toy_sa,
        sample_id_col = "FullRunName",
        color_by = "Condition"
    )

    expect_s3_class(density_plot, "ggplot")
    grouped_data <- density_plot$data[
        ,
        c("mean", "Type", ".pb_density_group"),
        drop = FALSE
    ]
    grouped_data <- grouped_data[
        order(grouped_data$.pb_density_group, grouped_data$mean),
        ,
        drop = FALSE
    ]
    expected <- data.frame(
        mean = c(7, 10, 5.5, 8.5, 12),
        Type = c(
            "Valid Value",
            "Valid Value",
            "Valid Value",
            "Valid Value",
            "Missing Value"
        ),
        .pb_density_group = c(
            rep("A", 2L),
            rep("B", 3L)
        ),
        stringsAsFactors = FALSE
    )

    expect_equal(grouped_data, expected, ignore_attr = TRUE)
    expect_equal(density_plot$labels$colour, "Condition")
    expect_equal(density_plot$labels$linetype, "Value Type")
})

test_that("plot_NA_density.default keeps colliding group keys distinct", {
    collision_annotation <- data.frame(
        Group = c(NA_character_, "<missing>", "other"),
        row.names = colnames(toy_matrix),
        stringsAsFactors = FALSE
    )
    density_data <- .pb_missing_density_df(
        data_matrix = toy_matrix,
        assay_nm = "",
        missing_label = "Missing Value",
        valid_label = "Valid Value",
        sample_annotation = collision_annotation,
        color_by = "Group"
    )
    group_keys <- unique(density_data$.pb_density_group)

    expect_false(anyNA(group_keys))
    expect_setequal(
        group_keys,
        c("<missing>", "<missing> #1", "other")
    )

    collision_scheme <- setNames(
        c("#111111", "#222222", "#333333"),
        c("<missing>", "<missing> #1", "other")
    )
    density_plot <- plot_NA_density(
        toy_matrix,
        sample_annotation = collision_annotation,
        color_by = "Group",
        color_scheme = collision_scheme
    )
    built <- suppressWarnings(ggplot_build(density_plot))

    expect_setequal(
        unique(built$data[[1L]]$colour),
        unname(collision_scheme)
    )
})

test_that("plot_NA_density.default forwards density parameters", {
    density_plot <- plot_NA_density(
        toy_matrix,
        adjust = 2,
        na.rm = FALSE
    )

    expect_equal(density_plot$layers[[1L]]$stat_params$adjust, 2)
    expect_false(density_plot$layers[[1L]]$stat_params$na.rm)
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

test_that("plot_NA_density.ProBatchFeatures supports grouped densities", {
    density_plot <- plot_NA_density(
        pbf_multi,
        pbf_name = c(toy_assay, toy_assay_alt),
        color_by = "Condition",
        nrow = 1,
        ncol = 2
    )

    expect_s3_class(density_plot, "ggplot")
    expect_setequal(
        unique(as.character(density_plot$data$pbf_name)),
        c(toy_assay, toy_assay_alt)
    )
    expect_setequal(
        unique(density_plot$data$.pb_density_group),
        c("A", "B")
    )
    expect_equal(density_plot$labels$colour, "Condition")
    expect_equal(density_plot$labels$linetype, "Value Type")
})

test_that("plot_NA_density.default applies grouped color schemes", {
    named_plot <- plot_NA_density(
        toy_matrix,
        sample_annotation = toy_sa,
        sample_id_col = "FullRunName",
        color_by = "Condition",
        color_scheme = c(B = "#222222", A = "#111111")
    )
    named_built <- suppressWarnings(ggplot_build(named_plot))

    expect_setequal(
        unique(named_built$data[[1L]]$colour),
        c("#111111", "#222222")
    )
    named_scale <- named_built$plot$scales$get_scales("colour")
    expect_equal(
        unname(named_scale$map(c("A", "B"))),
        c("#111111", "#222222")
    )

    brewer_plot <- plot_NA_density(
        toy_matrix,
        sample_annotation = toy_sa,
        sample_id_col = "FullRunName",
        color_by = "Condition",
        color_scheme = "brewer"
    )
    brewer_built <- suppressWarnings(ggplot_build(brewer_plot))

    expect_setequal(
        unique(brewer_built$data[[1L]]$colour),
        c("#E41A1C", "#377EB8")
    )

    override_plot <- plot_NA_density(
        toy_matrix,
        sample_annotation = toy_sa,
        sample_id_col = "FullRunName",
        color_by = "Condition",
        color_scheme = c(A = "#111111", B = "#222222"),
        col_vector = "#333333"
    )
    override_built <- suppressWarnings(ggplot_build(override_plot))

    expect_equal(unique(override_built$data[[1L]]$colour), "#333333")
})

test_that("plot_NA_density.ProBatchFeatures accepts color scheme lists", {
    density_plot <- plot_NA_density(
        pbf_multi,
        pbf_name = c(toy_assay, toy_assay_alt),
        color_by = "Condition",
        color_scheme = list(
            Condition = c(B = "#222222", A = "#111111")
        )
    )
    built <- suppressWarnings(ggplot_build(density_plot))

    expect_setequal(
        unique(built$data[[1L]]$colour),
        c("#111111", "#222222")
    )
})


test_that("plot_NA_frequency.default summarises observation counts", {
    frequency_plot <- plot_NA_frequency(toy_matrix)

    expect_s3_class(frequency_plot, "ggplot")
    freq_data <- frequency_plot$data
    freq_data <- freq_data[order(freq_data$valid_counts), , drop = FALSE]
    expected <- data.frame(
        valid_counts = c(2L, 3L),
        count = c(2L, 1L)
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
