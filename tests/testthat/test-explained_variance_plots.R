test_that("pvca_plot", {
    pb_test_load_example_data()

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

test_that("prepare_PVCA_df aggregates weights, classifies categories, and saves CSV", {
    pb_test_load_example_data()

    matrix_test <- na.omit(example_proteome_matrix)[1:80, ]
    tmp_dir <- tempfile("pvca_export")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

    pvca_res <- prepare_PVCA_df(
        matrix_test,
        example_sample_annotation,
        technical_factors = c("MS_batch", "digestion_batch"),
        biological_factors = c("Diet", "Sex", "Strain"),
        fill_the_missing = NULL,
        variance_threshold = 0.05,
        path_to_save_results = tmp_dir
    )

    expect_s3_class(pvca_res, "data.frame")
    expect_true(file.exists(file.path(tmp_dir, "PVCA_results_aggregated.csv")))
    expect_true(all(c("weights", "label", "category") %in% names(pvca_res)))
    expect_equal(tail(as.character(pvca_res$label), 1), "resid")
    expect_true(all(as.character(pvca_res$category) %in% c("biological", "technical", "biol:techn", "residual")))

    high_thr_res <- prepare_PVCA_df(
        matrix_test,
        example_sample_annotation,
        technical_factors = c("MS_batch", "digestion_batch"),
        biological_factors = c("Diet", "Sex", "Strain"),
        fill_the_missing = NULL,
        variance_threshold = 0.9
    )
    below_label <- "Below 90%"
    expect_true(below_label %in% as.character(high_thr_res$label))
    expect_equal(
        unique(as.character(high_thr_res$category[as.character(high_thr_res$label) == below_label])),
        "residual"
    )
})

test_that("prepare_PVCA_df categorizes residual labels consistently", {
    data_matrix <- matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        dimnames = list(paste0("feat", 1:2), paste0("Run", 1:2))
    )
    sample_annotation <- data.frame(
        FullRunName = colnames(data_matrix),
        stringsAsFactors = FALSE
    )

    pvca_stub <- function(...) {
        data.frame(
            feature_id = c("feat1", "feat1", "feat2", "feat2", "feat2"),
            label = c("Tech", "Bio", "Residuals", "Residual", "Below 5%"),
            weights = c(0.3, 0.25, 0.2, 0.1, 0.05),
            stringsAsFactors = FALSE
        )
    }

    vp_res <- with_mocked_bindings(
        prepare_PVCA_df.default(
            data_matrix,
            sample_annotation,
            technical_factors = "Tech",
            biological_factors = "Bio",
            fill_the_missing = NULL,
            variance_threshold = 0.05
        ),
        calculate_PVCA = pvca_stub
    )

    expect_true(all(vp_res$category[vp_res$label %in% c("Residuals", "Residual", "Below 5%")] == "residual"))
    expect_equal(vp_res$category[vp_res$label == "Tech"], "technical")
    expect_equal(vp_res$category[vp_res$label == "Bio"], "biological")
})

test_that("plot_PVCA legend follows axis order with stable colors", {
    vp_res <- data.frame(
        feature_id = rep("feat", 4),
        label = factor(
            c("TechVar", "BioVar", "Residuals", "Other"),
            levels = c("TechVar", "BioVar", "Residuals", "Other")
        ),
        weights = c(0.4, 0.35, 0.15, 0.1),
        category = c("technical", "biological", "residual", "biol:techn"),
        stringsAsFactors = FALSE
    )
    custom_colors <- c("#111111", "#222222", "#333333", "#444444")
    vp_plot <- plot_PVCA.df.default(vp_res, colors_for_bars = custom_colors)

    fill_scale <- vp_plot$scales$get_scales("fill")
    expect_equal(fill_scale$limits, levels(vp_plot$data$category))

    residual_color <- fill_scale$map("residual")
    expect_equal(residual_color, "#111111")
})

test_that("prepare_variance_partition_df keeps residual categories consistent", {
    sample_ann <- data.frame(FullRunName = paste0("Run", 1:2), stringsAsFactors = FALSE)
    data_matrix <- matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        dimnames = list(paste0("feat", 1:2), sample_ann$FullRunName)
    )
    stub_varpart <- function(...) {
        data.frame(
            feature_id = c("feat1", "feat1"),
            label = c("Residual", "Below 10%"),
            variance_explained = c(0.6, 0.05),
            stringsAsFactors = FALSE
        )
    }

    vp_res <- with_mocked_bindings(
        prepare_variance_partition_df.default(
            data_matrix,
            sample_ann,
            technical_factors = character(),
            biological_factors = character(),
            fill_the_missing = NULL
        ),
        calculate_variance_partition = stub_varpart
    )

    expect_equal(unique(as.character(vp_res$category)), "residual")
})

test_that("plot_PVCA ProBatchFeatures stacked bar sorts assays by requested factor", {
    pb_test_load_example_data()

    matrix_small <- na.omit(example_proteome_matrix)[1:40, 1:6]
    sample_ids <- colnames(matrix_small)
    sample_ann <- example_sample_annotation[match(sample_ids, example_sample_annotation$FullRunName), ]

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
    assays <- names(pbf)
    sort_label <- "Diet"

    stacked_plot <- suppressWarnings(plot_PVCA(
        pbf,
        sample_id_col = "FullRunName",
        technical_factors = c("MS_batch"),
        biological_factors = c("Diet", "Sex"),
        stacked_bar = TRUE,
        sort_stacked = sort_label,
        fill_the_missing = NULL
    ))

    expect_s3_class(stacked_plot, "ggplot")
    expect_setequal(as.character(unique(stacked_plot$data$assay)), assays)

    pvca_list <- prepare_PVCA_df(
        pbf,
        sample_id_col = "FullRunName",
        technical_factors = c("MS_batch"),
        biological_factors = c("Diet", "Sex"),
        fill_the_missing = NULL
    )
    expect_equal(length(pvca_list), length(assays))

    label_weights <- vapply(names(pvca_list), function(nm) {
        df <- pvca_list[[nm]]
        vals <- df$weights[df$label == sort_label]
        if (length(vals) == 0) {
            return(NA_real_)
        }
        sum(vals)
    }, numeric(1))
    sorted_assays <- names(sort(label_weights[!is.na(label_weights)], decreasing = TRUE))
    expected_order <- c(sorted_assays, setdiff(assays, sorted_assays))
    expect_equal(levels(stacked_plot$data$assay), rev(expected_order))
})

test_that("plot_PVCA stacked bar honors stacked_plot_title", {
    pb_test_load_example_data()

    matrix_small <- na.omit(example_proteome_matrix)[1:40, 1:6]
    sample_ids <- colnames(matrix_small)
    sample_ann <- example_sample_annotation[match(sample_ids, example_sample_annotation$FullRunName), ]

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

    stacked_title <- c("Stacked PVCA", "/path/to/data")
    stacked_plot <- suppressWarnings(plot_PVCA(
        pbf,
        sample_id_col = "FullRunName",
        technical_factors = c("MS_batch"),
        biological_factors = c("Diet", "Sex"),
        stacked_bar = TRUE,
        stacked_plot_title = stacked_title,
        fill_the_missing = NULL
    ))

    expect_equal(stacked_plot$labels$title, paste(stacked_title, collapse = "\n"))
})

test_that("prepare_variance_partition_df exports CSV and keeps categories", {
    skip_if_not_installed("variancePartition")

    pb_test_load_example_data()

    matrix_small <- na.omit(example_proteome_matrix)[1:30, 1:8]
    sample_ids <- colnames(matrix_small)
    sample_ann <- example_sample_annotation[match(sample_ids, example_sample_annotation$FullRunName), ]

    tmp_dir <- tempfile("variance_partition")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

    vp_res <- prepare_variance_partition_df(
        matrix_small,
        sample_ann,
        technical_factors = c("MS_batch"),
        biological_factors = c("Diet", "Sex"),
        fill_the_missing = NULL,
        variance_threshold = 0.3,
        path_to_save_results = tmp_dir
    )

    expect_s3_class(vp_res, "data.frame")
    expect_true(file.exists(file.path(tmp_dir, "variance_partition_results_long.csv")))
    expect_true(all(c("feature_id", "label", "variance_explained", "category") %in% names(vp_res)))
    expect_true(all(as.character(vp_res$category) %in% c("biological", "technical", "biol:techn", "residual")))

    below_label <- "Below 30%"
    expect_true(below_label %in% as.character(vp_res$label))
    expect_equal(
        unique(as.character(vp_res$category[as.character(vp_res$label) == below_label])),
        "residual"
    )
})

test_that("prepare_variance_partition_df can filter by abundance ranking", {
    sample_ann <- data.frame(FullRunName = paste0("Run", 1:3), stringsAsFactors = FALSE)
    data_matrix <- matrix(
        c(
            0, 0, 100,
            20, 20, 20,
            5, 5, 5
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(c("feat1", "feat2", "feat3"), sample_ann$FullRunName)
    )

    stub_varpart <- function(...) {
        data.frame(
            feature_id = rep(c("feat1", "feat2", "feat3"), each = 2),
            label = rep(c("Tech", "Bio"), times = 3),
            variance_explained = rep(c(0.6, 0.3), times = 3),
            stringsAsFactors = FALSE
        )
    }

    vp_mean <- with_mocked_bindings(
        prepare_variance_partition_df.default(
            data_matrix,
            sample_ann,
            technical_factors = "Tech",
            biological_factors = "Bio",
            fill_the_missing = NULL,
            abundance_top_n = 1,
            abundance_direction = "most",
            abundance_stat = "mean"
        ),
        calculate_variance_partition = stub_varpart
    )
    expect_equal(unique(vp_mean$feature_id), "feat1")

    vp_median <- with_mocked_bindings(
        prepare_variance_partition_df.default(
            data_matrix,
            sample_ann,
            technical_factors = "Tech",
            biological_factors = "Bio",
            fill_the_missing = NULL,
            abundance_top_n = 1,
            abundance_direction = "most",
            abundance_stat = "median"
        ),
        calculate_variance_partition = stub_varpart
    )
    expect_equal(unique(vp_median$feature_id), "feat2")
})

test_that("prepare_variance_partition_df validates abundance_top_n", {
    sample_ann <- data.frame(FullRunName = paste0("Run", 1:2), stringsAsFactors = FALSE)
    data_matrix <- matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        dimnames = list(c("feat1", "feat2"), sample_ann$FullRunName)
    )
    stub_varpart <- function(...) {
        data.frame(
            feature_id = rep(c("feat1", "feat2"), each = 2),
            label = rep(c("Tech", "Bio"), times = 2),
            variance_explained = rep(c(0.6, 0.3), times = 2),
            stringsAsFactors = FALSE
        )
    }
    invalid_vals <- list(0, -1, 1.5, Inf, NA_real_)
    for (val in invalid_vals) {
        expect_error(
            with_mocked_bindings(
                prepare_variance_partition_df.default(
                    data_matrix,
                    sample_ann,
                    technical_factors = "Tech",
                    biological_factors = "Bio",
                    fill_the_missing = NULL,
                    abundance_top_n = val
                ),
                calculate_variance_partition = stub_varpart
            ),
            "abundance_top_n"
        )
    }
})

test_that("plot_variance_partition.df adds medians and respects legend toggle", {
    skip_if_not_installed("variancePartition")

    pb_test_load_example_data()

    matrix_small <- na.omit(example_proteome_matrix)[1:30, 1:8]
    sample_ids <- colnames(matrix_small)
    sample_ann <- example_sample_annotation[match(sample_ids, example_sample_annotation$FullRunName), ]

    vp_res <- prepare_variance_partition_df(
        matrix_small,
        sample_ann,
        technical_factors = c("MS_batch"),
        biological_factors = c("Diet", "Sex"),
        fill_the_missing = NULL,
        variance_threshold = 0.1
    )

    vp_plot <- plot_variance_partition.df(
        vp_res,
        add_medians = TRUE,
        show_legend = FALSE,
        y_limits = c(0, 0.6)
    )

    expect_s3_class(vp_plot, "ggplot")
    expect_equal(vp_plot$theme$legend.position, "none")

    built <- ggplot2::ggplot_build(vp_plot)
    expect_equal(length(built$data), 2)
    expect_equal(vp_plot$coordinates$ylim, c(0, 0.6))
})

test_that("plot_variance_partition.df supports mean/median summaries", {
    vp_res <- data.frame(
        feature_id = rep(c("feat1", "feat2"), each = 2),
        label = rep(c("Tech", "Bio"), times = 2),
        variance_explained = c(0.4, 0.3, 0.5, 0.2),
        category = rep(c("technical", "biological"), times = 2),
        stringsAsFactors = FALSE
    )

    mean_plot <- plot_variance_partition.df.default(
        vp_res,
        summary_stat = "mean",
        add_medians = TRUE
    )
    expect_s3_class(mean_plot, "ggplot")
    mean_geoms <- vapply(mean_plot$layers, function(layer) class(layer$geom)[1], character(1))
    expect_equal(unname(mean_geoms), "GeomCol")
    expect_equal(mean_plot$labels$y, "Mean variance explained")

    median_plot <- plot_variance_partition.df.default(
        vp_res,
        summary_stat = "median",
        add_medians = TRUE
    )
    median_geoms <- vapply(median_plot$layers, function(layer) class(layer$geom)[1], character(1))
    expect_equal(unname(median_geoms), "GeomCol")
    expect_equal(median_plot$labels$y, "Median variance explained")
})

test_that("plot_variance_partition legend follows axis order with stable colors", {
    vp_res <- data.frame(
        feature_id = rep("feat", 4),
        label = factor(
            c("TechVar", "BioVar", "Residuals", "Below 5%"),
            levels = c("TechVar", "BioVar", "Residuals", "Below 5%")
        ),
        variance_explained = c(0.4, 0.35, 0.2, 0.1),
        category = c("technical", "biological", "residual", "residual"),
        stringsAsFactors = FALSE
    )
    custom_colors <- c("#111111", "#222222", "#333333", "#444444")
    vp_plot <- plot_variance_partition.df.default(vp_res, colors_for_boxes = custom_colors)

    fill_scale <- vp_plot$scales$get_scales("fill")
    expect_equal(fill_scale$limits, c("technical", "biological", "residual"))

    residual_color <- fill_scale$map("residual")
    expect_equal(residual_color, "#111111")
})

test_that("plot_variance_partition ProBatchFeatures arranges multiple assays", {
    skip_if_not_installed("variancePartition")
    skip_if_not_installed("gridExtra")

    pb_test_load_example_data()

    matrix_small <- na.omit(example_proteome_matrix)[1:35, 1:6]
    sample_ids <- colnames(matrix_small)
    sample_ann <- example_sample_annotation[match(sample_ids, example_sample_annotation$FullRunName), ]

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

    res <- suppressWarnings(plot_variance_partition(
        pbf,
        sample_id_col = "FullRunName",
        technical_factors = c("MS_batch"),
        biological_factors = c("Diet", "Sex"),
        fill_the_missing = NULL,
        variance_threshold = 0.1,
        add_medians = TRUE,
        show_legend = FALSE,
        return_gridExtra = TRUE
    ))

    expect_type(res, "list")
    expect_true(all(c("grob", "plots") %in% names(res)))
    expect_equal(length(res$plots), length(names(pbf)))
    expect_true(all(vapply(res$plots, inherits, logical(1), "ggplot")))
})

build_variance_partition_test_pbf <- function() {
    pb_test_load_example_data()

    matrix_small <- na.omit(example_proteome_matrix)[1:35, 1:6]
    sample_ids <- colnames(matrix_small)
    sample_ann <- example_sample_annotation[match(sample_ids, example_sample_annotation$FullRunName), ]

    pbf <- suppressMessages(ProBatchFeatures(
        data_matrix = matrix_small,
        sample_annotation = sample_ann,
        sample_id_col = "FullRunName",
        name = "feature::raw"
    ))
    suppressMessages(pb_transform(
        pbf,
        from = "feature::raw",
        steps = "log2",
        store_fast_steps = TRUE
    ))
}

test_that("plot_variance_partition shares y limits when organize_pbfs is TRUE", {
    skip_if_not_installed("variancePartition")
    skip_if_not_installed("gridExtra")

    pbf <- build_variance_partition_test_pbf()
    vp_list <- suppressWarnings(prepare_variance_partition_df(
        pbf,
        sample_id_col = "FullRunName",
        technical_factors = c("MS_batch"),
        biological_factors = c("Diet", "Sex"),
        fill_the_missing = NULL,
        variance_threshold = 0.1
    ))
    if (!is.list(vp_list)) {
        vp_list <- list(vp_list)
    }
    shared_range <- range(unlist(lapply(vp_list, function(df) df$variance_explained)), na.rm = TRUE)
    expect_true(all(is.finite(shared_range)))

    res <- suppressWarnings(plot_variance_partition(
        pbf,
        sample_id_col = "FullRunName",
        technical_factors = c("MS_batch"),
        biological_factors = c("Diet", "Sex"),
        fill_the_missing = NULL,
        variance_threshold = 0.1,
        organize_pbfs = TRUE,
        return_gridExtra = TRUE
    ))

    ylims <- lapply(res$plots, function(p) p$coordinates$ylim)
    expect_true(all(vapply(ylims, function(x) identical(x, shared_range), logical(1))))
})

test_that("explicit y limits override shared organize_pbfs limits", {
    skip_if_not_installed("variancePartition")
    skip_if_not_installed("gridExtra")

    pbf <- build_variance_partition_test_pbf()
    target_limits <- c(0, 0.5)
    res <- suppressWarnings(plot_variance_partition(
        pbf,
        sample_id_col = "FullRunName",
        technical_factors = c("MS_batch"),
        biological_factors = c("Diet", "Sex"),
        fill_the_missing = NULL,
        variance_threshold = 0.1,
        organize_pbfs = TRUE,
        y_limits = target_limits,
        return_gridExtra = TRUE
    ))

    ylims <- lapply(res$plots, function(p) p$coordinates$ylim)
    expect_true(all(vapply(ylims, function(x) identical(x, target_limits), logical(1))))
})
