
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

test_that("prepare_PVCA_df aggregates weights, classifies categories, and saves CSV", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

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

test_that("plot_PVCA ProBatchFeatures stacked bar sorts assays by requested factor", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

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

test_that("prepare_variance_partition_df exports CSV and keeps categories", {
    skip_if_not_installed("variancePartition")

    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

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

test_that("plot_variance_partition.df adds medians and respects legend toggle", {
    skip_if_not_installed("variancePartition")

    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

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

test_that("plot_variance_partition ProBatchFeatures arranges multiple assays", {
    skip_if_not_installed("variancePartition")
    skip_if_not_installed("gridExtra")

    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

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
