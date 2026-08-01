test_that("pvca_plot", {
    pb_test_load_example_data()

    matrix_test <- example_proteome_matrix[1:150, ]
    pb_test_expect_warnings(
        pvca <- plot_PVCA(
            matrix_test,
            example_sample_annotation,
            technical_factors = c("MS_batch", "digestion_batch"),
            biological_factors = c("Diet", "Sex", "Strain")
        ),
        c(
            "PVCA cannot operate with missing values in the matrix",
            "filling missing values with -1"
        ),
        fixed = TRUE
    )

    expect_equal(
        pvca$data$weights[1],
        0.39166175,
        tolerance = 3e-2,
        ignore_attr = TRUE
    )
    expect_equal(
        as.character(pvca$data$label[3]),
        "MS_batch",
        ignore_attr = TRUE
    )
    expect_equal(
        as.character(pvca$data$label[2]),
        "Sex:Strain",
        ignore_attr = TRUE
    )

    expect_equal(pvca$data$category[1], "biological")
    expect_equal(pvca$data$category[3], "technical")
})

test_that("calculate_PVCA rejects factors without sampled variation", {
    data_matrix <- matrix(
        seq_len(8),
        nrow = 2,
        dimnames = list(c("feature_1", "feature_2"), paste0("sample_", 1:4))
    )
    sample_annotation <- data.frame(
        FullRunName = colnames(data_matrix),
        constant_factor = rep("one_level", ncol(data_matrix)),
        stringsAsFactors = FALSE
    )

    expect_error(
        calculate_PVCA.default(
            data_matrix,
            sample_annotation,
            feature_id_col = NULL,
            factors_for_PVCA = "constant_factor"
        ),
        "^No PVCA factors with more than one sampled level\\.$"
    )
    expect_error(
        calculate_PVCA.default(
            data_matrix,
            sample_annotation["FullRunName"],
            feature_id_col = NULL,
            factors_for_PVCA = NULL
        ),
        "^No PVCA factors with more than one sampled level\\.$"
    )
})

test_that("calculate_PVCA ProBatchFeatures routes assays and preserves return shape", {
    pbf <- pb_test_make_pbf(n_rows = 8, n_cols = 4, add_log2 = TRUE)
    selected_assays <- rev(names(pbf))
    default_annotation <- as.data.frame(SummarizedExperiment::colData(pbf))

    annotations <- setNames(
        lapply(rev(selected_assays), function(assay_name) {
            annotation <- default_annotation
            annotation$route_marker <- assay_name
            annotation
        }),
        rev(selected_assays)
    )

    calls <- list()
    testthat::local_mocked_bindings(
        calculate_PVCA.default = function(
            data_matrix,
            sample_annotation,
            feature_id_col,
            sample_id_col,
            ...
        ) {
            calls[[length(calls) + 1L]] <<- list(
                data_matrix = data_matrix,
                sample_annotation = sample_annotation,
                dots = list(...)
            )
            data.frame(
                weights = nrow(data_matrix),
                label = sample_annotation$route_marker[[1]],
                stringsAsFactors = FALSE
            )
        },
        .package = "proBatch"
    )

    multi <- suppressMessages(calculate_PVCA(
        pbf,
        pbf_name = selected_assays,
        sample_annotation = annotations,
        pca_threshold = 0.75
    ))

    expect_type(multi, "list")
    expect_named(multi, selected_assays)
    expect_equal(
        unname(vapply(multi, function(result) result$label[[1]], character(1))),
        selected_assays
    )
    expect_length(calls, length(selected_assays))
    expect_equal(
        vapply(calls, function(call) call$dots$pca_threshold, numeric(1)),
        rep(0.75, length(selected_assays))
    )
    expect_equal(
        vapply(
            calls,
            function(call) call$sample_annotation$route_marker[[1]],
            character(1)
        ),
        selected_assays
    )
    for (i in seq_along(selected_assays)) {
        expect_equal(
            calls[[i]]$data_matrix,
            suppressMessages(pb_assay_matrix(pbf, selected_assays[[i]]))
        )
    }

    calls <- list()
    single_assay <- selected_assays[[1]]
    single <- suppressMessages(calculate_PVCA(
        pbf,
        pbf_name = single_assay,
        sample_annotation = annotations[[single_assay]],
        pca_threshold = 0.5
    ))

    expect_s3_class(single, "data.frame")
    expect_identical(single$label[[1]], single_assay)
    expect_length(calls, 1L)
    expect_identical(calls[[1]]$dots$pca_threshold, 0.5)
})

test_that("prepare_PVCA_df normalizes residual labels and optionally writes CSV", {
    calculated <- data.frame(
        weights = c(0.20, 0.15, 0.10, 0.05, 0.03, 0.18, 0.17, 0.12),
        label = c(
            "Residuals",
            "Residual",
            "resid",
            "RESIDUAL",
            "bElOw 1%",
            "batch",
            "biology",
            "batch:biology"
        ),
        stringsAsFactors = FALSE
    )
    testthat::local_mocked_bindings(
        calculate_PVCA = function(...) calculated,
        .package = "proBatch"
    )

    output_dir <- tempfile("prepared_pvca_")
    on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

    prepared <- prepare_PVCA_df.default(
        data_matrix = matrix(1, nrow = 1, ncol = 1),
        sample_annotation = data.frame(FullRunName = "sample_1"),
        feature_id_col = NULL,
        technical_factors = "batch",
        biological_factors = "biology",
        variance_threshold = 0.01,
        path_to_save_results = output_dir
    )

    residual_labels <- c(
        "Residuals",
        "Residual",
        "resid",
        "RESIDUAL",
        "bElOw 1%"
    )
    expect_true(all(
        prepared$category[prepared$label %in% residual_labels] == "residual"
    ))
    expect_identical(prepared$category[prepared$label == "batch"], "technical")
    expect_identical(
        prepared$category[prepared$label == "biology"],
        "biological"
    )
    expect_identical(
        prepared$category[prepared$label == "batch:biology"],
        "biol:techn"
    )

    output_file <- file.path(output_dir, "PVCA_results_aggregated.csv")
    expect_true(file.exists(output_file))
    saved <- utils::read.csv(output_file, stringsAsFactors = FALSE)
    expected <- as.data.frame(prepared, stringsAsFactors = FALSE)
    rownames(expected) <- NULL
    expect_equal(saved, expected)
})

test_that("plot_PVCA.df adds values to prepared data", {
    prepared <- data.frame(
        weights = c(0.2, 0.3, 0.5),
        label = c("Residual", "batch", "biology"),
        category = c("residual", "technical", "biological"),
        stringsAsFactors = FALSE
    )

    pvca_plot <- plot_PVCA.df(prepared, add_values = TRUE)

    expect_s3_class(pvca_plot, "ggplot")
    text_layers <- vapply(
        pvca_plot$layers,
        function(layer) inherits(layer$geom, "GeomText"),
        logical(1)
    )
    expect_true(any(text_layers))

    built <- ggplot2::ggplot_build(pvca_plot)
    expect_setequal(
        built$data[[which(text_layers)[1]]]$label,
        sprintf("%.2f", prepared$weights)
    )
})

test_that("plot_PVCA_stacked_from_saved rebuilds a core stacked summary", {
    pvca_dir <- tempfile("saved_pvca_")
    assay_a <- file.path(pvca_dir, "assay_a")
    assay_b <- file.path(pvca_dir, "assay_b")
    dir.create(assay_a, recursive = TRUE)
    dir.create(assay_b, recursive = TRUE)
    on.exit(unlink(pvca_dir, recursive = TRUE), add = TRUE)

    utils::write.csv(
        data.frame(
            label = c("Diet", "resid"),
            weights = c(0.7, 0.3),
            category = c("biological", "residual")
        ),
        file.path(assay_a, "PVCA_results_aggregated.csv"),
        row.names = FALSE
    )
    utils::write.csv(
        data.frame(
            label = c("Diet", "resid"),
            weights = c(0.2, 0.8),
            category = c("biological", "residual")
        ),
        file.path(assay_b, "PVCA_results_aggregated.csv"),
        row.names = FALSE
    )

    stacked <- plot_PVCA_stacked_from_saved(
        pvca_dir,
        sort_stacked = "Diet",
        stacked_plot_title = c("Saved", "PVCA"),
        category_order = c("biological", "residual")
    )

    expect_s3_class(stacked, "ggplot")
    expect_equal(stacked$labels$title, "Saved\nPVCA")
    expect_setequal(
        as.character(unique(stacked$data$assay)),
        c("assay_a", "assay_b")
    )
    expect_equal(levels(stacked$data$assay), c("assay_b", "assay_a"))
})
