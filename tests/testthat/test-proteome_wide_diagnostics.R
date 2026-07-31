test_that("hierarchical_clustering", {
    pb_test_load_example_data()

    matrix_test <- example_proteome_matrix[1:10, ]

    # exclude columns EarTag Strain Sex RunDate RunTime digestion_batch
    example_sample_annotation <- example_sample_annotation %>%
        select("FullRunName", "MS_batch", "Diet")

    color_list <- sample_annotation_to_colors(
        example_sample_annotation,
        sample_id_col = "FullRunName",
        factor_columns = c("MS_batch", "Diet")
    )

    pb_test_expect_warnings(
        hiearchical <- plot_hierarchical_clustering(matrix_test,
            sample_annotation = example_sample_annotation,
            factors_to_plot = c("MS_batch", "Diet"),
            distance = "euclidean", agglomeration = "complete",
            color_list = color_list,
            label_samples = FALSE
        ),
        c(
            "filling missing values with 0",
            "Hierarchical clustering cannot operate with missing values",
            "color list and sample annotation have different factors,"
        ),
        fixed = TRUE
    )

    expect_identical(hiearchical$mar[[1]], 1)
    expect_identical(hiearchical$mar[[2]], 5)
    expect_identical(hiearchical$mar[[3]], 0)
    expect_identical(hiearchical$mar[[4]], 1)
})


test_that("heatmap_plot", {
    pb_test_load_example_data()

    matrix_test <- example_proteome_matrix[1:20, ]
    example_sample_annotation <- example_sample_annotation %>%
        select("FullRunName", "MS_batch", "Sex", "digestion_batch", "Diet")


    color_list <- sample_annotation_to_colors(
        example_sample_annotation,
        sample_id_col = "FullRunName",
        factor_columns = c("MS_batch", "Sex", "digestion_batch", "Diet")
    )

    suppressWarnings(
        heatmap <- plot_heatmap_diagnostic(
            matrix_test,
            sample_annotation = example_sample_annotation,
            factors_to_plot = c("MS_batch", "Sex", "digestion_batch", "Diet"),
            cluster_cols = TRUE,
            show_rownames = TRUE, show_colnames = FALSE,
            color_list = color_list
        )
    )

    expect_equal(heatmap$tree_row$method, "complete")
    expect_equal(heatmap$tree_row$dist.method, "euclidean")

    expect_equal(heatmap$tree_row$labels[1], "10062_NVGVSFYADKPEVTQEQK_3")
    expect_equal(heatmap$tree_row$labels[2], "10063_NVGVSFYADKPEVTQEQKK_3")

    expect_equal(heatmap$gtable$layout$name[1], "col_tree")
    expect_equal(heatmap$gtable$layout$name[3], "matrix")
    expect_equal(heatmap$gtable$layout$name[5], "col_annotation")
    expect_equal(heatmap$gtable$layout$name[8], "legend")

    expect_equal(heatmap$gtable$layout$t, c(2, 4, 4, 4, 3, 3, 3, 3))
})


test_that("plot_PCA warns on missing values and honors x/y PC selection", {
    pb_test_load_example_data()

    expect_error(
        plot_PCA(
            example_proteome_matrix,
            example_sample_annotation,
            fill_the_missing = Inf
        ),
        "finite numeric scalar"
    )

    pb_test_expect_warnings(
        pca <- plot_PCA(
            data_matrix = example_proteome_matrix,
            sample_annotation = example_sample_annotation,
            color_by = "MS_batch", plot_title = "PCA colored by MS batch",
            fill_the_missing = -1
        ),
        c(
            "filling missing values with -1",
            "PCA cannot operate with missing values in the matrix",
            "color_scheme will be inferred automatically"
        ),
        fixed = TRUE
    )

    expect_equal(pca$labels$y, "PC2 (14.24%)")
    expect_equal(pca$labels$x, "PC1 (69.50%)")
    expect_equal(pca$labels$colour, "MS_batch")

    pca_xy <- suppressWarnings(plot_PCA(
        example_proteome_matrix, example_sample_annotation,
        color_by = "MS_batch",
        fill_the_missing = -1,
        x_nPC = 2,
        y_nPC = 3
    ))

    expect_match(pca_xy$labels$x, "^PC2 ")
    expect_match(pca_xy$labels$y, "^PC3 ")
})

test_that("plot_PCA builds numeric annotations with discrete or continuous colors", {
    pb_test_load_example_data()

    matrix_test <- example_proteome_matrix[1:20, 1:8]
    sample_ids <- colnames(matrix_test)
    sample_annotation <- example_sample_annotation[
        match(sample_ids, example_sample_annotation$FullRunName), ,
        drop = FALSE
    ]

    sample_annotation$numeric_group <- rep(1:2, length.out = nrow(sample_annotation))
    discrete_palette <- c("1" = "#1B9E77", "2" = "#D95F02")
    discrete_plot <- suppressWarnings(plot_PCA(
        matrix_test,
        sample_annotation,
        color_by = "numeric_group",
        color_scheme = discrete_palette,
        fill_the_missing = -1
    ))
    discrete_build <- ggplot2::ggplot_build(discrete_plot)

    expect_setequal(
        unique(discrete_build$data[[1]]$colour),
        unname(discrete_palette)
    )

    sample_annotation$numeric_gradient <- seq_len(nrow(sample_annotation))
    continuous_plot <- suppressWarnings(plot_PCA(
        matrix_test,
        sample_annotation,
        color_by = "numeric_gradient",
        color_scheme = c("#2166AC", "#B2182B"),
        fill_the_missing = -1
    ))
    continuous_build <- ggplot2::ggplot_build(continuous_plot)

    expect_gt(length(unique(continuous_build$data[[1]]$colour)), 2)
})

test_that("plot_PCA supports marginal density plots", {
    pb_test_load_example_data()

    pca <- suppressMessages(suppressWarnings(plot_PCA(
        example_proteome_matrix, example_sample_annotation,
        color_by = "MS_batch",
        fill_the_missing = -1,
        marginal_density = TRUE
    )))

    expect_true(inherits(pca, "ggplot") || grid::is.grob(pca))

    matrix_test <- stats::na.omit(example_proteome_matrix)[1:20, 1:8]
    annotation_test <- example_sample_annotation[
        match(
            colnames(matrix_test),
            example_sample_annotation$FullRunName
        ), ,
        drop = FALSE
    ]
    annotation_test$Dim1 <- factor(
        rep(c("first", "second"), length.out = nrow(annotation_test))
    )
    batch_levels <- unique(as.character(annotation_test$MS_batch))
    collision_colors <- list(
        MS_batch = stats::setNames(
            grDevices::hcl.colors(length(batch_levels), "Dark 3"),
            batch_levels
        ),
        Dim1 = c(first = "#1B9E77", second = "#D95F02")
    )
    collision_pca <- expect_no_warning(suppressMessages(plot_PCA(
        data_matrix = matrix_test,
        sample_annotation = annotation_test,
        color_by = "MS_batch",
        color_scheme = collision_colors,
        marginal_density = "Dim1"
    )))
    expect_true(inherits(collision_pca, "ggplot") ||
        grid::is.grob(collision_pca))
})

pb_test_embedding_fixture <- function(n_features = 12L, n_samples = 8L,
                                      with_missing = FALSE) {
    sample_ids <- paste0("sample_", seq_len(n_samples))
    feature_ids <- paste0("feature_", seq_len(n_features))
    values <- outer(
        seq_len(n_features),
        seq_len(n_samples),
        function(feature, sample) {
            sin(feature * sample / 3) +
                cos((feature + sample) / 4) +
                feature / 20
        }
    )
    dimnames(values) <- list(feature_ids, sample_ids)
    if (isTRUE(with_missing)) {
        values[1L, 2L] <- NA_real_
    }

    annotation <- data.frame(
        FullRunName = sample_ids,
        group = factor(rep(c("control", "treated"), length.out = n_samples)),
        score = seq_len(n_samples),
        marker = rep(c("circle", "triangle"), length.out = n_samples),
        stringsAsFactors = FALSE
    )
    annotation <- annotation[rev(seq_len(nrow(annotation))), , drop = FALSE]

    list(
        matrix = values,
        annotation = annotation,
        aligned_annotation = annotation[
            match(sample_ids, annotation$FullRunName), ,
            drop = FALSE
        ]
    )
}

pb_test_fake_embedding <- function(input, n_components = 2L) {
    components <- vapply(seq_len(n_components), function(component) {
        rowSums(input * component) + seq_len(nrow(input)) / (component + 1)
    }, numeric(nrow(input)))
    matrix(components, nrow = nrow(input), ncol = n_components)
}

test_that("matrix embeddings align annotations and forward backend arguments", {
    skip_if_not_installed("Rtsne")
    skip_if_not_installed("umap")

    fixture <- pb_test_embedding_fixture()
    captured <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        Rtsne = function(X, ...) {
            captured$tsne_input <- X
            captured$tsne_args <- list(...)
            list(Y = pb_test_fake_embedding(
                X,
                n_components = captured$tsne_args$dims
            ))
        },
        .package = "Rtsne"
    )
    testthat::local_mocked_bindings(
        umap = function(d, config, ...) {
            captured$umap_input <- d
            captured$umap_config <- config
            captured$umap_args <- list(...)
            list(layout = pb_test_fake_embedding(
                d,
                n_components = config$n_components
            ))
        },
        .package = "umap"
    )

    discrete_palette <- c(
        control = "#1B9E77",
        treated = "#D95F02"
    )
    tsne_plot <- plot_TSNE(
        data_matrix = fixture$matrix,
        sample_annotation = fixture$annotation,
        color_by = "group",
        shape_by = "marker",
        color_scheme = discrete_palette,
        perplexity = 2,
        max_iter = 250,
        random_seed = 71,
        theta = 0.17
    )
    tsne_build <- ggplot2::ggplot_build(tsne_plot)

    expect_s3_class(tsne_plot, "ggplot")
    expect_equal(captured$tsne_input, t(fixture$matrix))
    expect_identical(captured$tsne_args$theta, 0.17)
    expect_false("random_seed" %in% names(captured$tsne_args))
    expect_false("check_duplicates" %in% names(captured$tsne_args))
    expect_identical(tsne_plot$data$sample_id, colnames(fixture$matrix))
    expect_identical(
        tsne_plot$data$FullRunName,
        fixture$aligned_annotation$FullRunName
    )
    expect_identical(tsne_plot$data$score, fixture$aligned_annotation$score)
    expect_setequal(
        unique(tsne_build$data[[1L]]$colour),
        unname(discrete_palette)
    )
    expect_length(unique(tsne_build$data[[1L]]$shape), 2L)

    expect_warning(
        umap_plot <- plot_UMAP(
            data_matrix = fixture$matrix,
            sample_annotation = fixture$annotation,
            color_by = "score",
            shape_by = "marker",
            color_scheme = c("#2166AC", "#B2182B"),
            n_neighbors = 3,
            min_dist = 0.2,
            metric = "manhattan",
            spread = 1.4,
            learning_rate = 0.42,
            random_state = 19,
            preserve.seed = FALSE
        ),
        "very few values"
    )
    umap_build <- ggplot2::ggplot_build(umap_plot)

    expect_s3_class(umap_plot, "ggplot")
    expect_equal(captured$umap_input, t(fixture$matrix))
    expect_equal(captured$umap_config$n_neighbors, 3)
    expect_identical(captured$umap_config$min_dist, 0.2)
    expect_identical(captured$umap_config$metric, "manhattan")
    expect_identical(captured$umap_config$spread, 1.4)
    expect_identical(captured$umap_config$alpha, 0.42)
    expect_identical(captured$umap_args$preserve.seed, FALSE)
    expect_identical(umap_plot$data$sample_id, colnames(fixture$matrix))
    expect_identical(umap_plot$data$score, fixture$aligned_annotation$score)
    expect_gt(length(unique(umap_build$data[[1L]]$colour)), 2L)
    expect_length(unique(umap_build$data[[1L]]$shape), 2L)

    sample_name_collision_matrix <- fixture$matrix
    colnames(sample_name_collision_matrix)[1L] <- ".pb_feature_id"
    sample_name_collision_annotation <- fixture$annotation
    sample_name_collision_annotation$FullRunName[
        sample_name_collision_annotation$FullRunName == "sample_1"
    ] <- ".pb_feature_id"
    expect_no_warning(
        plot_TSNE(
            sample_name_collision_matrix,
            sample_name_collision_annotation,
            color_by = "group",
            color_scheme = discrete_palette,
            perplexity = 2,
            max_iter = 250,
            random_seed = 71
        )
    )
    expect_equal(captured$tsne_input, t(sample_name_collision_matrix))

    collision_annotation <- fixture$annotation
    names(collision_annotation)[
        names(collision_annotation) == "FullRunName"
    ] <- ".pb_feature_id"
    expect_no_warning(
        plot_TSNE(
            fixture$matrix,
            collision_annotation,
            sample_id_col = ".pb_feature_id",
            color_by = "group",
            color_scheme = discrete_palette,
            perplexity = 2,
            max_iter = 250,
            random_seed = 71
        )
    )
    expect_equal(captured$tsne_input, t(fixture$matrix))

    reserved_annotation <- fixture$annotation
    reserved_annotation$Dim1 <- reserved_annotation$group
    reserved_annotation$sample_id <- reserved_annotation$marker
    reserved_aligned <- reserved_annotation[
        match(
            colnames(fixture$matrix),
            reserved_annotation$FullRunName
        ), ,
        drop = FALSE
    ]
    reserved_plot <- expect_no_warning(plot_TSNE(
        data_matrix = fixture$matrix,
        sample_annotation = reserved_annotation,
        color_by = "Dim1",
        shape_by = "sample_id",
        color_scheme = discrete_palette,
        perplexity = 2,
        max_iter = 250,
        random_seed = 71
    ))
    expect_identical(
        as.character(reserved_plot$data$Dim1),
        as.character(reserved_aligned$Dim1)
    )
    expect_identical(
        as.character(reserved_plot$data$sample_id),
        as.character(reserved_aligned$sample_id)
    )
})

test_that("embedding defaults handle missing values consistently", {
    skip_if_not_installed("Rtsne")
    skip_if_not_installed("umap")

    expect_identical(eval(formals(plot_TSNE.default)$fill_the_missing), -1)
    expect_identical(eval(formals(plot_UMAP.default)$fill_the_missing), -1)

    fixture <- pb_test_embedding_fixture(with_missing = TRUE)
    captured <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        Rtsne = function(X, ...) {
            captured$tsne_inputs <- c(captured$tsne_inputs, list(X))
            list(Y = pb_test_fake_embedding(X))
        },
        .package = "Rtsne"
    )
    testthat::local_mocked_bindings(
        umap = function(d, config, ...) {
            captured$umap_inputs <- c(captured$umap_inputs, list(d))
            list(layout = pb_test_fake_embedding(d, config$n_components))
        },
        .package = "umap"
    )

    pb_test_expect_warnings(
        plot_TSNE(
            fixture$matrix,
            fixture$annotation,
            color_by = "group",
            color_scheme = c(control = "#1B9E77", treated = "#D95F02"),
            perplexity = 2,
            max_iter = 250,
            random_seed = 3
        ),
        c(
            "t-SNE cannot operate with missing values",
            "filling missing values with -1"
        ),
        fixed = TRUE
    )
    expect_false(anyNA(captured$tsne_inputs[[1L]]))
    expect_true(-1 %in% captured$tsne_inputs[[1L]])

    pb_test_expect_warnings(
        plot_TSNE(
            fixture$matrix,
            fixture$annotation,
            color_by = "group",
            color_scheme = c(control = "#1B9E77", treated = "#D95F02"),
            fill_the_missing = NULL,
            perplexity = 2,
            max_iter = 250,
            random_seed = 3
        ),
        c(
            "t-SNE cannot operate with missing values",
            "removed 1 rows and 0 columns"
        ),
        fixed = TRUE
    )
    expect_identical(
        ncol(captured$tsne_inputs[[2L]]),
        nrow(fixture$matrix) - 1L
    )

    pb_test_expect_warnings(
        plot_UMAP(
            fixture$matrix,
            fixture$annotation,
            color_by = "group",
            color_scheme = c(control = "#1B9E77", treated = "#D95F02"),
            n_neighbors = 3,
            random_state = 3
        ),
        c(
            "UMAP cannot operate with missing values",
            "filling missing values with -1"
        ),
        fixed = TRUE
    )
    expect_false(anyNA(captured$umap_inputs[[1L]]))
    expect_true(-1 %in% captured$umap_inputs[[1L]])

    pb_test_expect_warnings(
        plot_UMAP(
            fixture$matrix,
            fixture$annotation,
            color_by = "group",
            color_scheme = c(control = "#1B9E77", treated = "#D95F02"),
            fill_the_missing = FALSE,
            n_neighbors = 3,
            random_state = 3
        ),
        c(
            "UMAP cannot operate with missing values",
            "removed 1 rows and 0 columns"
        ),
        fixed = TRUE
    )
    expect_identical(
        ncol(captured$umap_inputs[[2L]]),
        nrow(fixture$matrix) - 1L
    )

    all_incomplete <- fixture$matrix
    all_incomplete[, 1L] <- NA_real_
    warnings_seen <- character()
    expect_error(
        withCallingHandlers(
            plot_TSNE(
                all_incomplete,
                fixture$annotation,
                color_by = "group",
                color_scheme = c(
                    control = "#1B9E77",
                    treated = "#D95F02"
                ),
                fill_the_missing = NULL,
                perplexity = 2,
                max_iter = 250,
                random_seed = 3
            ),
            warning = function(condition) {
                warnings_seen <<- c(
                    warnings_seen,
                    conditionMessage(condition)
                )
                invokeRestart("muffleWarning")
            }
        ),
        "discarded every feature"
    )
    expect_identical(
        warnings_seen,
        "t-SNE cannot operate with missing values in the matrix"
    )
})

test_that("embedding validation enforces dimensions, samples, and perplexity", {
    skip_if_not_installed("Rtsne")
    skip_if_not_installed("umap")

    fixture <- pb_test_embedding_fixture()
    captured <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        Rtsne = function(X, ...) {
            args <- list(...)
            captured$perplexity <- args$perplexity
            list(Y = pb_test_fake_embedding(X, args$dims))
        },
        .package = "Rtsne"
    )

    expect_error(
        plot_TSNE(
            fixture$matrix,
            fixture$annotation,
            tsne_dims = 1
        ),
        "tsne_dims.*2"
    )
    expect_error(
        plot_TSNE(
            fixture$matrix[, 1L, drop = FALSE],
            fixture$annotation[fixture$annotation$FullRunName == "sample_1", ],
            color_by = "group",
            color_scheme = c(control = "#1B9E77"),
            perplexity = 1
        ),
        "At least two samples"
    )
    expect_warning(
        plot_TSNE(
            fixture$matrix,
            fixture$annotation,
            color_by = "group",
            color_scheme = c(control = "#1B9E77", treated = "#D95F02"),
            perplexity = 999,
            max_iter = 250,
            random_seed = 1
        ),
        "adjusting"
    )
    expect_identical(captured$perplexity, 2)

    expect_error(
        plot_UMAP(
            fixture$matrix,
            fixture$annotation,
            n_components = 1
        ),
        "n_components.*2"
    )
    expect_error(
        plot_UMAP(
            fixture$matrix,
            fixture$annotation,
            color_by = "group",
            n_neighbors = ncol(fixture$matrix)
        ),
        "n_neighbors.*smaller than the number of samples"
    )

    duplicate_annotation <- rbind(
        fixture$annotation,
        fixture$annotation[1L, , drop = FALSE]
    )
    expect_error(
        plot_UMAP(
            fixture$matrix,
            duplicate_annotation,
            color_by = "group",
            n_neighbors = 3
        ),
        "uplicat"
    )
})

test_that("explicit embedding seeds reproduce backend coordinates", {
    skip_if_not_installed("Rtsne")
    skip_if_not_installed("umap")

    fixture <- pb_test_embedding_fixture(n_features = 10L, n_samples = 8L)
    common_tsne <- list(
        data_matrix = fixture$matrix,
        sample_annotation = fixture$annotation,
        color_by = "group",
        color_scheme = c(control = "#1B9E77", treated = "#D95F02"),
        perplexity = 2,
        max_iter = 250,
        random_seed = 2026
    )
    tsne_a <- do.call(plot_TSNE.default, common_tsne)
    tsne_b <- do.call(plot_TSNE.default, common_tsne)
    expect_equal(
        tsne_a$data[, c("Dim1", "Dim2")],
        tsne_b$data[, c("Dim1", "Dim2")],
        tolerance = 0
    )

    common_umap <- list(
        data_matrix = fixture$matrix,
        sample_annotation = fixture$annotation,
        color_by = "group",
        color_scheme = c(control = "#1B9E77", treated = "#D95F02"),
        n_neighbors = 3,
        random_state = 2026
    )
    umap_a <- do.call(plot_UMAP.default, common_umap)
    umap_b <- do.call(plot_UMAP.default, common_umap)
    expect_equal(
        umap_a$data[, c("Dim1", "Dim2")],
        umap_b$data[, c("Dim1", "Dim2")],
        tolerance = 0
    )
})

test_that("matrix embeddings return plotly objects only when requested", {
    skip_if_not_installed("Rtsne")
    skip_if_not_installed("umap")
    skip_if_not_installed("plotly")

    fixture <- pb_test_embedding_fixture()
    testthat::local_mocked_bindings(
        Rtsne = function(X, ...) list(Y = pb_test_fake_embedding(X)),
        .package = "Rtsne"
    )
    testthat::local_mocked_bindings(
        umap = function(d, config, ...) {
            list(layout = pb_test_fake_embedding(d, config$n_components))
        },
        .package = "umap"
    )

    tsne <- plot_TSNE(
        fixture$matrix,
        fixture$annotation,
        color_by = "group",
        shape_by = "marker",
        color_scheme = c(control = "#1B9E77", treated = "#D95F02"),
        perplexity = 2,
        use_plotlyrender = TRUE
    )
    umap <- plot_UMAP(
        fixture$matrix,
        fixture$annotation,
        color_by = "score",
        shape_by = "marker",
        color_scheme = c("#2166AC", "#B2182B"),
        n_neighbors = 3,
        use_plotlyrender = TRUE
    )

    expect_s3_class(tsne, "plotly")
    expect_s3_class(umap, "plotly")
    expect_gt(length(plotly::plotly_build(tsne)$x$data), 0L)
    expect_gt(length(plotly::plotly_build(umap)$x$data), 0L)

    reserved_annotation <- fixture$aligned_annotation
    reserved_annotation$.color_value <- factor(
        rep(c("round", "angular", "angular", "round"),
            length.out = nrow(reserved_annotation)
        )
    )
    reserved_annotation$.shape_value <- seq_len(nrow(reserved_annotation))
    reserved_annotation$.hover_text <- paste0(
        "annotation_", seq_len(nrow(reserved_annotation))
    )
    captured <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        plot_ly = function(...) {
            captured$args <- list(...)
            structure(list(), class = c("plotly", "htmlwidget"))
        },
        layout = function(p, ...) p,
        .package = "plotly"
    )
    proBatch:::.pb_create_embedding_plotly(
        embedding_matrix = pb_test_fake_embedding(t(fixture$matrix)),
        sample_ids = colnames(fixture$matrix),
        sample_annotation = reserved_annotation,
        sample_id_col = "FullRunName",
        color_by = "group",
        shape_by = ".color_value",
        color_scheme = c(
            control = "#1B9E77",
            treated = "#D95F02"
        ),
        point_size = 3,
        point_alpha = 0.8,
        plot_title = NULL,
        axis_labels = list(x = "Dim1", y = "Dim2", title = "Embedding"),
        plotly_param = list()
    )

    expect_identical(
        captured$args$data$.color_value,
        reserved_annotation$.color_value
    )
    expect_identical(
        captured$args$data$.shape_value,
        reserved_annotation$.shape_value
    )
    expect_identical(
        captured$args$data$.hover_text,
        reserved_annotation$.hover_text
    )
    internal_columns <- vapply(
        captured$args[c("color", "symbol", "text")],
        function(mapping) all.vars(mapping)[[1L]],
        character(1L)
    )
    expect_false(any(internal_columns %in%
        c(".color_value", ".shape_value", ".hover_text")))
    expect_identical(
        as.character(captured$args$data[[internal_columns[["symbol"]]]]),
        as.character(reserved_annotation$.color_value)
    )
})

test_that("ProBatchFeatures t-SNE preserves assay order and per-assay arguments", {
    skip_if_not_installed("Rtsne")
    skip_if_not_installed("gridExtra")

    pbf <- pb_test_make_pbf(
        n_rows = 12L,
        n_cols = 8L,
        add_log2 = TRUE,
        complete = TRUE
    )
    assays <- rev(names(pbf))
    color_annotation <- as.data.frame(
        SummarizedExperiment::colData(pbf)
    )[, c("FullRunName", "MS_batch"), drop = FALSE]
    batch_colors <- sample_annotation_to_colors(
        color_annotation,
        sample_id_col = "FullRunName",
        factor_columns = "MS_batch"
    )[["MS_batch"]]
    captured <- new.env(parent = emptyenv())
    captured$calls <- list()
    testthat::local_mocked_bindings(
        Rtsne = function(X, ...) {
            args <- list(...)
            captured$calls[[length(captured$calls) + 1L]] <- list(
                input = X,
                args = args
            )
            list(Y = pb_test_fake_embedding(X, args$dims))
        },
        .package = "Rtsne"
    )

    assay_args <- setNames(list(
        list(
            perplexity = 2,
            max_iter = 251,
            random_seed = 101,
            theta = 0.11
        ),
        list(
            perplexity = 2,
            max_iter = 252,
            random_seed = 202,
            theta = 0.22
        )
    ), assays)
    assay_titles <- setNames(c("Requested first", "Requested second"), assays)

    expect_no_warning(
        static <- plot_TSNE(
            data_matrix = pbf,
            pbf_name = assays,
            sample_id_col = "FullRunName",
            plot_title = assay_titles,
            color_by = "MS_batch",
            color_scheme = batch_colors,
            assay_args = assay_args,
            return_gridExtra = TRUE
        )
    )

    expect_type(static, "list")
    expect_named(static, c("grob", "plots"))
    expect_named(static$plots, assays)
    expect_true(all(vapply(static$plots, inherits, logical(1), "ggplot")))
    expect_identical(
        unname(vapply(
            static$plots,
            function(plot) plot$labels$title,
            character(1)
        )),
        unname(assay_titles)
    )
    expect_equal(
        vapply(captured$calls, function(call) call$args$theta, numeric(1)),
        c(0.11, 0.22)
    )
    expect_equal(
        vapply(captured$calls, function(call) call$args$max_iter, numeric(1)),
        c(251, 252)
    )

    expect_no_warning(
        single <- plot_TSNE(
            pbf,
            pbf_name = assays[[1L]],
            sample_id_col = "FullRunName",
            color_by = "MS_batch",
            color_scheme = batch_colors,
            assay_args = assay_args[assays[[1L]]]
        )
    )
    expect_s3_class(single, "ggplot")
})

test_that("ProBatchFeatures UMAP works without per-assay overrides", {
    skip_if_not_installed("umap")

    pbf <- pb_test_make_pbf(
        n_rows = 12L,
        n_cols = 8L,
        complete = TRUE
    )
    assay <- names(pbf)[[1L]]
    color_annotation <- as.data.frame(
        SummarizedExperiment::colData(pbf)
    )[, c("FullRunName", "MS_batch"), drop = FALSE]
    batch_colors <- sample_annotation_to_colors(
        color_annotation,
        sample_id_col = "FullRunName",
        factor_columns = "MS_batch"
    )[["MS_batch"]]
    testthat::local_mocked_bindings(
        umap = function(d, config, ...) {
            list(layout = pb_test_fake_embedding(d, config$n_components))
        },
        .package = "umap"
    )

    expect_no_warning(
        umap_plot <- plot_UMAP(
            pbf,
            pbf_name = assay,
            sample_id_col = "FullRunName",
            color_by = "MS_batch",
            color_scheme = batch_colors,
            n_neighbors = 3,
            random_state = 2026
        )
    )
    expect_s3_class(umap_plot, "ggplot")
})

test_that("ProBatchFeatures interactive embeddings return lists or subplots", {
    skip_if_not_installed("Rtsne")
    skip_if_not_installed("umap")
    skip_if_not_installed("plotly")

    pbf <- pb_test_make_pbf(
        n_rows = 12L,
        n_cols = 8L,
        add_log2 = TRUE,
        complete = TRUE
    )
    assays <- rev(names(pbf))
    testthat::local_mocked_bindings(
        Rtsne = function(X, ...) {
            args <- list(...)
            list(Y = pb_test_fake_embedding(X, args$dims))
        },
        .package = "Rtsne"
    )
    testthat::local_mocked_bindings(
        umap = function(d, config, ...) {
            list(layout = pb_test_fake_embedding(d, config$n_components))
        },
        .package = "umap"
    )

    tsne_args <- setNames(lapply(seq_along(assays), function(index) {
        list(
            perplexity = 2,
            max_iter = 250 + index,
            random_seed = 300 + index
        )
    }), assays)
    expect_no_warning(
        tsne_list <- plot_TSNE(
            pbf,
            pbf_name = assays,
            sample_id_col = "FullRunName",
            color_by = "MS_batch",
            assay_args = tsne_args,
            use_plotlyrender = TRUE
        )
    )

    expect_type(tsne_list, "list")
    expect_named(tsne_list, assays)
    expect_true(all(vapply(tsne_list, inherits, logical(1), "plotly")))

    umap_args <- setNames(lapply(seq_along(assays), function(index) {
        list(
            n_neighbors = 2 + index,
            random_state = 400 + index,
            preserve.seed = index == 1L
        )
    }), assays)
    expect_no_warning(
        umap_list <- plot_UMAP(
            data_matrix = pbf,
            pbf_name = assays,
            sample_id_col = "FullRunName",
            color_by = "MS_batch",
            assay_args = umap_args,
            use_plotlyrender = TRUE
        )
    )

    expect_type(umap_list, "list")
    expect_named(umap_list, assays)
    expect_true(all(vapply(umap_list, inherits, logical(1), "plotly")))

    expect_no_warning(
        subplot <- plot_UMAP(
            pbf,
            pbf_name = assays,
            sample_id_col = "FullRunName",
            color_by = "MS_batch",
            assay_args = umap_args,
            use_plotlyrender = TRUE,
            return_subplots = TRUE,
            subplot_ncol = 1L,
            share_axes = c(x = TRUE, y = FALSE)
        )
    )
    subplot_build <- plotly::plotly_build(subplot)

    expect_s3_class(subplot, "plotly")
    expect_gte(length(subplot_build$x$data), length(assays))
})

test_that("ProBatchFeatures embedding collection controls are validated", {
    skip_if_not_installed("Rtsne")

    pbf <- pb_test_make_pbf(
        n_rows = 12L,
        n_cols = 8L,
        add_log2 = TRUE,
        complete = TRUE
    )
    assays <- names(pbf)
    color_annotation <- as.data.frame(
        SummarizedExperiment::colData(pbf)
    )[, c("FullRunName", "MS_batch"), drop = FALSE]
    batch_colors <- sample_annotation_to_colors(
        color_annotation,
        sample_id_col = "FullRunName",
        factor_columns = "MS_batch"
    )[["MS_batch"]]
    backend_calls <- 0L
    testthat::local_mocked_bindings(
        Rtsne = function(X, ...) {
            backend_calls <<- backend_calls + 1L
            args <- list(...)
            list(Y = pb_test_fake_embedding(X, args$dims))
        },
        .package = "Rtsne"
    )

    expect_error(
        plot_TSNE(
            pbf,
            pbf_name = assays,
            sample_id_col = "FullRunName",
            color_by = "MS_batch",
            assay_args = list(
                not_an_assay = list(perplexity = 2)
            )
        ),
        "unknown assay"
    )

    reserved_args <- setNames(
        list(list(sample_id_col = "not_allowed")),
        assays[[1L]]
    )
    expect_error(
        plot_TSNE(
            pbf,
            pbf_name = assays,
            sample_id_col = "FullRunName",
            color_by = "MS_batch",
            assay_args = reserved_args
        ),
        "cannot override reserved"
    )

    valid_args <- setNames(lapply(assays, function(assay) {
        list(perplexity = 2, max_iter = 250, random_seed = 9)
    }), assays)
    expect_error(
        plot_TSNE(
            pbf,
            pbf_name = assays,
            sample_id_col = "FullRunName",
            color_by = "MS_batch",
            color_scheme = batch_colors,
            assay_args = valid_args,
            subplot_ncol = 0
        ),
        "subplot_ncol.*positive integer"
    )
    expect_error(
        plot_TSNE(
            pbf,
            pbf_name = assays,
            sample_id_col = "FullRunName",
            color_by = "MS_batch",
            color_scheme = batch_colors,
            assay_args = valid_args,
            share_axes = NA
        ),
        "share_axes.*logical"
    )
    expect_identical(backend_calls, 0L)
})

test_that("plot_PCA ProBatchFeatures arranges assays and preserves subset order", {
    skip_if_not_installed("gridExtra")
    pbf <- pb_test_make_pbf(n_rows = 40, n_cols = 6, add_log2 = TRUE)

    res <- suppressWarnings(plot_PCA(
        pbf,
        sample_id_col = "FullRunName",
        return_gridExtra = TRUE
    ))

    expect_type(res, "list")
    expect_named(res$plots, names(pbf))
    expect_equal(length(res$plots), length(names(pbf)))
    expect_true(all(vapply(res$plots, inherits, logical(1), "ggplot")))

    subset_assays <- rev(names(pbf))
    res_subset <- suppressWarnings(plot_PCA(
        pbf,
        pbf_name = subset_assays,
        sample_id_col = "FullRunName",
        return_gridExtra = TRUE
    ))

    expect_type(res_subset, "list")
    expect_equal(names(res_subset$plots), subset_assays)
    expect_true(all(vapply(res_subset$plots, inherits, logical(1), "ggplot")))
})

test_that("plot_heatmap_diagnostic ProBatchFeatures arranges assays and preserves subset order", {
    skip_if_not_installed("gridExtra")
    pbf <- pb_test_make_pbf(n_rows = 30, n_cols = 5, add_log2 = TRUE)

    res <- suppressWarnings(plot_heatmap_diagnostic(
        pbf,
        sample_id_col = "FullRunName",
        factors_to_plot = c("MS_batch"),
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        return_gridExtra = TRUE
    ))

    expect_type(res, "list")
    expect_equal(length(res$plots), length(names(pbf)))
    expect_true(all(vapply(res$plots, function(x) inherits(x, "pheatmap"), logical(1))))

    subset_assays <- names(pbf)[2:1]
    res_subset <- suppressWarnings(plot_heatmap_diagnostic(
        pbf,
        pbf_name = subset_assays,
        sample_id_col = "FullRunName",
        factors_to_plot = c("MS_batch"),
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        return_gridExtra = TRUE
    ))

    expect_type(res_subset, "list")
    expect_equal(names(res_subset$plots), subset_assays)
    expect_true(all(vapply(res_subset$plots, function(x) inherits(x, "pheatmap"), logical(1))))
})

test_that("plot_PCA ProBatchFeatures returns ggplot for single assay", {
    pbf <- pb_test_make_pbf(n_rows = 40, n_cols = 6, add_log2 = TRUE)

    single_assay <- names(pbf)[1]
    res <- suppressWarnings(plot_PCA(
        pbf,
        pbf_name = single_assay,
        sample_id_col = "FullRunName",
        x_nPC = 2,
        y_nPC = 3
    ))

    expect_s3_class(res, "ggplot")
    expect_match(res$labels$x, "^PC2 ")
    expect_match(res$labels$y, "^PC3 ")
})

test_that("plot_heatmap_diagnostic ProBatchFeatures single assay returns pheatmap", {
    pbf <- pb_test_make_pbf(n_rows = 30, n_cols = 5, add_log2 = FALSE)

    res <- suppressWarnings(plot_heatmap_diagnostic(
        pbf,
        pbf_name = names(pbf)[1],
        sample_id_col = "FullRunName",
        factors_to_plot = c("MS_batch"),
        cluster_rows = FALSE,
        cluster_cols = FALSE
    ))

    expect_true(inherits(res, "pheatmap"))
})

test_that("PVCA ownership and historical df registration stay consolidated", {
    check_root <- normalizePath(
        file.path(testthat::test_path(), "..", ".."),
        mustWork = TRUE
    )
    candidates <- c(
        check_root,
        file.path(check_root, "00_pkg_src", "proBatch")
    )
    valid <- vapply(
        candidates,
        function(candidate) {
            file.exists(file.path(candidate, "DESCRIPTION")) &&
                length(list.files(
                    file.path(candidate, "R"),
                    pattern = "[.]R$"
                )) > 0L
        },
        logical(1L)
    )
    if (!any(valid)) {
        skip("PVCA ownership checks require the proBatch source package")
    }

    source_root <- normalizePath(
        candidates[[which(valid)[[1L]]]],
        mustWork = TRUE
    )
    r_dir <- file.path(source_root, "R")
    definitions <- do.call(rbind, lapply(
        sort(list.files(r_dir, pattern = "[.]R$", full.names = TRUE)),
        function(file) {
            expressions <- parse(file = file)
            rows <- lapply(expressions, function(expression) {
                is_assignment <- is.call(expression) &&
                    identical(expression[[1L]], as.name("<-"))
                is_function <- is_assignment &&
                    is.call(expression[[3L]]) &&
                    identical(expression[[3L]][[1L]], as.name("function"))
                if (!is_function) {
                    return(NULL)
                }
                data.frame(
                    symbol = paste(deparse(expression[[2L]]), collapse = ""),
                    file = basename(file),
                    stringsAsFactors = FALSE
                )
            })
            do.call(rbind, rows)
        }
    ))

    public_symbols <- c(
        "calculate_PVCA",
        "calculate_PVCA.default",
        "calculate_PVCA.ProBatchFeatures",
        "plot_PVCA",
        "plot_PVCA.default",
        "plot_PVCA.ProBatchFeatures",
        "plot_PVCA.df",
        "prepare_PVCA_df",
        "prepare_PVCA_df.default",
        "prepare_PVCA_df.ProBatchFeatures"
    )
    owned <- definitions[definitions$symbol %in% public_symbols, ]

    expect_setequal(owned$symbol, public_symbols)
    expect_equal(nrow(owned), length(public_symbols))
    expect_true(all(owned$file == "proteome_wide_diagnostics.R"))
    expect_false(any(grepl("[.]df[.]default$", definitions$symbol)))

    source <- readLines(
        file.path(r_dir, "proteome_wide_diagnostics.R"),
        warn = FALSE
    )
    definition <- grep("^plot_PVCA[.]df <- function", source)
    expect_length(definition, 1L)
    cursor <- definition - 1L
    roxygen <- character()
    while (cursor > 0L && startsWith(source[[cursor]], "#'")) {
        roxygen <- c(source[[cursor]], roxygen)
        cursor <- cursor - 1L
    }

    expect_true(any(trimws(roxygen) == "#' @method plot_PVCA df"))
    expect_true(any(trimws(roxygen) == "#' @export"))
    expect_true(any(grepl(
        "@rawNamespace export(plot_PVCA.df)",
        roxygen,
        fixed = TRUE
    )))
    expect_true(grepl(
        "registered for the historical",
        paste(roxygen, collapse = " "),
        fixed = TRUE
    ))
})
