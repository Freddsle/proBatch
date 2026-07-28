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

    pb_test_expect_warnings(
        pca <- plot_PCA(
            example_proteome_matrix, example_sample_annotation,
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

test_that("plot_PCA supports marginal density plots", {
    pb_test_load_example_data()

    pca <- suppressMessages(suppressWarnings(plot_PCA(
        example_proteome_matrix, example_sample_annotation,
        color_by = "MS_batch",
        fill_the_missing = -1,
        marginal_density = TRUE
    )))

    expect_true(inherits(pca, "ggplot") || grid::is.grob(pca))
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
