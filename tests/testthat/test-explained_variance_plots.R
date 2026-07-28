test_that("pvca_plot", {
    pb_test_load_example_data()

    matrix_test <- example_proteome_matrix[1:150, ]
    pb_test_expect_warnings(
        pvca <- plot_PVCA(matrix_test, example_sample_annotation,
            technical_factors = c("MS_batch", "digestion_batch"),
            biological_factors = c("Diet", "Sex", "Strain")
        ),
        c(
            "PVCA cannot operate with missing values in the matrix",
            "filling missing values with -1"
        ),
        fixed = TRUE
    )

    expect_equal(pvca$data$weights[1], 0.39166175, tolerance = 3e-2, ignore_attr = TRUE)
    expect_equal(as.character(pvca$data$label[3]), "MS_batch", ignore_attr = TRUE)
    expect_equal(as.character(pvca$data$label[2]), "Sex:Strain", ignore_attr = TRUE)

    expect_equal(pvca$data$category[1], "biological")
    expect_equal(pvca$data$category[3], "technical")
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
    expect_setequal(as.character(unique(stacked$data$assay)), c("assay_a", "assay_b"))
    expect_equal(levels(stacked$data$assay), c("assay_b", "assay_a"))
})
