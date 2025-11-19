
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

