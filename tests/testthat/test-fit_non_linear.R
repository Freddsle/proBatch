test_that("fit works", {
    data(example_proteome, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    test_annotation <- example_sample_annotation[example_sample_annotation$MS_batch == "Batch_1", ]
    selected_files <- test_annotation$FullRunName

    df_selected <- example_proteome[example_proteome$peptide_group_label == example_proteome$peptide_group_label[1], ]
    df_selected <- df_selected[df_selected$FullRunName %in% selected_files, ]
    df_selected <- merge(df_selected, test_annotation, by = "FullRunName")

    fit_values <- fit_nonlinear(df_selected)

    # expect_length(fit_values, nrow(df_selected))
})

test_that("rle_func finds longest measured run", {
    df <- data.frame(order = 1:6, Intensity = c(1, 2, 3, NA, 4, 5))
    expect_equal(rle_func(df), 3)
 })

 test_that("rle_func handles custom columns and all NA", {
    df <- data.frame(idx = 6:1, value = c(NA, NA, 3, 4, 5, NA))
    expect_equal(rle_func(df, measure_col = "value", order_col = "idx"), 3)
    df_na <- data.frame(idx = 1:3, value = c(NA, NA, NA))
    expect_equal(rle_func(df_na, measure_col = "value", order_col = "idx"), 0)
 })

 # new tests for fit_nonlinear behaviour
 test_that("fit_nonlinear returns fitted curve with sufficient data", {
    set.seed(1)
    df <- data.frame(order = 1:10,
                     Intensity = sin(1:10),
                     m_score = rep(0, 10))
    vals <- fit_nonlinear(df, min_measurements = 5, qual_col = NULL, no_fit_imputed = FALSE)
    expect_length(vals, nrow(df))
    expect_true(any(!is.na(vals)))
 })

test_that("fit_nonlinear returns NA vector when insufficient data", {
    df <- data.frame(order = 1:4,
                     Intensity = c(1, NA, NA, 2),
                     m_score = c(0, 2, 2, 0))
    vals <- suppressWarnings(
        fit_nonlinear(df, min_measurements = 5, qual_col = NULL, no_fit_imputed = FALSE)
    )
    expect_true(all(is.na(vals)))
})