# Helper: small tolerance for numeric comparisons
tol <- 1e-8


test_that("log_transformed_matrix", {
    data(example_proteome_matrix, package = "proBatch")

    matrix_test <- example_proteome_matrix[1:10, ]

    # log2, offset 0.5
    log2_transformed_matrix <- log_transform_dm(matrix_test, log_base = 2, offset = 0.5)
    manual_log2 <- log(matrix_test + 0.5, base = 2)
    expect_equal(manual_log2, log2_transformed_matrix, ignore_attr = TRUE)

    # log10, offset 1
    log10_transformed_matrix <- log_transform_dm(matrix_test, log_base = 10, offset = 1)
    manual_log10 <- log(matrix_test + 1, base = 10)
    expect_equal(manual_log10, log10_transformed_matrix, ignore_attr = TRUE)

    # Warning when log_base NULL
    expect_warning(
        {
            log_warn <- log_transform_dm(matrix_test, log_base = NULL)
            expect_identical(log_warn, matrix_test)
        },
        regexp = "Log base is NULL"
    )
})


test_that("log_transform_dm: numeric matrix transforms correctly", {
    # Construct a small numeric matrix, including zeros and positive values
    mat <- matrix(c(0, 1, 2.5, 10), nrow = 2, ncol = 2)

    # Case 1: log2 with offset 0.5
    offset <- 0.5
    base <- 2
    expected_log2 <- log(mat + offset, base = base)
    result_log2 <- log_transform_dm(mat, log_base = base, offset = offset)
    expect_equal(result_log2, expected_log2)

    # Case 2: log10 with offset 1
    offset2 <- 1
    base2 <- 10
    expected_log10 <- log(mat + offset2, base = base2)
    result_log10 <- log_transform_dm(mat, log_base = base2, offset = offset2)
    expect_equal(result_log10, expected_log10)

    # Round-trip: log → unlog approximates original
    unlogged2 <- unlog_dm(result_log2, log_base = base, offset = offset)
    # Because of floating arithmetic, allow small tolerance
    expect_equal(unlogged2, mat, tolerance = tol)

    unlogged10 <- unlog_dm(result_log10, log_base = base2, offset = offset2)
    expect_equal(unlogged10, mat, tolerance = tol)
})

test_that("log_transform_dm: error on non-numeric matrix", {
    mat_char <- matrix(c("a", "b", "c", "d"), nrow = 2)
    expect_error(log_transform_dm(mat_char, log_base = 2, offset = 1),
        regexp = "data_matrix must be numeric"
    )
})

test_that("log_transform_dm and unlog_dm: warning when log_base is NULL", {
    mat <- matrix(1:4, nrow = 2)
    expect_warning(
        {
            out <- log_transform_dm(mat, log_base = NULL, offset = 1)
            # Should return original matrix
            expect_identical(out, mat)
        },
        regexp = "Log base is NULL"
    )

    expect_warning(
        {
            out2 <- unlog_dm(mat, log_base = NULL, offset = 1)
            expect_identical(out2, mat)
        },
        regexp = "Log base is NULL"
    )
})

test_that("log_transform_df: data frame long format transforms correctly", {
    # Create a simple tibble with an "Intensity" column
    df_long <- tibble(
        Sample = rep(c("A", "B"), each = 3),
        Intensity = c(0, 1, 4, 0.2, 2.5, 10)
    )

    # Case: log_base = 2, offset = 0.5
    base <- 2
    offset <- 0.5
    df_logged <- log_transform_df(df_long, log_base = base, offset = offset, measure_col = "Intensity")

    # Expect a new column "beforeLog_Intensity" equal to original
    expect_true("beforeLog_Intensity" %in% names(df_logged))
    expect_equal(df_logged$beforeLog_Intensity, df_long$Intensity)

    # Expect transformed column equals log(Intensity + offset, base = 2)
    expect_equal(df_logged$Intensity, log(df_long$Intensity + offset, base = base))

    # Round-trip via unlog_df
    df_unlogged <- unlog_df(df_logged, log_base = base, offset = offset, measure_col = "Intensity")
    # Should add "beforeUnLog_Intensity"
    expect_true("beforeUnLog_Intensity" %in% names(df_unlogged))
    # The “beforeUnLog_Intensity” should equal the logged values
    expect_equal(df_unlogged$beforeUnLog_Intensity, df_logged$Intensity)
    # The final Intensity should approximate original df_long$Intensity
    expect_equal(df_unlogged$Intensity, df_long$Intensity, tolerance = tol)
})

test_that("unlog_df: handles negative or zero after unlog (informative check)", {
    # Construct a logged value that yields small or negative on unlog
    # For example, if logged is very negative: base^(-10) - offset → nearly zero minus offset → negative
    df_long <- tibble(
        Sample = "X",
        Intensity = -10 # logged value
    )
    base <- 2
    offset <- 1
    df_unlogged <- unlog_df(df_long, log_base = base, offset = offset, measure_col = "Intensity")
    # Compute manually: 2^(-10) - 1
    expected_val <- 2^(-10) - 1
    expect_equal(df_unlogged$Intensity, expected_val)
    # This may be negative; test that function does not error
    expect_true(is.numeric(df_unlogged$Intensity))
})

test_that("log_transform_df produces expected values", {
    data(example_proteome, package = "proBatch")

    df_test <- example_proteome[1:3, ]
    log_df <- log_transform_df(df_test, log_base = 2, offset = 0.5)

    expected <- log2(df_test$Intensity + 0.5)

    expect_true("beforeLog_Intensity" %in% names(log_df))
    expect_equal(log_df$Intensity, expected, ignore_attr = TRUE)
    expect_equal(log_df$beforeLog_Intensity, df_test$Intensity, ignore_attr = TRUE)

    expect_warning(log_transform_df(df_test, log_base = NULL))
})

test_that("unlog functions invert log transformations", {
    data(example_proteome, package = "proBatch")
    data(example_proteome_matrix, package = "proBatch")

    df_test <- example_proteome[1:3, ]
    log_df <- log_transform_df(df_test, log_base = 2, offset = 0)
    unlog_df_res <- unlog_df(log_df, log_base = 2, offset = 0)

    expect_true("beforeUnLog_Intensity" %in% names(unlog_df_res))
    expect_equal(unlog_df_res$Intensity, df_test$Intensity, ignore_attr = TRUE)

    mat_test <- example_proteome_matrix[1:5, 1:4]
    log_mat <- log_transform_dm(mat_test, log_base = 2, offset = 0)
    unlog_mat <- unlog_dm(log_mat, log_base = 2, offset = 0)

    expect_equal(unlog_mat, mat_test, ignore_attr = TRUE)

    expect_warning(unlog_df(log_df, log_base = NULL))
    expect_warning(unlog_dm(log_mat, log_base = NULL))
})
