test_that("log_transformed_matrix", {
  data(example_proteome_matrix, package = "proBatch")

  matrix_test <- example_proteome_matrix[1:10, ]
  log2_transformed_matrix <- log_transform_dm(matrix_test, log_base = 2, offset = 0.5)
  log10_transformed_matrix <- log_transform_dm(matrix_test, log_base = 10, offset = 1)

  log2_matrix <- log2(matrix_test + 0.5)
  log10_matrix <- log10(matrix_test + 1)

  expect_equal(log2_matrix[1:10], log2_transformed_matrix[1:10], ignore_attr = TRUE)
  expect_equal(log10_matrix[1:10], log10_transformed_matrix[1:10], ignore_attr = TRUE)

  expect_warning(log_warn <- log_transform_dm(matrix_test, log_base = NULL))
})
