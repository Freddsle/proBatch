pb_test_load_example_data <- function() {
    data(
        list = c("example_proteome_matrix", "example_sample_annotation"),
        package = "proBatch",
        envir = parent.frame()
    )
    invisible(TRUE)
}

pb_test_make_pbf <- function(n_rows = 30, n_cols = 6, add_log2 = FALSE) {
    pb_test_load_example_data()

    matrix_small <- example_proteome_matrix[1:n_rows, 1:n_cols]
    sample_ids <- colnames(matrix_small)
    sample_ann <- example_sample_annotation[match(sample_ids, example_sample_annotation$FullRunName), ]

    pbf <- suppressMessages(ProBatchFeatures(
        data_matrix = matrix_small,
        sample_annotation = sample_ann,
        sample_id_col = "FullRunName",
        name = "feature::raw"
    ))

    if (isTRUE(add_log2)) {
        pbf <- suppressMessages(pb_transform(
            pbf,
            from = "feature::raw",
            steps = "log2",
            store_fast_steps = TRUE
        ))
    }

    pbf
}
