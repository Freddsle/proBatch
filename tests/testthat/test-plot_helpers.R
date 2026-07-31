test_that(".pb_split_arg_by_assay keeps atomic vectors intact", {
    assays <- paste0("assay", 1:4)
    vec <- c("tech1", "tech2", "tech3")
    res <- proBatch:::`.pb_split_arg_by_assay`(vec, assays)

    expect_equal(length(res), length(assays))
    expect_true(all(vapply(res, identical, logical(1), vec)))
})

test_that(".pb_split_arg_by_assay can map atomic vectors per assay", {
    assays <- paste0("assay", 1:4)
    vec <- c("title1", "title2", "title3")
    res <- proBatch:::`.pb_split_arg_by_assay`(
        vec,
        assays,
        atomic_vector_mode = "per_assay"
    )

    expect_equal(res[[1]], "title1")
    expect_equal(res[[2]], "title2")
    expect_equal(res[[3]], "title3")
    expect_equal(res[[4]], "title3")
    expect_equal(
        proBatch:::`.pb_resolve_titles`(assays, vec),
        c("title1", "title2", "title3", "title3")
    )
})

test_that(".pb_split_arg_by_assay honours named vectors", {
    assays <- c("assay1", "assay2", "assay3")
    vec <- c(assay1 = "x", assay3 = "z")
    res <- proBatch:::`.pb_split_arg_by_assay`(vec, assays)

    expect_equal(res[[1]], "x")
    expect_equal(res[[2]], "z")
    expect_equal(res[[3]], "z")
})

test_that(".pb_prepare_embedding_assay_args supplies empty defaults", {
    assays <- c("assay1", "assay2")

    empty <- proBatch:::`.pb_prepare_embedding_assay_args`(NULL, assays)
    expect_identical(empty, setNames(list(list(), list()), assays))

    partial <- proBatch:::`.pb_prepare_embedding_assay_args`(
        list(assay1 = list(random_state = 7L)),
        assays
    )
    expect_identical(partial[["assay1"]], list(random_state = 7L))
    expect_identical(partial[["assay2"]], list())
})

test_that(".pb_resolve_assay_for_input validates NA assay names", {
    object <- structure(list(), class = "ProBatchFeatures")

    expect_error(
        proBatch:::`.pb_resolve_assay_for_input`(object, pbf_name = NA_character_),
        "`pbf_name` must contain exactly one non-empty assay name."
    )
})

test_that(".pb_resolve_assays_for_input validates and deduplicates selections", {
    object <- structure(list(assay1 = 1, assay2 = 2), class = "ProBatchFeatures")

    expect_equal(
        proBatch:::`.pb_resolve_assays_for_input`(
            object = object,
            pbf_name = c("assay1", "assay1", "assay2"),
            deduplicate = TRUE
        ),
        c("assay1", "assay2")
    )
    expect_error(
        proBatch:::`.pb_resolve_assays_for_input`(
            object = object,
            pbf_name = c("assay1", "", NA_character_)
        ),
        "non-missing, non-empty assay names"
    )
    expect_error(
        proBatch:::`.pb_resolve_assays_for_input`(
            object = object,
            pbf_name = character()
        ),
        "non-missing, non-empty assay names"
    )

    expect_equal(
        proBatch:::`.pb_resolve_assays_for_input`(
            object = object,
            pbf_name = NULL,
            default = "all",
            deduplicate = TRUE
        ),
        c("assay1", "assay2")
    )
})

test_that(".pb_handle_missing_wrapper preserves diagnostic legacy policies", {
    symmetric <- diag(3)
    dimnames(symmetric) <- list(paste0("f", 1:3), paste0("f", 1:3))
    symmetric[1, 3] <- NA_real_
    symmetric[3, 1] <- NA_real_

    removed <- pb_test_expect_warnings(
        suppressMessages(proBatch:::.pb_handle_missing_wrapper(
            data_matrix = symmetric,
            warning_message = "diagnostic requires complete data",
            fill_the_missing = NULL
        )),
        c(
            "diagnostic requires complete data",
            "removed 2 rows and 2 columns"
        ),
        fixed = TRUE
    )
    expect_identical(dim(removed), c(1L, 1L))
    expect_identical(rownames(removed), "f2")
    expect_identical(colnames(removed), "f2")

    filled <- pb_test_expect_warnings(
        suppressMessages(proBatch:::.pb_handle_missing_wrapper(
            data_matrix = matrix(c(1, NA_real_), nrow = 1),
            warning_message = "diagnostic requires complete data",
            fill_the_missing = 0
        )),
        c(
            "diagnostic requires complete data",
            "filling missing values with 0"
        ),
        fixed = TRUE
    )
    expect_identical(unname(filled), matrix(c(1, 0), nrow = 1))
})

test_that(".pb_align_matrix_and_annotation merges check arguments", {
    data_matrix <- matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s2", "s1"))
    )
    annotation <- data.frame(
        FullRunName = c("s1", "s2"),
        Batch = c("A", "B"),
        stringsAsFactors = FALSE
    )

    aligned <- proBatch:::.pb_align_matrix_and_annotation(
        data_matrix,
        annotation,
        sample_id_col = "FullRunName",
        check_args = list(batch_col = "Batch"),
        allow_partial_annotation = FALSE
    )

    expect_equal(aligned$data_matrix, data_matrix)
    expect_identical(aligned$sample_ids, c("s2", "s1"))
    expect_identical(
        as.character(aligned$sample_annotation$FullRunName),
        aligned$sample_ids
    )
})
