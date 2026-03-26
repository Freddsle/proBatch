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
})

test_that(".pb_split_arg_by_assay honours named vectors", {
    assays <- c("assay1", "assay2", "assay3")
    vec <- c(assay1 = "x", assay3 = "z")
    res <- proBatch:::`.pb_split_arg_by_assay`(vec, assays)

    expect_equal(res[[1]], "x")
    expect_equal(res[[2]], "z")
    expect_equal(res[[3]], "z")
})

test_that(".pb_resolve_assay_for_input validates NA assay names", {
    object <- structure(list(), class = "ProBatchFeatures")

    expect_error(
        proBatch:::`.pb_resolve_assay_for_input`(object, pbf_name = NA_character_),
        "`pbf_name` must contain exactly one non-empty assay name."
    )
})

test_that(".pb_resolve_assays_for_input supports filtering and deduplication", {
    object <- structure(list(assay1 = 1, assay2 = 2), class = "ProBatchFeatures")

    expect_equal(
        proBatch:::`.pb_resolve_assays_for_input`(
            object = object,
            pbf_name = c("assay1", "", NA, "assay1", "assay2"),
            deduplicate = TRUE
        ),
        c("assay1", "assay2")
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
